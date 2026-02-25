use extendr_api::prelude::*;
use rust_deseq2_core::prelude::*;
use rust_deseq2_core::{deseq_results, run_deseq};
use ndarray::ShapeBuilder;

#[derive(Debug)]
struct DdsWithDesignInfo {
    dds: DESeqDataSet,
    design_info: Option<DesignInfo>,
}

fn arrayview2_to_rmatrix_colmajor(a: ndarray::ArrayView2<'_, f64>) -> RMatrix<f64> {
    let (nrows, ncols) = a.dim();
    let mut v = Vec::with_capacity(nrows * ncols);
    for c in 0..ncols {
        for r in 0..nrows {
            v.push(a[(r, c)]);
        }
    }
    let v = std::rc::Rc::new(v);
    let v2 = v.clone();
    RMatrix::new_matrix(nrows, ncols, move |r, c| v2[c * nrows + r])
}

fn other_err<T: std::fmt::Display>(e: T) -> extendr_api::Error {
    extendr_api::Error::Other(e.to_string())
}

#[extendr]
/// @export
fn deseq_dataset_from_matrix(
    count_data: Vec<f64>,
    nrows: usize,
    ncols: usize,
    gene_ids: Vec<String>,
    sample_ids: Vec<String>,
    condition: Vec<String>,
    design_variable: String,

) -> extendr_api::Result<ExternalPtr<DdsWithDesignInfo>> {
    let counts_arr = ndarray::Array2::from_shape_vec((nrows, ncols).f(), count_data)
        .map_err(other_err)?;
    let counts = CountMatrix::new(counts_arr, gene_ids, sample_ids.clone()).map_err(other_err)?;

    let mut md = SampleMetadata::new(sample_ids);
    md.add_condition(&design_variable, condition).map_err(other_err)?;

    let dds = DESeqDataSet::new(counts, md, &design_variable).map_err(other_err)?;
    Ok(ExternalPtr::new(DdsWithDesignInfo {
        dds,
        design_info: None,
    }))
}

#[extendr]
/// @export
fn deseq(mut dds: ExternalPtr<DdsWithDesignInfo>) -> extendr_api::Result<ExternalPtr<DdsWithDesignInfo>> {
    let design_info = run_deseq(&mut dds.dds).map_err(other_err)?;
    dds.design_info = Some(design_info);
    Ok(dds)
}

#[extendr]
/// @export
fn results(
    dds: ExternalPtr<DdsWithDesignInfo>,
    numerator: String,
    denominator: String,
    alpha: f64,

) -> extendr_api::Result<List> {
    let design_info = dds
        .design_info
        .as_ref()
        .ok_or_else(|| extendr_api::Error::Other("DESeq has not been run yet. Call deseq(dds) before results(dds, ...)".to_string()))?;

    let res = deseq_results(&dds.dds, design_info, &numerator, &denominator, alpha).map_err(other_err)?;

    Ok(list!(
        log2FoldChange = res.log2_fold_changes,
        lfcSE = res.lfc_se,
        stat = res.stat,
        pvalue = res.pvalues,
        padj = res.padj,
        baseMean = res.base_means
    ))
}

#[extendr]
/// @export
fn counts(dds: ExternalPtr<DdsWithDesignInfo>, normalized: bool) -> RMatrix<f64> {
    if normalized {
        if let Some(nc) = dds.dds.normalized_counts() {
            arrayview2_to_rmatrix_colmajor(nc.view())
        } else {
            arrayview2_to_rmatrix_colmajor(dds.dds.counts().counts())
        }
    } else {
        arrayview2_to_rmatrix_colmajor(dds.dds.counts().counts())
    }
}

#[extendr]
/// @export
fn vst(
    dds: ExternalPtr<DdsWithDesignInfo>,
    blind: bool,
    nsub: usize,
) -> extendr_api::Result<ExternalPtr<VstData>> {
    use rust_deseq2_core::transform::{vst, VstMethod};
    
    // Use proper VST from underlying crate
    let vst_result = vst(&dds.dds, VstMethod::Parametric, blind)
        .map_err(|e| extendr_api::Error::Other(e.to_string()))?;
    
    Ok(ExternalPtr::new(VstData {
        data: vst_result.data,
        gene_ids: vst_result.gene_ids,
        sample_ids: vst_result.sample_ids,
    }))
}

#[extendr]
/// @export
fn assay(vst_data: ExternalPtr<VstData>) -> RMatrix<f64> {
    arrayview2_to_rmatrix_colmajor(vst_data.data.view())
}

// Structure to hold VST transformed data
#[derive(Debug)]
struct VstData {
    data: ndarray::Array2<f64>,
    gene_ids: Vec<String>,
    sample_ids: Vec<String>,
}

#[extendr]
/// @export
fn dispersion_estimates(dds: ExternalPtr<DdsWithDesignInfo>) -> Robj {
    let dispersions = dds.dds.dispersions().unwrap();
    // Convert to R vector
    let disp_vec: Vec<f64> = dispersions.to_vec();
    extendr_api::Robj::from(disp_vec)
}

#[extendr]
/// @export
fn dispersion_function(dds: ExternalPtr<DdsWithDesignInfo>) -> List {
    let (asympt_disp, extra_pois) = dds.dds.dispersion_function().unwrap();
    // Return as list to match standard DESeq2 interface
    list!(
        asympt_disp = asympt_disp,
        extra_pois = extra_pois
    )
}

#[extendr]
/// @export
fn size_factors(dds: ExternalPtr<DdsWithDesignInfo>) -> RMatrix<f64> {
    let sf = dds.dds.size_factors().unwrap();
    // Convert 1D array to 2D matrix
    let sf_matrix = sf.view().into_shape((sf.len(), 1)).unwrap();
    arrayview2_to_rmatrix_colmajor(sf_matrix)
}

#[extendr]
/// @export
fn base_mean(dds: ExternalPtr<DdsWithDesignInfo>) -> Robj {
    // Calculate base mean from normalized counts
    let nc = dds.dds.normalized_counts().unwrap();
    // Calculate mean on Axis(1) to get mean of each gene (row)
    let base_means = nc.view().map_axis(ndarray::Axis(1), |view| view.mean().unwrap_or(0.0));
    // Convert to R vector
    let base_vec: Vec<f64> = base_means.to_vec();
    extendr_api::Robj::from(base_vec)
}

extendr_module! {
    mod RustDESeq2;
    fn deseq_dataset_from_matrix;
    fn deseq;
    fn results;
    fn counts;
    fn vst;
    fn assay;
    fn dispersion_estimates;
    fn dispersion_function;
    fn size_factors;
    fn base_mean;
}
