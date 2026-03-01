use extendr_api::prelude::*;
use rust_deseq2_core::{deseq_results, run_deseq};
use ndarray::{Array1, Array2};
use rust_deseq2_core::prelude::*;

// Struct to mimic R DESeqDataSet
#[derive(Debug)]
struct RustDESeqDataSet {
    dds: DESeqDataSet,
    design_info: Option<DesignInfo>,
}

// Struct to mimic R DESeqTransform
#[derive(Debug)]
struct RustDESeqTransform {
    vst: VstResult,
}

fn arrayview_to_rmatrix(
    a: ndarray::ArrayView2<'_, f64>,
    rownames: Vec<String>,
    colnames: Vec<String>
) -> extendr_api::Result<RMatrix<f64>> {
    let (nrows, ncols) = a.dim();
    
    // 1. Create the matrix directly from the view.
    // This closure runs inside R's memory allocation step.
    // Note: R matrices are column-major, so we map (r, c) carefully.
    let mut matrix = RMatrix::new_matrix(nrows, ncols, |r, c| {
        a[(r, c)]
    });

    // 2. Convert Vec<String> to R objects (Character Vectors)
    let r_rows = Robj::from(rownames);
    let r_cols = Robj::from(colnames);
 
    // 3. Set dimnames (must be a list of 2)
    matrix.set_dimnames(list!(r_rows, r_cols));

    Ok(matrix)
}


#[extendr]
/// @export
fn deseq_dataset_from_matrix(
    count_data: RMatrix<i32>,
    gene_ids: Vec<String>,
    sample_ids: Vec<String>,
    condition: Vec<String>,
    design_variable: String,
) -> extendr_api::Result<ExternalPtr<RustDESeqDataSet>> {
    let [nrows, ncols] = *count_data.dim();

    // Early validation — prevents most misuse
    if gene_ids.len() != nrows {
        return Err(extendr_api::Error::Other(format!(
            "gene_ids length {} ≠ nrows {}", gene_ids.len(), nrows
        )));
    }
    if sample_ids.len() != ncols {
        return Err(extendr_api::Error::Other(format!(
            "sample_ids length {} ≠ ncols {}", sample_ids.len(), ncols
        )));
    }
    if condition.len() != ncols {
        return Err(extendr_api::Error::Other(format!(
            "condition length {} ≠ number of samples {}", condition.len(), ncols
        )));
    }
    if design_variable.is_empty() {
        return Err(extendr_api::Error::Other("design_variable cannot be empty".to_string()));
    }

    // Build Array2 directly by indexing into the RMatrix
    let counts_arr: Array2<f64> = Array2::from_shape_fn((nrows, ncols), |(row, col)| {
        count_data[[row, col]] as f64
    });

    let counts = CountMatrix::new(counts_arr, gene_ids, sample_ids.clone())
        .map_err(|e| extendr_api::Error::Other(format!("CountMatrix creation failed: {}", e)))?;

    let mut md = SampleMetadata::new(sample_ids);
    md.add_condition(&design_variable, condition)
        .map_err(|e| extendr_api::Error::Other(format!("Failed to add condition '{}': {}", design_variable, e)))?;

    let dds = DESeqDataSet::new(counts, md, &design_variable)
        .map_err(|e| extendr_api::Error::Other(format!("DESeqDataSet creation failed: {}", e)))?;

    Ok(ExternalPtr::new(RustDESeqDataSet {
        dds,
        design_info: None,
    }))
}

#[extendr]
/// @export
fn deseq(mut dds: ExternalPtr<RustDESeqDataSet>) -> extendr_api::Result<ExternalPtr<RustDESeqDataSet>> {
    let design_info = run_deseq(&mut dds.dds)
        .map_err(|e| extendr_api::Error::Other(format!("Failed to run DESeq: {}", e)))?;
    
    // modify the struct field
    dds.design_info = Some(design_info);
    
    Ok(dds)
}


#[extendr]
/// @export
fn results(
    dds: ExternalPtr<RustDESeqDataSet>,
    numerator: String,
    denominator: String,
    alpha: f64,
) -> extendr_api::Result<Robj> {
    // Check if design info is available
    let design_info = match &dds.design_info {
        Some(info) => info,
        None => {
            return Err(extendr_api::Error::Other(
                "Design info not found. Please run deseq() first to generate design information.".to_string()
            ));
        }
    };
    
    let res = match deseq_results(
        &dds.dds,
        design_info,
        &numerator,
        &denominator,
        alpha,
    ) {
        Ok(res) => res,
        Err(e) => return Err(extendr_api::Error::Other(format!("Failed to get results: {}", e))),
    };

    // Extract gene IDs to use as row names
    let gene_ids = dds.dds.counts().gene_ids();

    // Build the data frame
    let mut df = data_frame!(
        baseMean = res.base_means,
        log2FoldChange = res.log2_fold_changes,
        lfcSE = res.lfc_se,
        stat = res.stat,
        pvalue = res.pvalues,
        padj = res.padj
    );

    // Set the gene IDs as row names so it matches DESeq2 style
    df.set_attrib("row.names", gene_ids.into_robj())?;

    Ok(df)
}

#[extendr]
/// @export
fn counts(dds: ExternalPtr<RustDESeqDataSet>, normalized: bool) -> extendr_api::Result<RMatrix<f64>> {
    let gene_ids = dds.dds.counts().gene_ids().to_vec();
    let sample_ids = dds.dds.counts().sample_ids().to_vec();
    if normalized {
        arrayview_to_rmatrix(dds.dds.normalized_counts().unwrap().view(), gene_ids, sample_ids)
    } else {
        arrayview_to_rmatrix(dds.dds.counts().counts(), gene_ids, sample_ids)
    }
}

#[extendr]
/// @export
fn vst_transform(
    dds: ExternalPtr<RustDESeqDataSet>,
    blind: bool,
) -> extendr_api::Result<ExternalPtr<RustDESeqTransform>> {
    use rust_deseq2_core::transform::{vst, VstMethod};
    
    // Use proper VST from underlying crate
    let vst_result = match vst(&dds.dds, VstMethod::Parametric, blind) {
        Ok(result) => result,
        Err(e) => return Err(extendr_api::Error::Other(format!("Failed to run VST: {}", e))),
    };
    
    Ok(ExternalPtr::new(RustDESeqTransform {
        vst: vst_result,
    }))
}

#[extendr]
/// @export
fn assay(obj: Robj, assay_name: Option<String>) -> extendr_api::Result<RMatrix<f64>> {
    // 1. Try to treat it as a RustDESeqDataSet
    if let Ok(dds_ptr) = <ExternalPtr<RustDESeqDataSet>>::try_from(obj.clone()) {
        let gene_ids = dds_ptr.dds.counts().gene_ids().to_vec();
        let sample_ids = dds_ptr.dds.counts().sample_ids().to_vec();

        return match assay_name.as_deref() {
            Some("normalized") => {
                arrayview_to_rmatrix(dds_ptr.dds.normalized_counts().unwrap().view(), gene_ids, sample_ids)
            },
            Some("beta") | Some("betaSE") => {
                let design_info = match dds_ptr.design_info.as_ref() {
                    Some(info) => info,
                    None => {
                        return Err(extendr_api::Error::Other("Design info not found. Run deseq() first.".to_string()));
                    }
                };

                let mut colnames = vec!["(Intercept)".to_string()];
                let design_var = dds_ptr.dds.design_variable();
                colnames.extend(design_info.levels.iter().skip(1).map(|l| format!("{}{}", design_var, l)));

                let data = if assay_name.as_deref() == Some("beta") {
                    dds_ptr.dds.coefficients().unwrap()
                } else {
                    dds_ptr.dds.standard_errors().unwrap()
                };
                
                arrayview_to_rmatrix(data.view(), gene_ids, colnames)
            },
            None | Some("counts") => {
                arrayview_to_rmatrix(dds_ptr.dds.counts().counts(), gene_ids, sample_ids)
            },
            _ => Err(extendr_api::Error::Other(format!("Unknown assay: {:?}", assay_name))),
        };
    }

    // 2. Try to treat it as a RustDESeqTransform (VST)
    if let Ok(vst_ptr) = <ExternalPtr<RustDESeqTransform>>::try_from(obj) {
        let rownames = vst_ptr.vst.gene_ids.clone();
        let colnames = vst_ptr.vst.sample_ids.clone();
        return arrayview_to_rmatrix(vst_ptr.vst.data.view(), rownames, colnames);
    }

    Err(extendr_api::Error::Other("Object is neither a DESeqDataSet nor a DESeqTransform".to_string()))
}

#[extendr]
/// @export
fn dispersions(dds: ExternalPtr<RustDESeqDataSet>) -> extendr_api::Result<Robj> {
    let dispersions: &Array1<f64> = dds.dds.dispersions()
        .ok_or_else(|| extendr_api::Error::Other(
            "Dispersions not found. Run deseq() or estimateDispersions first.".to_string()
        ))?;

    let disp_vec: Vec<f64> = dispersions.to_vec();

    let gene_ids: Vec<String> = dds.dds.counts().gene_ids().to_vec();

    if disp_vec.len() != gene_ids.len() {
        return Err(extendr_api::Error::Other("Length mismatch between dispersions and gene IDs".to_string()));
    }

    let mut robj = disp_vec.into_robj();
    robj.set_names(&gene_ids)?;
    Ok(robj)
}

#[extendr]
/// @export
fn dispersion_function(dds: ExternalPtr<RustDESeqDataSet>) -> extendr_api::Result<Robj> {
    if let Some((asympt_disp, extra_pois)) = dds.dds.dispersion_function() {
        Ok(list!(
            fitType               = "parametric",
            asymptotic_dispersion = asympt_disp,
            extra_poisson         = extra_pois
        ).into_robj())
    } else {
        // Returns R's NULL — matches DESeq2 behavior when not fitted
        Ok(extendr_api::Robj::from(extendr_api::NULL))
    }
}

#[extendr]
/// @export
fn size_factors(dds: ExternalPtr<RustDESeqDataSet>) -> extendr_api::Result<Robj> {
    let sf: Vec<f64> = dds.dds.size_factors()
        .ok_or_else(|| extendr_api::Error::Other(
            "Size factors not found. Run estimateSizeFactors() or deseq() first.".to_string()
        ))?
        .to_vec();  // assuming size_factors() returns &[f64] or ArrayView1 → .to_vec()

    let sample_ids: Vec<String> = dds.dds.counts().sample_ids().to_vec();

    if sf.len() != sample_ids.len() {
        return Err(extendr_api::Error::Other("Length mismatch between size factors and samples".to_string()));
    }

    let mut robj = sf.into_robj();
    robj.set_names(&sample_ids)?;
    Ok(robj)
}

//TODO rowMeans(normalized_counts())
#[extendr]
/// @export
fn base_mean(dds: ExternalPtr<RustDESeqDataSet>) -> extendr_api::Result<Robj> {
    let nc = dds.dds.normalized_counts()
        .ok_or_else(|| extendr_api::Error::Other("Normalized counts not available. Run deseq() first.".to_string()))?;

    let n_samples = dds.dds.n_samples() as f64;
    let base_means: Vec<f64> = nc
        .axis_iter(ndarray::Axis(0))
        .map(|row| row.sum() / n_samples)
        .collect();

    let gene_ids: Vec<String> = dds.dds.counts().gene_ids().to_vec();

    if base_means.len() != gene_ids.len() {
        return Err(extendr_api::Error::Other("Length mismatch between base means and gene IDs".to_string()));
    }

    let mut robj = base_means.into_robj();
    robj.set_names(&gene_ids)?;
    Ok(robj)
}

#[extendr]
/// @export
fn design_matrix(dds: ExternalPtr<RustDESeqDataSet>) -> extendr_api::Result<RMatrix<f64>> {
    // Check if design info is available (ensures deseq() was run)
    match &dds.design_info {
        Some(design_info) => {
            // Get design matrix from DESeqDataSet
            let matrix = dds.dds.design_matrix().unwrap();
            let sample_ids = dds.dds.counts().sample_ids().to_vec();
            let coef_names = design_info.coef_names.clone();
            
            // Use sample IDs as rownames and coef_names from design_info as colnames
            arrayview_to_rmatrix(matrix.view(), sample_ids, coef_names)
        },
        None => {
            Err(extendr_api::Error::Other(
                "DESeq has not been run yet — design matrix and coefficient names are not available. Call deseq() first.".to_string()
            ))
        }
    }
}

#[extendr]
/// @export
fn coefficients(dds: ExternalPtr<RustDESeqDataSet>) -> extendr_api::Result<RMatrix<f64>> {
    let design_info = dds.design_info.as_ref().ok_or_else(|| {
        extendr_api::Error::Other("Please run deseq() first.".to_string())
    })?;

    // get coefficients as Array2
    let coeffs = dds.dds.coefficients().ok_or_else(|| {
        extendr_api::Error::Other("Coefficients not found.".to_string())
    })?;

    // use helper to return a 2D Matrix (Genes x Coefficients)
    let gene_ids = dds.dds.counts().gene_ids().to_vec();
    let coef_names = design_info.coef_names.clone();
    
    arrayview_to_rmatrix(coeffs.view(), gene_ids, coef_names)
}

#[extendr]
/// @export
fn coefficient_se(dds: ExternalPtr<RustDESeqDataSet>) -> extendr_api::Result<RMatrix<f64>> {
    let design_info = dds.design_info.as_ref()
        .ok_or_else(|| extendr_api::Error::Other("Please run deseq() first.".to_string()))?;

    let se_flat = dds.dds.standard_errors()
        .ok_or_else(|| extendr_api::Error::Other("Standard errors not available.".to_string()))?;

    let n_genes = dds.dds.counts().n_genes() as usize;
    let n_coefs = design_info.coef_names.len();

    if se_flat.len() != n_genes * n_coefs {
        return Err(extendr_api::Error::Other(format!(
            "Dimension mismatch: got {} SE values, expected {} × {} = {}",
            se_flat.len(), n_genes, n_coefs, n_genes * n_coefs
        )));
    }

    // Reshape to (n_genes, n_coefs)
    let se_2d = se_flat.to_shape((n_genes, n_coefs))
        .map_err(|e| extendr_api::Error::Other(format!("Cannot reshape: {}", e)))?;

    let gene_ids = dds.dds.counts().gene_ids().to_vec();
    let coef_names = design_info.coef_names.clone();

    arrayview_to_rmatrix(se_2d.view(), gene_ids, coef_names)
}

extendr_module! {
    mod RustDESeq2;
    fn deseq_dataset_from_matrix;
    fn deseq;
    fn results;
    fn counts;
    fn vst_transform;
    fn assay;
    fn dispersions;
    fn dispersion_function;
    fn size_factors;
    fn base_mean;
    fn design_matrix;
    fn coefficients;
    fn coefficient_se;
}
