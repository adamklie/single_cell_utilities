# Function to transfer labels from reference data
label_transfer = function(igvf_object, ref_object, ref_object_label, predictions_label){
    DefaultAssay(igvf_object) <- "SCT"
    DefaultAssay(ref_object) <- "SCT"

    transfer_anchors <- FindTransferAnchors(
        reference = ref_object,
        query = igvf_object,
        reference.assay = "SCT",
        features = rownames(igvf_object)[rownames(igvf_object) %in% rownames(ref_object)],
        normalization.method = "SCT",
        reference.reduction = "pca",
        recompute.residuals = FALSE,
        dims = 1:50,
        verbose = F)
    
    # gen_celltype
    predictions <- TransferData(
        anchorset = transfer_anchors, 
        refdata = ref_object@meta.data[,colnames(ref_object@meta.data) %in% ref_object_label], 
        weight.reduction = igvf_object[['pca']],
        dims = 1:50,
        verbose = F)

    igvf_object <- AddMetaData(
        object = igvf_object,
        metadata = predictions)
    
    igvf_object@meta.data[[predictions_label]] = igvf_object$predicted.id
         
    out = igvf_object

    }
