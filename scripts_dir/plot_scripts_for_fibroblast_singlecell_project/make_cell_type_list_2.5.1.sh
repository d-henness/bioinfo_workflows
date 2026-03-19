for cell in Fibroblasts:breast Fibroblasts:foreskin Tissue_stem_cells:BM_MSC MSC Tissue_stem_cells:iliac_MSC Chondrocytes:MSC-derived Tissue_stem_cells:BM_MSC:BMP2 Tissue_stem_cells:BM_MSC:TGFb3 Osteoblasts Osteoblasts:BMP2 Smooth_muscle_cells:bronchial Smooth_muscle_cells:bronchial:vit_D Smooth_muscle_cells:vascular:IL-17 iPS_cells:adipose_stem_cells iPS_cells:CRL2097_foreskin iPS_cells:fibroblasts iPS_cells:PDB_fibroblasts iPS_cells:skin_fibroblast Tissue_stem_cells:adipose-derived_MSC_AM3;do
    echo "${cell}_normal_skin"
    echo "${cell}_scar"
    echo "${cell}_keloid"
    echo "${cell}_ssc"
done
