#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: spatial_10x_visium/mouse_brain_coronal_section1
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1/CytAssist_FFPE_Mouse_Brain_Rep1_spatial.tar.gz"
    dataset_name: Mouse Brain Coronal Section 1 (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/mouse-brain-coronal-section-1-ffpe-2-standard"
    dataset_summary: Gene expression library of Mouse Brain (CytAssist FFPE) using the Mouse Whole Transcriptome Probe Set
    dataset_description: "FFPE Mouse Brain tissue blocks sectioned as described in Visium CytAssist Spatial Gene Expression for FFPE - Tissue Preparation Guide Demonstrated Protocol. The H&E stained glass slide with tissue section was processed via Visium CytAssist instrument to transfer analytes to a Visium CytAssist Spatial Gene Expression slide. The probe extension and library construction steps follow the standard Visium for FFPE workflow outside of the instrument. The H&E image was acquired using Olympus VS200 Slide Scanning Microscope. Sequencing depth was 53,497 reads per spot. Sequencing configuration: 28bp read 1 (16bp Visium spatial barcode, 12bp UMI), 90bp read 2 (transcript), 10bp i7 sample barcode and 10bp i5 sample barcode. Key metrics include: 2,310 spots detected under tissue; 6,736 median genes per spot; 24,862 median UMI counts per spot."
    dataset_organism: Mus musculus

  - id: spatial_10x_visium/human_colorectal_cancer
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Colorectal_Cancer/CytAssist_11mm_FFPE_Human_Colorectal_Cancer_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Colorectal_Cancer/CytAssist_11mm_FFPE_Human_Colorectal_Cancer_spatial.tar.gz"
    dataset_name: Human Colorectal Cancer (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/human-colorectal-cancer-11-mm-capture-area-ffpe-2-standard"
    dataset_summary: Gene expression library of Human Colorectal Cancer (CytAssist FFPE) using the Human Whole Transcriptome Probe Set
    dataset_description: "The tissue was sectioned as described in the Visium CytAssist Spatial Gene Expression for FFPE Tissue Preparation Guide (CG000518). Tissue section of 5 µm was placed on a standard glass slide, then stained following the Deparaffinization, H&E Staining, Imaging & Decrosslinking Demonstrated Protocol (CG000520). The glass slide with tissue section was processed via Visium CytAssist instrument to transfer analytes to a Visium CytAssist Spatial Gene Expression Slide v2, with 11 mm capture areas following the Visium CytAssist Spatial Gene Expression Reagent Kits User Guide (CG000495)."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_heart
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Human_Heart/V1_Human_Heart_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Human_Heart/V1_Human_Heart_spatial.tar.gz"
    dataset_name: Human Heart
    dataset_url: "https://www.10xgenomics.com/datasets/human-heart-1-standard-1-0-0"
    dataset_summary: V1_Human_Heart
    dataset_description: "10x Genomics obtained fresh frozen human heart tissue from BioIVT Asterand. The tissue was embedded and cryosectioned as described in Visium Spatial Protocols - Tissue Preparation Guide Demonstrated Protocol (CG000240). Tissue sections of 10 µm thickness were placed on Visium Gene Expression Slides."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/mouse_embryo
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_11mm_FFPE_Mouse_Embryo/CytAssist_11mm_FFPE_Mouse_Embryo_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_11mm_FFPE_Mouse_Embryo/CytAssist_11mm_FFPE_Mouse_Embryo_spatial.tar.gz"
    dataset_name: Mouse Embryo (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/visium-cytassist-mouse-embryo-11-mm-capture-area-ffpe-2-standard"
    dataset_summary: Gene expression library of Mouse Embryo (CytAssist FFPE) using the Mouse Whole Transcriptome Probe Set
    dataset_description: "The tissue was sectioned as described in Visium CytAssist Spatial Gene Expression for FFPE Tissue Preparation Guide Demonstrated Protocol CG000518. Tissue sections of 5 µm was placed on a standard glass slide, and H&E-stained following deparaffinization. Sections were coverslipped with 85% glycerol, imaged, decoverslipped, followed by dehydration & decrosslinking (Demonstrated Protocol CG000520). The glass slide with the tissue section was processed with the Visium CytAssist instrument to transfer analytes to a Visium CytAssist Spatial Gene Expression slide (11 mm Capture Area). The probe extension and library construction steps follow the standard Visium for FFPE workflow outside of the instrument."
    dataset_organism: Mus musculus

  - id: spatial_10x_visium/mouse_olfactory_bulb 
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz"
    dataset_name: Adult Mouse Olfactory Bulb (Fresh Frozen)
    dataset_url: "https://www.10xgenomics.com/datasets/adult-mouse-olfactory-bulb-1-standard-1"
    dataset_summary: 10X Genomics obtained fresh frozen mouse olfactory bulb tissue from BioIVT.
    dataset_description: "The tissue was embedded and cryosectioned as described in Visium Spatial Protocols Tissue Preparation Guide (Demonstrated Protocol CG000240). Tissue sections of 10µm were placed on Visium Gene Expression slides, then fixed and stained following Methanol Fixation, H&E Staining & Imaging for Visium Spatial Protocols (CG000160)."
    dataset_organism: Mus musculus

  - id: spatial_10x_visium/human_breast_cancer_1
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer/Parent_Visium_Human_BreastCancer_spatial.tar.gz"
    dataset_name: Adult Human Breast Cancer
    dataset_url: "https://www.10xgenomics.com/datasets/human-breast-cancer-whole-transcriptome-analysis-1-standard-1-2-0"
    dataset_summary: Whole transcriptome analysis, Adult Human Breast Cancer (Visium)
    dataset_description: "10X Genomics obtained fresh frozen human Invasive Lobular Carcinoma breast tissue from BioIVT Asterand. The tissue was embedded and cryosectioned as described in Visium Spatial Protocols Tissue Preparation Guide Demonstrated Protocol (CG000240). Tissue sections of 10µm were placed on Visium Gene Expression slides and fixed and stained following Methanol Fixation, H&E Staining & Imaging for Visium Spatial Protocols (CG000160)."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_lymph_node 
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_spatial.tar.gz"
    dataset_name: Human Lymph Node (Fresh Frozen)
    dataset_url: "https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-0-0"
    dataset_summary: Whole transcriptome analysis, Human Lymph Node
    dataset_description: "10x Genomics obtained fresh frozen human lymph node from BioIVT Asterand. The tissue was embedded and cryosectioned as described in Visium Spatial Protocols - Tissue Preparation Guide Demonstrated Protocol (CG000240). Tissue sections of 10 µm thickness were placed on Visium Gene Expression Slides."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_normal_prostate 
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Normal_Prostate/Visium_FFPE_Human_Normal_Prostate_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Normal_Prostate/Visium_FFPE_Human_Normal_Prostate_spatial.tar.gz"
    dataset_name: Human Normal Prostate (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/normal-human-prostate-ffpe-1-standard-1-3-0"
    dataset_summary: Gene expression library of Human Normal Prostate (Visium FFPE) using the Human Whole Transcriptome Probe Set
    dataset_description: "10x Genomics obtained FFPE human prostate tissue from Indivumed Human Tissue Specimens. The tissue was sectioned as described in Visium Spatial Gene Expression for FFPE – Tissue Preparation Guide Demonstrated Protocol (CG000408). Tissue sections of 5 µm were placed on Visium Gene Expression slides, then stained following Deparaffinization, H&E Staining, Imaging & Decrosslinking Demonstrated Protocol (CG000409)."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_prostate_cancer
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Prostate_IF/Visium_FFPE_Human_Prostate_IF_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Prostate_IF/Visium_FFPE_Human_Prostate_IF_spatial.tar.gz"
    dataset_name: Human Prostate Cancer (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/human-prostate-cancer-adjacent-normal-section-with-if-staining-ffpe-1-standard"
    dataset_summary: Gene expression library of Human Prostate Cancer (Visium FFPE) with an IF image using the Human Whole Transcriptome Probe Set
    dataset_description: "10x Genomics obtained FFPE human prostate tissue from Indivumed Human Tissue Specimens. Original diagnosis with adenocarcinoma. The tissue was sectioned as described in Visium Spatial Gene Expression for FFPE Tissue Preparation Guide Demonstrated Protocol (CG000408). Tissue sections of 10 µm were placed on Visium Gene Expression slides, then stained following Deparaffinization, Decrosslinking, Immunofluorescence Staining & Imaging Demonstrated Protocol (CG000410)."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_cerebellum
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Cerebellum/Parent_Visium_Human_Cerebellum_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_Cerebellum/Parent_Visium_Human_Cerebellum_spatial.tar.gz"
    dataset_name: Adult Human Cerebellum (Visium)
    dataset_url: "https://www.10xgenomics.com/datasets/human-cerebellum-whole-transcriptome-analysis-1-standard-1-2-0"
    dataset_summary: Human Cerebellum: Whole Transcriptome Analysis
    dataset_description: "10X Genomics obtained fresh frozen human cerebellum tissue from BioIVT Asterand. The tissue was embedded and cryosectioned as described in Visium Spatial Protocols Tissue Preparation Guide (Demonstrated Protocol CG000240). Tissue sections of 10µm were placed on Visium Gene Expression slides and fixed and stained following Methanol Fixation, H&E Staining & Imaging for Visium Spatial Protocols (CG000160)."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/mouse_kidney_v1
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Kidney/V1_Mouse_Kidney_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Kidney/V1_Mouse_Kidney_spatial.tar.gz"
    dataset_name: Mouse Kidney Section (Coronal)
    dataset_url: "https://www.10xgenomics.com/datasets/mouse-kidney-section-coronal-1-standard-1-1-0"
    dataset_summary: Mouse Kidney: Whole Transcriptome Analysis
    dataset_description: "10x Genomics obtained fresh frozen mouse kidney tissue from BioIVT Asterand. The tissue was embedded and cryosectioned as described in Visium Spatial Protocols - Tissue Preparation Guide Demonstrated Protocol (CG000240). Tissue sections of 10 µm thickness from a slice of the coronal plane were placed on Visium Gene Expression slides, then stained following the Methanol Fixation, H&E Staining & Imaging Demonstrated Protocol (CG000160)."
    dataset_organism: Mus musculus

  - id: spatial_10x_visium/human_lung_cancer
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Lung_Cancer/CytAssist_11mm_FFPE_Human_Lung_Cancer_spatial.tar.gz"
    dataset_name: Human Lung Cancer (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/human-lung-cancer-11-mm-capture-area-ffpe-2-standard"
    dataset_summary: Gene expression library of Human Lung Cancer (CytAssist FFPE) using the Human Whole Transcriptome Probe Set
    dataset_description: "10x Genomics obtained FFPE human lung cancer tissue from Avaden Biosciences. The tissue was sectioned as described in the Visium CytAssist Spatial Gene Expression for FFPE Tissue Preparation Guide (CG000518). Tissue section of 5 µm was placed on a standard glass slide, then stained following the Deparaffinization, H&E Staining, Imaging & Decrosslinking Demonstrated Protocol (CG000520). The glass slide with tissue section was processed via Visium CytAssist instrument to transfer analytes to a Visium CytAssist Spatial Gene Expression Slide v2, with 11 mm capture areas following the Visium CytAssist Spatial Gene Expression Reagent Kits User Guide (CG000495)."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_glioblastoma
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_spatial.tar.gz"
    dataset_name: Human Glioblastoma (FFPE)
    dataset_url: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_web_summary.html"
    dataset_summary: Gene expression library of Human Glioblastoma (CytAssist FFPE) using the Human Whole Transcriptome Probe Set
    dataset_description: ""
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_kidney
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Kidney/CytAssist_11mm_FFPE_Human_Kidney_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Kidney/CytAssist_11mm_FFPE_Human_Kidney_spatial.tar.gz"
    dataset_name: Human Kidney (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/human-kidney-11-mm-capture-area-ffpe-2-standard"
    dataset_summary: Gene expression library of Human Kidney (CytAssist FFPE) using the Human Whole Transcriptome Probe Set
    dataset_description: "10x Genomics obtained FFPE human kidney tissue from Avaden Biosciences. The tissue was sectioned as described in the Visium CytAssist Spatial Gene Expression for FFPE – Tissue Preparation Guide (CG000518). Tissue section of 5 µm was placed on a standard glass slide, then stained following the Deparaffinization, H&E Staining, Imaging & Decrosslinking Demonstrated Protocol (CG000520). The glass slide with tissue section was processed via Visium CytAssist instrument to transfer analytes to a Visium CytAssist Spatial Gene Expression Slide v2, with 11 mm capture areas following the Visium CytAssist Spatial Gene Expression Reagent Kits User Guide (CG000495)."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_intestinal_cancer
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Intestinal_Cancer/Visium_FFPE_Human_Intestinal_Cancer_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Intestinal_Cancer/Visium_FFPE_Human_Intestinal_Cancer_spatial.tar.gz"
    dataset_name: Human Intestine Cancer (FPPE)
    dataset_url: "https://www.10xgenomics.com/datasets/human-intestine-cancer-1-standard"
    dataset_summary: Gene expression library of Human Intestinal Cancer (Visium FFPE) using the Human Whole Transcriptome Probe Set
    dataset_description: "5 µm section from Human Intestinal Cancer. FFPE tissue purchased from BioIVT Asterand Human Tissue Specimens. Libraries were prepared following the Visium Spatial Gene Expression Reagent Kits for FFPE User Guide (CG000407 Rev A)."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_skin_melanoma
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz"
    dataset_name: Human Melanoma, IF Stained (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/human-melanoma-if-stained-ffpe-2-standard"
    dataset_summary: Gene expression library of Human Skin Melanoma (CytAssist FFPE) using the Human Whole Transcriptome Probe Set
    dataset_description: "10x Genomics obtained FFPE Human Melanoma tissue blocks from Avaden Biosciences. The tissue was sectioned as described in Visium CytAssist Spatial Gene Expression for FFPE Tissue Preparation Guide Demonstrated Protocol (CG000518). Tissue sections of 5 µm was placed on a standard glass slide, deparaffinized followed by immunofluorescence (IF) staining. Sections were coverslipped with 85% glycerol, imaged, decoverslipped, followed by dehydration & decrosslinking Demonstrated Protocol (CG000519). The glass slide with tissue section was processed via Visium CytAssist instrument to transfer analytes to a Visium CytAssist Spatial Gene Expression slide. The probe extension and library construction steps follow the standard Visium for FFPE workflow outside of the instrument."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_cervical_cancer
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Cervical_Cancer/Visium_FFPE_Human_Cervical_Cancer_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Cervical_Cancer/Visium_FFPE_Human_Cervical_Cancer_spatial.tar.gz"
    dataset_name: Human Cervical Cancer (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/human-cervical-cancer-1-standard"
    dataset_summary: Gene expression library of Human Cervical Cancer (Visium FFPE) using the Human Whole Transcriptome Probe Set
    dataset_description: "5 µm section from squamous cell carcinoma of human cervical cancer. FFPE tissue purchased from Discovery Life Sciences."
    dataset_organism: Homo sapiens

  - id: spatial_10x_visium/human_breast_cancer_2
    input_expression: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
    input_spatial: "https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Breast_Cancer/Visium_FFPE_Human_Breast_Cancer_spatial.tar.gz"
    dataset_name: Human Breast Cancer: Ductal Carcinoma In Situ, Invasive Carcinoma (FFPE)
    dataset_url: "https://www.10xgenomics.com/datasets/human-breast-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0"
    dataset_summary: Gene expression library of Human Breast Cancer (Visium FFPE) using the Human Whole Transcriptome Probe Set
    dataset_description: "10x Genomics obtained FFPE human breast tissue from BioIVT Asterand Human Tissue Specimens. The tissue was annotated with Ductal Carcinoma In Situ, Invasive Carcinoma. The tissue was sectioned as described in Visium Spatial Gene Expression for FFPE – Tissue Preparation Guide Demonstrated Protocol (CG000408). Tissue sections of 5 µm were placed on Visium Gene Expression slides, then stained following Deparaffinization, H&E Staining, Imaging & Decrosslinking Demonstrated Protocol (CG000409)."
    dataset_organism: Homo sapiens

normalization_methods: [log_cp10k]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
publish_dir: s3://openproblems-data/resources/datasets
spot_filter_min_genes: 200
gene_filter_min_spots: 50
remove_mitochondrial: true
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
  withLabel: highmem {
    memory = '350GB'
  }
  withName: '.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision integration_build \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_10x_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config 