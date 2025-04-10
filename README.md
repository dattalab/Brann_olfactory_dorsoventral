# Brann_olfactory_dorsoventral

Code to replicate analyses by Brann et al.

# Installation

## Requirements 

1. Create a new conda environment. `conda create -n dv_score python=3.9`
2. Activate that env `conda activate dv_score`.
3. Clone and enter this repo: `git clone git@github.com:dattalab/Brann_olfactory_dorsoventral.git && cd Brann_olfactory_dorsoventral`
4. Install the code in this directory via `pip install -e .`
5. To install the specific versions of packages used when this repo was created do `pip install -r requirements.txt`. The additional requirements for running the notebooks in this repo are to `pip install numpy seaborn scikit-learn jupyter notebook`. The scripts also rely on `pip install pysam scanpy`.


## Data
1. Processed data is available on the NCBI GEO at Accession number [GSE173947](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173947) and raw fastq files can be found on the SRA (accession SRP318630).
2. Data were preprocessed by running the Nextflow [pipeline](scripts/nextflow_pipeline.nf) in the scripts folder.
3. Additional instructions for how to work with the raw data and to work with the olfactory gene expression programs (GEPs) can be found in the follow GitHub repo: [Tsukahara_Brann_OSN](https://github.com/dattalab/Tsukahara_Brann_OSN).

# Examples
Code to generate key results, focusing on the dorsoventral (DV) score.

1. Open a new jupyter notebook with `jupyter notebook`.
2. Run the [notebooks](./notebooks). 
3. Additional stand-alone [scripts](./scripts) demonstrate the [Nextflow](https://nextflow.io/) pipelines that were used to uniformly preprocess scRNA-seq data, as well the [scVI](https://scvi-tools.org/) and scANVI models that were used for data integration.

# Contact
For more details, please post an issue here or contact the authors.