# DeepMASS2 Documentation
### Oct 16th, 2023

**Package:**
	Version: 1.0.0, released at Oct 16th, 2023    
	Maintainer: Ji Hongchao *(jihongchao@caas.cn)*     

## Introduction

DeepMASS is an innovative software tool offering a powerful solution for annotating and 
discovering metabolites within complex biological systems. Its foundation lies in a sophisticated 
deep-learning-based semantic similarity model, which seamlessly connects mass spectra to 
structurally related compounds, effectively mapping the chemical space of the unknown. 
DeepMASS maximizes the utility of mass spectrometry big data, positioning itself for 
further development as data scales continue to expand.


## Installation

**System Recommended:**    
Operating Systems: 

    - Windows 11
    - MacOS
 
Recommended Hardware:

    - Intel Core i5 or greater
    - 16 GB RAM or more
    - 5 GB hard drive space

**Please follow the following installation steps:**

1. Install [Anaconda](https://www.anaconda.com/)  or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)   
2. Create a new conda environment and activate:

        conda create -n deepmass python=3.8.13
        conda activate deepmass

3. Clone the repository and enter:

        git clone https://github.com/hcji/DeepMASS2_GUI.git
        cd DeepMASS2_GUI
    
4. Install dependency (for *MacOS*, some dependency may install with conda manually):

        pip install -r requirements.txt
        
5. Download the [dependent data](https://github.com/hcji/DeepMASS2_GUI/releases/tag/v0.99.0).     

    1) put the following files into *data* folder:
    
                DeepMassStructureDB-v1.0.csv
                references_index_negative_spec2vec.bin
                references_index_positive_spec2vec.bin
                references_spectrums_negative.pickle
                references_spectrums_positive.pickle
    
    2) put the following files into *model* folder:
    
                Ms2Vec_allGNPSnegative.hdf5
                Ms2Vec_allGNPSnegative.hdf5.syn1neg.npy
                Ms2Vec_allGNPSnegative.hdf5.wv.vectors.npy
                Ms2Vec_allGNPSpositive.hdf5
                Ms2Vec_allGNPSpositive.hdf5.syn1neg.npy
                Ms2Vec_allGNPSpositive.hdf5.wv.vectors.npy

Please note, these dependent data are introduced as public version of the published paper. It 
means they are based on the GNPS dataset only. If you test on CASMI dataset with these version, 
there may be some difference compared with the reported results in the paper. If you have the 
accessibility of NIST 20, please referred the introduction of the *Advanced usage* part, and 
re-train the model.

6. Run DeepMASS

        python DeepMASS2.py


## Quick start

1. DeepMASS may need some time for auto-loading the dependent data. Please wait until the 
buttons become active.

<div align="center">
<img src="https://github.com/hcji/DeepMASS2_GUI/blob/main/document/imgs/sceenshot_1.png" width="50%">
</div>

2. Click **Open** button, select a mgf file containing one or multiple MS/MS spectra. 
See [this](https://github.com/hcji/DeepMASS2_GUI/blob/main/example/all_casmi.mgf) for an example. 
Except *SMILES* and *INCHIKEY* lines, the other meta information is necessary for each spectrum.

<div align="center">
<img src="https://github.com/hcji/DeepMASS2_GUI/blob/main/document/imgs/sceenshot_2.png" width="50%">
</div>

3. Click **Run DeepMASS** for annotating with DeepMASS algorithm, or click **Run MatchMS** for library matching. 
Wait for the progress bar to finish.

<div align="center">
<img src="https://github.com/hcji/DeepMASS2_GUI/blob/main/document/imgs/sceenshot_3.png" width="50%">
</div>

4. Click **Save** button, select the folder path to save the annotation results.

<div align="center">
<img src="https://github.com/hcji/DeepMASS2_GUI/blob/main/document/imgs/sceenshot_4.png" width="50%">
</div>


## Advanced usage

1. Constructing *mgf* file from MS-DIAL results.

    1) Process your ms files of DDA/DIA mode metabolomic study following the MS-DIAL [tutorial](https://mtbinfo-team.github.io/mtbinfo.github.io/MS-DIAL/tutorial).
    2) Export the alignment result with *txt* format. Refer the [tutorial-section 5-6-(B)](https://mtbinfo-team.github.io/mtbinfo.github.io/MS-DIAL/tutorial#section-5-6).
    3) Refer the scripts [here](https://github.com/hcji/DeepMASS2_Data_Processing/blob/master/Scripts/test_data_collection/processing_tomato.py).

2. Training models with NIST 20 spectra.

    1) Use [LIB2NIST](https://chemdata.nist.gov/mass-spc/ms-search/Library_conversion_tool.html) tool to export NIST 20 database to *mgf* format.
    2) Refer the scripts [here](https://github.com/hcji/DeepMASS2_Data_Processing/blob/master/Scripts/training_data_collection/clean_nist.py),
       and transform the data into DeepMASS required format.
    3) Refer the scripts [here](https://github.com/hcji/DeepMASS2_Data_Processing/blob/master/Scripts/training_models/train_ms2vec.py),
       and train your *ms2vec* model.
    4) Refer the scripts [here](https://github.com/hcji/DeepMASS2_Data_Processing/blob/master/Scripts/training_models/vectorize_reference_by_ms2vec.py),
       and build index for the spectra of NIST 20.
    5) Copy all the generated files into corresponding folder of DeepMASS.

## Reference

Comming soon ...


