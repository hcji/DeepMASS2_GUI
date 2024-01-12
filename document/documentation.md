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

2. Click **Open** button, select an **mgf** file containing one or multiple MS/MS spectra. 
See [this](https://github.com/hcji/DeepMASS2_GUI/blob/main/example/all_casmi.mgf) for an example. 

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


## Input file format

DeepMASS accept **mgf** file containing one or multiple MS/MS spectra. The Mascot Generic Format (MGF) is a 
standard file format used to store mass spectrometry (MS) data, particularly tandem mass spectrometry (MS/MS) 
data. MGF files are commonly employed for exchanging spectral information between different mass spectrometry 
software tools. Below is a brief description of the typical format of an MGF file:

1. Header Information:

  The file begins with optional header information, which may include metadata about the MS data, such as 
instrument parameters, acquisition settings, and other experimental details.  

2. BEGIN IONS:

  Marks the beginning of the MS/MS data section.  

3. Spectrum Information:

  Each MS/MS spectrum is represented by a set of parameters and a list of mass-to-charge ratio (m/z) and intensity pairs.
Common parameters include precursor ion charge, precursor ion m/z, compound name, challenge information, molecular formula, 
SMILES notation, InChIKey, adduct information, ionization mode, parent mass, and others.  

4. Peak Data:

  The mass spectrum is described by a list of m/z and intensity pairs. Each line typically contains two values separated by a space: the m/z value and the intensity value.
The intensity values indicate the abundance or intensity of ions at specific m/z values.  

5. END IONS:

  Marks the end of the MS/MS data section for a specific spectrum.  

6. Repeating Sections:

  The file may contain multiple MS/MS spectra, each preceded by the "BEGIN IONS" header and followed by the "END IONS" footer.      

**These spectrum information are essential for DeepMASS inputs:**    

  **CHARGE:** Indicates the number of charge of the precursor ion.   
  **PRECURSOR_MZ:** Specifies the precursor ion's mass-to-charge ratio.   
  **COMPOUND_NAME:** Not the true name of the compound, but a identifier (in fact the true name is unknown), for example, "unknown_1".    
  **IONMODE:** Specifies the ionization mode of the precursor ion.    
  **ADDUCT:** Specifies the type of the adduct of the precursor ion. 

**At least one of the following information is necessary for DeepMASS inputs:**  

  **FORMULA:** Represents the molecular formula of the compound.    
  **PARENT_MASS**: Specifies the monoisotopic mass of the compound.
 
Here's a simplified example of an MGF file:

        BEGIN IONS
        CHARGE=1-
        PRECURSOR_MZ=194.04588176800002
        COMPOUND_NAME=challenge_255
        CHALLENGE=CASMI2016_Challenge-006
        FORMULA=C9H9NO4
        SMILES=CC(=O)NC1=CC=C(O)C(=C1)C(O)=O
        INCHIKEY=GEFDRROBUCULOD-UHFFFAOYSA-N
        ADDUCT=[M-H]-
        IONMODE=negative
        PARENT_MASS=195.053157768
        132.0454 0.006008828855732822 
        134.0245 0.0022295306804225633 
        149.0483 0.00297960480724686 
        150.056 1.0 
        194.0459 0.36915136378467533 
        END IONS
        
        BEGIN IONS
        CHARGE=1+
        PRECURSOR_MZ=223.970551912
        COMPOUND_NAME=challenge_368
        CHALLENGE=CASMI2016_Challenge-119
        FORMULA=C9H6BrNO
        SMILES=BrC1=CC2=C(NC(=O)C=C2)C=C1
        INCHIKEY=YLAFBGATSQRSTB-UHFFFAOYSA-N
        ADDUCT=[M+H]+
        IONMODE=positive
        PARENT_MASS=222.963275912
        127.0418 0.01388754827761238 
        145.0524 0.0058134492435322465 
        205.9602 0.1043324271433068 
        223.9708 1.0 
        END IONS


## Constructing *mgf* file from MS-DIAL

  1) Process your ms files of DDA/DIA mode metabolomic study following the MS-DIAL [tutorial](https://mtbinfo-team.github.io/mtbinfo.github.io/MS-DIAL/tutorial).
  2) Export the alignment result with *txt* format. Refer the [tutorial-section 5-6-(B)](https://mtbinfo-team.github.io/mtbinfo.github.io/MS-DIAL/tutorial#section-5-6).
  3) Refer the scripts [here](https://github.com/hcji/DeepMASS2_Data_Processing/blob/master/Scripts/test_data_collection/processing_tomato.py).

## Training models with NIST 20.

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


