# DeepMASS2 Documentation
*Version 1.0.0*  
*Released on October 16th, 2023*  
*Maintainer: Ji Hongchao*  
*Email: [jihongchao@caas.cn](mailto:jihongchao@caas.cn)*  

## Introduction

DeepMASS is an innovative software tool providing a robust solution for annotating and 
discovering metabolites within complex biological systems. Its foundation lies in a 
sophisticated deep-learning-based semantic similarity model, seamlessly connecting mass 
spectra to structurally related compounds. This connection effectively maps the chemical 
space of the unknown metabolites. DeepMASS maximizes the utility of mass spectrometry big 
data, positioning itself for further development as data scales continue to expand.

## Installation

**System Requirements:**  
**Operating Systems:**
- Windows 11
- MacOS

**Recommended Hardware:**
- Intel Core i5 or greater
- 16 GB RAM or more
- 5 GB hard drive space

**Installation Steps:**

1. Install [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).
2. Create a new conda environment and activate:

    ```bash
    conda create -n deepmass python=3.8.13
    conda activate deepmass
    ```

3. Clone the repository and enter:

    ```bash
    git clone https://github.com/hcji/DeepMASS2_GUI.git
    cd DeepMASS2_GUI
    ```

4. Install dependencies (for *MacOS*, some dependencies may need manual installation with conda):

    ```bash
    pip install -r requirements.txt
    ```

5. Download the [dependent data](https://github.com/hcji/DeepMASS2_GUI/releases/tag/v0.99.0).

    - Put the following files into the *data* folder:

        ```
        DeepMassStructureDB-v1.0.csv
        references_index_negative_spec2vec.bin
        references_index_positive_spec2vec.bin
        references_spectrums_negative.pickle
        references_spectrums_positive.pickle
        ```

    - Put the following files into the *model* folder:

        ```
        Ms2Vec_allGNPSnegative.hdf5
        Ms2Vec_allGNPSnegative.hdf5.syn1neg.npy
        Ms2Vec_allGNPSnegative.hdf5.wv.vectors.npy
        Ms2Vec_allGNPSpositive.hdf5
        Ms2Vec_allGNPSpositive.hdf5.syn1neg.npy
        Ms2Vec_allGNPSpositive.hdf5.wv.vectors.npy
        ```

    Please note that these dependent data are based on the GNPS dataset. If testing on the CASMI dataset, 
    there may be differences compared to reported results in the paper. For NIST 20 access, 
    refer to the *Training models with NIST 20* section and re-train the model.

6. Run DeepMASS:

    ```bash
    python DeepMASS2.py
    ```


## Quick Start

1. DeepMASS may require some time for auto-loading dependent data. Please wait until the buttons become active.

<div align="center">
<img src="https://github.com/hcji/DeepMASS2_GUI/blob/main/document/imgs/sceenshot_1.png" width="50%">
</div>

2. Click the **Open** button and select an **mgf** file containing one or multiple MS/MS spectra. [Example here](https://github.com/hcji/DeepMASS2_GUI/blob/main/example/all_casmi.mgf)

<div align="center">
<img src="https://github.com/hcji/DeepMASS2_GUI/blob/main/document/imgs/sceenshot_2.png" width="50%">
</div>

3. Click **Run DeepMASS** for annotating with the DeepMASS algorithm or **Run MatchMS** for library matching. Wait for the progress bar to finish.

<div align="center">
<img src="https://github.com/hcji/DeepMASS2_GUI/blob/main/document/imgs/sceenshot_3.png" width="50%">
</div>

4. Click the **Save** button and select the folder path to save the annotation results.

<div align="center">
<img src="https://github.com/hcji/DeepMASS2_GUI/blob/main/document/imgs/sceenshot_4.png" width="50%">
</div>


## Input File Format

DeepMASS accepts **mgf** files containing one or multiple MS/MS spectra. The Mascot Generic Format (MGF) 
is a standard file format used to store mass spectrometry (MS) data, particularly tandem mass spectrometry (MS/MS) data.

- **Header Information:** Optional metadata about the MS data, including instrument parameters, acquisition settings, and other experimental details.
- **BEGIN IONS:** Marks the beginning of the MS/MS data section.
- **Spectrum Information:** Parameters and a list of mass-to-charge ratio (m/z) and intensity pairs.
- **Peak Data:** Describes the mass spectrum with m/z and intensity pairs.
- **END IONS:** Marks the end of the MS/MS data section for a specific spectrum.

Essential DeepMASS input information includes:
- **CHARGE:** Number of charges of the precursor ion.
- **PRECURSOR_MZ:** Precursor ion's mass-to-charge ratio.
- **COMPOUND_NAME:** Identifier for the compound.
- **IONMODE:** Ionization mode of the precursor ion.
- **ADDUCT:** Type of adduct of the precursor ion.

At least one of the following is necessary for DeepMASS inputs:
- **FORMULA:** Molecular formula of the compound.
- **PARENT_MASS:** Monoisotopic mass of the compound.
 

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

[MS-DIAL](http://prime.psc.riken.jp/compms/msdial/main.html) is open-source software for processing and analyzing mass spectrometry (MS) data. You can get 
Follow [tutorial](https://github.com/systemsomicslab/mtbinfo.github.io/blob/master/MS-DIAL/tutorial.md). 
to process DDA/DIA mode metabolomic study MS files and export the alignment result. An example output can be found [here](https://github.com/hcji/DeepMASS2_GUI/blob/main/example/MS_DIAL_example.txt).  

Use the following code to transform the obtained csv or txt file to an mgf file, which can be directly upload in DeepMASS GUI:


        from matchms.exporting import save_as_mgf
        from core.external import load_MS_DIAL_Alginment
        
        spectrums = load_MS_DIAL_Alginment('example/MS_DIAL_example.txt')
        save_as_mgf(spectrums, 'example/MS_DIAL_example.mgf')


## Training models with NIST 20.

  1) Use [LIB2NIST](https://chemdata.nist.gov/mass-spc/ms-search/Library_conversion_tool.html) to export the NIST 20 database with *msp* and *mol* format respectively. 
  Then you will get an *msp* file and a folder of *mol* files. The *msp* file includes the information of peak list of ms/ms and the NIST index of each molecule, while 
  the *mol* file includes the structural information of molecules.
  
  2) Refer the scripts [here] to correlate the ms/ms and the chemical structures. 
  
  
  
  
  2) Refer the scripts [here](https://github.com/hcji/DeepMASS2_Data_Processing/blob/master/Scripts/training_data_collection/clean_nist.py),
       and transform the data into DeepMASS required format.
  3) Refer the scripts [here](https://github.com/hcji/DeepMASS2_Data_Processing/blob/master/Scripts/training_models/train_ms2vec.py),
       and train your *ms2vec* model.
  4) Refer the scripts [here](https://github.com/hcji/DeepMASS2_Data_Processing/blob/master/Scripts/training_models/vectorize_reference_by_ms2vec.py),
       and build index for the spectra of NIST 20.
  5) Copy all the generated files into corresponding folder of DeepMASS.

## Reference

Comming soon ...


