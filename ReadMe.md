# DeepMASS2

DeepMASS2 is a cross-platform GUI software tool, which enables deep-learning based metabolite annotation
 via semantic similarity analysis of mass spectral language. This approach enables the prediction 
 of structurally related metabolites for the unknown compounds. By considering the chemical space, these 
 structurally related metabolites provide valuable information about the potential location of the unknown 
 metabolites and assist in ranking candidates obtained from molecular structure databases. 

https://github.com/hcji/DeepMASS2_GUI/assets/17610691/16589373-3c8f-4e85-9e7b-41310dfb34b2

## News
- [10/2024] Using DeepMASS2, we made a web UI interface, check it out! [Website](http://218.245.102.112/)
- [04/2025] We've enhanced DeepMASS with a distributed search feature. Curious how it works? See the guide: [Distributed_ReadMe.md](./Distributed_ReadMe.md)

## Data Input Specifications

To ensure **DeepMASS2** accurately identifies metabolites and correctly names output files, your input data must include specific metadata tags. While various formats are supported, the following specifications use the **.mgf** format as a reference.

For a practical example of a compatible `.mgf` file structure, please refer to the [CASMI Example Dataset](https://github.com/hcji/DeepMASS2_Data_Processing/blob/master/Example/CASMI/all_casmi.mgf).

**Mandatory Metadata Tags**

1. **Precursor m/z** - **Required**  
   This tag specifies the precursor ion mass-to-charge ratio ($m/z$). This is a fundamental requirement for the search engine to filter candidates within the structural databases.  
   - **Example**: `PRECURSOR_MZ=517.22098`

2. **Ion Mode** - **Required**  
   This specifies the polarity of the data, ensuring DeepMASS2 utilizes the correct deep-learning model (Positive vs. Negative) and reference libraries.  
   - **Positive mode**: `IONMODE=positive`
   - **Negative mode**: `IONMODE=negative`

**Optional & Recommended Tags**

3. **Automatic Naming** - **Recommended**  
   DeepMASS2 uses this tag to define the output filename for the exported semantic similarity analysis. If provided, the results will be saved as `<COMPOUND_NAME>.csv`.  
   - **Example**: `COMPOUND_NAME=challenge_0` will export results to `challenge_0.csv`.

4. **Molecular Formula** - **Optional**  
   Adding the molecular formula helps the semantic similarity engine constrain potential chemical space, significantly improving the ranking accuracy of structurally related metabolites.  
   - **Example**: `FORMULA=C25H38O9`
   
## Installation
Please follow the following installation steps:

1. Install [Anaconda](https://www.anaconda.com/)  or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)   

2. Download the [lastest released version](https://github.com/hcji/DeepMASS2_GUI/tags) of source codes, unzip, enter:

        # download and unzip, then,
        cd DeepMASS2_GUI

3. For the installation of dependencies
   
   Use the following step for installation:

        conda env create -f environment.yml
        conda activate deepmass
        
   Or follow the steps below:

        (1) Create a new conda environment and activate:
                conda create -n deepmass python=3.8.13
                conda activate deepmass

        (2) Install dependency (note, for *MacOS* some dependency may install with conda manually):

                pip install -r requirements.txt
        
4. Download the [dependent data](https://github.com/hcji/DeepMASS2_GUI/releases/tag/v0.99.1).    
    1) put the following files into *data* folder:
    
                DeepMassStructureDB-v1.1.csv
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

5. Run DeepMASS

        python DeepMASS2.py


## Release

* [Version 0.99.0](https://github.com/hcji/DeepMASS2_GUI/releases/tag/v0.99.0)
* [Version 0.99.1](https://github.com/hcji/DeepMASS2_GUI/releases/tag/v0.99.1)

## Documentation

For the details on how to use DeepMASS, please check [Ducomentation](https://hcji.github.io/DeepMASS2_GUI/).

## Citation

In preparation
        
## Contact

Ji Hongchao   
E-mail: ji.hongchao@foxmail.com    
<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0002-7364-0741" href="https://orcid.org/0000-0002-7364-0741" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">https://orcid.org/0000-0002-7364-0741</a></div>
    
WeChat public account: Chemocoder    
<img align="center" src="https://github.com/hcji/hcji/blob/main/img/qrcode.jpg" width="20%"/>
