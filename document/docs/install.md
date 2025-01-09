## System Requirements

**Operating Systems:**

- Windows 11
- MacOS

**Recommended Hardware:**

- Intel Core i5 or greater
- 16 GB RAM or more
- 5 GB hard drive space


## Installation Steps

1. Install [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).
2. Create a new conda environment and activate:

    ```
    conda create -n deepmass python=3.8.13
    conda activate deepmass
    ```

3. Clone the repository and enter:

    ```
    git clone https://github.com/hcji/DeepMASS2_GUI.git
    cd DeepMASS2_GUI
    ```

4. Install dependencies (for *MacOS*, some dependencies may need manual installation with conda):

    ```
    pip install -r requirements.txt
    ```

5. Download the [dependent data](https://github.com/hcji/DeepMASS2_GUI/releases/tag/v0.99.2).

    - Put the following files into the *data* folder:

        ```
        DeepMassStructureDB-v1.1.csv
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

    Please note that these dependent data are based on the GNPS dataset only. If you have licence of NIST software,
    please refer the next section: **Training models with NIST data** 

6. Run DeepMASS:

    ```
    python DeepMASS2.py
    ```

## Training models with NIST data
The [NIST MS/MS Library](https://www.nist.gov/programs-projects/nist20-updates-nist-tandem-and-electron-ionization-spectral-libraries) is not available for download directly, 
but can be purchased from an authorized distributor and exported using the instructions below.


### Exporting NIST data
*Note: this step requires a Windows System or Virtual Machine.*

The spectra and associated compounds can be exported to MSP/MOL format using the free [lib2nist software](https://chemdata.nist.gov/mass-spc/ms-search/Library_conversion_tool.html). 
The resulting export will contain a single MSP file with all of the mass spectra, 
and multiple MOL files which include the molecular structure information (linked to the spectra by ID). 
The screenshot below indicates appropriate lib2nist export settings.

[![nist-export.png](https://i.postimg.cc/c40KbnpK/nist-export.png)](https://postimg.cc/QVYxBHxs)

After exporting the files, create a directory on your primary computer and save them there. 
If done correctly, inside "nist_20" there should be a single .MSP file with all the spectra, 
`hr_nist_msms.MSP`, and a directory of .MOL files, `hr_nist_msms.MOL`.


### Preprocessing NIST data
Refer the [Python script](https://github.com/hcji/DeepMASS2_GUI/blob/main/scripts/Processing_NIST20_data.ipynb) 
processes mass spectrometry (MS/MS) data and corresponding chemical structures 
to extract and clean metadata, filter the data based on adduct types, and organize the data for further analysis. 


### Training the Word2Vec model
Refer the [Python script](https://github.com/hcji/DeepMASS2_GUI/blob/main/scripts/Traning_Spec2vec_Model.ipynb) 
processes mass spectrometry (MS/MS) data to train word2vec models for positive and negative ion mode 
spectra using the spec2vec approach.This script ultimately creates two spec2vec models: one for positive 
ion mode spectra and one for negative ion mode spectra. 


### Generating spectral representation index
Refer the [Python script](https://github.com/hcji/DeepMASS2_GUI/blob/main/scripts/Vectorize_Spectra.ipynb) 
processes mass spectrometry (MS/MS) data to compute vector representations of spectra using a pre-trained 
word2vec model, and then indexes these vectors for fast similarity searching using the HNSW 
algorithm. It handles both positive and negative ion mode spectra. This script ultimately creates and saves 
HNSW indices for fast similarity searching of positive and negative ion mode spectra using the spec2vec vector 
representations.

### Replacing the data
Copy all the generated files into corresponding folder of DeepMASS, refer the installation step 5.
