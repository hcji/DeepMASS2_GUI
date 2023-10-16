# DeepMASS2

DeepMASS2 is a cross-platform GUI software tool, which enables deep-learning based metabolite annotation
 via semantic similarity analysis of mass spectral language. This approach enables the prediction 
 of structurally related metabolites for the unknown compounds. By considering the chemical space, these 
 structurally related metabolites provide valuable information about the potential location of the unknown 
 metabolites and assist in ranking candidates obtained from molecular structure databases. 


## Installation
Please follow the following installation steps:

1. Install [Anaconda](https://www.anaconda.com/)  or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)   
2. Create a new conda environment and activate:

        conda create -n deepmass python=3.8.13
        conda activate deepmass

3. Clone the repository and enter:

        git clone https://github.com/hcji/DeepMASS2_GUI.git
        cd DeepMASS2_GUI
    
4. Install dependency (note, for *MacOS* some dependency may install with conda manually):

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

6. Run DeepMASS

        python DeepMASS2.py


## Release

* [Version 0.99.0](https://github.com/hcji/DeepMASS2_GUI/releases/tag/v0.99.0)

## Documentation

For the details on how to use DeepMASS, please check [Ducomentation](https://github.com/hcji/DeepMASS2_GUI/blob/main/document/documentation.md).

## Citation

In preparation
        
## Contact

Ji Hongchao   
E-mail: ji.hongchao@foxmail.com    
<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0002-7364-0741" href="https://orcid.org/0000-0002-7364-0741" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">https://orcid.org/0000-0002-7364-0741</a></div>
    
WeChat public account: Chemocoder    
<img align="center" src="https://github.com/hcji/hcji/blob/main/img/qrcode.jpg" width="20%"/>
