### Input File Format

#### MGF File

DeepMASS accepts **mgf** files containing one or multiple MS/MS spectra. The Mascot Generic Format (MGF) 
is a standard file format used to store mass spectrometry (MS) data, particularly tandem mass spectrometry (MS/MS) data.    

Structure of MGF Files:    

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


#### MSP File

DeepMASS also accepts **msp** files, which contain one or multiple MS/MS spectra.    

Structure of MSP Files:    

- **Name**: A descriptor or identifier for the spectrum, often the compound name.
- **Num Peaks**: The number of mass-to-charge ratio (m/z) and intensity pairs that follow.
- **Peak Data**: Describes the mass spectrum with m/z and intensity pairs listed.

Essential DeepMASS input information includes:    

- **CHARGE:** Number of charges of the precursor ion.
- **PRECURSOR_MZ:** Precursor ion's mass-to-charge ratio.
- **COMPOUND_NAME:** Identifier for the compound.
- **IONMODE:** Ionization mode of the precursor ion.
- **ADDUCT:** Type of adduct of the precursor ion.

At least one of the following is necessary for DeepMASS inputs:    

- **FORMULA:** Molecular formula of the compound.
- **PARENT_MASS:** Monoisotopic mass of the compound.

Here's a simplified example of an MSP file with multiple spectra:   

        Name: challenge_255
        PrecursorMZ: 194.04588176800002
        ExactMass: 195.053157768
        Formula: C9H9NO4
        Ion_mode: negative
        Adduct: [M-H]-
        Num Peaks: 5
        132.0454 0.0060
        134.0245 0.0022
        149.0483 0.0030
        150.0560 1.0
        194.0459 0.3692
        
        Name: challenge_368
        PrecursorMZ: 223.970551912
        ExactMass: 222.963275912
        Formula: C9H6BrNO
        Ion_mode: positive
        Adduct: [M+H]+
        Num Peaks: 4
        127.0418 0.0139
        145.0524 0.0058
        205.9602 0.1043
        223.9708 1.0


#### MAT File

DeepMASS accepts **mat** files which contain detailed information about MS/MS spectra. The file format is 
the same as *MS-Finder* used, which also contains isotopical information.    

Here's a simplified example of an MAT file with multiple spectra:  

        NAME: L-Glutamic acid
        SCANNUMBER: 1228
        RETENTIONTIME: 5.8845
        PRECURSORMZ: 148.0424
        PRECURSORTYPE: [M+H]+
        IONMODE: Positive
        SPECTRUMTYPE: Profile
        INTENSITY: 10818.72
        FORMULA: C5H9NO4
        ONTOLOGY: 
        INCHIKEY: WHUUTDBJXJRKMK-VKHMYHEASA-N
        SMILES: 
        MSTYPE: MS1
        Num Peaks: 5
        148.04244	17967
        149.04424	1255
        150.03773	1096
        151.04014	113
        152.07201	68
        MSTYPE: MS2
        Num Peaks: 7
        44.97939	304
        61.00939	32
        73.78966	70
        74.00583	7225
        76.82613	13
        87.01311	80
        120.01151	988


### Output File Format

The output is CSV file contains information about ranked annotated chemical structures. 
The columns provide details about each compound, including its title, molecular formula, 
structure, and scores related to formula and structure predictions.

Columns Description:   

- **Title**: The name of the compound.
- **MolecularFormula**: The molecular formula of the compound, indicating the number and type of atoms present.
- **CanonicalSMILES**: The Canonical Simplified Molecular Input Line Entry System (SMILES) string, representing the structure of the compound.
- **InChIKey**: The International Chemical Identifier Key, a unique alphanumeric string to identify the compound.
- **Database IDs**: Identifiers from various databases where the compound is registered.
- **Formula Score**: A score indicating the accuracy or confidence in the predicted molecular formula.
- **Structure Score**: A score indicating the accuracy or confidence in the predicted molecular structure.
- **Consensus Score**: An overall score combining the formula and structure scores, reflecting the overall confidence in the compound identification.

Here's a simplified example of an output file:  

|  | Title                                | MolecularFormula | CanonicalSMILES     | InChIKey                     | Database IDs                               | Formula Score | Structure Score     | Consensus Score       |
|--|--------------------------------------|------------------|---------------------|------------------------------|--------------------------------------------|----------------|----------------------|------------------------|
| 0 | "Acetic acid N',N'-dimethylhydrazide" | C4H10N2O          | CC(=O)NN(C)C         | SLIKWWJXVUHCPJ-UHFFFAOYSA-N  | NIST:CAS6233041                            | 1.0            | 0.9999999991779194   | 0.9999999994245434     |
| 1 | N-(2-Aminoethyl)acetamide            | C4H10N2O          | CC(=O)NCCN           | DAKZISABEDGGSV-UHFFFAOYSA-N  | PMHub:MS000036761                          | 1.0            | 0.957735969543063    | 0.970415178680144      |
| 2 | N-Propylurea                         | C4H10N2O          | CCCNC(N)=O           | ZQZJKHIIQFPZCS-UHFFFAOYSA-N  | BloodExp:CID12303;;HMDB0255240             | 1.0            | 0.940423010788881    | 0.9582961075522167     |
| 3 | "2-Amino-N,N-dimethylacetamide"      | C4H10N2O          | CN(C)C(=O)CN         | KNVRBEGQERGQRP-UHFFFAOYSA-N  | NIST:CAS1857198                            | 1.0            | 0.8933594423863824   | 0.9253516096704677     |
| 4 | 4-Aminobutanamide                    | C4H10N2O          | NCCCC(N)=O           | WCVPFJVXEXJFLB-UHFFFAOYSA-N  | NIST:CAS3251089                            | 1.0            | 0.8222707059586819   | 0.8755894941710773     |
| 5 | "1,1,3-Trimethylurea"                | C4H10N2O          | CNC(=O)N(C)C         | COSWCAGTKRUTQV-UHFFFAOYSA-N  | NIST:CAS632144                             | 1.0            | 0.7851222117039512   | 0.8495855481927659     |
| 6 | Methylpropylnitrosamine              | C4H10N2O          | CCCN(C)N=O           | ITBDKUCVKYSWMF-UHFFFAOYSA-N  | BloodExp:CID13545                          | 1.0            | 0.7422341116397019   | 0.8195638781477912     |
| 7 | (2S)-2-(Methylamino)propanamide      | C4H10N2O          | CN[C@@H](C)C(N)=O    | QKNFFJHHPCWXTH-VKHMYHEASA-N  | NIST:CAS55988120                           | 1.0            | 0.7230225309588504   | 0.8061157716711953     |
| 8 | N-nitrosodiethylamine                | C4H10N2O          | CCN(CC)N=O           | WBNQDOYYEUMPFS-UHFFFAOYSA-N  | CHEBI:34873;;BloodExp:CID5921;;C14422;;PMHub:MS000007343;;DTXSID2021028 | 1.0 | 0.7174501691978581 | 0.8022151184385007     |
| 9 | ISOPROPYLUREA                        | C4H10N2O          | CC(C)NC(=N)O         | LZMATGARSSLFMQ-UHFFFAOYSA-N  | PMHub:MS000231594;;NIST:CAS691601          | 1.0            | 0.6049047550645834   | 0.7234333285452084     |


### Export MS/MS file from MS-DIAL

To export an MGF (Mascot Generic Format) file from MS-DIAL processing results, follow these steps:

[![ms-dial-1.png](https://i.postimg.cc/CKVG1pyc/ms-dial-1.png)](https://postimg.cc/RWX6s8Zt)

A) Peak list export
B) Alignment result export
C) Molecular spectrum networking export
D) Copy screenshot to clipboard (emf)
E) Parameter export (Tab-delimited text) F) Export as lipoquality database format
G) Export normalization result

[![ms-dial-2.png](https://i.postimg.cc/c6yHmZBj/ms-dial-2.png)](https://postimg.cc/GTKdLWNJ)

A) Peak list export: You can get the peak list information of each sample including retention time, m/z, MS/MS spectra information, and so on. Available formats are MSP, MGF or Text.

Step1. Choose an export folder path.
Step2. Choose files which you want to export and click button "Add ->".
Step3. Select export format.
Step4. Click the export button.

[![ms-dial-3.png](https://i.postimg.cc/4NVbTxqY/ms-dial-3.png)](https://postimg.cc/SJQz6ky4)

B) Alignment result export: You can get data matrix or spectral information.

Step1. Choose an export folder path.
Step2. Choose an alignment file which you want to export.
Step3. Select export format if you want to export the representative spectra.
Step4. Click the export button.


### Export MS/MS file from XCMS

To export MS/MS data into an MGF (Mascot Generic Format) file from XCMS result, we can use the combineSpectra function. 
For more details see the documentation of the consensusSpectrum function in the MSnbase R package. Refer this script to define 
a function exportSpectraToMGF to convert the spectra data into MGF format.
 
        exportSpectraToMGF <- function(spectra, file) {
            mgf_data <- lapply(spectra, function(sp) {
                list(
                    TITLE = paste("Scan", sp@scanIndex),
                    RTINSECONDS = sp@rtime,
                    PEPMASS = c(sp@precursorMz, sp@precursorIntensity),
                    CHARGE = sp@precursorCharge,
                    MZ = sp@mz,
                    INTENSITY = sp@intensity
                )
            })
            con <- file(file, "w")
            for (spectrum in mgf_data) {
                cat("BEGIN IONS\n", file = con)
                cat(paste("TITLE=", spectrum$TITLE, sep = ""), "\n", file = con)
                cat(paste("RTINSECONDS=", spectrum$RTINSECONDS, sep = ""), "\n", file = con)
                cat(paste("PEPMASS=", paste(spectrum$PEPMASS, collapse = " "), sep = ""), "\n", file = con)
                cat(paste("CHARGE=", spectrum$CHARGE, sep = ""), "\n", file = con)
                for (i in seq_along(spectrum$MZ)) {
                    cat(paste(spectrum$MZ[i], spectrum$INTENSITY[i], sep = " "), "\n", file = con)
                }
                cat("END IONS\n\n", file = con)
            }
            close(con)
        }
