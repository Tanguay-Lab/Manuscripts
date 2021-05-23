## SUPPLEMENTAL INFORMATION
### The multi-dimensional embryonic zebrafish platform predicts flame retardant bioactivity

Lisa Truong<sup>1</sup>, Skylar Marvel<sup>2</sup>, David M. Reif<sup>2</sup>, Dennis Thomas<sup>3</sup>, Paritosh Pande<sup>4</sup>, Subham Dasgupta<sup>1</sup>, Michael T. Simonich<sup>1</sup>, Katrina M. Watersa<sup>3</sup>, and Robyn L. Tanguay<sup>1</sup>

  1. Department of Environmental and Molecular Toxicology, the Sinnhuber Aquatic Research Laboratory, and the Environmental Health Sciences Center at Oregon State University, Corvallis, OR, USA
  2. Bioinformatics Research Center, Department of Biological Sciences, North Carolina State University, Raleigh, NC, USA
  3. Biological Sciences Division, Pacific Northwest Laboratory, 902 Battelle Boulevard, P.O. Box 999, Richland, WA 99352 USA


### Table S1 – Chemical Data and Raw Data for Figures 2 and 3
* Table S1A - ChemProps and class – FRC properties, class, and sources for the 61 FRCs.
* Table S1B- LEL hit table – summary hit table used to generate the heatmap (expressed in uM)
* Table S1C- BMD – BMD values for each FRC and the 22 endpoints. 10000[SM1] indicates no bioactivity
* Table S1D - ToxCast Analysis – ToxCast data for the 45 FRCs
* Table S1E-Comparison of 12 FRCs – in vitro and in vivo LEC, POD or BMD values (log transformed)

### Table S2 – Raw Zebrafish Data 
* my_key: 
* chemical.ID (Tanguay Lab IDs used throughout)
* preferred_name – preferred name of FRC
* casrn – cas number
* pd_sample – EPA NCCT ID
* gsid – EPA ID
* morph – developmental exposure raw morphology data 
* behav_24 – EPR raw data (bin every second)
* behav_5d – LPR raw data (bin every 6 seconds)

### Table S3 – Model Input and Results
* Table S3A - Model input (LEC) – The data matrix input into the model with LEC responses
* Table S3B - Model input (BMD) - The data matrix input into the model with BMD responses
* Table S3C -  MI (BMD, LEC) – The data matrix input into the model with BMD and LEC responses
* Table S3D -  MI phys chem – The data matrix input into the model with the 13 physico-chemical properties
* Table S3E - Model Output LEC – Summary Results, and model input for the training and test set using physical properties and LEC values
* Table S3F – MO Phys Chem - Summary Results, and model input for the training and test set using the 13 physico-chemical properties
* Table S3G – MO BMD (Morph, AUC) - Summary Results, and model input for the training and test set using physical properties and BMD values for morphology and AUC calculation for LPR
* Table S3H – MO BMD (Morph, MOV) - Summary Results, and model input for the training and test set using BMD values for morphology and MOV calculation for LPR
* Table S3I – MO BMD (all) - Summary Results, and model input for the training and test set using physical properties and BMD values for morphology and both LPR endpoints
* Table S3J – MO (all) - Summary Results, and model input for the training and test set using all endpoints 
* Table S3K – MO BMD (morph only) - Summary Results, and model input for the training and test set using physical properties and BMD for zebrafish morphology only

### Table S4 – Data Dictionary

### R Code and Related Files
* SI1 – Noyes Hit Call Comparison R Script
* SI2 – Classification Model R script
