# Deep autoencoder-based behavioral pattern recognition outperforms standard statistical methods in high-dimensional zebrafish studies 
Deep autoencoders utilizing raw behavioral tracking data from zebrafish larvae outperforms traditional statistical methodologies resulting in a comprehensive evaluation of behavioral data following toxicant exposure. Our models can accurately distinguish between normal and abnormal behavior with near complete overlap with existing statistical approaches, with many chemicals detectable at lower concentrations than with traditional statistical tests – a critical finding relevant to protecting public health. 

![Autoencoder Model Design](https://github.com/ajgreen4/auto-behavior-zt/blob/main/Autoencoder%20Model%20Design.png)


This repository contains ten Jupyter Notebook [files](https://github.com/Tanguay-Lab/Manuscripts/tree/main/Green_et_al_(2022)_Manuscript/Files):
 - Data preprocessing notebooks
    - etho_data_Step_1_QC.ipynb
        - Removes data for dead or abnormal larvae as well as incomplete data
    - etho_data_Step_2_standardize_tabular_data.ipynb
        - 5dpf movement dataset with all timepoints in the second light cycle using QC'd data from step 1
    - etho_data_Step_3_control_threshold.ipynb
        - Standardized data from step 2 and remove outliers
    - etho_data_Step_4_tabular_data.ipynb
        - Create 3D numpy data tables and split movement data and experimental details
 - Training behavioral autoencoder notebooks
    - Train models using dark phase data
        - etho_data_Step_5_train_dark_low_zf_Autoencoder.ipynb
            - model uses hypoactive control data
        - etho_data_Step_5_train_dark_medium_zf_Autoencoder.ipynb
            - model uses normal control data
        - etho_data_Step_5_train_dark_high_zf_Autoencoder.ipynb
            - model uses hyperactive control data
    - Train models using light phase data
        - etho_data_Step_5_train_dark_low_zf_Autoencoder.ipynb
            - model uses hypoactive control data
        - etho_data_Step_5_train_dark_medium_zf_Autoencoder.ipynb
            - model uses normal control data
        - etho_data_Step_5_train_dark_high_zf_Autoencoder.ipynb
            - model uses hyperactive control data
 - Local routines
    - tkmetrics.py
    - Green_scripts_v1.py
- high_dark_sample_z_data.zip
    - Sample training data that may be used to run the etho_data_Step_5_train_dark_high_zf_Autoencoder.ipynb notebook.
    - Complete dataset available apon request.

 These files should be run in the Singularity container environment build for this project to ensure all dependencies are present and that the code runs error free. The container can be found at: https://cloud.sylabs.io/library/ajgreen/default/ml-nvdashboard-container

This container can be run with port forwarding on a remote high-performance computing cluster using the following code slurm:

```basg
srun -p gpu -w node94 -J jupAuto -c 4 --mem=50000 --pty bash

# get tunneling info
XDG_RUNTIME_DIR=""
port=8194
# port=$(shuf -i8000-9999 -n1)
# GPU node
node=--IP Address--
user=$(whoami)
cluster=--cluster address--

# print tunneling instructions jupyter-log
echo -e "
MacOS or linux terminal command to create your ssh tunnel
ssh -N -L ${port}:${node}:${port} ${user}@${cluster}

Windows MobaXterm info for Tools->ModaSSHTunnel
Forwarded port: ${port}
Remote server: ${node}
Remote port: ${port}
SSH server: ${cluster}
SSH login: ${user}
SSH port: 22

Use a Browser on your local machine to go to:
http://localhost:${port}  (prefix w/ https:// if using password)

DON'T USE ADDRESS BELOW.
DO USE TOKEN BELOW

"

nvidia-modprobe -u -c=0
singularity exec --nv tensorflow2-gpu-r-jupyter-nvdashboard.sif jupyter lab --no-browser --port=${port} --ip=${node}
```

