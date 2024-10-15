# aind-spot-spectral-unmixing

This approach has been developed by Michael Xie and Tim Wang at Janelia, and modified by Hannah Schryver and Carson Berry at the Allen Institute. 

It is part of the HCR processing pipeline implemented in Code Ocean. 

The goal is to assign correct identities to all detected spots, and output a cell-by-gene table. 

The inputs to this capsule are the outputs from spot detection (https://github.com/AllenNeuralDynamics/aind-z1-spot-detection) and multichannel spot evaluation (https://github.com/AllenNeuralDynamics/aind-spot-get-multichannel). 

It works in a few steps: 

    1. Load the data 
        For each channel: 
            A. Detected spots
            B. Intensity information of all other channels at detected spot locations 
            
    2. Calculate Intensities and apply percentile threshold filtering 
    
    3. Calculate multidimensional "dye-lines" (aka ratios) from the intensity of each spots in each channel. 
    These are used later in assessing how the likelihood that a given spot really belongs to a given channel. 
    
    4. Calculate the distance every spot is from each dye-line and produce some other stats
    
    5. Apply some basic QC filters to the spots (roundness, distance to other spots, intensity)
    
    6. Find the closest dye-line for each spot and eliminate any spots that might be double counts. 
    
    7. Generate cell-by-gene table!

