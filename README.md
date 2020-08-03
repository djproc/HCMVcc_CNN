# Overview:

#TODO: add overview figure

This ia a several step image analysis pipeline:
1. Identidy nuclei in dataset (use either CellProfiler or StarDist). 
2. Extract RGB images of single nuclei (only using Hoechst and HCMV gB channels)
3. Use a pre-trained convolutional neural network to clasify images.
4. Sort images based on classification confiedence into folders for manual check of classification accuracy.
5. Choose folders to include in the analysis, then collect line scans (150 pixel wide) across nucleus for all cells.
6. Determine mean intensity and SEM over three biologically independent replicates. Vizualize data.  

# Getting started:

#TODO: Add description of setting up environment (or add a docker hub link). 

Open the GenerateIDs_Drawlinescans file and follow the instructions there 
