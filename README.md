# Meta-analysis
A meta-analysis of curated breast cancer datasets

# Steps to follow

It is assumed that the user has Docker installed. Installation steps for Docker can be found at "https://docs.docker.com/get-docker/".

1) Create a directory, eg "Meta_Analysis" 

2) Download folders (1_prepare_data, 2_metaanalysis, 3_cross_validation as well as the "Data" folder) into this directory

3) Each step should be run in sequence, i.e 1_prepare_data should be run first, followed by 2_metaanalysis, followed by 3_cross-validation 

4) Approximate runtimes are included in the ".sh" files for each folder (Using a server running a 20 core Intel(R) Xeon(R) W-2155 CPU with 256GB of RAM Step 3 took about 3.5 days).  

5) To run each step change into tha subdirectory, e.g for step 1 change into 1_prepare_data.

6) Each step is run by executing the run_*.sh script in that step, e.g for step 2 "bash run_data_prep.sh"

