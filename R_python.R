# This script enables the use of Python in Rstudio using the reticulate package

# install reticulate and load it (it needs to be loaded to access miniconda)
install.packages("reticulate")
library(reticulate)

# install miniconda to your required path
# you only need to do this once (the first time you want to install miniconda)
install_miniconda(path = "/inwosu/miniconda", update = T)

# reticulate creates a virtual environment as part of the miniconda installation. 
# Check that this is available. it may be in a different location from "~/miniconda/bin/conda"
conda_list(conda = "/inwosu/miniconda/bin/conda")

# If the default virtual environment (r-reticulate) is listed from the above command, 
# bind the environment to the current R session.
use_condaenv(condaenv = "r-reticulate", conda = "/inwosu/miniconda/bin/conda")

# If the default virtual environment above is not available, or you would like to create a new one
# you only need to do this once.
conda_create(envname = "myenv", conda = "/inwosu/miniconda/bin/conda")

# load the environment with ("use_condaenv")
use_condaenv(condaenv = "myenv", conda = "/inwosu/miniconda/bin/conda") 

#install python libraries to work with.
py_install(packages = c("pandas", "scikit-learn"))


