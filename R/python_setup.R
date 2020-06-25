install.packages("reticulate")
library(reticulate)
Sys.which("python") # find out what version of python R has found
py_discover_config() # same thing basically
use_python("/usr/bin/python3") # force it to use python3
py_config() # check that it worked 
py_install("pandas") # get this error: "could not find a Python environment for /usr/local/bin/python3". need to create an environment with either virtualenv or conda

# sounds like changes to Python packages, etc. can either be made globally for the entire computer using Terminal or Virtualcode, or locally using Conda. Let's use Conda for now, it seems more bulletproof / allows R customization?

conda_create("r-reticulate") # I think this created an environment in "/Users/afh/opt/anaconda3/envs/r-reticulate/bin/python" that can be accessed from any R window? 

# RUN IN TERMINAL: conda activate r-reticulate
# now in terminal instaed of saying (python) (base) it says (r-reticulate) (python)

use_condaenv("r-reticulate") # confirm that we want to use the conda installations in r-reticulate

# none of these instalations worked for me
# kept getting an error about a file environments.txt that exists already: 
# Verifying transaction: / WARNING conda.core.path_actions:verify(963): Unable to create environments file. Path not writable.
# environment location: /Users/afh/.conda/environments.txt
# 
# done
# Executing transaction: - WARNING conda.core.envs_manager:register_env(52): Unable to register environment. Path not writable or missing.
# environment location: /Users/afh/opt/anaconda3
# registry file: /Users/afh/.conda/environments.txt
conda_install("r-reticulate", "numpy")
conda_install("r-reticulate", "scipy")
conda_install("r-reticulate", "pandas")
conda_install("r-reticulate", "matplotlib")

# install via terminal instead: conda install SciPy
# that didn't work either 

pd <- import("pandas")
np <- import("numpy")
plt <- import("matplotlib.pyplot")
sci <- import("scipy")
# all of these except numpy return: 
# Error in py_module_import(module, convert = convert) : 
#   ModuleNotFoundError: No module named 'scipy'