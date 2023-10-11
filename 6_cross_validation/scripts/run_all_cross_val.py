# this script just runs all the cross validation scripts at once

from my_libs import *

result_dir = "/Data/cross_validation_results"

if not os.path.exists(result_dir):
    os.mkdir(result_dir)
    
from race_cross_val import *
from HER2_cross_val import *
from PR_cross_val import *
from ER_cross_val import *












# set directory
# script_path = "/home/inwosu/Meta_Analysis"
# os.chdir(script_path)

# use this for testing
# num_genes = [5, 10]
# num_cv_runs = 3

# this snipppet prints the present working directory 
# a = os.getcwd()
# print(a)

# print(datetime.now()-start)