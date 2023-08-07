from datetime import datetime
start = datetime.now()

import os, sys

# Data Processing
import pandas as pd
import numpy as np

from numpy import mean
# from numpy import std
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict

# set directory
script_path = "/home/inwosu/Meta_Analysis"
os.chdir(script_path)

# number of genes to use in cross validation
num_genes = [5, 10, 50, 100, 500]
# num_genes = [5, 10,]

# number of cross validation runs
num_cv_runs = 100 

# define methods
kfold = StratifiedKFold(n_splits = 5) # test diff values 
cv_method = kfold  

# some other options we may consider
# # kfold = KFold(n_splits = 3)
# (cv_method = LeaveOneOut()

# create model
model = RandomForestClassifier(n_estimators = 100, random_state = 1) # look at the default parameters

accuracy_file = open("accuracy_file_ER.txt", "w")
accuracy_file_results = open("accuracy_file_ER_result.txt", "w")

# define directories
meta_results_dir = "Data/ER_meta_results"
cross_val_data_dir = "Data/cross_validation_ER"

cv_mean_scores = []
results_file = []

for num in num_genes:

    # This outer loop determines the accuracy from list of significant genes identified by the metaanalysis
    for filename in os.listdir(meta_results_dir):    
        
        # get file name
        dataset_id = filename.replace("meta_results_without_", "")
        meta_path_name = os.path.join(meta_results_dir, "meta_results_without_" + dataset_id)
        
        # read in meta results dataframe
        meta_results = pd.read_table(meta_path_name, sep ='\t') 
        gene = meta_results['Gene']

        # select top "n" genes
        gene = gene[:num]
        
        # read in gene expression data
        cross_val_path_name = os.path.join(cross_val_data_dir, dataset_id)
        cross_val_df = pd.read_table(cross_val_path_name, sep ='\t')

        # split matrix using only top "n" genes from meta results
        X = cross_val_df[gene]
        y = cross_val_df['ER_Status'].map({'negative':0, 'positive':1})

        gene_names = X.keys()

        # evaluate model
        scores = cross_val_score(model, X, y, scoring = 'balanced_accuracy', cv = cv_method) #roc_auc
        meta_mean_score = scores.mean()
        meta_std_dev = scores.std()
        # new_scores = cross_val_predict(model, X, y, cv = cv_method, method = "predict_proba") 

        # report performance        
        meta_accuracy_value = 'Meta genes accuracy: %.3f (%.3f)' % (meta_mean_score, meta_std_dev)    
        

        # accuracy_file.writelines(dataset_id + "\t" + gene_names + "\t" + meta_accuracy_value + "\n")  
        accuracy_file.writelines(meta_path_name + "\n")
        accuracy_file.writelines(cross_val_path_name + "\n")
        accuracy_file.writelines(gene_names + " ")
        accuracy_file.writelines("\n" + meta_accuracy_value + "\n" + "\n")
        

        # This inner loop performs cross-validation and gets mean scores from a random set of genes in each dataset 
        mean_total_score = []
        total_scores = []  
        for i in range(num_cv_runs):              
            
            genes_only = cross_val_df.drop('ER_Status', axis = 1)
            
            # make sure these genes are not part of the meta analysis genes to ensure they are truly random
            new_df = genes_only.drop(columns = gene, axis = 1)
            
            # get random genes 
            X_random = new_df.sample(n = num, axis = 'columns', random_state = i)        
            y_random = cross_val_df['ER_Status']      

            # evaluate model
            random_scores = cross_val_score(model, X_random, y_random, scoring = 'balanced_accuracy', cv = cv_method) 
            random_mean_score = random_scores.mean()
            random_std_dev = random_scores.std()
            cv_mean_scores.append(random_mean_score)
            total_scores.append(random_mean_score)

            results_file.append({'Dataset_ID': dataset_id, 'num_genes': num, 'mean_score': random_mean_score})

            # report performance  
            random_accuracy_value = 'Random genes accuracy: %.3f (%.3f)' % (random_mean_score, random_std_dev)        
                
            random_gene_names = X_random.keys()

            accuracy_file.writelines("Random genes" + "\n")
            accuracy_file.writelines(random_gene_names + " ")
            accuracy_file.writelines("\n" + random_accuracy_value + "\n" + "\n")
            
            # accuracy_file.writelines(dataset_id + "\t" + gene_names + "\t" + meta_accuracy_value + "\t" + random_gene_names + "\t" + random_accuracy_value + "\t" + str(num) + "\n") 
        
        mean_total_score.append(total_scores)
        accuracy_file_results.writelines(dataset_id + "\t" + str(num) + "\t" + "Meta_genes" + "\t" + str(meta_mean_score) + "\n")
        accuracy_file_results.writelines(dataset_id + "\t" + str(num) + "\t" + "Random_genes" + "\t" + str(mean(mean_total_score)) + "\n")  

    # Create a DataFrame to store the mean scores
    results_df = pd.DataFrame(results_file)
    results_df.to_csv("results_df.tsv", sep='\t', index = False)


accuracy_file.close()
accuracy_file_results.close()


print(datetime.now()-start )



