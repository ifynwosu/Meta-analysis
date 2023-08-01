
import os, sys

# Data Processing
import pandas as pd
import numpy as np

from numpy import mean
from numpy import std
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_predict

# set directory
script_path = "/home/inwosu/Meta_Analysis"
os.chdir(script_path)

# number of genes to use in cross validation
n = 5

# define methods
kfold = StratifiedKFold(n_splits = 5) # test diff values 
cv_method = kfold
# kfold = KFold(n_splits = 3)
# cv_method = LeaveOneOut()

# create model
model = RandomForestClassifier(n_estimators = 100, random_state = 1) # look at the default parameters
# model = LogisticRegression(max_iter = 1000, random_state = 1)

accuracy_file = open("accuracy_file_ER.txt", "w")


# define directories
meta_results_dir = "Data/ER_meta_results"
cross_val_data_dir = "Data/cross_validation_ER"
# filename = os.listdir(meta_results_dir)
# print(filename)
# sys.exit()

for filename in os.listdir(meta_results_dir):

    # This first section determines the accuracy from list of significant genes identified by the metaanalysis
    # get file name
    filename = filename.replace("meta_results_without_", "")
    meta_path_name = os.path.join(meta_results_dir, "meta_results_without_" + filename)
    
    # read in data
    meta_results = pd.read_table(meta_path_name, sep ='\t') 
    gene = meta_results['Gene']

    # select top "n" genes
    gene = gene[:n]
    
    # read in data
    cross_val_path_name = os.path.join(cross_val_data_dir, filename)
    print(cross_val_path_name)
    cross_val_df = pd.read_table(cross_val_path_name, sep ='\t')

    # split matrix using only top "n" genes
    X = cross_val_df[gene]
    y = cross_val_df['ER_Status'].map({'negative':0, 'positive':1})

    gene_names = X.keys()

    # evaluate model
    # scores = cross_val_score(model, X, y, scoring = 'balanced_accuracy', cv = cv_method)
    scores = cross_val_score(model, X, y, scoring = 'roc_auc', cv = cv_method)
    print(scores)

    new_scores = cross_val_predict(model, X, y, cv = cv_method, method = "predict_proba")
    print(new_scores)
    # sys.exit()

    # report performance
    print('Meta genes accuracy: %.3f (%.3f)' % (mean(scores), std(scores)))
    meta_accuracy_value = 'Meta genes accuracy: %.3f (%.3f)' % (mean(scores), std(scores))

    accuracy_file.writelines(meta_path_name + "\n")
    accuracy_file.writelines(cross_val_path_name + "\n")
    accuracy_file.writelines(gene_names + " ")
    accuracy_file.writelines("\n")
    accuracy_file.writelines(meta_accuracy_value + "\n" + "\n")

# This second section determines the accuracy from a random list of genes in each dataset
    for i in range(5):  
    
        genes_only = cross_val_df.drop('ER_Status', axis = 1)
        
        # make sure thesse gene are not part of the meta analysis genes
        new_df = genes_only.drop(columns = gene, axis = 1)
        
        X_random = new_df.sample(n = n, axis = 'columns', random_state = i) # write selected genes to file        
        y_random = cross_val_df['ER_Status']      

        # evaluate model
        # random_scores = cross_val_score(model, X_random, y_random, scoring = 'balanced_accuracy', cv = cv_method) 
        random_scores = cross_val_score(model, X_random, y_random, scoring = 'roc_auc', cv = cv_method) 

        # # report performance
        print('Random genes accuracy: %.3f (%.3f)' % (mean(random_scores), std(random_scores)))
        random_accuracy_value = 'Random genes accuracy: %.3f (%.3f)' % (mean(random_scores), std(random_scores))
        
        random_gene_names = X_random.keys()
        
        accuracy_file.writelines("Random genes" + "\n")
        accuracy_file.writelines(random_gene_names + " ")
        accuracy_file.writelines("\n")
        accuracy_file.writelines(random_accuracy_value + "\n" + "\n")
        print(X_random.keys())
        # sys.exit()


accuracy_file.close()






