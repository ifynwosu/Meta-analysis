import os, sys

# Data Processing
import pandas as pd
import numpy as np

from numpy import mean
from numpy import std
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold

# look at pandas read_table doc for stripping
# set directory
script_path = "/home/inwosu/Meta_Analysis"
os.chdir(script_path)

# number of genes to use in cross validation
n = 100

# define methods
# kfold = KFold(n_splits = 3)
kfold = StratifiedKFold(n_splits = 3) # test diff values 
cv_method = kfold
# cv_method = LeaveOneOut()

# create model
model = RandomForestClassifier(random_state = 1) # look at the default parameters

# This first section determines the accuracy from list of significant genes identified by the metaanalysis

# read in data 
# meta_results = pd.read_table("/inwosu/Meta_Analysis/Data/meta_results/meta_results_without_E_TABM_158.tsv", sep ='\t') 
meta_results = pd.read_table("Data/meta_results/meta_results_without_GSE5847.tsv", sep ='\t')  
# meta_results.columns = meta_results.columns.str.strip()
gene_column = meta_results['Gene']
print(gene_column)

# select top n genes
# gene_column = gene_column[:n]
# gene_column = gene_column.str.strip()

# cross_val_df = pd.read_table("/inwosu/Meta_Analysis/Data/cross_validation_full/E_TABM_158.tsv")
# cross_val_df = pd.read_table("Data/cross_validation_full/GSE5847.tsv")
# cross_val_df.columns = cross_val_df.columns.str.strip()

# split matrix
# X = cross_val_df[gene_column]
# y = cross_val_df['race'].str.strip().map({'Black':0, 'White':1})
# print(X)
# print(y)

# sys.exit()
# evaluate model
# scores = cross_val_score(model, X, y, scoring = 'balanced_accuracy', cv = cv_method)

# report performance
# print('Meta_genes Accuracy: %.3f (%.3f)' % (mean(scores), std(scores)))

# accuracy_file = open("accuracy_file.txt", "w")

# read in dataset  
# cross_val_df = pd.read_table("Data/cross_validation/E_TABM_158.tsv")
# cross_val_df = pd.read_table("Data/cross_validation_full/E_TABM_158.tsv") 
# cross_val_df.columns = cross_val_df.columns.str.strip()
# cross_val_df['race'] = cross_val_df['race'].str.strip().map({'Black':0, 'White':1})

# for i in range(10): 
    
#     y = cross_val_df['race']
#     new_cross_df = cross_val_df.drop('race', axis = 1)
#     X = new_cross_df.sample(n = n, axis = 'columns', random_state = i)

#     # predict_df = pd.merge(race_col, random_df, how = "inner", left_index = True, right_index = True)

#     # X = predict_df.drop('race', axis = 1)
#     # y = predict_df['race']

#     # # evaluate model
#     scores = cross_val_score(model, X, y, scoring = 'balanced_accuracy', cv = cv_method)

#     # # report performance
#     print('Accuracy: %.3f (%.3f)' % (mean(scores), std(scores)))
#     accuracy_value = ['Accuracy: %.3f (%.3f)' % (mean(scores), std(scores)), "\n"]
    # accuracy_file.writelines(accuracy_value)
    
# accuracy_file.close()










# y_true, y_pred = list(), list()
# for train_ix, test_ix in cv.split(X):
#   X_train = X[train_ix]
#   X_test = X[test_ix]
#   y_train, y_test = y[train_ix], y[test_ix]
#   model = RandomForestClassifier(random_state=1)
#   model.fit(X_train, y_train)
#   yhat = model.predict(X_test)
# #   y_true.append(y_test[0])
#   y_pred.append(yhat[0])
#     #break

# acc = accuracy_score(y_true, y_pred)
# print('Accuracy: %.3f' % acc)




# save predictions and actual values + auroc diff file
# randomly select x(100 - 500) number of genes (100 times)
# set random seed before for loop (set to 0)
# make a histogram of accuracy scores for randomly selected genes
# plot a vertical red line that shows where the accuracy score is for the metaanalysis selected genes
