import os

# Data Processing
import pandas as pd
import numpy as np

# Modelling
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from scipy.stats import randint

# Tree Visualisation
# from sklearn.tree import export_graphviz
# from IPython.display import Image
# import graphviz

script_path = "/inwosu/Meta_Analysis"
os.chdir(script_path)

df = pd.read_table("Data/cross_validation/GSE62944_Tumor.tsv") 
#feature selection

df['race'] = df['race'].map({'Black':0,'White':1})

# # Splitting arrays or matrices into random train and test subsets
# # 70 % training dataset and 30 % test datasets
X = df.drop('race', axis=1)
y = df['race']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.30)   #use leave one out

# # creating a RF classifier
rf = RandomForestClassifier(n_estimators = 100)

# # Training the model on the training dataset
# # fit function is used to train the model using the training sets as parameters
rf.fit(X_train, y_train)

# # performing predictions on the test dataset
y_pred = rf.predict(X_test)
y_pred_prob = rf.predict_proba(X_test)

# print(y_pred)
# print(y_pred_prob)

print(len(y_pred))
print(len(y_pred_prob))
# df['predicted_race'] = y_pred

y_pred_df = pd.DataFrame(data = y_pred, columns = ['y_pred'], index = X_test.index.copy())
df_out = pd.merge(y_pred_df, df, how = 'left', left_index = True, right_index = True)
print(df_out.head)


# # # metrics are used to find accuracy or error
# from sklearn import metrics  
# # print()

# # # using metrics module for accuracy calculation
# print("ACCURACY OF THE MODEL: ", metrics.accuracy_score(y_test, y_pred))


# scores = cross_val_score(model, X, y, scoring='roc_auc', cv=cv, n_jobs=-1)











# save predictions to file
# 2 folders, cross_validation predictions, cross_validation metrices

# data_list = list()
# # directory = 'cross_validation'
# directory = 'meta_results'

# for filename in os.listdir(directory):
#     file = os.path.join(directory, filename)
#     df = pd.read_table(file)
#     data_list.append(df)

# with open("file1.txt", 'w') as output:
#     for row in data_list:
#         output.write(str(row) + '\n')


# # Initialize your model
# model = RandomForestClassifier()

# # Initialize the hold-one-out cross-validation object
# loo = LeaveOneOut()

# # Initialize a list to store the results
# results = []

# # Iterate over the datasets
# for i, dataset in enumerate(data_list):
#     X = dataset  # Assuming your dataset contains features only
#     y = labels[i]  # Assuming you have corresponding labels for each dataset
    
#     # Perform cross-validation
#     scores = cross_val_score(model, X, y, cv=loo)
    
#     # Store the mean accuracy of the model for this dataset
#     results.append(scores.mean())

# # Calculate the average performance across all datasets
# average_performance = sum(results) / len(results)
