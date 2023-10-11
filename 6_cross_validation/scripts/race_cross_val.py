
from my_libs import *

# number of genes to use in cross validation
num_genes = [5, 10, 20, 50, 100, 200, 500]

# number of cross validation runs
num_cv_runs = 100 

# define methods
kfold = StratifiedKFold(n_splits = 5) # test diff values 
cv_method = kfold

# some other options we may consider
# kfold = KFold(n_splits = 3)
# (cv_method = LeaveOneOut())

# create model
model = RandomForestClassifier(n_estimators = 100, random_state = 1) # look at the default parameters

from run_all_cross_val import result_dir
accuracy_file = os.path.join(result_dir, "accuracy_file_race.txt")
accuracy_file_results = os.path.join(result_dir, "accuracy_file_result_race.txt")

accuracy_file = open(accuracy_file, "w")
accuracy_file_results = open(accuracy_file_results, "w")

# define directories
meta_results_dir = "/Data/results/race/"
cross_val_data_dir = "/Data/cross_validation_data/race/"

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
        
        # split matrix using only top "n" genes  from meta results
        X = cross_val_df[gene]
        y = cross_val_df['race'].map({'Black':0, 'White':1})

        gene_names = X.keys()

        # evaluate model
        scores = cross_val_score(model, X, y, scoring = 'balanced_accuracy', cv = cv_method)
        meta_mean_score = scores.mean()
        meta_std_dev = scores.std()
        # new_scores = cross_val_predict(model, X, y, cv = cv_method, method = "predict_proba") 

        # report performance   
        meta_accuracy_value = 'Meta genes accuracy: %.3f (%.3f)' % (meta_mean_score, meta_std_dev)

        #write results to file
        # accuracy_file.writelines(dataset_id + "\t" + gene_names + "\t" + meta_accuracy_value + "\n")  
        accuracy_file.writelines(meta_path_name + "\n")
        accuracy_file.writelines(cross_val_path_name + "\n")
        accuracy_file.writelines(gene_names + " ")
        accuracy_file.writelines("\n" + meta_accuracy_value + "\n" + "\n")

        # This inner loop performs cross-validation and gets mean scores from a random set of genes in each dataset 
        mean_total_score = []
        total_scores = []  
        
        for i in range(num_cv_runs):          
            genes_only = cross_val_df.drop('race', axis = 1)

            # make sure these gene are not part of the meta analysis genes
            new_df = genes_only.drop(columns = gene, axis = 1)

            # get random genes 
            X_random = new_df.sample(n = num, axis = 'columns', random_state = i)     
            y_random = cross_val_df['race']  

            # evaluate model
            random_scores = cross_val_score(model, X_random, y_random, scoring = 'balanced_accuracy', cv = cv_method)
            random_mean_score = random_scores.mean()
            random_std_dev = random_scores.std()
            # cv_mean_scores.append(random_mean_score)
            total_scores.append(random_mean_score)

            #report performance
            print('Random genes accuracy: %.3f (%.3f)' % (random_mean_score, random_std_dev))        
            
            random_gene_names = X_random.keys()            
            random_accuracy_value = 'Random genes accuracy: %.3f (%.3f)' % (random_mean_score, random_std_dev)
            
            accuracy_file.writelines("Random genes" + "\n")
            accuracy_file.writelines(random_gene_names + " ")
            accuracy_file.writelines("\n" + random_accuracy_value + "\n" + "\n")
            # accuracy_file.writelines(dataset_id + "\t" + gene_names + "\t" + meta_accuracy_value + "\t" + random_gene_names + "\t" + random_accuracy_value + "\t" + str(num) + "\n") 
            
            results_file.append({'Dataset_ID': dataset_id, 'num_genes': num, 'mean_score': random_mean_score})

        mean_total_score.append(total_scores)
        accuracy_file_results.writelines(dataset_id + "\t" + str(num) + "\t" + "Meta genes" + "\t" + str(meta_mean_score) + "\n")
        accuracy_file_results.writelines(dataset_id + "\t" + str(num) + "\t" + "Random genes" + "\t" + str(mean(mean_total_score)) + "\n")  
   
    # Create a DataFrame to store the mean scores
    results_df = pd.DataFrame(results_file)
    final_results = os.path.join(result_dir, "race_results_df.tsv")
    results_df.to_csv(final_results, sep='\t', index = False)  
       
accuracy_file.close()
accuracy_file_results.close()











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
# make a histogram of accuracy scores for randomly selected genes
# plot a vertical red line that shows where the accuracy score is for the metaanalysis selected genes
