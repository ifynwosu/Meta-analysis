from my_libs import *

result_dir = "/Data/cross_validation_results"

if not os.path.exists(result_dir):
    os.mkdir(result_dir)

# number of genes to use in cross validation
# num_genes = [5, 10, 20]                           # values to use when testing to reduce run time
num_genes = [5, 10, 20, 50, 100, 200, 500, 1000]

# number of cross validation runs
# num_cv_runs = 50                                  # values to use when testing
num_cv_runs = 100

# define methods
kfold = StratifiedKFold(n_splits = 5)
cv_method = kfold  

# create model
model = RandomForestClassifier(n_estimators = 100, random_state = 1) 

def calculate_accuracy(accuracy_file, accuracy_file_results, meta_results_dir, cross_val_data_dir, variable, val_neg, val_pos):
    accuracy_file_path = os.path.join(result_dir, accuracy_file)
    accuracy_file_results_path = os.path.join(result_dir, accuracy_file_results)

    accuracy_file = open(accuracy_file_path, "w")
    accuracy_file_results = open(accuracy_file_results_path, "w")

    results_file = []

    for num in num_genes:

        # This outer loop determines the accuracy from list of significant genes identified by the metaanalysis
        for filename in os.listdir(meta_results_dir):

            # get file name 
            dataset_id = filename.replace("meta_results_without_", "")
            meta_path_name = os.path.join(meta_results_dir, "meta_results_without_" + dataset_id)

            # read in meta results dataframe
            meta_results = pd.read_table(meta_path_name, sep='\t')
            gene = meta_results['Gene']

            # select top "n" genes
            gene = gene[:num]

            # read in gene expression data
            cross_val_path_name = os.path.join(cross_val_data_dir, dataset_id)
            cross_val_df = pd.read_table(cross_val_path_name, sep='\t')

            # split matrix using only top "n" genes from meta results
            X = cross_val_df[gene]
            y = cross_val_df[variable].map({val_neg: 0, val_pos: 1})
            gene_names = X.keys()

            # evaluate model
            scores = cross_val_score(model, X, y, scoring='balanced_accuracy', cv=cv_method)
            meta_mean_score = scores.mean()
            meta_std_dev = scores.std()

            # report performance
            meta_accuracy_value = 'Meta genes accuracy: %.3f (%.3f)' % (meta_mean_score, meta_std_dev)

            # write results to file
            accuracy_file.writelines(meta_path_name + "\n")
            accuracy_file.writelines(cross_val_path_name + "\n")
            accuracy_file.writelines(gene_names + " ")
            accuracy_file.writelines("\n" + meta_accuracy_value + "\n" + "\n")

            # This inner loop performs cross-validation and gets mean scores from a random set of genes in each dataset
            mean_total_score = []
            total_scores = []

            for i in range(num_cv_runs):
                genes_only = cross_val_df.drop(variable, axis=1)

                # make sure these genes are not part of the meta analysis genes to ensure they are truly random
                new_df = genes_only.drop(columns=gene, axis=1)

                # if number of random genes to use in cross validation is greater than number of genes in input dataframe, skip
                # if num > new_df.shape[1]:
                #     continue

                X_random = new_df.sample(n=num, axis='columns', random_state=i)
                y_random = cross_val_df[variable]

                # evaluate model
                random_scores = cross_val_score(model, X_random, y_random, scoring='balanced_accuracy', cv=cv_method)
                random_mean_score = random_scores.mean()
                random_std_dev = random_scores.std()
                total_scores.append(random_mean_score)

                # report performance
                print(dataset_id, variable, " number of genes:", num, " iteration:", i, ' Random genes accuracy: %.3f (%.3f)' % (random_mean_score, random_std_dev))                   

                random_gene_names = X_random.keys()
                random_accuracy_value = 'Random genes accuracy: %.3f (%.3f)' % (random_mean_score, random_std_dev)

                accuracy_file.writelines("Random genes" + "\n")
                accuracy_file.writelines(random_gene_names + " ")
                accuracy_file.writelines("\n" + random_accuracy_value + "\n" + "\n")

                results_file.append({'Dataset_ID': dataset_id, 'num_genes': num, 'mean_score': random_mean_score})

            mean_total_score.append(total_scores)
            accuracy_file_results.writelines(dataset_id + "\t" + str(num) + "\t" + "Meta-analysis genes" + "\t" + str(meta_mean_score) + "\n")
            accuracy_file_results.writelines(dataset_id + "\t" + str(num) + "\t" + "Random genes" + "\t" + str(mean(mean_total_score)) + "\n")

        # Create a DataFrame to store the mean scores
        results_df = pd.DataFrame(results_file)
        final_results = os.path.join(result_dir, (variable + ".tsv"))
        results_df.to_csv(final_results, sep='\t', index=False)

    accuracy_file.close()
    accuracy_file_results.close()


# Define the parameters for each call
calls = [    
    {   "accuracy_file": "accuracy_file_ER_status.txt",
        "accuracy_file_results": "accuracy_file_result_ER_status.txt",
        "meta_results_dir": "/Data/metaanalysis_results/ER_status/",
        "cross_val_data_dir": "/Data/cross_validation_data/ER_status/",
        "variable": "ER_status",
        "val_neg": "negative",
        "val_pos": "positive"    
    },
    {   "accuracy_file": "accuracy_file_PR_status.txt",
        "accuracy_file_results": "accuracy_file_result_PR_status.txt",
        "meta_results_dir": "/Data/metaanalysis_results/PR_status/",
        "cross_val_data_dir": "/Data/cross_validation_data/PR_status/",
        "variable": "PR_status",
        "val_neg": "negative",
        "val_pos": "positive"
    },
    {   "accuracy_file": "accuracy_file_HER2_status.txt",
        "accuracy_file_results": "accuracy_file_result_HER2_status.txt",
        "meta_results_dir": "/Data/metaanalysis_results/HER2_status/",
        "cross_val_data_dir": "/Data/cross_validation_data/HER2_status/",
        "variable": "HER2_status",
        "val_neg": "negative",
        "val_pos": "positive"
    },
    {   "accuracy_file": "accuracy_file_race.txt",
        "accuracy_file_results": "accuracy_file_result_race.txt",
        "meta_results_dir": "/Data/metaanalysis_results/race/",
        "cross_val_data_dir": "/Data/cross_validation_data/race/",
        "variable": "race",
        "val_neg": "Black",
        "val_pos": "White" 
    },
    {   "accuracy_file": "accuracy_file_tri_neg_status.txt",
        "accuracy_file_results": "accuracy_file_result_tri_neg_status.txt",
        "meta_results_dir": "/Data/metaanalysis_results/tri_neg_status/",
        "cross_val_data_dir": "/Data/cross_validation_data/tri_neg_status/",
        "variable": "tri_neg_status",
        "val_neg": "tri_neg",
        "val_pos": "non_tri_neg"
    },    
    {   "accuracy_file": "accuracy_file_race_tri_neg.txt",
        "accuracy_file_results": "accuracy_file_result_race_tri_neg.txt",
        "meta_results_dir": "/Data/metaanalysis_results/race_tri_neg/",
        "cross_val_data_dir": "/Data/cross_validation_data/race_tri_neg/",
        "variable": "race",
        "val_neg": "Black",
        "val_pos": "White"    
    }    
]

# Call the function for each set of parameters
for call_params in calls:
    calculate_accuracy(**call_params)

