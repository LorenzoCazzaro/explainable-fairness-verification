from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import Dataset

import numpy as np
import pandas as pd
import json

#part of the code is inspired by the released code of "Fairness-Aware Training of Decision Trees by Abstract Interpretation" by Francesco Ranzato, Caterina Urban, and Marco Zanella (https://github.com/fatt21/fatt)

def preproc_adult(feats_to_remove = []): #45222 88
	train_file = "../adult/adult.data"
	test_file = "../adult/adult.test"

	column_names = [
	'age', 'workclass', 'fnlwgt', 'education', 'education_num', 'marital_status', 'occupation', 'relationship',
	'race', 'sex', 'capital_gain', 'capital_loss', 'hours_per_week', 'native_country', 'income'
	]
	
	training_set = pd.read_csv(train_file, sep=',', header=None, names=column_names)
	test_set = pd.read_csv(test_file, sep=',', header=None, names=column_names)

	training_set = training_set.applymap(lambda x: x.strip() if isinstance(x, str) else x)
	test_set = test_set.applymap(lambda x: x.strip() if isinstance(x, str) else x)
	training_set.replace(to_replace='?', value=np.nan, inplace=True)
	test_set.replace(to_replace='?', value=np.nan, inplace=True)
	training_set.dropna(axis=0, inplace=True)
	test_set.dropna(axis=0, inplace=True)
	training_set = training_set.drop("education", 1)
	test_set = test_set.drop("education", 1)

	binary = ['sex']
	categorical = ['workclass', 'marital_status', 'occupation', 'relationship', 'race', 'native_country']
	for c in feats_to_remove:
		if c in categorical:
			categorical.remove(c)
		column_names.remove(c)

	training_set = training_set.drop(feats_to_remove, 1)
	test_set = test_set.drop(feats_to_remove, 1)
		
	training_set.replace({'<=50K': 0, '>50K': 1}, inplace=True)
	test_set.replace({'<=50K.': 0, '>50K.': 1}, inplace=True)

	numerical = Dataset.numericalColumns(training_set, binary, categorical, ['income'])
	training_set = Dataset.normalizeColumnNames(training_set)
	training_set = Dataset.normalizeColumnValues(training_set, binary + categorical)
	training_set = Dataset.oneHotEncoding(training_set, binary, categorical)
	training_set[numerical] = training_set[numerical].values.astype(np.float32)
	training_set = Dataset.selectLabelColumn(training_set, 'income')
	training_set.loc[training_set["income"]==0, 'income'] = -1
	
	test_set = Dataset.normalizeColumnNames(test_set)
	test_set = Dataset.normalizeColumnValues(test_set, binary + categorical)
	test_set = Dataset.oneHotEncoding(test_set, binary, categorical)
	test_set[numerical] = test_set[numerical].values.astype(np.float32)
	test_set = Dataset.selectLabelColumn(test_set, 'income')
	test_set.loc[test_set["income"] == 0, "income"] = -1

	binary = ['sex_male']

	if "native_country" in column_names:
		test_set.insert(loc=training_set.columns.get_loc('native_country=holand_netherlands'), column='native_country=holand_netherlands', value=0)
	
	dataset = pd.concat([training_set, test_set], sort=False)

	norm_info = []
	maxs = dataset.max(axis=0)
	mins = dataset.min(axis=0)
	avgs = dataset.mean(axis=0)
	stds = dataset.std(axis=0)
	for i in range(dataset.shape[1]):
		norm_info += [[float(mins[i]), float(maxs[i]), avgs[i], stds[i]]]
	norm_info = norm_info[1:]

	training_set = Dataset.normalize(training_set, numerical)
	test_set = Dataset.normalize(test_set, numerical)
	dataset = pd.concat([training_set, test_set], sort=False)

	categorical_columns_indices = []
	categorical_columns_names = []
	for categorical_column in categorical:
		indices = [i-1 for i, x in enumerate(dataset.columns) if categorical_column in x]
		column_names = [name for name in dataset.columns if categorical_column in name]
		categorical_columns_indices += [indices]
		categorical_columns_names += [column_names]

	numerical_binary_columns = numerical + binary
	numerical_binary_columns_indices = [i-1 for i, x in enumerate(dataset.columns) if x in numerical_binary_columns]
	numerical_binary_columns_names = [name for name in dataset.columns if name in numerical_binary_columns]

	suffix_filename = "" if len(feats_to_remove) == 0 else "-"
	for c in feats_to_remove:
		suffix_filename += (c + "_")
	if len(suffix_filename) != 0:
		suffix_filename = suffix_filename[:-1]

	Dataset.save(dataset, "../adult/dataset" + suffix_filename + ".csv")
	Dataset.save(training_set, "../adult/training-set" + suffix_filename + ".csv")
	Dataset.save(test_set, "../adult/test-set" + suffix_filename + ".csv")
	Dataset.exportColumns(dataset, "../adult/columns" + suffix_filename + ".csv")
	dataset_column_names_file = open("../adult/ad_column_names" + suffix_filename + ".json", "w")
	dataset_column_names_file.write(json.dumps(dataset.columns.tolist()[1:]))
	dataset_column_names_file.close()
	dataset_categorical_column_names_file = open("../adult/ad_categorical_column_names" + suffix_filename + ".json", "w")
	dataset_categorical_column_names_file.write(json.dumps(categorical_columns_names))
	dataset_categorical_column_names_file.close()
	dataset_categorical_column_index_file = open("../adult/ad_categorical_column_index" + suffix_filename + ".json", "w")
	dataset_categorical_column_index_file.write(json.dumps(categorical_columns_indices))
	dataset_categorical_column_index_file.close()
	dataset_numerical_binary_column_names_file = open("../adult/ad_numerical_binary_column_names" + suffix_filename + ".json", "w")
	dataset_numerical_binary_column_names_file.write(json.dumps(numerical_binary_columns_names))
	dataset_numerical_binary_column_names_file.close()
	dataset_numerical_binary_column_index_file = open("../adult/ad_numerical_binary_column_index" + suffix_filename + ".json", "w")
	dataset_numerical_binary_column_index_file.write(json.dumps(numerical_binary_columns_indices))
	dataset_numerical_binary_column_index_file.close()
	dataset_lims_file = open("../adult/ad_normalization_info" + suffix_filename + ".json", "w")
	dataset_lims_file.write(json.dumps(norm_info))
	dataset_lims_file.close()

def load_entire_adult(feats_to_remove = []):
	suffix_filename = "" if len(feats_to_remove) == 0 else "-"
	for c in feats_to_remove:
		suffix_filename += (c + "_")
	if len(suffix_filename) != 0:
		suffix_filename = suffix_filename[:-1]

	column_names = pd.read_csv("../adult/columns" + suffix_filename + ".csv", sep=',', header=None).to_numpy()
	dataset = pd.read_csv("../adult/dataset" + suffix_filename + ".csv", sep=',', header=None).to_numpy()
	dataset_lims_file = open("../adult/ad_normalization_info" + suffix_filename + ".json", "r")
	norm_info = json.load(dataset_lims_file)
	return dataset[:, 1:], dataset[:, 0], norm_info, column_names[0][1:]

def load_splitted_adult():
	column_names = pd.read_csv("../adult/columns.csv", sep=',', header=None).to_numpy()
	training_set = pd.read_csv("../adult/training-set.csv", sep=',', header=None).to_numpy()
	test_set = pd.read_csv("../adult/test-set.csv", sep=',', header=None).to_numpy()
	return training_set, test_set, column_names[0][1:]

def normalize(X):
    X = X - X.min(axis=0)
    X = X/X.max(axis=0)
    X = np.nan_to_num(X)
    return X

def str_to_json(ensamble, dim, lab=[-1,1]):
    n_trees = len(ensamble)
    i = 0
    str_tree = "\t" + str(dim) + ",\n\t["
    for l in lab:
        str_tree += str(l) + ", "
    
    str_tree = str_tree[:-2] + "],\n\t[\n"
    
    for dict_tree in ensamble:
        app = np.array(dict_tree['feature'])
        app+=1
        app[app<0] = 0
        dict_tree['feature'] = app
        
        app = np.array(dict_tree['children_left'])
        app+=1
        app[app<0] = 0
        dict_tree['children_left'] = app
        
        app = np.array(dict_tree['children_right'])
        app+=1
        app[app<0] = 0
        dict_tree['children_right'] = app
        
        str_f = '\t\t\t[' + ', '.join([str(v) for v in dict_tree['feature']]) + '],\n'
        str_t = '\t\t\t[' + ', '.join([str(v) for v in dict_tree['threshold']]) + '],\n'
        str_l = '\t\t\t[' + ', '.join([str(v) for v in dict_tree['children_left']]) + '],\n'
        str_r = '\t\t\t[' + ', '.join([str(v) for v in dict_tree['children_right']]) + '],\n'
        str_p = '\t\t\t[' + ', '.join([str(v) for v in dict_tree['prediction']]) + ']\n'
        
        i+=1
        if i < n_trees:
            str_tree += '\t\t[\n' + str_f + str_t + str_l + str_r + str_p + '\t\t],\n'
        else:
            str_tree += '\t\t[\n' + str_f + str_t + str_l + str_r + str_p + '\t\t]'
        
    str_ensamble = '[\n' + str_tree + '\n\t]\n]'
    #print(str_ensamble)
    return str_ensamble

def _rec_(tree, dict_tree, idx_n):
    #print(tree)
    if 'split_index' in tree:
        dict_tree['feature'].append(tree['split_feature'])
        dict_tree['threshold'].append(tree['threshold'])
        dict_tree['children_left'].append(idx_n+1)
        dict_tree['children_right'].append(idx_n+2)
        dict_tree['prediction'].append(0)
        
        idx_n = _rec_(tree['left_child'], dict_tree, idx_n + 2)
        idx_n = _rec_(tree['right_child'], dict_tree, idx_n)
        return idx_n;
    else:
        dict_tree['feature'].append(-2)
        dict_tree['threshold'].append(-2)
        dict_tree['children_left'].append(-2)
        dict_tree['children_right'].append(-2)
        dict_tree['prediction'].append(tree['leaf_value'])
        return idx_n + 1

# prende in unput una lista di alberi
def sklearn_to_json(trees, dim, filename, mode=0):
    ensamble = []
    classes = np.asarray([-1, 1])
    n_trees = len(trees)
    i = 0
    str_tree = ""
    for tr in trees:
        dict_tree = {}
        n_nodes = len(tr.tree_.value)
        n_labels = len(tr.tree_.value[0][0])   
        value = tr.tree_.value.reshape(n_nodes, n_labels)
        prediction = classes[np.argmax(value, axis=1)]
        
        dict_tree['feature'] = tr.tree_.feature
        dict_tree['threshold'] = tr.tree_.threshold
        dict_tree['children_left'] = tr.tree_.children_left
        dict_tree['children_right'] = tr.tree_.children_right
        dict_tree['prediction'] = prediction
        ensamble.append(dict_tree)

    if mode == 0:
        with open(filename, 'w') as f:
            f.write(str_to_json(ensamble, dim))
    else:
        print(str_to_json(ensamble, dim))

# saves dataset in json format
def dataset_json(data, labels, filename):
    lbs = labels.tolist()
    dt = data.tolist()
    l = [lbs, dt]
    with open(filename, 'w') as f:
        f.write(json.dumps(l))

def rec_visit(node, dict_tree):
    input_id = len(dict_tree["feature"])

    if node!=None:
        print("FILL TREE ####", input_id, "#############")
        if node.is_leaf():
            dict_tree["feature"].append(-1)
            dict_tree["threshold"].append(-2)
            p = node.prediction
            p = -1 if p == 0 else 1
            dict_tree["prediction"].append(p)

            return input_id, [-1], [-1]
        else:
            
            dict_tree["feature"].append(node.best_split_feature_id)
            dict_tree["threshold"].append(node.best_split_feature_value)
            dict_tree["prediction"].append(0)

            node_id_l, l1l, l1r = rec_visit(node.left, dict_tree)
            node_id_r, l2l, l2r = rec_visit(node.right, dict_tree)

            return input_id, [node_id_l] + l1l + l2l, [node_id_r] + l1r + l2r

# prende in unput una lista di alberi
def sklearn_to_dict(trees, dim, filename):
    ensamble = []
    classes = np.asarray([-1, 1])

    for tr in trees:
        dict_tree = {}
        n_nodes = len(tr.tree_.value)
        n_labels = len(tr.tree_.value[0][0])   
        value = np.asarray(tr.tree_.value.reshape(n_nodes, n_labels), dtype=int)
        prediction = classes[np.argmax(value, axis=1)]

        dict_tree['feature'] = tr.tree_.feature
        dict_tree['threshold'] = tr.tree_.threshold
        dict_tree['children_left'] = tr.tree_.children_left
        dict_tree['children_right'] = tr.tree_.children_right
        dict_tree['prediction'] = prediction
        dict_tree['value'] = value
        ensamble.append(dict_tree)

    with open(filename, 'w') as f:
        f.write(foo(ensamble, dim))