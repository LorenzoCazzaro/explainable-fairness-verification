from misc import *
import json
from sklearn.model_selection import GridSearchCV
from joblib import dump
import sys

seed = int(sys.argv[3]) #7
np.random.seed(seed)

dataset_name = sys.argv[1]
also_type = int(sys.argv[2])

if also_type == 3:
    feats_to_remove = sys.argv[4].strip().split("-")
    print(feats_to_remove)

if dataset_name == "ad":
    X, y, norm_info, column_names = load_entire_adult(feats_to_remove) if also_type == 3 else load_entire_adult()
    name = 'adult'

path = '../' + name

n_ist, n_feat = X.shape
dim = n_feat
print('original', X.shape)

X, idx = np.unique(X, axis=0, return_index=True)
y = y[idx]
print(y)
print('unique', X.shape)
print(np.unique(y))
print("Class distr (positive over total)", np.sum(y==1) / len(y))

#dataset are pre-processed before executing this
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=seed, stratify=y) #stratify added now

np.random.seed(seed)

if also_type == 0: 
    path += '/models/'
    n_estimators = int(sys.argv[4])
    max_depth = int(sys.argv[5])
    filename = 'rf_' + dataset_name + '_' + str(n_estimators) + '_' + str(max_depth) + "_" + str(seed)
    clf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed, max_depth=max_depth)
    clf.fit(X_train, y_train)
    print("TREE SCORE: ", clf.score(X_test, y_test))
    sklearn_to_json(clf.estimators_, dim, path + filename + ".json")

if also_type == 1:
    print(X_test.shape)
    print(path + '/test_' + name + '_orig.json')
    dataset_json(X_test, y_test, path + '/test_' + name + '_orig.json')
    dataset_json(X_train, y_train, path + '/train_' + name + '_orig.json')