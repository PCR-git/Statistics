"""
 Stat 202A 2019 Fall - Homework 07
 Author: Peter Racioppo
 Date: 11/17/2019

 INSTRUCTIONS: Please fill in the corresponding function. Do not change function names, 
 function inputs or outputs. Do not write anything outside the function.
"""

import numpy as np
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import KFold, GridSearchCV
from matplotlib import pyplot as plt
import xgboost as xgb

cancer = load_breast_cancer()
X = cancer.data
y = cancer.target

#############################################################################################################
# TODO (1) Perform 5-fold validation for cancer data.
# You may consider use KFold function in sklearn.model_selection
# Print the mean and std of the 5-fold validation accuracy
#############################################################################################################

def XGB(X,y,max_depth,min_child_weight):
    from sklearn.model_selection import KFold
#     from sklearn.tree import DecisionTreeClassifier
    kf = KFold(5,True,1)
    kf.get_n_splits(X)
    for train_index, test_index in kf.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

#     clf = DecisionTreeClassifier(criterion="entropy", max_depth=max_depth, min_weight_fraction_leaf = min_child_weight)
    clf = xgb.XGBClassifier(max_depth=max_depth, min_weight_child_weight = min_child_weight)
    from sklearn.model_selection import cross_val_score
    all_accuracies = cross_val_score(estimator=clf, X=X_train, y=y_train, cv=5)
    return all_accuracies.mean(), all_accuracies.std()

max_depth = 3
min_child_weight = 0.1
mean, std = XGB(X,y,max_depth,min_child_weight)
print("PART 1:")
print("mean =",round(mean,3))
print("std =",round(std,3))

#############################################################################################################
# TODO (2) Perform Grid Search for parameter max_depth in [3,5,7] and min_child_weight in [0.1, 1, 5]
# For each combination use 5-Fold validation
# You may consider use GridSearchCV function in sklearn.modelselection
# Print the grid search mean test score for each parameter combination (use cv_results_)
#############################################################################################################

from sklearn.model_selection import GroupKFold
# from sklearn.tree import DecisionTreeClassifier

def GridXGB(X,y,max_depth,min_child_weight):
    param_grid = dict(max_depth=max_depth, min_child_weight=min_child_weight)
#     clf = DecisionTreeClassifier(criterion="entropy", max_depth=max_depth, min_weight_fraction_leaf = min_child_weight)
    clf = xgb.XGBClassifier(max_depth=max_depth, min_child_weight = min_child_weight) 

    kf = KFold(5,True,1)
    kf.get_n_splits(X)
    grid_search = GridSearchCV(clf, param_grid, cv=kf)
    grid_result = grid_search.fit(X, y)
    result_dict = grid_result.cv_results_
#     print(result_dict)
    mean = result_dict.get('mean_test_score')
    std = result_dict.get('std_test_score')
    y_pred = grid_search.predict(X)
    bst = grid_result.best_estimator_
    return mean,std,y_pred,bst

max_depth = [3,5,7]
min_child_weight = [0.1,1,5]
mean, std, y_pred, bst = GridXGB(X,y,max_depth,min_child_weight)
print("PART 2:")
print("mean =",np.round(mean,3))
print("std =",np.round(std,3))

opt_max_depth = max_depth[int(np.floor(np.argmax(mean)/3))]
opt_min_child_weight = min_child_weight[np.argmax(mean)%3]

opt_max_depth = bst.max_depth
opt_min_child_weight = bst.min_child_weight

print("optimal max depth =",opt_max_depth)
print("optimal min child weight =",opt_min_child_weight)

width = np.shape(X[1])[0]
i = np.random.randint(width)
j = np.random.randint(width)
print('dim1=',i)
print('dim2=',j)
plt.scatter(X[:,i],X[:,j],y,color='blue',marker='o',linewidth=4)
plt.scatter(X[:,i],X[:,j],y_pred,color='red',marker='x',linewidth=4)
plt.grid()
plt.show()

#############################################################################################################
# TODO (3) Plot the feature importance of the best model
# You may fit a new xgboost model using all the data and then plot the importance using xgb.plot_importance()
#############################################################################################################

def XGB_importances(X,y,max_depth,min_child_weight):
    clf = xgb.XGBClassifier(max_depth=max_depth, min_weight_fraction_leaf = min_child_weight)
    clf_r = clf.fit(X,y)
    importances = clf_r.feature_importances_
    xgb.plot_importance(clf_r)
    plt.show()
    return importances

max_depth = opt_max_depth
min_child_weight = opt_min_child_weight
importances = XGB_importances(X,y,max_depth,min_child_weight)
# print(importances)
