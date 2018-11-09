
# coding: utf-8

# In[24]:


import sys, os
from sklearn.externals import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn import model_selection
from scipy.stats import randint as sp_randint
from sklearn import svm
from time import time
from operator import itemgetter
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn import preprocessing
from sklearn.model_selection import train_test_split


# In[3]:


trainingSetDir = "/media/kevin/TOSHIBA_EXT/TTT_RecombinationGenomeScans/ml_project/feature_vecs_h12"
classifierPickleFileName = "addH12.p"
statsToUse = "all"
print(os.listdir(trainingSetDir))


# In[6]:



classList = []
trainingData = []
labelToClassName = {}
headerH = {}

for trainingSetFileName in os.listdir(trainingSetDir):
    classList.append(trainingSetFileName.split(".fvec")[0])
    trainingSetFile = open(trainingSetDir + "/" + trainingSetFileName)
    currTrainingData = trainingSetFile.readlines()
    trainingSetFile.close()

    trainingData += currTrainingData[1:]#append all training data from the current set (minus the header)

    currLabelH = {}
    for example in currTrainingData[1:]:
        currLabelH[example.split("\t")[0]] = 1
    assert len(currLabelH) == 1, "Length: %d  label: %s" %(len(currLabelH), currLabelH)
    labelToClassName[list(currLabelH.keys())[0]] = trainingSetFileName.split(".fvec")[0]
    
    header = currTrainingData[0].strip().split("\t")
    headerH[currTrainingData[0].strip()] = 1
    assert header[0] == "classLabel"
    statIndices = []
    if "all" in statsToUse:
        statIndices = range(1, len(header))
    else:
        for i in range(1, len(header)):
            if header[i] in statsToUse or header[i].split("_win")[0] in statsToUse:
                statIndices.append(i)
assert len(headerH) == 1

sys.stderr.write("using these features: %s (indices: %s)\n" %(str(statsToUse), str(statIndices)))
XH = {}
for i in range(len(trainingData)):
    trainingData[i] = trainingData[i].strip().split("\t")
    currVector = []
    if not "nan" in trainingData[i]:
        for j in statIndices:
            try:
                currVector.append(float(trainingData[i][j]))
            except:
                print("Invalid data at coordinates: " + str(i) + "," + str(j))
                surr = trainingData[i-1 : i + 5]
                for site in surr: print(site)
                print(trainingData[i])
                print(header[j])
        assert len(currVector) == len(statIndices),         "length of current vector: %s doesn't match length of stat indices: %s" %(len(currVector), len(statIndices))
        if trainingData[i][0] not in XH:
            XH[trainingData[i][0]] = []
        XH[trainingData[i][0]].append(currVector)


# In[8]:


X = []
y = []
for classLabel in sorted(XH.keys()):
    print('{:24} : {:>6}'.format(classLabel, str(len(XH[classLabel]))))
    random.shuffle(XH[classLabel])
    for i in range(1000):
        try:
            currVector = XH[classLabel][i]
        except IndexError:
            break
        X.append(currVector)
        y.append(classLabel)
        
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
sys.stderr.write("Training set size after split: %s\n" %(len(y_train)))
sys.stderr.write("Testing set size: %s\n" %len(y_test))


# In[9]:


print(XH.keys())
labels = XH.keys()


# In[47]:


# Utility function to report best scores
def report(grid_scores, n_top=3):
    scores = pd.DataFrame(grid_scores)
    top_scores = scores.sort_values(by = ['mean_test_score'], ascending = False).iloc[:n_top]
    
    for i in range(len(top_scores)):
        score = top_scores.iloc[i]
        print("Model with rank: {0}".format(i + 1))
        print("Mean test score: {0:.3f} (std: {1:.3f})".format(
              score.mean_test_score,
              score.std_test_score))
        print("Parameters: {0}".format(score.params))
        print("")


# In[35]:


sys.stderr.write("Checking accuracy when distinguishing among all %s classes\n" %(len(XH.keys())))

maxMaxFeatures = len(X[0])
param_grid_forest = {"max_depth": [3, 10, None],
              "max_features": [1, 3, int(maxMaxFeatures**0.5), maxMaxFeatures],
              "min_samples_split": [2, 3, 10],
              "min_samples_leaf": [1, 3, 10],
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}

clf, mlType, paramGrid = ExtraTreesClassifier(n_estimators=100), "extraTreesClassifier", param_grid_forest

heatmap = []
sys.stderr.write("Training %s\n" %(mlType))
grid_search = GridSearchCV(clf, param_grid=param_grid_forest, cv=10,n_jobs=-1, return_train_score = False)
start = time()
grid_search.fit(X_train, y_train)
sys.stderr.write("GridSearchCV took %.2f seconds for %d candidate parameter settings.\n"
      % (time() - start, len(grid_search.cv_results_)))
print("Results for %s" %(mlType))
report(grid_search.cv_results_)
joblib.dump((X_test, y_test, grid_search, list(XH.keys())), classifierPickleFileName)


# In[48]:


report(grid_search.cv_results_)


# In[ ]:


plotReportForComparison("QTNdistinct.p")

