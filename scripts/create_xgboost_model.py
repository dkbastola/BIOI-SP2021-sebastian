#!/usr/bin/python3
import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from itertools import cycle
import matplotlib.pyplot as plt


import matplotlib.pyplot as plt

import xgboost as xgb


# Training dataframe
tdf = pd.read_csv('training_set.csv')
mito_df = pd.read_csv('testing_set.csv')

subset = ['alphaproteobacteria', 'betaproteobacteria', 'deltaproteobacteria', 'epsilonproteobacteria']


tdf = tdf[tdf['label'].isin(subset)]
tdf = tdf.groupby('label')
tdf = pd.DataFrame(tdf.apply(lambda x: x.sample(tdf.size().min()).reset_index(drop=True)))

print(tdf.label.value_counts())
del tdf['species']
del mito_df['species']

# Independent variables
X = tdf.iloc[:, 1:64]
X_mito = mito_df.iloc[:, 1:64]

# Dependent variables
y = tdf.select_dtypes(include=[object])
print(y.label.unique())

# Label encoder
le = preprocessing.LabelEncoder()
y = y.apply(le.fit_transform)
print(y.label.unique())

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.10)

xgb_model = xgb.XGBClassifier(objective="multi:softprob", random_state=42)
xgb_model.fit(X_train, y_train.values.ravel())

predictions = xgb_model.predict(X_test)
y_pred = xgb_model.predict(X_mito)
mito_df['predicted'] = y_pred
mito_df.to_csv('xgboost_predicted.csv')

print(confusion_matrix(y_test, predictions))
print(classification_report(y_test, predictions))

# Mitochondrial predictions
print('mitochondrial predictions: ', y_pred)
y_score = xgb_model.fit(X_train, y_train).predict_proba(X_test)

print('AUC ',roc_auc_score(y_test, y_score, multi_class='ovr'))

varimp= xgb_model.feature_importances_.argsort()

print('Variable importance: ', varimp)

