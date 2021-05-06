#!/usr/bin/python3
import numpy as np
import pandas as pd

from sklearn import preprocessing, metrics
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from itertools import cycle
import matplotlib.pyplot as plt


import matplotlib.pyplot as plt  

# Training dataframe
tdf = pd.read_csv('training_set.csv')
mito_df = pd.read_csv('testing_set.csv')
mito_df2 = pd.read_csv('testing_set.csv')

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

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20)

# Feature scaling
scaler = StandardScaler()
scaler.fit(X_train)

X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)

# Training
mlp = MLPClassifier(hidden_layer_sizes=(10, 10, 10), max_iter=1000, random_state=1)
mlp.fit(X_train, y_train.values.ravel())

# predictions
predictions = mlp.predict(X_test)
predictions_mito = mlp.predict(X_mito)
mito_df2['predicted'] = predictions_mito
mito_df2.to_csv('nnet_predicted.csv', index=False)



print(confusion_matrix(y_test,predictions))
print(classification_report(y_test,predictions))

print(predictions_mito)

# ROC curve
y_score = mlp.predict_proba(X_test)
n_classes = 3


y_score = mlp.predict_proba(X_test)
print(roc_auc_score(y_test, y_score, multi_class='ovr'))

# ROC

y = tdf.select_dtypes(include=[object])
y = label_binarize(y, classes=['alphaproteobacteria', 'betaproteobacteria', 'deltaproteobacteria', 'epsilonproteobacteria'])
classifier = OneVsRestClassifier(MLPClassifier(hidden_layer_sizes=(10, 10, 10), max_iter=1000, random_state=1))

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20)

y_score = classifier.fit(X_train, y_train).predict_proba(X_test)
n_classes = y.shape[1]

fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

colors = cycle(['blue', 'red', 'green', 'yellow'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=1,
             label='ROC curve of class {0} (area = {1:0.2f})'
             ''.format(subset[i], roc_auc[i]))

print(roc_auc_score(y_test, y_score, multi_class='ovr'))
plt.plot([0, 1], [0, 1], 'k--', lw=1)
plt.xlim([-0.05, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Artificial Neural Network model (AUC=0.9580)')
plt.legend(loc="lower right")
plt.show()

