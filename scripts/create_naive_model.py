#!/usr/bin/python3
import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score

from sklearn.naive_bayes import GaussianNB
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from itertools import cycle
import matplotlib.pyplot as plt

# Training dataframe
tdf = pd.read_csv('training_set.csv')
mito_df = pd.read_csv('testing_set.csv')

# Training dataframe
tdf = pd.read_csv('training_set.csv')
mito_df = pd.read_csv('testing_set.csv')

subset = ['alphaproteobacteria', 'betaproteobacteria', 'deltaproteobacteria', 'epsilonproteobacteria']

# Make all samples equal in size
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
gnb = GaussianNB().fit(X_train, y_train.values.ravel())
gnb_predictions = gnb.predict(X_test)

# creating a confusion matrix
cm = confusion_matrix(y_test, gnb_predictions)
cr = classification_report(y_test,gnb_predictions)
print(cm)
print(cr)
predictions_mito = gnb.predict(X_mito)
print('mitochondrial predictions: ',predictions_mito)

# ROC/AUC calculations
y = tdf.select_dtypes(include=[object])
y = label_binarize(y, classes=['alphaproteobacteria', 'betaproteobacteria', 'deltaproteobacteria', 'epsilonproteobacteria'])
classifier = OneVsRestClassifier(GaussianNB())

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20)

y_score = classifier.fit(X_train, y_train).predict_proba(X_test)
n_classes = y.shape[1]

# Create ROC plot
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

print('AUC', roc_auc_score(y_test, y_score, multi_class='ovr'))
plt.plot([0, 1], [0, 1], 'k--', lw=1)
plt.xlim([-0.05, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Naive Bayes model (AUC=0.8383)')
plt.legend(loc="lower right")
plt.show()
