r"""°°°
pycaret based mechine learning , integrated multi algorithm to predict morphology features extracted by cellprofiler<br>
kernel:'pycaret'
°°°"""
# |%%--%%| <nXVe1j5o6P|bVxiQFUk5s>

import os 
import pandas as pd

# |%%--%%| <bVxiQFUk5s|2CQXQNBJq5>

ucec = pd.read_csv('/data/dk/ucec/cf_out/Image.csv',delimiter=',')
ucec.select_dtypes(include='object')
columns_to_keep = [col for col in ucec.columns if col == 'FileName_myimage' or ucec[col].dtype != 'object']
ucec2 = ucec[columns_to_keep]
ucec2.loc[:,'case']  = ucec2['FileName_myimage'].str.split('_tile').str[0]
ucec3 = ucec2.groupby('case').mean(numeric_only = True)
refer_data = pd.read_csv('/data/dk/ucec/DataSheet_2.csv')
ucec3 = ucec3[refer_data.drop(columns=['case']).columns]
ucec3 = ucec3.reset_index()
ucec3['case_submitter_id'] = ucec3['case'].str.extract(r'^(.*?-.*?-.*?)-',expand = False)
ucec_meta = pd.read_csv('/data/dk/ucec/clinical.tsv', delimiter='\t')[['case_submitter_id','figo_stage']]
data_ucec = ucec3.merge(ucec_meta,on='case_submitter_id',how='inner')
replacements = {
    r'Stage IV.*': "Stage_IV",
    r'Stage III.*': "Stage_III",
    r'Stage II.*': "Stage_II",
    r'Stage I.*': "Stage_I"
}
for pattern ,replace in replacements.items():
    data_ucec['figo_stage'] = data_ucec['figo_stage'].str.replace(pattern,replace)
data_ucec.drop(columns=['case','case_submitter_id'],inplace=True)

# |%%--%%| <2CQXQNBJq5|9jTdIhYTfo>

tma = pd.read_csv('/data/dk/endometrium/cf_out/Image.csv',delimiter=',')
tma.select_dtypes(include='object')
columns_to_keep = [col for col in tma.columns if col == 'FileName_myimage' or tma[col].dtype != 'object']
filtered_df = tma[columns_to_keep]
filtered_df.loc[:,'case']  = filtered_df['FileName_myimage'].str.split('_tile').str[0]
filtered_df2 = filtered_df.groupby('case').mean(numeric_only = True)
# refer_data = pd.read_csv('/data/dk/ucec/cf_out/DataSheet_2.csv')
tma3 = filtered_df2[refer_data.drop(columns=['case']).columns]
endo_meta = pd.read_csv('/data/dk/ucec/endo_meta.csv', delimiter=',').iloc[:,[0,15]]
endo_meta.rename(columns={'序号': 'case_submitter_id', 'TNM': 'figo_stage'}, inplace=True)
tma3 = tma3.reset_index()
tma3['case_submitter_id'] = tma3['case'].str.extract('([A-Z])(\d+)').apply(lambda x: f"{x[0]}{x[1].zfill(2)}", axis=1)
data_tma3 = tma3.merge(endo_meta,on='case_submitter_id',how='inner')
replacements = {
    r'T1.*': "Stage_I",
    r'T4.*': "Stage_IV",
    r'T3.*': "Stage_III"
}
for pattern ,replace in replacements.items():
    data_tma3['figo_stage'] = data_tma3['figo_stage'].str.replace(pattern,replace)
data_tma3.drop(columns=['case','case_submitter_id'],inplace=True)

# |%%--%%| <9jTdIhYTfo|bHcakn8m1q>

data_tma3.drop(list(data_tma3.filter(regex = 'Threshold')), axis = 1, inplace = True)
data_ucec.drop(list(data_ucec.filter(regex = 'Threshold')), axis = 1, inplace = True)


# |%%--%%| <bHcakn8m1q|4Xrt7v8ZoJ>


import pycaret
from pycaret.classification import *


# |%%--%%| <4Xrt7v8ZoJ|21qVvyO0qO>

from pycaret.classification import ClassificationExperiment
exp = ClassificationExperiment()

# |%%--%%| <21qVvyO0qO|domK13gAMx>

# init setup on exp
setup(data_ucec, target = 'figo_stage', train_size = 0.8, normalize = False, 
            test_data = data_tma3,session_id = 123,index=False, 
            fix_imbalance = True)

# |%%--%%| <domK13gAMx|EkXayZMCKx>

best = compare_models(n_select=3)

# |%%--%%| <EkXayZMCKx|xYQaQiyOme>

dt = create_model('lr')
tuned_dt = tune_model(dt, optimize = 'AUC')

# |%%--%%| <xYQaQiyOme|FabGQks7GF>

import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from itertools import cycle
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer

# |%%--%%| <FabGQks7GF|bIYVP1JunB>

X_train = data_ucec.drop('figo_stage', axis=1).values
X_test = data_tma3.drop('figo_stage', axis=1).values
y_train = data_ucec['figo_stage'].values
y_test = data_tma3['figo_stage'].values


# |%%--%%| <bIYVP1JunB|zCzJ6tDZsc>

from sklearn.linear_model import LogisticRegression

y_score = tuned_dt.predict_proba(X_test)

# |%%--%%| <zCzJ6tDZsc|WieIw6fTLf>

from sklearn.preprocessing import LabelBinarizer

label_binarizer = LabelBinarizer().fit(y_train)
y_onehot_test = label_binarizer.transform(y_test)
y_onehot_test.shape  

# |%%--%%| <WieIw6fTLf|9evYMTRczB>

from itertools import cycle

fig, ax = plt.subplots(figsize=(6, 6))

plt.plot(
    fpr["micro"],
    tpr["micro"],
    label=f"micro-average ROC curve (AUC = {roc_auc['micro']:.2f})",
    color="deeppink",
    linestyle=":",
    linewidth=4,
)

plt.plot(
    fpr["macro"],
    tpr["macro"],
    label=f"macro-average ROC curve (AUC = {roc_auc['macro']:.2f})",
    color="navy",
    linestyle=":",
    linewidth=4,
)

colors = cycle(["aqua", "darkorange", "cornflowerblue"])
for class_id, color in zip(range(n_classes), colors):
    RocCurveDisplay.from_predictions(
        y_onehot_test[:, class_id],
        y_score[:, class_id],
        name=f"ROC curve for {target_names[class_id]}",
        color=color,
        ax=ax
    )

plt.axis("square")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Extension of Receiver Operating Characteristic\nto One-vs-Rest multiclass")
plt.legend()
plt.show()

# |%%--%%| <9evYMTRczB|rjAB8mnz4e>

ktarget_names = data_ucec['figo_stage'].unique()
target_names


# |%%--%%| <rjAB8mnz4e|Zl5jC1BnXw>


