import glob
from scipy import interp
import  numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve, auc, RocCurveDisplay
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelBinarizer 
import os
os.chdir('/data/dk/kirc/fig')
files = glob.glob('/data/dk/kirc/CLAM/eval_results/EVAL_kirc70_val/fold*')
lb = LabelBinarizer()
fprs, tprs = [], []
i = 0
for file in files:
# Parse the table pd.read_csv(file)
    i += 1
    y_true = pd.read_csv(file)['Y']
    y_pred = pd.read_csv(file)['Y_hat']
    y_pred_probs = pd.read_csv(file)[['p_0','p_1','p_2', 'p_3']].values
    ytruehot = lb.fit_transform(y_true) 
    ytruehot
    fpr, tpr, _= roc_curve(ytruehot.ravel(), y_pred_probs.ravel())
    print(len(fpr))
    fprs.append(fpr)
    tprs.append(tpr)
    roc_auc = auc(fpr,tpr)
    plt.plot(fpr, tpr, label='Fold {}'.format(i) + ' (AUC = %0.1f)' % roc_auc, color='grey')



#the fprs and tprs which roc_curve returned arrays of different length for each fold , 
#To handle this situation, one can calculate the average ROC curve by interpolating 
#the different-length arrays to a common set of thresholds. 
mean_fpr = np.linspace(0, 1, 100)
mean_tpr = np.zeros_like(mean_fpr)


for i in range(len(fprs)):
    fpr = fprs[i]
    tpr = tprs[i]
    mean_tpr += interp(mean_fpr, fpr, tpr)

mean_tpr /= len(fprs)
mean_tpr[0] = 0.0
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)

plt.plot(mean_fpr, mean_tpr, label=f'Mean ROC (AUC = {mean_auc:.2f})', color='red')
plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc='lower right')
plt.show()
plt.savefig('./auc.png')









##plot confusion matrix


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from matplotlib import ticker

for i in range(0, 10):
    print(i)
    file_name = f'/data/dk/kirc/CLAM/eval_results/EVAL_kirc70_val/fold_{i}.csv'
    files = pd.read_csv(file_name)
    true_labels = files['Y']
    predicted_labels = files['Y_hat']

    # Calculate the confusion matrix
    cm = confusion_matrix(true_labels, predicted_labels)

    # Create a list of class labels
    class_labels = ['Class 0', 'Class 1', 'Class 2', 'Class 3']

    # Create a heatmap of the confusion matrix
    ax = sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=class_labels, yticklabels=class_labels, cbar=True)

    # Add labels, title, and axis ticks
    plt.xlabel('Predicted Labels')
    plt.ylabel('True Labels')
    plt.title('Confusion Matrix')

    # Set the x-axis and y-axis ticks as integers
    plt.xticks(np.arange(len(class_labels)), np.arange(len(class_labels)))
    plt.yticks(np.arange(len(class_labels)), np.arange(len(class_labels)))

    # Set the color bar ticks as integers
    cbar = ax.collections[0].colorbar
    cbar.set_ticks(np.arange(np.max(cm) + 1))
    cbar.set_ticklabels(np.arange(np.max(cm) + 1))

    # Remove the legend
    plt.legend().remove()

    # Display the plot
    plt.show()









##### plot precision recall curve
#################################
import pandas as pd
from sklearn.metrics import precision_recall_curve
from sklearn.preprocessing import LabelBinarizer 
from sklearn.metrics import average_precision_score
data = pd.read_csv('/data/dk/kirc/CLAM/eval_results/EVAL_kirc70_val/fold_1.csv')

Y_test = data['Y']
lb = LabelBinarizer()
yhot = lb.fit_transform(Y_test)
y_score = data[['p_0','p_1','p_2', 'p_3']].values
precision = dict()
recall = dict()
average_precision = dict()
for i in range(3):
    precision[i], recall[i], _ = precision_recall_curve(yhot[:, i],
                                                        y_score[:, i])

    average_precision[i] = average_precision_score(yhot[:, i], 
                                                   y_score[:, i])


#"macro-average": quantifying score on all classes jointly
precision["macro"], recall["macro"], _ = precision_recall_curve(yhot.ravel(), y_score.ravel())

average_precision["macro"] = average_precision_score(yhot, y_score,                                                     average="macro")

print('Average precision score, macro-averaged over all classes: {0:0.2f}'.format(average_precision["macro"]))
plt.figure()
plt.step(recall['macro'], precision['macro'], where='post')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title('Average precision score, macro-averaged over all classes: AP={0:0.3f}'.format(average_precision["macro"]))
plt.show()
