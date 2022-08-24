from __future__ import absolute_import, division, print_function, unicode_literals

# standard python
import numpy as np
from sklearn.metrics import roc_auc_score, confusion_matrix, cohen_kappa_score
import pandas as pd

# plotting, especially for jupyter notebooks
from matplotlib import pyplot as plt
import matplotlib as mpl
import seaborn as sns; sns.set()


def display_conf_matrix(y_pred, y_true, Gmodelname=None, Display=False,  Save_path=None):
    """Create confusion matrix and display model fit metrics.
    
    :parameter y_pred: predicted list of abnormal (1) and normal (0) values.
            :type: np.array
    :parameter y_true: impirical list of abnormal (1) and normal (0) values.
            :type: np.array
    :parameter Gmodelname: image output file name containing the type of model.
            :type: string
    :parameter Display: Display confusion matrix.
            :type: Boolean
    :parameter Save_path: Path were the confusion matrix should be saved to file (only when Display=True).
            :type: str
            
    :returns: NN_kappa_score, NN_auroc, SE, SP, PPV
    :rtype: int
        NN_kappa_score -- Kappa Statistic.
        NN_auroc -- Area under the Receiver Operating Curve.
        SE -- Sensitivity
        SP --  Specificity
        PPV -- Positive Predictive Value
    """
    
    y_pred = y_pred
    y_true = y_true
    conf_matrix = confusion_matrix(y_true, y_pred)
    tp = conf_matrix[1,1]
    tn = conf_matrix[0,0]
    conf_matrix[1,1] = tn
    conf_matrix[0,0] = tp
    # print(conf_matrix)
    NN_kappa_score = round(cohen_kappa_score(y_true, y_pred),3)
    NN_auroc = round(roc_auc_score(y_true, y_pred),4)
    

    mpl.style.use('seaborn')
    cm = pd.DataFrame(conf_matrix, 
                        index = ['Abnormal', 'Normal'],
                        columns = ['Abnormal', 'Normal'])    
    cm_rsum = np.sum(conf_matrix, axis=1, keepdims=True)
    cm_csum = np.sum(conf_matrix, axis=0, keepdims=True)

    cm_perc = cm / cm_csum.astype(float) * 100

    cm_csum = pd.DataFrame(cm_csum, 
                               index = ['True'],
                               columns = ['Abnormal', 'Normal'])
    cm_rsum = pd.DataFrame(cm_rsum, 
                               index = ['Abnormal', 'Normal'],
                               columns = ['Predicted'])

    cm = pd.concat([cm, cm_csum], axis=0)
    cm = pd.concat([cm, cm_rsum], axis=1, sort=False)
    cm.at['True','Predicted'] = cm.iloc[0, 0] + cm.iloc[1, 1]
    annot = np.empty_like(cm).astype(str)
    cm_max = np.sum(conf_matrix)
    
    cmp = cm_perc.copy()
    cm_perc = np.zeros([cmp.shape[0]+1,cmp.shape[1]+1])
    cm_perc[0:2,0:2] = cmp  
    cm_perc = pd.DataFrame(cm_perc)

    nrows, ncols = cm.shape
    for i in range(nrows):
        for j in range(ncols):
            c = cm.iloc[i, j]
            if i < 2 and j < 2:
                p = cm_perc.iloc[i, j]     
                if c <= 0:
                    annot[i, j] = '{}\n0.0%'.format(c)
                else:
                    annot[i, j] = '{}\n{:2.1f}%'.format(c,p)
                cm_perc.iloc[i, j] = 0
            elif i == 2 and j == 2:
                p = c / cm_max.astype(float) * 100
                annot[i, j] = 'ACC\n{:2.1f}%'.format(p)
                cm_perc.iloc[i, j] = p
            elif i == 2:
                p = cm.iloc[j,j] / c * 100
                if j == 0:
                    annot[i, j] = 'SE\n{:2.1f}%'.format(p)
                    SE = p
                    cm_perc.iloc[i, j] = p
                if j == 1:
                    annot[i, j] = 'SP\n{:2.1f}%'.format(p)
                    SP = p
                    cm_perc.iloc[i, j] = p
            elif j == 2:
                p = cm.iloc[i,i] / c * 100
                if i == 0:
                    annot[i, j] = 'PPV\n{:2.1f}%'.format(p)
                    PPV = p
                    cm_perc.iloc[i, j] = p
                if i == 1:
                    annot[i, j] = 'NPV\n{:2.1f}%'.format(p)
                    cm_perc.iloc[i, j] = p
    cm_perc.index.name = 'Predicted'
    cm_perc.columns.name = 'True'
    
    if Display:
        fig = plt.figure()
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_aspect(1)
        cmap = sns.cubehelix_palette(light=0.85, dark=0.15, as_cmap=True)
#         res = sns.heatmap(cm_perc, annot=annot, vmin=0.0, vmax=cm_max, fmt='', ax=ax, cmap=cmap)        
        res = sns.heatmap(cm_perc, annot=annot, vmin=0.0, vmax=100, fmt='', ax=ax, cmap=cmap)
        plt.yticks([0.5,1.5], ['Abnormal', 'Normal'],va='center')
        plt.xticks([0.5,1.5], ['Abnormal', 'Normal'],va='center')

        if Gmodelname is not None:
            if Save_path is None:
                imageOut = '/home2/ajgreen4/Read-Across_w_GAN/imageOut/'
            else:
                imageOut = Save_path
            ## Is this a regression generator?
            if (Gmodelname.find('Go-ZT')!=-1):
                plt.title('Go-ZT Confusion Matrix')
                cm_file = imageOut+Gmodelname+'-confusion_matrix.png'
            elif (Gmodelname.find('GAN-ZT')!=-1):
                plt.title('GAN-ZT Confusion Matrix')
                cm_file = imageOut+Gmodelname+'-confusion_matrix.png' 
            else:
                plt.title(Gmodelname+' \nConfusion Matrix')
                cm_file = imageOut+Gmodelname+'-confusion_matrix.png'
                
            if Save_path is not None:
                plt.savefig(cm_file, dpi=600, bbox_inches='tight')  
        plt.pause(0.5)
        print('\n\nKappa: ',NN_kappa_score)       
        print('Auroc: ',NN_auroc)
        
    return NN_kappa_score, NN_auroc, SE, SP, PPV

def display_random_conf_matrix(y_pred, y_true, Display=0, Save_path=None):
    """Create confusion matrix and display model fit metrics.
    
    :parameter y_pred: predicted list of abnormal (1) and normal (0) values.
            :type: np.array
    :parameter y_true: impirical list of abnormal (1) and normal (0) values.
            :type: np.array
    :parameter Display: Display confusion matrix.
            :type: Boolean
    :parameter Save: Path were the confusion matrix should be saved to file.
            :type: str
    """
    
    y_pred = y_pred[:,-1]
    y_true = y_true[:,-1]
    conf_matrix = confusion_matrix(y_true, y_pred)
    conf_matrix = np.round(conf_matrix/Display,0)  
    tp = conf_matrix[1,1]
    tn = conf_matrix[0,0]
    conf_matrix[1,1] = tn
    conf_matrix[0,0] = tp
    # print(conf_matrix)
    NN_kappa_score = round(cohen_kappa_score(y_true, y_pred),3)
    NN_auroc = round(roc_auc_score(y_true, y_pred),4)
    

    mpl.style.use('seaborn')
    cm = pd.DataFrame(conf_matrix, 
                        index = ['Active', 'Inactive'],
                        columns = ['Active', 'Inactive'])    
    cm_rsum = np.sum(conf_matrix, axis=1, keepdims=True)
    cm_csum = np.sum(conf_matrix, axis=0, keepdims=True)

    cm_perc = cm / cm_csum.astype(float) * 100

    cm_csum = pd.DataFrame(cm_csum, 
                               index = ['Real Toxicity Data'],
                               columns = ['Active', 'Inactive'])
    cm_rsum = pd.DataFrame(cm_rsum, 
                               index = ['Active', 'Inactive'],
                               columns = ['Generated Toxicity Data'])

    cm = pd.concat([cm, cm_csum], axis=0)
    cm = pd.concat([cm, cm_rsum], axis=1, sort=False)
    cm.at['Real Toxicity Data','Generated Toxicity Data'] = cm.iloc[0, 0] + cm.iloc[1, 1]
    annot = np.empty_like(cm).astype(str)
    cm_max = np.sum(conf_matrix)    

    nrows, ncols = cm.shape
    for i in range(nrows):
        for j in range(ncols):
            c = cm.iloc[i, j]
            if i < 2 and j < 2:
                p = cm_perc.iloc[i, j]
                if c <= 0:
                    annot[i, j] = '{}\n0.0%'.format(c)
                else:
                    annot[i, j] = '{}\n{:2.1f}%'.format(c,p)
            elif i == 2 and j == 2:
                p = c / cm_max.astype(float) * 100
                annot[i, j] = 'ACC\n{:2.1f}%'.format(p)
            elif i == 2:
                p = cm.iloc[j,j] / c * 100
                if j == 0:
                    annot[i, j] = 'SE\n{:2.1f}%'.format(p)
                    SE = p
                if j == 1:
                    annot[i, j] = 'SP\n{:2.1f}%'.format(p)
                    SP = p
            elif j == 2:
                p = cm.iloc[i,i] / c * 100
                if i == 0:
                    annot[i, j] = 'PPV\n{:2.1f}%'.format(p)
                    PPV = p
                if i == 1:
                    annot[i, j] = 'FOR\n{:2.1f}%'.format(p)
    cm.index.name = 'Generated Toxicity Data'
    cm.columns.name = 'Real Toxicity Data'
    
    if Display:
        fig = plt.figure()
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_aspect(1)
        cmap = sns.cubehelix_palette(light=0.85, dark=0.15, as_cmap=True)
        res = sns.heatmap(cm, annot=annot, vmin=0.0, vmax=cm_max, fmt='', ax=ax, cmap=cmap)
        plt.yticks([0.5,1.5], ['Active', 'Inactive'],va='center')
        plt.xticks([0.5,1.5], ['Active', 'Inactive'],va='center')
        plt.title('Random Confusion Matrix')
        plt.pause(0.5)    
        print('\n\nKappa: ',NN_kappa_score)
        print('Auroc: ',NN_auroc)   
        
    if Save_path is not None:
        imageOut = Save_path
        ## Is this a regression generator?
        cm_file = imageOut+'random_confusion_matrix.png'              
        plt.savefig(cm_file, dpi=600, bbox_inches='tight' )   
    
#     return NN_kappa_score, NN_auroc, SE, SP, PPV