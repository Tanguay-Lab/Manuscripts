################################################################################
######################### MISCELLANEOUS SCRIPTS ################################
################################################################################

"""
A variety of smaller scripts for data input, data wrangling and transformation,
data plotting and data analysis. Scripts were usually tested interactively in JupyterLab notebook
before being transferred here
"""

################################################################################
################################ DEPENDENCIES ##################################
################################################################################

import os, math, datetime, random, itertools
from collections import Counter
import numpy as np
import pandas as pd
import pickle
import matplotlib as mpl
import seaborn as sns; sns.set()
import scipy.stats

from matplotlib import pyplot as plt
from scipy.optimize import fmin
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import ks_2samp
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.metrics import roc_auc_score, confusion_matrix, cohen_kappa_score
from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LogisticRegression 
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import MinMaxScaler 

from tensorflow.keras.layers import Dense, Flatten, Reshape, Input, InputLayer, BatchNormalization, Dropout
from tensorflow.keras.models import Model, Sequential 
from tensorflow.keras import regularizers 


################################################################################
############################## DATA INPUT-OUTPUT################################
################################################################################

def save_obj(obj, path, name):
    
    """
    Saves Dictionary or model to pickle.
    ----------
    obj: [Dictionary] or [model].
    path: path to saving location.
    name: name of saved file. Format: /path/name.
    """

    # remove potential extensions
    if os.path.isfile(path+name):
        name = ".".join(name.split(".")[:-1])
    
    with open(path + name + '.pkl', 'wb') as file:
        pickle.dump(obj, file, pickle.HIGHEST_PROTOCOL)
    print('File saved: ', path, name, '.pkl')

################################################################################

def load_obj(path, name):
    
    """
    loads Dictionary or model from pickle.
    ----------
    path: path to load location.
    name: name of saved file. Format: /path/name.
    ----------
    returns [Dictionary] or [model]
    """
    
    # remove potential extensions
    if os.path.isfile(path+name):
        name = ".".join(name.split(".")[:-1])
    
    with open(path + name + '.pkl', 'rb') as file:
        return pickle.load(file)
    
################################################################################
################# DATA TRANSFORMATION AND FEATURE SELECTION ####################
################################################################################

def agg_by_zf(movement_df, begin, end, verbose=0):
    """Aggregate movement data by individual fish (plate, chemical, well, concentration).

    :parameter: movement_df: zebrafish behaviour data.
    :type: movement_df: DataFrame
    :parameter: begin:  Starting frame.
    :type begin: int
    :parameter: end:  Ending frame (exclude).
    :type end: int
    :parameter verbose: Print info or not
    :type verbose: Boolean
    
    :returns: Movement data aggregated by plate, chemical, well, and concentration.
    :rtype: DataFrame

    """
    
    df = movement_df.copy()
    
    # Subset data
    df = df.loc[(df['Trial time'] >= begin) & (df['Trial time'] < end)]
    
    # Aggregate movement data and reindex
    df1 = df.groupby(by=['PlateID', 'ChemID', 'Well', 'Conc'], 
                     sort=False)['Distance moved'].sum().reset_index()

    if verbose:
        print("Start time: ", begin, " Stop time: ", end, "\n")
        print('Total number of minutes: ', str(int((end/60))-int((begin/60))))
        print(df1.head())
    
    return df1

################################################################################

def agg_by_time(movement_df, begin, end, Nsec, agg_group, verbose=0):
    """Aggregate movement data by plate, chemical, well, concentration and second.

    :parameter: movement_df: zebrafish behaviour data.
    :type: movement_df: DataFrame
    :parameter: begin:  Starting frame.
    :type begin: int
    :parameter: end:  Ending frame (exclude).
    :type end: int
    :parameter: Nsec: number of seconds to bin
    :type: Nsec: int
    :parameter agg_group: Columns to aggregate by e.g ['PlateID', 'ChemID', 'Well', 'Conc']
    :type agg_group: str array
    :parameter verbose: Print info or not
    :type verbose: Boolean
    
    :returns: Movement data aggregated by columns provided and seconds.
                Labeled by the end of the bin (i.e Nsec = 30; labels: 30, 60, 90, ...)
    :rtype: DataFrame

    """

    df = movement_df.copy()
    agg_group.append('second')
    
    # Subset data
    df = df.loc[(df['Trial time'] >= begin) & (df['Trial time'] < end)]

    df['second'] = np.zeros(df.shape[0])
    for second in range(int((begin/Nsec)), int((end/Nsec))):
        # experimental time in seconds
        start = (second*Nsec)
        stop = (second*Nsec)+Nsec

        # Using trial time create a new column with second data
        df.loc[(df['Trial time'] >= start) & (df['Trial time'] < stop), 'second'] = stop
    
    # Aggregate movement data and reindex
    df1 = df.groupby(by=agg_group, 
                     sort=False)['Distance moved'].sum().reset_index()
    df1 = df.groupby(by=agg_group, 
                         sort=False)['Distance moved'].sum().reset_index()
    
    #rename specific column names
    df1.rename(columns = {'second':'Trial time'}, inplace = True)
    

    if verbose:
        print("Start time: ", begin, " Stop time: ", end, "\n")
        print('Total number of minutes: ', str(int((end))-int((begin))))
        print(df1.head())
    
    return df1

################################################################################

def zf_behavior_KStest(zf_movement_data, aggregation=0, verbose=0):
    """Finds the mean and KS test p-value for zebrafish behaviour during acclimation, light, and dark periods.

    :parameter: zf_movement_data: zebrafish behaviour data.
    :type: zf_movement_data: DataFrame
    :parameter: aggregation: Aggregate movement data into time bins. Frame (0), second (1), minute (2), fish (3).
    :type: aggregation: int
    :parameter verbose: Print info or not
    :type verbose: Boolean
    
    :returns: Mean distance moved and KS results per chemical, concentration, and period.
    :rtype: DataFrame

    """
    
    ks_results = pd.DataFrame([], columns=['ChemID', 'Conc', 'L_pvalue', 'D_pvalue',
                                          'LC_mean', 'DC_mean', 'LT_mean', 'DT_mean'])
    zf_movement = zf_movement_data.copy()
    
    for chem in zf_movement['ChemID'].unique():
        # retreive control data for each experimental period per chemical and calcualate the means
        L_control_df = zf_movement.loc[(zf_movement['ChemID'] == chem) & 
                               (zf_movement['Trial time'] >= 720) & (zf_movement['Trial time'] < 900) &
                               (zf_movement['Conc'] == 0.0)]
        LC_mean = L_control_df['Distance moved'].mean()
        D_control_df = zf_movement.loc[(zf_movement['ChemID'] == chem) & 
                               (zf_movement['Trial time'] >= 900) & (zf_movement['Trial time'] < 1080) &
                               (zf_movement['Conc'] == 0.0)]
        DC_mean = D_control_df['Distance moved'].mean()
        for conc in zf_movement.loc[(zf_movement['ChemID'] == chem)].loc[:,'Conc'].unique():               
            if conc != 0.0:
                # retreive treated data for each experimental period per chemical and calcualate the means
                L_treated_df = zf_movement.loc[(zf_movement['ChemID'] == chem) & 
                                       (zf_movement['Trial time'] >= 720) & (zf_movement['Trial time'] < 900) &
                                       (zf_movement['Conc'] == conc)]
                LT_mean = L_treated_df['Distance moved'].mean()
                D_treated_df = zf_movement.loc[(zf_movement['ChemID'] == chem) & 
                                       (zf_movement['Trial time'] >= 900) & (zf_movement['Trial time'] < 1080) &
                                       (zf_movement['Conc'] == conc)] 
                DT_mean = D_treated_df['Distance moved'].mean()

                # Aggregate time
                if aggregation == 1:
                    # run KS test
                    [L_test, L_pvalue] = ks_2samp(agg_by_time(L_control_df, 720, 1080, 1, 
                                                              ['PlateID', 'ChemID', 'Well', 'Conc'], verbose)['Distance moved'],
                                                  agg_by_time(L_treated_df, 720, 1080, 1, 
                                                              ['PlateID', 'ChemID', 'Well', 'Conc'], verbose)['Distance moved'])
                    [D_test, D_pvalue] = ks_2samp(agg_by_time(D_control_df, 720, 1080, 1, 
                                                              ['PlateID', 'ChemID', 'Well', 'Conc'], verbose)['Distance moved'],
                                                  agg_by_time(D_treated_df, 720, 1080, 1, 
                                                              ['PlateID', 'ChemID', 'Well', 'Conc'], verbose)['Distance moved'])
                elif aggregation == 2:
                    # run KS test
                    [L_test, L_pvalue] = ks_2samp(agg_by_min(L_control_df, 720, 1080, 60, 
                                                              ['PlateID', 'ChemID', 'Well', 'Conc'], verbose)['Distance moved'],
                                                  agg_by_min(L_treated_df, 720, 1080, 60, 
                                                              ['PlateID', 'ChemID', 'Well', 'Conc'], verbose)['Distance moved'])
                    [D_test, D_pvalue] = ks_2samp(agg_by_min(D_control_df, 720, 1080, 60, 
                                                              ['PlateID', 'ChemID', 'Well', 'Conc'], verbose)['Distance moved'],
                                                  agg_by_min(D_treated_df, 720, 1080, 60, 
                                                              ['PlateID', 'ChemID', 'Well', 'Conc'], verbose)['Distance moved'])
                elif aggregation == 3:
                    # run KS test
                    [L_test, L_pvalue] = ks_2samp(agg_by_zf(L_control_df, 720, 1080, verbose)['Distance moved'],
                                                  agg_by_zf(L_treated_df, 720, 1080, verbose)['Distance moved'])
                    [D_test, D_pvalue] = ks_2samp(agg_by_zf(D_control_df, 720, 1080, verbose)['Distance moved'],
                                                  agg_by_zf(D_treated_df, 720, 1080, verbose)['Distance moved'])
                else:
                    # run KS test
                    [L_test, L_pvalue] = ks_2samp(L_control_df['Distance moved'],L_treated_df['Distance moved'])
                    [D_test, D_pvalue] = ks_2samp(D_control_df['Distance moved'],D_treated_df['Distance moved'])

                if verbose:
                    print("Controls: ", L_control_df.shape, " Treated: ", L_treated_df.shape)
                ks_results = ks_results.append({'ChemID':chem, 'Conc':conc, 
                                                'L_pvalue':L_pvalue, 'D_pvalue':D_pvalue,
                                                'LC_mean':LC_mean, 'DC_mean':DC_mean, 
                                                'LT_mean':LT_mean, 'DT_mean':DT_mean},
                                                         ignore_index=True)
    return ks_results

################################################################################

def zf_Control_behavior_KStest(zf_movement_1, zf_movement_2, aggregation=0, verbose=0):
    """Finds the mean and KS test p-value for zebrafish behaviour during second light, and dark periods 
       between dataframes provided.

    :parameter: zf_movement_1: zebrafish behaviour data.
    :type: zf_movement_1: DataFrame
    :parameter: zf_movement_2: zebrafish behaviour data.
    :type: zf_movement_2: DataFrame
    :parameter: aggregation: Aggregate movement data into time bins. Frame (0), second (1), minute (2), fish (3).
    :type: aggregation: int
    :parameter verbose: Print info or not
    :type verbose: Boolean
    
    :returns: Mean distance moved and KS results per chemical, concentration, and period.
    :rtype: DataFrame

    """
    
    ks_results = pd.DataFrame([], columns=['L_pvalue', 'D_pvalue',
                                           'LC1_mean', 'DC1_mean', 
                                           'LC2_mean', 'DC2_mean'])
    zf_movement_1_data = zf_movement_1.copy()
    zf_movement_2_data = zf_movement_2.copy()
    
    # retreive control data for each experimental period per chemical and calcualate the means
    L_control_df = zf_movement_1_data.loc[(zf_movement_1_data['Trial time'] >= 720) & 
                                          (zf_movement_1_data['Trial time'] < 900)]
    LC1_mean = L_control_df['Distance moved'].mean()
    D_control_df = zf_movement_1_data.loc[(zf_movement_1_data['Trial time'] >= 900) & 
                                          (zf_movement_1_data['Trial time'] < 1080)]
    DC1_mean = D_control_df['Distance moved'].mean()

    # retreive treated data for each experimental period per chemical and calcualate the means
    L_treated_df = zf_movement_2_data.loc[(zf_movement_2_data['Trial time'] >= 720) & 
                                          (zf_movement_2_data['Trial time'] < 900)]
    LC2_mean = L_treated_df['Distance moved'].mean()
    D_treated_df = zf_movement_2_data.loc[(zf_movement_2_data['Trial time'] >= 900) & 
                                          (zf_movement_2_data['Trial time'] < 1080)] 
    DC2_mean = D_treated_df['Distance moved'].mean()

    # Aggregate time
    if aggregation == 1:
        # run KS test
        [L_test, L_pvalue] = ks_2samp(agg_by_sec(L_control_df, 720, 1080, verbose)['Distance moved'],
                                      agg_by_sec(L_treated_df, 720, 1080, verbose)['Distance moved'])
        [D_test, D_pvalue] = ks_2samp(agg_by_sec(D_control_df, 720, 1080, verbose)['Distance moved'],
                                      agg_by_sec(D_treated_df, 720, 1080, verbose)['Distance moved'])
    elif aggregation == 2:
        # run KS test
        [L_test, L_pvalue] = ks_2samp(agg_by_min(L_control_df, 720, 1080, verbose)['Distance moved'],
                                      agg_by_min(L_treated_df, 720, 1080, verbose)['Distance moved'])
        [D_test, D_pvalue] = ks_2samp(agg_by_min(D_control_df, 720, 1080, verbose)['Distance moved'],
                                      agg_by_min(D_treated_df, 720, 1080, verbose)['Distance moved'])
    elif aggregation == 3:
        # run KS test
        [L_test, L_pvalue] = ks_2samp(agg_by_zf(L_control_df, 720, 1080, verbose)['Distance moved'],
                                      agg_by_zf(L_treated_df, 720, 1080, verbose)['Distance moved'])
        [D_test, D_pvalue] = ks_2samp(agg_by_zf(D_control_df, 720, 1080, verbose)['Distance moved'],
                                      agg_by_zf(D_treated_df, 720, 1080, verbose)['Distance moved'])
    else:
        # run KS test
        [A_test, A_pvalue] = ks_2samp(A_control_df['Distance moved'],A_treated_df['Distance moved'])
        [L_test, L_pvalue] = ks_2samp(L_control_df['Distance moved'],L_treated_df['Distance moved'])
        [D_test, D_pvalue] = ks_2samp(D_control_df['Distance moved'],D_treated_df['Distance moved'])

    if verbose:
        print("Control 1: ", L_control_df.shape, " Control 2: ", L_treated_df.shape)
    ks_results = ks_results.append({'L_pvalue':L_pvalue, 'D_pvalue':D_pvalue,
                                    'LC1_mean':LC1_mean, 'DC1_mean':DC1_mean, 
                                    'LC2_mean':LC2_mean, 'DC2_mean':DC2_mean},
                                             ignore_index=True)
    return ks_results

################################################################################

def diff_behavior(movement, results_df, activity_state, stringency = 1):
    """Finds the fish with a given activity state in the light and dark.

    :parameter movement: zebrafish movement data.
    :type movement: DataFrame
    :parameter results_df: significant treatments.
    :type results_df: DataFrame
    :parameter activity_state: zebrafish activity state (hyperactive, hypoactive).
    :type activity_state: str
    :parameter stringency: How strict should be activity call be? Significant activity in either light or dark (0),
                            (1) Significant activity in both light and dark.
    :type stringency: int

    :returns: Mean distance moved for each fish per chemical, and concentration along with corresponding control mean.
    :rtype: DataFrame

    """

    # Initialize output variable
    hyperactive_zfs= pd.DataFrame([], columns=['PlateID', 'ChemID', 'Well', 'Conc', 'meanDM', 'C_meanDM'])  
    if stringency:
        hyperactive_chems= results_df.loc[(results_df['L_Activity'] == activity_state) &
                                               (results_df['D_Activity'] == activity_state)].loc[:,['ChemID', 'Conc']]
    else:
        hyperactive_chems= results_df.loc[(results_df['L_Activity'] == activity_state) |
                                               (results_df['D_Activity'] == activity_state)].loc[:,['ChemID', 'Conc']]
    for i in range(len(hyperactive_chems)):
        chem = hyperactive_chems.iloc[i,0]
        conc = hyperactive_chems.iloc[i,1]
        hyperactive_zf = movement.loc[(movement['ChemID'] == chem) & (movement['Conc'] == conc)]
        control_zf_mean = movement.loc[(movement['ChemID'] == chem) & 
                                           (movement['Conc'] == 0.0)].loc[:,'Distance moved'].mean()
        for plate in hyperactive_zf['PlateID'].unique():
            for well in hyperactive_zf.loc[(hyperactive_zf['PlateID'] == plate)].loc[:,'Well'].unique():
                hyperactive_zf_mean = hyperactive_zf.loc[(hyperactive_zf['PlateID'] == plate) & 
                                                      (hyperactive_zf['Well'] == well)].loc[:,'Distance moved'].mean()
                hyperactive_zfs = hyperactive_zfs.append({'PlateID': plate, 'ChemID':chem, 'Well': well, 'Conc':conc, 
                                                          'meanDM':hyperactive_zf_mean, 'C_meanDM':control_zf_mean},
                                                             ignore_index=True)
    return hyperactive_zfs

################################################################################

def phase_diff_behavior(movement, results_df, activity_state, phase = None):
    """Finds the fish with a given activity state in the light and dark.

    :parameter movement: zebrafish movement data.
    :type movement: DataFrame
    :parameter results_df: significant treatments.
    :type results_df: DataFrame
    :parameter activity_state: zebrafish activity state (hyperactive, hypoactive).
    :type activity_state: str
    :parameter Phase: 'Light' or 'Dark' phase of experiment
    :type stringency: str

    :returns: Mean distance moved for each fish per chemical, and concentration along with corresponding control mean.
    :rtype: DataFrame

    """

    # Initialize output variable
    active_zfs= pd.DataFrame([], columns=['PlateID', 'ChemID', 'Well', 'Conc', 'meanDM', 'C_meanDM'])  
    if phase == 'Light':
        active_chems= results_df.loc[(results_df['L_Activity'] == activity_state)].loc[:,['ChemID', 'Conc']]
    elif phase == 'Dark':
        active_chems= results_df.loc[(results_df['D_Activity'] == activity_state)].loc[:,['ChemID', 'Conc']]
    else:
        active_chems= results_df.loc[(results_df['L_Activity'] == activity_state) &
                                               (results_df['D_Activity'] == activity_state)].loc[:,['ChemID', 'Conc']]
    for i in range(len(active_chems)):
        chem = active_chems.iloc[i,0]
        conc = active_chems.iloc[i,1]
        active_zf = movement.loc[(movement['ChemID'] == chem) & (movement['Conc'] == conc)]
        control_zf_mean = movement.loc[(movement['ChemID'] == chem) & 
                                           (movement['Conc'] == 0.0)].loc[:,'Distance moved'].mean()
        for plate in active_zf['PlateID'].unique():
            for well in active_zf.loc[(active_zf['PlateID'] == plate)].loc[:,'Well'].unique():
                active_zf_mean = active_zf.loc[(active_zf['PlateID'] == plate) & 
                                                      (active_zf['Well'] == well)].loc[:,'Distance moved'].mean()
                active_zfs = active_zfs.append({'PlateID': plate, 'ChemID':chem, 'Well': well, 'Conc':conc, 
                                                          'meanDM':active_zf_mean, 'C_meanDM':control_zf_mean},
                                                             ignore_index=True)
    return active_zfs

################################################################################

def find_abnormal(zf_activity, ouput_df, state='Normal', phase='None', verbose=0):
    """Assign ClassID to plateTable and movement output for heatmap and ML analysis.

    :parameter zf_activity: zebrafish activity data file (diff_behavior() function output).
    :type zf_activity: DataFrame
    :parameter ouput_df: DataFrame to which ClassID will be added
    :type ouput_df: DataFrame
    :parameter state: String that indicates the activity state (hyper- vs hypoactive)
    :type state: Str
    :parameter Phase: 'Light' or 'Dark' phase of experiment
    :type stringency: str
    :parameter verbose: Print info or not
    :type verbose: Boolean
    
    :returns: ouput_df with ClassID added
    :rtype: DataFrame

    """    
    active_zfs = pd.DataFrame([], columns=['PlateID', 'ChemID', 'Well', 'Conc'])  
    i = None
    for i in range(round((len(zf_activity)*0.3))):
        PlateID = zf_activity.loc[i, 'PlateID']
        ChemID = zf_activity.loc[i, 'ChemID']
        Well = zf_activity.loc[i, 'Well']
        Conc = zf_activity.loc[i, 'Conc']
        ouput_df.loc[(ouput_df['PlateID'] == PlateID) & (ouput_df['ChemID'] == ChemID) & 
                     (ouput_df['Well'] == Well) & (ouput_df['Conc'] == Conc), 'ClassID_'+phase] = 1
        ouput_df.loc[(ouput_df['PlateID'] == PlateID) & (ouput_df['ChemID'] == ChemID) & 
                     (ouput_df['Well'] == Well) & (ouput_df['Conc'] == Conc), 'Activity_'+phase] = state

        if verbose:
            active_zfs = active_zfs.append(zf_activity.loc[(zf_activity['PlateID'] == PlateID) & 
                                                           (zf_activity['ChemID'] == ChemID) & 
                                                           (zf_activity['Well'] == Well) &
                                                           (zf_activity['Conc'] == Conc)])
#             print("Plate: ", PlateID, " Chem: ", ChemID, " Well: ", Well, " Conc: ", Conc)
    if i == None:
        i = 0
    else:
        i+=1
    if verbose:
        print(active_zfs.head(active_zfs.shape[0]))
    
    return ouput_df, i

################################################################################
############################### DATA ANALYSIS ##################################
################################################################################

def fdr(p_vals, rate):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    
    :parameter p_vals: An array of p-values
    :type pending: np.array  
    
    :parameter rate: Acceptable false dicovery rate
    :type pending: float     
   
    :returns: BH adjusted p-value
    :rtype: np.array
    """
    
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > rate] = 1
    return fdr

################################################################################
############################### DATA PLOTTING ##################################
################################################################################

def tsne_plot(x, y,  Save_path=None): 
    
    """Perform t-distributed Stochastic Neighbor Embedding and plot data.
    
    :parameter x: features.
            :type: np.array
    :parameter y: class labels.
            :type: np.array
    :parameter Save: name and path were the tSNE plot should be saved to file.
            :type: str
    """

    import seaborn as sns 
    from sklearn.manifold import TSNE
    import matplotlib.pyplot as plt 
    
    # Setting the plotting background 
    sns.set(style ="whitegrid") 
    
    tsne = TSNE(n_components = 2, random_state = 0) 
    
    # Reducing the dimensionality of the data 
    X_transformed = tsne.fit_transform(x)
    
    plt.figure(figsize =(12, 8)) 
    
    # Building the scatter plot 
    plt.scatter(X_transformed[np.where(y == 0), 0], 
                X_transformed[np.where(y == 0), 1], 
                marker ='o', color ='y', linewidth ='1', 
                alpha = 0.8, label ='Normal') 
    plt.scatter(X_transformed[np.where(y == 1), 0], 
                X_transformed[np.where(y == 1), 1], 
                marker ='o', color ='k', linewidth ='1', 
                alpha = 0.8, label ='Abnormal') 

    # Specifying the location of the legend 
    plt.legend(loc ='best') 
    
    if Save_path is not None:            
        plt.savefig(Save_path, dpi=600, bbox_inches='tight' )   
        
    # Plotting the reduced data 
    plt.show() 

################################################################################

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

################################################################################

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

################################################################################
############################### MISCELLANEOUS ##################################
################################################################################

def round_sig(x, sig=2):
    """
    Returns a value rounded to a specific number of significant figures.

    ----------
    x: number
    sig: number of significant figures
    ----------
    """
    return round(x, sig-int(math.floor(math.log10(abs(x))))-1)

    
################################################################################
############################## BUILD AUTOENCODER ################################
################################################################################

def build_autoencoder(data_shape, parameters):
    """Construct an autoencoder using the given inputs
    
    :parameter data_shape: The shape of the input data i.e (64, 64)
    :type data_shape: tuple
    :parameter parameters: network parameters [HP_EncoderLayers, HP_DecoderLayers, HP_Kern, HP_EDrop,
                                                HP_DDrop, HP_OPTIMIZER, HP_ACTIVATION, HP_L2, HP_L1]
    :type data_shape: narray
    
    :returns: The autiencoder and encoder
    :rtype: Sequential
    
    """
    scaler = 15
    neurons = int(np.prod(data_shape))
    # Building the Input Layer 
    input_layer = Input(shape = data_shape) 

    # Building the Encoder network 
    encoder = Flatten(name='Hidden1')(input_layer)
    for i in range(parameters['HP_EncoderLayers']):
        layer_name = 'Hidden'+str(i+2)
        encoder = Dense(neurons//(scaler*(i+1)), activation=parameters['HP_ACTIVATION'], 
                        activity_regularizer = regularizers.l1_l2(l1=parameters['HP_L1'], l2=parameters['HP_L2']), 
                        kernel_initializer = parameters['HP_Kern'], name=layer_name)(encoder)
        if parameters['HP_EDrop'] < 1:
            encoder = Dropout(parameters['HP_EDrop'])(encoder)

    Dneurons = neurons//(scaler*parameters['HP_DecoderLayers'])
    for i in range(parameters['HP_DecoderLayers']):
        if i == 0:
            decoded = Dense(Dneurons*(i+1), activation=parameters['HP_ACTIVATION'])(encoder) 
        elif (i+1) < parameters['HP_DecoderLayers']:
            decoded = Dense(Dneurons*(i+1), activation=parameters['HP_ACTIVATION'])(decoded)
        else:
            decoded = Dense(neurons)(decoded) 
        if parameters['HP_DDrop'] < 1:
            decoded = Dropout(parameters['HP_DDrop'])(decoded)
    # Building the Output Layer 
    output_layer = Reshape(data_shape)(decoded)

    # Defining the parameters of the Auto-encoder network 
    autoencoder = Model(input_layer, output_layer) 
    autoencoder.compile(optimizer=parameters['HP_OPTIMIZER'], loss ="mse") 


    hidden_representation = Sequential() 
    for i in range(parameters['HP_EncoderLayers']):
        hidden_representation.add(autoencoder.layers[i]) 

    return autoencoder, hidden_representation

################################################################################

def build_encoder(data_shape, parameters):
    """Construct an encoder using the given inputs
    
    :parameter data_shape: The shape of the input data i.e (64, 64)
    :type data_shape: tuple
    :parameter parameters: network parameters [HP_EncoderLayers, HP_Kern, HP_EDrop, 
                                                HP_OPTIMIZER, HP_ACTIVATION, HP_L2, HP_L1]
    :type data_shape: narray
    
    :returns: The encoder
    :rtype: Sequential
    
    """
    scaler = 15
    neurons = int(np.prod(data_shape))
    # Building the Input Layer 
    input_layer = Input(shape = data_shape) 

    # Building the Encoder network 
    encoder = Flatten(name='Hidden1')(input_layer)
    for i in range(parameters['HP_EncoderLayers']):
        layer_name = 'Hidden'+str(i+2)
        encoder = Dense(neurons//(scaler*(i+1)), activation=parameters['HP_ACTIVATION'], 
                        activity_regularizer = regularizers.l1_l2(l1=parameters['HP_L1'], l2=parameters['HP_L2']), 
                        kernel_initializer = parameters['HP_Kern'], name=layer_name)(encoder)
        if parameters['HP_EDrop'] < 1:
            encoder = Dropout(parameters['HP_EDrop'])(encoder)

    encoder = Dense(1, activation='softmax')(encoder)
    encoder = Model(input_layer, encoder) 
    encoder.compile(optimizer=parameters['HP_OPTIMIZER'], loss ="mse") 

    return encoder