#!/usr/bin/env python
# coding: utf-8

# In[1]:


## classification of transitory myeloid cells 
## Logistic regression, XGBoost, RandomForest 

#import cuml.accel
runs = 5 

import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder

#from sklearn.preprocessing import OneHotEncoder
#from sklearn.preprocessing import LabelBinarizer

print('Reading data ...')
data=pd.read_csv("../MLdata_TNK/Merged_DomTNK_data_AllSamples_LUAD.csv")
#print(data)

data_cleaned = data.dropna(subset=['label']) 


# In[2]:


print('Encoding classes ...')

mapping = {'NoTransMyl': 0, 'WithTransMyl': 1}
data_cleaned.loc[:,'status_encoded'] = data_cleaned.loc[:,'label'].map(mapping)
#print(data_cleaned)
numClasses = 2

#le = LabelEncoder()
#data_cleaned['status_encoded'] = le.fit_transform(data_cleaned['label']) 

#mapping = {label: index for index, label in enumerate(le.classes_)}
#numClasses = len(le.classes_)
#mapping = {label: index for index, label in enumerate(le.classes_)}
#print("Encoded Mapping:", mapping)

print(data_cleaned['status_encoded'].value_counts()) 

mapping = None 


# In[3]:


## Load
from sklearn.feature_selection import SelectKBest, f_classif, RFECV, SelectFromModel, VarianceThreshold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.inspection import permutation_importance
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.preprocessing import QuantileTransformer
from sklearn.pipeline import Pipeline
import xgboost as xgb

from sklearn.metrics import classification_report, balanced_accuracy_score,accuracy_score, confusion_matrix, f1_score, precision_score, recall_score, roc_auc_score
import matplotlib.pyplot as plt
import statistics
import shap

scaler = StandardScaler() 
scaler.set_output(transform="pandas")
normalizer = QuantileTransformer(output_distribution='normal', n_quantiles=4)


# In[ ]:


## =========== XGBoost ================

#clf = RandomForestClassifier(class_weight="balanced")
lrm = LogisticRegression(solver='lbfgs', C=0.1, max_iter=1000)
rfecv = RFECV(estimator=lrm, step=1, cv=StratifiedKFold(5), scoring='roc_auc')

pipeline = Pipeline(
    steps=[
        ("scaling",scaler),
        ("norm",normalizer),
        ("variance", VarianceThreshold()),
        ("feature_sel1",SelectKBest(score_func=f_classif, k=200)),
        ("feature_sel2",SelectFromModel(LogisticRegression(solver='lbfgs', C=0.1, max_iter=1000))),
        ("feature_sel3",rfecv),
        ("xgb_model", xgb.XGBClassifier(tree_method="hist"))
    ]
)

param_grid = {
    "xgb_model__max_depth": [ 4, 6, 8 ], 
    "xgb_model__n_estimators": [ 50, 100, 150], 
    "xgb_model__learning_rate": [ 0.1, 0.05]
}

wrname = '../results_MergedTNK/Results_MLrun_Merged_LUAD_run_XGB.txt'
wrf = open(wrname,'w')

acc_scores = []
acc_scores_bal = []
f1_scores = []
prec_scores = []
rec_scores = []
#roc_scores = []

nat_imp = {}
perm_imp = {}
shap_imp = {}


print('Sampling data ...')

for k in range(0,runs):
    print(k)
    wrf.write('Round: {}\n'.format(k))

    df1 = data_cleaned.loc[data_cleaned['status_encoded'] == 0]
    df2 = data_cleaned.loc[data_cleaned['status_encoded'] == 1]

    min_val = min(len(df1),len(df2))
    min_val=200

    df1_sample = df1.sample(n=min_val)
    df1 = None 
    df2_sample = df2.sample(n=min_val)
    df2 = None 

    final_df = pd.concat([df1_sample, df2_sample], axis=0)
    df1_sample = None 
    df2_sample = None 


    print(final_df['status_encoded'].value_counts()) 

    xnames=final_df.columns[1:len(final_df.columns)-2]
    #print(xnames)
    X=final_df[xnames]
    y=final_df.status_encoded

    xnames = None 
    final_df = None  

    train_X, test_X, train_y, test_y = train_test_split(X, y, stratify=y, train_size=0.75)

    #sample_weights = compute_sample_weight(
    #    class_weight='balanced',
    #    y = train_y #provide your own target name
    #)

    #train_X_embedded, test_X_embedded = feature_selection1(train_X_scaled,train_y,test_X_scaled)
    #train_X_embedded, test_X_embedded = feature_selection2(train_X_scaled,train_y,test_X_scaled)

    grid_search = GridSearchCV(estimator=pipeline, param_grid=param_grid, cv=5, scoring='accuracy', verbose=1, n_jobs=-1)
    grid_search.fit(train_X, train_y) #,xgb_model__sample_weight=sample_weights) 

    best_pipeline = grid_search.best_estimator_
    best_params = grid_search.best_params_

    print(f"Best hyperparameters: {best_params}")
    wrf.write('  Best hyperparameters: {}\n'.format(best_params))
    #print(f"Best score (roc_auc): {grid_search.best_score_}")

    #pred_values=best_model.predict(test_X)
    pred_values = grid_search.predict(test_X)

    accr = accuracy_score(test_y,pred_values)
    accr_bal = balanced_accuracy_score(test_y,pred_values)
    f1s = f1_score(test_y,pred_values,average='weighted',zero_division = np.nan)
    prec = precision_score(test_y,pred_values,average='weighted',zero_division = np.nan)
    rec = recall_score(test_y,pred_values,average='weighted',zero_division = np.nan)
    #roc = roc_auc_score(test_y,pred_values,average='macro',multi_class='ovo')

    print(confusion_matrix(test_y, pred_values))
    wrf.write('  {}\n'.format(confusion_matrix(test_y, pred_values)))

    acc_scores.append(accr)
    acc_scores_bal.append(accr_bal)
    f1_scores.append(f1s)
    prec_scores.append(prec)
    rec_scores.append(rec)
    #roc_scores.append(roc)

    accr = None
    accr_bal = None 
    f1s = None 
    prec = None 
    rec = None 
    roc = None 
    #best_model = None 
    #best_params = None 


    ## === Feature importance === 

    col_names = np.array(X.columns)
    mask1 = best_pipeline.named_steps['variance'].get_support()
    mask2 = best_pipeline.named_steps['feature_sel1'].get_support()
    mask3 = best_pipeline.named_steps['feature_sel2'].get_support()
    mask4 = best_pipeline.named_steps['feature_sel3'].get_support()
    features_mask1 = col_names[mask1]
    features_mask2 = features_mask1[mask2]
    features_mask3 = features_mask2[mask3]
    feature_names = features_mask3[mask4]
    ## print("Selected Features:", feature_names)

    ## Native method
    print('\tNative importance')
    best_model = best_pipeline.named_steps['xgb_model']
    importance = best_model.feature_importances_

    #print(importance)
    for lk,feature in enumerate(feature_names): 
        if feature in nat_imp.keys(): 
            tmp = nat_imp[feature]
            tmp.append(importance[lk])
            nat_imp[feature] = tmp
            tmp = None 
        else:
            nat_imp[feature] = [importance[lk]]

    importance = None 
    lk = None 
    feature = None 

    ## Permutaiton method 
    print('\tPermutation importance')
    test_X_select = test_X[feature_names]
    result = permutation_importance(best_model, test_X_select.values, test_y, n_repeats=10, random_state=42)
    importance = result.importances_mean
    for lk,feature in enumerate(feature_names): 
        if feature in perm_imp.keys(): 
            tmp = perm_imp[feature]
            tmp.append(importance[lk])
            perm_imp[feature] = tmp
            tmp = None 
        else:
            perm_imp[feature] = [importance[lk]]    

    importance = None 
    lk = None 
    feature = None 
    result = None 

   ## SHAP 

    print('\tSHAP importance')
    explainer = shap.Explainer(best_model.predict_proba,test_X_select.values)
    shap_values = explainer(test_X_select.values)
    shap_df1 = pd.DataFrame(shap_values.values[:, :, 0], columns=test_X_select.columns)
    shap_df2 = pd.DataFrame(shap_values.values[:, :, 1], columns=test_X_select.columns)

    ## shap_df_sum = shap_df1.abs() + shap_df2.abs()

    importance_class0 = shap_df1.abs().mean()
    importance_class1 = shap_df2.abs().mean()

    ##print(type(importance))
    ##importance_df = pd.DataFrame({
    ##    'Feature': feature_names,
    ##    'Importance': shap_df_sum.abs().mean()
    ##})

    for lk,feature in enumerate(feature_names): 
        if feature in shap_imp.keys(): 
            tmp = shap_imp[feature]
            tmp.append([importance_class0.iloc[lk],importance_class1.iloc[lk]])
            shap_imp[feature] = tmp
            tmp = None 
        else:
            shap_imp[feature] = [[importance_class0.iloc[lk],importance_class1.iloc[lk]]]   

    shap_df1 = None 
    shap_df2 = None 
    shap_values = None 
    explainer = None 
    lk = None 
    feature = None 
    importance_class0 = None
    importance_class1 = None


    X = None 
    y = None 
    best_model = None 
    best_params = None
    train_X = None 
    test_X = None 
    train_y = None 
    test_y = None 
    test_X_select = None 



print(f'average_accuracy: {statistics.mean(acc_scores):.2f}, sd_accuracy: {statistics.pstdev(acc_scores):.2f}\n')
print(f'average_accuracy_balanced: {statistics.mean(acc_scores_bal):.2f}, sd_accuracy: {statistics.pstdev(acc_scores_bal):.2f}\n')
print(f'average_f1-score: {statistics.mean(f1_scores):.2f}, sd_f1-score: {statistics.pstdev(f1_scores):.2f}\n')
print(f'average_precision: {statistics.mean(prec_scores):.2f}, sd_precision: {statistics.pstdev(prec_scores):.2f}\n')
print(f'average_recall: {statistics.mean(rec_scores):.2f}, sd_Recall: {statistics.stdev(rec_scores):.2f}\n')
#print(f'average_ROC-AUC_score: {statistics.mean(roc_scores):.2f}, sd_ROC-AUC_score: {statistics.stdev(roc_scores):.2f}\n')
wrf.write('\taverage_accuracy: {}, sd_accuracy: {}\n'.format(round(statistics.mean(acc_scores),2), round(statistics.pstdev(acc_scores),2)))
wrf.write('\taverage_accuracy_balanced: {}, sd_accuracy_balanced: {}\n'.format(round(statistics.mean(acc_scores_bal),2), round(statistics.pstdev(acc_scores_bal),2)))
wrf.write('\taverage_f1-score: {}, sd_f1-score: {}\n'.format(round(statistics.mean(f1_scores),2), round(statistics.pstdev(f1_scores),2)))
wrf.write('\taverage_precision: {}, sd_precision: {}\n'.format(round(statistics.mean(prec_scores),2), round(statistics.pstdev(prec_scores),2)))
wrf.write('\taverage_recall: {}, sd_Recall: {}\n'.format(round(statistics.mean(rec_scores),2), round(statistics.pstdev(rec_scores),2)))

wrf.write("\n ===== Top Features by Native importance method =====\n")
print("Top Features by Native importance method:")

avg_imp = []
sd_imp = []

for key in nat_imp.keys():
    #print(key,nat_imp[key])
    if len(nat_imp[key]) < runs:
        xtlist = [0]*(runs-len(nat_imp[key]))
        tmp = nat_imp[key] + xtlist 
        nat_imp[key] = tmp 
        tmp = None 
        xtlist = None 

    #print(key,nat_imp[key])    
    avg = round(statistics.mean(nat_imp[key]),4)
    ##sd = round(statistics.stdev(nat_imp[key]),2)
    series = pd.Series(nat_imp[key]) 
    sd = round(series.std(ddof=0),4)
    series = None 

    avg_imp.append(avg)
    sd_imp.append(sd)

data = {
    'Features': nat_imp.keys(),
    'Mean_Imp': avg_imp, 
    'SD_Imp': sd_imp, 
}

feature_importance = pd.DataFrame(data)
feature_importance = feature_importance.sort_values(by='Mean_Imp', ascending=False)
print(feature_importance.head(20))
wrf.write('{}\n'.format(feature_importance))

feature_importance.to_csv('../results_MergedTNK/LUAD_TNK_FeatureImp_XGB_Native.tsv',sep='\t',index=False)

avg_imp = None 
sd_imp = None 
nat_imp = None 
data = None 
feature_importance = None 


wrf.write("\n====== Top Features by Permutation Importance ======= \n")
print("\nTop Features by Permutation Importance:")

avg_imp = []
sd_imp = []

for key in perm_imp.keys():
    #print(key,nat_imp[key])
    if len(perm_imp[key]) < runs:
        xtlist = [0]*(runs-len(perm_imp[key]))
        tmp = perm_imp[key] + xtlist 
        perm_imp[key] = tmp 
        tmp = None 
        xtlist = None 

    #print(key,nat_imp[key])    
    avg = round(statistics.mean(perm_imp[key]),4)
    series = pd.Series(perm_imp[key]) 
    sd = round(series.std(ddof=0),4)
    series = None 

    avg_imp.append(avg)
    sd_imp.append(sd)

data = {
    'Features': perm_imp.keys(),
    'Mean_Imp': avg_imp, 
    'SD_Imp': sd_imp, 
}

feature_importance = pd.DataFrame(data)
feature_importance = feature_importance.sort_values(by='Mean_Imp', ascending=False)
print(feature_importance.head(20))
wrf.write('{}\n'.format(feature_importance))

feature_importance.to_csv('../results_MergedTNK/LUAD_TNK_FeatureImp_XGB_Perm.tsv',sep='\t',index=False)

avg_imp = None 
sd_imp = None 
perm_imp = None 
data = None 
feature_importance = None 


print("Top Features by SHAP importance method:")

avg_imp = []
sd_imp = []

for key in shap_imp.keys():
   #print(key,shap_imp[key])
   if len(shap_imp[key]) < runs:
       xtlist = [[0]*numClasses]*(runs-len(shap_imp[key]))
       tmp = shap_imp[key] + xtlist 
       shap_imp[key] = tmp 
       tmp = None 
       xtlist = None 

   #print(key,shap_imp[key])    

   keyavg = []
   keysd = []
   for i in range(0,numClasses):
       vals = []
       for j in range(0,runs):
           #print(shap_imp[key][j][i])
           vals.append(shap_imp[key][j][i])

       avg = round(statistics.mean(vals),4)
       sd = round(statistics.stdev(vals),4)
       keyavg.append(avg)
       keysd.append(sd)
       avg = None 
       sd = None 
   ## print(keyavg,keysd)

   avg_imp.append(keyavg)
   sd_imp.append (keysd)



data = {
   'Features': shap_imp.keys(),
   'Mean_Imp': avg_imp, 
   'SD_Imp': sd_imp, 
}

feature_importance = pd.DataFrame(data)
feature_importance = feature_importance.sort_values(by='Mean_Imp', ascending=False)
print(feature_importance.head(20))

wrf.write('{}\n'.format(feature_importance))

feature_importance.to_csv('../results_MergedTNK/LUAD_TNK_FeatureImp_XGB_SHAP.tsv',sep='\t',index=False)

avg_imp = None 
sd_imp = None 
shap_imp = None 
data = None 
feature_importance = None 

wrf.close()



# In[ ]:


## ===== Random Forest ============= 

#clf = RandomForestClassifier(class_weight="balanced")
lrm = LogisticRegression(solver='lbfgs', C=0.1, max_iter=1000)
rfecv = RFECV(estimator=lrm, step=1, cv=StratifiedKFold(5), scoring='roc_auc')

pipeline = Pipeline(
    steps=[
        ("scaling",scaler),
        ("norm",normalizer),
        ("variance", VarianceThreshold()),
        ("feature_sel1",SelectKBest(score_func=f_classif, k=200)),
        ("feature_sel2",SelectFromModel(LogisticRegression(solver='lbfgs', C=0.1, max_iter=1000))),
        ("feature_sel3",rfecv),
        ("rf_model", RandomForestClassifier(n_estimators=100, max_depth=3))
    ]
)


param_grid = {
    "rf_model__max_depth": [ 4, 6, 8 ], 
    "rf_model__n_estimators": [ 50, 100, 150], 
}

wrname = '../results_MergedTNK/Results_MLrun_Merged_LUAD_run_RF.txt'
wrf = open(wrname,'w')

acc_scores = []
acc_scores_bal = []
f1_scores = []
prec_scores = []
rec_scores = []
#roc_scores = []

nat_imp = {}
perm_imp = {}
shap_imp = {}


print('Sampling data ...')

for k in range(0,runs):
    print(k)
    wrf.write('Round: {}\n'.format(k))

    df1 = data_cleaned.loc[data_cleaned['status_encoded'] == 0]
    df2 = data_cleaned.loc[data_cleaned['status_encoded'] == 1]

    min_val = min(len(df1),len(df2))

    df1_sample = df1.sample(n=min_val)
    df1 = None 
    df2_sample = df2.sample(n=min_val)
    df2 = None 

    final_df = pd.concat([df1_sample, df2_sample], axis=0)
    df1_sample = None 
    df2_sample = None 


    print(final_df['status_encoded'].value_counts()) 

    xnames=final_df.columns[1:len(final_df.columns)-2]
    #print(xnames)
    X=final_df[xnames]
    y=final_df.status_encoded


    xnames = None 
    final_df = None  

    train_X, test_X, train_y, test_y = train_test_split(X, y, stratify=y, train_size=0.75)

    #sample_weights = compute_sample_weight(
    #    class_weight='balanced',
    #    y = train_y #provide your own target name
    #)

    grid_search = GridSearchCV(estimator=pipeline, param_grid=param_grid, cv=5, scoring='accuracy', verbose=1, n_jobs=-1)
    grid_search.fit(train_X, train_y) #,rf_model__sample_weight=sample_weights) 

    best_pipeline = grid_search.best_estimator_
    best_params = grid_search.best_params_

    print(f"Best hyperparameters: {best_params}")
    wrf.write('  Best hyperparameters: {}\n'.format(best_params))
    #print(f"Best score (roc_auc): {grid_search.best_score_}")

    #pred_values=best_model.predict(test_X)
    pred_values = grid_search.predict(test_X)

    accr = accuracy_score(test_y,pred_values)
    accr_bal = balanced_accuracy_score(test_y,pred_values)
    f1s = f1_score(test_y,pred_values,average='weighted',zero_division = np.nan)
    prec = precision_score(test_y,pred_values,average='weighted',zero_division = np.nan)
    rec = recall_score(test_y,pred_values,average='weighted',zero_division = np.nan)
    #roc = roc_auc_score(test_y,pred_values,average='macro',multi_class='ovo')

    print(confusion_matrix(test_y, pred_values))
    wrf.write('  {}\n'.format(confusion_matrix(test_y, pred_values)))

    acc_scores.append(accr)
    acc_scores_bal.append(accr_bal)
    f1_scores.append(f1s)
    prec_scores.append(prec)
    rec_scores.append(rec)
    #roc_scores.append(roc)

    accr = None
    accr_bal = None 
    f1s = None 
    prec = None 
    rec = None 
    roc = None 
    #best_model = None 
    #best_params = None 


    ## === Feature importance === 

    col_names = np.array(X.columns)
    mask1 = best_pipeline.named_steps['variance'].get_support()
    mask2 = best_pipeline.named_steps['feature_sel1'].get_support()
    mask3 = best_pipeline.named_steps['feature_sel2'].get_support()
    mask4 = best_pipeline.named_steps['feature_sel3'].get_support()
    features_mask1 = col_names[mask1]
    features_mask2 = features_mask1[mask2]
    features_mask3 = features_mask2[mask3]
    feature_names = features_mask3[mask4]
    ## print("Selected Features:", feature_names)

    ## Native method
    print('\tNative importance')
    best_model = best_pipeline.named_steps['rf_model']
    importance = best_model.feature_importances_

    #print(importance)
    for lk,feature in enumerate(feature_names): 
        if feature in nat_imp.keys(): 
            tmp = nat_imp[feature]
            tmp.append(importance[lk])
            nat_imp[feature] = tmp
            tmp = None 
        else:
            nat_imp[feature] = [importance[lk]]

    importance = None 
    lk = None 
    feature = None 

    ## Permutaiton method 
    print('\tPermutation importance')
    test_X_select = test_X[feature_names]

    result = permutation_importance(best_model, test_X_select.values, test_y, n_repeats=10, random_state=42)
    importance = result.importances_mean
    for lk,feature in enumerate(feature_names): 
        if feature in perm_imp.keys(): 
            tmp = perm_imp[feature]
            tmp.append(importance[lk])
            perm_imp[feature] = tmp
            tmp = None 
        else:
            perm_imp[feature] = [importance[lk]]    

    importance = None 
    lk = None 
    feature = None 
    result = None 

   ## SHAP 

    print('\tSHAP importance')
    explainer = shap.Explainer(best_model.predict_proba,test_X_select.values)
    shap_values = explainer(test_X_select.values)
    shap_df1 = pd.DataFrame(shap_values.values[:, :, 0], columns=test_X_select.columns)
    shap_df2 = pd.DataFrame(shap_values.values[:, :, 1], columns=test_X_select.columns)

    ## shap_df_sum = shap_df1.abs() + shap_df2.abs()

    importance_class0 = shap_df1.abs().mean()
    importance_class1 = shap_df2.abs().mean()

    ##print(type(importance))
    ##importance_df = pd.DataFrame({
    ##    'Feature': feature_names,
    ##    'Importance': shap_df_sum.abs().mean()
    ##})

    for lk,feature in enumerate(feature_names): 
        if feature in shap_imp.keys(): 
            tmp = shap_imp[feature]
            tmp.append([importance_class0.iloc[lk],importance_class1.iloc[lk]])
            shap_imp[feature] = tmp
            tmp = None 
        else:
            shap_imp[feature] = [[importance_class0.iloc[lk],importance_class1.iloc[lk]]]   

    shap_df1 = None 
    shap_df2 = None 
    shap_values = None 
    explainer = None 
    lk = None 
    feature = None 
    importance_class0 = None
    importance_class1 = None


    X = None 
    y = None 
    best_model = None 
    best_params = None
    train_X = None 
    test_X = None 
    train_y = None 
    test_y = None 
    test_X_select = None 



print(f'average_accuracy: {statistics.mean(acc_scores):.2f}, sd_accuracy: {statistics.pstdev(acc_scores):.2f}\n')
print(f'average_accuracy_balanced: {statistics.mean(acc_scores_bal):.2f}, sd_accuracy: {statistics.pstdev(acc_scores_bal):.2f}\n')
print(f'average_f1-score: {statistics.mean(f1_scores):.2f}, sd_f1-score: {statistics.pstdev(f1_scores):.2f}\n')
print(f'average_precision: {statistics.mean(prec_scores):.2f}, sd_precision: {statistics.pstdev(prec_scores):.2f}\n')
print(f'average_recall: {statistics.mean(rec_scores):.2f}, sd_Recall: {statistics.stdev(rec_scores):.2f}\n')
#print(f'average_ROC-AUC_score: {statistics.mean(roc_scores):.2f}, sd_ROC-AUC_score: {statistics.stdev(roc_scores):.2f}\n')
wrf.write('\taverage_accuracy: {}, sd_accuracy: {}\n'.format(round(statistics.mean(acc_scores),2), round(statistics.pstdev(acc_scores),2)))
wrf.write('\taverage_accuracy_balanced: {}, sd_accuracy_balanced: {}\n'.format(round(statistics.mean(acc_scores_bal),2), round(statistics.pstdev(acc_scores_bal),2)))
wrf.write('\taverage_f1-score: {}, sd_f1-score: {}\n'.format(round(statistics.mean(f1_scores),2), round(statistics.pstdev(f1_scores),2)))
wrf.write('\taverage_precision: {}, sd_precision: {}\n'.format(round(statistics.mean(prec_scores),2), round(statistics.pstdev(prec_scores),2)))
wrf.write('\taverage_recall: {}, sd_Recall: {}\n'.format(round(statistics.mean(rec_scores),2), round(statistics.pstdev(rec_scores),2)))

wrf.write("\n ===== Top Features by Native importance method =====\n")
print("Top Features by Native importance method:")

avg_imp = []
sd_imp = []

for key in nat_imp.keys():
    #print(key,nat_imp[key])
    if len(nat_imp[key]) < runs:
        xtlist = [0]*(runs-len(nat_imp[key]))
        tmp = nat_imp[key] + xtlist 
        nat_imp[key] = tmp 
        tmp = None 
        xtlist = None 

    #print(key,nat_imp[key])    
    avg = round(statistics.mean(nat_imp[key]),4)
    ##sd = round(statistics.stdev(nat_imp[key]),2)
    series = pd.Series(nat_imp[key]) 
    sd = round(series.std(ddof=0),4)
    series = None 

    avg_imp.append(avg)
    sd_imp.append(sd)

data = {
    'Features': nat_imp.keys(),
    'Mean_Imp': avg_imp, 
    'SD_Imp': sd_imp, 
}

feature_importance = pd.DataFrame(data)
feature_importance = feature_importance.sort_values(by='Mean_Imp', ascending=False)
print(feature_importance.head(20))
wrf.write('{}\n'.format(feature_importance))

feature_importance.to_csv('../results_MergedTNK/LUAD_TNK_FeatureImp_RF_Native.tsv',sep='\t',index=False)

avg_imp = None 
sd_imp = None 
nat_imp = None 
data = None 
feature_importance = None 


wrf.write("\n====== Top Features by Permutation Importance ======= \n")
print("\nTop Features by Permutation Importance:")

avg_imp = []
sd_imp = []

for key in perm_imp.keys():
    #print(key,nat_imp[key])
    if len(perm_imp[key]) < runs:
        xtlist = [0]*(runs-len(perm_imp[key]))
        tmp = perm_imp[key] + xtlist 
        perm_imp[key] = tmp 
        tmp = None 
        xtlist = None 

    #print(key,nat_imp[key])    
    avg = round(statistics.mean(perm_imp[key]),4)
    series = pd.Series(perm_imp[key]) 
    sd = round(series.std(ddof=0),4)
    series = None 

    avg_imp.append(avg)
    sd_imp.append(sd)

data = {
    'Features': perm_imp.keys(),
    'Mean_Imp': avg_imp, 
    'SD_Imp': sd_imp, 
}

feature_importance = pd.DataFrame(data)
feature_importance = feature_importance.sort_values(by='Mean_Imp', ascending=False)
print(feature_importance.head(20))
wrf.write('{}\n'.format(feature_importance))

feature_importance.to_csv('../results_MergedTNK/LUAD_TNK_FeatureImp_RF_Perm.tsv',sep='\t',index=False)

avg_imp = None 
sd_imp = None 
perm_imp = None 
data = None 
feature_importance = None 


print("Top Features by SHAP importance method:")

avg_imp = []
sd_imp = []

for key in shap_imp.keys():
   #print(key,shap_imp[key])
   if len(shap_imp[key]) < runs:
       xtlist = [[0]*numClasses]*(runs-len(shap_imp[key]))
       tmp = shap_imp[key] + xtlist 
       shap_imp[key] = tmp 
       tmp = None 
       xtlist = None 

   #print(key,shap_imp[key])    

   keyavg = []
   keysd = []
   for i in range(0,numClasses):
       vals = []
       for j in range(0,runs):
           #print(shap_imp[key][j][i])
           vals.append(shap_imp[key][j][i])

       avg = round(statistics.mean(vals),4)
       sd = round(statistics.stdev(vals),4)
       keyavg.append(avg)
       keysd.append(sd)
       avg = None 
       sd = None 
   ## print(keyavg,keysd)

   avg_imp.append(keyavg)
   sd_imp.append (keysd)



data = {
   'Features': shap_imp.keys(),
   'Mean_Imp': avg_imp, 
   'SD_Imp': sd_imp, 
}

feature_importance = pd.DataFrame(data)
feature_importance = feature_importance.sort_values(by='Mean_Imp', ascending=False)
print(feature_importance.head(20))

wrf.write('{}\n'.format(feature_importance))

feature_importance.to_csv('../results_MergedTNK/LUAD_TNK_FeatureImp_RF_SHAP.tsv',sep='\t',index=False)

avg_imp = None 
sd_imp = None 
shap_imp = None 
data = None 
feature_importance = None 

wrf.close()




# In[ ]:


## ======= Logistic regression ==========

lr_model=LogisticRegression(solver='saga',max_iter=10000)

#clf = RandomForestClassifier(class_weight="balanced")
lrm = LogisticRegression(solver='lbfgs', C=0.1, max_iter=1000)
rfecv = RFECV(estimator=lrm, step=1, cv=StratifiedKFold(5), scoring='roc_auc')

pipeline = Pipeline(
    steps=[
        ("scaling",scaler),
        ("norm",normalizer),
        ("variance", VarianceThreshold()),
        ("feature_sel1",SelectKBest(score_func=f_classif, k=200)),
        ("feature_sel2",SelectFromModel(LogisticRegression(solver='lbfgs', C=0.1, max_iter=1000))),
        ("feature_sel3",rfecv),
        ("lr_model", LogisticRegression(solver='saga',max_iter=10000))
    ]
)

param_grid = {
    "lr_model__C": [ 0.1, 0.5, 1], 
    "lr_model__l1_ratio": [ 0, 0.5, 1]
}

wrname = '../results_MergedTNK/Results_MLrun_Merged_LUAD_run_LogReg.txt'
wrf = open(wrname,'w')

acc_scores = []
acc_scores_bal = []
f1_scores = []
prec_scores = []
rec_scores = []
#roc_scores = []

nat_imp = {}
perm_imp = {}
shap_imp = {}


print('Sampling data ...')

for k in range(0,runs):
    print(k)
    wrf.write('Round: {}\n'.format(k))

    df1 = data_cleaned.loc[data_cleaned['status_encoded'] == 0]
    df2 = data_cleaned.loc[data_cleaned['status_encoded'] == 1]

    min_val = min(len(df1),len(df2))

    df1_sample = df1.sample(n=min_val)
    df1 = None 
    df2_sample = df2.sample(n=min_val)
    df2 = None 

    final_df = pd.concat([df1_sample, df2_sample], axis=0)
    df1_sample = None 
    df2_sample = None 


    print(final_df['status_encoded'].value_counts()) 

    xnames=final_df.columns[1:len(final_df.columns)-2]
    #print(xnames)
    X=final_df[xnames]
    y=final_df.status_encoded

    xnames = None 
    final_df = None  

    train_X, test_X, train_y, test_y = train_test_split(X, y, stratify=y, train_size=0.75)

    #sample_weights = compute_sample_weight(
    #    class_weight='balanced',
    #    y = train_y #provide your own target name
    #)

    #train_X_embedded, test_X_embedded = feature_selection1(train_X_scaled,train_y,test_X_scaled)
    #train_X_embedded, test_X_embedded = feature_selection2(train_X_scaled,train_y,test_X_scaled)

    grid_search = GridSearchCV(estimator=pipeline, param_grid=param_grid, cv=5, scoring='accuracy', verbose=1, n_jobs=-1)
    grid_search.fit(train_X, train_y) #,xgb_model__sample_weight=sample_weights) 

    best_pipeline = grid_search.best_estimator_
    best_params = grid_search.best_params_

    print(f"Best hyperparameters: {best_params}")
    wrf.write('  Best hyperparameters: {}\n'.format(best_params))
    #print(f"Best score (roc_auc): {grid_search.best_score_}")

    #pred_values=best_model.predict(test_X)
    pred_values = grid_search.predict(test_X)

    accr = accuracy_score(test_y,pred_values)
    accr_bal = balanced_accuracy_score(test_y,pred_values)
    f1s = f1_score(test_y,pred_values,average='weighted',zero_division = np.nan)
    prec = precision_score(test_y,pred_values,average='weighted',zero_division = np.nan)
    rec = recall_score(test_y,pred_values,average='weighted',zero_division = np.nan)
    #roc = roc_auc_score(test_y,pred_values,average='macro',multi_class='ovo')

    print(confusion_matrix(test_y, pred_values))
    wrf.write('  {}\n'.format(confusion_matrix(test_y, pred_values)))

    acc_scores.append(accr)
    acc_scores_bal.append(accr_bal)
    f1_scores.append(f1s)
    prec_scores.append(prec)
    rec_scores.append(rec)
    #roc_scores.append(roc)

    accr = None
    accr_bal = None 
    f1s = None 
    prec = None 
    rec = None 
    roc = None 
    #best_model = None 
    #best_params = None 


    ## === Feature importance === 

    best_model = best_pipeline.named_steps['lr_model']

    col_names = np.array(X.columns)
    mask1 = best_pipeline.named_steps['variance'].get_support()
    mask2 = best_pipeline.named_steps['feature_sel1'].get_support()
    mask3 = best_pipeline.named_steps['feature_sel2'].get_support()
    mask4 = best_pipeline.named_steps['feature_sel3'].get_support()
    features_mask1 = col_names[mask1]
    features_mask2 = features_mask1[mask2]
    features_mask3 = features_mask2[mask3]
    feature_names = features_mask3[mask4]
    ## print("Selected Features:", feature_names)

    ## Permutaiton method 
    print('\tPermutation importance')
    test_X_select = test_X[feature_names]
    result = permutation_importance(best_model, test_X_select.values, test_y, n_repeats=10, random_state=42)
    importance = result.importances_mean
    for lk,feature in enumerate(feature_names): 
        if feature in perm_imp.keys(): 
            tmp = perm_imp[feature]
            tmp.append(importance[lk])
            perm_imp[feature] = tmp
            tmp = None 
        else:
            perm_imp[feature] = [importance[lk]]    

    importance = None 
    lk = None 
    feature = None 
    result = None 

    X = None 
    y = None 
    best_model = None 
    best_params = None
    train_X = None 
    test_X = None 
    train_y = None 
    test_y = None 
    test_X_select = None 



print(f'average_accuracy: {statistics.mean(acc_scores):.2f}, sd_accuracy: {statistics.pstdev(acc_scores):.2f}\n')
print(f'average_accuracy_balanced: {statistics.mean(acc_scores_bal):.2f}, sd_accuracy: {statistics.pstdev(acc_scores_bal):.2f}\n')
print(f'average_f1-score: {statistics.mean(f1_scores):.2f}, sd_f1-score: {statistics.pstdev(f1_scores):.2f}\n')
print(f'average_precision: {statistics.mean(prec_scores):.2f}, sd_precision: {statistics.pstdev(prec_scores):.2f}\n')
print(f'average_recall: {statistics.mean(rec_scores):.2f}, sd_Recall: {statistics.stdev(rec_scores):.2f}\n')
#print(f'average_ROC-AUC_score: {statistics.mean(roc_scores):.2f}, sd_ROC-AUC_score: {statistics.stdev(roc_scores):.2f}\n')
wrf.write('\taverage_accuracy: {}, sd_accuracy: {}\n'.format(round(statistics.mean(acc_scores),2), round(statistics.pstdev(acc_scores),2)))
wrf.write('\taverage_accuracy_balanced: {}, sd_accuracy_balanced: {}\n'.format(round(statistics.mean(acc_scores_bal),2), round(statistics.pstdev(acc_scores_bal),2)))
wrf.write('\taverage_f1-score: {}, sd_f1-score: {}\n'.format(round(statistics.mean(f1_scores),2), round(statistics.pstdev(f1_scores),2)))
wrf.write('\taverage_precision: {}, sd_precision: {}\n'.format(round(statistics.mean(prec_scores),2), round(statistics.pstdev(prec_scores),2)))
wrf.write('\taverage_recall: {}, sd_Recall: {}\n'.format(round(statistics.mean(rec_scores),2), round(statistics.pstdev(rec_scores),2)))

wrf.write("\n====== Top Features by Permutation Importance ======= \n")
print("\nTop Features by Permutation Importance:")

avg_imp = []
sd_imp = []

for key in perm_imp.keys():
    #print(key,nat_imp[key])
    if len(perm_imp[key]) < runs:
        xtlist = [0]*(runs-len(perm_imp[key]))
        tmp = perm_imp[key] + xtlist 
        perm_imp[key] = tmp 
        tmp = None 
        xtlist = None 

    #print(key,nat_imp[key])    
    avg = round(statistics.mean(perm_imp[key]),4)
    series = pd.Series(perm_imp[key]) 
    sd = round(series.std(ddof=0),4)
    series = None 

    avg_imp.append(avg)
    sd_imp.append(sd)

data = {
    'Features': perm_imp.keys(),
    'Mean_Imp': avg_imp, 
    'SD_Imp': sd_imp, 
}

feature_importance = pd.DataFrame(data)
feature_importance = feature_importance.sort_values(by='Mean_Imp', ascending=False)
print(feature_importance.head(20))
wrf.write('{}\n'.format(feature_importance))

feature_importance.to_csv('../results_MergedTNK/LUAD_TNK_FeatureImp_LogReg_Perm.tsv',sep='\t',index=False)

avg_imp = None 
sd_imp = None 
perm_imp = None 
data = None 
feature_importance = None 

wrf.close()


# In[ ]:


data_cleaned = None 


# In[ ]:




