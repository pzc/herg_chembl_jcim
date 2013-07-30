from rdkit import Chem
from scipy import interp
import numpy as np
from sklearn import cross_validation
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.cross_validation import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_score,recall_score
from sklearn import preprocessing
import os,sys,csv,cPickle,gzip
from time import time
t0 = time()
inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]
csv_file        = open(sys.argv[1]+".CV.100est.csv","wb")
csv_file_writer = csv.writer(csv_file,delimiter=";",quotechar=' ')
title_line = ["accuracy","MCC","precision","recall","f1","auc","assaytype","threshold","randomseed"]
csv_file_writer.writerow(title_line)
directory = os.getcwd().split("/")[-2:]
dir_string  = ';'.join(directory)
#
hERG_TL_list = []
property_list_list = []
for cpd in cpds:
    property_list = []
    property_name_list = []
    prop_name = cpd.GetPropNames()
    for property in prop_name:
        if property not in ['hERG_TL','value']:
            property_list.append(float(cpd.GetProp(property)))
            property_name_list.append(property)
        elif property == 'hERG_TL':
            hERG_TL_list.append(int(cpd.GetProp(property)))
        else:
            pass
    property_list_list.append(property_list)
dataDescrs_array = np.asarray(property_list_list)
dataActs_array   = np.array(hERG_TL_list)

for randomseedcounter in range(1,11):
    print "seed", randomseedcounter
    X_train,X_test,y_train,y_test = cross_validation.train_test_split(dataDescrs_array,dataActs_array,test_size=.4,random_state=randomseedcounter)
    ##
    # RF
    clf_RF     = RandomForestClassifier(compute_importances=True,n_estimators=100,random_state=randomseedcounter)
    clf_RF     = clf_RF.fit(X_train,y_train)
    cv_counter = 5
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.zero_one_score)
    accuracy_CV = round(scores.mean(),3)
    accuracy_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.matthews_corrcoef)
    MCC_CV = round(scores.mean(),3)
    MCC_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.f1_score)
    scores_rounded = [round(x,3) for x in scores]
    f1_CV = round(scores.mean(),3)
    f1_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.precision_score)
    scores_rounded = [round(x,3) for x in scores]
    precision_CV = round(scores.mean(),3)
    precision_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.recall_score)
    scores_rounded = [round(x,3) for x in scores]
    recall_CV = round(scores.mean(),3)
    recall_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.auc_score)
    scores_rounded = [round(x,3) for x in scores]
    auc_CV = round(scores.mean(),3)
    auc_std_CV = round(scores.std(),3)
    result_string_cut = [str(accuracy_CV)+"_"+str(accuracy_std_CV),
                         str(MCC_CV)+"_"+str(MCC_std_CV),
                         str(precision_CV)+"_"+str(precision_std_CV),
                         str(recall_CV)+"_"+str(recall_std_CV),
                         str(f1_CV)+"_"+str(f1_std_CV),
                         str(auc_CV)+"_"+str(auc_std_CV),
                         dir_string,randomseedcounter]
    csv_file_writer.writerow(result_string_cut)
t1 = time()-t0
print "done in %1.2f seconds" %(t1)
