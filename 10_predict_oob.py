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
import matplotlib
import matplotlib.pyplot as plt

from pylab import *

inF = "./bindingdata.ready2model.sdf.gz.2class_hERGTL.CUT_10000nM.sdf.gz.removedSDtags.sdf.sdf.gz.descr.sdf.sdf.gz"
#inF = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

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

# Den Zaehler benoetige ich, wenn ich mehrere train/test splits durchfuehre 
randomseedcounter = 1
X_train,X_test,y_train,y_test = train_test_split(\
dataDescrs_array,dataActs_array,test_size=.4,random_state=randomseedcounter)
##
# RF
clf_RF     = RandomForestClassifier(compute_importances=True,random_state=0,oob_score=True,n_estimators=100)
clf_RF     = clf_RF.fit(X_train,y_train)
#    print "oob_decision_function_", len(clf_RF.oob_decision_function_)
<
nbins = 10
binhits = [0 for i in range(nbins)]
bin_actives = [0 for i in range(nbins)]
for score,cl in zip(clf_RF.oob_decision_function_,y_train):
  bin_idx = int(nbins * score[1]-0.000000001)
  binhits[bin_idx] += 1
  if cl == 1: bin_actives[bin_idx] += 1

y_frac_active = [float(i)/j for i,j in zip(bin_actives,binhits)]

y_frac_inactive = [float(j-i)/j for i,j in zip(bin_actives,binhits)]

x_active   = [float(i)/10 for i in range(10)]
x_inactive = [float(i)/10 for i in range(9,-1,-1)] #range(9,-1,-1)

print "train set: %i inactives || %i actives \n" %(np.bincount(y_train)[0],np.bincount(y_train)[1])

print "bin_actives %s  %s" %(bin_actives,len(bin_actives))
print "bin_hits    %s  %s" %(binhits,len(binhits))
print "x_inactive ", x_inactive
print "x_active   ", x_active
#
#
print "\n\n"
y_frac_active_fake   = [float(i) for i,j in zip(bin_actives,binhits)]
y_frac_inactive_fake = [float(j-i) for i,j in zip(bin_actives,binhits)]
print "y_frac_active_fake   ", y_frac_active_fake
print "y_frac_inactive_fake ", y_frac_inactive_fake

fig = plt.figure(figsize=(10,5))
#http://matplotlib.1069221.n5.nabble.com/title-when-using-subplot-td22151.html
fig.text(.5, .95, 'train set - out of bag probabilties', \
         horizontalalignment='center', fontsize=14)
ax2 = fig.add_subplot(121)
ax2.set_xlim([0,1])
ax2.set_ylim([0,1])
plt.scatter(x_inactive,y_frac_inactive,marker="o",s=40)
ax2.set_xlabel("relative inactives count   [x_inactive]",fontsize=10)
ax2.set_ylabel("fraction of true_inactives [y_frac_inactive]",fontsize=10)
ax1 = fig.add_subplot(122)
ax1.set_xlim([0,1])
ax1.set_ylim([0,1])
plt.scatter(x_active,y_frac_active,marker="o",s=40)
ax1.set_xlabel("relative actives count   [x_active]",fontsize=10)
ax1.set_ylabel("fraction of true_actives [y_frac_active]",fontsize=10)
plt.show()
