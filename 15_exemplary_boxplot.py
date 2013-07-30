from rdkit import Chem
import cPickle,gzip
import numpy as np
from matplotlib import pyplot
from pylab import *
from scipy import stats

filename_inactives = sys.argv[1]
filename_actives = sys.argv[2]
cpds_actives = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(filename_actives)) if x is not None]
cpds_inactives = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(filename_inactives)) if x is not None]

hERG_TL_list_actives = []
property_list_list_actives = []
cpd_name_list_actives = []
prop_dict_actives_all = {}
for cpd in cpds_actives:
    property_list_actives = []
    property_name_list_actives = []
    prop_name_actives = cpd.GetPropNames()
    cpd_name_list_actives.append(cpd.GetProp("_Name"))
    for property in prop_name_actives:
            if property not in prop_dict_actives_all:
                prop_dict_actives_all[property] = []
                prop_dict_actives_all[property].append(float(cpd.GetProp(property)))
            else:
                prop_dict_actives_all[property].append(float(cpd.GetProp(property)))
                property_name_list_actives.append(property)
#
hERG_TL_list_inactives = []
property_list_list_inactives = []
cpd_name_list_inactives = []
prop_dict_inactives_all = {}
for cpd in cpds_inactives:
    property_list_inactives = []
    property_name_list_inactives = []
    prop_name_inactives = cpd.GetPropNames()
    cpd_name_list_inactives.append(cpd.GetProp("_Name"))
    for property in prop_name_inactives:
            if property not in prop_dict_inactives_all:
                list_element = cpd.GetProp(property)
                prop_dict_inactives_all[property] = []
                prop_dict_inactives_all[property].append(float(cpd.GetProp(property)))
            else:
                prop_dict_inactives_all[property].append(float(cpd.GetProp(property)))
                property_name_list_inactives.append(property)
intersectionlist = ['NumHAcceptors']
counter = 0
data = []
label_list = []
for relevant_descr in intersectionlist:
    data.append([prop_dict_actives_all[relevant_descr],prop_dict_inactives_all[relevant_descr]])
    label_list.append(["actives","inactives"])
    counter += 1
    
for ii,x in enumerate(data):
    print intersectionlist[ii]
    medianlist = []
    for i,value in enumerate(x):
        print label_list[ii][i]
        print "%1.1f // %1.1f // %1.1f" %(round(stats.scoreatpercentile(x[i],25),1),round(median(x[i]),1),round(stats.scoreatpercentile(x[i],75),1))


fig = plt.figure(figsize=(10,10))

for x,i in enumerate(intersectionlist):
    plotname = "ax"+str(x)
    subplotnumbering = int("11"+str(x+1))
    plotname = fig.add_subplot(subplotnumbering)
    xtickNames = plt.setp(plotname, xticklabels=label_list[x])
    plt.setp(xtickNames, rotation=45, fontsize=24)
    bp = plt.boxplot(data[x], notch=0, sym='+', vert=1, whis=1.5)
    plotname.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.5)
    plotname.set_title(intersectionlist[x],fontsize=24)
    plotname.tick_params(axis='both', which='major', labelsize=20)

show()
