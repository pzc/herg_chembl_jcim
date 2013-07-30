from rdkit import Chem
import cPickle,gzip
from pylab import *

# CAVE
# run this on uniq.func.combined/func.uniq.dupl.sdf.gz

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]


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

intersectionlist = ['fr_COO']


counter = 0
data = []
label_list = []
for relevant_descr in intersectionlist:
    data.append([prop_dict_actives_all[relevant_descr],prop_dict_inactives_all[relevant_descr]])
    label_list.append(["actives","inactives"])
    counter += 1
fig = plt.figure(figsize=(10,10))

master_dict_list = [prop_dict_inactives_all,prop_dict_actives_all]

title_list = ["inactives","actives"]

for x,i in enumerate(master_dict_list):
    plotname = "ax"+str(x)
    subplotnumbering = int("12"+str(x+1))
    plotname = fig.add_subplot(subplotnumbering)
    plotname.hist(master_dict_list[x]['fr_COO'])
    plotname.tick_params(axis='both', which='major', labelsize=14)
    plotname.set_title(title_list[x],fontsize=30)
    plotname.set_xlabel("COO count",fontsize=20)
    plotname.set_ylabel("compound count",fontsize=20)
show()
