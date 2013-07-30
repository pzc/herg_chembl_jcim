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


hERG_TL_list = []
property_list_list = []
cpd_name_list = []
for cpd in cpds:
    property_list = []
    property_name_list = []
    prop_name = cpd.GetPropNames()
    cpd_name_list.append(cpd.GetProp("_Name"))
    for property in prop_name:
        if property in ['value']:
            pIC50 = -log10(float(cpd.GetProp(property))*1e-9)
        else:
            pass
    property_list_list.append(pIC50)

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(111)
plt.hist(property_list_list)
#ax1.set_title("pIC50 distribution",fontsize=32)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlabel("pIC50",fontsize=24)
ax1.set_ylabel("compound count",fontsize=24)
show()
