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

all_struct_dict = {}
for cpd in cpds:
    cansmi = cpd.GetProp('cansmirdkit')
    try:
        all_struct_dict[cansmi].append(cpd)
    except:
        all_struct_dict[cansmi] = [cpd]

all_struct_list = []
for ent in all_struct_dict:
    for blub in all_struct_dict[ent]:
        all_struct_list.append(len(all_struct_dict[ent]))
        

fig = plt.figure(figsize=(10,8))
ax1 = fig.add_subplot(111)
plt.xticks(range(1,20))
plt.hist(all_struct_list,bins=20)

ax1.set_title("measurements per compound",fontsize=32)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlabel("compound count",fontsize=24)
ax1.set_ylabel("measurement count",fontsize=24)
show()
