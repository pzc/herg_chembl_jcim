from rdkit import Chem
import cPickle,gzip
from pylab import *

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

all_struct_dict = {}
all_struct_list = []
for cpd in cpds:
    chemblassayid = cpd.GetProp("assay__chemblid")
    try:
        all_struct_dict[chemblassayid].append(cpd)
    except:
        all_struct_dict[chemblassayid] = [cpd]
    all_struct_list.append(chemblassayid)

new_dict = {}
new_list = []
for ent in all_struct_dict:
    new_dict[ent] = len(all_struct_dict[ent])
    new_list.append(len(all_struct_dict[ent]))

'''
# binding assay
fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(111)
plt.hist(new_list,bins=20)
ax1.set_title("measurements \n per CHEMBL assay ID",fontsize=32)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlabel("number of measurements",fontsize=24)
ax1.set_ylabel("number of compounds",fontsize=24)
plt.show()

'''
    
'''
# functional assay
fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(121)
ax1.set_title("measurements",fontsize=32)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlabel("number of",fontsize=24)
ax1.set_ylabel("number of compounds",fontsize=24)
bins1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
#lt.yticks([30,40])
#plt.yticks([5,10,15,20,25,30,35,40])
plt.hist(new_list,bins1)
ax2 = fig.add_subplot(122)
ax2.set_title("per CHEMBL assay ID",fontsize=32)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.set_xlabel("measurements",fontsize=24)
bins2 = [600,610,620,630,640,650,660,670,680,690,700]
plt.hist(new_list,bins2)
'''
# functional assay
fig = plt.figure(figsize=(10,10))

###ax1 = fig.add_subplot(111)
###bins1=[1,2,3,4,5,6,7,8,9,10]
###ax1.hist(new_list,bins1,rwidth=0.05,align='left')


ax2 = fig.add_subplot(111)
bins2=[673,675]
ax2.hist(new_list,bins2,rwidth=0.05)

plt.show()

print new_list
