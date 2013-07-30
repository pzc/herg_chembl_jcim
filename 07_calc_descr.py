import sys,gzip,os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

outF     = inF+".descr.sdf"
cpds_out = Chem.SDWriter(outF)
nms=[x[0] for x in Descriptors._descList]
nms.remove('MolecularFormula')
print len(Descriptors._descList)

calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
for i in range(len(cpds)):
    descrs = calc.CalcDescriptors(cpds[i])
    for x in range(len(descrs)):
        cpds[i].SetProp(str(nms[x]),str(descrs[x]))
    cpds_out.write(cpds[i])
cpds_out.close()
f_in = open(outF, 'rb')
file_name_gz = outF+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
os.remove(outF)
