import gzip,sys,os
from rdkit import Chem

inF1  = sys.argv[1]
inF2  = sys.argv[2]

if inF1.endswith(".sdf.gz") or inF1.endswith(".sd.gz"):
    cpds1 = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
else:
    cpds1 = [x for x in  Chem.SDMolSupplier(inF1) if x is not None]
if inF2.endswith(".sdf.gz") or inF2.endswith(".sd.gz"):
    cpds2 = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[2])) if x is not None]
else:
    cpds2 = [x for x in  Chem.SDMolSupplier(inF2) if x is not None]


file_out = inF1.split(".")[0]+".ready2model"
test_cpd_out = Chem.SDWriter(file_out)

for cpd1 in cpds1:
    test_cpd_out.write(cpd1)
for cpd2 in cpds2:
    test_cpd_out.write(cpd2)

test_cpd_out.close()

f_in = open(file_out, 'rb')
file_name_gz = file_out+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
f_out.close()
f_in.close()
os.remove(file_out)
