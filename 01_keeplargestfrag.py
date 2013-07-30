import gzip,sys,os
from rdkit import Chem
#zims = [x for x in Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]


#file_name = sys.argv[1]+".onlylargestfrag.sdf.gz"
#test_output = gzip.open(file_name,'w+')
#test_cpd_out = Chem.SDWriter(test_output)
#test_cpd_out = Chem.ForwardSDMolSupplier(test_output)

file_name    = sys.argv[1]+".onlylargestfrag.sdf"
test_cpd_out = Chem.SDWriter(file_name)

for cpd in cpds:
    fragmented = Chem.GetMolFrags(cpd,asMols=True)
    list_cpds_fragsize = []
    for x in fragmented:
        list_cpds_fragsize.append(x.GetNumAtoms())
    largest_frag_index = list_cpds_fragsize.index(max(list_cpds_fragsize))
    largest_frag = fragmented[largest_frag_index]
    test_cpd_out.write(largest_frag)
test_cpd_out.close()

#http://stackoverflow.com/questions/8156707/gzip-a-file-in-python
f_in = open(file_name, 'rb')
file_name_gz = file_name+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
os.remove(file_name)
#f_out.flush()


#http://stackoverflow.com/questions/8156707/gzip-a-file-in-python
# another method
##file_name_gz = file_name+".sdf.gz"
##with open(file_name, 'rb') as orig_file:
##    with gzip.open(file_name_gz, 'wb') as zipped_file:
##        zipped_file.writelines(orig_file)
##    zipped_file.close()
##orig_file.close()

#f_in.flush()


#test_cpd_out.flush()
#test_output.flush()
#test_cpd_out=None
#test_output=None
