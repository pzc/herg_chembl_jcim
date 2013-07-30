from rdkit import Chem
import sys,gzip,os

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    #test_output = gzip.open(inF,'w+')
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

uniqF = inF+".uniq.sdf"
duplF = inF+".dupl.sdf"
uniqF_SD = Chem.SDWriter(uniqF)
duplF_SD = Chem.SDWriter(duplF)

# key: canonical_smiles == value: mol_object, occurence
all_struct_dict = {}
all_struct_list = []
for cpd in cpds:
    Chem.RemoveHs(cpd)
    cansmi = Chem.MolToSmiles(cpd,canonical=True)
    try:
        all_struct_dict[cansmi].append(cpd)
    except:
        all_struct_dict[cansmi] = [cpd]
    all_struct_list.append(cansmi)
    
for ent in all_struct_dict:
    if len(all_struct_dict[ent]) == 1:
        can_smi = ent
        all_struct_dict[ent][0].SetProp('cansmirdkit',can_smi)
        uniqF_SD.write(all_struct_dict[ent][0])
    else:
        for singleentries in all_struct_dict[ent]:
            can_smi = ent
            singleentries.SetProp('cansmirdkit',can_smi)
            duplF_SD.write(singleentries)
uniqF_SD.close()
duplF_SD.close()

# uniq-gz            
f_in = open(uniqF, 'rb')
file_name_gz = uniqF+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
f_out.close()
f_in.close()
os.remove(uniqF)

# dupl-gz
f_in = open(duplF, 'rb')
file_name_gz = duplF+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
f_out.close()
f_in.close()
os.remove(duplF)
