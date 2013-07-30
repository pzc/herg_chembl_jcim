from rdkit import Chem
from rdkit.Chem import AllChem
import sys,numpy,gzip,os

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

file_out = inF+".avg.sdf"
cpds_out = Chem.SDWriter(file_out)

IC50_dict = {}
for cpd in cpds:
    cansmi = str(cpd.GetProp("cansmirdkit"))
    IC50_dict[cansmi]={}
    
def get_mean_IC50(mol_list):
    IC50 = 0
    IC50_avg = 0
    for bla in mol_list:
        try:
            IC50 +=  float(bla.GetProp("value"))
        except:
            print "no IC50 reported",bla.GetProp("_Name")
    IC50_avg = IC50 / len(mol_list)
    return IC50_avg

def get_stddev_IC50(mol_list):
    IC50_list = []
    for bla in mol_list:
        try:
            IC50_list.append(round(float(bla.GetProp("value")),2))
        except:
            print "no IC50 reported",bla.GetProp("_Name")
    IC50_stddev = numpy.std(IC50_list,ddof=1)
    # http://stackoverflow.com/questions/4575645/different-standard-deviation-for-same-input-from-wolfram-and-numpy
    return IC50_stddev,IC50_list
    
for cpd in cpds:
    cansmi = str(cpd.GetProp("cansmirdkit"))
    try:
        IC50_dict[cansmi].append(cpd)
    except:
        IC50_dict[cansmi] = [cpd]
        
for ent in IC50_dict:
    #print ent,IC50_dict[ent]
    #for x in IC50_dict[ent]:
    #    print x.GetProp("_Name")
    IC50_avg = str(get_mean_IC50(IC50_dict[ent]))
    #IC50_dict[ent][0].SetProp("value_avg",IC50_avg)
    IC50_stddev = get_stddev_IC50(IC50_dict[ent])[0]
    IC50_dict[ent][0].SetProp("value_stddev",str(IC50_stddev))
    IC50_list   = get_stddev_IC50(IC50_dict[ent])[1]
    IC50_dict[ent][0].SetProp("value_avg",str(IC50_list))
    minimumvalue = float(IC50_avg)-3*float(IC50_stddev)
    maximumvalue = float(IC50_avg)+3*float(IC50_stddev)
    IC50_dict[ent][0].SetProp("value",IC50_avg)
    if round(IC50_stddev,1) == 0.0:
        cpds_out.write(IC50_dict[ent][0])
        #print "perfect!", IC50_list, IC50_avg,IC50_stddev
    elif IC50_stddev > float(IC50_avg):
        runawaylist =[]
        for x in IC50_dict[ent]:
            runawaylist.append(x.GetProp("_Name"))
        print "stddev larger than mean", runawaylist, IC50_list, IC50_avg,IC50_stddev
        #print "stddev larger than mean", IC50_list, IC50_avg,IC50_stddev    
    elif numpy.min(IC50_list) < minimumvalue or numpy.max(IC50_list) > maximumvalue:
        pass
        #print "LARGER than 3fold_sigma", IC50_list, IC50_avg,IC50_stddev,minimumvalue,maximumvalue
    else:
        pass
        #print "deviation from mean in reasonable range. " ##,IC50_list, IC50_avg,IC50_stddev,minimumvalue,maximumvalue
        cpds_out.write(IC50_dict[ent][0])

cpds_out.close()

f_in = open(file_out, 'rb')
file_name_gz = file_out+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
f_out.close()
f_in.close()
os.remove(file_out)
