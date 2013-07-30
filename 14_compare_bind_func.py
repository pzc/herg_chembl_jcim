from rdkit import Chem
import sys,gzip,os



bind_inF  = sys.argv[1]
func_inF  = sys.argv[2]
outF      = Chem.SDWriter("merged.sdf")

if bind_inF.endswith(".sdf.gz") or bind_inF.endswith(".sd.gz"):
    bind_cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(bind_inF)) if x is not None]
else:
    bind_cpds = [x for x in  Chem.SDMolSupplier(bind_inF) if x is not None]

if func_inF.endswith(".sdf.gz") or func_inF.endswith(".sd.gz"):
    func_cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(func_inF)) if x is not None]
else:
    func_cpds = [x for x in  Chem.SDMolSupplier(func_inF) if x is not None]


value_bind_dict = {}
for bind_cpd in bind_cpds:
    cansmi = bind_cpd.GetProp("cansmirdkit")
    value_bind = bind_cpd.GetProp("value")
    bind_cpd.SetProp("value_bind",value_bind)
    bind_cpd.ClearProp("value")
    hERG_TL_bind = bind_cpd.GetProp("hERG_TL")
    bind_cpd.SetProp("hERG_TL_bind",hERG_TL_bind)
    bind_cpd.ClearProp("hERG_TL")
    try:
        value_avg_bind = bind_cpd.GetProp("value_avg")
        bind_cpd.SetProp("value_avg_bind",value_avg_bind)
        bind_cpd.ClearProp("value_avg")
    except:
        pass
    try:
        value_stddev_bind = bind_cpd.GetProp("value_stddev")
        bind_cpd.SetProp("value_stddev_bind",value_stddev_bind)
        bind_cpd.ClearProp("value_stddev")
    except:
        pass

    value_bind_dict[cansmi] = bind_cpd
    
value_func_dict = {}
for func_cpd in func_cpds:
    cansmi = func_cpd.GetProp("cansmirdkit")
    if value_bind_dict.has_key(cansmi):
        value_func_dict[cansmi] = value_bind_dict[cansmi]
        value_func = func_cpd.GetProp("value")

        hERG_TL_func = func_cpd.GetProp("hERG_TL")
        value_bind_dict[cansmi].SetProp("hERG_TL_func",hERG_TL_func)
        value_bind_dict[cansmi].SetProp("value_func",value_func)
        
        #func_cpd.SetProp("value_func",value_func)
        #func_cpd.ClearProp("value")
        try:
            value_avg_func = bind_cpd.GetProp("value_avg")
            func_cpd.SetProp("value_avg_func",value_avg_func)
            func_cpd.ClearProp("value_avg")
        except:
            pass
        try:
            value_stddev_func = func_cpd.GetProp("value_stddev")
            func_cpd.SetProp("value_stddev_func",value_stddev_func)
            func_cpd.ClearProp("value_stddev")
        except:
            pass

print len(value_func_dict)

for x in value_func_dict:
    outF.write(value_func_dict[x])
outF.close()
