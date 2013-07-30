import sys,gzip,os
from rdkit import Chem

inF       = sys.argv[1]
SD_tag    = sys.argv[2]
cutoff    = sys.argv[3]
file01_wr = inF+".2class_hERGTL.CUT_%snM.sdf" %(cutoff)
ms_wr01   = Chem.SDWriter(file01_wr)

def main(argv=[__name__]):
    if len(argv) != 4:
        print (" \n\nUSAGE: %s <SD file> <IC50 flag> <cutoff value in nM>\n" % argv[0])
        sys.exit()
    if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
        cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
    else:
        cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]
    for mol in cpds:
        # highly active cpds -> class1 (smaller than 1 mikroMol)
        if float(mol.GetProp(SD_tag))> float(float(cutoff)):
            mol.SetProp('hERG_TL','0')
        else:
            mol.SetProp('hERG_TL','1')
        ms_wr01.write(mol)

    ms_wr01.close()
    f_in = open(file01_wr, 'rb')
    file_name_gz = file01_wr+".gz"
    f_out = gzip.open(file_name_gz, 'wb')
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()
    os.remove(file01_wr)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
