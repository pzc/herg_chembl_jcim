import sys,gzip,os
from rdkit import Chem

inF            = sys.argv[1]
SDFilename_out = inF+'.removedSDtags.sdf'
file_SDtags =  sys.argv[2]

def main(argv=[__name__]):
    if len(argv) != 3:
        print (" \n\nUSAGE: %s <SD file> <TXT file with SDtags2remove, line-wise!>\n" % argv[0])
    if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
        SDFile = Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1]))
    else:
        SDFile = Chem.SDMolSupplier(inF)    
    SDTags2remove = open(sys.argv[2],'r').readlines()
    SDFileout = Chem.SDWriter(SDFilename_out)
    SDTags2remove_rstripped = []
    for ent in SDTags2remove:
        SDTags2remove_rstripped.append(ent.rstrip())
    for mol in SDFile:
        for tags2remove in SDTags2remove_rstripped:
            if mol.HasProp(tags2remove):
                mol.ClearProp(tags2remove)
        SDFileout.write(mol)
    SDFileout.close()
    f_in = open(SDFilename_out, 'rb')
    file_name_gz = SDFilename_out+".sdf.gz"
    f_out = gzip.open(file_name_gz, 'wb')
    f_out.writelines(f_in)
    os.remove(SDFilename_out)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
