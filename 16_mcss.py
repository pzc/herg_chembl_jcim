import sys
import re
from rdkit import Chem
import hacked_fmcs
from time import time

import argparse
parser = argparse.ArgumentParser("Find fragments")
parser.add_argument("filename",
                    help="SMILES filename")
parser.add_argument("--threshold", "-t", type=float, default=0.5,
                    help="minimum match threshold (default: 0.5)")
parser.add_argument("--min-num-atoms", type=int, default=2,
                    help="minimum number of atom for the match (default: 2)")
parser.add_argument("--min-num-bonds", type=int, default=1,
                    help="minimum number of bonds for the match (default: 1)")
parser.add_argument("--complete-rings-only", action="store_true",
                    help="don't match part ring in a fragment")

args = parser.parse_args()
if args.min_num_atoms < 2:
    parser.error("--min-num-atoms must be at least 2")
if args.min_num_bonds < 1:
    parser.error("--min-num-bonds must be at least 2")
if args.threshold <= 0.0:
    parser.error("--threshold must be positive")

table = Chem.GetPeriodicTable()
_atomic_num_to_symbol = dict((str(i), table.GetElementSymbol(i)) for i in range(105))
atom_pat = re.compile("\[#(\d+)\]")
def atomic_num_to_symbol(m):
    return _atomic_num_to_symbol[m.group(1)]
    
def smarts_to_unique_smiles_fragment(smarts):
    fragment_smiles = atom_pat.sub(atomic_num_to_symbol, smarts)
    fragment_smiles = fragment_smiles.replace("!@", "").replace("@", "")
    mol = Chem.MolFromSmiles(fragment_smiles, sanitize=False)
    # Correct the aromaticity information
    for bond in mol.GetBonds():
        if bond.GetIsAromatic():
            bond.GetBeginAtom().SetIsAromatic(1)
            bond.GetEndAtom().SetIsAromatic(1)
    return Chem.MolToSmiles(mol)

## Hack in code to check the results of the complete_rings_only check

only_has_complete_rings = {}
old_check_complete_rings_only = hacked_fmcs.check_complete_rings_only

def new_check_complete_rings_only(smarts, subgraph, enumeration_mol):
    result = old_check_complete_rings_only(smarts, subgraph, enumeration_mol)
    only_has_complete_rings[smarts] = result
    return result

hacked_fmcs.check_complete_rings_only = new_check_complete_rings_only

hacked_fmcs._maximize_options[("all-matches", False)] = (
    hacked_fmcs.no_pruning, hacked_fmcs.SingleBestAtoms)
hacked_fmcs._maximize_options[("all-matches", True)] = (
    hacked_fmcs.no_pruning, hacked_fmcs.SingleBestAtomsCompleteRingsOnly)

##
# PC
##
t0 = time()

##### Read in the input file

mols = []
for lineno, line in enumerate(open(args.filename)):
    fields = line.split()
    if not fields:
        sys.stderr.write("Missing SMILES on line %d" % (lineno+1,))
        continue
    smiles = fields[0]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.stderr.write("Could not process SMILES %s on line %d\n" % (smiles, lineno+1))
        continue
    mols.append(mol)

mcs_result = hacked_fmcs.fmcs(mols,
                              threshold=args.threshold,
                              min_num_atoms=args.min_num_atoms,
                              complete_rings_only=args.complete_rings_only,
                              # fmcs doesn't implement min_num_bonds; use post-processing
                              maximize="all-matches")

# Convert to unique fragment SMILES
matches = []
unique_fragment_smiles = set()
for smarts, is_match in mcs_result.matches.items():
    if not is_match:  # remember, this is hacked code, and a bit ugly!
        continue
    if args.complete_rings_only:
        if not only_has_complete_rings.get(smarts, False):
            continue

    fragment_smiles = smarts_to_unique_smiles_fragment(smarts)
    if fragment_smiles in unique_fragment_smiles:
        continue
    unique_fragment_smiles.add(fragment_smiles)
    pattern = Chem.MolFromSmarts(smarts)
    if pattern.GetNumBonds() < args.min_num_bonds:
        continue
    matches.append( (pattern, smarts, fragment_smiles) )

def order_by_number_of_atoms(term):
    return term[0].GetNumAtoms()
    
matches.sort(key=order_by_number_of_atoms)


print "frequency #atoms #bonds SMILES SMARTS"
for pattern, smarts, fragment_smiles in matches:
    count = sum(1 for mol in mols if mol.HasSubstructMatch(pattern))
    print count, pattern.GetNumAtoms(), pattern.GetNumBonds(), fragment_smiles, smarts

print "done in %0.3fs" % (time() - t0)
