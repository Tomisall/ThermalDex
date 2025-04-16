#from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Mol, MolFromSmiles, MolFromSmarts, rdmolfiles, rdDetermineBonds, GetBonds
from rdkit import Chem

mol = Chem.MolFromSmiles("CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C")
#print(mol)
#test = rdDetermineBonds.DetermineBonds(mol)
#print(test)

def get_atom_types(molecule):
    atom_types = [atom.GetAtomicNum() for atom in molecule.GetAtoms()]
    return atom_types

def get_bond_info(molecule):
    bonds = []
    for bond in molecule.GetBonds():
        atom_types = get_atom_types(molecule)
        bond_type = bond.GetBondTypeAsDouble()
        bonds.append({
            'atom1_idx': atom_types[bond.GetBeginAtomIdx()],
            'atom2_idx': atom_types[bond.GetEndAtomIdx()],
            'bond_type': bond_type
        })
    return bonds

result = get_bond_info(mol)
print(result)