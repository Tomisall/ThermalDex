from dataclasses import dataclass, field
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import Mol
import csv

@dataclass
class PSMolecule:
	smiles: str = None
	name: str = None
	HEG: int = None
	mp: float = None
	MW: float = None
	TE: float = None
	proj: str = None

	def toIterable(self):
		return iter(
			[
				self.smiles,
				self.name,
				self.HEG,
				self.mp,
				self.MW,
				self.TE,
				self.proj
			]
	)

	def toHeader(self):
		return [
			"Molecule",
			"Name",
			"Number of High Energy Groups",
			"Melting point",
			"Molecular Weight",
			"Thermal Event",
			"Project",
			]





entryOne = PSMolecule(smiles="O[C@@H]([C@H](O)C(O)=O)C(O)=O |r|", name="TestMol", HEG=1, mp=136.7, TE=242.8, proj="testProject") #"N1N=NN=C1C1=CC=CC=C1 |c:1,3,8,10,t:6|", name="TestMol", HEG=1, mp=136.7, TE=242.8, proj="testProject")
RDMol = Chem.MolFromSmiles(entryOne.smiles)
#cmpdMW = Descriptors.ExactMolWt(RDMol) # gives exact weight
cmpdMW = Descriptors.MolWt(RDMol)

print(RDMol)

#entryOne = PSMolecule(smiles="N1N=NN=C1C1=CC=CC=C1 |c:1,3,8,10,t:6|", HEG=1, mp=136.7, MW=cmpdMW, TE=242.8)
print(entryOne)

listforDB = [entryOne]

def writetoCSV(PSDBList: list):
	with open("PSdatabase.csv", "a") as DB:
		writer = csv.writer(DB)
		writer.writerow(PSMolecule().toHeader())
		for item in PSDBList:
			writer.writerow(item.toIterable())

writetoCSV(listforDB)

testSubstructure = Chem.MolFromSmiles("N1C=NN=N1 |c:1,3|")
matchList = Mol.GetSubstructMatches(RDMol, testSubstructure)

print(matchList)
print(len(matchList))


fullMatch = 0
with open("HighEnergyGroups.csv", "r") as HEGroups:
	for line in HEGroups:
		print(line)
		HeSubstructure = Chem.MolFromSmiles(line)
		fullmatchList = Mol.GetSubstructMatches(RDMol, HeSubstructure)
		fullMatch += len(fullmatchList)

print(fullMatch)

#Draw.ShowMol(RDMol)
#Draw.ShowMol(testSubstructure)
