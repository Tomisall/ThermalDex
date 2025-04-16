import pandas as pd
import os  
from dataclasses import asdict
from thermDex.thermDexMolecule import thermalDexMolecule


def cleanMolDataFrame(molecule):
    dictMol = asdict(molecule)
    for key in dictMol:
        dictMol[key] = str(dictMol[key])
    print(dictMol)
    molData = pd.DataFrame.from_dict([dictMol])
    print(molData)
    selectedMolData = molData.drop(columns=['mol','molIMG','molPixmap','mwStr','obStr','isStr','epStr','eleList','Catoms','Hatoms','Oatoms'])
    selectedMolData.index = selectedMolData['SMILES'].values
    print(selectedMolData)
    return selectedMolData

def createDatabase(molecule):
    dictMol = asdict(molecule)
    print(dictMol)
    for key in dictMol:
        dictMol[key] = str(dictMol[key])
    print(dictMol)
    molData = pd.DataFrame.from_dict([dictMol])
    print(molData)
    selectedMolData = molData.drop(columns=['mol','molIMG','molPixmap','mwStr','obStr','isStr','epStr','eleList','Catoms','Hatoms','Oatoms'])
    print(selectedMolData)
    selectedMolData.to_csv('./_core/altDB.csv', index=False)



#molecule = thermalDexMolecule(SMILES='CC1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]', name='TNT', Q_dsc=770.26, Qunits='J g⁻¹', onsetT=90.14, #initT=74.62, proj='PDF_Test')
#molecule.genCoreValues()
#molecule.genAdditionalValues()

#createDatabase(molecule)