from PyQt5.QtGui import QPixmap
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Mol, MolFromSmiles, MolFromSmarts, rdmolfiles
from dataclasses import dataclass, field
from pubchempy import get_compounds
import cirpy
from io import BytesIO
from numpy import log10
import re

highEnergyGroups = ".\\_core\\HighEnergyGroups.csv"
expEnergyGroups = ".\\_core\\ExplosiveGroups.csv"

@dataclass
class thermalDexMolecule:
    # Class for keeping track of Molecules in ThermalDex

    # User Provided Values
    SMILES: str
    name: str = ''
    mp: float = None
    mpEnd: float = None
    Q_dsc: float = None
    Qunits: str = 'J g<sup>-1</sup>'
    onsetT: float = None
    initT: float = None
    proj: str = ''

    # Calculated Values
    MW: float = None
    molForm: str = None
    eleComp: str = None
    HEG: int = None
    HEG_list: list = field(default_factory=list)
    EFG: int = None
    EFG_list: list = field(default_factory=list)
    RoS_val: float = None
    RoS_des: str = None
    OB_val: float = None
    OB_des: str = None
    IS_val: float = None
    IS_des: str = None
    EP_val: float = None
    EP_des: str = None
    Td24: float = None
    oreo: int = 0
    oreoSmallScale_val: int = None
    oreoSmallScale_des: str = None
    oreoTensScale_val: int = None
    oreoTensScale_des: str = None
    oreoHundsScale_val: int = None
    oreoHundsScale_des: str = None
    oreoLargeScale_val: int = None
    oreoLargeScale_des: str = None
    IUPAC_name: str = ''
    CAS_no: str = ''
 
    # Internal Processing Values
    mol: any = None
    molPixmap: any = None
    mwStr: str = None
    obStr: str = None
    isStr: str = None
    epStr: str = None
    eleList: any = None
    Catoms: int = None
    Hatoms: int = None
    Oatoms: int = None
    yoshidaMethod: str = 'Pfizer'

    def genMol(self):
        RDmol = MolFromSmiles(self.SMILES)
        self.mol = RDmol
        return RDmol

    def molToQPixmap(self):
        # Generate a molecular drawing as a PNG image
        img = Draw.MolToImage(self.mol)

        # Convert the image to a byte array
        byte_array = BytesIO()
        img.save(byte_array, format='PNG')

        # Convert the byte array to a QPixmap and display it
        pixmap = QPixmap()
        pixmap.loadFromData(byte_array.getvalue())
        self.molPixmap = pixmap
        return pixmap

    def mwFromMol(self):
        cmpdMW = Descriptors.MolWt(self.mol)
        self.MW = cmpdMW
        self.mwStr = "{:.2f}".format(cmpdMW)

    def HEGFromMol(self):
        fullMatch = 0
        with open(highEnergyGroups, "r") as HEGroups: 
            for line in HEGroups:
                HeSubstructure = MolFromSmiles(line)
                fullmatchList = Mol.GetSubstructMatches(self.mol, HeSubstructure)
                if len(fullmatchList) > 0:
                    print('High Energy Group Found: ' + line[:-1])
                    self.HEG_list += [line[:-1]]
                fullMatch += len(fullmatchList)

        self.HEG = fullMatch

    def EFGFromMol(self):
        expMatch= 0
        with open(expEnergyGroups, "r") as expGroups: 
            for line in expGroups:
                expSubstructure = MolFromSmiles(line)
                expmatchList = Mol.GetSubstructMatches(self.mol, expSubstructure)
                if len(expmatchList) > 0:
                    print('Explosive Group Found: ' + line[:-1])
                    self.EFG_list += [line[:-1]]
                expMatch += len(expmatchList)

        self.EFG = expMatch
        if self.EFG >= 1:
            self.oreo += 8
        elif self.EFG == 0:
            self.oreo += 1
        else:
            print("Error: O.R.E.O.S. EFG number gives unexpected value: " + EFG)   

    def eleCompFromMol(self):
        self.molForm = rdMolDescriptors.CalcMolFormula(self.mol)
        eleComp = ""
        eleCompList = []
        eleList = []
        niceList = []
        atomPattern = r"Cl\d*|H\d*|O\d*|N\d*|Si\d*|S\d*|F\d*|Cs\d*|Br\d*|I\d*|B\d*|Al\d*|Na\d*|K\d*|Mg\d*|Zn\d*|Ti\d*|Pd\d*|C\d*"
        match = re.findall(atomPattern, self.molForm, re.I)
        if match:
            items = match
        for ele in match:
            ment = re.findall(r"([a-z]+)([0-9]+)?", ele, re.I)
            matchedEleComp = ment[0]
            eleList += [matchedEleComp]
            
        for compostion in eleList:
            eleCompList += compostion[::-1]
            niceList += [('').join(compostion[::-1])]

        self.eleComp = (', ').join(niceList)
        self.eleList = eleList
        return eleList

    def CHOFromEleComp(self):
        carbonAtoms = [ele[1] for ele in self.eleList if "C" in ele[0]]
        hydrogenAtoms = [ele[1] for ele in self.eleList if "H" in ele[0]]
        oxygenAtoms = [ele[1] for ele in self.eleList if "O" in ele[0]]
     
        # Get Number of Carbon Atoms
        if carbonAtoms == []:
            Catoms = 0
        elif carbonAtoms[0] == '':
            Catoms = 1
        else:
            Catoms = int(carbonAtoms[0])

        # Get Number of Hydrogen Atoms
        if hydrogenAtoms == []:
            Hatoms = 0
        elif hydrogenAtoms[0] == '':
            Hatoms = 1
        else:
            Hatoms = int(hydrogenAtoms[0])
        
        # Get Number of Oxygen Atoms
        if oxygenAtoms == []:
            Oatoms = 0
        elif oxygenAtoms[0] == '':
            Oatoms = 1
        else:
            Oatoms = int(oxygenAtoms[0])

        self.Catoms = Catoms
        self.Hatoms = Hatoms
        self.Oatoms = Oatoms
        return Catoms, Hatoms, Oatoms

    def RoSFromEleComp(self): 
        ruleSixCrit = 6*self.HEG - self.Catoms
        if ruleSixCrit > 0:
            ruleSix = " <b>(Explosive)</b>"
            self.oreo += 8
        else:
            ruleSix = " <b>(Not Explosive)</b>"
            self.oreo += 2 
        self.RoS_val = ruleSixCrit
        self.RoS_des = ruleSix

    def OBFromEleComp(self):
        oxygenBalance = (-1600*((2*self.Catoms)+(self.Hatoms/2)-self.Oatoms))/self.MW
        if oxygenBalance > 160:
            obRisk = "<b>(Low Risk)</b>"
            self.oreo += 2
        elif oxygenBalance > 80 and oxygenBalance <= 160:
            obRisk = "<b>(Medium Risk)</b>"
            self.oreo += 4
        elif oxygenBalance >= -120 and oxygenBalance <= 80:
            obRisk = "<b>(High Risk)</b>"
            self.oreo += 8
        elif oxygenBalance >= -240 and oxygenBalance < -120:
            obRisk = "<b>(Medium Risk)</b>"
            self.oreo += 4
        elif oxygenBalance < -240:
            obRisk = "<b>(Low Risk)</b>"
            self.oreo += 2
        self.OB_val = oxygenBalance
        self.OB_des = obRisk
        self.obStr = "{:.2f}".format(oxygenBalance)

    def oreoOnsetTadjustment(self):
        if self.onsetT != 'nan' and self.onsetT != '' and self.onsetT != None:
            print('Onset Temperature given as ' + str(self.onsetT))
            self.onsetT = float(self.onsetT)
            if self.onsetT <= 125:
                self.oreo += 8
            elif self.onsetT in range(126,200):
                self.oreo += 4
            elif self.onsetT in range(200,300):
                self.oreo += 2
            elif self.onsetT >= 300:
                self.oreo += 1
        else:
            self.onsetT = None
             
    def oreoSafeScaleCal(self):
        self.oreoSmallScale_val = self.oreo + 1
        self.oreoTensScale_val = self.oreo + 2
        self.oreoHundsScale_val = self.oreo + 4
        self.oreoLargeScale_val = self.oreo + 8

        scaleList = [self.oreoSmallScale_val, self.oreoTensScale_val, self.oreoHundsScale_val, self.oreoLargeScale_val]
        hazardList = []

        if self.onsetT == None:
            hazardValuesRangeList = [15, 22]
        elif self.onsetT != None:
            hazardValuesRangeList = [18, 28]

        for scale in scaleList:
            if scale < hazardValuesRangeList[0]:
                oreosHazard = "Low Hazard"
                hazardList.append(oreosHazard)
            elif scale in range(hazardValuesRangeList[0], hazardValuesRangeList[1]):  
                # Note: Remember the last integer in range isnt included thats why I can get away with a list of 2 items. 
                oreosHazard = "Medium Hazard"
                hazardList.append(oreosHazard)
            elif scale >= hazardValuesRangeList[1]:
                oreosHazard = "High Hazard"
                hazardList.append(oreosHazard)
            else:
                print('Error: ' + str(scale))
        print(hazardList)
        self.oreoSmallScale_des,self.oreoTensScale_des,self.oreoHundsScale_des,self.oreoLargeScale_des = hazardList


    def Td24FromThermalProps(self):
        # Estimation of Maximum Recommended Process Temperature To Avoid Hazardous Thermal Decomposition
        if self.initT != 'nan' and self.initT != '' and self.initT != None:
            print('Tinit given as ' + str(self.initT))
            d24Temp = 0.7*self.initT - 46
            print('T_D24 =' + str(d24Temp))

            self.Td24 = d24Temp

    def yoshidaFromThermalProps(self):
        # Determin if Original Yoshida calculation is to be used or if Pfizer modification is to be used.
        if self.yoshidaMethod == 'Yoshida':
            T_yoshida = self.onsetT
            ISConstant = 0.72
            EPConstant = 0.38
        elif self.yoshidaMethod == 'Pfizer':
            T_yoshida = self.initT
            ISConstant = 0.54
            EPConstant = 0.285

        # Perform the calculation
        if T_yoshida != 'nan' and T_yoshida != '' and T_yoshida != None and self.Q_dsc != 'nan' and self.Q_dsc != '' and self.Q_dsc != None:
                # Convert Qdsc to cal/g if needed
                print('Qdsc given as ' + str(self.Q_dsc) + ' ' + self.Qunits)
                self.Q_dsc = float(self.Q_dsc)
                if self.Qunits == 'J g<sup>-1</sup>':
                    calQ = self.Q_dsc/4.184
                elif self.Qunits == 'cal g<sup>-1</sup>':
                    calQ == self.Q_dsc
             
                # Impact Sensitvity (IS)
                impactSens = log10(calQ) - ISConstant*(log10(abs(self.onsetT-25))) - 0.98
                print('IS = ' + str(impactSens))
             
                # Explosive Propagation (EP)
                exProp = log10(calQ) - EPConstant*(log10(abs(self.onsetT-25))) - 1.67
                print('EP = ' + str(exProp))
          
                # Warning Text
                if impactSens >= 0:
                    impact = ' <b>(Impact Sensitive)</b>'
                elif impactSens < 0:
                    impact = ' <b>(Not Impact Sensitive)</b>'
                if exProp >= 0:
                    explos = ' <b>(Propagates)</b>'
                elif exProp < 0:
                    explos = ' <b>(Should Not Propagate)</b>'

                # Save to thermalDexMolecule
                self.IS_val = impactSens
                self.IS_des = impact
                self.EP_val = exProp
                self.EP_des = explos
                self.isStr = "{:.2f}".format(impactSens)
                self.epStr = "{:.2f}".format(exProp)


    def IUPACNameFromSMILES(self):
        foundNames = get_compounds(self.SMILES, namespace='smiles')
        IUPACName = foundNames[0].iupac_name
        print(IUPACName)
        self.IUPAC_name = IUPACName

    def genCoreValues(self):
        self.genMol()
        self.mwFromMol()
        self.HEGFromMol()
        self.EFGFromMol()
        self.eleCompFromMol()
        self.CHOFromEleComp()
        self.RoSFromEleComp()
        self.OBFromEleComp()
        self.oreoOnsetTadjustment()
        self.oreoSafeScaleCal()
        self.yoshidaFromThermalProps()
        self.Td24FromThermalProps()

    def genAdditionalValues(self):
        self.IUPACNameFromSMILES()

    def genAllValues(self):
        self.genCoreValues()
        self.genAdditionalValues()

    def genCoreValuesFromMol(self):
        self.mwFromMol()
        self.HEGFromMol()
        self.EFGFromMol()
        self.eleCompFromMol()
        self.CHOFromEleComp()
        self.RoSFromEleComp()
        self.OBFromEleComp()
        self.oreoOnsetTadjustment()
        self.oreoSafeScaleCal()
        self.yoshidaFromThermalProps()
        self.Td24FromThermalProps()


nicatinamide = thermalDexMolecule(SMILES="ClC(=O)C1=CC=NC=C1 |c:5,7,t:3|",onsetT=80.63,Q_dsc=670,initT=93.72,yoshidaMethod='Yoshida')
nicatinamide.genAllValues()
print('\n')
print(nicatinamide)
print('\n')