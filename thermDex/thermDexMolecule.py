import sys
import os
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Mol, MolFromSmiles, MolFromSmarts, rdmolfiles
from io import BytesIO
import base64
from numpy import log10
from pubchempy import get_compounds
from dataclasses import dataclass, field, asdict
import pandas as pd
import re
from datetime import datetime
from PyQt5.QtGui import QPixmap
import random
import string

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
    hammerDrop: str = 'Not Performed'
    friction: str = 'Not Performed'

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
    molIMG: any = None
    mwStr: str = None
    obStr: str = None
    isStr: str = None
    epStr: str = None
    eleList: any = None
    Catoms: int = None
    Hatoms: int = None
    Oatoms: int = None
    yoshidaMethod: str = 'Pfizer'
    dataFolder: str = ''
    noDSCPeak: str = ''
    measuredTD24: bool = False
    empiricalTD24: float = None

    def genMol(self):
        RDmol = MolFromSmiles(self.SMILES)
        self.mol = RDmol
        return RDmol

    def molToQPixmap(self):
        # Generate a molecular drawing as a PNG image
        opts = Draw.MolDrawOptions()
        opts.bondLineWidth = 2.
        img = Draw.MolToImage(self.mol, size=(550, 275), options=opts)

        # Convert the image to a byte array
        byte_array = BytesIO()
        img.save(byte_array, format='PNG')

        # Convert the byte array to a QPixmap and display it
        pixmap = QPixmap()
        pixmap.loadFromData(byte_array.getvalue())
        self.molPixmap = pixmap
        return pixmap

    def molToIMG(self):
        # Generate a molecular drawing as a PNG image
        opts = Draw.MolDrawOptions()
        opts.bondLineWidth = 5.
        img = Draw.MolToImage(self.mol, size=(1200, 600), options=opts)

        # Convert the image to a byte array
        #byte_array = BytesIO()
        #img.save(byte_array, format='PNG')

        # Convert the byte array to a QPixmap and display it
        # pixmap = QPixmap()
        # pixmap.loadFromData(byte_array.getvalue())
        self.molIMG = img
        return img

    def molToBytes(self):
        image = self.molToIMG()
        byte_array = BytesIO()
        image.save(byte_array, format='PNG')
        return base64.b64encode(byte_array.getvalue()).decode('utf-8')

    def mwFromMol(self):
        cmpdMW = Descriptors.MolWt(self.mol)
        self.MW = cmpdMW
        self.mwStr = "{:.2f}".format(cmpdMW)

    def HEGFromMol(self, highEnergyGroups):
        fullMatch = 0
        with open(highEnergyGroups, "r") as HEGroups: 
            for line in HEGroups:
                HeSubstructure = MolFromSmiles(line)
                fullmatchList = Mol.GetSubstructMatches(self.mol, HeSubstructure)
                if len(fullmatchList) > 0:
                    print('High Energy Group Found: ' + line[:-1])
                    for match in range(0, len(fullmatchList)):
                        self.HEG_list += [line[:-1]]
                fullMatch += len(fullmatchList)

        nitroCount = self.HEG_list.count('C[N+](=O)[O-]')
        nitroCorrection = -2*nitroCount
        if nitroCount >= 1:
            for nitro in range(0, nitroCount):
                self.HEG_list.remove('CNO')
                self.HEG_list.remove('CN=O')

        alkylOximeCount = self.HEG_list.count('C=NOC')
        alkylOximeCorrection = -2*alkylOximeCount
        if alkylOximeCount >= 1:
            for alkylOxime in range(0, alkylOximeCount):
                self.HEG_list.remove('ON=C')
                self.HEG_list.remove('CON')

        self.HEG = fullMatch + nitroCorrection + alkylOximeCorrection

    def EFGFromMol(self, expEnergyGroups):
        expMatch= 0
        with open(expEnergyGroups, "r") as expGroups: 
            for line in expGroups:
                expSubstructure = MolFromSmiles(line)
                expmatchList = Mol.GetSubstructMatches(self.mol, expSubstructure)
                if len(expmatchList) > 0:
                    print('Explosive Group Found: ' + line[:-1])
                    for match in range(0, len(expmatchList)):
                        self.EFG_list += [line[:-1]]
                expMatch += len(expmatchList)

        nitroCount = self.EFG_list.count('C[N+](=O)[O-]')
        nitroCorrection = -2*nitroCount
        if nitroCount >= 1:
            for nitro in range(0, nitroCount):
                self.EFG_list.remove('CNO')
                self.EFG_list.remove('CN=O')

        alkylOximeCount = self.EFG_list.count('C=NOC')
        alkylOximeCorrection = -2*alkylOximeCount
        if alkylOximeCount >= 1:
            for alkylOxime in range(0, alkylOximeCount):
                self.EFG_list.remove('ON=C')
                self.EFG_list.remove('CON')

        self.EFG = expMatch + nitroCorrection + alkylOximeCorrection

        if self.EFG >= 1:
            self.oreo += 8
        elif self.EFG == 0:
            self.oreo += 1
        else:
            print("Error: O.R.E.O.S. EFG number gives unexpected value: " + EFG)

        print(f'| OREOS | EFG: {self.EFG} -> OREOS SubTotal: {self.oreo}')

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

    def makeFolderForMolData(self):
        if self.dataFolder == '' or self.dataFolder == None:
            print('\nMaking Folder...')
            now = datetime.now()
            neatNow = now.strftime("%d-%b-%Y_%H-%M-%S")
            uniqueID = ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))
            folderName = f'./_core/UserAddedData/{str(self.molForm)}_{str(neatNow)}_{uniqueID}'
            os.mkdir(folderName)
            self.dataFolder = folderName
            print(f'{os.path.abspath(folderName)}')
            return folderName
        else:
            print('Folder provided, using that.')
            pass

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
        print(f'| OREOS | RoS: {self.RoS_des} -> OREOS SubTotal: {self.oreo}')

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
        print(f'| OREOS | OB: {self.OB_des} -> OREOS SubTotal: {self.oreo}')

    def oreoOnsetTadjustment(self):
        if self.onsetT != 'nan' and self.onsetT != '' and self.onsetT != None and self.noDSCPeak == '':
            print('Onset Temperature given as ' + str(self.onsetT))
            self.onsetT = float(self.onsetT)
            if self.onsetT <= 125:
                self.oreo += 8
            elif 125 < self.onsetT <= 200:
                self.oreo += 4
            elif 200 < self.onsetT <= 300:
                self.oreo += 2
            elif self.onsetT > 300:
                self.oreo += 1
        elif self.noDSCPeak == 'True':
            self.onsetT = None
        else:
            self.onsetT = None
        print(f'| OREOS | Tonset: {self.onsetT} -> OREOS SubTotal: {self.oreo}')
             
    def oreoSafeScaleCal(self):
        self.oreoSmallScale_val = self.oreo + 1
        self.oreoTensScale_val = self.oreo + 2
        self.oreoHundsScale_val = self.oreo + 4
        self.oreoLargeScale_val = self.oreo + 8

        scaleList = [self.oreoSmallScale_val, self.oreoTensScale_val, self.oreoHundsScale_val, self.oreoLargeScale_val]
        hazardList = []

        if self.onsetT == None and self.noDSCPeak == '':
            hazardValuesRangeList = [15, 22]
        elif self.onsetT != None:
            hazardValuesRangeList = [18, 28]
        elif self.noDSCPeak == 'True':
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
        if self.initT != 'nan' and self.initT != '' and self.initT != None and self.noDSCPeak == '':
            print('Tinit given as ' + str(self.initT))
            self.initT = float(self.initT)
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
        if T_yoshida != 'nan' and T_yoshida != '' and T_yoshida != None and self.Q_dsc != 'nan' and self.Q_dsc != '' and self.Q_dsc != None and self.noDSCPeak == '':
                # Convert Qdsc to cal/g if needed
                print('Qdsc given as ' + str(self.Q_dsc) + ' ' + self.Qunits)
                self.Q_dsc = float(self.Q_dsc)
                if self.Qunits == 'J g<sup>-1</sup>' or self.Qunits == 'J g⁻¹':
                    calQ = self.Q_dsc/4.184
                elif self.Qunits == 'cal g<sup>-1</sup>' or self.Qunits == 'cal g⁻¹':
                    calQ = self.Q_dsc

                print(f'Q_DSC expressed as cal/g for the calculation (calQ var): {calQ}')
             
                # Impact Sensitvity (IS)
                impactSens = log10(calQ) - (ISConstant*(log10(abs(T_yoshida-25)))) - 0.98
                print('IS = ' + str(impactSens))
                print()
             
                # Explosive Propagation (EP)
                exProp = log10(calQ) - EPConstant*(log10(abs(T_yoshida-25))) - 1.67
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

    def genCoreValues(self, highEnergyGroups, expEnergyGroups):
        self.genMol()
        self.mwFromMol()
        self.HEGFromMol(highEnergyGroups)
        self.EFGFromMol(expEnergyGroups)
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
        self.molToQPixmap()
        self.genAdditionalValues()
