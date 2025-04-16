#import sys
#from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Mol, MolFromSmiles, MolFromSmarts, rdmolfiles
#from io import BytesIO
#from numpy import log10
#from pubchempy import get_compounds
#from dataclasses import dataclass, field, asdict
#import re
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene, QFrame, QTableWidget, QTableWidgetItem, QTabWidget, QGraphicsPixmapItem,  QMessageBox, QComboBox #, QTableView, QToolTip
from PyQt5.QtGui import QPixmap, QColor, QIcon, QRegExpValidator #QValidator #, QCursor
from PyQt5.QtCore import Qt, QRegExp, pyqtSignal #, QCoreApplication
#from thermDex.thermDexReport import *
from thermDex.thermDexMolecule import * #thermalDexMolecule
from thermDex.thermDexHTMLRep import *
from thermDex.attachedFileManager import *
from thermDex.Section import Section
from thermDex.thermDexPlots import *
from thermDex.thermDexPandasTools import *
import pyperclip
import configparser
from numpy import log10
from contextlib import redirect_stdout
from os import path
from shutil import copy2

versionNumber = "1.0.3b"

try:
    import pyi_splash
    pyi_splash.close()
    window.raise_()
except:
    pass


def importConfig():
    conf = open('.\\_core\\ThermalDex.config', 'r')
    confCounter = 0
    for line in conf:
        #print(confCounter)
        if confCounter == 4:
           defaultDB = line.strip("\n")
           confCounter += 1
        elif confCounter == 8:
           highEnergyGroups = line.strip("\n")
           confCounter += 1
        elif confCounter == 12:
           expEnergyGroups = line.strip("\n")
           confCounter += 1
        else:
           confCounter += 1

    #print(defaultDB)
    return defaultDB, highEnergyGroups, expEnergyGroups

def altImportConfig():
    config = configparser.ConfigParser()
    config.read('./_core/ThermalDex.ini')
    defaultDB = config.get('Database', 'defaultDB')
    highEnergyGroups = config.get('Lists of Groups', 'highenergygroups')
    expEnergyGroups = config.get('Lists of Groups', 'expenergygroups')
    yoshidaMethod = config.get('Calculations', 'yoshidacalcs')
    qdscUnits = config.get('Default Values', 'qdscunits')
    ambertd24limit = config.get('Warnings', 'ambertd24limit')
    redtd24limit = config.get('Warnings', 'redtd24limit')
    oreohazardlimit = config.get('Warnings', 'oreohazardlimit')
    oreohazardwarningscale = config.get('Warnings', 'oreohazardwarningscale')

    return defaultDB, highEnergyGroups, expEnergyGroups, yoshidaMethod, qdscUnits, ambertd24limit, redtd24limit, oreohazardlimit, oreohazardwarningscale

class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)

class QVLine(QFrame):
    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QFrame.VLine)
        self.setFrameShadow(QFrame.Sunken)

class ClickableLineEdit(QLineEdit):
    clicked = pyqtSignal() # signal when the text entry is left clicked

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton: self.clicked.emit()
        else: super().mousePressEvent(event)

class MolDrawer(QWidget):
    def __init__(self):
        super().__init__()

        self.init_ui()

    def init_ui(self):
        self.df = None
        #current_index = 0
        self.current_index = 0
        self.result_smiles = None
        self.selectedDatabase = None
        self.error_flag = None
        # Set up the main layout
        layout = QVBoxLayout()

        # Create a tab widget
        self.tab_widget = QTabWidget()

        # Tab for molecule rendering
        self.molecule_tab = QWidget()
        molecule_layout = QVBoxLayout()

        # Display area for the molecular drawing
        self.mol_display = QGraphicsView(self)
        top_info_sublayout = QHBoxLayout()
        top_info_sublayout.addWidget(QLabel('Molecule:'))
        top_info_sublayout.addStretch()
        self.approval_needed = QLabel("<h3 style='color: red;'>Seek Approval Before Use</h3>")
        top_info_sublayout.addWidget(self.approval_needed)
        self.approval_needed.hide()
        molecule_layout.addLayout(top_info_sublayout)
        molecule_layout.addWidget(self.mol_display)

        ResultsContainLayout = QHBoxLayout()
        ResultsLeftLayout = QVBoxLayout()
        ResultsRightLayout = QVBoxLayout()

	#Add labels for calculated values
        self.mwText = 'MW: '
        self.HEGText = 'Number of High Energy Groups:'
        self.EFGText = 'Number of Explosive Functional Groups:'
        self.eleText = 'Elemental Composition: '
        self.RoSText = 'Rule of Six: '
        self.obText = 'Oxygen Balance: '
        self.ISText = 'Yoshida Impact Sensitivity: '
        self.EPText = 'Yoshida Explosive Propagation: '
        self.Td24Text = 'T<sub>D24</sub>: '

        self.mwLabel = QLabel(self.mwText)
        self.mwLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsLeftLayout.addWidget(self.mwLabel)
        self.HEGlabel = QLabel(self.HEGText)
        self.HEGlabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsLeftLayout.addWidget(self.HEGlabel)
        self.EFGlabel = QLabel(self.EFGText)
        self.EFGlabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsLeftLayout.addWidget(self.EFGlabel)
        self.eleLabel = QLabel(self.eleText)
        self.eleLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsLeftLayout.addWidget(self.eleLabel)
        self.RoSLabel = QLabel(self.RoSText)
        self.RoSLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsRightLayout.addWidget(self.RoSLabel)
        self.obLabel = QLabel(self.obText)
        self.obLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsRightLayout.addWidget(self.obLabel)
        self.ISLabel = QLabel(self.ISText)
        self.ISLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsRightLayout.addWidget(self.ISLabel)
        self.EPLabel = QLabel(self.EPText)
        self.EPLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsRightLayout.addWidget(self.EPLabel)
        self.Td24Label = QLabel(self.Td24Text)
        self.Td24Label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        ResultsRightLayout.addWidget(self.Td24Label)

        ResultsContainLayout.addWidget(QVLine())
        ResultsContainLayout.addLayout(ResultsLeftLayout)
        ResultsContainLayout.addWidget(QVLine())
        ResultsContainLayout.addLayout(ResultsRightLayout)
        ResultsContainLayout.addWidget(QVLine())
        molecule_layout.addLayout(ResultsContainLayout)
        molecule_layout.addWidget(QHLine())

        self.tableLabel = QLabel('O.R.E.O.S. Assessment of Hazard by Scale:')
        molecule_layout.addWidget(self.tableLabel)
        tableLayout = QHBoxLayout()
        self.table = QTableWidget(1, 4)
        self.table.setHorizontalHeaderLabels(['<5 g', '5 to <100 g', '100 to 500 g', '>500 g'])
        self.table.verticalHeader().setVisible(False)
        #self.table.setStyleSheet("QTableWidget::item { border-bottom: 2px solid black; }")
        self.table.setMaximumHeight(53)
        self.table.setMaximumWidth(402)
        self.table.setMinimumHeight(53)
        self.table.setMinimumWidth(402)
        #self.table.setAlignment(Qt.AlignVCenter)
        tableLayout.addWidget(self.table)

        molecule_layout.addLayout(tableLayout)
        molecule_layout.addWidget(QHLine())

        # Input field for SMILES string
        #self.smiles_input = QLineEdit(self)
        self.smiles_input = ClickableLineEdit(self)
        self.smiles_input.clicked.connect(self.mrvToSMILES)
        molecule_layout.addWidget(QLabel('Enter SMILES String:'))
        molecule_layout.addWidget(self.smiles_input)

        InputContainLayout = QHBoxLayout()
        InputLeftLayout = QVBoxLayout()
        InputRightLayout = QVBoxLayout()
        numValidator = QRegExpValidator(QRegExp(r'[-]?\d+[.]?\d*'))
        posNumValidator = QRegExpValidator(QRegExp(r'\d+[.]?\d*'))

        # Input field for Name string
        self.name_input = QLineEdit(self)
        InputLeftLayout.addWidget(QLabel('Name:'))
        nameUnitsSubLayout = QHBoxLayout()
        nameUnitsSubLayout.addWidget(self.name_input)
        nameUnitLabel = QLabel('°C    ')
        nameUnitLabel.setStyleSheet('color: white')
        nameUnitsSubLayout.addWidget(nameUnitLabel)
        InputLeftLayout.addLayout(nameUnitsSubLayout)

        # Input field for mp string
        self.mp_input = QLineEdit(self)
        self.mp_input.setMaximumWidth(110)
        self.mp_input.setValidator(numValidator)
        self.mpEnd_input = QLineEdit(self)
        self.mpEnd_input.setMaximumWidth(110)
        self.mpEnd_input.setValidator(numValidator)
        InputLeftLayout.addWidget(QLabel('m.p.:'))
        mpUnitsSubLayout = QHBoxLayout()
        mpUnitsSubLayout.addWidget(self.mp_input)
        mpUnitsSubLayout.addWidget(QLabel('  to  '))
        mpUnitsSubLayout.addWidget(self.mpEnd_input)
        mpUnitsSubLayout.addWidget(QLabel('°C    '))
        InputLeftLayout.addLayout(mpUnitsSubLayout)

        # Input field for Q_DSC string
        self.Qdsc_input = QLineEdit(self)
        self.Qdsc_input.setValidator(posNumValidator)
        InputRightLayout.addWidget(QLabel('Q<sub>DSC</sub>:'))
        QUnitsSubLayout = QHBoxLayout()
        QUnitsSubLayout.addWidget(self.Qdsc_input)
        #QUnitsSubLayout.addWidget(QLabel('cal g<sup>-1</sup> '))
        self.QUnitsSelection = QComboBox(self)
        self.QUnitsSelection.addItems(['J g⁻¹', 'cal g⁻¹'])
        self.QUnitsSelection.setCurrentIndex(int(qdscUnits))
        QUnitsSubLayout.addWidget(self.QUnitsSelection)
        InputRightLayout.addLayout(QUnitsSubLayout)

        # Input field for Onset string
        self.TE_input = QLineEdit(self)
        self.TE_input.setValidator(numValidator)
        InputRightLayout.addWidget(QLabel('Onset Temperature:'))
        TEUnitsSubLayout = QHBoxLayout()
        TEUnitsSubLayout.addWidget(self.TE_input)
        TEUnitsSubLayout.addWidget(QLabel('°C       '))
        InputRightLayout.addLayout(TEUnitsSubLayout)

        # Input field for Init string
        self.Tinit_input = QLineEdit(self)
        self.Tinit_input.setValidator(numValidator)
        InputRightLayout.addWidget(QLabel('Initiation Temperature:'))
        TinitUnitsSubLayout = QHBoxLayout()
        TinitUnitsSubLayout.addWidget(self.Tinit_input)
        TinitUnitsSubLayout.addWidget(QLabel('°C       '))
        InputRightLayout.addLayout(TinitUnitsSubLayout)

        # Input field for proj string
        self.proj_input = QLineEdit(self)
        InputLeftLayout.addWidget(QLabel('Project:'))
        projUnitsSubLayout = QHBoxLayout()
        projUnitsSubLayout.addWidget(self.proj_input)
        projUnitLabel = QLabel('°C    ')
        projUnitLabel.setStyleSheet('color: white')
        projUnitsSubLayout.addWidget(projUnitLabel)
        InputLeftLayout.addLayout(projUnitsSubLayout)

        # Input field for Hammer Drop Test
        hamSubLayout = QHBoxLayout()
        InputLeftLayout.addWidget(QLabel('Hammer Drop Test:'))
        self.hamSelection = QComboBox(self)
        self.hamSelection.addItems(['Not Performed', 'Positive', 'Negative'])
        hamSubLayout.addWidget(self.hamSelection)
        hamUnitLabel = QLabel('°C    ')
        hamUnitLabel.setStyleSheet('color: white')
        hamSubLayout.addWidget(hamUnitLabel)
        InputLeftLayout.addLayout(hamSubLayout)

        # Input field for Friction Test
        fricSubLayout = QHBoxLayout()
        InputRightLayout.addWidget(QLabel('Friction Test:'))
        self.fricSelection = QComboBox(self)
        self.fricSelection.addItems(['Not Performed', 'Positive', 'Negative'])
        fricSubLayout.addWidget(self.fricSelection)
        fricUnitLabel = QLabel('°C    ')
        fricUnitLabel.setStyleSheet('color: white')
        fricSubLayout.addWidget(fricUnitLabel)
        InputRightLayout.addLayout(fricSubLayout)

        InputContainLayout.addLayout(InputLeftLayout)
        #ResultsContainLayout.addWidget(QVLine())
        InputContainLayout.addLayout(InputRightLayout)
        #ResultsContainLayout.addWidget(QVLine())
        molecule_layout.addLayout(InputContainLayout)

        # Attach Data Files
        filesSubLayout = QHBoxLayout()
        self.attachedFilesLabel = QLabel('Attached Files:')
        self.attachedFilesLabel.resize(120, 120)
        filesSubLayout.addWidget(self.attachedFilesLabel)
        self.filesCount = QLabel('0 Attached Files')
        self.filesCount.resize(90, 120)
        filesSubLayout.addWidget(self.filesCount)
        self.attach_button = QPushButton('Add/View Files', self)
        self.attach_button.clicked.connect(self.openFileManager) 
        self.attach_button.setMaximumWidth(140)
        filesSubLayout.addWidget(self.attach_button)
        attachSpacerLabel = QLabel('hidden spacer')
        attachSpacerLabel.setStyleSheet('color: white')
        filesSubLayout.addWidget(attachSpacerLabel)
        attachSpacerLabelTwo = QLabel('hidden spacer')
        attachSpacerLabelTwo.setStyleSheet('color: white')
        filesSubLayout.addWidget(attachSpacerLabelTwo)
        molecule_layout.addSpacing(15)
        molecule_layout.addLayout(filesSubLayout)

        # Buttons
        buttonSpacerLabel = QLabel('hidden spacer')
        buttonSpacerLabel.setStyleSheet('color: white')
        molecule_layout.addWidget(buttonSpacerLabel)
        buttonContainLayout = QHBoxLayout()
        render_button = QPushButton('Evaluate Molecule', self)
        render_button.clicked.connect(self.render_molecule)
        render_button.setMaximumWidth(180)
        buttonContainLayout.addWidget(render_button)
        clear_button = QPushButton('Clear Sheet', self)
        clear_button.clicked.connect(self.resetToDefaultState)
        clear_button.setMaximumWidth(180)
        buttonContainLayout.addWidget(clear_button)
        #molecule_layout.addLayout(buttonContainLayout)
        msg_button = QPushButton('Save to PDF', self)
        msg_button.clicked.connect(self.createReport)
        msg_button.setMaximumWidth(180)
        buttonContainLayout.addWidget(msg_button)
        ''' Button for error testing/Debugging.
        error_testing_button = QPushButton('SecretErrorButton')
        error_testing_button.clicked.connect(self.errorTestingTool)
        buttonContainLayout.addWidget(error_testing_button)       
        '''
        molecule_layout.addLayout(buttonContainLayout)

        self.molecule_tab.setLayout(molecule_layout)
        self.tab_widget.addTab(self.molecule_tab, "Add")

        # Hide File Handling Widgets 
        self.attachedFilesLabel.hide()
        self.filesCount.hide()
        self.attach_button.hide()

        # Tab for Search
        search_tab = QWidget()
        search_layout = QVBoxLayout()

        # Entry widget for searching
        lbl_search = QLabel('Search:')
        entry_search = ClickableLineEdit(self) #QLineEdit()
        entry_search.clicked.connect(self.mrvToSMILES)
        self.searchTypeSelection = QComboBox(self)
        self.searchTypeSelection.addItems(['SMILES', 'Name', 'Project', 'MW', 'Qdsc', 'Tinit', 'Tonset', 'Td24', 'O.R.E.O.S. at >500 g'])
        smilesSearchList = ['Substructure', 'Exact']
        nameProjSearchList = ['is', 'contains']
        valuesSearchList = ['=', '<', '>', '</=', '>/=']
        self.listOfSearchTypes = [smilesSearchList, nameProjSearchList, nameProjSearchList, valuesSearchList, valuesSearchList, valuesSearchList, valuesSearchList, valuesSearchList, valuesSearchList]
        self.searchSubType = QComboBox(self)
        self.searchSubType.addItems(smilesSearchList)
        self.searchTypeSelection.currentIndexChanged.connect(self.set_sub_search)

        result_label = QLabel('click search')
        counter_label = QLabel('none')

        def search_database():
             self.readDatabase = pd.read_csv(defaultDB) #, index_col=0) #, encoding='mbcs')

             if entry_search.text() == '' or entry_search.text() == None:    
                 self.selectedDatabase = self.readDatabase
                 show_result(self, self.readDatabase, True)

             elif self.searchTypeSelection.currentText() == 'SMILES':
                 if self.searchSubType.currentText() == 'Exact':
                     searchSMILES = entry_search.text()
                     row_index = self.readDatabase.index[self.readDatabase['SMILES'] == searchSMILES].tolist()
                     foundDataFrame = self.readDatabase.iloc[row_index]
                     foundDataFrame.reset_index(drop=True)
                     self.selectedDatabase = foundDataFrame
                     show_result(self, foundDataFrame, True)

                 elif self.searchSubType.currentText() == 'Substructure':
                    try:
                        searchTest = MolFromSmiles(entry_search.text())
                        #checkifrealbyMW = Descriptors.MolWt(searchTest)
                        if entry_search.text() != '' and entry_search.text() != None and searchTest is not None:
                            smilesList = self.readDatabase['SMILES'].tolist()  #.index.values
                            #print(smilesList)
                            foundList = []
                            for smile in smilesList:
                                searchStructure = MolFromSmiles(smile)
                                fullmatchList = Mol.GetSubstructMatches(searchStructure, searchTest)
                                if len(fullmatchList) > 0:
                                    print('Substructure Match Found: ' + smile)
                                    foundList += [smile]
                            #print('\n\n')
                            print(foundList)
                            indexList = []
                            for foundMatch in foundList:
                                row_index = self.readDatabase.index[self.readDatabase['SMILES'] == foundMatch].tolist()
                                indexList += row_index
                            #print(indexList)
                            #print('\n\n')
                            foundDataFrame = self.readDatabase.iloc[indexList]
                            print(foundDataFrame)
                            print('\n\n')
                            foundDataFrame.reset_index(drop=True)
                            self.selectedDatabase = foundDataFrame
                            show_result(self, foundDataFrame, True)
                    except:
                        errorInfo = "Enter Valid SMILES"
                        self.interactiveErrorMessage(errorInfo)

             elif self.searchTypeSelection.currentText() == 'Name':
                 if self.searchSubType.currentText() == 'is':
                     searchName = entry_search.text()
                     row_index = self.readDatabase.index[self.readDatabase['name'] == searchName].tolist()
                     foundDataFrame = self.readDatabase.iloc[row_index]
                     foundDataFrame.reset_index(drop=True)
                     self.selectedDatabase = foundDataFrame
                     show_result(self, foundDataFrame, True)

                 elif self.searchSubType.currentText() == 'contains':
                     searchName = entry_search.text()
                     nameList = self.readDatabase['name'].tolist()  #.index.values
                     cleanNameList = []
                     foundList = []
                     for name in nameList:
                         if name == name: # a check to remove nan values
                             cleanNameList += [name]
                     for name in cleanNameList:
                         print(name)
                         if searchName.lower() in name.lower():
                             foundList += [name]
                     #print('\n\n')
                     #print(foundList)
                     foundList = list(dict.fromkeys(foundList))
                     indexList = []
                     for foundMatch in foundList:
                         row_index = self.readDatabase.index[self.readDatabase['name'] == foundMatch].tolist()
                         indexList += row_index
                     foundDataFrame = self.readDatabase.iloc[indexList]
                     print(foundDataFrame)
                     print('\n\n')
                     foundDataFrame.reset_index(drop=True)
                     self.selectedDatabase = foundDataFrame
                     show_result(self, foundDataFrame, True)

             elif self.searchTypeSelection.currentText() == 'Project':
                 if self.searchSubType.currentText() == 'is':
                     searchName = entry_search.text()
                     row_index = self.readDatabase.index[self.readDatabase['proj'] == searchName].tolist()
                     foundDataFrame = self.readDatabase.iloc[row_index]
                     foundDataFrame.reset_index(drop=True)
                     self.selectedDatabase = foundDataFrame
                     show_result(self, foundDataFrame, True)

                 elif self.searchSubType.currentText() == 'contains':
                     searchName = entry_search.text()
                     nameList = self.readDatabase['proj'].tolist()  #.index.values
                     cleanNameList = []
                     foundList = []
                     for name in nameList:
                         if name == name: # a check to remove nan values
                             cleanNameList += [name]
                     for name in cleanNameList:
                         print(name)
                         if searchName.lower() in name.lower():
                             foundList += [name]
                     #print('\n\n')
                     #print(foundList)
                     foundList = list(dict.fromkeys(foundList))
                     indexList = []
                     for foundMatch in foundList:
                         row_index = self.readDatabase.index[self.readDatabase['proj'] == foundMatch].tolist()
                         indexList += row_index
                     foundDataFrame = self.readDatabase.iloc[indexList]
                     print(foundDataFrame)
                     print('\n\n')
                     foundDataFrame.reset_index(drop=True)
                     self.selectedDatabase = foundDataFrame
                     show_result(self, foundDataFrame, True)

             elif self.searchTypeSelection.currentText() == 'MW':
                 if self.searchSubType.currentText() == '=':
                     try:
                        searchName = float(entry_search.text())
                        row_index = self.readDatabase.index[self.readDatabase['MW'] == searchName].tolist()
                        foundDataFrame = self.readDatabase.iloc[row_index]
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo) 

                 elif self.searchSubType.currentText() == '<':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['MW'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value < float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['MW'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '>':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['MW'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value > float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['MW'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '</=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['MW'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value <= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['MW'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)           

                 elif self.searchSubType.currentText() == '>/=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['MW'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value >= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['MW'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)  

             elif self.searchTypeSelection.currentText() == 'Qdsc':
                 if self.searchSubType.currentText() == '=':
                     try:
                        searchName = float(entry_search.text())
                        row_index = self.readDatabase.index[self.readDatabase['Q_dsc'] == searchName].tolist()
                        foundDataFrame = self.readDatabase.iloc[row_index]
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo) 

                 elif self.searchSubType.currentText() == '<':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['Q_dsc'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value < float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['Q_dsc'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '>':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['Q_dsc'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value > float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['Q_dsc'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '</=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['Q_dsc'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value <= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['Q_dsc'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)           

                 elif self.searchSubType.currentText() == '>/=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['Q_dsc'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value >= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['Q_dsc'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)  

             elif self.searchTypeSelection.currentText() == 'Tinit':
                 if self.searchSubType.currentText() == '=':
                     try:
                        searchName = float(entry_search.text())
                        row_index = self.readDatabase.index[self.readDatabase['initT'] == searchName].tolist()
                        foundDataFrame = self.readDatabase.iloc[row_index]
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo) 

                 elif self.searchSubType.currentText() == '<':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['initT'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value < float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['initT'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '>':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['initT'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value > float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['initT'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        foundList = list(dict.fromkeys(foundList))
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '</=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['initT'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value <= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['initT'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        foundList = list(dict.fromkeys(foundList))
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)           

                 elif self.searchSubType.currentText() == '>/=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['initT'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value >= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['initT'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)               

             elif self.searchTypeSelection.currentText() == 'Tonset':
                 if self.searchSubType.currentText() == '=':
                     try:
                        searchName = float(entry_search.text())
                        row_index = self.readDatabase.index[self.readDatabase['onsetT'] == searchName].tolist()
                        foundDataFrame = self.readDatabase.iloc[row_index]
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo) 

                 elif self.searchSubType.currentText() == '<':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['onsetT'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value < float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['onsetT'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '>':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['onsetT'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value > float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['onsetT'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '</=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['onsetT'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value <= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['onsetT'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)           

                 elif self.searchSubType.currentText() == '>/=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['onsetT'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value >= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['onsetT'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)     

             elif self.searchTypeSelection.currentText() == 'Td24':
                 if self.searchSubType.currentText() == '=':
                     try:
                        searchName = float(entry_search.text())
                        row_index = self.readDatabase.index[self.readDatabase['Td24'] == searchName].tolist()
                        foundDataFrame = self.readDatabase.iloc[row_index]
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo) 

                 elif self.searchSubType.currentText() == '<':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['Td24'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value < float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['Td24'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '>':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['Td24'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value > float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['Td24'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '</=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['Td24'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value <= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['Td24'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)           

                 elif self.searchSubType.currentText() == '>/=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['Td24'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value >= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['Td24'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)     

             elif self.searchTypeSelection.currentText() == 'O.R.E.O.S. at >500 g':
                 if self.searchSubType.currentText() == '=':
                     try:
                        searchName = float(entry_search.text())
                        row_index = self.readDatabase.index[self.readDatabase['oreoLargeScale_val'] == searchName].tolist()
                        foundDataFrame = self.readDatabase.iloc[row_index]
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo) 

                 elif self.searchSubType.currentText() == '<':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['oreoLargeScale_val'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value < float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['oreoLargeScale_val'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '>':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['oreoLargeScale_val'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value > float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        foundList = list(dict.fromkeys(foundList))
                        indexList = []
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['oreoLargeScale_val'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

                 elif self.searchSubType.currentText() == '</=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['oreoLargeScale_val'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value <= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        indexList = []
                        foundList = list(dict.fromkeys(foundList))
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['oreoLargeScale_val'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)
                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)           

                 elif self.searchSubType.currentText() == '>/=':
                     try:
                        searchName = entry_search.text()
                        valueList = self.readDatabase['oreoLargeScale_val'].tolist()  #.index.values
                        foundList = []
                        for value in valueList:
                            print(value)
                            if value >= float(searchName):
                                foundList += [value]
                        print('\n\n')
                        print(foundList)
                        indexList = []
                        foundList = list(dict.fromkeys(foundList))
                        for foundMatch in foundList:
                            row_index = self.readDatabase.index[self.readDatabase['oreoLargeScale_val'] == foundMatch].tolist()
                            indexList += row_index
                        foundDataFrame = self.readDatabase.iloc[indexList]
                        print(foundDataFrame)
                        print('\n\n')
                        foundDataFrame.reset_index(drop=True)
                        self.selectedDatabase = foundDataFrame
                        show_result(self, foundDataFrame, True)

                     except:
                        errorInfo = "Enter a Valid number"
                        self.interactiveErrorMessage(errorInfo)   

             else:
                 errorInfo = "How!? How have you gotten to this error message... I... Ugh... Report this to the developer okay?"
                 self.interactiveErrorMessage(errorInfo)   

        def show_result(self, Database, resetIndex):
             #print(self)
             layout = self.layout()
             if Database.empty:
                  errorInfo = "No matches found. Try a diffrent search?"
                  self.interactiveErrorMessage(errorInfo)
                  counter_label.setText("")
                  result_label.setText("")
                  self.mol_result_display.hide()
                  self.molLabel.hide()
                  make_plot_label.hide()
                  make_plot_button.hide()
                  self.select_x_values.hide()
                  vs_label.hide()
                  self.select_y_values.hide()
                  try:
                      export_button.hide()
                      edit_button.hide()
                      del_button.hide()
                      prev_button.hide()
                      next_button.hide()
                      results_table_Label.hide()
                      results_table.hide()
                      results_approval_warning.hide()
                  except:
                      print('Not shown')
             if self.error_flag is not None:
                  self.error_message.setText('')
                  layout.removeWidget(self.error_message)
                  self.error_flag = None
             if resetIndex == True:
                  self.current_index = 0             
             if Database is not None and not Database.empty:
                  self.mol_result_display.show()
                  self.molLabel.show()
                  current_row = Database.iloc[self.current_index]
                  print(current_row)
                  print('\n\n')
                  dictRow = current_row.to_dict()
                  print(dictRow)
                  print('\n\n')
                  readMolecule = thermalDexMolecule(**dictRow)
                  print(readMolecule.MW)
                  self.result_smiles = current_row['SMILES']

                  if current_row['Td24'] != None and current_row['Td24'] != '':
                    d24Str = "{:.1f}".format(current_row['Td24'])
                    if int(ambertd24limit) >= current_row['Td24'] > int(redtd24limit):
                        formTd24 = f"T<sub>D24</sub>: <b style='color: orange;'> {d24Str} °C</b>"
                    elif current_row['Td24'] <= int(redtd24limit):
                        formTd24 = f"<h3 style='color: red;'>T<sub>D24</sub>: {d24Str} °C<br>Approval Needed Before Use</h3>"

                    else:
                        formTd24 = 'T<sub>D24</sub>: <b>' + d24Str + ' °C</b>'
                  else:
                    formTd24 = f"T<sub>D24</sub>: current_row['Td24']"

                  if self.check_if_oreos_need_approval(readMolecule) == 'Show Approval Message':
                      oreoApprovalWarning = "<h3 style='color: red;'>Approval Needed Before Use</h3>"
                  else:
                      oreoApprovalWarning = ''
                  results_approval_warning.setText(oreoApprovalWarning)

                  result_text = f"SMILES: {current_row['SMILES']}<br>Name: {current_row['name']}<br>High Energy Groups: {current_row['HEG']}<br>Explosive Functional Groups: {current_row['EFG']}<br>mp: {current_row['mp']} to {current_row['mpEnd']} °C<br>MW: {'{:.2f}'.format(current_row['MW'])} g mol<sup>-1</sup><br>Q<sub>DSC</sub>: {current_row['Q_dsc']} {current_row['Qunits']}<br>T<sub>onset</sub>: {current_row['onsetT']} °C<br>T<sub>init</sub>: {current_row['initT']} °C<br>Oxygen Balance: {'{:.2f}'.format(current_row['OB_val'])} {current_row['OB_des']}<br>Rule of Six: {current_row['RoS_val']} {current_row['RoS_des']}<br>Impact Sensitivity: {'{:.2f}'.format(current_row['IS_val'])} {current_row['IS_des']}<br>Explosive Propagation: {'{:.2f}'.format(current_row['EP_val'])} {current_row['EP_des']}<br>Yosida Calculation Method Used: {current_row['yoshidaMethod']}<br>{formTd24}<br>Hammer Drop Test: {current_row['hammerDrop']}<br>Friction Test: {current_row['friction']}<br>Project: {current_row['proj']}"
                  result_label.setText(result_text)
                  result_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
                  result_label.setWordWrap(True)

                  results_table.clearContents()

                  oreoSmall = current_row['oreoSmallScale_des']
                  oreoTens = current_row['oreoTensScale_des']
                  oreoHunds = current_row['oreoHundsScale_des']
                  oreoLarge = current_row['oreoLargeScale_des']

                  smallEntry = QTableWidgetItem(oreoSmall)
                  results_table.setItem(0, 0, smallEntry)
                  classColor = self.getColorForValue(oreoSmall)
                  smallEntry.setBackground(classColor)

                  tensEntry = QTableWidgetItem(oreoTens)
                  results_table.setItem(0, 1, tensEntry)
                  classColor = self.getColorForValue(oreoTens)
                  tensEntry.setBackground(classColor)

                  hundsEntry = QTableWidgetItem(oreoHunds)
                  results_table.setItem(0, 2, hundsEntry)
                  classColor = self.getColorForValue(oreoHunds)
                  hundsEntry.setBackground(classColor)

                  largeEntry = QTableWidgetItem(oreoLarge)
                  results_table.setItem(0, 3, largeEntry)
                  classColor = self.getColorForValue(oreoLarge)
                  largeEntry.setBackground(classColor)

                  counter_label.setText(f"Result {self.current_index + 1} of {len(Database)}")
                  #search_layout.insertWidget(4, edit_button) #  .addWidget(edit_button)
                  #search_layout.insertWidget(5, del_button)
                  edit_button.show()
                  del_button.show()
                  search_layout.insertWidget(5, result_label)
                  search_layout.insertWidget(6, results_table_Label)
                  search_layout.insertWidget(7, results_table)
                  search_layout.insertWidget(8, results_approval_warning)
                  search_layout.addWidget(counter_label)
                  prev_next_layout = QHBoxLayout()
                  prev_next_layout.addWidget(prev_button)
                  prev_next_layout.addWidget(next_button)
                  search_layout.addLayout(prev_next_layout)
                  export_button.show()
                  edit_button.show()
                  del_button.show()
                  make_plot_label.show()
                  make_plot_button.show()
                  self.plot_database = Database
                  self.plot_current_value = dictRow
                  self.select_x_values.show()
                  vs_label.show()
                  self.select_y_values.show()
                  results_table_Label.show()
                  results_table.show()
                  prev_button.show()
                  next_button.show()

                  if self.result_smiles is not None:
                      mol = MolFromSmiles(self.result_smiles)

                  else: 
                      mol = None

                  if mol is not None:

                      # Generate a molecular drawing as a PNG image
                      opts = Draw.MolDrawOptions()
                      opts.bondLineWidth = 2.
                      img = Draw.MolToImage(mol, size=(400, 150), options=opts)

                      # Convert the image to a byte array
                      byte_array = BytesIO()
                      img.save(byte_array, format='PNG')

                      # Convert the byte array to a QPixmap and display it
                      pixmap = QPixmap()
                      pixmap.loadFromData(byte_array.getvalue())
                      scene = QGraphicsScene()
                      scene.addPixmap(pixmap)
                      #scene.setSceneRect(0,50,500,300)
                      self.mol_result_display.setScene(scene)

        def prev_result(self):
             if self.current_index > 0:
                  #search_layout.removeWidget(self.mol_result_display)
                  #search_layout.removeWidget(self.molLabel)
                  self.current_index -= 1
                  show_result(self, self.selectedDatabase, False)

        def next_result(self, Database):
             if Database is not None:
                  if self.current_index < len(Database) - 1:
                       #search_layout.removeWidget(self.mol_result_display)
                       #search_layout.removeWidget(self.molLabel)
                       self.current_index += 1
                       show_result(self, self.selectedDatabase, False)

        # Search Buttons & display area for the molecular drawing
        self.mol_result_display = QGraphicsView(self)
        #self.mol_result_display.scale(0.7,0.7)
        self.molLabel = QLabel('Molecule:')
        search_layout.addWidget(self.molLabel)
        search_layout.addWidget(self.mol_result_display)
        self.mol_result_display.hide()
        self.molLabel.hide()

        btn_search = QPushButton('Search', clicked=search_database)
        edit_button = QPushButton('Edit')
        del_button = QPushButton('Delete')
        export_button = QPushButton('Export Results')
        edit_button.hide()
        del_button.hide()
        export_button.hide()
        eddel_sublayout = QHBoxLayout()
        eddel_sublayout.addWidget(export_button)
        eddel_sublayout.addWidget(del_button)
        eddel_sublayout.addWidget(edit_button)
        prev_button = QPushButton('Previous')
        next_button = QPushButton('Next')
        prev_button.clicked.connect(lambda: prev_result(self))
        next_button.clicked.connect(lambda: next_result(self, self.selectedDatabase))
        edit_button.clicked.connect(self.changeTabForEditing)
        del_button.clicked.connect(lambda: subFunctionForDeletion(self.plot_current_value,defaultDB))
        export_button.clicked.connect(lambda: self.exportSearchResults(self.selectedDatabase))

        results_table_Label = QLabel('O.R.E.O.S. Assessment of Hazard by Scale:')
        results_table = QTableWidget(1, 4)
        results_table.setHorizontalHeaderLabels(['<5 g', '5 to <100 g', '100 to 500 g', '>500 g'])
        results_table.verticalHeader().setVisible(False)
        results_table.setMaximumHeight(53)
        results_table.setMaximumWidth(402)
        results_table.setMinimumHeight(53)
        results_table.setMinimumWidth(402)
        results_approval_warning = QLabel('')

        sub_search_layout = QHBoxLayout()
        search_layout.addWidget(lbl_search)
        sub_search_layout.addWidget(self.searchTypeSelection)
        sub_search_layout.addWidget(self.searchSubType)
        sub_search_layout.addWidget(entry_search)
        sub_search_layout.addWidget(btn_search)
        search_layout.addLayout(sub_search_layout)
        search_layout.addLayout(eddel_sublayout)
        search_layout.addStretch()

        make_plot_label = QLabel('For search results: ')
        make_plot_button = QPushButton('Plot')
        self.select_x_values = QComboBox(self)
        vs_label = QLabel('<em>vs</em>')
        self.select_y_values = QComboBox(self)
        selectable_xy_values = ['Qdsc', 'Tonset', 'Tinit', 'Td24', 'IS', 'EP', 'OB', 'RoS', 'OREOS', 'MW', 'Log(Qdsc)', 'Log(Tonset-25)', 'Log(Tinit-25)']
        self.select_x_values.addItems(selectable_xy_values)
        self.select_y_values.addItems(selectable_xy_values)
        make_plot_layout = QHBoxLayout()
        make_plot_layout.addWidget(make_plot_label)
        make_plot_layout.addWidget(make_plot_button)
        make_plot_layout.addWidget(self.select_x_values)
        make_plot_layout.addWidget(vs_label)
        make_plot_layout.addWidget(self.select_y_values)
        make_plot_layout.addStretch()
        search_layout.addLayout(make_plot_layout)
        make_plot_label.hide()
        make_plot_button.hide()
        self.select_x_values.hide()
        vs_label.hide()
        self.select_y_values.hide()
        make_plot_button.clicked.connect(lambda: self.plotSearchResults(self.plot_database, self.plot_current_value))

        def subFunctionForDeletion(Entry, Database):
            self.delCurrentEntry(Entry, Database)
            search_database()

        search_tab.setLayout(search_layout)
        self.tab_widget.addTab(search_tab, "Search")

        # Tab for Settings
        settings_tab = QWidget()
        settings_layout = QVBoxLayout()
        self.config = configparser.ConfigParser()
        self.config.read('./_core/ThermalDex.ini')

        settings_intro = QLabel('<h1>ThermalDex Settings</h1><p>This pane contains the settings for this application. Change relavent settings as needed.</p>')
        settings_layout.addWidget(settings_intro)
        #settings_layout.addWidget(self.default_file_label)
        #settings_layout.addWidget(self.default_file_input)
        #settings_layout.addWidget(self.default_file_button)
        #settings_layout.addWidget(self.save_button)

        # Defaults
        defaults_section = Section("Defaults", 100, self)
        defaults_layout = QVBoxLayout()
        defaults_layout.addWidget(QLabel("Select ThermalDex Defaults", defaults_section))
        #any_layout.addWidget(QPushButton("Button in Section", section))
        QDSCDefaultsLayout = QHBoxLayout()
        QDSCDefaultsLabel = QLabel('Select Q<sub>DSC</sub> Units: ')
        self.QDSCDefaultsCombo = QComboBox()
        self.QDSCDefaultsCombo.addItems(['J g⁻¹', 'cal g⁻¹'])
        QDSCApplyButton = QPushButton('Apply')
        QDSCApplyButton.clicked.connect(self.save_units_settings)
        QDSCDefaultsLayout.addWidget(QDSCDefaultsLabel)
        QDSCDefaultsLayout.addWidget(self.QDSCDefaultsCombo)
        QDSCDefaultsLayout.addWidget(QDSCApplyButton)
        defaults_layout.addLayout(QDSCDefaultsLayout)

        yoshidaDefaultsLayout = QHBoxLayout()
        yoshidaDefaultsLabel = QLabel('Yoshia Method: ')
        self.yoshidaDefaultsCombo = QComboBox()
        self.yoshidaDefaultsCombo.addItems(['Pfizer', 'Yoshida'])
        yoshidaApplyButton = QPushButton('Apply')
        yoshidaApplyButton.clicked.connect(self.save_yoshidacal_settings)
        yoshidaDefaultsLayout.addWidget(yoshidaDefaultsLabel)
        yoshidaDefaultsLayout.addWidget(self.yoshidaDefaultsCombo)
        yoshidaDefaultsLayout.addWidget(yoshidaApplyButton)
        defaults_layout.addLayout(yoshidaDefaultsLayout)

        defaults_section.setContentLayout(defaults_layout)

        # Warnings
        warnings_section = Section("Warnings", 100, self)
        warnings_layout = QVBoxLayout()
        warnings_layout.addWidget(QLabel("Set Warning Behaviour for ThermalDex", warnings_section))
        #any_layout.addWidget(QPushButton("Button in Section", section))
        Td24WarnLayout = QHBoxLayout()
        Td24WarnLabel = QLabel('Select T<sub>D24</sub> temperature warning limits:')
        Td24WarnUpperLabel = QLabel('Upper Limit (°C) = ')
        self.amberTd24Combo = QComboBox()
        self.amberTd24Combo.addItems([str(n) for n in range(1, 500)])
        self.amberTd24Combo.setMaximumWidth(75)
        self.amberTd24Combo.setCurrentText(ambertd24limit)
        Td24WarnLowerLabel = QLabel('Lower Limit (°C) = ')
        self.amberTd24LowerCombo = QComboBox()
        self.amberTd24LowerCombo.addItems([str(n) for n in range(1, 300)])
        self.amberTd24LowerCombo.setMaximumWidth(75)
        self.amberTd24LowerCombo.setCurrentText(redtd24limit)
        Td24WarnApplyButton = QPushButton('Apply')
        Td24WarnApplyButton.setMaximumWidth(100)
        Td24WarnApplyButton.clicked.connect(self.save_td24_warn_settings)
        warnings_layout.addWidget(Td24WarnLabel)
        warnings_layout.addSpacing(10)
        Td24WarnLayout.addWidget(Td24WarnUpperLabel)
        Td24WarnLayout.addWidget(self.amberTd24Combo)
        Td24WarnLayout.addSpacing(10)
        Td24WarnLayout.addWidget(Td24WarnLowerLabel)
        Td24WarnLayout.addWidget(self.amberTd24LowerCombo)
        Td24WarnLayout.addStretch()
        Td24WarnLayout.addWidget(Td24WarnApplyButton)
        warnings_layout.addLayout(Td24WarnLayout)
        warnings_layout.addSpacing(10)


        oreoWarnLayout = QHBoxLayout()
        oreoWarnLabel = QLabel('Select O.R.E.O. warning behavour:')
        oreoWarnFirstLabel = QLabel('Warn if scale = ')
        #oreoWarnLabel.setWordWrap(True)
        oreoWarnFirstLabel.setMaximumWidth(75)
        #oreoWarnLabel.setFixedHeight(40)
        self.oreoWarnScaleCombo = QComboBox()
        self.oreoWarnScaleCombo.addItems(['<5 g', '5 to <100 g', '100 to 500 g', '>500 g'])
        self.oreoWarnScaleCombo.setCurrentText(oreohazardwarningscale)
        self.oreoWarnScaleCombo.setMaximumWidth(90)
        oreoWarnContLabel = QLabel('has hazard rating of =<br>(or greater)')
        oreoWarnContLabel.setWordWrap(True)
        oreoWarnContLabel.setFixedHeight(40)
        self.oreoWarnHazCombo = QComboBox()
        self.oreoWarnHazCombo.addItems(["Low Hazard","Medium Hazard","High Hazard"])
        self.oreoWarnHazCombo.setCurrentText(oreohazardlimit)
        self.oreoWarnHazCombo.setMaximumWidth(100)
        oreoWarnApplyButton = QPushButton('Apply')
        oreoWarnApplyButton.setMaximumWidth(100)
        oreoWarnApplyButton.clicked.connect(self.save_oreo_warn_settings)
        warnings_layout.addWidget(oreoWarnLabel)
        oreoWarnLayout.addWidget(oreoWarnFirstLabel)
        oreoWarnLayout.addWidget(self.oreoWarnScaleCombo)
        oreoWarnLayout.addSpacing(10)
        oreoWarnLayout.addWidget(oreoWarnContLabel)
        oreoWarnLayout.addWidget(self.oreoWarnHazCombo)
        oreoWarnLayout.addStretch()
        oreoWarnLayout.addWidget(oreoWarnApplyButton)
        warnings_layout.addLayout(oreoWarnLayout)
        warnings_layout.addSpacing(10)

        warnings_section.setContentLayout(warnings_layout)       


        # Core Settings
        core_section = Section("Core Settings", 100, self)
        core_layout = QVBoxLayout()
        core_layout.addWidget(QLabel("Core Settings for ThermalDex (be careful...)", core_section))
        
        databaseCoreLayout = QHBoxLayout()
        databaseCoreLabel = QLabel("Compound Database:")
        self.databaseCoreInput = QLineEdit(self.config.get('Database', 'defaultDB'))
        databaseCoreSelectButton = QPushButton("Browse")
        databaseCoreSelectButton.clicked.connect(self.select_database)
        databaseCoreApplyButton = QPushButton("Apply")
        databaseCoreApplyButton.clicked.connect(self.save_database_settings)
        databaseCoreLayout.addWidget(databaseCoreLabel)
        databaseCoreLayout.addSpacing(14)
        databaseCoreLayout.addWidget(self.databaseCoreInput)
        databaseCoreLayout.addWidget(databaseCoreSelectButton)
        databaseCoreLayout.addWidget(databaseCoreApplyButton)
        core_layout.addLayout(databaseCoreLayout)

        hegCoreLayout = QHBoxLayout()
        hegCoreLabel = QLabel("High Energy Groups List:")
        self.hegCoreInput = QLineEdit(self.config.get('Lists of Groups', 'highenergygroups'))
        hegCoreSelectButton = QPushButton("Browse")
        hegCoreSelectButton.clicked.connect(self.select_heglist)
        hegCoreApplyButton = QPushButton("Apply")
        hegCoreApplyButton.clicked.connect(self.save_heglist_settings)
        hegCoreLayout.addWidget(hegCoreLabel)
        hegCoreLayout.addWidget(self.hegCoreInput)
        hegCoreLayout.addWidget(hegCoreSelectButton)
        hegCoreLayout.addWidget(hegCoreApplyButton)
        core_layout.addLayout(hegCoreLayout)

        efgCoreLayout = QHBoxLayout()
        efgCoreLabel = QLabel("Explosive Groups List:")
        self.efgCoreInput = QLineEdit(self.config.get('Lists of Groups', 'expenergygroups'))
        efgCoreSelectButton = QPushButton("Browse")
        efgCoreSelectButton.clicked.connect(self.select_efglist)
        efgCoreApplyButton = QPushButton("Apply")
        efgCoreApplyButton.clicked.connect(self.save_efglist_settings)
        efgCoreLayout.addWidget(efgCoreLabel)
        efgCoreLayout.addSpacing(14)
        efgCoreLayout.addWidget(self.efgCoreInput)
        efgCoreLayout.addWidget(efgCoreSelectButton)
        efgCoreLayout.addWidget(efgCoreApplyButton)
        core_layout.addLayout(efgCoreLayout)

        core_section.setContentLayout(core_layout)

        # Apperance Settings
        apperance_section = Section("Apperance Settings", 100, self)
        apperance_layout = QVBoxLayout()
        apperance_layout.addWidget(QLabel("Apperance Settings for ThermalDex (currently inactive)", core_section))
        
        palletApperLayout = QHBoxLayout()
        palletApperLabel = QLabel('Molecule Drawing Pallet:')
        palletApperCombo = QComboBox()
        palletApperCombo.addItems(['Colorful', 'Black and White'])
        palletApperButton = QPushButton('Apply')
        palletApperLayout.addWidget(palletApperLabel)
        palletApperLayout.addWidget(palletApperCombo)
        palletApperLayout.addWidget(palletApperButton)
        apperance_layout.addLayout(palletApperLayout)

        themeApperLayout = QHBoxLayout()
        themeApperLabel = QLabel('ThermalDex Theme:')
        themeApperCombo = QComboBox()
        themeApperCombo.addItems(['Light Mode', 'Dark Mode'])
        themeApperButton = QPushButton('Apply')
        themeApperLayout.addWidget(themeApperLabel)
        themeApperLayout.addWidget(themeApperCombo)
        themeApperLayout.addWidget(themeApperButton)
        apperance_layout.addLayout(themeApperLayout)

        fullscreenApperLayout = QHBoxLayout()
        fullscreenApperLabel = QLabel('Full Screen on Start-up?:')
        fullscreenApperCombo = QComboBox()
        fullscreenApperCombo.addItems(['True', 'False'])
        fullscreenApperButton = QPushButton('Apply')
        fullscreenApperLayout.addWidget(fullscreenApperLabel)
        fullscreenApperLayout.addWidget(fullscreenApperCombo)
        fullscreenApperLayout.addWidget(fullscreenApperButton)
        apperance_layout.addLayout(fullscreenApperLayout)
        
        apperance_section.setContentLayout(apperance_layout)       

        settings_layout.addWidget(defaults_section)
        settings_layout.addWidget(warnings_section)
        settings_layout.addWidget(core_section)
        settings_layout.addWidget(apperance_section)

        # User Reload of Config
        reload_sublayout = QHBoxLayout()
        reload_config_label = QLabel('Configuration changes will take effect the next time ThermalDex is opened.<br> To force changes to take effect now, click:')
        reload_config_button = QPushButton('Reload Config')
        reload_config_button.clicked.connect(self.reload_config_func)
        reload_sublayout.addWidget(reload_config_label)
        reload_sublayout.addWidget(reload_config_button)
        settings_layout.addLayout(reload_sublayout)

        # Import Database -> a merge of current database with existing file
        import_sublayout = QHBoxLayout()
        import_heading = QLabel('<h2>Bulk Import</h2>')
        import_db_label = QLabel('To bulk import molecules into ThermalDex click:')
        import_explain_label = QLabel('<p>Note these compounds will be added directly to your database. You will be warned for any molecule which would overwrite an existing entry in you database. Please ensure the correct .csv format is used. To import molecules with no other imformation availible, simlpy create a .csv with a column heading of "SMILES" and a new line between the SMILES of the molecules.</p>')
        #import_explain_label.setMaximumWidth(400)
        import_explain_label.setWordWrap(True)
        import_db_button = QPushButton('Import Database File')
        import_db_button.setMaximumWidth(190)
        import_db_button.clicked.connect(self.import_external_database)
        import_sublayout.addWidget(import_db_label)
        import_sublayout.addWidget(import_db_button)
        settings_layout.addSpacing(40)
        settings_layout.addWidget(import_heading)
        settings_layout.addLayout(import_sublayout)
        settings_layout.addWidget(import_explain_label)

        settings_layout.addStretch()

        settings_tab.setLayout(settings_layout)
        self.tab_widget.addTab(settings_tab, "Settings")

        # Tab for "About" information

        #def hover(url):
        #        QToolTip.showText(QCursor.pos(), titles.get(url, url))
        #    else:
        #        QToolTip.hideText()

        about_tab = QWidget()
        about_layout = QVBoxLayout()
        about_title = QLabel("<b>About ThermalDex</b>\n\n")
        about_blank = QLabel("\nVersion: " + versionNumber + "\n")
        about_text = QLabel("\n\nThis is a simple tool for assessing and recording the potential thermal hazards assoicated with a molecule. It uses the <b>'O.R.E.O.S.'</b> assement scale and other ideas that can be read about in <a href=\"https://pubs.acs.org/doi/10.1021/acs.oprd.0c00467\"><em>Org. Process Res. Dev.</em> 2021, 25, 2, 212–224</a> by Jeffrey B. Sperry et. al.")
        iconLabel = QLabel()
        iconImage = QPixmap(".\\_core\\ThermalDexIcon.jpg")
        scaledIcon = iconImage.scaled(400, 400, Qt.KeepAspectRatio)
        iconLabel.setText("test") #.setPixmap(scaledIcon)
        
        scene = QGraphicsScene()
        pic = QGraphicsPixmapItem()
        pic.setPixmap(scaledIcon) #QPixmap.fromImage(iconImage))
        scene.setSceneRect(0, 0, 400, 400)
        scene.addItem(pic)
        view = QGraphicsView()
        view.setScene(scene)
        view.setStyleSheet("border: 0px")
        altNames_text = QLabel("The alternative name for this project is <b>C.O.O.K.I.E.S.</b> which stands for <em>C</em>alculation <em>O</em>f <em>O</em>.R.E.O.S., <em>K</em>illing <em>I</em>nelegant <em>E</em>xcel <em>S</em>olutions. This would be the offical relase name, but it felt a bit glib.")
        moreabout_text = QLabel("This tool has been developed by Tom Sheridan and Matt Mulheir using PyQt5 and RDKit. You can contribute to this project on <a href=\"https://github.com/Tomisall/ProcessSafteyDB\">GitHub</a>.")
        #about_text.setAlignment(Qt.AlignCenter)
        about_layout.addWidget(about_title)
        about_layout.addWidget(about_text)
        about_layout.addWidget(about_blank)
        about_layout.addWidget(view) #iconLabel)
        about_layout.addWidget(altNames_text)
        about_layout.addWidget(moreabout_text)
        about_text.setTextFormat(Qt.RichText)
        about_text.setOpenExternalLinks(True)
        about_text.setTextInteractionFlags(Qt.LinksAccessibleByMouse)
        #about_text.linkHovered.connect(hover)
        about_text.setWordWrap(True) 
        altNames_text.setWordWrap(True) 
        moreabout_text.setTextInteractionFlags(Qt.LinksAccessibleByMouse)
        moreabout_text.setOpenExternalLinks(True)
        moreabout_text.setWordWrap(True)

        #abouttestLabel = QLabel("<a href=\"http://www.google.com\">'Click this link to go to Google'</a>")
        #about_layout.addWidget(abouttestLabel)
        #abouttestLabel.setOpenExternalLinks(True)
        #abouttestLabel.setTextInteractionFlags(Qt.LinksAccessibleByMouse)

        about_tab.setLayout(about_layout)
        self.tab_widget.addTab(about_tab, "About")
        layout.addWidget(self.tab_widget)



        # Set up the main window
        self.setLayout(layout)
        self.setGeometry(100, 100, 600, 895)
        self.setWindowTitle('ThermalDex')

    def reload_config_func(self):
        #sys.exit(app.exec_())
        #QCoreApplication.quit()
        global defaultDB
        global highEnergyGroups
        global expEnergyGroups
        global yoshidaMethod
        global qdscUnits
        global ambertd24limit
        global redtd24limit
        global oreohazardlimit
        global oreohazardwarningscale

        defaultDB, highEnergyGroups, expEnergyGroups, yoshidaMethod, qdscUnits, ambertd24limit, redtd24limit, oreohazardlimit, oreohazardwarningscale = altImportConfig()
        QMessageBox.information(self, "Config Settings", "Settings have been loaded successfully.")

    def exportSearchResults(self, Database):
        savefile, _ = QFileDialog.getSaveFileName(self, "Export Search Results as CSV", "ThermalDexResults.csv", "CSV Files (*.csv)")
        if savefile:
            Database.to_csv(savefile, index=False)

    def delCurrentEntry(self, currentResult, Database):
        errorInfo = "Really? Delete this entry? Are you sure?"
        userInteract = self.interactiveErrorMessage(errorInfo) 
        if userInteract == QMessageBox.Yes:
            print(f'dict Results: {currentResult}')
            storedData = pd.read_csv(Database, index_col=0)
            removeDataTranspose = pd.DataFrame.from_dict(currentResult, orient='index')
            removeDataIndexed = removeDataTranspose.T
            removeData = removeDataIndexed.set_index('SMILES')
            print(f'\n\n\nDatabase:\n{storedData}\n\nValue to Delete:\n{removeData}')

            reducedData = storedData.drop([currentResult['SMILES']])
            print(f'\n\nDatabase after drop would be:\n{reducedData}')
            reducedData['SMILES'] = reducedData.index
            reducedData = reducedData[ ['SMILES'] + [ col for col in reducedData.columns if col != 'SMILES' ] ]
            print(reducedData)
            reducedData.to_csv(Database, index=False)


    def plotSearchResults(self, Results, currentResult):
        print(f'\n\nDatabase:\n{Results}\n\nCurrent:\n{currentResult}')

        selectable_xy_values = ['Qdsc', 'Tonset', 'Tinit', 'Td24', 'IS', 'EP', 'OB', 'RoS', 'OREOS', 'MW', 'Log(Qdsc)', 'Log(Tonset-25)', 'Log(Tinit-25)']
        databaseHeadings = ['Q_dsc', 'onsetT', 'initT', 'Td24', 'IS_val', 'EP_val', 'OB_val', 'RoS_val', 'oreoLargeScale_val', 'MW']
 
        print(f'\n\nxvalues = {self.select_x_values.currentText()}')
        print(f'yvalues = {self.select_y_values.currentText()}')

        if self.select_x_values.currentText() == 'Log(Qdsc)':
            rawListValues = Results['Q_dsc'].tolist()
            x_values = log10(rawListValues)
            current_x = log10(currentResult['Q_dsc'])
        elif self.select_x_values.currentText() == 'Log(Tonset-25)':
            rawListValues = Results['onsetT'].tolist()
            subtractedValues = [x - 25 for x in rawListValues]
            x_values = log10(subtractedValues)
            current_x = log10(currentResult['onsetT']-25)
        elif self.select_x_values.currentText() == 'Log(Tinit-25)':
            rawListValues = Results['initT'].tolist()
            subtractedValues = [x - 25 for x in rawListValues]
            x_values = log10(subtractedValues)
            current_x = log10(currentResult['initT']-25)
        else:
            selected_x_index = selectable_xy_values.index(self.select_x_values.currentText())
            databaseCol = databaseHeadings[selected_x_index]
            x_values = Results[databaseCol].tolist()
            current_x = currentResult[databaseCol]

        if self.select_y_values.currentText() == 'Log(Qdsc)':
            rawListValues = Results['Q_dsc'].tolist()
            y_values = log10(rawListValues)
            current_y = log10(currentResult['Q_dsc'])
        elif self.select_y_values.currentText() == 'Log(Tonset-25)':
            rawListValues = Results['onsetT'].tolist()
            subtractedValues = [y - 25 for y in rawListValues]
            y_values = log10(subtractedValues)
            current_y = log10(currentResult['onsetT']-25)
        elif self.select_y_values.currentText() == 'Log(Tinit-25)':
            rawListValues = Results['initT'].tolist()
            subtractedValues = [y - 25 for y in rawListValues]
            y_values = log10(subtractedValues)
            current_y = log10(currentResult['initT']-25)
        else:
            selected_y_index = selectable_xy_values.index(self.select_y_values.currentText())
            databaseCol = databaseHeadings[selected_y_index]
            y_values = Results[databaseCol].tolist()
            current_y = currentResult[databaseCol]

        scatterSelection(x_values,y_values,self.select_x_values.currentText(),self.select_y_values.currentText(),current_x,current_y)


    def set_sub_search(self, subIndex):
        self.searchSubType.clear()
        self.searchSubType.addItems(self.listOfSearchTypes[subIndex])

    def import_external_database(self):
        options = QFileDialog.Options()
        import_file, _ = QFileDialog.getOpenFileName(self, "Select Database to Import:", "", "CSV Files (*.csv)", options=options)
        if import_file:
                  importedDB = pd.read_csv(import_file)
                  for index, row in importedDB.iterrows():
                    dictRow = row.to_dict()
                    readMolecule = thermalDexMolecule(**dictRow)
                    self.writeToDatabase(readMolecule, defaultDB)
                
                  QMessageBox.information(self, "Database Import", "Import has completed.")
                  
            

    def select_database(self):
        options = QFileDialog.Options()
        default_file, _ = QFileDialog.getOpenFileName(self, "Select Database to Use:", "", "CSV Files (*.csv)", options=options)
        if default_file:
            rel_file = path.relpath(default_file)
            self.databaseCoreInput.setText(rel_file)

    def save_oreo_warn_settings(self):
        self.config.set('Warnings', 'oreohazardwarningscale', self.oreoWarnScaleCombo.currentText())
        self.config.set('Warnings', 'oreohazardlimit', self.oreoWarnHazCombo.currentText())

        with open('./_core/ThermalDex.ini', 'w') as configfile:
            self.config.write(configfile)

    def save_td24_warn_settings(self):
        if int(self.amberTd24Combo.currentText()) > int(self.amberTd24LowerCombo.currentText()):
            self.config.set('Warnings', 'ambertd24limit', self.amberTd24Combo.currentText())
            self.config.set('Warnings', 'redtd24limit', self.amberTd24LowerCombo.currentText())

            with open('./_core/ThermalDex.ini', 'w') as configfile:
                self.config.write(configfile)

            QMessageBox.information(self, "Settings Saved", "Settings have been saved.")
        else:
            self.interactiveErrorMessage('T<sub>D24</sub> lower limit must be smaller than the upper limit!')

    def save_yoshidacal_settings(self):
        self.config.set('Calculations', 'yoshidacalcs', self.yoshidaDefaultsCombo.currentText())

        with open('./_core/ThermalDex.ini', 'w') as configfile:
            self.config.write(configfile)

        QMessageBox.information(self, "Settings Saved", "Settings have been saved.")

    def save_units_settings(self):
        self.config.set('Default Values', 'qdscunits', str(self.QDSCDefaultsCombo.currentIndex()))

        with open('./_core/ThermalDex.ini', 'w') as configfile:
            self.config.write(configfile)

        QMessageBox.information(self, "Settings Saved", "Settings have been saved.")

    def save_database_settings(self):
        self.config.set('Database', 'defaultDB', self.databaseCoreInput.text())

        with open('./_core/ThermalDex.ini', 'w') as configfile:
            self.config.write(configfile)
            
        QMessageBox.information(self, "Settings Saved", "Settings have been saved.")

    def select_heglist(self):
        options = QFileDialog.Options()
        default_file, _ = QFileDialog.getOpenFileName(self, "Select List to Use:", "", "CSV Files (*.csv)", options=options)
        if default_file:
            rel_file = path.relpath(default_file)
            self.hegCoreInput.setText(rel_file)

    def save_heglist_settings(self):
        self.config.set('Lists of Groups', 'highenergygroups', self.hegCoreInput.text())

        with open('./_core/ThermalDex.ini', 'w') as configfile:
            self.config.write(configfile)
            
        QMessageBox.information(self, "Settings Saved", "Settings have been saved.")

    def select_efglist(self):
        options = QFileDialog.Options()
        default_file, _ = QFileDialog.getOpenFileName(self, "Select List to Use:", "", "CSV Files (*.csv)", options=options)
        if default_file:
            rel_file = path.relpath(default_file)
            self.efgCoreInput.setText(rel_file)

    def save_efglist_settings(self):
        self.config.set('Lists of Groups', 'expenergygroups', self.efgCoreInput.text())

        with open('./_core/ThermalDex.ini', 'w') as configfile:
            self.config.write(configfile)
            
        QMessageBox.information(self, "Settings Saved", "Settings have been saved.")


    def changeTabForEditing(self):
        #try:
        self.approval_needed.hide()
        editDB = self.selectedDatabase.fillna('')
        current_row = editDB.iloc[self.current_index]
        dictRow = current_row.to_dict()
        print(dictRow)
        print('\n\n')
        readMolecule = thermalDexMolecule(**dictRow)
        readMolecule.genMol()
        # Make Pixmap Image to Display.
        pixmap = readMolecule.molToQPixmap()
        scaledPixmap = pixmap #.scaled(550, 275, Qt.KeepAspectRatio)
        scene = QGraphicsScene()
        #scene.setSceneRect(0, 0, 400, 400)
        scene.addPixmap(scaledPixmap) #pixmap)
        self.mol_display.setScene(scene)
        readMolecule.molPixmap = None
        self.smiles_input.setText(readMolecule.SMILES)
        self.name_input.setText(readMolecule.name)
        self.mp_input.setText(str(readMolecule.mp))
        self.mpEnd_input.setText(str(readMolecule.mpEnd))
        self.Qdsc_input.setText(str(readMolecule.Q_dsc))
        comboIndex = self.QUnitsSelection.findText(readMolecule.Qunits)
        self.QUnitsSelection.setCurrentIndex(comboIndex)
        self.TE_input.setText(str(readMolecule.onsetT))
        self.Tinit_input.setText(str(readMolecule.initT))
        self.proj_input.setText(str(readMolecule.proj))
        niceMWStr = "{:.2f}".format(readMolecule.MW)
        self.mwLabel.setText('MW: ' + niceMWStr) #str(readMolecule.MW))
        self.HEGlabel.setText('Number of High Energy Groups: ' + str(readMolecule.HEG))
        self.EFGlabel.setText('Number of Explosive Functional Groups: ' + str(readMolecule.EFG)) 
        self.eleLabel.setText('Elemental Composition: ' + readMolecule.eleComp)
        self.RoSLabel.setText('Rule of Six: ' + str(readMolecule.RoS_val) + readMolecule.RoS_des)
        self.obLabel.setText('Oxygen Balance: ' + '{:.2f}'.format(readMolecule.OB_val) + ' ' + readMolecule.OB_des)
        self.table.clearContents()
        hazardList = [readMolecule.oreoSmallScale_des, readMolecule.oreoTensScale_des, readMolecule.oreoHundsScale_des, readMolecule.oreoLargeScale_des]
        # Add values to cells
        for i, hazardClass in enumerate(hazardList):
            item = QTableWidgetItem(hazardClass)
            self.table.setItem(0, i, item)

            # Color code cells based on values
            classColor = self.getColorForValue(hazardClass)
            print(classColor)
            item.setBackground(classColor)

        if readMolecule.hammerDrop == 'Positive':
            self.hamSelection.setCurrentIndex(1)

        elif readMolecule.hammerDrop == 'Negative':
            self.hamSelection.setCurrentIndex(2)

        else:
            self.hamSelection.setCurrentIndex(0)

        if readMolecule.friction == 'Positive':
            self.fricSelection.setCurrentIndex(1)

        elif readMolecule.friction == 'Negative':
            self.fricSelection.setCurrentIndex(2)

        else:
            self.fricSelection.setCurrentIndex(0)

        fileCounter = self.countFiles(defaultDB)
        self.filesCount.setText(f"{str(fileCounter)} Attached Files")
        self.attachedFilesLabel.show()
        self.filesCount.show()
        self.attach_button.show()

        if self.check_if_oreos_need_approval(readMolecule) == 'Show Approval Message':
            self.approval_needed.show()
    
        try:
            isStr = "{:.2f}".format(readMolecule.IS_val)
        except:
            isStr = ''
        self.ISLabel.setText('Yoshida Impact Sensitivity: ' + isStr + readMolecule.IS_des)
        try:
            epStr = "{:.2f}".format(readMolecule.EP_val)
        except:
            epStr = ''
        self.EPLabel.setText('Yoshida Explosive Propagation: ' + epStr + readMolecule.EP_des)
        try:
            d24Str = "{:.1f}".format(readMolecule.Td24)
            if int(ambertd24limit) >= readMolecule.Td24 > int(redtd24limit):
                self.Td24Label.setText(f"T<sub>D24</sub>: <b style='color: orange;'> {d24Str} °C</b>")
            elif readMolecule.Td24 <= int(redtd24limit):
                self.Td24Label.setText(f"<h3 style='color: red;'>T<sub>D24</sub>: {d24Str} °C</h3>")
                self.approval_needed.show()
            else:
                self.Td24Label.setText('T<sub>D24</sub>: <b>' + d24Str + ' °C</b>')       
        except:
            pass
        self.tab_widget.setCurrentWidget(self.molecule_tab)   #0) #self.tab_widget.findChild(QWidget, "Add"))

        #except:
        #    window.showErrorMessage("Editing current selection. Has a molecule been found?")

    def mrvToSMILES(self):
        try:
            copyMolXML = pyperclip.paste()
            copyMolReal = rdmolfiles.MolFromMrvBlock(copyMolXML)
            smilesFromMarvin = rdmolfiles.MolToSmiles(copyMolReal)
            print(smilesFromMarvin)
            pyperclip.copy(smilesFromMarvin)

        except:
            print("No mrv XML found.")

    def openFileManager(self):
        self.fileWindow = FileManagementWindow(self.smiles_input.text(), defaultDB, window)
        self.fileWindow.show()
        self.fileWindow.raise_()
        self.fileWindow.activateWindow()

    def findFolder(self, database):
        try:
            storedData = pd.read_csv(database)
            row_index = storedData.index[storedData['SMILES'] == self.smiles_input.text()].tolist()
            folderInfo = storedData['dataFolder'][row_index[0]]
            return folderInfo

        except:
            return ''
     

    def countFiles(self, database):
        try:
            folderInfo = self.findFolder(database)
            files = [f for f in listdir(folderInfo) if isfile(join(folderInfo, f))]
            fileCount = len(files)
            return fileCount
        except:
            print('No Folder Given.')
            return 0

    def showErrorMessage(self, errorCode):
        self.msg = QMessageBox()
        self.msg.setIcon(QMessageBox.Warning)
        self.msg.setText("An error has occured")
        self.msg.setInformativeText("The source of the problem seems to be with: " + errorCode + "\nTry again and contact developer if the problem persists.\n\n Select 'OK' to export an error log or 'Cancel' to dismiss this message.")
        self.msg.setWindowTitle("ThermalDex Error")
        self.msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        returnValue = self.msg.exec()

        if returnValue == QMessageBox.Ok:
            self.exportErrorLog()

    def exportErrorLog(self):
        saveLog, _ = QFileDialog.getSaveFileName(self, "Export Error Log", "ThermalDex.log", "Text Files (*.log)")
        if saveLog:
            copy2('./_core/ThermalDex.log',saveLog)

    def interactiveErrorMessage(self, errorInfo):
        self.interactMsg = QMessageBox()
        self.interactMsg.setIcon(QMessageBox.Information)
        self.interactMsg.setText("Action Needed")
        self.interactMsg.setInformativeText(errorInfo)
        self.interactMsg.setWindowTitle("ThermalDex - Info Box")
        self.interactMsg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        returnValue = self.interactMsg.exec()
        return returnValue

    def errorTestingTool(self):
        self.showErrorMessage("This is an Alternative test method for error handling")

    def clearTheCalcdValues(self):
        scene = QGraphicsScene()
        self.approval_needed.hide()
        self.mol_display.setScene(scene)
        self.mwLabel.setText(self.mwText)
        self.HEGlabel.setText(self.HEGText)
        self.EFGlabel.setText(self.EFGText)
        self.eleLabel.setText(self.eleText)
        self.RoSLabel.setText(self.RoSText)
        self.obLabel.setText(self.obText)
        self.ISLabel.setText(self.ISText)
        self.EPLabel.setText(self.EPText)
        self.Td24Label.setText(self.Td24Text)
        clearEntry = ['', '', '', '']
        for i, entry in enumerate(clearEntry):
            clear = QTableWidgetItem(entry)
            self.table.setItem(0, i, clear)
            clear.setBackground(QColor(255, 255, 255))

    def clearUserValues(self):
        self.smiles_input.setText('')
        self.name_input.setText('')
        self.mp_input.setText('')
        self.mpEnd_input.setText('')
        self.Qdsc_input.setText('')
        self.QUnitsSelection.setCurrentIndex(int(qdscUnits))
        self.TE_input.setText('')
        self.Tinit_input.setText('')
        self.proj_input.setText('')
        self.hamSelection.setCurrentIndex(0)
        self.fricSelection.setCurrentIndex(0)
        self.attachedFilesLabel.hide()
        self.filesCount.hide()
        self.attach_button.hide()

    def resetToDefaultState(self):
        self.clearTheCalcdValues()
        self.clearUserValues()
        layout = self.layout()
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)

    def genCoreValuesFromMol(self, molecule):
        try:
            molecule.mwFromMol()
        except:
            window.showErrorMessage("Calculating MW from RDChem Mol.")

        try:
            molecule.HEGFromMol(highEnergyGroups)
        except:
            window.showErrorMessage("Determining High Energy Groups from RDChem Mol.")

        try:
           molecule.EFGFromMol(expEnergyGroups)
        except:
            window.showErrorMessage("Determining Explosive Groups from RDChem Mol.")

        try:
           molecule.eleCompFromMol()
        except:
            window.showErrorMessage("Determining Elemental Compostion from RDChem Mol.")
        
        try:
           molecule.CHOFromEleComp()
        except:
            window.showErrorMessage("Determining CHO from Elemental Compostion.")

        try:        
           molecule.RoSFromEleComp()
        except:
            window.showErrorMessage("Calculating Rule of Six from Elemental Compostion.")

        try:
           molecule.OBFromEleComp()
        except:
            window.showErrorMessage("Calculating Oxygen Balance from Elemental Compostion.")

        try:
           molecule.oreoOnsetTadjustment()
        except:
            window.showErrorMessage("Adjusting O.R.E.O.S. Calculation for Onset Temperature")

        try:
           molecule.oreoSafeScaleCal()
        except:
            window.showErrorMessage("Determining O.R.E.O.S. Hazard by Scale.")

        try:
           molecule.Td24FromThermalProps()
        except:
            window.showErrorMessage("Calculating T<sub>D24</sub> from Thermal Properties.")

        try:
           molecule.yoshidaFromThermalProps()
        except:
            window.showErrorMessage("Calculating Yoshida values from Thermal Properties.")

        try:
            molecule.makeFolderForMolData()
        except:
            window.showErrorMessage("Making Folder to Hold Molecule Data Files.")

    def createReport(self):
        try:
            moleculeData = self.render_molecule()
            #img = moleculeData.molToIMG()
            #results = asdict(moleculeData)
            #create_pdf(results["name"], results, img) #results["molIMG"]) 
            imageData = moleculeData.molToBytes()
            dataURL = 'data:image/png;base64,' + imageData
            mdReportCreation(moleculeData, dataURL, int(ambertd24limit), int(redtd24limit))
            if self.error_flag is not None:
                self.error_message.setText('')
                layout.removeWidget(self.error_message)
        except:
            window.showErrorMessage("Generating Memo PDF from given values.")

    def writeToDatabase(self, molecule, Database):
        #molecule.genAdditionalValues()
        selectedMolData = cleanMolDataFrame(molecule)
        try:
            storedData = pd.read_csv(Database, index_col=0)
            checkData = pd.read_csv(Database)
        except:
            data = {'SMILES': []}
            storedData = pd.DataFrame(data)
            checkData = pd.DataFrame(data)
        
        print('\n\n\n')
        if selectedMolData['SMILES'][0] in storedData.index:
            print('found')
            userInteract = self.interactiveErrorMessage(f'{selectedMolData["SMILES"][0]}\n\nMolecule Already in Database. Would you like to overwrite it?')
            if userInteract == QMessageBox.Yes:
                storedData.update(selectedMolData)
                outputData = storedData
                row_index = checkData.index[checkData['SMILES'] == selectedMolData['SMILES'][0]].tolist()
                print(row_index)
                folderInfo = checkData['dataFolder'][row_index[0]]
                print(folderInfo)
                output_index = outputData.index[checkData['SMILES'] == selectedMolData['SMILES'][0]].tolist()
                outputData['dataFolder'][output_index[0]] = folderInfo

            elif userInteract == QMessageBox.No:
                outputData = storedData

            else:
                outputData = storedData

        else:
            outputData = pd.concat([storedData, selectedMolData])

        outputData['SMILES'] = outputData.index
        outputData = outputData[ ['SMILES'] + [ col for col in outputData.columns if col != 'SMILES' ] ]
        print(outputData)
        outputData.to_csv(Database, index=False)

    def getColorForValue(self, hazardClass):
        # Color-coding logic
        if hazardClass == 'High Hazard':
            return QColor(255, 0, 0)  # Red
        elif hazardClass == 'Medium Hazard':
            return QColor(255, 255, 0)  # Yellow
        elif hazardClass == 'Low Hazard':
            return QColor(0, 255, 0)  # Green
        else:
            return QColor(0, 0, 255)  # Blue

    def genMoleculeFromUserInput(self):
        # Get the SMILES string from the input field
        smiles = self.smiles_input.text()
        name = self.name_input.text()
        mp = self.mp_input.text()
        mpEnd = self.mpEnd_input.text()
        Qdsc = self.Qdsc_input.text()
        QUnits = self.QUnitsSelection.currentText()
        TE = self.TE_input.text()
        Tinit = self.Tinit_input.text()
        proj = self.proj_input.text()
        hammerDrop = self.hamSelection.currentText()
        friction = self.fricSelection.currentText()

        dataFolder = self.findFolder(defaultDB)

        # Create an RDKit molecule from the SMILES string
        addedMolecule = thermalDexMolecule(SMILES=smiles, name=name, mp=mp, mpEnd=mpEnd, Q_dsc=Qdsc, Qunits=QUnits, onsetT=TE, initT=Tinit, proj=proj, hammerDrop=hammerDrop, friction=friction, dataFolder=dataFolder, yoshidaMethod=yoshidaMethod)
        return addedMolecule
    
    def check_if_oreos_need_approval(self, molecule):
            
            print(f'oreohazardlimit: {oreohazardlimit}')
            if oreohazardlimit == 'High Hazard':
                checkOreoHazards = ['High Hazard']
            
            elif oreohazardlimit == 'Medium Hazard':
                checkOreoHazards = ['High Hazard', 'Medium Hazard']

            elif oreohazardlimit == 'Low Hazard':
                checkOreoHazards = ['High Hazard', 'Medium Hazard', 'Low Hazard']

            else:
                self.showErrorMessage('Error Reading Config Warning Settings.')
            
            print(f'oreohazardwarningscale: {oreohazardwarningscale}')
            if oreohazardwarningscale == '<5 g':
                for hazardRating in checkOreoHazards:
                    if molecule.oreoSmallScale_des == hazardRating:
                        #self.approval_needed.show()
                        return('Show Approval Message')

            elif oreohazardwarningscale == '5 to <100 g':
                for hazardRating in checkOreoHazards:
                    if molecule.oreoTensScale_des == hazardRating:
                        #self.approval_needed.show()
                        return('Show Approval Message')

            elif oreohazardwarningscale == '100 to 500 g':
                for hazardRating in checkOreoHazards:
                    if molecule.oreoHundsScale_des == hazardRating:
                        #self.approval_needed.show()
                        return('Show Approval Message')

            elif oreohazardwarningscale == '>500 g':
                for hazardRating in checkOreoHazards:
                    if molecule.oreoLargeScale_des == hazardRating:
                        #self.approval_needed.show()
                        return('Show Approval Message')

            else:
                self.showErrorMessage('Error Reading Config Warning Settings.')


    def checkIfSMILESAreValid(self, molecule):
        if molecule.SMILES != '' and molecule.SMILES is not None:
            molecule.mol = MolFromSmiles(molecule.SMILES)

        else:
            molecule.mol = None

    def displayTheMolecule(self, molecule, display):
        # Make Pixmap Image to Display.
        pixmap = molecule.molToQPixmap()
        scene = QGraphicsScene()
        scene.addPixmap(pixmap)
        display.setScene(scene)
        molecule.molPixmap = None

    def render_molecule(self):
        layout = self.layout()
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)

        # Generate a thermalDexMolecule from the user input.
        addedMolecule = self.genMoleculeFromUserInput()

        # Check that the user provided SMILES are valid.
        self.checkIfSMILESAreValid(addedMolecule)
        
        if addedMolecule.mol is not None:
            # Make Pixmap Image to Display.
            self.clearTheCalcdValues()
            self.displayTheMolecule(addedMolecule, self.mol_display)

            # Calculate Core Properties
            self.genCoreValuesFromMol(addedMolecule)

            # Format and Display Properties
            self.mwLabel.setText('MW: ' + addedMolecule.mwStr)
            self.HEGlabel.setText('Number of High Energy Groups: ' + str(addedMolecule.HEG))
            self.EFGlabel.setText('Number of Explosive Functional Groups: ' + str(addedMolecule.EFG)) 
            self.eleLabel.setText('Elemental Composition: ' + addedMolecule.eleComp)
            self.RoSLabel.setText('Rule of Six: ' + str(addedMolecule.RoS_val) + addedMolecule.RoS_des)
            self.obLabel.setText('Oxygen Balance: ' + addedMolecule.obStr + ' ' + addedMolecule.OB_des)
            self.table.clearContents()
            hazardList = [addedMolecule.oreoSmallScale_des, addedMolecule.oreoTensScale_des, addedMolecule.oreoHundsScale_des, addedMolecule.oreoLargeScale_des]
            # Add values to cells
            for i, hazardClass in enumerate(hazardList):
                item = QTableWidgetItem(hazardClass)
                self.table.setItem(0, i, item)

                # Color code cells based on values
                classColor = self.getColorForValue(hazardClass)
                print(classColor)
                item.setBackground(classColor)

            if self.check_if_oreos_need_approval(addedMolecule) == 'Show Approval Message':
                self.approval_needed.show()

            if addedMolecule.isStr != None:
                self.ISLabel.setText('Yoshida Impact Sensitivity: ' + addedMolecule.isStr + addedMolecule.IS_des)
            if addedMolecule.epStr != None:
                self.EPLabel.setText('Yoshida Explosive Propagation: ' + addedMolecule.epStr + addedMolecule.EP_des)
            if addedMolecule.Td24 != '' and addedMolecule.Td24 != 'nan' and addedMolecule.Td24 != None:
                d24Str = "{:.1f}".format(addedMolecule.Td24)
                if int(ambertd24limit) >= addedMolecule.Td24 > int(redtd24limit):
                    self.Td24Label.setText(f"T<sub>D24</sub>: <b style='color: orange;'> {d24Str} °C</b>")
                elif addedMolecule.Td24 <= int(redtd24limit):
                    self.Td24Label.setText(f"<h3 style='color: red;'>T<sub>D24</sub>: {d24Str} °C</h3>")
                    self.approval_needed.show()
                else:
                    self.Td24Label.setText('T<sub>D24</sub>: <b>' + d24Str + ' °C</b>')
                
            if addedMolecule.onsetT != 'nan' and addedMolecule.onsetT != '' and addedMolecule.onsetT != None and addedMolecule.onsetT <= 24.99:
                self.error_message = QLabel('Sub-Ambient Onset Temperature? Are you sure? Yoshida values will not be accurate.')
                layout.addWidget(self.error_message)
                self.error_flag = 201
            if addedMolecule.initT != 'nan' and addedMolecule.initT != '' and addedMolecule.initT != None and addedMolecule.initT <= 24.99:
                self.error_message = QLabel('Sub-Ambient Initiation Temperature? Are you sure? Yoshida values will not be accurate.')
                layout.addWidget(self.error_message)
                self.error_flag = 202

            #createDatabase(addedMolecule)
            self.writeToDatabase(addedMolecule, defaultDB)
            print('\n')
            print(addedMolecule)
            print('\n')
            fileCounter = self.countFiles(defaultDB)
            self.filesCount.setText(f"{str(fileCounter)} Attached Files")
            self.attachedFilesLabel.show()
            self.filesCount.show()
            self.attach_button.show()
            tempFiles = [f for f in listdir('./_core/UserAddedData/temp') if isfile(join('./_core/UserAddedData/temp', f))]
            for file in tempFiles:
                try:
                    remove(f'./_core/UserAddedData/temp/{file}')
                except:
                    pass
            return addedMolecule


        else:
            # Display an error message if the SMILES string is invalid
            self.error_message = QLabel('Invalid SMILES string. Please enter a valid SMILES.')
            #layout = self.layout()
            layout.addWidget(self.error_message)
            self.error_flag = 100

if __name__ == '__main__':
    with open('./_core/ThermalDex.log', 'w', encoding='utf-8') as logFile:
        with redirect_stdout(logFile):
            defaultDB, highEnergyGroups, expEnergyGroups, yoshidaMethod, qdscUnits, ambertd24limit, redtd24limit, oreohazardlimit, oreohazardwarningscale = altImportConfig()
            app = QApplication(sys.argv)
            app.setWindowIcon(QIcon('.\\_core\\ThermalDexIcon.ico'))
            window = MolDrawer()
            window.show()
            window.raise_()
            window.activateWindow()
            sys.exit(app.exec_())
