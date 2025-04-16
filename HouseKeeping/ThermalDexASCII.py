import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene, QFrame, QTableWidget, QTableWidgetItem, QTabWidget, QToolTip, QGraphicsPixmapItem
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Mol
from PyQt5.QtGui import QPixmap, QColor, QIcon, QCursor
from PyQt5.QtCore import Qt
from io import BytesIO
import pandas as pd
import re

try:
    import pyi_splash
    pyi_splash.close()
except:
    pass

def importConfig():
    conf = open('ThermalDex.config', 'r')
    confCounter = 0
    for line in conf:
        #print(confCounter)
        if confCounter == 4:
           defaultDB = line.strip("\n")
           confCounter += 1
        elif confCounter == 8:
           highEnergyGroups = line.strip("\n")
           confCounter += 1
        else:
           confCounter += 1

    #print(defaultDB)
    return defaultDB, highEnergyGroups

class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)

class MolDrawer(QWidget):
    def __init__(self):
        super().__init__()

        self.init_ui()

    def init_ui(self):
        self.df = None
        #current_index = 0
        self.current_index = 0
        self.result_smiles = None
        self.error_flag = None
        # Set up the main layout
        layout = QVBoxLayout()

        # Create a tab widget
        tab_widget = QTabWidget()

        # Tab for molecule rendering
        molecule_tab = QWidget()
        molecule_layout = QVBoxLayout()

        # Display area for the molecular drawing
        self.mol_display = QGraphicsView(self)
        molecule_layout.addWidget(QLabel('Molecule:'))
        molecule_layout.addWidget(self.mol_display)

	#Add labels for calculated values
        self.mwLabel = QLabel('MW: ')
        self.mwLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        molecule_layout.addWidget(self.mwLabel)
        self.HEGlabel = QLabel('Number of High Energy Groups:')
        self.HEGlabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        molecule_layout.addWidget(self.HEGlabel)
        self.EFGlabel = QLabel('Number of Explosive Functional Groups:')
        self.EFGlabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        molecule_layout.addWidget(self.EFGlabel)
        self.eleLabel = QLabel('Elemental Composition: ')
        self.eleLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        molecule_layout.addWidget(self.eleLabel)
        self.RoSLabel = QLabel('Rule of Six: ')
        self.RoSLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        molecule_layout.addWidget(self.RoSLabel)
        self.obLabel = QLabel('Oxygen Balance: ')
        self.obLabel.setTextInteractionFlags(Qt.TextSelectableByMouse)
        molecule_layout.addWidget(self.obLabel)

        molecule_layout.addWidget(QHLine())

        self.tableLabel = QLabel('O.R.E.O.S. Assessment of Hazard by Scale:')
        molecule_layout.addWidget(self.tableLabel)
        self.table = QTableWidget(1, 4)
        self.table.setHorizontalHeaderLabels(['<5 g', '5 to <100 g', '100 to 500 g', '>500 g'])
        self.table.verticalHeader().setVisible(False)
        #self.table.setStyleSheet("QTableWidget::item { border-bottom: 2px solid black; }")
        self.table.setMaximumHeight(52)
        molecule_layout.addWidget(self.table)


        molecule_layout.addWidget(QHLine())

        # Input field for SMILES string
        self.smiles_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('Enter SMILES String:'))
        molecule_layout.addWidget(self.smiles_input)

        # Input field for Name string
        self.name_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('Name:'))
        molecule_layout.addWidget(self.name_input)

        # Input field for mp string
        self.mp_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('m.p.:'))
        molecule_layout.addWidget(self.mp_input)

        # Input field for Onset string
        self.TE_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('Onset Temperature:'))
        molecule_layout.addWidget(self.TE_input)

        # Input field for proj string
        self.proj_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('Project:'))
        molecule_layout.addWidget(self.proj_input)

        # Button to render the molecule
        render_button = QPushButton('Evaluate Molecule', self)
        render_button.clicked.connect(self.render_molecule)
        molecule_layout.addWidget(render_button)

        molecule_tab.setLayout(molecule_layout)
        tab_widget.addTab(molecule_tab, "Add")


        # Tab for Search
        search_tab = QWidget()
        search_layout = QVBoxLayout()

        # Entry widget for searching
        lbl_search = QLabel('Search:')
        entry_search = QLineEdit()
        result_label = QLabel('click search')
        counter_label = QLabel('none')

        def search_database():
        #    query = entry_search.text().lower()
        #    results = [entry for entry in self.database if query in entry['Name'].lower() or query in entry['Formula'].lower()]

        #    list_search_results.clear()
        #    for result in results:
        #        list_search_results.addItem(f"{result['Name']} ({result['Formula']})")

             self.df = pd.read_csv(defaultDB, encoding='mbcs')
             show_result(self)

        def show_result(self):
             #print(self)
             layout = self.layout()
             if self.error_flag is not None:
                  self.error_message.setText('')
                  layout.removeWidget(self.error_message)
                  self.error_flag = None
             if self.df is not None and not self.df.empty:
                  #print(self.current_index)
                  current_row = self.df.iloc[self.current_index]
                  #print(self.result_smiles)
                  self.result_smiles = current_row['Molecule']
                  #print(self.result_smiles)
                  result_text = f"SMILES: {current_row['Molecule']}\nName: {current_row['Name']}\nHEG: {current_row['Number of High Energy Groups']}\nmp: {current_row['Melting point']}\nMW: {current_row['Molecular Weight']}\nTE: {current_row['Thermal Event']}\nProject: {current_row['Project']}"
                  result_label.setText(result_text)
                  result_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
                  result_label.setWordWrap(True) 
                  counter_label.setText(f"Result {self.current_index + 1} of {len(self.df)}")
                  search_layout.addWidget(prev_button)
                  search_layout.addWidget(next_button)
                  #print(self.current_index)

                  if self.result_smiles is not None:
                      mol = Chem.MolFromSmiles(self.result_smiles)

                  else: 
                      mol = None

                  if mol is not None:
                      # Generate a molecular drawing as a PNG image
                      img = Draw.MolToImage(mol)

                      # Convert the image to a byte array
                      byte_array = BytesIO()
                      img.save(byte_array, format='PNG')

                      # Convert the byte array to a QPixmap and display it
                      pixmap = QPixmap()
                      pixmap.loadFromData(byte_array.getvalue())
                      scene = QGraphicsScene()
                      scene.addPixmap(pixmap)
                      self.mol_result_display.setScene(scene)

        def prev_result(self):
             if self.current_index > 0:
                  #search_layout.removeWidget(self.mol_result_display)
                  #search_layout.removeWidget(self.molLabel)
                  self.current_index -= 1
                  show_result(self)

        def next_result(self):
             if self.df is not None:
                  if self.current_index < len(self.df) - 1:
                       #search_layout.removeWidget(self.mol_result_display)
                       #search_layout.removeWidget(self.molLabel)
                       self.current_index += 1
                       show_result(self)


        # Search Buttons & display area for the molecular drawing
        self.mol_result_display = QGraphicsView(self)
        self.molLabel = QLabel('Molecule:')
        search_layout.addWidget(self.molLabel)
        search_layout.addWidget(self.mol_result_display)

        btn_search = QPushButton('Search', clicked=search_database)
        prev_button = QPushButton('Previous')
        next_button = QPushButton('Next')
        prev_button.clicked.connect(lambda: prev_result(self))
        next_button.clicked.connect(lambda: next_result(self))

        search_layout.addWidget(lbl_search)
        search_layout.addWidget(entry_search)
        search_layout.addWidget(result_label)
        search_layout.addWidget(counter_label)
        search_layout.addWidget(btn_search)

        search_tab.setLayout(search_layout)
        tab_widget.addTab(search_tab, "Search")


        # Tab for "About" information

        def hover(url):
            if url:
                QToolTip.showText(QCursor.pos(), titles.get(url, url))
            else:
                QToolTip.hideText()

        about_tab = QWidget()
        about_layout = QVBoxLayout()
        about_title = QLabel("<b>About ThermalDex</b>\n\n")
        about_blank = QLabel("\nVersion: 0.1.0  (This is currently an alpha build)\n")
        about_text = QLabel("\n\nThis is a simple tool for assessing and recording the potential thermal hazards assoicated with a molecule. It uses the <b>'O.R.E.O.S.'</b> assement scale and other ideas that can be read about in <a href=\"https://pubs.acs.org/doi/10.1021/acs.oprd.0c00467\"><em>Org. Process Res. Dev.</em> 2021, 25, 2, 212?224</a> by Jeffrey B. Sperry et. al.")
        iconLabel = QLabel()
        iconImage = QPixmap("ThermalDexIcon.jpg")
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
        #about_text.setAlignment(QtCore.Qt.AlignCenter)
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
        tab_widget.addTab(about_tab, "About")

        layout.addWidget(tab_widget)



        # Set up the main window
        self.setLayout(layout)
        self.setGeometry(100, 100, 448, 850)
        self.setWindowTitle('ThermalDex')

    def getColorForValue(self, hazardClass):
        # Example color-coding logic
        if hazardClass == 'High Hazard':
            return QColor(255, 0, 0)  # Red
        elif hazardClass == 'Medium Hazard':
            return QColor(255, 255, 0)  # Yellow
        elif hazardClass == 'Low Hazard':
            return QColor(0, 255, 0)  # Green
        else:
            return QColor(0, 0, 255)  # Blue

    def render_molecule(self):
        layout = self.layout()
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)
        # Get the SMILES string from the input field
        smiles = self.smiles_input.text()
        name = self.name_input.text()
        mp = self.mp_input.text()
        TE = self.TE_input.text()
        proj = self.proj_input.text()
        writeSmiles = '"'+smiles+'"' #repr(str(smiles))
        writeName = '"'+name+'"'
        writemp = '"'+mp+'"'
        writeTE = '"'+TE+'"'
        writeProj = '"'+proj+'"'

        # Create an RDKit molecule from the SMILES string
        if smiles != '' and smiles is not None:
            mol = Chem.MolFromSmiles(smiles)

        else:
            mol = None
        
        test = 'hi'

        if mol is not None:
            # Generate a molecular drawing as a PNG image
            img = Draw.MolToImage(mol)

            # Convert the image to a byte array
            byte_array = BytesIO()
            img.save(byte_array, format='PNG')

            # Convert the byte array to a QPixmap and display it
            pixmap = QPixmap()
            pixmap.loadFromData(byte_array.getvalue())
            scene = QGraphicsScene()
            scene.addPixmap(pixmap)
            self.mol_display.setScene(scene)
            cmpdMW = Descriptors.MolWt(mol)
            mwStr = "{:.2f}".format(cmpdMW)
            self.mwLabel.setText('MW: ' + mwStr)
            fullMatch = 0
            with open("HighEnergyGroups.csv", "r") as HEGroups:
               for line in HEGroups:
                    #print(line)
                    HeSubstructure = Chem.MolFromSmiles(line)
                    fullmatchList = Mol.GetSubstructMatches(mol, HeSubstructure)
                    if len(fullmatchList) > 0:
                        print('High Energy Group Found: ' + line)
                    fullMatch += len(fullmatchList)

            HEG = str(fullMatch)
            expMatch = 0
            with open("ExplosiveGroups.csv", "r") as expGroups:
               for line in expGroups:
                    #print(line)
                    expSubstructure = Chem.MolFromSmiles(line)
                    expmatchList = Mol.GetSubstructMatches(mol, expSubstructure)
                    #print('\n')
                    #print(expmatchList)
                    #print('\n')
                    if len(expmatchList) > 0:
                        print('Explosive Group Found: ' + line)
                    expMatch += len(expmatchList)

            EFG = str(expMatch)

            print('Explosive Groups =' + ' ' + EFG)

            self.HEGlabel.setText('Number of High Energy Groups: ' + HEG)
            self.EFGlabel.setText('Number of Explosive Functional Groups: ' + EFG)
            oreo = 0 # Initiallise O.R.E.O.S. calculation. Lack of 's' in var name is deliberate, as Scale is factored in later.        
            # O.R.E.O.S. check for EFG:
            if int(EFG) >= 1:
                oreo += 8
            elif int(EFG) == 0:
                oreo += 1
            else:
                print("Error: O.R.E.O.S. EFG number gives unexpected value: " + EFG)
            #carbonAtoms = len(Mol.GetSubstructMatches(mol, Chem.MolFromSmarts("[C]")))
            #print(carbonAtoms)
            chemForm = rdMolDescriptors.CalcMolFormula(mol)
            print(chemForm)
            eleComp = ""
            eleCompList = []
            eleList = []
            niceList = []
            #match = re.findall(r"[A-Z][a-z]?\d*|\((?:[^()]*(?:\(.*\))?[^()]*)+\)\d+", chemForm, re.I) #re.findall(r"([a-z]+)([0-9]+)", chemForm, re.I) #match
            #match = re.findall(r"[A-Z][a-z]?\d?", chemForm, re.I) #re.findall(r"([a-z]+)([0-9]+)", chemForm, re.I) #match
            #match = re.findall(r"[A-Za-z]\d?", chemForm, re.I)
            #match = re.findall(r"[A-Z][0-9]+?[a-z]?[0-9]+?", chemForm, re.I)
            #match = re.findall(r"[A-Z][a-z]\d+?|[A-Z]\d+?", chemForm, re.I)
            atomPattern = r"Cl\d*|H\d*|O\d*|N\d*|Si\d*|S\d*|F\d*|Cs\d*|Br\d*|I\d*|B\d*|Al\d*|Na\d*|K\d*|Mg\d*|Zn\d*|Ti\d*|Pd\d*|C\d*"
            match = re.findall(atomPattern, chemForm, re.I)
            if match:
               items = match #match.groups()
            for ele in match:
                print(ele)
                ment = re.findall(r"([a-z]+)([0-9]+)?", ele, re.I)
                #print(ment)
                matchedEleComp = ment[0]
                #print(matchedEleComp)
                eleList += [matchedEleComp]
            
            for compostion in eleList:
               #print(compostion)
               eleCompList += compostion[::-1]
               niceList += [('').join(compostion[::-1])]
   
            eleComp = (', ').join(niceList)
            #print(eleCompList)
            print(eleList)
            print(eleComp)
            self.eleLabel.setText('Elemental Composition: ' + eleComp)
            carbonAtoms = [ele[1] for ele in eleList if "C" in ele[0]]
            hydrogenAtoms = [ele[1] for ele in eleList if "H" in ele[0]]
            oxygenAtoms = [ele[1] for ele in eleList if "O" in ele[0]]
            print(oxygenAtoms)
            if carbonAtoms == []:
                Catoms = 0
            elif carbonAtoms[0] == '':
                Catoms = 1
            else:
                Catoms = int(carbonAtoms[0])
            if hydrogenAtoms == []:
                Hatoms = 0
            elif hydrogenAtoms[0] == '':
                Hatoms = 1
            else:
                Hatoms = int(hydrogenAtoms[0])
            if oxygenAtoms == []:
                Oatoms = 0
            elif oxygenAtoms[0] == '':
                Oatoms = 1
            else:
                Oatoms = int(oxygenAtoms[0])
            #print(int(carbonAtoms[0]))
            print(Catoms)
            ruleSixCrit = 6*int(HEG) - Catoms #int(carbonAtoms[0])
            print(ruleSixCrit)
            if ruleSixCrit > 0:
                ruleSix = "Explosive"
                oreo += 8
            else:
                ruleSix = "Not Explosive"
                oreo += 2
            self.RoSLabel.setText('Rule of Six: ' + ruleSix)
            #oxygenBalance = (-1600*((2*int(carbonAtoms[0]))+(int(hydrogenAtoms[0])/2)-int(oxygenAtoms[0])))/cmpdMW
            oxygenBalance = (-1600*((2*Catoms)+(Hatoms/2)-Oatoms))/cmpdMW
            print(oxygenBalance)
            obStr = "{:.2f}".format(oxygenBalance)
            if oxygenBalance > 160:
                obRisk = "(Low Risk)"
                oreo += 2
            elif oxygenBalance > 80 and oxygenBalance <= 160:
                obRisk = "(Medium Risk)"
                oreo += 4
            elif oxygenBalance >= -120 and oxygenBalance <= 80:
                obRisk = "(High Risk)"
                oreo += 8
            elif oxygenBalance >= -240 and oxygenBalance < -120:
                obRisk = "(Medium Risk)"
                oreo += 4
            elif oxygenBalance < -240:
                obRisk = "(Low Risk)"
                oreo += 2
            self.obLabel.setText('Oxygen Balance: ' + obStr + ' ' + obRisk)
            
            # Calculation of O.R.E.O.S. Safe Scale
            # Onset temp adjustment
            if TE != 'nan' and TE != '':
                onsetTemp = int(TE)
                print('Onset Temperature given as' + TE)
                if onsetTemp <= 125:
                    oreo += 8
                elif onsetTemp in range(126,200):
                    oreo += 4
                elif onsetTemp in range(200,300):
                    oreo += 2
                elif onsetTemp >= 300:
                    oreo += 1
            else:
                onsetTemp = None

            print(oreo)
            print(onsetTemp)

            largeScaleSafety = oreo + 8
            hundsScaleSafety = oreo + 4
            tensScaleSafety = oreo + 2
            smallScaleSafety = oreo + 1

            scaleList = [smallScaleSafety, tensScaleSafety, hundsScaleSafety, largeScaleSafety]
            hazardList = []
            hazardValuesRangeList = []

            if onsetTemp == None:
                hazardValuesRangeList = [15, 22]
            if onsetTemp != None:
                hazardValuesRangeList = [18, 28]

            for scale in scaleList:
                if scale < hazardValuesRangeList[0]:
                    oreosHazard = "Low Hazard"
                    hazardList.append(oreosHazard)
                elif scale in range(hazardValuesRangeList[0], hazardValuesRangeList[1]):  # Note: Remember the last integer in range isnt included thats why I can get away with a list of 2 items. 
                    oreosHazard = "Medium Hazard"
                    hazardList.append(oreosHazard)
                elif scale >= hazardValuesRangeList[1]:
                    oreosHazard = "High Hazard"
                    hazardList.append(oreosHazard)
                else:
                    print('Error: ' + str(scale))
            print(hazardList)
            self.table.clearContents()
            # Add values to cells


            for i, hazardClass in enumerate(hazardList):
                item = QTableWidgetItem(hazardClass)
                self.table.setItem(0, i, item)

                # Color code cells based on values
                classColor = self.getColorForValue(hazardClass)
                print(classColor)
                item.setBackground(classColor)

            addMol = open(defaultDB, 'a')
            addMol.write(writeSmiles + ',' + writeName + ',' + HEG + ',' + writemp + ',' + mwStr + ',' + writeTE + ',' + writeProj + '\n')


        else:
            # Display an error message if the SMILES string is invalid
            self.error_message = QLabel('Invalid SMILES string. Please enter a valid SMILES.')
            #layout = self.layout()
            layout.addWidget(self.error_message)
            self.error_flag = 100

if __name__ == '__main__':
    defaultDB, highEnergyGroups = importConfig()
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('ThermalDexIcon.ico'))
    window = MolDrawer()
    window.show()
    sys.exit(app.exec_())
