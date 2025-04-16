import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene, QFrame, QTableWidget, QTableWidgetItem, QTabWidget, QGraphicsPixmapItem,  QMessageBox, QComboBox #, QTableView, QToolTip
from PyQt5.QtGui import QPixmap, QColor, QIcon, QRegExpValidator #QValidator #, QCursor
from PyQt5.QtCore import Qt, QRegExp, pyqtSignal
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Mol, MolFromSmiles, MolFromSmarts, rdmolfiles
import pandas as pd
import re
import pyperclip

versionNumber = "0.6.4"

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

        InputContainLayout.addLayout(InputLeftLayout)
        #ResultsContainLayout.addWidget(QVLine())
        InputContainLayout.addLayout(InputRightLayout)
        #ResultsContainLayout.addWidget(QVLine())
        molecule_layout.addLayout(InputContainLayout)

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
        molecule_layout.addLayout(buttonContainLayout)

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
                      mol = MolFromSmiles(self.result_smiles)

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

        #def hover(url):
        #    if url:
        #        QToolTip.showText(QCursor.pos(), titles.get(url, url))
        #    else:
        #        QToolTip.hideText()

        about_tab = QWidget()
        about_layout = QVBoxLayout()
        about_title = QLabel("<b>About ThermalDex</b>\n\n")
        about_blank = QLabel("\nVersion: " + versionNumber + " (This is currently an alpha build)\n")
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
        tab_widget.addTab(about_tab, "About")

        layout.addWidget(tab_widget)



        # Set up the main window
        self.setLayout(layout)
        self.setGeometry(100, 100, 600, 825)
        self.setWindowTitle('ThermalDex')

    def mrvToSMILES(self):
        try:
            copyMolXML = pyperclip.paste()
            copyMolReal = rdmolfiles.MolFromMrvBlock(copyMolXML)
            smilesFromMarvin = rdmolfiles.MolToSmiles(copyMolReal)
            print(smilesFromMarvin)
            pyperclip.copy(smilesFromMarvin)

        except:
            print("No mrv XML found.")

    def showErrorMessage(self, errorCode):
        self.msg = QMessageBox()
        self.msg.setIcon(QMessageBox.Warning)
        self.msg.setText("An error has occured")
        self.msg.setInformativeText("The source of the problem seems to be with: " + errorCode + "\nTry again and contact developer if the problem persists.")
        self.msg.setWindowTitle("ThermalDex Error")
        self.msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        returnValue = self.msg.exec()

    def errorTestingTool(self):
        self.showErrorMessage("This is an Alternative test method for error handling")

    def clearTheCalcdValues(self):
        scene = QGraphicsScene()
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
        self.TE_input.setText('')
        self.Tinit_input.setText('')
        self.proj_input.setText('')

    def resetToDefaultState(self):
        self.clearTheCalcdValues()
        self.clearUserValues()
        layout = self.layout()
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)

    def createReport(self):
        moleculeData = self.render_molecule()
        img = moleculeData.molToIMG()
        results = asdict(moleculeData)
        create_pdf(results["name"], results, img) #results["molIMG"]) 
        if self.error_flag is not None:
            self.error_message.setText('')
            layout.removeWidget(self.error_message)

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
        Qdsc = self.Qdsc_input.text()
        QUnits = self.QUnitsSelection.currentText()
        TE = self.TE_input.text()
        Tinit = self.Tinit_input.text()
        proj = self.proj_input.text()

        # To be removed and replaced by call to thermalDexMolecule, see later.
        writeSmiles = '"'+smiles+'"' #repr(str(smiles))
        writeName = '"'+name+'"'
        writemp = '"'+mp+'"'
        writeTE = '"'+TE+'"'
        writeProj = '"'+proj+'"'

        # Create an RDKit molecule from the SMILES string
        addedMolecule = thermalDexMolecule(SMILES=smiles, name=name, mp=mp, Q_dsc=Qdsc, Qunits=QUnits, onsetT=TE, initT=Tinit, proj=proj)
        if smiles != '' and smiles is not None:
            addedMolecule.mol = MolFromSmiles(smiles)

        else:
            addedMolecule.mol = None
        
        if addedMolecule.mol is not None:
            # Make Pixmap Image to Display.
            pixmap = addedMolecule.molToQPixmap()
            scaledPixmap = pixmap #.scaled(550, 275, Qt.KeepAspectRatio)
            scene = QGraphicsScene()
            #scene.setSceneRect(0, 0, 400, 400)
            scene.addPixmap(scaledPixmap) #pixmap)
            self.mol_display.setScene(scene)
            addedMolecule.molPixmap = None

            # Calculate Core Properties
            addedMolecule.genCoreValuesFromMol()

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

            if addedMolecule.isStr != None:
                self.ISLabel.setText('Yoshida Impact Sensitivity: ' + addedMolecule.isStr + addedMolecule.IS_des)
            if addedMolecule.epStr != None:
                self.EPLabel.setText('Yoshida Explosive Propagation: ' + addedMolecule.epStr + addedMolecule.EP_des)
            if addedMolecule.Td24 != '' and addedMolecule.Td24 != 'nan' and addedMolecule.Td24 != None:
                d24Str = "{:.1f}".format(addedMolecule.Td24)
                self.Td24Label.setText('T<sub>D24</sub>: ' + '<b>' + d24Str + ' °' + 'C' + '</b>')
            if addedMolecule.onsetT != 'nan' and addedMolecule.onsetT != '' and addedMolecule.onsetT != None and addedMolecule.onsetT <= 24.99:
                self.error_message = QLabel('Sub-Ambient Onset Temperature? Are you sure? Yoshida values will not be accurate.')
                layout.addWidget(self.error_message)
                self.error_flag = 201
            if addedMolecule.initT != 'nan' and addedMolecule.initT != '' and addedMolecule.initT != None and addedMolecule.initT <= 24.99:
                self.error_message = QLabel('Sub-Ambient Initiation Temperature? Are you sure? Yoshida values will not be accurate.')
                layout.addWidget(self.error_message)
                self.error_flag = 202

            addMol = open(defaultDB, 'a')
            #addMol.write(writeSmiles + ',' + writeName + ',' + HEG + ',' + writemp + ',' + mwStr + ',' + writeTE + ',' + writeProj + '\n')
            print('\n')
            print(addedMolecule)
            print('\n')
            return addedMolecule


        else:
            # Display an error message if the SMILES string is invalid
            self.error_message = QLabel('Invalid SMILES string. Please enter a valid SMILES.')
            #layout = self.layout()
            layout.addWidget(self.error_message)
            self.error_flag = 100