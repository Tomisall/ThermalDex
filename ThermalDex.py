from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QShortcut, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene, QFrame, QTableWidget, QTableWidgetItem, QTabWidget, QGraphicsPixmapItem,  QMessageBox, QComboBox, QSpacerItem, QCheckBox, QSizePolicy
from PyQt5.QtGui import QPixmap, QColor, QIcon, QRegExpValidator, QFont, QPalette, QKeySequence
from PyQt5.QtCore import Qt, QRegExp, pyqtSignal, QEvent
from thermDex.thermDexMolecule import *
from thermDex.thermDexHTMLRep import *
from thermDex.attachedFileManager import *
from thermDex.Section import Section
from thermDex.thermDexPlots import *
from thermDex.thermDexPandasTools import *
from thermDex.thermDexTextEdit import *
from thermDex.thermDexTable import *
import pyperclip
import configparser
import sqlite3
from numpy import log10
from contextlib import redirect_stdout
from os import path, environ
from shutil import copy2

versionNumber = "1.4.1a"

try:
    import pyi_splash
    pyi_splash.close()
    window.raise_()
except:
    pass

os.environ["QT_ENABLE_HIGHDPI_SCALING"]   = "1"
os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"
os.environ["QT_SCALE_FACTOR"]             = "1"

if hasattr(Qt, 'AA_DisableHighDpiScaling'): #'AA_EnableHighDpiScaling'
    QApplication.setAttribute(Qt.AA_DisableHighDpiScaling, True)

if hasattr(Qt, 'AA_UseHighDpiPixmaps'):
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

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

class RichTextPushButton(QPushButton):
    def __init__(self, parent=None, text=None):
        if parent is not None:
            super().__init__(parent)
        else:
            super().__init__()
        self.__lbl = QLabel(self)
        if text is not None:
            self.__lbl.setText(text)
        self.__lyt = QHBoxLayout()
        self.__lyt.setContentsMargins(0, 0, 0, 0)
        self.__lyt.setSpacing(0)
        self.setLayout(self.__lyt)
        self.__lbl.setAttribute(Qt.WA_TranslucentBackground)
        self.__lbl.setAttribute(Qt.WA_TransparentForMouseEvents)
        self.__lbl.setSizePolicy(
            QSizePolicy.Expanding,
            QSizePolicy.Expanding,
        )
        self.__lbl.setTextFormat(Qt.RichText)
        self.__lyt.addWidget(self.__lbl)
        return

    def setText(self, text):
        self.__lbl.setText(text)
        self.updateGeometry()
        return

    def sizeHint(self):
        s = QPushButton.sizeHint(self)
        w = self.__lbl.sizeHint()
        s.setWidth(w.width())
        s.setHeight(w.height())
        return s

class Td24OverrideWindow(QWidget):
    submitClicked = pyqtSignal(bool,float)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Override Value")
        self.setGeometry(100, 100, 250, 100)
        layout = QVBoxLayout()
        self.setLayout(layout)

        self.measured = False
        self.td24_measured_float = None

        self.override_value = QLineEdit(self)
        layout.addWidget(QLabel('Enter Measured T<sub>D24</sub> value:'))
        layout.addWidget(self.override_value)
        self.overide_button = RichTextPushButton()
        self.overide_button.setText('&nbsp;&nbsp;Record T<sub>D24</sub>')
        self.overide_button.setStyleSheet("text-align:center;")
        self.overide_button.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        button_centre = QHBoxLayout()
        button_centre.addWidget(self.overide_button)
        self.overide_button.clicked.connect(self.override_button_click)
        layout.addLayout(button_centre)

    def interactiveErrorMessage(self, errorInfo):
        self.interactMsg = QMessageBox()
        self.interactMsg.setIcon(QMessageBox.Information)
        self.interactMsg.setText("Action Needed")
        self.interactMsg.setInformativeText(errorInfo)
        self.interactMsg.setWindowTitle("ThermalDex - Info Box")
        self.interactMsg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        returnValue = self.interactMsg.exec()
        return returnValue

    def open_window(self):
        self.newWindow = Td24OverrideWindow(self)
        self.newWindow.show()

    def override_button_click(self):
        try:
            print(f'{self.override_value.text()=}')
            float_of_empirical_td24 = float(self.override_value.text())
            print(f'{float_of_empirical_td24=}')
            self.measured = True
            self.td24_measured_float = float_of_empirical_td24
            self.submitClicked.emit(self.measured,self.td24_measured_float)
        except Exception as e:
            errorInfo = f"Enter Valid Number Please:\n\nDetailed info: {e}"
            self.interactiveErrorMessage(errorInfo)

        self.close()

class CommentsBox(QWidget):
    submitClicked = pyqtSignal()

    def __init__(self, comments_location):
        super().__init__()
        self.setWindowTitle("Comments Box")
        self.setGeometry(100, 100, 250, 100)
        self.comments_location = comments_location
        layout = QVBoxLayout()
        self.setLayout(layout)
        self.text_edit = TextEdit()
        if os.path.isfile(f'{comments_location}/comments.html'):
            with open(f'{comments_location}/comments.html') as file:
                self.text_edit.setHtml(file.read())

        

        # Shortcuts
        save_shortcut = QKeySequence(Qt.CTRL + Qt.Key_S)



        bold_shortcut = QKeySequence(Qt.CTRL + Qt.Key_B)
        italic_shortcut = QKeySequence(Qt.CTRL + Qt.Key_I)
        underline_shortcut = QKeySequence(Qt.CTRL + Qt.Key_U)
        table_shortcut = QKeySequence(Qt.CTRL + Qt.Key_T)

        # Tool bar
        tool_bar = QToolBar('Writing Tools')

        tb_save = QAction(icon=QIcon('./_core/app_icons/floppy-disk.png'),text='Save',parent=self)
        tb_save.setShortcut(save_shortcut)
        tb_save.triggered.connect(self.save_body_html)

        tb_bold = QAction(icon=QIcon('./_core/app_icons/bold.png'),text='Bold',parent=self)
        tb_bold.setShortcut(bold_shortcut)
        tb_bold.triggered.connect(self.set_text_bold)

        tb_italic = QAction(icon=QIcon('./_core/app_icons/italic.png'),text='Italic',parent=self)
        tb_italic.setShortcut(italic_shortcut)
        tb_italic.triggered.connect(self.set_text_italic)

        tb_underline = QAction(icon=QIcon('./_core/app_icons/underline.png'),text='Underline',parent=self)
        tb_underline.setShortcut(underline_shortcut)
        tb_underline.triggered.connect(self.set_text_underline)

        tb_colour_picker = QAction(icon=QIcon('./_core/app_icons/color-wheel.png'),text='Color Picker',parent=self)
        #tb_colour_picker.setShortcut(xxxxx_shortcut)
        tb_colour_picker.triggered.connect(self.set_text_colour)

        tb_font_picker = QAction(icon=QIcon('./_core/app_icons/text.png'),text='Font Picker',parent=self)
        #tb_font_picker.setShortcut(xxxxx_shortcut)
        tb_font_picker.triggered.connect(self.set_text_font)

        tb_table = QAction(icon=QIcon('./_core/app_icons/table.png'),text='Insert Table',parent=self)
        tb_table.setStatusTip("Insert table")
        tb_table.setShortcut(table_shortcut)
        tb_table.triggered.connect(Table(parent=self.text_edit).show)

        tool_bar.addAction(tb_save)
        tool_bar.addAction(tb_bold)
        tool_bar.addAction(tb_italic)
        tool_bar.addAction(tb_underline)
        tool_bar.addAction(tb_colour_picker)
        tool_bar.addAction(tb_font_picker)
        tool_bar.addAction(tb_table)
        layout.addWidget(tool_bar)
        layout.addWidget(self.text_edit)

    def save_body_html(self):
        output_html = self.text_edit.toHtml()
        #data_img_html = local_parse_HTML_for_imgs(output_html)
        #body_cont_html = local_get_body_contents(data_img_html)
        #Html_no_newlines = body_cont_html.replace('\n', '').replace('\r', '')
        with open(f'{self.comments_location}/comments.html', "w") as text_file:
            text_file.write(output_html)
        
        self.submitClicked.emit()
    
    def set_text_bold(self):
        current_font_weight = self.text_edit.fontWeight()
        print(current_font_weight)
        if current_font_weight < 75:
            self.text_edit.setFontWeight(75)#QFont.Bold)
        else:
            self.text_edit.setFontWeight(50)#QFont.Normal)

    def set_text_italic(self):
        current_font_italic = self.text_edit.fontItalic()
        print(current_font_italic)
        if current_font_italic == True:
            self.text_edit.setFontItalic(False)
        else:
            self.text_edit.setFontItalic(True)

    def set_text_underline(self):
        current_font_underline = self.text_edit.fontUnderline()
        print(current_font_underline)
        if current_font_underline == True:
            self.text_edit.setFontUnderline(False)
        else:
            self.text_edit.setFontUnderline(True)

    def set_text_colour(self):
        dialog = ColorPickerDialog(color=self.text_edit.textColor(), orientation='horizontal')#QColor(0, 0, 0), orientation='horizontal')
        reply = dialog.exec()
        if reply == QDialog.Accepted: 
            color = dialog.getColor() # return type is QColor
            self.text_edit.setTextColor(color)

    def set_text_font(self):
        font, valid = QFontDialog.getFont()
        if valid:
            self.text_edit.setFont(font)

class MolDrawer(QWidget):
    def __init__(self, parent=None):
        super(MolDrawer, self).__init__(parent)
        app.installEventFilter(self)
        self.activateWindow()
        self.raise_()
        self.init_ui()

    def eventFilter(self, obj, event):
        if event.type() == QEvent.WindowStateChange:
            if event.oldState() == Qt.WindowNoState or self.windowState() == Qt.WindowMaximized:
                print("WindowMaximized")
                self.give_max_layout()
            else:
                self.give_normal_layout()
                #self.setGeometry(100, 100, 633, 950)
                self.resize(633, 950)
        return super(MolDrawer, self).eventFilter(obj, event)


    def init_ui(self):
        self.df = None
        #current_index = 0
        self.current_index = 0
        self.result_smiles = None
        self.selectedDatabase = None
        self.error_flag = None
        self.DSC_run_flag = False
        # Set up the main layout

        layout = QVBoxLayout()
        
        # Set up the main window
        self.setLayout(layout)
        self.setGeometry(100, 100, 600, 895)
        self.setWindowTitle('ThermalDex')
        #self.setMinimumSize(450,810)
        self.resize(int(screen.size().width()*0.33), int(screen.size().height()*0.85))
        #self.setFixedSize(int(screen.size().width()*0.33), int(screen.size().height()*0.85))

        # Create a tab widget
        self.tab_widget = QTabWidget()

        # Tab for molecule rendering
        self.molecule_tab = QWidget()
        self.mol_overall_layout = QHBoxLayout()
        self.mol_left_layout = QVBoxLayout()
        self.molecule_layout = QVBoxLayout()

        # Display area for the molecular drawing
        self.mol_display = QGraphicsView(self)
        self.top_info_sublayout = QHBoxLayout()
        self.mol_view_title_label = QLabel('Molecule:')
        self.top_info_sublayout.addWidget(self.mol_view_title_label)
        self.top_info_sublayout.addStretch()
        self.approval_needed = QLabel("<h3 style='color: red;'>Seek Approval Before Use</h3>")
        self.top_info_sublayout.addWidget(self.approval_needed)
        self.approval_needed.hide()



        self.ResultsContainLayout = QHBoxLayout()
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

        self.overrideOn = False
        self.overideValue = None
        self.measuredTd24 = QShortcut(QKeySequence('Ctrl+O'), self)
        self.measuredTd24.activated.connect(self.openTd24Override)
        self.openCommentsBox = QShortcut(QKeySequence('Ctrl+K'), self)
        self.openCommentsBox.activated.connect(self.openComments)

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

        self.ResultsContainLayout.addWidget(QVLine())
        self.ResultsContainLayout.addLayout(ResultsLeftLayout)
        self.ResultsContainLayout.addWidget(QVLine())
        self.ResultsContainLayout.addLayout(ResultsRightLayout)
        self.ResultsContainLayout.addWidget(QVLine())


        #self.molecule_layout.addLayout(self.ResultsContainLayout)
        #self.molecule_layout.addWidget(QHLine())


        self.tableLabel = QLabel('O.R.E.O.S. Assessment of Hazard by Scale:')
        #self.molecule_layout.addWidget(self.tableLabel)
        self.tableLayout = QHBoxLayout()
        self.table = QTableWidget(1, 4)
        self.table.setHorizontalHeaderLabels(['<5 g', '5 to <100 g', '100 to 500 g', '>500 g'])
        self.table.verticalHeader().setVisible(False)
        #self.table.setStyleSheet("QTableWidget::item { border-bottom: 2px solid black; }")
        #self.table.setMaximumHeight(53)
        #self.table.setMaximumWidth(402)
        #self.table.setMinimumHeight(53)
        #self.table.setMinimumWidth(402)
        if dpi >= 150:
            self.table.setMaximumHeight(int(screen.size().height()*0.069)) #53) self.resize(int(screen.size().width()*0.33), int(screen.size().height()*0.85))
            self.table.setMaximumWidth(int(screen.size().width()*0.314))
            self.table.setMinimumHeight(int(screen.size().height()*0.069)) #53) self.resize(int(screen.size().width()*0.33), int(screen.size().height()*0.85))
            self.table.setMinimumWidth(int(screen.size().width()*0.314)) 

        else:
            self.table.setMaximumHeight(int(screen.size().height()*0.0495)) #53) self.resize(int(screen.size().width()*0.33), int(screen.size().height()*0.85))
            self.table.setMaximumWidth(int(screen.size().width()*0.2097))
            self.table.setMinimumHeight(int(screen.size().height()*0.0495)) #53) self.resize(int(screen.size().width()*0.33), int(screen.size().height()*0.85))
            self.table.setMinimumWidth(int(screen.size().width()*0.2097))       
       #self.table.setAlignment(Qt.AlignVCenter)
        self.tableLayout.addWidget(self.table)

        #self.molecule_layout.addLayout(tableLayout)
        #self.molecule_layout.addWidget(QHLine())



        # Input field for SMILES string
        #self.smiles_input = QLineEdit(self)
        self.smiles_input = ClickableLineEdit(self)
        self.smiles_input.clicked.connect(self.mrvToSMILES)
        self.smiles_input_label = QLabel('Enter SMILES String:')
        #self.molecule_layout.addWidget(self.smiles_input_label)
        #self.molecule_layout.addWidget(self.smiles_input)

        self.InputContainLayout = QHBoxLayout()
        InputLeftLayout = QVBoxLayout()
        InputRightLayout = QVBoxLayout()
        numValidator = QRegExpValidator(QRegExp(r'[-]?\d+[.]?\d*'))
        posNumValidator = QRegExpValidator(QRegExp(r'\d+[.]?\d*'))

        # Input field for Name string
        self.name_input = QLineEdit(self)
        InputLeftLayout.addWidget(QLabel('Name:'))
        nameUnitsSubLayout = QHBoxLayout()
        nameUnitsSubLayout.addWidget(self.name_input)
        nameUnitsSubLayout.addSpacing(37)
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
        mpUnitsSubLayout.addStretch()
        InputLeftLayout.addLayout(mpUnitsSubLayout)

        # Input field for Q_DSC string
        self.Qdsc_input = QLineEdit(self)
        self.Qdsc_input.setValidator(posNumValidator)
        QLabelSubLayout = QHBoxLayout()
        QLabelSubLayout.addWidget(QLabel('Q<sub>DSC</sub>:'))
        # DSC run check
        self.dsc_check = QCheckBox('No Peak')
        self.dsc_check.setChecked(False)
        self.dsc_check.stateChanged.connect(self.check_if_dsc_performed)
        #QLabelSubLayout.addSpacing(100)
        QLabelSubLayout.addStretch()
        QLabelSubLayout.addWidget(self.dsc_check)
        #QLabelSubLayout.addStretch()
        QUnitsSubLayout = QHBoxLayout()
        QUnitsSubLayout.addWidget(self.Qdsc_input)
        #QUnitsSubLayout.addWidget(QLabel('cal g<sup>-1</sup> '))
        self.QUnitsSelection = QComboBox(self)
        self.QUnitsSelection.addItems(['J g⁻¹', 'cal g⁻¹'])
        self.QUnitsSelection.setCurrentIndex(int(qdscUnits))
        QUnitsSubLayout.addWidget(self.QUnitsSelection)
        InputRightLayout.addLayout(QLabelSubLayout)
        InputRightLayout.addLayout(QUnitsSubLayout)

        # Input field for Onset string
        self.TE_input = QLineEdit(self)
        self.TE_input.setValidator(numValidator)
        InputRightLayout.addWidget(QLabel('Onset Temperature:'))
        TEUnitsSubLayout = QHBoxLayout()
        TEUnitsSubLayout.addWidget(self.TE_input)
        TEUnitsSubLayout.addWidget(QLabel('°C'))
        TEUnitsSubLayout.addSpacing(45)
        InputRightLayout.addLayout(TEUnitsSubLayout)

        # Input field for Init string
        self.Tinit_input = QLineEdit(self)
        self.Tinit_input.setValidator(numValidator)
        InputRightLayout.addWidget(QLabel('Initiation Temperature:'))
        TinitUnitsSubLayout = QHBoxLayout()
        TinitUnitsSubLayout.addWidget(self.Tinit_input)
        TinitUnitsSubLayout.addWidget(QLabel('°C'))
        TinitUnitsSubLayout.addSpacing(45)
        InputRightLayout.addLayout(TinitUnitsSubLayout)

        # Input field for proj string
        self.proj_input = QLineEdit(self)
        InputLeftLayout.addWidget(QLabel('Project:'))
        projUnitsSubLayout = QHBoxLayout()
        projUnitsSubLayout.addWidget(self.proj_input)
        projUnitsSubLayout.addSpacing(37)
        InputLeftLayout.addLayout(projUnitsSubLayout)

        # Input field for Hammer Drop Test
        hamSubLayout = QHBoxLayout()
        InputLeftLayout.addWidget(QLabel('Hammer Drop Test:'))
        self.hamSelection = QComboBox(self)
        self.hamSelection.addItems(['Not Performed', 'Positive', 'Negative'])
        hamSubLayout.addWidget(self.hamSelection)
        hamSubLayout.addSpacing(37)
        InputLeftLayout.addLayout(hamSubLayout)

        # Input field for Friction Test
        fricSubLayout = QHBoxLayout()
        InputRightLayout.addWidget(QLabel('Friction Test:'))
        self.fricSelection = QComboBox(self)
        self.fricSelection.addItems(['Not Performed', 'Positive', 'Negative'])
        fricSubLayout.addWidget(self.fricSelection)
        fricSubLayout.addSpacing(62)
        InputRightLayout.addLayout(fricSubLayout)

        self.InputContainLayout.addLayout(InputLeftLayout)
        #ResultsContainLayout.addWidget(QVLine())
        self.InputContainLayout.addLayout(InputRightLayout)
        #ResultsContainLayout.addWidget(QVLine())
        #self.molecule_layout.addLayout(self.InputContainLayout)


        # Attach Data Files
        self.filesSubLayout = QHBoxLayout()
        self.attachedFilesLabel = QLabel('Attached Files:')
        self.attachedFilesLabel.resize(120, 120)
        self.filesSubLayout.addWidget(self.attachedFilesLabel)
        self.filesCount = QLabel('0 Attached Files')
        self.filesCount.resize(90, 120)
        self.filesSubLayout.addWidget(self.filesCount)
        self.attach_button = QPushButton('Add/View Files', self)
        self.attach_button.clicked.connect(self.openFileManager) 
        self.attach_button.setMaximumWidth(140)
        self.filesSubLayout.addWidget(self.attach_button)
        attachSpacerLabel = QLabel('hidden spacer')
        attachSpacerLabel.setStyleSheet('color: white')
        self.filesSubLayout.addWidget(attachSpacerLabel)
        attachSpacerLabelTwo = QLabel('hidden spacer')
        attachSpacerLabelTwo.setStyleSheet('color: white')
        self.filesSubLayout.addWidget(attachSpacerLabelTwo)
        #self.molecule_layout.addSpacing(15)
        self.file_spacer = QSpacerItem(15, 15)
        #self.molecule_layout.addLayout(self.filesSubLayout)


        # Buttons
        buttonSpacerLabel = QLabel('hidden spacer')
        buttonSpacerLabel.setStyleSheet('color: white')
        self.buttonContainLayout = QHBoxLayout()
        render_button = QPushButton('Evaluate Molecule', self)
        render_button.clicked.connect(self.render_molecule)
        render_button.setMaximumWidth(180)
        self.buttonContainLayout.addWidget(render_button)
        clear_button = QPushButton('Clear Sheet', self)
        clear_button.clicked.connect(self.resetToDefaultState)
        clear_button.setMaximumWidth(180)
        self.buttonContainLayout.addWidget(clear_button)
        #molecule_layout.addLayout(buttonContainLayout)
        msg_button = QPushButton('Save to PDF', self)
        msg_button.clicked.connect(self.createReport)
        msg_button.setMaximumWidth(180)
        self.buttonContainLayout.addWidget(msg_button)
        ''' Button for error testing/Debugging.
        error_testing_button = QPushButton('SecretErrorButton')
        error_testing_button.clicked.connect(self.errorTestingTool)
        buttonContainLayout.addWidget(error_testing_button)       
        '''
        #self.molecule_layout.addStretch(5)
        #self.molecule_layout.addLayout(self.buttonContainLayout)

        self.layout_spacer = QSpacerItem(20, 80)
        self.hline = QHLine()
        self.hlineTwo = QHLine()
        self.hlineThree = QHLine()
        self.strechBox = QVBoxLayout()
        self.strechBox.addStretch()
        self.strechBoxLabel = QLabel('')
        self.strechBox.addWidget(self.strechBoxLabel)
        self.lhs_title = QLabel('Title')
        self.lhs_title.setStyleSheet('color: white')

        self.molecule_layout.addLayout(self.ResultsContainLayout)
        self.molecule_layout.addWidget(self.hline)
        self.molecule_layout.addWidget(self.tableLabel)
        self.molecule_layout.addLayout(self.tableLayout)
        self.molecule_layout.addWidget(self.hlineThree)
        self.molecule_layout.addWidget(self.smiles_input_label)
        self.molecule_layout.addWidget(self.smiles_input)
        self.molecule_layout.addLayout(self.InputContainLayout)
        self.molecule_layout.addSpacerItem(self.file_spacer)
        self.molecule_layout.addLayout(self.filesSubLayout)
        self.molecule_layout.addWidget(buttonSpacerLabel)
        self.molecule_layout.addLayout(self.buttonContainLayout)

        #self.mol_left_layout_check()
        self.mol_overall_layout.addLayout(self.mol_left_layout)
        self.mol_overall_layout.addLayout(self.molecule_layout)
        self.molecule_tab.setLayout(self.mol_overall_layout) #self.mol_left_layout_check())
        self.tab_widget.addTab(self.molecule_tab, "Add")

        # Hide File Handling Widgets 
        self.attachedFilesLabel.hide()
        self.filesCount.hide()
        self.attach_button.hide()

        # Tab for Search
        search_tab = QWidget()
        self.search_layout = QVBoxLayout()
        self.left_search_layout = QVBoxLayout()
        self.overall_search_layout = QHBoxLayout()
        self.search_widget_index = 0

        # Entry widget for searching
        lbl_search = QLabel('Search:')
        self.entry_search = ClickableLineEdit(self) #QLineEdit()
        self.entry_search.clicked.connect(self.mrvToSMILES)
        self.entry_search.returnPressed.connect(lambda: self.search_database(self.entry_search.text())) 
        self.searchTypeSelection = QComboBox(self)
        self.searchTypeSelection.addItems(['SMILES', 'Name', 'Project', 'MW', 'Qdsc', 'Tinit', 'Tonset', 'Td24', 'O.R.E.O.S. at >500 g', 'Rule of Six', 'Oxygen Balance'])
        smilesSearchList = ['Substructure', 'Exact']
        nameProjSearchList = ['is', 'contains', 'is not', 'does not contain']
        valuesSearchList = ['=', '<', '>', '</=', '>/=']
        self.listOfSearchTypes = [smilesSearchList, nameProjSearchList, nameProjSearchList, valuesSearchList, valuesSearchList, valuesSearchList, valuesSearchList, valuesSearchList, valuesSearchList, valuesSearchList, valuesSearchList]
        self.searchSubType = QComboBox(self)
        self.searchSubType.addItems(smilesSearchList)
        self.searchTypeSelection.currentIndexChanged.connect(self.set_sub_search)
        self.result_label = QLabel('click search')
        self.counter_label = QLabel('none')

        # Search Buttons & display area for the molecular drawing
        self.mol_result_display = QGraphicsView(self)
        #self.mol_result_display.scale(0.7,0.7)
        self.molLabel = QLabel('Molecule:')
        #self.search_layout.addWidget(self.molLabel)
        #self.search_layout.addWidget(self.mol_result_display)
        self.mol_result_display.hide()
        self.molLabel.hide()
        self.give_normal_layout()        

        btn_search = QPushButton('Search')
        self.edit_button = QPushButton('Edit')
        self.del_button = QPushButton('Delete')
        self.export_button = QPushButton('Export Results')
        self.gen_multi_rep_button = QPushButton('Generate Report')
        self.edit_button.hide()
        self.del_button.hide()
        self.export_button.hide()
        self.gen_multi_rep_button.hide()
        eddel_sublayout = QHBoxLayout()
        eddel_sublayout.addWidget(self.export_button)
        eddel_sublayout.addWidget(self.gen_multi_rep_button)
        eddel_sublayout.addWidget(self.del_button)
        eddel_sublayout.addWidget(self.edit_button)
        self.prev_button = QPushButton('Previous')
        self.next_button = QPushButton('Next')
        self.first_button = QPushButton('First')
        self.last_button = QPushButton('Last')
        btn_search.clicked.connect(lambda: self.search_database(self.entry_search.text()))
        self.prev_button.clicked.connect(self.prev_result)
        self.next_button.clicked.connect(lambda: self.next_result(self.selectedDatabase))
        self.first_button.clicked.connect(self.first_result)
        self.last_button.clicked.connect(lambda: self.last_result(self.selectedDatabase))
        self.edit_button.clicked.connect(self.changeTabForEditing)
        self.del_button.clicked.connect(lambda: self.subFunctionForDeletion(self.plot_current_value,defaultDB))
        self.export_button.clicked.connect(lambda: self.exportSearchResults(self.selectedDatabase))
        self.gen_multi_rep_button.clicked.connect(lambda: self.createMultiReport(self.selectedDatabase))

        self.results_table_Label = QLabel('O.R.E.O.S. Assessment of Hazard by Scale:')
        self.results_table = QTableWidget(1, 4)
        self.results_table.setHorizontalHeaderLabels(['<5 g', '5 to <100 g', '100 to 500 g', '>500 g'])
        self.results_table.verticalHeader().setVisible(False)

        if dpi >= 150:
            self.results_table.setMaximumHeight(int(screen.size().height()*0.069)) #53) self.resize(int(screen.size().width()*0.33), int(screen.size().height()*0.85))
            self.results_table.setMaximumWidth(int(screen.size().width()*0.314))
            self.results_table.setMinimumHeight(int(screen.size().height()*0.069)) #53) self.resize(int(screen.size().width()*0.33), int(screen.size().height()*0.85))
            self.results_table.setMinimumWidth(int(screen.size().width()*0.314)) 

        else:
            self.results_table.setMaximumHeight(int(screen.size().height()*0.0495)) #53) self.resize(int(screen.size().width()*0.33), int(screen.size().height()*0.85))
            self.results_table.setMaximumWidth(int(screen.size().width()*0.2097))
            self.results_table.setMinimumHeight(int(screen.size().height()*0.0495)) #53) self.resize(int(screen.size().width()*0.33), int(screen.size().height()*0.85))
            self.results_table.setMinimumWidth(int(screen.size().width()*0.2097))         
        #self.results_table.setMinimumHeight(53)
        #self.results_table.setMinimumWidth(402)
        self.results_approval_warning = QLabel('')

        sub_search_layout = QHBoxLayout()
        self.search_layout.addWidget(lbl_search)
        sub_search_layout.addWidget(self.searchTypeSelection)
        sub_search_layout.addWidget(self.searchSubType)
        sub_search_layout.addWidget(self.entry_search)
        sub_search_layout.addWidget(btn_search)
        self.search_layout.addLayout(sub_search_layout)
        self.search_layout.addLayout(eddel_sublayout)
        self.search_layout.addStretch()

        self.make_plot_label = QLabel('For search results: ')
        self.make_plot_button = QPushButton('Plot')
        self.select_x_values = QComboBox(self)
        self.vs_label = QLabel('<em>vs</em>')
        self.select_y_values = QComboBox(self)
        selectable_xy_values = ['Qdsc', 'Tonset', 'Tinit', 'Td24', 'IS', 'EP', 'OB', 'RoS', 'OREOS', 'MW', 'Log(Qdsc)', 'Log(Tonset-25)', 'Log(Tinit-25)']
        self.select_x_values.addItems(selectable_xy_values)
        self.select_y_values.addItems(selectable_xy_values)
        make_plot_layout = QHBoxLayout()
        make_plot_layout.addWidget(self.make_plot_label)
        make_plot_layout.addWidget(self.make_plot_button)
        make_plot_layout.addWidget(self.select_x_values)
        make_plot_layout.addWidget(self.vs_label)
        make_plot_layout.addWidget(self.select_y_values)
        make_plot_layout.addStretch()
        self.search_layout.addLayout(make_plot_layout)
        self.make_plot_label.hide()
        self.make_plot_button.hide()
        self.select_x_values.hide()
        self.vs_label.hide()
        self.select_y_values.hide()
        self.make_plot_button.clicked.connect(lambda: self.plotSearchResults(self.plot_database, self.plot_current_value))

        self.overall_search_layout.addLayout(self.left_search_layout)
        self.overall_search_layout.addLayout(self.search_layout)
        search_tab.setLayout(self.overall_search_layout)
        self.tab_widget.addTab(search_tab, "Search")

        # Tab for Settings
        settings_tab = QWidget()
        settings_layout = QVBoxLayout()
        self.config = configparser.ConfigParser()
        self.config.read('./_core/ThermalDex.ini')
        settings_intro = QLabel('<h1>ThermalDex Settings</h1><p>This pane contains the settings for this application. Change relavent settings as needed.</p>')
        settings_layout.addWidget(settings_intro)

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
        self.oreoWarnHazCombo.setMaximumWidth(110)
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
        import_db_button = QPushButton('Bulk Import .csv')
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

    def give_max_layout(self):
        # Molecule Tab
        self.molecule_layout.removeItem(self.top_info_sublayout)
        self.molecule_layout.removeWidget(self.mol_display)
 
        self.mol_left_layout.insertItem(0, self.top_info_sublayout)
        self.mol_left_layout.insertWidget(1, self.mol_display)
        self.molecule_layout.insertWidget(0, self.lhs_title)
        self.lhs_title.show()
        self.molecule_layout.insertWidget(1, self.hlineTwo)
        self.hlineTwo.show()
        self.molecule_layout.insertItem(20, self.strechBox) #10 for buttons on bottom

        # Search Tab
        self.search_widget_index = 2
        self.search_layout.removeWidget(self.molLabel)
        self.search_layout.removeWidget(self.mol_result_display)
        self.left_search_layout.insertWidget(0, self.molLabel)
        self.left_search_layout.insertWidget(1, self.mol_result_display)        

    def give_normal_layout(self):
        try:
            # Molecule Tab
            self.mol_left_layout.removeItem(self.top_info_sublayout)
            self.mol_left_layout.removeWidget(self.mol_display)

            self.molecule_layout.removeWidget(self.lhs_title)
            self.lhs_title.hide()
            self.molecule_layout.removeWidget(self.hlineTwo)
            self.hlineTwo.hide()
            self.molecule_layout.removeItem(self.strechBox)

            # Search Tab
            self.left_search_layout.removeWidget(self.molLabel)
            self.left_search_layout.removeWidget(self.mol_result_display)

        except:
            print('Layout adjustment not needed.')

        # Molecule Tab
        self.molecule_layout.insertLayout(0, self.top_info_sublayout)
        self.molecule_layout.insertWidget(1, self.mol_display)

        # Search Tab
        self.search_widget_index = 0
        self.search_layout.insertWidget(0, self.molLabel)
        self.search_layout.insertWidget(1, self.mol_result_display)


    def mol_left_layout_check(self):
        if self.mol_left_layout.isEmpty():
            #masterLayout = self.molecule_layout
            pass
        else:
            self.molecule_layout.addStretch()
            #self.mol_overall_layout.addLayout(self.mol_left_layout)
            #self.mol_overall_layout.addLayout(self.molecule_layout)
            #masterLayout = self.mol_overall_layout
        #return masterLayout

    def lazy_compare(self,substring,string):
        for char in substring:
            if char in string:
                continue
            else:
                return False
        return True

    def search_database(self, search_text):
            # This needs to be rewitten as a switch

            self.readDatabase = pd.read_csv(defaultDB) #, index_col=0) #, encoding='mbcs')

            if search_text == '' or search_text == None:    
                self.selectedDatabase = self.readDatabase
                self.show_result(self.selectedDatabase, True)

            elif self.searchTypeSelection.currentText() == 'SMILES':
                if self.searchSubType.currentText() == 'Exact':
                    # Note: I should think about 'canonicalising' the SMILES on input to make comparisons easier and prevent duplications. I should think about this more: canon_SMILES = Chem.CanonSmiles(Molecule).
                    searchSMILES = search_text
                    row_index = self.readDatabase.index[self.readDatabase['SMILES'] == searchSMILES].tolist()
                    foundDataFrame = self.readDatabase.iloc[row_index]
                    foundDataFrame.reset_index(drop=True)
                    self.selectedDatabase = foundDataFrame
                    self.show_result(self.selectedDatabase, True)

                elif self.searchSubType.currentText() == 'Substructure':
                    try:
                        searchTest = MolFromSmiles(search_text)
                        if search_text != '' and search_text != None and searchTest != None:
                            smilesList = self.readDatabase['SMILES'].tolist()  
                            foundList = []
                            for smile in smilesList:
                                if len(smile) < len(search_text):
                                    continue

                                elif search_text in smile:
                                    foundList += [smile]
                                    continue

                                elif not self.lazy_compare(search_text,smile):
                                    continue

                                else:
                                    try:
                                        searchStructure = MolFromSmiles(smile)
                                        fullmatchList = Mol.GetSubstructMatches(searchStructure, searchTest)
                                        if len(fullmatchList) > 0:
                                            print('Substructure Match Found: ' + smile)
                                            foundList += [smile]
                                    except:
                                        print(f'Bad SMILES found in Database: {smile}')

                            print(foundList)
                            indexList = []
                            for foundMatch in foundList:
                                row_index = self.readDatabase.index[self.readDatabase['SMILES'] == foundMatch].tolist()
                                indexList += row_index
                            foundDataFrame = self.readDatabase.iloc[indexList]
                            print(foundDataFrame)
                            print('\n\n')
                            foundDataFrame.reset_index(drop=True)
                            self.selectedDatabase = foundDataFrame
                            self.show_result(self.selectedDatabase, True)
                    except Exception as e:
                        errorInfo = f"Enter Valid SMILES:\n\nDetailed info: {e}"
                        self.interactiveErrorMessage(errorInfo)

            elif self.searchTypeSelection.currentText() == 'Name':
                self.text_search_type('name',search_text)
                self.show_result(self.selectedDatabase, True)

            elif self.searchTypeSelection.currentText() == 'Project':
                self.text_search_type('proj',search_text)
                self.show_result(self.selectedDatabase, True)              

            elif self.searchTypeSelection.currentText() == 'MW':
                self.value_search_type('MW',search_text)
                self.show_result(self.selectedDatabase, True)

            elif self.searchTypeSelection.currentText() == 'Qdsc':
                self.value_search_type('Q_dsc',search_text)
                self.show_result(self.selectedDatabase, True)
                
            elif self.searchTypeSelection.currentText() == 'Tinit':
                self.value_search_type('initT',search_text)
                self.show_result(self.selectedDatabase, True)   

            elif self.searchTypeSelection.currentText() == 'Tonset':
                self.value_search_type('onsetT',search_text)
                self.show_result(self.selectedDatabase, True)               

            elif self.searchTypeSelection.currentText() == 'Td24':
                self.value_search_type('Td24',search_text)
                self.show_result(self.selectedDatabase, True)                   

            elif self.searchTypeSelection.currentText() == 'O.R.E.O.S. at >500 g':
                self.value_search_type('oreoLargeScale_val',search_text)
                self.show_result(self.selectedDatabase, True)

            elif self.searchTypeSelection.currentText() == 'Rule of Six':
                self.value_search_type('RoS_val',search_text)
                self.show_result(self.selectedDatabase, True)

            elif self.searchTypeSelection.currentText() == 'Oxygen Balance':
                self.value_search_type('OB_val',search_text)
                self.show_result(self.selectedDatabase, True)

            else:
                errorInfo = "How!? How have you gotten to this error message... I... Ugh... Report this to the developer okay?"
                self.interactiveErrorMessage(errorInfo)   

    def show_result(self, Database, resetIndex):
            #print(f'show_result run with args: Database={Database} and resetIndex={resetIndex}')
            print(f'\nshow_result run, resetIndex = {resetIndex}')
            layout = self.layout()
            if Database.empty:
                print('Database empty')
                errorInfo = "No matches found. Try a diffrent search?"
                self.interactiveErrorMessage(errorInfo)
                self.counter_label.setText("")
                self.result_label.setText("")
                self.mol_result_display.hide()
                self.molLabel.hide()
                self.make_plot_label.hide()
                self.make_plot_button.hide()
                self.select_x_values.hide()
                self.vs_label.hide()
                self.select_y_values.hide()
                try:
                    self.export_button.hide()
                    self.gen_multi_rep_button.hide()
                    self.edit_button.hide()
                    self.del_button.hide()
                    self.prev_button.hide()
                    self.next_button.hide()
                    self.first_button.hide()
                    self.last_button.hide()
                    self.results_table_Label.hide()
                    self.results_table.hide()
                    self.results_approval_warning.hide()
                except:
                    print('Not shown')
            if self.error_flag is not None:
                self.error_message.setText('')
                layout.removeWidget(self.error_message)
                self.error_flag = None
            if resetIndex == True:
                self.current_index = 0             
            if Database is not None and not Database.empty:
                print(f'Database not empty, current index = {self.current_index}\n')
                """try:
                    self.search_layout.removeWidget(self.results_table_Label)
                    self.search_layout.removeWidget(self.result_label)
                    self.search_layout.removeWidget(self.results_table)
                    self.search_layout.removeWidget(self.results_approval_warning)
                    self.search_layout.removeWidget(self.counter_label)
                    #self.search_layout.removaddLayout(prev_next_layout)
                except:
                    print('Not in layout for valid dataframe') """
                self.mol_result_display.show()
                self.molLabel.show()
                current_row = Database.iloc[self.current_index]
                dictRow = current_row.to_dict()
                print(f'The current row (as a dict) is: {dictRow}')
                readMolecule = thermalDexMolecule(**dictRow)
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
                self.results_approval_warning.setText(oreoApprovalWarning)


                try:
                    MW_as_fmtd_float = '{:.2f}'.format(float(current_row['MW']))
                except (TypeError,ValueError) as e:
                    MW_as_fmtd_float = ''
                    print(f"{e}\n\n{current_row['MW']}")
                try:
                    OB_val_as_fmtd_float = '{:.2f}'.format(float(current_row['OB_val']))
                except (TypeError,ValueError) as e:
                    OB_val_as_fmtd_float = ''
                    print(f"{e}\n\n{current_row['OB_val']}")
                try:
                    IS_val_as_fmtd_float = '{:.2f}'.format(float(current_row['IS_val']))
                except (TypeError,ValueError) as e:
                    IS_val_as_fmtd_float = ''
                    print(f"{e}\n\n{current_row['IS_val']}")
                try:
                    EP_val_as_fmtd_float = '{:.2f}'.format(float(current_row['EP_val']))
                except (TypeError,ValueError) as e:
                    EP_val_as_fmtd_float = ''
                    print(f"{e}\n\n{current_row['EP_val']}")

                result_text = f'''SMILES: {str(current_row['SMILES'])}<br>
                                  Name: {current_row['name']}<br>
                                  High Energy Groups: {current_row['HEG']}<br>
                                  Explosive Functional Groups: {current_row['EFG']}<br>
                                  mp: {current_row['mp']} to {current_row['mpEnd']} °C<br>
                                  MW: {MW_as_fmtd_float} g mol<sup>-1</sup><br>
                                  Q<sub>DSC</sub>: {current_row['Q_dsc'] if not pd.isna(current_row['Q_dsc']) else ' N/A '} {current_row['Qunits']}<br>
                                  T<sub>onset</sub>: {current_row['onsetT'] if not pd.isna(current_row['onsetT']) else ' N/A '} °C<br>
                                  T<sub>init</sub>: {current_row['initT'] if not pd.isna(current_row['initT']) else ' N/A '} °C<br>
                                  Oxygen Balance: {OB_val_as_fmtd_float} {current_row['OB_des']}<br>
                                  Rule of Six: {current_row['RoS_val']} {current_row['RoS_des']}<br>
                                  Impact Sensitivity: {IS_val_as_fmtd_float} {current_row['IS_des']}<br>
                                  Explosive Propagation: {EP_val_as_fmtd_float} {current_row['EP_des']}<br>
                                  Yosida Calculation Method Used: {current_row['yoshidaMethod']}<br>{formTd24}<br>
                                  Hammer Drop Test: {current_row['hammerDrop']}<br>Friction Test: {current_row['friction']}<br>
                                  Project: {current_row['proj']}'''
                
                self.result_label.setText(result_text)
                self.result_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
                self.result_label.setWordWrap(True)
                self.results_table.clearContents()
                oreoSmall = current_row['oreoSmallScale_des']
                oreoTens = current_row['oreoTensScale_des']
                oreoHunds = current_row['oreoHundsScale_des']
                oreoLarge = current_row['oreoLargeScale_des']

                smallEntry = QTableWidgetItem(oreoSmall)
                self.results_table.setItem(0, 0, smallEntry)
                classColor = self.getColorForValue(oreoSmall)
                smallEntry.setBackground(classColor)

                tensEntry = QTableWidgetItem(oreoTens)
                self.results_table.setItem(0, 1, tensEntry)
                classColor = self.getColorForValue(oreoTens)
                tensEntry.setBackground(classColor)

                hundsEntry = QTableWidgetItem(oreoHunds)
                self.results_table.setItem(0, 2, hundsEntry)
                classColor = self.getColorForValue(oreoHunds)
                hundsEntry.setBackground(classColor)

                largeEntry = QTableWidgetItem(oreoLarge)
                self.results_table.setItem(0, 3, largeEntry)
                classColor = self.getColorForValue(oreoLarge)
                largeEntry.setBackground(classColor)

                self.counter_label.setText(f"Result {self.current_index + 1} of {len(Database)}")
                self.edit_button.show()
                self.del_button.show()
                self.search_layout.insertWidget(5-self.search_widget_index, self.result_label)#5
                self.search_layout.insertWidget(6-self.search_widget_index, self.results_table_Label)#6
                self.search_layout.insertWidget(7-self.search_widget_index, self.results_table)#7
                self.search_layout.insertWidget(8-self.search_widget_index, self.results_approval_warning)#8
                self.search_layout.addWidget(self.counter_label)
                prev_next_layout = QHBoxLayout()
                prev_next_layout.addWidget(self.first_button)
                prev_next_layout.addWidget(self.prev_button)
                prev_next_layout.addWidget(self.next_button)
                prev_next_layout.addWidget(self.last_button)
                self.search_layout.addLayout(prev_next_layout)
                self.export_button.show()
                self.gen_multi_rep_button.show()
                self.edit_button.show()
                self.del_button.show()
                self.make_plot_label.show()
                self.make_plot_button.show()
                self.plot_database = Database
                self.plot_current_value = dictRow
                self.select_x_values.show()
                self.vs_label.show()
                self.select_y_values.show()
                self.results_table_Label.show()
                self.results_table.show()
                self.first_button.show()
                self.prev_button.show()
                self.next_button.show()
                self.last_button.show()

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
            self.show_result(self.selectedDatabase, False)

    def next_result(self, Database):
        if Database is not None:
            if self.current_index < len(Database) - 1:
                #search_layout.removeWidget(self.mol_result_display)
                #search_layout.removeWidget(self.molLabel)
                self.current_index += 1
                self.show_result(self.selectedDatabase, False)

    def first_result(self):
        self.current_index = 0
        self.show_result(self.selectedDatabase, False)

    def last_result(self, Database):
        self.current_index = len(Database)-1
        self.show_result(self.selectedDatabase, False)

    def subFunctionForDeletion(self, Entry, Database):
        self.delCurrentEntry(Entry, Database)
        self.search_database(self.entry_search.text())

    def exportSearchResults(self, Database):
        savefile, _ = QFileDialog.getSaveFileName(self, "Export Search Results as CSV", "ThermalDexResults.csv", "CSV Files (*.csv)")
        if savefile:
            Database.to_csv(savefile, index=False)

    def createMultiReport(self, Database):
        string_of_report = '''
                            <table>
                                <thead>
                                    <tr>
                                        <th>Structure</th>
                                        <th>O.R.E.O. Risk <5 g</th>
                                        <th>O.R.E.O. Risk 5 to 100 g</th>
                                        <th>O.R.E.O. Risk 100 to 500 g</th>
                                        <th>O.R.E.O. Risk >500 g</th>
                                        <th>No. of High Energy Groups</th>
                                        <th>No. of Explosive Groups</th>
                                    </tr>
                                </thead>
                                <tbody>
                            ''' # removed <th>Name</th> and <th>SMILES</th> for clarity on page
        
        string_dsc_report = '''
                            <table>
                                <thead>
                                    <tr>
                                        <th>Structure</th>
                                        <th>Imact Sensitivity</th>
                                        <th>Explosive Propagation</th>
                                        <th>T<sub>D24 (°C)</sub></th>
                                    </tr>
                                </thead>
                                <tbody>
                            ''' 
        
        for index, row in Database.iterrows():
            dictRow = row.to_dict()
            print(dictRow)
            repMolecule = thermalDexMolecule(**dictRow)
            repMolecule.genMol()
            imageData = repMolecule.molToBytes()
            dataURL = 'data:image/png;base64,' + imageData
            string_of_report += f'''
                                    <tr>
                                        <td><img src="{dataURL}" alt="HTML image test" style="object-fit: cover;"/></td>
                                        <td>{repMolecule.oreoSmallScale_des}</td>
                                        <td>{repMolecule.oreoTensScale_des}</td>
                                        <td>{repMolecule.oreoHundsScale_des}</td>
                                        <td>{repMolecule.oreoLargeScale_des}</td>
                                        <td>{repMolecule.HEG}</td>
                                        <td>{repMolecule.EFG}</td>
                                    </tr>
                                ''' # removed  <td>{repMolecule.name}</td> and <td>{repMolecule.SMILES}</td> for clarity
            
            string_dsc_report += f'''
                                    <tr>
                                        <td><img src="{dataURL}" alt="HTML image test" style="object-fit: cover;"/></td>
                                        <td>{repMolecule.IS_des if not pd.isna(repMolecule.IS_des) and repMolecule.IS_des != None and repMolecule.IS_des != 'None' else 'N/A'}</td>
                                        <td>{repMolecule.EP_des if not pd.isna(repMolecule.EP_des) and repMolecule.EP_des != None and repMolecule.EP_des != 'None' else 'N/A'}</td>
                                        <td>{round(repMolecule.Td24,1) if isinstance(repMolecule.Td24,float) and not pd.isna(repMolecule.Td24) else 'N/A'}</td>
                                    </tr>
                                ''' 
            
        string_of_report += f'''
                                </tbody>
                            </table>
                            '''  
                    
        string_dsc_report += f'''
                                </tbody>
                            </table>
                            '''  
        multiReportCreation(string_of_report,string_dsc_report)


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
        imported_df = pd.DataFrame()
        if import_file:
                  overrideDialouge = QMessageBox()
                  overrideDialouge.setIcon(QMessageBox.Information)
                  overrideDialouge.setWindowTitle("ThermalDex - Info Box")
                  overrideDialouge.setText("Run With Override Check?")
                  overrideDialouge.setInformativeText("By default ThermalDex checks each compound to ensure existing compounds in the database are not overwritten. This may results in multiple warning dialouges. Would you like to leave these checks in place? (recommend 'yes' unless you have a reason not to).")
                  overrideDialouge.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
                  returnValue = overrideDialouge.exec()
                  if returnValue == QMessageBox.Yes:
                      import_override_protection = True
                  else:
                      import_override_protection = False
                  importedDB = pd.read_csv(import_file)
                  invalid_SMILES = []
                  for index, row in importedDB.iterrows():
                    dictRow = row.to_dict()
                    readMolecule = thermalDexMolecule(**dictRow)
                    readMolecule.genMol()
                    invalid_SMILES.append(self.checkIfSMILESAreValid(readMolecule))
                    if readMolecule.mol is not None:
                        # Calculate Core Properties
                        self.genCoreValuesFromMol(readMolecule)
                    sql_data = self.writeToDatabase(readMolecule, defaultDB, sqlflag = False, override_protection = import_override_protection)
                    imported_df = pd.concat([imported_df,cleanMolDataFrame(readMolecule)])

                  invalid_SMILES = filter(None, invalid_SMILES)
                  print(f'The following SMILES were found to be invalid: {invalid_SMILES}')

                  self.sqlite_db_implementation(sql_data)
                  self.import_qmsg = QMessageBox()
                  self.import_qmsg.setIcon(QMessageBox.Information)
                  self.import_qmsg.setWindowTitle("ThermalDex - Info Box")
                  self.import_qmsg.setText("Import Complete")
                  self.import_qmsg.setInformativeText("Import has completed. Do you want to generate a summary pdf?")
                  self.import_qmsg.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
                  returnValue = self.import_qmsg.exec()
                  if returnValue == QMessageBox.Yes:
                      self.createMultiReport(imported_df)


    def text_search_type(self, search_column, search_text):
        if self.searchSubType.currentText() == 'is':
            searchName = search_text
            row_index = self.readDatabase.index[self.readDatabase[search_column] == searchName].tolist()
            foundDataFrame = self.readDatabase.iloc[row_index]
            foundDataFrame.reset_index(drop=True)
            self.selectedDatabase = foundDataFrame
            #show_result(self, foundDataFrame, True)

        elif self.searchSubType.currentText() == 'contains':
            searchName = search_text
            nameList = self.readDatabase[search_column].tolist()
            cleanNameList = []
            foundList = []
            for name in nameList:
                if name == name: # a check to remove nan values
                    cleanNameList += [name]
            for name in cleanNameList:
                print(name)
                if searchName.lower() in name.lower():
                    foundList += [name]
            foundList = list(dict.fromkeys(foundList))
            indexList = []
            for foundMatch in foundList:
                row_index = self.readDatabase.index[self.readDatabase[search_column] == foundMatch].tolist()
                indexList += row_index
            foundDataFrame = self.readDatabase.iloc[indexList]
            #print(foundDataFrame)
            #print('\n\n')
            foundDataFrame.reset_index(drop=True)
            self.selectedDatabase = foundDataFrame
            #show_result(self, foundDataFrame, True)

        elif self.searchSubType.currentText() == 'is not':
            searchName = search_text
            row_index = self.readDatabase.index[self.readDatabase[search_column] != searchName].tolist()
            foundDataFrame = self.readDatabase.iloc[row_index]
            foundDataFrame.reset_index(drop=True)
            self.selectedDatabase = foundDataFrame

        elif self.searchSubType.currentText() == 'does not contain':
            searchName = search_text
            nameList = self.readDatabase[search_column].tolist()
            cleanNameList = []
            foundList = []
            for name in nameList:
                if name == name: # check for nan values!
                    if searchName.lower() in name.lower():
                        foundList += [name]
            foundList = list(dict.fromkeys(foundList))
            indexList = []
            if len(foundList) == 0:
                foundDataFrame = self.readDatabase
            else:    
                for foundMatch in foundList:
                    row_index = self.readDatabase.index[self.readDatabase[search_column] != foundMatch].tolist()
                    indexList += row_index
                foundDataFrame = self.readDatabase.iloc[indexList]
                foundDataFrame.reset_index(drop=True)
            self.selectedDatabase = foundDataFrame


    def value_search_operator_selection(self, valueList, searchName):
        foundList = []
        if self.searchSubType.currentText() == '<':
            for value in valueList:
                #print(value)
                if value < float(searchName):
                    foundList += [value]

        elif self.searchSubType.currentText() == '>':
            for value in valueList:
                #print(value)
                if value > float(searchName):
                    foundList += [value]

        elif self.searchSubType.currentText() == '</=':
            for value in valueList:
                #print(value)
                if value <= float(searchName):
                    foundList += [value]

        elif self.searchSubType.currentText() == '>/=':
            for value in valueList:
                #print(value)
                if value >= float(searchName):
                    foundList += [value]
        return foundList

    def value_search_type(self, search_column, search_text):
        if self.searchSubType.currentText() == '=':
            try:
                searchName = float(search_text)
                row_index = self.readDatabase.index[self.readDatabase[search_column] == searchName].tolist()
                foundDataFrame = self.readDatabase.iloc[row_index]
                foundDataFrame.reset_index(drop=True)
                self.selectedDatabase = foundDataFrame
                #show_result(self, foundDataFrame, True)
            except:
                errorInfo = "Enter a Valid number"
                self.interactiveErrorMessage(errorInfo) 

        else:
            try:
                searchName = search_text
                foundList = self.value_search_operator_selection(self.readDatabase[search_column].tolist(),searchName)
                print(foundList)
                foundList = list(dict.fromkeys(foundList))
                indexList = []
                for foundMatch in foundList:
                    row_index = self.readDatabase.index[self.readDatabase[search_column] == foundMatch].tolist()
                    indexList += row_index
                foundDataFrame = self.readDatabase.iloc[indexList]
                #print(foundDataFrame)
                #print('\n\n')
                foundDataFrame.reset_index(drop=True)
                self.selectedDatabase = foundDataFrame
                #show_result(self, foundDataFrame, True)
            except:
                errorInfo = "Enter a Valid number"
                self.interactiveErrorMessage(errorInfo)               

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

    def check_if_dsc_performed(self):
        if self.dsc_check.isChecked():
            self.Qdsc_input.setReadOnly(True)
            self.TE_input.setReadOnly(True)
            self.Tinit_input.setReadOnly(True)
            self.Qdsc_input.setStyleSheet("QLineEdit"
                                "{"
                                "background : lightgrey;"
                                "}") 
            self.TE_input.setStyleSheet("QLineEdit"
                                "{"
                                "background : lightgrey;"
                                "}")
            self.Tinit_input.setStyleSheet("QLineEdit"
                                "{"
                                "background : lightgrey;"
                                "}")
            self.DSC_run_flag = True
        else:          
            self.Qdsc_input.setReadOnly(False)
            self.TE_input.setReadOnly(False)
            self.Tinit_input.setReadOnly(False)
            self.Qdsc_input.setStyleSheet("QLineEdit"
                                "{"
                                "background : white;"
                                "}") 
            self.TE_input.setStyleSheet("QLineEdit"
                                "{"
                                "background : white;"
                                "}")
            self.Tinit_input.setStyleSheet("QLineEdit"
                                "{"
                                "background : white;"
                                "}")
  

    def changeTabForEditing(self):
        #try:
        self.resetToDefaultState()
        self.approval_needed.hide()
        self.lhs_title.setText("Title")
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
        if readMolecule.noDSCPeak == '':
            self.dsc_check.setChecked(False)
        else:
            self.dsc_check.setChecked(True)
        self.check_if_dsc_performed()
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

        #fileCounter = self.countFiles(defaultDB)
        #self.filesCount.setText(f"{str(fileCounter)} Attached Files")
        self.fileCounterUpdate()
        self.attachedFilesLabel.show()
        self.filesCount.show()
        self.attach_button.show()

        if self.check_if_oreos_need_approval(readMolecule) == 'Show Approval Message':
            self.approval_needed.show()
            self.lhs_title.setText("<h3>Title</h3>")
    
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
            if readMolecule.measuredTD24:
                d24Str = "{:.1f}".format(readMolecule.empiricalTD24)
                if int(ambertd24limit) >= readMolecule.empiricalTD24 > int(redtd24limit):
                    self.Td24Label.setText(f"T<sub>D24</sub>: <b style='color: orange;'> {d24Str} °C</b> (measured)")
                elif readMolecule.empiricalTD24 <= int(redtd24limit):
                    self.Td24Label.setText(f"<h3 style='color: red;'>T<sub>D24</sub>: {d24Str} °C</h3> (measured)")
                    self.approval_needed.show()
                    self.lhs_title.setText("<h3>Title</h3>")
                else:
                    self.Td24Label.setText('T<sub>D24</sub>: <b>' + d24Str + ' °C</b> (measured)')  

            else:
                d24Str = "{:.1f}".format(readMolecule.Td24)
                if int(ambertd24limit) >= readMolecule.Td24 > int(redtd24limit):
                    self.Td24Label.setText(f"T<sub>D24</sub>: <b style='color: orange;'> {d24Str} °C</b>")
                elif readMolecule.Td24 <= int(redtd24limit):
                    self.Td24Label.setText(f"<h3 style='color: red;'>T<sub>D24</sub>: {d24Str} °C</h3>")
                    self.approval_needed.show()
                    self.lhs_title.setText("<h3>Title</h3>")
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

    def openComments(self):
        comment_location = self.set_comment_location(database=defaultDB)
        if comment_location == None:
            errorInfo = f"Please Enter A Valid SMILES Before Commenting."
            self.interactiveErrorMessage(errorInfo)
        else:
            self.commentsWindow = CommentsBox(comments_location=comment_location)
            self.commentsWindow.submitClicked.connect(self.fileCounterUpdate)
            #self.commentsWindow.exec_()
            self.commentsWindow.show()
            self.commentsWindow.raise_()
            self.commentsWindow.activateWindow()

    def fileCounterUpdate(self) -> None:
        fileCounter = self.countFiles(defaultDB)
        self.filesCount.setText(f"{str(fileCounter)} Attached Files")
        return None

    def openTd24Override(self):
        self.overrideWindow = Td24OverrideWindow()
        self.overrideWindow.submitClicked.connect(self.on_sub_window_confirm)
        #self.overrideWindow.exec_()
        self.overrideWindow.show()
        self.overrideWindow.raise_()
        self.overrideWindow.activateWindow()

    def on_sub_window_confirm(self,onBool,onVal):
        self.overrideOn = onBool
        self.overideValue = onVal
        print(f'{onVal=}')
        self.render_molecule()

    def set_comment_location(self, database) -> str | None:
        folder = self.findFolder(database=database)
        if folder != '' and folder != 'nan':
            return folder
        
        else:
            newMolecule = self.render_molecule()
            if newMolecule != None:
                newfolder = self.findFolder(database=database)
                return newfolder
            else:
                return None

    def findFolder(self, database) -> str:
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
        self.overrideOn = False
        self.overrideValue = None

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
        #try:
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
        #except:
         #   window.showErrorMessage("Generating Memo PDF from given values.")
            
    def sqlite_db_implementation(self, pandas_data):
        sqlite_db_connection = sqlite3.connect('./_core/ThermalDex.db')
        sqlite_db_cursor = sqlite_db_connection.cursor()
        now = datetime.now()
        neatNow = now.strftime("%d-%b-%Y_%H-%M-%S")
        #column_headers = tuple(pandas_data.columns.values)
        #sqlite_db_cursor.execute(f"CREATE TABLE {neatNow}{column_headers}")
        pandas_data.to_sql(name=neatNow, con=sqlite_db_connection)

    def writeToDatabase(self, molecule, Database, sqlflag=True, override_protection=True):
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

        if selectedMolData['SMILES'][0] in storedData.index and override_protection == True:
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

        elif selectedMolData['SMILES'][0] in storedData.index and override_protection == False:
            storedData.update(selectedMolData)
            outputData = storedData
            row_index = checkData.index[checkData['SMILES'] == selectedMolData['SMILES'][0]].tolist()
            print(row_index)
            folderInfo = checkData['dataFolder'][row_index[0]]
            print(folderInfo)
            output_index = outputData.index[checkData['SMILES'] == selectedMolData['SMILES'][0]].tolist()
            outputData['dataFolder'][output_index[0]] = folderInfo

        else:
            outputData = pd.concat([storedData, selectedMolData])

        outputData['SMILES'] = outputData.index
        outputData = outputData[ ['SMILES'] + [ col for col in outputData.columns if col != 'SMILES' ] ]
        print(outputData)
        outputData.to_csv(Database, index=False)
        sql_data = outputData.drop('SMILES', axis=1)

        if sqlflag == True:
            self.sqlite_db_implementation(sql_data)
        
        return sql_data

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
        
    def get_dsc_checkstate(self):
        if self.dsc_check.isChecked():
            return 'True'
        else:
            return ''

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
        noDSCPeak = self.get_dsc_checkstate()
        measuredTD24 = self.overrideOn
        empiricalTD24 = self.overideValue
        print(f'{measuredTD24=}\n{empiricalTD24=}\n')

        # Create an RDKit molecule from the SMILES string
        addedMolecule = thermalDexMolecule(SMILES=smiles, name=name, mp=mp, mpEnd=mpEnd, Q_dsc=Qdsc, Qunits=QUnits, onsetT=TE, initT=Tinit, proj=proj, hammerDrop=hammerDrop, friction=friction, dataFolder=dataFolder, yoshidaMethod=yoshidaMethod, noDSCPeak=noDSCPeak, measuredTD24=measuredTD24, empiricalTD24=empiricalTD24)
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
            try:
                realtest = molecule.molToBytes()
            except:
                molecule.mol = None
                self.showErrorMessage(f'Invalid SMILES detected: {molecule.SMILES}')
                print(f'Invalid SMILES: {molecule.SMILES}')
                return molecule.SMILES

        else:
            molecule.mol = None

    def displayTheMolecule(self, molecule, display):
        # Make Pixmap Image to Display.
        pixmap = molecule.molToQPixmap()
        scene = QGraphicsScene()
        scene.addPixmap(pixmap)
        display.setScene(scene)
        molecule.molPixmap = None

    def render_molecule(self) -> thermalDexMolecule | None:
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
                self.lhs_title.setText("<h3>Title</h3>")
            else:
                self.lhs_title.setText("Title")

            if addedMolecule.isStr != None:
                self.ISLabel.setText('Yoshida Impact Sensitivity: ' + addedMolecule.isStr + addedMolecule.IS_des)
            if addedMolecule.epStr != None:
                self.EPLabel.setText('Yoshida Explosive Propagation: ' + addedMolecule.epStr + addedMolecule.EP_des)
            if addedMolecule.Td24 != '' and addedMolecule.Td24 != 'nan' and addedMolecule.Td24 != None and not addedMolecule.measuredTD24:
                d24Str = "{:.1f}".format(addedMolecule.Td24)
                if int(ambertd24limit) >= addedMolecule.Td24 > int(redtd24limit):
                    self.Td24Label.setText(f"T<sub>D24</sub>: <b style='color: orange;'> {d24Str} °C</b>")
                elif addedMolecule.Td24 <= int(redtd24limit):
                    self.Td24Label.setText(f"<h3 style='color: red;'>T<sub>D24</sub>: {d24Str} °C</h3>")
                    self.approval_needed.show()
                    self.lhs_title.setText("<h3>Title</h3>")
                else:
                    self.Td24Label.setText('T<sub>D24</sub>: <b>' + d24Str + ' °C</b>')
                    self.lhs_title.setText("Title")
                
            elif addedMolecule.measuredTD24:
                print(f'{addedMolecule.empiricalTD24=}')
                d24Str = "{:.1f}".format(addedMolecule.empiricalTD24)
                if int(ambertd24limit) >= addedMolecule.empiricalTD24 > int(redtd24limit):
                    self.Td24Label.setText(f"T<sub>D24</sub>: <b style='color: orange;'> {d24Str} °C</b> (measured)")
                elif addedMolecule.empiricalTD24 <= int(redtd24limit):
                    self.Td24Label.setText(f"<h3 style='color: red;'>T<sub>D24</sub>: {d24Str} °C</h3> (measured)")
                    self.approval_needed.show()
                    self.lhs_title.setText("<h3>Title</h3>")
                else:
                    self.Td24Label.setText('T<sub>D24</sub>: <b>' + d24Str + ' °C</b> (measured)')
                    self.lhs_title.setText("Title")

            if addedMolecule.noDSCPeak == '' and addedMolecule.onsetT != 'nan' and addedMolecule.onsetT != '' and addedMolecule.onsetT != None and addedMolecule.onsetT <= 24.99:
                self.error_message = QLabel('Sub-Ambient Onset Temperature? Are you sure? Yoshida values will not be accurate.')
                layout.addWidget(self.error_message)
                self.error_flag = 201
            if addedMolecule.noDSCPeak == '' and addedMolecule.initT != 'nan' and addedMolecule.initT != '' and addedMolecule.initT != None and addedMolecule.initT <= 24.99:
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
            return None

    def recal_db(self):
        db_as_df = pd.read_csv(defaultDB)
        for index, row in db_as_df.iterrows():
            dictRow = row.to_dict()
            readMolecule = thermalDexMolecule(**dictRow)
            readMolecule.genMol()
            sql_data = self.writeToDatabase(readMolecule, defaultDB, sqlflag = False, override_protection = False)

            self.sqlite_db_implementation(sql_data)
            QMessageBox.information(self, "Database Recalculation", "Recaluation of the entire database has been completed.")


if __name__ == '__main__':
    #with open('./_core/ThermalDex.log', 'w', encoding='utf-8') as logFile:
     #   with redirect_stdout(logFile):
    defaultDB, highEnergyGroups, expEnergyGroups, yoshidaMethod, qdscUnits, ambertd24limit, redtd24limit, oreohazardlimit, oreohazardwarningscale = altImportConfig()
    font = QFont("Segoe UI")
    font.setPixelSize(11)
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('.\\_core\\ThermalDexIcon.ico'))
    app.setFont(font, "QLabel")
    app.setFont(font, "QLineEdit")
    app.setFont(font, "QComboBox")
    app.setFont(font, "QTableWidget")
    app.setFont(font, "QTableWidget.horizontalHeader")
    app.setFont(font, "QTabWidget")
    app.setFont(font, "QWidget")
    app.setFont(font, "QPushButton")
    screen = app.primaryScreen()
    dpi = screen.physicalDotsPerInch()
    print(f'Screen DPI: {dpi}')
    window = MolDrawer()
    #window.layout().setSizeConstraint(QLayout.SetFixedSize
    #window.showMaximized()
    window.show()
    window.raise_()
    window.activateWindow()
    sys.exit(app.exec_())
