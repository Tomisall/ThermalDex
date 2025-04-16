import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene, QFrame, QTableWidget, QTableWidgetItem, QTabWidget, QGraphicsPixmapItem,  QMessageBox #, QTableView, QToolTip
from PyQt5.QtGui import QPixmap, QColor, QIcon, QRegExpValidator #QValidator #, QCursor
from PyQt5.QtCore import Qt, QRegExp, pyqtSignal
from rdkit.Chem import Draw, rdmolfiles, Mol
import pyperclip


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
        # Set up the main layout
        layout = QVBoxLayout()
        self.mwLabelTest = QLabel("Paste from .mrv here!")
        layout.addWidget(self.mwLabelTest)
        self.mwLabel = ClickableLineEdit(self)
        self.mwLabel.clicked.connect(self.mrvToSMILES)
        layout.addWidget(self.mwLabel)
        self.othermwLabel = ClickableLineEdit(self)
        layout.addWidget(self.othermwLabel)
        self.setLayout(layout)
        self.setGeometry(100, 100, 250, 250)
        self.setWindowTitle('SMILES Gen')

    def mrvToSMILES(self):
        try:
            copyMolXML = pyperclip.paste()
            copyMolReal = rdmolfiles.MolFromMrvBlock(copyMolXML)
            #Draw.ShowMol(copyMolReal)
            smilesFromMarvin = rdmolfiles.MolToSmiles(copyMolReal)
            print(smilesFromMarvin)
            pyperclip.copy(smilesFromMarvin)

        except:
            print("Poo")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('.\\_core\\ThermalDexIcon.ico'))
    window = MolDrawer()
    window.show()
    sys.exit(app.exec_())