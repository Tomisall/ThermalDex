import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene, QFrame, QTableWidget, QTableWidgetItem, QTabWidget, QGraphicsPixmapItem,  QMessageBox, QComboBox #, QTableView, QToolTip
from PyQt5.QtGui import QPixmap, QColor, QIcon, QRegExpValidator #QValidator #, QCursor
from PyQt5.QtCore import Qt, QRegExp, pyqtSignal
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Mol, MolFromSmiles, MolFromSmarts, rdmolfiles
from thermDex.thermDexQt import *
from io import BytesIO
from numpy import log10
from pubchempy import get_compounds
from dataclasses import dataclass, field, asdict
from thermDex.thermDexMolecule import * #thermalDexMolecule
from thermDex.thermDexReport import *
import pandas as pd
import re
import pyperclip

try:
    import pyi_splash
    pyi_splash.close()
except:
    pass

if __name__ == '__main__':
    defaultDB, highEnergyGroups, expEnergyGroups = importConfig()
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon('.\\_core\\ThermalDexIcon.ico'))
    window = MolDrawer()
    window.show()
    sys.exit(app.exec_())
