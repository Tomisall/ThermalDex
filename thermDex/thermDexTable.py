from PyQt5.QtWidgets import (QSpinBox, QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QListWidget, QTextEdit, QLabel, QListWidgetItem, QPushButton, QShortcut,QToolBar,QAction,QDialog,QFontDialog, QGridLayout, QMessageBox, QMenu, QTableView)
from PyQt5.QtGui import QFont, QFontDatabase, QImage, QImageReader, QTextDocument, QPixmap, QPainter, QPen, QKeySequence, QIcon, QTextFormat, QTextCharFormat, QTextTable, QContextMenuEvent, QTextTableFormat
from PyQt5.QtCore import Qt, QUrl, QFile, QMimeData, QTranslator, QSize, QPoint, QAbstractTableModel, QVariant, QModelIndex, pyqtProperty, pyqtSlot
import pandas as pd

class Table(QDialog):
    def __init__(self,parent = None):
        QDialog.__init__(self, parent)
 
        self.parent = parent
 
        self.initUI()
 
    def initUI(self):
 
        # Rows
        rowsLabel = QLabel("Rows: ",self)
 
        self.rows = QSpinBox(self)
 
        # Columns
        colsLabel = QLabel("Columns",self)
 
        self.cols = QSpinBox(self)
 
        # Cell spacing (distance between cells)
        spaceLabel = QLabel("Cell spacing",self)
 
        self.space = QSpinBox(self)
 
        # Cell padding (distance between cell and inner text)
        padLabel = QLabel("Cell padding",self)
 
        self.pad = QSpinBox(self)
 
        self.pad.setValue(10)
 
        # Button
        insertButton = QPushButton("Insert",self)
        insertButton.clicked.connect(self.insert)
 
        # Layout
        layout = QGridLayout()
 
        layout.addWidget(rowsLabel,0,0)
        layout.addWidget(self.rows,0,1)
 
        layout.addWidget(colsLabel,1,0)
        layout.addWidget(self.cols,1,1)
 
        layout.addWidget(padLabel,2,0)
        layout.addWidget(self.pad,2,1)
 
        layout.addWidget(spaceLabel,3,0)
        layout.addWidget(self.space,3,1)
 
        layout.addWidget(insertButton,4,0,1,2)
 
        self.setWindowTitle("Insert Table")
        self.setGeometry(300,300,200,100)
        self.setLayout(layout)
 
    def insert(self):
 
        cursor = self.parent.textCursor()#self.parent.text.textCursor()
 
        # Get the configurations
        rows = self.rows.value()
 
        cols = self.cols.value()
 
        if not rows or not cols:
 
            popup = QMessageBox(QMessageBox.Warning,
                                      "Parameter error",
                                      "Row and column numbers may not be zero!",
                                      QMessageBox.Ok,
                                      self)
            popup.show()
 
        else:
 
            padding = self.pad.value()
 
            space = self.space.value()
 
            # Set the padding and spacing
            fmt = QTextTableFormat()
 
            fmt.setCellPadding(padding)
 
            fmt.setCellSpacing(space)
 
            # Inser the new table
            cursor.insertTable(rows,cols,fmt)
 
            self.close()


class DataFrameModel(QAbstractTableModel):
    DtypeRole = Qt.UserRole + 1000
    ValueRole = Qt.UserRole + 1001

    def __init__(self, df=pd.DataFrame(), parent=None):
        super(DataFrameModel, self).__init__(parent)
        self._dataframe = df

    def setDataFrame(self, dataframe):
        self.beginResetModel()
        self._dataframe = dataframe.copy()
        self.endResetModel()

    def dataFrame(self):
        return self._dataframe

    dataFrame = pyqtProperty(pd.DataFrame, fget=dataFrame, fset=setDataFrame)

    pyqtSlot(int, Qt.Orientation, result=str)
    def headerData(self, section: int, orientation: Qt.Orientation, role: int = Qt.DisplayRole):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self._dataframe.columns[section]
            else:
                return str(self._dataframe['Reactant Name'][section])
        return QVariant()

    def rowCount(self, parent=QModelIndex()):
        if parent.isValid():
            return 0
        return len(self._dataframe.index)

    def columnCount(self, parent=QModelIndex()):
        if parent.isValid():
            return 0
        return self._dataframe.columns.size

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid() or not (0 <= index.row() < self.rowCount() \
            and 0 <= index.column() < self.columnCount()):
            return QVariant()
        row = self._dataframe.index[index.row()]
        col = self._dataframe.columns[index.column()]
        dt = self._dataframe[col].dtype

        val = self._dataframe.iloc[row][col]
        if role == Qt.DisplayRole:
            return str(val)
        elif role == DataFrameModel.ValueRole:
            return val
        if role == DataFrameModel.DtypeRole:
            return dt
        return QVariant()

    def roleNames(self):
        roles = {
            Qt.DisplayRole: b'display',
            DataFrameModel.DtypeRole: b'dtype',
            DataFrameModel.ValueRole: b'value'
        }
        return roles