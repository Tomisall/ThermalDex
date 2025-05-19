from PyQt5.QtWidgets import (QCheckBox, QLineEdit, QSpinBox, QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QListWidget, QTextEdit, QLabel, QListWidgetItem, QPushButton, QShortcut,QToolBar,QAction,QDialog,QFontDialog, QGridLayout, QMessageBox, QMenu, QTableView)
from PyQt5.QtGui import QFont, QFontDatabase, QImage, QImageReader, QTextDocument, QPixmap, QPainter, QPen, QKeySequence, QIcon, QTextFormat, QTextCharFormat, QTextTable, QContextMenuEvent
from PyQt5.QtCore import Qt, QUrl, QFile, QMimeData, QTranslator, QSize, QPoint, QAbstractTableModel, QVariant, QModelIndex, pyqtProperty, pyqtSlot
import os
import sys
from rdkit import Chem
from rdkit.Chem import Draw, rdmolfiles, rdChemReactions
import xml.etree.ElementTree as ET
import tempfile
import shutil
from pyqt_color_picker import ColorPickerDialog

class TextEdit(QTextEdit):
    def __init__(self):
        super(QTextEdit, self).__init__()
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.context)

    def canInsertFromMimeData(self, source: QMimeData) -> bool:
        # Check if the mime data contains an image, URLs, or fallback to default behavior
        return source.hasImage() or source.hasFormat('application/x-qt-windows-mime;value="MRV"') or source.hasUrls() or source.hasFormat("text/plain") or source.hasFormat("chemical/x-cml") or source.hasFormat("chemical/x-mdl-molfile") or source.hasFormat("chemical/x-daylight-smiles") or super().canInsertFromMimeData(source)

    def insertFromMimeData(self, source: QMimeData):
        print("Available MIME formats:", source.formats()) 
        if source.hasImage():
            # Handle dropped image
            image = source.imageData()
            self._saveImageWithTempName(image=image)
        elif source.hasFormat('application/x-qt-windows-mime;value="MRV"'):
            print('MRV found')
            xml_data = source.text()#source.data("application/xml").data().decode("utf-8")
            self._insertMoleculeFromXML(xml_data)
        elif source.hasUrls():
            # Handle dropped URLs
            for url in source.urls():
                local_file = url.toLocalFile()
                if not os.path.exists(local_file):
                    continue
                # Check if the URL points to an image
                suffix = os.path.splitext(local_file)[1].lower()[1:].encode('utf-8')
                if suffix in QImageReader.supportedImageFormats():
                    backend_file = self._copyFileToBackend(local_file)
                    self.dropImage(QUrl(backend_file), QImage(backend_file))
                else:
                    self.dropTextFile(url)
        elif source.hasFormat("application/xml"):
            print('application/xml found')
            xml_data = source.data("application/xml").data().decode("utf-8")
            self._insertMoleculeFromXML(xml_data)
        elif source.hasFormat("text/xml"):
            print('text/xml found')
            xml_data = source.data("text/xml").data().decode("utf-8")
            self._insertMoleculeFromXML(xml_data)
        elif source.hasFormat("text/plain"):
            print('text found')
            text = source.text()
            if text.startswith("<?xml"):
                print('xml found')
                #xml_data = source.data("text/xml").data().decode("utf-8")
                self._insertMoleculeFromXML(text)
            else:
                super().insertFromMimeData(source)
        elif source.hasFormat("chemical/x-cml"):
            print('chem x-cml')
            cml_data = source.data("chemical/x-cml").data().decode("utf-8")
            self._insertMoleculeFromCML(cml_data)
        elif source.hasFormat("chemical/x-mdl-molfile"):
            print('chem mol')
            mol_data = source.data("chemical/x-mdl-molfile").data().decode("utf-8")
            self._insertMolecule(mol_data, format="mol")
        elif source.hasFormat("chemical/x-daylight-smiles"):
            print('chem smiles')
            mol_data = source.data("chemical/x-daylight-smiles").data().decode("utf-8")
            self._insertMolecule(mol_data, format="smiles")
        else:
            print(source.formats)
            # Default behavior
            super().insertFromMimeData(source)

    def dropImage(self, url: QUrl, image):
        if not image.isNull():
            # Add the image as a resource to the document and insert it
            self.document().addResource(QTextDocument.ImageResource, url, image)
            self.textCursor().insertImage(url.toString())

    def dropTextFile(self, url: QUrl):
        # Read and insert the content of the text file
        file_path = url.toLocalFile()
        try:
            with open(file_path, 'r', encoding='utf-8') as file:
                content = file.read()
                self.textCursor().insertText(content)
        except:
            self.textCursor().insertText(str(file_path))

    def _insertMoleculeFromXML(self, xml_data: str):
        try:
            if rdChemReactions.MrvBlockIsReaction(xml_data):
                print('Rxn found')
                rxn = rdChemReactions.ReactionFromMrvBlock(xml_data)
                self._renderReactionToImage(rxn=rxn)
            else:
                mol = rdmolfiles.MolFromMrvBlock(xml_data)
                self._renderMoleculeToImage(mol=mol)
        except:
            return

    def _insertMoleculeFromCML(self, cml_data: str):
        # Parse the CML (XML) data
        try:
            mol_block = self._extractMolBlockFromCML(cml_data)
            if mol_block:
                self._insertMolecule(mol_block, format="mol")
            else:
                self.textCursor().insertText("[Invalid CML data]")
        except Exception as e:
            self.textCursor().insertText(f"[Error parsing CML: {e}]")

    def _extractMolBlockFromCML(self, cml_data: str) -> str:
        # Parse the CML data and extract the MolBlock (MDL MOL format)
        root = ET.fromstring(cml_data)
        for molecule in root.findall(".//molecule"):
            mol_block = molecule.find("molfile")
            if mol_block is not None:
                return mol_block.text.strip()
        return None
    
    def _insertMolecule(self, mol_data: str, format: str):
        # Parse molecule using RDKit
        if format == "mol":
            mol = Chem.MolFromMolBlock(mol_data)
            self._renderMolecule(mol=mol)
        elif format == "smiles":
            mol = Chem.MolFromSmiles(mol_data)
            self._renderMolecule(mol=mol)
        else:
            return

        if mol is None:
            self.textCursor().insertText("[Invalid Molecule]")
            return

    def _renderMoleculeToImage(self, mol) -> QImage:
        # Use RDKit to render the molecule to a PIL image
        pil_img = Draw.MolToImage(mol)
        self._saveImageWithTempName(image=pil_img)
        return
    
    def _saveImageWithTempName(self,image):
        temp_name = next(tempfile._get_candidate_names())
        image.save(f'./backend/downloads/{temp_name}.png', format="PNG")
        qimage = QImage(f'./backend/downloads/{temp_name}.png')
        file_url = f'./backend/downloads/{temp_name}.png'
        self.dropImage(url=QUrl(file_url),image=qimage)
        return
    
    def _renderReactionToImage(self, rxn) -> QImage:
        # Use RDKit to render the molecule to a PIL image
        pil_img = Draw.ReactionToImage(rxn)
        self._saveImageWithTempName(image=pil_img)
        return
    
    def _copyFileToBackend(self, file_path: str) -> str:
        temp_name = next(tempfile._get_candidate_names())
        destination_file = f'./backend/downloads/{temp_name}{os.path.splitext(file_path)[1].lower()}'
        shutil.copyfile(file_path,destination_file)
        return destination_file
    
    def context(self,pos):
 
        # Grab the cursor
        cursor = self.textCursor()
 
        # Grab the current table, if there is one
        table = cursor.currentTable()
 
        # Above will return 0 if there is no current table, in which case
        # we call the normal context menu. If there is a table, we create
        # our own context menu specific to table interaction
        if table:
 
            menu = QMenu(self)
 
            appendRowAction = QAction("Append row",self)
            appendRowAction.triggered.connect(lambda: table.appendRows(1))
 
            appendColAction = QAction("Append column",self)
            appendColAction.triggered.connect(lambda: table.appendColumns(1))
 
 
            removeRowAction = QAction("Remove row",self)
            removeRowAction.triggered.connect(self.removeRow)
 
            removeColAction = QAction("Remove column",self)
            removeColAction.triggered.connect(self.removeCol)
 
 
            insertRowAction = QAction("Insert row",self)
            insertRowAction.triggered.connect(self.insertRow)
 
            insertColAction = QAction("Insert column",self)
            insertColAction.triggered.connect(self.insertCol)
 
 
            mergeAction = QAction("Merge cells",self)
            mergeAction.triggered.connect(lambda: table.mergeCells(cursor))
 
            # Only allow merging if there is a selection
            if not cursor.hasSelection():
                mergeAction.setEnabled(False)
 
 
            splitAction = QAction("Split cells",self)
 
            cell = table.cellAt(cursor)
 
            # Only allow splitting if the current cell is larger
            # than a normal cell
            if cell.rowSpan() > 1 or cell.columnSpan() > 1:
 
                splitAction.triggered.connect(lambda: table.splitCell(cell.row(),cell.column(),1,1))
 
            else:
                splitAction.setEnabled(False)
 
 
            menu.addAction(appendRowAction)
            menu.addAction(appendColAction)
 
            menu.addSeparator()
 
            menu.addAction(removeRowAction)
            menu.addAction(removeColAction)
 
            menu.addSeparator()
 
            menu.addAction(insertRowAction)
            menu.addAction(insertColAction)
 
            menu.addSeparator()
 
            menu.addAction(mergeAction)
            menu.addAction(splitAction)
 
            # Convert the widget coordinates into global coordinates
            pos = self.mapToGlobal(pos)
 
            # Add pixels for the tool and formatbars, which are not included
            # in mapToGlobal(), but only if the two are currently visible and
            # not toggled by the user
 
            #if self.toolbar.isVisible():
            #  pos.setY(pos.y() + 45)
 
            #if self.formatbar.isVisible():
            #    pos.setY(pos.y() + 45)
 
            # Move the menu to the new position
            menu.move(pos)
 
            menu.show()
 
        else:
 
            event = QContextMenuEvent(QContextMenuEvent.Mouse,QPoint())
 
            self.contextMenuEvent(event)
 
    def removeRow(self):
 
        # Grab the cursor
        cursor = self.textCursor()
 
        # Grab the current table (we assume there is one, since
        # this is checked before calling)
        table = cursor.currentTable()
 
        # Get the current cell
        cell = table.cellAt(cursor)
 
        # Delete the cell's row
        table.removeRows(cell.row(),1)
 
    def removeCol(self):
 
        # Grab the cursor
        cursor = self.textCursor()
 
        # Grab the current table (we assume there is one, since
        # this is checked before calling)
        table = cursor.currentTable()
 
        # Get the current cell
        cell = table.cellAt(cursor)
 
        # Delete the cell's column
        table.removeColumns(cell.column(),1)
 
    def insertRow(self):
 
        # Grab the cursor
        cursor = self.textCursor()
 
        # Grab the current table (we assume there is one, since
        # this is checked before calling)
        table = cursor.currentTable()
 
        # Get the current cell
        cell = table.cellAt(cursor)
 
        # Insert a new row at the cell's position
        table.insertRows(cell.row(),1)
 
    def insertCol(self):
 
        # Grab the cursor
        cursor = self.textCursor()
 
        # Grab the current table (we assume there is one, since
        # this is checked before calling)
        table = cursor.currentTable()
 
        # Get the current cell
        cell = table.cellAt(cursor)
 
        # Insert a new row at the cell's position
        table.insertColumns(cell.column(),1)