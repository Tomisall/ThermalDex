import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, \
    QGraphicsView, QGraphicsScene, QTabWidget

from rdkit import Chem
from rdkit.Chem import Draw
from PyQt5.QtGui import QPixmap
from io import BytesIO


class MolDrawer(QWidget):
    def __init__(self):
        super().__init__()

        self.init_ui()

    def init_ui(self):
        # Set up the main layout
        layout = QVBoxLayout()

        # Create a tab widget
        tab_widget = QTabWidget()

        # Tab for molecule rendering
        molecule_tab = QWidget()
        molecule_layout = QVBoxLayout()

        self.smiles_input = QLineEdit(self)
        molecule_layout.addWidget(QLabel('Enter SMILES String:'))
        molecule_layout.addWidget(self.smiles_input)

        render_button = QPushButton('Render Molecule', self)
        render_button.clicked.connect(self.render_molecule)
        molecule_layout.addWidget(render_button)

        self.mol_display = QGraphicsView(self)
        molecule_layout.addWidget(QLabel('Molecular Drawing:'))
        molecule_layout.addWidget(self.mol_display)

        molecule_tab.setLayout(molecule_layout)
        tab_widget.addTab(molecule_tab, "Molecule")

        # Tab for "About" information
        about_tab = QWidget()
        about_layout = QVBoxLayout()

        about_text = QLabel("This is a simple molecule drawer using PyQt5 and RDKit.")
        #about_text.setAlignment(QtCore.Qt.AlignCenter)
        about_layout.addWidget(about_text)

        about_tab.setLayout(about_layout)
        tab_widget.addTab(about_tab, "About")

        layout.addWidget(tab_widget)

        # Set up the main window
        self.setLayout(layout)
        self.setGeometry(100, 100, 600, 400)
        self.setWindowTitle('Molecule Drawer')

    def render_molecule(self):
        smiles = self.smiles_input.text()

        mol = Chem.MolFromSmiles(smiles)

        if mol is not None:
            img = Draw.MolToImage(mol)
            byte_array = BytesIO()
            img.save(byte_array, format='PNG')

            pixmap = QPixmap()
            pixmap.loadFromData(byte_array.getvalue())
            scene = QGraphicsScene()
            scene.addPixmap(pixmap)
            self.mol_display.setScene(scene)
        else:
            error_message = QLabel('Invalid SMILES string. Please enter a valid SMILES.')
            layout = self.layout()
            layout.addWidget(error_message)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MolDrawer()
    window.show()
    sys.exit(app.exec_())
