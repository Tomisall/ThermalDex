import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene
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

        # Input field for SMILES string
        self.smiles_input = QLineEdit(self)
        layout.addWidget(QLabel('Enter SMILES String:'))
        layout.addWidget(self.smiles_input)

        # Button to render the molecule
        render_button = QPushButton('Render Molecule', self)
        render_button.clicked.connect(self.render_molecule)
        layout.addWidget(render_button)

        # Display area for the molecular drawing
        self.mol_display = QGraphicsView(self)
        layout.addWidget(QLabel('Molecular Drawing:'))
        layout.addWidget(self.mol_display)

        # Set up the main window
        self.setLayout(layout)
        self.setGeometry(100, 100, 600, 400)
        self.setWindowTitle('Molecule Drawer')

    def render_molecule(self):
        # Get the SMILES string from the input field
        smiles = self.smiles_input.text()

        # Create an RDKit molecule from the SMILES string
        mol = Chem.MolFromSmiles(smiles)

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
        else:
            # Display an error message if the SMILES string is invalid
            error_message = QLabel('Invalid SMILES string. Please enter a valid SMILES.')
            layout = self.layout()
            layout.addWidget(error_message)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MolDrawer()
    window.show()
    sys.exit(app.exec_())
