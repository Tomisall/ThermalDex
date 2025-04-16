import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QGraphicsView, QGraphicsScene, QFrame
from rdkit import Chem
from rdkit.Chem import Draw
from PyQt5.QtGui import QPixmap
from io import BytesIO

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
        # Set up the main layout
        layout = QVBoxLayout()

        # Display area for the molecular drawing
        self.mol_display = QGraphicsView(self)
        layout.addWidget(QLabel('Molecule:'))
        layout.addWidget(self.mol_display)

	#Add labels for calculated values
        test = 'MW: '
        self.mwLabel = QLabel(test)
        layout.addWidget(self.mwLabel)
        layout.addWidget(QLabel('Number of High Energy Groups:'))

        layout.addWidget(QHLine())

        # Input field for SMILES string
        self.smiles_input = QLineEdit(self)
        layout.addWidget(QLabel('Enter SMILES String:'))
        layout.addWidget(self.smiles_input)

        # Input field for Name string
        self.name_input = QLineEdit(self)
        layout.addWidget(QLabel('Name:'))
        layout.addWidget(self.name_input)

        # Input field for mp string
        self.mp_input = QLineEdit(self)
        layout.addWidget(QLabel('m.p.:'))
        layout.addWidget(self.mp_input)

        # Input field for TE string
        self.TE_input = QLineEdit(self)
        layout.addWidget(QLabel('Thermal Event:'))
        layout.addWidget(self.TE_input)

        # Input field for proj string
        self.proj_input = QLineEdit(self)
        layout.addWidget(QLabel('Project:'))
        layout.addWidget(self.proj_input)

        # Button to render the molecule
        render_button = QPushButton('Render Molecule', self)
        render_button.clicked.connect(self.render_molecule)
        layout.addWidget(render_button)



        # Set up the main window
        self.setLayout(layout)
        self.setGeometry(100, 100, 400, 650)
        self.setWindowTitle('Molecule Drawer')

    def render_molecule(self):
        # Get the SMILES string from the input field
        smiles = self.smiles_input.text()

        # Create an RDKit molecule from the SMILES string
        mol = Chem.MolFromSmiles(smiles)
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
            self.mwLabel.setText(test) 

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
