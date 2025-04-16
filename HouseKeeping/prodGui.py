import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QListWidget, QTabWidget, QFileDialog

class ChemicalInventoryApp(QMainWindow):
    def __init__(self):
        super(ChemicalInventoryApp, self).__init__()

        # Database
        self.database = []

        # Create tab widget
        self.tab_widget = QTabWidget()

        # Search Database Tab
        self.create_search_database_tab()

        # Add to Database Tab
        self.create_add_to_database_tab()

        # Settings Tab
        self.create_settings_tab()

        # About Tab
        self.create_about_tab()

        self.setCentralWidget(self.tab_widget)
        self.setWindowTitle("Chemical Inventory System")

    def create_search_database_tab(self):
        tab = QWidget()
        self.tab_widget.addTab(tab, 'Search Database')

        # Layout for the "Search Database" tab
        layout = QVBoxLayout(tab)

        # Entry widget for searching
        lbl_search = QLabel('Search:')
        entry_search = QLineEdit()
        layout.addWidget(lbl_search)
        layout.addWidget(entry_search)

        # List widget to display search results
        #list_search_results = QListWidget()
        #layout.addWidget(list_search_results)

        #def search_database():
        #    query = entry_search.text().lower()
        #    results = [entry for entry in self.database if query in entry['Name'].lower() or query in entry['Formula'].lower()]

        #    list_search_results.clear()
        #    for result in results:
        #        list_search_results.addItem(f"{result['Name']} ({result['Formula']})")

        # Button to search
        btn_search = QPushButton('Search') #, clicked=search_database)
        layout.addWidget(btn_search)

    def create_add_to_database_tab(self):
        tab = QWidget()
        self.tab_widget.addTab(tab, 'Add to Database')

        # Layout for the "Add to Database" tab
        layout = QHBoxLayout(tab)

        # Molecule View
        list_molecule_view = QListWidget()
        layout.addWidget(list_molecule_view)

        # Entry widgets for adding to the database
        entry_layout = QVBoxLayout()

        lbl_smiles = QLabel('SMILES:')
        lbl_name = QLabel('Name:')
        lbl_mw = QLabel('MW:')
        lbl_mp = QLabel('mp:')
        lbl_high_energy_groups = QLabel('Number of High Energy Groups:')
        lbl_thermal_event = QLabel('Thermal Event:')
        lbl_project = QLabel('Project:')

        entry_smiles = QLineEdit()
        entry_name = QLineEdit()
        entry_mw = QLineEdit()
        entry_mp = QLineEdit()
        entry_high_energy_groups = QLineEdit()
        entry_thermal_event = QLineEdit()
        entry_project = QLineEdit()

        entry_layout.addWidget(lbl_smiles)
        entry_layout.addWidget(entry_smiles)
        entry_layout.addWidget(lbl_name)
        entry_layout.addWidget(entry_name)
        entry_layout.addWidget(lbl_mw)
        entry_layout.addWidget(entry_mw)
        entry_layout.addWidget(lbl_mp)
        entry_layout.addWidget(entry_mp)
        entry_layout.addWidget(lbl_high_energy_groups)
        entry_layout.addWidget(entry_high_energy_groups)
        entry_layout.addWidget(lbl_thermal_event)
        entry_layout.addWidget(entry_thermal_event)
        entry_layout.addWidget(lbl_project)
        entry_layout.addWidget(entry_project)

        # Add to database function
        def add_to_database():
            smiles = entry_smiles.text()
            name = entry_name.text()
            mw = entry_mw.text()
            mp = entry_mp.text()
            high_energy_groups = entry_high_energy_groups.text()
            thermal_event = entry_thermal_event.text()
            project = entry_project.text()

            if name and smiles:
                self.database.append({
                    'SMILES': smiles,
                    'Name': name,
                    'MW': mw,
                    'mp': mp,
                    'Number of High Energy Groups': high_energy_groups,
                    'Thermal Event': thermal_event,
                    'Project': project
                })
                entry_smiles.clear()
                entry_name.clear()
                entry_mw.clear()
                entry_mp.clear()
                entry_high_energy_groups.clear()
                entry_thermal_event.clear()
                entry_project.clear()
                self.update_molecule_view(list_molecule_view)

        # Button to add to database
        btn_add = QPushButton('Add to Database', clicked=add_to_database)
        entry_layout.addWidget(btn_add)

        layout.addLayout(entry_layout)

    def update_molecule_view(self, list_molecule_view):
        # Update list widget in the "Add to Database" tab
        list_molecule_view.clear()
        for entry in self.database:
            list_molecule_view.addItem(f"{entry['Name']} ({entry['SMILES']})")

    def create_settings_tab(self):
        tab = QWidget()
        self.tab_widget.addTab(tab, 'Settings')

        # Layout for the "Settings" tab
        layout = QVBoxLayout(tab)

        # Button to select a database file
        btn_select_file = QPushButton('Select Database File', clicked=self.select_database_file)
        layout.addWidget(btn_select_file)

    def select_database_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly

        file_name, _ = QFileDialog.getOpenFileName(self, "Select Database File", "", "All Files (*);;Text Files (*.txt);;CSV Files (*.csv)", options=options)

        if file_name:
            # Load the database from the selected file (You can implement this part based on your file format)
            print(f"Selected Database File: {file_name}")

    def create_about_tab(self):
        tab = QWidget()
        self.tab_widget.addTab(tab, 'About')

        # Layout for the "About" tab
        layout = QVBoxLayout(tab)

        # Text label for information about the application
        about_text = "Chemical Inventory System\nVersion 1.0\n\nDeveloped by Your Name"
        lbl_about = QLabel(about_text)
        layout.addWidget(lbl_about)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainWin = ChemicalInventoryApp()
    mainWin.show()
    sys.exit(app.exec_())
