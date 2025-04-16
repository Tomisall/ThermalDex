import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QListWidget, QTabWidget

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
        list_search_results = QListWidget()
        layout.addWidget(list_search_results)

        def search_database():
            query = entry_search.text().lower()
            results = [entry for entry in self.database if query in entry['Name'].lower() or query in entry['Formula'].lower()]

            list_search_results.clear()
            for result in results:
                list_search_results.addItem(f"{result['Name']} ({result['Formula']})")

        # Button to search
        btn_search = QPushButton('Search', clicked=search_database)
        layout.addWidget(btn_search)

    def create_add_to_database_tab(self):
        tab = QWidget()
        self.tab_widget.addTab(tab, 'Add to Database')

        # Layout for the "Add to Database" tab
        layout = QVBoxLayout(tab)

        # Entry widgets for adding to the database
        lbl_name = QLabel('Chemical Name:')
        lbl_formula = QLabel('Chemical Formula:')
        entry_name = QLineEdit()
        entry_formula = QLineEdit()

        layout.addWidget(lbl_name)
        layout.addWidget(entry_name)
        layout.addWidget(lbl_formula)
        layout.addWidget(entry_formula)

        # Add to database function
        def add_to_database():
            name = entry_name.text()
            formula = entry_formula.text()

            if name and formula:
                self.database.append({'Name': name, 'Formula': formula})
                entry_name.clear()
                entry_formula.clear()

        # Button to add to database
        btn_add = QPushButton('Add to Database', clicked=add_to_database)
        layout.addWidget(btn_add)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    mainWin = ChemicalInventoryApp()
    mainWin.show()
    sys.exit(app.exec_())
