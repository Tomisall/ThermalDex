import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, QPushButton, QVBoxLayout, QWidget, QMessageBox, QFileDialog
import configparser

class SettingsWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Settings")
        self.setGeometry(100, 100, 400, 200)

        self.config = configparser.ConfigParser()
        self.config.read('./_core/ThermalDex.ini')

        self.default_file_label = QLabel("Default CSV File:")
        self.default_file_input = QLineEdit(self.config.get('Database', 'defaultDB'))
        self.default_file_button = QPushButton("Browse")
        self.default_file_button.clicked.connect(self.select_default_file)

        self.save_button = QPushButton("Save")
        self.save_button.clicked.connect(self.save_settings)

        layout = QVBoxLayout()
        layout.addWidget(self.default_file_label)
        layout.addWidget(self.default_file_input)
        layout.addWidget(self.default_file_button)
        layout.addWidget(self.save_button)

        self.setLayout(layout)

    def select_default_file(self):
        options = QFileDialog.Options()
        default_file, _ = QFileDialog.getOpenFileName(self, "Select Database to Use:", "", "CSV Files (*.csv)", options=options)
        if default_file:
            self.default_file_input.setText(default_file)

    def save_settings(self):
        self.config.set('Database', 'defaultDB', self.default_file_input.text())

        with open('./_core/ThermalDex.ini', 'w') as configfile:
            self.config.write(configfile)

        QMessageBox.information(self, "Settings Saved", "Settings have been saved.")

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Main Window")
        self.setGeometry(200, 200, 400, 200)

        self.settings_window = SettingsWindow()

        self.settings_button = QPushButton("Open Settings")
        self.settings_button.clicked.connect(self.open_settings)

        layout = QVBoxLayout()
        layout.addWidget(self.settings_button)

        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

    def open_settings(self):
        self.settings_window.show()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
