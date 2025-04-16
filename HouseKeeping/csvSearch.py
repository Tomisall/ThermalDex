import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QPushButton, QFileDialog
import pandas as pd

class CSVViewer(QWidget):
    def __init__(self):
        super().__init__()

        self.init_ui()

    def init_ui(self):
        self.df = None
        self.current_index = 0

        self.layout = QVBoxLayout()

        self.result_label = QLabel()
        self.counter_label = QLabel()

        self.layout.addWidget(self.result_label)
        self.layout.addWidget(self.counter_label)

        self.prev_button = QPushButton('Previous')
        self.next_button = QPushButton('Next')

        self.prev_button.clicked.connect(self.prev_result)
        self.next_button.clicked.connect(self.next_result)

        self.layout.addWidget(self.prev_button)
        self.layout.addWidget(self.next_button)

        self.setLayout(self.layout)

        self.load_csv()

        self.show_result()

    def load_csv(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(self, "Open CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)

        if file_name:
            self.df = pd.read_csv(file_name)

    def show_result(self):
        if self.df is not None and not self.df.empty:
            current_row = self.df.iloc[self.current_index]
            result_text = f"SMILES: {current_row['SMILES']}\nName: {current_row['Name']}\nMW: {current_row['MW']}\nmp: {current_row['mp']}\nProject: {current_row['Project']}"
            self.result_label.setText(result_text)
            self.counter_label.setText(f"Result {self.current_index + 1} of {len(self.df)}")

    def prev_result(self):
        if self.current_index > 0:
            self.current_index -= 1
            self.show_result()

    def next_result(self):
        if self.df is not None and self.current_index < len(self.df) - 1:
            self.current_index += 1
            self.show_result()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    viewer = CSVViewer()
    viewer.setGeometry(100, 100, 800, 600)
    viewer.setWindowTitle('CSV Viewer')
    viewer.show()
    sys.exit(app.exec_())
