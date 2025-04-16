import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem, QPushButton, QLabel
from PyQt5.QtGui import QColor

class MyWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        self.setWindowTitle('Table Example')

        layout = QVBoxLayout()

        # Create a table with 1 row and 3 columns
        title_label = QLabel('Assessment of Hazard by Scale:')
        layout.addWidget(title_label)
        self.table = QTableWidget(1, 4)
        self.table.setHorizontalHeaderLabels(['<5 g', '5 to <100 g', '100 to 500 g', '>500 g'])
        self.table.verticalHeader().setVisible(False)
        #self.table.setStyleSheet("QTableWidget::item { border-bottom: 2px solid black; }")
        self.table.setMaximumHeight(52)
        layout.addWidget(self.table)


        # Create a button
        self.button = QPushButton('Fill Table')
        self.button.clicked.connect(self.fillTable)
        layout.addWidget(self.button)

        self.setLayout(layout)

    def fillTable(self):
        # Clear previous contents
        self.table.clearContents()

        # Add values to cells
        values = [10, 5, 20, 66]  # You can replace this with your own values
        for i, value in enumerate(values):
            item = QTableWidgetItem(str(value))
            self.table.setItem(0, i, item)

            # Color code cells based on values
            color = self.getColorForValue(value)
            item.setBackground(color)

    def getColorForValue(self, value):
        # Example color-coding logic
        if value < 10:
            return QColor(255, 0, 0)  # Red
        elif value < 20:
            return QColor(255, 255, 0)  # Yellow
        else:
            return QColor(0, 255, 0)  # Green


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MyWindow()
    window.setGeometry(100, 100, 400, 200)
    window.show()
    sys.exit(app.exec_())
