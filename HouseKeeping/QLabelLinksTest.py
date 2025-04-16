import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QVBoxLayout

def window():
   app = QApplication(sys.argv)
   win = QWidget() 
	
   l1 = QLabel()
   l2 = QLabel()
   l3 = QLabel()
   l4 = QLabel()
	
   l1.setText("Hello World")
   l4.setText("<a href=\"http://www.google.com\">'Click this link to go to Google'</a>")
   l2.setText("<a href=\"http://www.google.com\">'Click this link to go to Google'</a>")
	
   l1.setAlignment(Qt.AlignCenter)
   l3.setAlignment(Qt.AlignCenter)
   l4.setAlignment(Qt.AlignRight)
   l3.setPixmap(QPixmap("ThermalDexIcon.jpg"))
	
   vbox = QVBoxLayout()
   vbox.addWidget(l1)
   vbox.addStretch()
   vbox.addWidget(l2)
   vbox.addStretch()
   vbox.addWidget(l3)
   vbox.addStretch()
   vbox.addWidget(l4)
	
   l1.setOpenExternalLinks(True)
   l4.linkActivated.connect(clicked)
   l2.linkHovered.connect(hovered)
   l1.setTextInteractionFlags(Qt.TextSelectableByMouse)
   win.setLayout(vbox)
	
   win.setWindowTitle("QLabel Demo")
   win.show()
   sys.exit(app.exec_())
	
def hovered():
   print("hovering")
def clicked():
   print("clicked")
	
if __name__ == '__main__':
   window()