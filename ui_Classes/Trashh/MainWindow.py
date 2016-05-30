# This is to deal with path issues for the sake of project organizaiton
import sys
sys.path.append("..")

from PySide import QtCore, QtGui
from PySide.QtCore import QSize
from PySide.QtGui import QMainWindow
from ui_Files.ui_MainApp import Ui_MainApp

class MainApp(QMainWindow, Ui_MainApp):
    def __init__(self):
        super(MainApp, self).__init__()
        # self.ui = Ui_MainWindow()
        self.setupUi(self)


def main():
    app = QtGui.QApplication(sys.argv)
    mySW = MainApp()
    mySW.setWindowTitle("ARGO Homepage")
    mySW.show()
    sys.exit(app.exec_())
    return

if __name__ == "__main__":
    main()

