''' 
MainApp.py
Callum Colvine - Research Assistant
Callum.Colvine@dfo-mpo.gc.ca

MainApp is the widget allowing access to TimeSeriesApp and CirculationApp.

Following Pep 8 formatting with the following exceptions:
- There is no spacing between a docstring and a function

MainApp hosts the GUI window surrounding the TimeSeries and Circulation 
Applications.  
'''

# This is to deal with path issues for the sake of project organizaiton
import sys
import os
sys.path.append("..")
# For development, make this your path to the project ui_Classes folder
# sys.path.append('C:\Users\ColvineC\IOS_DFO\ARGO-Float-Software')
sys.path.append(os.path.abspath(os.curdir))
# ToDo: os.path.join(os.getcwd()) <--- Try this instead

from PySide import QtCore, QtGui
from TimeSeriesApp import TimeSeriesApp
from ui_Files.ui_MainApp import Ui_MainApp
from PySide.QtGui import QMainWindow


class MainApp(QMainWindow, Ui_MainApp):

    def __init__(self):
        # hasInit
        super(MainApp, self).__init__()
        self.setupUi(self)
        self.setupSignals()
        # self.timeSeriesApp.setupUi(self)
        return

    def setupSignals(self):
        self.timeSeriesButton.clicked.connect(self.timeSeriesButtonClicked)
        self.exitButton.clicked.connect(self.close)
        self.backToProgramListButton.clicked.connect(
            self.backToProgramListButtonClicked)
        self.backToProgramListButton2.clicked.connect(
            self.backToProgramListButtonClicked)
        self.circulationButton.clicked.connect(self.circulationButtonClicked)

        return

    def circulationButtonClicked(self):
        # sys.path.append(self.argoPathLineEdit.text())
        self.listAllPages.setCurrentWidget(self.circulationPage)
        self.circulationApp.experimentSelected()
        self.circulationApp.argoPath = self.argoPathLineEdit.text()
        self.circulationApp.path0 = self.argoDataLineEdit.text() + \
                "argo_mirror\pacific_ocean\\"
        self.circulationApp.outPath = self.circulationApp.argoPath + \
                "argo_out_TEST\\Circulation\\"
        self.circulationApp.readEVals()
        
        return

    def timeSeriesButtonClicked(self):
        # sys.path.append(self.argoPathLineEdit.text())
        self.listAllPages.setCurrentWidget(self.timeSeriesPage)
        self.timeSeriesApp.experimentSelected()
        self.timeSeriesApp.localSoftwarePath = self.argoPathLineEdit.text()
        self.timeSeriesApp.drive = self.argoDataLineEdit.text()
        self.timeSeriesApp.driveEdited()
        return

    def backToProgramListButtonClicked(self):
        self.listAllPages.setCurrentWidget(self.programListPage)
        return

def main():
    app = QtGui.QApplication(sys.argv)
    mySW = MainApp()
    mySW.setWindowTitle("ARGO Homepage")
    mySW.show()
    sys.exit(app.exec_())
    return

if __name__ == "__main__":
    main()
