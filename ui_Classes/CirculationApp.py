# UI Imports
from PySide import QtCore, QtGui
from PySide.QtGui import QWidget
# Calculations imports
import numpy as np

from ui_Files.ui_circulationapp import Ui_CirculationApp



class CirculationApp(QWidget, Ui_CirculationApp):
    def __init__(self, parent):
        super(CirculationApp, self).__init__()
        self.setupUi(self)
        self.initAllClassVariables()
        self.setupSignals()
        return

    def initAllClassVariables(self):

        return

    def setupSignals(self):
        
        return