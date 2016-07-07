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
        self.plotCentre = None
        self.rangeCentre = None
        self.firstLatitude = self.firstLatitudeBox.value()
        self.secondLatitude = self.secondLatitudeBox.value()
        self.firstLongitude = self.firstLongitudeBox.value()
        self.secondLongitude = self.secondLongitudeBox.value()
        self.pressureCutOff = self.pressureCutOffBox.value()
        self.simplifyList = self.simplifyListCheckBox.isChecked()
        self.maxInterpDepth = self.maxInterpDepthBox.value()
        self.pressureStepSize = self.stepSizeBox.value()
        self.dynHeightsAtP = self.dynHeightAtPBox.value()
        self.relativeToPref = self.relativeToPrefBox.value()
        self.totalModes = self.totalModesBox.value()
        self.generateArray = self.generateArrayCheckBox.isChecked()
        self.pauseOn20 = self.pauseOn20thCheckBox.isChecked()
        self.entryString = self.entryStringLineEdit.text()
        return

    def setupSignals(self):
        self.firstLatitudeBox.editingFinished.connect(
            self.firstLatitudeBoxEditingFinished)
        self.secondLatitudeBox.editingFinished.connect(
            self.secondLatitudeBoxEditingFinished)
        self.firstLongitudeBox.editingFinished.connect(
            self.firstLongitudeBoxEditingFinished)
        self.secondLongitudeBox.editingFinished.connect(
            self.secondLongitudeBoxEditingFinished)
        self.pressureCutOffBox.editingFinished.connect(
            self.pressureCutOffBoxEditingFinished)
        self.maxInterpDepthBox.editingFinished.connect(
            self.maxInterpDepthBoxEditingFinished)
        self.stepSizeBox.editingFinished.connect(
            self.stepSizeBoxEditingFinished)
        self.simplifyListCheckBox.stateChanged.connect(
            self.simplifyListCheckBoxStateChanged)
        self.dynHeightAtPBox.editingFinished.connect(
            self.dynHeightAtPBoxEditingFinished)
        self.relativeToPrefBox.editingFinished.connect(
            self.relativeToPrefBoxEditingFinished)
        self.totalModesBox.editingFinished.connect(
            self.totalModesBoxEditingFinished)
        self.generateArrayCheckBox.stateChanged.connect(
            self.generateArrayCheckBoxStateChanged)
        self.pauseOn20thCheckBox.stateChanged.connect(
            self.pauseOn20thCheckBoxStateChanged)
        self.entryStringLineEdit.editingFinished.connect(
            self.entryStringLineEditEditingFinished)
        return


    def firstLatitudeBoxEditingFinished(self):
        self.firstLatitude = self.firstLatitudeBox.value()
        return

    def secondLatitudeBoxEditingFinished(self):
        self.secondLatitude = self.secondLatitudeBox.value()
        return

    def firstLongitudeBoxEditingFinished(self):
        self.firstLongitude = self.firstLongitudeBox.value()
        return

    def secondLongitudeBoxEditingFinished(self):
        self.secondLongitude = self.secondLongitudeBox.value()
        return

    def pressureCutOffBoxEditingFinished(self):
        newP = self.pressureCutOffBox.value()
        if newP > 2000:
            style = "QLabel { color : red; }"
            self.pressureCutOffWarningLabel.setStyleSheet(style)
            out = "Warning: unlikely to find many usable floats past -2000m"
            self.pressureCutOffWarningLabel.setText(out)
        else:
            self.pressureCutOffWarningLabel.setText("")
        self.pressureCutOff = newP
        return

    def maxInterpDepthBoxEditingFinished(self):
        self.maxInterpDepth = self.maxInterpDepthBox.value()
        self.updateNPress()
        return

    def stepSizeBoxEditingFinished(self):
        self.stepSize = self.stepSizeBox.value()
        self.updateNPress()
        return

    def simplifyListCheckBoxStateChanged(self):
        self.simplifyList = self.simplifyListCheckBox.isChecked()
        return

    def relativeToPrefBoxEditingFinished(self):
        self.relativeToPref = self.relativeToPrefBox.value()
        return

    def totalModesBoxEditingFinished(self):
        self.totalModes = self.totalModesBox.value()
        return

    def generateArrayCheckBoxStateChanged(self):
        self.generateArray = self.generateArrayCheckBox.isChecked()
        return

    def pauseOn20thCheckBoxStateChanged(self):
        self.pauseOn20 = self.pauseOn20thCheckBox.isChecked()
        return

    def entryStringLineEditEditingFinished(self):
        self.entryString = self.entryStringLineEdit.text()
        return

    def dynHeightAtPBoxEditingFinished(self):
        self.dynHeightAtPBox = self.dynHeightAtPBox.value()
        return
