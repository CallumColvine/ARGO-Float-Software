# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\ColvineC\IOS_DFO\ARGO-Float-Software\ui_Files\circulationapp.ui'
#
# Created: Mon Jul 04 12:17:23 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_CirculationApp(object):
    def setupUi(self, CirculationApp):
        CirculationApp.setObjectName("CirculationApp")
        CirculationApp.resize(594, 416)
        self.verticalLayout = QtGui.QVBoxLayout(CirculationApp)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtGui.QLabel(CirculationApp)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)

        self.retranslateUi(CirculationApp)
        QtCore.QMetaObject.connectSlotsByName(CirculationApp)

    def retranslateUi(self, CirculationApp):
        CirculationApp.setWindowTitle(QtGui.QApplication.translate("CirculationApp", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("CirculationApp", "Under Development...", None, QtGui.QApplication.UnicodeUTF8))

