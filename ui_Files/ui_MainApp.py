# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\ColvineC\IOS_DFO\ARGO-Float-Software\ui_Files\MainApp.ui'
#
# Created: Tue Jul 19 13:47:03 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_MainApp(object):
    def setupUi(self, MainApp):
        MainApp.setObjectName("MainApp")
        MainApp.setEnabled(True)
        MainApp.resize(943, 806)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainApp.sizePolicy().hasHeightForWidth())
        MainApp.setSizePolicy(sizePolicy)
        MainApp.setAutoFillBackground(True)
        self.centralwidget = QtGui.QWidget(MainApp)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label = QtGui.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setUnderline(True)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.horizontalLayout_2.addWidget(self.label)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.settingsButton = QtGui.QPushButton(self.centralwidget)
        self.settingsButton.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.settingsButton.setFont(font)
        self.settingsButton.setObjectName("settingsButton")
        self.horizontalLayout_2.addWidget(self.settingsButton)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.listAllPages = QtGui.QStackedWidget(self.centralwidget)
        self.listAllPages.setFrameShape(QtGui.QFrame.StyledPanel)
        self.listAllPages.setFrameShadow(QtGui.QFrame.Plain)
        self.listAllPages.setObjectName("listAllPages")
        self.programListPage = QtGui.QWidget()
        self.programListPage.setObjectName("programListPage")
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.programListPage)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_2 = QtGui.QLabel(self.programListPage)
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setUnderline(True)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_3.addWidget(self.label_2)
        spacerItem1 = QtGui.QSpacerItem(600, 20, QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem1)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        spacerItem2 = QtGui.QSpacerItem(20, 26, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout_2.addItem(spacerItem2)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.timeSeriesButton = QtGui.QPushButton(self.programListPage)
        self.timeSeriesButton.setMinimumSize(QtCore.QSize(0, 0))
        self.timeSeriesButton.setMaximumSize(QtCore.QSize(99, 16777215))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.timeSeriesButton.setFont(font)
        self.timeSeriesButton.setObjectName("timeSeriesButton")
        self.horizontalLayout.addWidget(self.timeSeriesButton)
        self.label_19 = QtGui.QLabel(self.programListPage)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_19.setFont(font)
        self.label_19.setWordWrap(True)
        self.label_19.setObjectName("label_19")
        self.horizontalLayout.addWidget(self.label_19)
        spacerItem3 = QtGui.QSpacerItem(100, 20, QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem3)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        spacerItem4 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Maximum)
        self.verticalLayout_2.addItem(spacerItem4)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.circulationButton = QtGui.QPushButton(self.programListPage)
        self.circulationButton.setMinimumSize(QtCore.QSize(0, 0))
        self.circulationButton.setMaximumSize(QtCore.QSize(99, 16777215))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.circulationButton.setFont(font)
        self.circulationButton.setObjectName("circulationButton")
        self.horizontalLayout_5.addWidget(self.circulationButton)
        self.label_20 = QtGui.QLabel(self.programListPage)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_20.setFont(font)
        self.label_20.setWordWrap(True)
        self.label_20.setObjectName("label_20")
        self.horizontalLayout_5.addWidget(self.label_20)
        spacerItem5 = QtGui.QSpacerItem(100, 20, QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_5.addItem(spacerItem5)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        spacerItem6 = QtGui.QSpacerItem(20, 330, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.MinimumExpanding)
        self.verticalLayout_2.addItem(spacerItem6)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        spacerItem7 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem7)
        self.exitButton = QtGui.QPushButton(self.programListPage)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.exitButton.setFont(font)
        self.exitButton.setObjectName("exitButton")
        self.horizontalLayout_4.addWidget(self.exitButton)
        self.verticalLayout_2.addLayout(self.horizontalLayout_4)
        self.listAllPages.addWidget(self.programListPage)
        self.timeSeriesPage = QtGui.QWidget()
        self.timeSeriesPage.setObjectName("timeSeriesPage")
        self.verticalLayout_6 = QtGui.QVBoxLayout(self.timeSeriesPage)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.horizontalLayout_18 = QtGui.QHBoxLayout()
        self.horizontalLayout_18.setObjectName("horizontalLayout_18")
        self.label_4 = QtGui.QLabel(self.timeSeriesPage)
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setUnderline(True)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_18.addWidget(self.label_4)
        spacerItem8 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_18.addItem(spacerItem8)
        self.verticalLayout_6.addLayout(self.horizontalLayout_18)
        self.timeSeriesApp = TimeSeriesApp(self.timeSeriesPage)
        self.timeSeriesApp.setMinimumSize(QtCore.QSize(400, 600))
        self.timeSeriesApp.setObjectName("timeSeriesApp")
        self.verticalLayout_6.addWidget(self.timeSeriesApp)
        self.horizontalLayout_19 = QtGui.QHBoxLayout()
        self.horizontalLayout_19.setObjectName("horizontalLayout_19")
        self.backToProgramListButton = QtGui.QPushButton(self.timeSeriesPage)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.backToProgramListButton.setFont(font)
        self.backToProgramListButton.setObjectName("backToProgramListButton")
        self.horizontalLayout_19.addWidget(self.backToProgramListButton)
        spacerItem9 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_19.addItem(spacerItem9)
        self.verticalLayout_6.addLayout(self.horizontalLayout_19)
        self.listAllPages.addWidget(self.timeSeriesPage)
        self.circulationPage = QtGui.QWidget()
        self.circulationPage.setObjectName("circulationPage")
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.circulationPage)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.horizontalLayout_20 = QtGui.QHBoxLayout()
        self.horizontalLayout_20.setObjectName("horizontalLayout_20")
        self.label_5 = QtGui.QLabel(self.circulationPage)
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setUnderline(True)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_20.addWidget(self.label_5)
        spacerItem10 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_20.addItem(spacerItem10)
        self.verticalLayout_3.addLayout(self.horizontalLayout_20)
        self.circulationApp = CirculationApp(self.circulationPage)
        self.circulationApp.setMinimumSize(QtCore.QSize(0, 600))
        self.circulationApp.setObjectName("circulationApp")
        self.verticalLayout_3.addWidget(self.circulationApp)
        self.horizontalLayout_21 = QtGui.QHBoxLayout()
        self.horizontalLayout_21.setObjectName("horizontalLayout_21")
        self.backToProgramListButton2 = QtGui.QPushButton(self.circulationPage)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.backToProgramListButton2.setFont(font)
        self.backToProgramListButton2.setObjectName("backToProgramListButton2")
        self.horizontalLayout_21.addWidget(self.backToProgramListButton2)
        spacerItem11 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_21.addItem(spacerItem11)
        self.verticalLayout_3.addLayout(self.horizontalLayout_21)
        self.listAllPages.addWidget(self.circulationPage)
        self.verticalLayout.addWidget(self.listAllPages)
        MainApp.setCentralWidget(self.centralwidget)
        self.statusbar = QtGui.QStatusBar(MainApp)
        self.statusbar.setObjectName("statusbar")
        MainApp.setStatusBar(self.statusbar)

        self.retranslateUi(MainApp)
        self.listAllPages.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainApp)

    def retranslateUi(self, MainApp):
        MainApp.setWindowTitle(QtGui.QApplication.translate("MainApp", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainApp", "ARGO Programs", None, QtGui.QApplication.UnicodeUTF8))
        self.settingsButton.setText(QtGui.QApplication.translate("MainApp", "Settings", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainApp", "Program List", None, QtGui.QApplication.UnicodeUTF8))
        self.timeSeriesButton.setText(QtGui.QApplication.translate("MainApp", "Time Series", None, QtGui.QApplication.UnicodeUTF8))
        self.label_19.setText(QtGui.QApplication.translate("MainApp", "Time Series is a program that interpolates data at a specific location and time, based on the data available from ARGO floats within the selected nearby area", None, QtGui.QApplication.UnicodeUTF8))
        self.circulationButton.setText(QtGui.QApplication.translate("MainApp", "Circulation", None, QtGui.QApplication.UnicodeUTF8))
        self.label_20.setText(QtGui.QApplication.translate("MainApp", "This thing does something about water circulation somewhere on the coast of somewhere?", None, QtGui.QApplication.UnicodeUTF8))
        self.exitButton.setText(QtGui.QApplication.translate("MainApp", "Exit", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("MainApp", "Time Series", None, QtGui.QApplication.UnicodeUTF8))
        self.backToProgramListButton.setText(QtGui.QApplication.translate("MainApp", "Back to Program List", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("MainApp", "Circulation", None, QtGui.QApplication.UnicodeUTF8))
        self.backToProgramListButton2.setText(QtGui.QApplication.translate("MainApp", "Back to Program List", None, QtGui.QApplication.UnicodeUTF8))

from ui_Classes.CirculationApp import CirculationApp
from ui_Classes.TimeSeriesApp import TimeSeriesApp
