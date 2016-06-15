# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainApp.ui'
#
# Created: Tue Jun 14 14:40:04 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_MainApp(object):
    def setupUi(self, MainApp):
        MainApp.setObjectName("MainApp")
        MainApp.setEnabled(True)
        MainApp.resize(925, 765)
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
        font = QtGui.QFont()
        font.setPointSize(12)
        self.timeSeriesButton.setFont(font)
        self.timeSeriesButton.setObjectName("timeSeriesButton")
        self.horizontalLayout.addWidget(self.timeSeriesButton)
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem3)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        spacerItem4 = QtGui.QSpacerItem(20, 330, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.MinimumExpanding)
        self.verticalLayout_2.addItem(spacerItem4)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        spacerItem5 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem5)
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
        spacerItem6 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_18.addItem(spacerItem6)
        self.verticalLayout_6.addLayout(self.horizontalLayout_18)
        self.timeSeriesStackedWidget = QtGui.QStackedWidget(self.timeSeriesPage)
        self.timeSeriesStackedWidget.setFrameShape(QtGui.QFrame.StyledPanel)
        self.timeSeriesStackedWidget.setObjectName("timeSeriesStackedWidget")
        self.settingsPage = QtGui.QWidget()
        self.settingsPage.setObjectName("settingsPage")
        self.verticalLayout_5 = QtGui.QVBoxLayout(self.settingsPage)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.verticalLayout_5.addLayout(self.horizontalLayout_6)
        self.defaultSettingsCheckBox = QtGui.QCheckBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(75)
        font.setBold(True)
        self.defaultSettingsCheckBox.setFont(font)
        self.defaultSettingsCheckBox.setObjectName("defaultSettingsCheckBox")
        self.verticalLayout_5.addWidget(self.defaultSettingsCheckBox)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.floatRadiusBox = QtGui.QSpinBox(self.settingsPage)
        self.floatRadiusBox.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.floatRadiusBox.setFont(font)
        self.floatRadiusBox.setMaximum(9999)
        self.floatRadiusBox.setProperty("value", 50)
        self.floatRadiusBox.setObjectName("floatRadiusBox")
        self.horizontalLayout_5.addWidget(self.floatRadiusBox)
        self.label_3 = QtGui.QLabel(self.settingsPage)
        self.label_3.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_5.addWidget(self.label_3)
        spacerItem7 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_5.addItem(spacerItem7)
        self.verticalLayout_5.addLayout(self.horizontalLayout_5)
        self.label_5 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_5.addWidget(self.label_5)
        self.horizontalLayout_16 = QtGui.QHBoxLayout()
        self.horizontalLayout_16.setObjectName("horizontalLayout_16")
        spacerItem8 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_16.addItem(spacerItem8)
        self.startRangeDateEdit = QtGui.QDateEdit(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.startRangeDateEdit.setFont(font)
        self.startRangeDateEdit.setCalendarPopup(True)
        self.startRangeDateEdit.setDate(QtCore.QDate(2016, 5, 18))
        self.startRangeDateEdit.setObjectName("startRangeDateEdit")
        self.horizontalLayout_16.addWidget(self.startRangeDateEdit)
        self.label_15 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_15.setFont(font)
        self.label_15.setObjectName("label_15")
        self.horizontalLayout_16.addWidget(self.label_15)
        self.endRangeDateEdit = QtGui.QDateEdit(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.endRangeDateEdit.setFont(font)
        self.endRangeDateEdit.setCalendarPopup(True)
        self.endRangeDateEdit.setDate(QtCore.QDate(2016, 5, 23))
        self.endRangeDateEdit.setObjectName("endRangeDateEdit")
        self.horizontalLayout_16.addWidget(self.endRangeDateEdit)
        spacerItem9 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_16.addItem(spacerItem9)
        self.verticalLayout_5.addLayout(self.horizontalLayout_16)
        self.horizontalLayout_17 = QtGui.QHBoxLayout()
        self.horizontalLayout_17.setObjectName("horizontalLayout_17")
        spacerItem10 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_17.addItem(spacerItem10)
        self.currentDayLabel = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(50)
        font.setBold(False)
        self.currentDayLabel.setFont(font)
        self.currentDayLabel.setObjectName("currentDayLabel")
        self.horizontalLayout_17.addWidget(self.currentDayLabel)
        self.currentDayDateEdit = QtGui.QDateEdit(self.settingsPage)
        self.currentDayDateEdit.setEnabled(True)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.currentDayDateEdit.setFont(font)
        self.currentDayDateEdit.setReadOnly(True)
        self.currentDayDateEdit.setButtonSymbols(QtGui.QAbstractSpinBox.NoButtons)
        self.currentDayDateEdit.setObjectName("currentDayDateEdit")
        self.horizontalLayout_17.addWidget(self.currentDayDateEdit)
        spacerItem11 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_17.addItem(spacerItem11)
        self.currentDayRangeLabel = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.currentDayRangeLabel.setFont(font)
        self.currentDayRangeLabel.setObjectName("currentDayRangeLabel")
        self.horizontalLayout_17.addWidget(self.currentDayRangeLabel)
        self.latestDayDateEdit = QtGui.QDateEdit(self.settingsPage)
        self.latestDayDateEdit.setEnabled(True)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.latestDayDateEdit.setFont(font)
        self.latestDayDateEdit.setReadOnly(True)
        self.latestDayDateEdit.setButtonSymbols(QtGui.QAbstractSpinBox.NoButtons)
        self.latestDayDateEdit.setDate(QtCore.QDate(2020, 1, 1))
        self.latestDayDateEdit.setObjectName("latestDayDateEdit")
        self.horizontalLayout_17.addWidget(self.latestDayDateEdit)
        spacerItem12 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_17.addItem(spacerItem12)
        self.currentDayRangeLabel_2 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.currentDayRangeLabel_2.setFont(font)
        self.currentDayRangeLabel_2.setObjectName("currentDayRangeLabel_2")
        self.horizontalLayout_17.addWidget(self.currentDayRangeLabel_2)
        self.earliestDayDateEdit = QtGui.QDateEdit(self.settingsPage)
        self.earliestDayDateEdit.setEnabled(True)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.earliestDayDateEdit.setFont(font)
        self.earliestDayDateEdit.setReadOnly(True)
        self.earliestDayDateEdit.setButtonSymbols(QtGui.QAbstractSpinBox.NoButtons)
        self.earliestDayDateEdit.setDate(QtCore.QDate(2001, 9, 7))
        self.earliestDayDateEdit.setObjectName("earliestDayDateEdit")
        self.horizontalLayout_17.addWidget(self.earliestDayDateEdit)
        spacerItem13 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_17.addItem(spacerItem13)
        self.verticalLayout_5.addLayout(self.horizontalLayout_17)
        self.horizontalLayout_8 = QtGui.QHBoxLayout()
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        self.dayStepSizeBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.dayStepSizeBox.setFont(font)
        self.dayStepSizeBox.setMaximum(9999)
        self.dayStepSizeBox.setProperty("value", 5)
        self.dayStepSizeBox.setObjectName("dayStepSizeBox")
        self.horizontalLayout_8.addWidget(self.dayStepSizeBox)
        self.label_6 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.horizontalLayout_8.addWidget(self.label_6)
        spacerItem14 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem14)
        self.verticalLayout_5.addLayout(self.horizontalLayout_8)
        self.horizontalLayout_9 = QtGui.QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.sampleWindowBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.sampleWindowBox.setFont(font)
        self.sampleWindowBox.setMaximum(9999)
        self.sampleWindowBox.setProperty("value", 10)
        self.sampleWindowBox.setObjectName("sampleWindowBox")
        self.horizontalLayout_9.addWidget(self.sampleWindowBox)
        self.label_7 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.horizontalLayout_9.addWidget(self.label_7)
        spacerItem15 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_9.addItem(spacerItem15)
        self.verticalLayout_5.addLayout(self.horizontalLayout_9)
        self.horizontalLayout_10 = QtGui.QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.firstLatitudeBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.firstLatitudeBox.setFont(font)
        self.firstLatitudeBox.setMaximum(9999)
        self.firstLatitudeBox.setProperty("value", 40)
        self.firstLatitudeBox.setObjectName("firstLatitudeBox")
        self.horizontalLayout_10.addWidget(self.firstLatitudeBox)
        self.label_8 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.horizontalLayout_10.addWidget(self.label_8)
        spacerItem16 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem16)
        self.secondLatitudeBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.secondLatitudeBox.setFont(font)
        self.secondLatitudeBox.setMaximum(9999)
        self.secondLatitudeBox.setProperty("value", 65)
        self.secondLatitudeBox.setObjectName("secondLatitudeBox")
        self.horizontalLayout_10.addWidget(self.secondLatitudeBox)
        self.label_16 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_16.setFont(font)
        self.label_16.setObjectName("label_16")
        self.horizontalLayout_10.addWidget(self.label_16)
        spacerItem17 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem17)
        self.verticalLayout_5.addLayout(self.horizontalLayout_10)
        self.horizontalLayout_11 = QtGui.QHBoxLayout()
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        self.firstLongitudeBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.firstLongitudeBox.setFont(font)
        self.firstLongitudeBox.setMaximum(9999)
        self.firstLongitudeBox.setProperty("value", 200)
        self.firstLongitudeBox.setObjectName("firstLongitudeBox")
        self.horizontalLayout_11.addWidget(self.firstLongitudeBox)
        self.label_9 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.horizontalLayout_11.addWidget(self.label_9)
        spacerItem18 = QtGui.QSpacerItem(28, 20, QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_11.addItem(spacerItem18)
        self.secondLongitudeBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.secondLongitudeBox.setFont(font)
        self.secondLongitudeBox.setMaximum(9999)
        self.secondLongitudeBox.setProperty("value", 230)
        self.secondLongitudeBox.setObjectName("secondLongitudeBox")
        self.horizontalLayout_11.addWidget(self.secondLongitudeBox)
        self.label_17 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_17.setFont(font)
        self.label_17.setObjectName("label_17")
        self.horizontalLayout_11.addWidget(self.label_17)
        spacerItem19 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_11.addItem(spacerItem19)
        self.verticalLayout_5.addLayout(self.horizontalLayout_11)
        self.horizontalLayout_12 = QtGui.QHBoxLayout()
        self.horizontalLayout_12.setObjectName("horizontalLayout_12")
        self.pressureCutOffBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.pressureCutOffBox.setFont(font)
        self.pressureCutOffBox.setMaximum(99999)
        self.pressureCutOffBox.setProperty("value", 970)
        self.pressureCutOffBox.setObjectName("pressureCutOffBox")
        self.horizontalLayout_12.addWidget(self.pressureCutOffBox)
        self.label_10 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.horizontalLayout_12.addWidget(self.label_10)
        spacerItem20 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_12.addItem(spacerItem20)
        self.verticalLayout_5.addLayout(self.horizontalLayout_12)
        self.horizontalLayout_13 = QtGui.QHBoxLayout()
        self.horizontalLayout_13.setObjectName("horizontalLayout_13")
        self.maxInterpDepthBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.maxInterpDepthBox.setFont(font)
        self.maxInterpDepthBox.setMaximum(99999)
        self.maxInterpDepthBox.setProperty("value", 1000)
        self.maxInterpDepthBox.setObjectName("maxInterpDepthBox")
        self.horizontalLayout_13.addWidget(self.maxInterpDepthBox)
        self.label_11 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.horizontalLayout_13.addWidget(self.label_11)
        spacerItem21 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_13.addItem(spacerItem21)
        self.stepSizeBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.stepSizeBox.setFont(font)
        self.stepSizeBox.setMaximum(9999)
        self.stepSizeBox.setProperty("value", 5)
        self.stepSizeBox.setObjectName("stepSizeBox")
        self.horizontalLayout_13.addWidget(self.stepSizeBox)
        self.label_12 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.horizontalLayout_13.addWidget(self.label_12)
        spacerItem22 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_13.addItem(spacerItem22)
        self.verticalLayout_5.addLayout(self.horizontalLayout_13)
        spacerItem23 = QtGui.QSpacerItem(20, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_5.addItem(spacerItem23)
        self.label_13 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setWeight(75)
        font.setBold(True)
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.verticalLayout_5.addWidget(self.label_13)
        self.horizontalLayout_15 = QtGui.QHBoxLayout()
        self.horizontalLayout_15.setObjectName("horizontalLayout_15")
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.tempCheckBox = QtGui.QCheckBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.tempCheckBox.setFont(font)
        self.tempCheckBox.setChecked(True)
        self.tempCheckBox.setObjectName("tempCheckBox")
        self.verticalLayout_3.addWidget(self.tempCheckBox)
        self.salinityCheckBox = QtGui.QCheckBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.salinityCheckBox.setFont(font)
        self.salinityCheckBox.setChecked(True)
        self.salinityCheckBox.setObjectName("salinityCheckBox")
        self.verticalLayout_3.addWidget(self.salinityCheckBox)
        self.sigmaTCheckBox = QtGui.QCheckBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.sigmaTCheckBox.setFont(font)
        self.sigmaTCheckBox.setChecked(True)
        self.sigmaTCheckBox.setObjectName("sigmaTCheckBox")
        self.verticalLayout_3.addWidget(self.sigmaTCheckBox)
        self.spicinessCheckBox = QtGui.QCheckBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.spicinessCheckBox.setFont(font)
        self.spicinessCheckBox.setChecked(True)
        self.spicinessCheckBox.setObjectName("spicinessCheckBox")
        self.verticalLayout_3.addWidget(self.spicinessCheckBox)
        self.dynamicHeightCheckBox = QtGui.QCheckBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.dynamicHeightCheckBox.setFont(font)
        self.dynamicHeightCheckBox.setChecked(True)
        self.dynamicHeightCheckBox.setObjectName("dynamicHeightCheckBox")
        self.verticalLayout_3.addWidget(self.dynamicHeightCheckBox)
        self.horizontalLayout_15.addLayout(self.verticalLayout_3)
        spacerItem24 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_15.addItem(spacerItem24)
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.horizontalLayout_14 = QtGui.QHBoxLayout()
        self.horizontalLayout_14.setObjectName("horizontalLayout_14")
        self.latitudeDesiredBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.latitudeDesiredBox.setFont(font)
        self.latitudeDesiredBox.setProperty("value", 50)
        self.latitudeDesiredBox.setObjectName("latitudeDesiredBox")
        self.horizontalLayout_14.addWidget(self.latitudeDesiredBox)
        self.longitudeDesiredBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.longitudeDesiredBox.setFont(font)
        self.longitudeDesiredBox.setMaximum(300)
        self.longitudeDesiredBox.setProperty("value", 215)
        self.longitudeDesiredBox.setObjectName("longitudeDesiredBox")
        self.horizontalLayout_14.addWidget(self.longitudeDesiredBox)
        self.label_14 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_14.setFont(font)
        self.label_14.setObjectName("label_14")
        self.horizontalLayout_14.addWidget(self.label_14)
        spacerItem25 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_14.addItem(spacerItem25)
        self.verticalLayout_4.addLayout(self.horizontalLayout_14)
        self.appendCheckBox = QtGui.QCheckBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.appendCheckBox.setFont(font)
        self.appendCheckBox.setChecked(True)
        self.appendCheckBox.setObjectName("appendCheckBox")
        self.verticalLayout_4.addWidget(self.appendCheckBox)
        self.verboseCheckBox = QtGui.QCheckBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.verboseCheckBox.setFont(font)
        self.verboseCheckBox.setObjectName("verboseCheckBox")
        self.verticalLayout_4.addWidget(self.verboseCheckBox)
        spacerItem26 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_4.addItem(spacerItem26)
        self.horizontalLayout_15.addLayout(self.verticalLayout_4)
        spacerItem27 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_15.addItem(spacerItem27)
        self.verticalLayout_5.addLayout(self.horizontalLayout_15)
        spacerItem28 = QtGui.QSpacerItem(20, 50, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_5.addItem(spacerItem28)
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        spacerItem29 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_7.addItem(spacerItem29)
        self.backButton_2 = QtGui.QPushButton(self.settingsPage)
        self.backButton_2.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.backButton_2.setFont(font)
        self.backButton_2.setObjectName("backButton_2")
        self.horizontalLayout_7.addWidget(self.backButton_2)
        self.nextButton = QtGui.QPushButton(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.nextButton.setFont(font)
        self.nextButton.setObjectName("nextButton")
        self.horizontalLayout_7.addWidget(self.nextButton)
        self.verticalLayout_5.addLayout(self.horizontalLayout_7)
        self.timeSeriesStackedWidget.addWidget(self.settingsPage)
        self.calculationsPage = QtGui.QWidget()
        self.calculationsPage.setObjectName("calculationsPage")
        self.verticalLayout_7 = QtGui.QVBoxLayout(self.calculationsPage)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.horizontalLayout_21 = QtGui.QHBoxLayout()
        self.horizontalLayout_21.setObjectName("horizontalLayout_21")
        self.outputFilesLabel = QtGui.QLabel(self.calculationsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.outputFilesLabel.setFont(font)
        self.outputFilesLabel.setObjectName("outputFilesLabel")
        self.horizontalLayout_21.addWidget(self.outputFilesLabel)
        spacerItem30 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_21.addItem(spacerItem30)
        self.label_19 = QtGui.QLabel(self.calculationsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_19.setFont(font)
        self.label_19.setObjectName("label_19")
        self.horizontalLayout_21.addWidget(self.label_19)
        self.verticalLayout_7.addLayout(self.horizontalLayout_21)
        self.horizontalLayout_20 = QtGui.QHBoxLayout()
        self.horizontalLayout_20.setObjectName("horizontalLayout_20")
        self.label_21 = QtGui.QLabel(self.calculationsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_21.setFont(font)
        self.label_21.setObjectName("label_21")
        self.horizontalLayout_20.addWidget(self.label_21)
        spacerItem31 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_20.addItem(spacerItem31)
        self.label_20 = QtGui.QLabel(self.calculationsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_20.setFont(font)
        self.label_20.setObjectName("label_20")
        self.horizontalLayout_20.addWidget(self.label_20)
        self.verticalLayout_7.addLayout(self.horizontalLayout_20)
        self.horizontalLayout_23 = QtGui.QHBoxLayout()
        self.horizontalLayout_23.setObjectName("horizontalLayout_23")
        spacerItem32 = QtGui.QSpacerItem(200, 20, QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_23.addItem(spacerItem32)
        self.verticalLayout_7.addLayout(self.horizontalLayout_23)
        self.horizontalLayout_22 = QtGui.QHBoxLayout()
        self.horizontalLayout_22.setObjectName("horizontalLayout_22")
        spacerItem33 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_22.addItem(spacerItem33)
        self.backToSettingsButton = QtGui.QPushButton(self.calculationsPage)
        self.backToSettingsButton.setEnabled(True)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.backToSettingsButton.setFont(font)
        self.backToSettingsButton.setObjectName("backToSettingsButton")
        self.horizontalLayout_22.addWidget(self.backToSettingsButton)
        self.finishedButton = QtGui.QPushButton(self.calculationsPage)
        self.finishedButton.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.finishedButton.setFont(font)
        self.finishedButton.setObjectName("finishedButton")
        self.horizontalLayout_22.addWidget(self.finishedButton)
        self.verticalLayout_7.addLayout(self.horizontalLayout_22)
        self.timeSeriesStackedWidget.addWidget(self.calculationsPage)
        self.verticalLayout_6.addWidget(self.timeSeriesStackedWidget)
        self.horizontalLayout_19 = QtGui.QHBoxLayout()
        self.horizontalLayout_19.setObjectName("horizontalLayout_19")
        self.backToProgramListButton = QtGui.QPushButton(self.timeSeriesPage)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.backToProgramListButton.setFont(font)
        self.backToProgramListButton.setObjectName("backToProgramListButton")
        self.horizontalLayout_19.addWidget(self.backToProgramListButton)
        spacerItem34 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_19.addItem(spacerItem34)
        self.verticalLayout_6.addLayout(self.horizontalLayout_19)
        self.listAllPages.addWidget(self.timeSeriesPage)
        self.verticalLayout.addWidget(self.listAllPages)
        MainApp.setCentralWidget(self.centralwidget)
        self.statusbar = QtGui.QStatusBar(MainApp)
        self.statusbar.setObjectName("statusbar")
        MainApp.setStatusBar(self.statusbar)

        self.retranslateUi(MainApp)
        self.listAllPages.setCurrentIndex(1)
        self.timeSeriesStackedWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainApp)

    def retranslateUi(self, MainApp):
        MainApp.setWindowTitle(QtGui.QApplication.translate("MainApp", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainApp", "ARGO Programs", None, QtGui.QApplication.UnicodeUTF8))
        self.settingsButton.setText(QtGui.QApplication.translate("MainApp", "Settings", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainApp", "Program List", None, QtGui.QApplication.UnicodeUTF8))
        self.timeSeriesButton.setText(QtGui.QApplication.translate("MainApp", "Time Series", None, QtGui.QApplication.UnicodeUTF8))
        self.exitButton.setText(QtGui.QApplication.translate("MainApp", "Exit", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("MainApp", "Time Series", None, QtGui.QApplication.UnicodeUTF8))
        self.defaultSettingsCheckBox.setText(QtGui.QApplication.translate("MainApp", "Default Settings", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainApp", "List floats within radius (Km)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("MainApp", "Day Range (Day/Month/Year)", None, QtGui.QApplication.UnicodeUTF8))
        self.startRangeDateEdit.setDisplayFormat(QtGui.QApplication.translate("MainApp", "d/M/yyyy", None, QtGui.QApplication.UnicodeUTF8))
        self.label_15.setText(QtGui.QApplication.translate("MainApp", "to", None, QtGui.QApplication.UnicodeUTF8))
        self.endRangeDateEdit.setDisplayFormat(QtGui.QApplication.translate("MainApp", "d/M/yyyy", None, QtGui.QApplication.UnicodeUTF8))
        self.currentDayLabel.setText(QtGui.QApplication.translate("MainApp", "Current Day:", None, QtGui.QApplication.UnicodeUTF8))
        self.currentDayDateEdit.setDisplayFormat(QtGui.QApplication.translate("MainApp", "d/M/yyyy", None, QtGui.QApplication.UnicodeUTF8))
        self.currentDayRangeLabel.setText(QtGui.QApplication.translate("MainApp", "Latest Recommended Day: ", None, QtGui.QApplication.UnicodeUTF8))
        self.latestDayDateEdit.setDisplayFormat(QtGui.QApplication.translate("MainApp", "d/M/yyyy", None, QtGui.QApplication.UnicodeUTF8))
        self.currentDayRangeLabel_2.setText(QtGui.QApplication.translate("MainApp", "Earliest Recommended Day: ", None, QtGui.QApplication.UnicodeUTF8))
        self.earliestDayDateEdit.setDisplayFormat(QtGui.QApplication.translate("MainApp", "d/M/yyyy", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("MainApp", "Day Step Size", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("MainApp", "Sample Window", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("MainApp", "First Latitude Point", None, QtGui.QApplication.UnicodeUTF8))
        self.label_16.setText(QtGui.QApplication.translate("MainApp", "Second Latitude Point", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setText(QtGui.QApplication.translate("MainApp", "First Longitude Point", None, QtGui.QApplication.UnicodeUTF8))
        self.label_17.setText(QtGui.QApplication.translate("MainApp", "Second Longitude Point", None, QtGui.QApplication.UnicodeUTF8))
        self.label_10.setText(QtGui.QApplication.translate("MainApp", "Pressure Cut-Off (m)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setText(QtGui.QApplication.translate("MainApp", "Interpolate Vertically to P-max  (m)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_12.setText(QtGui.QApplication.translate("MainApp", "Step size for depth/pressure (m)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_13.setText(QtGui.QApplication.translate("MainApp", "Desired outputs:", None, QtGui.QApplication.UnicodeUTF8))
        self.tempCheckBox.setText(QtGui.QApplication.translate("MainApp", "Temperature", None, QtGui.QApplication.UnicodeUTF8))
        self.salinityCheckBox.setText(QtGui.QApplication.translate("MainApp", "Salinity", None, QtGui.QApplication.UnicodeUTF8))
        self.sigmaTCheckBox.setText(QtGui.QApplication.translate("MainApp", "Sigma-T", None, QtGui.QApplication.UnicodeUTF8))
        self.spicinessCheckBox.setText(QtGui.QApplication.translate("MainApp", "Spiciness", None, QtGui.QApplication.UnicodeUTF8))
        self.dynamicHeightCheckBox.setText(QtGui.QApplication.translate("MainApp", "Dynamic Height", None, QtGui.QApplication.UnicodeUTF8))
        self.label_14.setText(QtGui.QApplication.translate("MainApp", "Latitude and Longtitude of desired station", None, QtGui.QApplication.UnicodeUTF8))
        self.appendCheckBox.setText(QtGui.QApplication.translate("MainApp", "Append to existing files (pick yes)", None, QtGui.QApplication.UnicodeUTF8))
        self.verboseCheckBox.setText(QtGui.QApplication.translate("MainApp", "Verbose", None, QtGui.QApplication.UnicodeUTF8))
        self.backButton_2.setText(QtGui.QApplication.translate("MainApp", "Back", None, QtGui.QApplication.UnicodeUTF8))
        self.nextButton.setText(QtGui.QApplication.translate("MainApp", "Next", None, QtGui.QApplication.UnicodeUTF8))
        self.outputFilesLabel.setText(QtGui.QApplication.translate("MainApp", "Output Files:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_19.setText(QtGui.QApplication.translate("MainApp", "Output files location", None, QtGui.QApplication.UnicodeUTF8))
        self.label_21.setText(QtGui.QApplication.translate("MainApp", "Write to file as settings for next time", None, QtGui.QApplication.UnicodeUTF8))
        self.label_20.setText(QtGui.QApplication.translate("MainApp", "Station to Interpolate location", None, QtGui.QApplication.UnicodeUTF8))
        self.backToSettingsButton.setText(QtGui.QApplication.translate("MainApp", "Back", None, QtGui.QApplication.UnicodeUTF8))
        self.finishedButton.setText(QtGui.QApplication.translate("MainApp", "Finished", None, QtGui.QApplication.UnicodeUTF8))
        self.backToProgramListButton.setText(QtGui.QApplication.translate("MainApp", "Back to Program List", None, QtGui.QApplication.UnicodeUTF8))

