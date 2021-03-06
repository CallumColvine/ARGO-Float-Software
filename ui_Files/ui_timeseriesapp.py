# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\rosst\Documents\Argo\Software\ui_Files\timeseriesapp.ui'
#
# Created: Fri Aug 26 11:04:42 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_TimeSeriesApp(object):
    def setupUi(self, TimeSeriesApp):
        TimeSeriesApp.setObjectName("TimeSeriesApp")
        TimeSeriesApp.resize(818, 595)
        self.horizontalLayout = QtGui.QHBoxLayout(TimeSeriesApp)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.timeSeriesStackedWidget = QtGui.QStackedWidget(TimeSeriesApp)
        self.timeSeriesStackedWidget.setFrameShape(QtGui.QFrame.StyledPanel)
        self.timeSeriesStackedWidget.setObjectName("timeSeriesStackedWidget")
        self.settingsPage = QtGui.QWidget()
        self.settingsPage.setObjectName("settingsPage")
        self.verticalLayout_5 = QtGui.QVBoxLayout(self.settingsPage)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.verticalLayout_5.addLayout(self.horizontalLayout_6)
        self.label_5 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_5.addWidget(self.label_5)
        self.horizontalLayout_16 = QtGui.QHBoxLayout()
        self.horizontalLayout_16.setObjectName("horizontalLayout_16")
        spacerItem = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_16.addItem(spacerItem)
        self.startRangeDateEdit = QtGui.QDateEdit(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.startRangeDateEdit.setFont(font)
        self.startRangeDateEdit.setCalendarPopup(True)
        self.startRangeDateEdit.setDate(QtCore.QDate(2016, 1, 2))
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
        self.endRangeDateEdit.setDate(QtCore.QDate(2016, 1, 22))
        self.endRangeDateEdit.setObjectName("endRangeDateEdit")
        self.horizontalLayout_16.addWidget(self.endRangeDateEdit)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_16.addItem(spacerItem1)
        self.verticalLayout_5.addLayout(self.horizontalLayout_16)
        self.horizontalLayout_17 = QtGui.QHBoxLayout()
        self.horizontalLayout_17.setObjectName("horizontalLayout_17")
        spacerItem2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_17.addItem(spacerItem2)
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
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_17.addItem(spacerItem3)
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
        spacerItem4 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_17.addItem(spacerItem4)
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
        spacerItem5 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem5)
        self.verticalLayout_5.addLayout(self.horizontalLayout_8)
        self.horizontalLayout_9 = QtGui.QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.sampleWindowBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.sampleWindowBox.setFont(font)
        self.sampleWindowBox.setMinimum(1)
        self.sampleWindowBox.setMaximum(30)
        self.sampleWindowBox.setProperty("value", 10)
        self.sampleWindowBox.setObjectName("sampleWindowBox")
        self.horizontalLayout_9.addWidget(self.sampleWindowBox)
        self.label_7 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.horizontalLayout_9.addWidget(self.label_7)
        spacerItem6 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_9.addItem(spacerItem6)
        self.verticalLayout_5.addLayout(self.horizontalLayout_9)
        self.horizontalLayout_10 = QtGui.QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.firstLatitudeBox = QtGui.QDoubleSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.firstLatitudeBox.setFont(font)
        self.firstLatitudeBox.setDecimals(1)
        self.firstLatitudeBox.setMinimum(-90.0)
        self.firstLatitudeBox.setMaximum(90.0)
        self.firstLatitudeBox.setSingleStep(0.1)
        self.firstLatitudeBox.setProperty("value", 40.0)
        self.firstLatitudeBox.setObjectName("firstLatitudeBox")
        self.horizontalLayout_10.addWidget(self.firstLatitudeBox)
        self.label_8 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.horizontalLayout_10.addWidget(self.label_8)
        spacerItem7 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem7)
        self.secondLatitudeBox = QtGui.QDoubleSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.secondLatitudeBox.setFont(font)
        self.secondLatitudeBox.setDecimals(1)
        self.secondLatitudeBox.setMinimum(-90.0)
        self.secondLatitudeBox.setMaximum(90.0)
        self.secondLatitudeBox.setSingleStep(0.1)
        self.secondLatitudeBox.setProperty("value", 60.0)
        self.secondLatitudeBox.setObjectName("secondLatitudeBox")
        self.horizontalLayout_10.addWidget(self.secondLatitudeBox)
        self.label_16 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_16.setFont(font)
        self.label_16.setObjectName("label_16")
        self.horizontalLayout_10.addWidget(self.label_16)
        spacerItem8 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem8)
        self.verticalLayout_5.addLayout(self.horizontalLayout_10)
        self.horizontalLayout_11 = QtGui.QHBoxLayout()
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        self.firstLongitudeBox = QtGui.QDoubleSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.firstLongitudeBox.setFont(font)
        self.firstLongitudeBox.setDecimals(1)
        self.firstLongitudeBox.setMinimum(0.0)
        self.firstLongitudeBox.setMaximum(360.0)
        self.firstLongitudeBox.setSingleStep(0.1)
        self.firstLongitudeBox.setProperty("value", 200.0)
        self.firstLongitudeBox.setObjectName("firstLongitudeBox")
        self.horizontalLayout_11.addWidget(self.firstLongitudeBox)
        self.label_9 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.horizontalLayout_11.addWidget(self.label_9)
        spacerItem9 = QtGui.QSpacerItem(30, 20, QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_11.addItem(spacerItem9)
        self.secondLongitudeBox = QtGui.QDoubleSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.secondLongitudeBox.setFont(font)
        self.secondLongitudeBox.setDecimals(1)
        self.secondLongitudeBox.setMaximum(360.0)
        self.secondLongitudeBox.setSingleStep(0.1)
        self.secondLongitudeBox.setProperty("value", 230.0)
        self.secondLongitudeBox.setObjectName("secondLongitudeBox")
        self.horizontalLayout_11.addWidget(self.secondLongitudeBox)
        self.label_17 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_17.setFont(font)
        self.label_17.setObjectName("label_17")
        self.horizontalLayout_11.addWidget(self.label_17)
        spacerItem10 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_11.addItem(spacerItem10)
        self.verticalLayout_5.addLayout(self.horizontalLayout_11)
        self.horizontalLayout_12 = QtGui.QHBoxLayout()
        self.horizontalLayout_12.setObjectName("horizontalLayout_12")
        self.pressureCutOffBox = QtGui.QSpinBox(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.pressureCutOffBox.setFont(font)
        self.pressureCutOffBox.setMaximum(2100)
        self.pressureCutOffBox.setProperty("value", 970)
        self.pressureCutOffBox.setObjectName("pressureCutOffBox")
        self.horizontalLayout_12.addWidget(self.pressureCutOffBox)
        self.label_10 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.horizontalLayout_12.addWidget(self.label_10)
        self.pressureCutOffWarningLabel = QtGui.QLabel(self.settingsPage)
        self.pressureCutOffWarningLabel.setText("")
        self.pressureCutOffWarningLabel.setObjectName("pressureCutOffWarningLabel")
        self.horizontalLayout_12.addWidget(self.pressureCutOffWarningLabel)
        spacerItem11 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_12.addItem(spacerItem11)
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
        spacerItem12 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_13.addItem(spacerItem12)
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
        spacerItem13 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_13.addItem(spacerItem13)
        self.verticalLayout_5.addLayout(self.horizontalLayout_13)
        spacerItem14 = QtGui.QSpacerItem(20, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_5.addItem(spacerItem14)
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
        spacerItem15 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_15.addItem(spacerItem15)
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.horizontalLayout_14 = QtGui.QHBoxLayout()
        self.horizontalLayout_14.setObjectName("horizontalLayout_14")
        self.latitudeDesiredBox = QtGui.QDoubleSpinBox(self.settingsPage)
        self.latitudeDesiredBox.setEnabled(True)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.latitudeDesiredBox.setFont(font)
        self.latitudeDesiredBox.setReadOnly(True)
        self.latitudeDesiredBox.setButtonSymbols(QtGui.QAbstractSpinBox.NoButtons)
        self.latitudeDesiredBox.setMaximum(180.0)
        self.latitudeDesiredBox.setProperty("value", 52.5)
        self.latitudeDesiredBox.setObjectName("latitudeDesiredBox")
        self.horizontalLayout_14.addWidget(self.latitudeDesiredBox)
        self.longitudeDesiredBox = QtGui.QDoubleSpinBox(self.settingsPage)
        self.longitudeDesiredBox.setEnabled(True)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.longitudeDesiredBox.setFont(font)
        self.longitudeDesiredBox.setReadOnly(True)
        self.longitudeDesiredBox.setButtonSymbols(QtGui.QAbstractSpinBox.NoButtons)
        self.longitudeDesiredBox.setMaximum(360.0)
        self.longitudeDesiredBox.setProperty("value", 215.0)
        self.longitudeDesiredBox.setObjectName("longitudeDesiredBox")
        self.horizontalLayout_14.addWidget(self.longitudeDesiredBox)
        self.label_14 = QtGui.QLabel(self.settingsPage)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_14.setFont(font)
        self.label_14.setObjectName("label_14")
        self.horizontalLayout_14.addWidget(self.label_14)
        spacerItem16 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_14.addItem(spacerItem16)
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
        spacerItem17 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_4.addItem(spacerItem17)
        self.horizontalLayout_15.addLayout(self.verticalLayout_4)
        spacerItem18 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_15.addItem(spacerItem18)
        self.verticalLayout_5.addLayout(self.horizontalLayout_15)
        spacerItem19 = QtGui.QSpacerItem(20, 50, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_5.addItem(spacerItem19)
        self.horizontalLayout_25 = QtGui.QHBoxLayout()
        self.horizontalLayout_25.setObjectName("horizontalLayout_25")
        spacerItem20 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_25.addItem(spacerItem20)
        self.progressLabel = QtGui.QLabel(self.settingsPage)
        self.progressLabel.setObjectName("progressLabel")
        self.horizontalLayout_25.addWidget(self.progressLabel)
        self.verticalLayout_5.addLayout(self.horizontalLayout_25)
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        spacerItem21 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_7.addItem(spacerItem21)
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
        self.pleaseWaitPage = QtGui.QWidget()
        self.pleaseWaitPage.setObjectName("pleaseWaitPage")
        self.verticalLayout_8 = QtGui.QVBoxLayout(self.pleaseWaitPage)
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.label_18 = QtGui.QLabel(self.pleaseWaitPage)
        self.label_18.setObjectName("label_18")
        self.verticalLayout_8.addWidget(self.label_18)
        self.calculationProgressBar = QtGui.QProgressBar(self.pleaseWaitPage)
        self.calculationProgressBar.setProperty("value", 24)
        self.calculationProgressBar.setObjectName("calculationProgressBar")
        self.verticalLayout_8.addWidget(self.calculationProgressBar)
        self.timeSeriesStackedWidget.addWidget(self.pleaseWaitPage)
        self.calculationsPage = QtGui.QWidget()
        self.calculationsPage.setObjectName("calculationsPage")
        self.verticalLayout_7 = QtGui.QVBoxLayout(self.calculationsPage)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.outputFilesLabel_2 = QtGui.QLabel(self.calculationsPage)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.outputFilesLabel_2.setFont(font)
        self.outputFilesLabel_2.setObjectName("outputFilesLabel_2")
        self.verticalLayout_7.addWidget(self.outputFilesLabel_2)
        self.horizontalLayout_20 = QtGui.QHBoxLayout()
        self.horizontalLayout_20.setObjectName("horizontalLayout_20")
        self.plotTemperatureButton = QtGui.QPushButton(self.calculationsPage)
        self.plotTemperatureButton.setMinimumSize(QtCore.QSize(99, 0))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.plotTemperatureButton.setFont(font)
        self.plotTemperatureButton.setObjectName("plotTemperatureButton")
        self.horizontalLayout_20.addWidget(self.plotTemperatureButton)
        self.plotSalinityButton = QtGui.QPushButton(self.calculationsPage)
        self.plotSalinityButton.setMinimumSize(QtCore.QSize(99, 0))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.plotSalinityButton.setFont(font)
        self.plotSalinityButton.setObjectName("plotSalinityButton")
        self.horizontalLayout_20.addWidget(self.plotSalinityButton)
        spacerItem22 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_20.addItem(spacerItem22)
        self.verticalLayout_7.addLayout(self.horizontalLayout_20)
        self.horizontalLayout_23 = QtGui.QHBoxLayout()
        self.horizontalLayout_23.setObjectName("horizontalLayout_23")
        self.plotSigmaTButton = QtGui.QPushButton(self.calculationsPage)
        self.plotSigmaTButton.setMinimumSize(QtCore.QSize(99, 0))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.plotSigmaTButton.setFont(font)
        self.plotSigmaTButton.setObjectName("plotSigmaTButton")
        self.horizontalLayout_23.addWidget(self.plotSigmaTButton)
        self.plotSpicinessButton = QtGui.QPushButton(self.calculationsPage)
        self.plotSpicinessButton.setMinimumSize(QtCore.QSize(99, 0))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.plotSpicinessButton.setFont(font)
        self.plotSpicinessButton.setObjectName("plotSpicinessButton")
        self.horizontalLayout_23.addWidget(self.plotSpicinessButton)
        spacerItem23 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_23.addItem(spacerItem23)
        self.verticalLayout_7.addLayout(self.horizontalLayout_23)
        spacerItem24 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout_7.addItem(spacerItem24)
        self.horizontalLayout_24 = QtGui.QHBoxLayout()
        self.horizontalLayout_24.setObjectName("horizontalLayout_24")
        self.outputFilesLabel = QtGui.QLabel(self.calculationsPage)
        font = QtGui.QFont()
        font.setPointSize(11)
        self.outputFilesLabel.setFont(font)
        self.outputFilesLabel.setObjectName("outputFilesLabel")
        self.horizontalLayout_24.addWidget(self.outputFilesLabel)
        spacerItem25 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_24.addItem(spacerItem25)
        self.verticalLayout_7.addLayout(self.horizontalLayout_24)
        spacerItem26 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_7.addItem(spacerItem26)
        self.horizontalLayout_22 = QtGui.QHBoxLayout()
        self.horizontalLayout_22.setObjectName("horizontalLayout_22")
        spacerItem27 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_22.addItem(spacerItem27)
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
        self.horizontalLayout.addWidget(self.timeSeriesStackedWidget)

        self.retranslateUi(TimeSeriesApp)
        self.timeSeriesStackedWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(TimeSeriesApp)

    def retranslateUi(self, TimeSeriesApp):
        TimeSeriesApp.setWindowTitle(QtGui.QApplication.translate("TimeSeriesApp", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("TimeSeriesApp", "Date Range (Day/Month/Year)", None, QtGui.QApplication.UnicodeUTF8))
        self.startRangeDateEdit.setDisplayFormat(QtGui.QApplication.translate("TimeSeriesApp", "d/M/yyyy", None, QtGui.QApplication.UnicodeUTF8))
        self.label_15.setText(QtGui.QApplication.translate("TimeSeriesApp", "to", None, QtGui.QApplication.UnicodeUTF8))
        self.endRangeDateEdit.setDisplayFormat(QtGui.QApplication.translate("TimeSeriesApp", "d/M/yyyy", None, QtGui.QApplication.UnicodeUTF8))
        self.currentDayLabel.setText(QtGui.QApplication.translate("TimeSeriesApp", "Current Date:", None, QtGui.QApplication.UnicodeUTF8))
        self.currentDayDateEdit.setDisplayFormat(QtGui.QApplication.translate("TimeSeriesApp", "d/M/yyyy", None, QtGui.QApplication.UnicodeUTF8))
        self.currentDayRangeLabel_2.setText(QtGui.QApplication.translate("TimeSeriesApp", "Earliest Recommended Date: ", None, QtGui.QApplication.UnicodeUTF8))
        self.earliestDayDateEdit.setDisplayFormat(QtGui.QApplication.translate("TimeSeriesApp", "d/M/yyyy", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("TimeSeriesApp", "Date Step Size", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("TimeSeriesApp", "Range of days used on both sides of target date to interpolate float (Max 30)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("TimeSeriesApp", "Start Latitude", None, QtGui.QApplication.UnicodeUTF8))
        self.label_16.setText(QtGui.QApplication.translate("TimeSeriesApp", "End Latitude", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setText(QtGui.QApplication.translate("TimeSeriesApp", "Start Longitude", None, QtGui.QApplication.UnicodeUTF8))
        self.label_17.setText(QtGui.QApplication.translate("TimeSeriesApp", "End Longitude", None, QtGui.QApplication.UnicodeUTF8))
        self.label_10.setText(QtGui.QApplication.translate("TimeSeriesApp", "Minimum float depth during drift (m)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setText(QtGui.QApplication.translate("TimeSeriesApp", "Interpolate to depth (m)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_12.setText(QtGui.QApplication.translate("TimeSeriesApp", "Step size for depth/pressure (m)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_13.setText(QtGui.QApplication.translate("TimeSeriesApp", "Desired outputs:", None, QtGui.QApplication.UnicodeUTF8))
        self.tempCheckBox.setText(QtGui.QApplication.translate("TimeSeriesApp", "Temperature", None, QtGui.QApplication.UnicodeUTF8))
        self.salinityCheckBox.setText(QtGui.QApplication.translate("TimeSeriesApp", "Salinity", None, QtGui.QApplication.UnicodeUTF8))
        self.sigmaTCheckBox.setText(QtGui.QApplication.translate("TimeSeriesApp", "Sigma-T", None, QtGui.QApplication.UnicodeUTF8))
        self.spicinessCheckBox.setText(QtGui.QApplication.translate("TimeSeriesApp", "Spiciness", None, QtGui.QApplication.UnicodeUTF8))
        self.dynamicHeightCheckBox.setText(QtGui.QApplication.translate("TimeSeriesApp", "Dynamic Height", None, QtGui.QApplication.UnicodeUTF8))
        self.label_14.setText(QtGui.QApplication.translate("TimeSeriesApp", "Latitude and Longtitude of desired station", None, QtGui.QApplication.UnicodeUTF8))
        self.appendCheckBox.setText(QtGui.QApplication.translate("TimeSeriesApp", "Append to existing files", None, QtGui.QApplication.UnicodeUTF8))
        self.verboseCheckBox.setText(QtGui.QApplication.translate("TimeSeriesApp", "Verbose", None, QtGui.QApplication.UnicodeUTF8))
        self.progressLabel.setText(QtGui.QApplication.translate("TimeSeriesApp", "Waiting for Settings", None, QtGui.QApplication.UnicodeUTF8))
        self.backButton_2.setText(QtGui.QApplication.translate("TimeSeriesApp", "Back", None, QtGui.QApplication.UnicodeUTF8))
        self.nextButton.setText(QtGui.QApplication.translate("TimeSeriesApp", "Next", None, QtGui.QApplication.UnicodeUTF8))
        self.label_18.setText(QtGui.QApplication.translate("TimeSeriesApp", "Calculations in progress, please do not interrupt", None, QtGui.QApplication.UnicodeUTF8))
        self.outputFilesLabel_2.setText(QtGui.QApplication.translate("TimeSeriesApp", "Data available to plot:", None, QtGui.QApplication.UnicodeUTF8))
        self.plotTemperatureButton.setText(QtGui.QApplication.translate("TimeSeriesApp", "Temperature", None, QtGui.QApplication.UnicodeUTF8))
        self.plotSalinityButton.setText(QtGui.QApplication.translate("TimeSeriesApp", "Salinity", None, QtGui.QApplication.UnicodeUTF8))
        self.plotSigmaTButton.setText(QtGui.QApplication.translate("TimeSeriesApp", "Sigma-T", None, QtGui.QApplication.UnicodeUTF8))
        self.plotSpicinessButton.setText(QtGui.QApplication.translate("TimeSeriesApp", "Spiciness", None, QtGui.QApplication.UnicodeUTF8))
        self.outputFilesLabel.setText(QtGui.QApplication.translate("TimeSeriesApp", "Output Files:", None, QtGui.QApplication.UnicodeUTF8))
        self.backToSettingsButton.setText(QtGui.QApplication.translate("TimeSeriesApp", "Back", None, QtGui.QApplication.UnicodeUTF8))
        self.finishedButton.setText(QtGui.QApplication.translate("TimeSeriesApp", "Finished", None, QtGui.QApplication.UnicodeUTF8))

