''' 
TimeSeriesApp.py
Callum Colvine - Research Assistant
Callum.Colvine@dfo-mpo.gc.ca
CallumColvine@gmail.com

TimeSeries has the functionality from Howard Freeland's TimeSeries which was 
written in HT Basic.

TimeSeries reads ARGO data and outputs interpolated float data into several 
TS_*.csv files. 

Potential upgrades to TimeSeries if there is time in the futre:
- Eliminate all the individual data arrays, and use a dictionary instead
- Concurrently handle mutiple tasks to improve runtime
- Examine loop structure in an attempt to improve runtime

Following Pep 8 formatting with the following exceptions:
- There is no spacing between a docstring and a function
'''

# This is to deal with path issues for the sake of project organizaiton
# import sys
# sys.path.append("..")


# Utility File imports
from utilities.ARGO_Utilities import formatToDateTime, julianToDate, \
    dateToJulian, lastRun, todayInDate, getSigmaT, getSvanom, getSpiciness, \
    getProfile, checkPressureMonotonic, removeIndexFromPTS
# UI imports
from PySide import QtCore, QtGui
from PySide.QtGui import QWidget
from PySide.QtCore import QSize, QDate
# Calculations imports
import numpy as np
import os
import time
import csv
from datetime import datetime
import matplotlib.pyplot as plt
from math import ceil
# I/O imports
import configparser
from ui_Files.ui_timeseriesapp import Ui_TimeSeriesApp


# TESTING is used to define if the project should limit which files it takes in
# assuming the user does not have access to all ARGO files
TESTING = False
# SAVELOCALLY will allow the user to control if output is saved onto the network
# or if it is saved to a local directory. Useful for testing.
SAVELOCALLY = True


class TimeSeriesApp(QWidget, Ui_TimeSeriesApp):

    def __init__(self, parent):
        super(TimeSeriesApp, self).__init__()
        self.setupUi(self)
        return

    def experimentSelected(self):
        self.loadOldSettings()
        self.initAllClassVariables()
        self.lastRunCalls()
        self.setupSignals()
        self.userDefinedSettings()
        self.timeSeriesStackedWidget.setCurrentWidget(self.settingsPage)   
        self.progressLabel.setText("Waiting for Settings")     
        return
        
    ''' This method holds references to all class variables. For anyone unclear,
    Python allows any method inside the class to access these variables. They 
    work similarly to "Global" variables, but in the class scope.'''
    def initAllClassVariables(self):
        # Declare input parameter variables
        self.xCoord = 0
        self.startDate = self.startRangeDateEdit.date().getDate()
        self.endDate = self.endRangeDateEdit.date().getDate()
        self.firstDayTuple = formatToDateTime(self.startDate[0],
                                              self.startDate[1],
                                              self.startDate[2])
        self.lastDayTuple = formatToDateTime(self.endDate[0],
                                             self.endDate[1],
                                             self.endDate[2])
        self.dayStepSize = self.dayStepSizeBox.value()
        self.sampleWindow = self.sampleWindowBox.value()
        self.firstLatitude = self.firstLatitudeBox.value()
        self.secondLatitude = self.secondLatitudeBox.value()
        self.firstLongitude = self.firstLongitudeBox.value()
        self.secondLongitude = self.secondLongitudeBox.value()
        self.pressureCutOff = self.pressureCutOffBox.value()
        self.maxInterpDepth = self.maxInterpDepthBox.value()

        self.stepSize = self.stepSizeBox.value()
        self.nPress = 0
        self.updateNPress()
        self.temp = self.tempCheckBox.isChecked()
        self.salinity = self.salinityCheckBox.isChecked()
        self.sigmaT = self.sigmaTCheckBox.isChecked()
        self.spiciness = self.spicinessCheckBox.isChecked()
        self.dynamicHeight = self.dynamicHeightCheckBox.isChecked()

        self.updateDesiredLonLat()

        self.append = True
        self.verbose = False
        self.weight = 0

        self.tempCSV = None
        self.salinityCSV = None
        self.sigmaTCSV = None
        self.spicinessCSV = None
        self.dynamicHeightCSV = None
        self.hgtCSV = None
        # self.cLimCSV = None

        self.yearMonthDayPath0 = ''
        self.sCY = 111.2  # ! sCY = km/degree of latitude
        self.sCX0 = self.sCY * np.cos(self.latitudeDesired)
        self.rho0 = 300  # ! Decay scale of weighting function

        # Pdel never used
        self.arSva = np.zeros((1000))
        self.dynH = np.zeros((1000))
        # self.dynH[:] = 0

        if TESTING:
            self.curDay = 18
            self.curMonth = 05
            self.curYear = 2016
        else:
            self.curDay, self.curMonth, self.curYear = todayInDate()
        self.currentDayDateEdit.setDate(QDate(self.curYear,
                                              self.curMonth, 
                                              self.curDay))
        self.numFloats = 0

        self.Te = np.empty((500, 500))
        self.Sa = np.empty((500, 500))
        self.Lat = np.empty((500))
        self.Lon = np.empty((500))
        self.floatUsable = np.empty((500))

        self.Te[:] = 999.9
        self.Sa[:] = 999.9
        self.Lat[:] = 999.9
        self.Lon[:] = 999.9
        self.floatUsable[:] = False

        self.sigRefSigT = np.empty((71))
        self.sigRefTemp = np.empty((71))
        self.sigRefSal = np.empty((71))
        # 110! Dimensions of data to be read from files
        self.P = np.empty((2000))
        self.T = np.empty((2000))
        self.S = np.empty((2000))
        self.St = np.empty((2000))
        self.Sp = np.empty((2000))
        # 140 ! ALPHA PEN 1
        # 150 Tabrow=11
        self.closestDist = 9999.9
        self.closestFloatNum = ''

        # Max is 100, indicating completed
        self.progress = 0
            
        self.drive = "P:\\"
        self.path0 = self.drive + "argo_mirror\\pacific_ocean\\"
        self.outPath = self.drive + "argo_out_TEST\\TimeSeries\\"
        self.sigPath = self.drive + "projects\\Sigma_Climate\\"
        self.localSoftwarePath = ''

        # self.timeSeriesStackedWidget.setCurrentWidget(self.settingsPage)
        return

    def updateDesiredLonLat(self):
        self.latitudeDesired = \
            abs((self.secondLatitude + self.firstLatitude) / 2.0)
        self.latitudeDesiredBox.setValue(self.latitudeDesired)
        self.longitudeDesired = \
            abs((self.secondLongitude + self.firstLongitude) / 2.0)
        self.longitudeDesiredBox.setValue(self.longitudeDesired)
        return

    def lastRunCalls(self):
        if SAVELOCALLY:
            self.outPath = self.localSoftwarePath + "argo_out_TEST\\TimeSeries\\"
            self.sigPath = self.localSoftwarePath + "projects\\Sigma_Climate\\"
            lastRun(self.localSoftwarePath, 2)
        else:
            lastRun(self.drive, 2)
        # Made this file empty here when he writes 3x71 empties
        # if problems: change this maybe
        csvF = open((self.sigPath + "Mp26_i.csv"), 'w')
        csvF.close()
        return

    ''' Methods called by the user interacting with the GUI '''
    def setupSignals(self):
        # TimeSeries stuff
        self.setupInputParameterSignals()
        self.setupResultsPageSignals()
        self.nextButton.clicked.connect(self.nextButtonClicked)
        return

    def setupResultsPageSignals(self):
        self.plotTemperatureButton.clicked.connect(self.plotTempButtonClicked)
        self.plotSalinityButton.clicked.connect(self.plotSalinityButtonClicked)
        self.plotSigmaTButton.clicked.connect(self.plotSigmaTButtonClicked)
        self.plotSpicinessButton.clicked.connect(self.plotSpicinessButtonClicked)
        return

    def setupInputParameterSignals(self):
        self.dayStepSizeBox.editingFinished.connect(
            self.dayStepSizeBoxEditingFinished)
        self.sampleWindowBox.editingFinished.connect(
            self.sampleWindowBoxEditingFinished)
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
        self.tempCheckBox.stateChanged.connect(
            self.tempCheckBoxStateChanged)
        self.salinityCheckBox.stateChanged.connect(
            self.salinityCheckBoxStateChanged)
        self.sigmaTCheckBox.stateChanged.connect(
            self.sigmaTCheckBoxStateChanged)
        self.spicinessCheckBox.stateChanged.connect(
            self.spicinessCheckBoxStateChanged)
        self.dynamicHeightCheckBox.stateChanged.connect(
            self.dynamicHeightCheckBoxStateChanged)
        self.latitudeDesiredBox.editingFinished.connect(
            self.latitudeDesiredBoxEditingFinished)
        self.longitudeDesiredBox.editingFinished.connect(
            self.longitudeDesiredBoxEditingFinished)
        self.appendCheckBox.stateChanged.connect(
            self.appendCheckBoxStateChanged)
        self.verboseCheckBox.stateChanged.connect(
            self.verboseCheckBoxStateChanged)
        self.backToSettingsButton.clicked.connect(
            self.backToSettingsButtonClicked)
        self.startRangeDateEdit.dateChanged.connect(
            self.startDateChanged)
        self.endRangeDateEdit.dateChanged.connect(
            self.endDateChanged)
        return

    def driveEdited(self):
        self.path0 = self.drive + "argo_mirror\\pacific_ocean\\"
        self.outPath = self.drive + "argo_out_TEST\\TimeSeries\\"
        self.sigPath = self.drive + "projects\\Sigma_Climate\\"
        if SAVELOCALLY:
            self.outPath = self.localSoftwarePath + "argo_out_TEST\\TimeSeries\\"
            self.sigPath = self.localSoftwarePath + "projects\\Sigma_Climate\\"
            lastRun(self.localSoftwarePath, 2)
        else:
            lastRun(self.drive, 2)
        return

    ''' The following functions are called by button presses.
    They will not receive individual documentation'''
    # Start of parameter boxes
    def startDateChanged(self, date):
        self.startDate = date.getDate()
        self.firstDayTuple = formatToDateTime(self.startDate[0],
                                              self.startDate[1],
                                              self.startDate[2])
        return

    def endDateChanged(self, date):
        self.endDate = date.getDate()
        self.lastDayTuple = formatToDateTime(self.endDate[0],
                                             self.endDate[1],
                                             self.endDate[2])
        return

    # def defaultSettingsCheckBoxStateChanged(self, boxInput):
    #     if boxInput == 2:
    #         self.defaultOptions()
    #     else:
    #         self.userDefinedSettings()
    #     return

    def dayRangeBoxEditingFinished(self):
        # self.dayRange = self.dayRangeBox.value()
        return

    def dayStepSizeBoxEditingFinished(self):
        self.dayStepSize = self.dayStepSizeBox.value()
        return

    def sampleWindowBoxEditingFinished(self):
        self.sampleWindow = self.sampleWindowBox.value()
        return

    def firstLatitudeBoxEditingFinished(self):
        self.firstLatitude = self.firstLatitudeBox.value()
        self.updateDesiredLonLat()
        return

    def secondLatitudeBoxEditingFinished(self):
        self.secondLatitude = self.secondLatitudeBox.value()
        self.updateDesiredLonLat()
        return

    def firstLongitudeBoxEditingFinished(self):
        self.firstLongitude = self.firstLongitudeBox.value()
        self.updateDesiredLonLat()
        return

    def secondLongitudeBoxEditingFinished(self):
        self.secondLongitude = self.secondLongitudeBox.value()
        self.updateDesiredLonLat()
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

    def tempCheckBoxStateChanged(self, boxInput):
        self.temp = self.tempCheckBox.isChecked()
        return

    def salinityCheckBoxStateChanged(self):
        self.salinity = self.salinityCheckBox.isChecked()
        return

    def sigmaTCheckBoxStateChanged(self):
        self.sigmaT = self.sigmaTCheckBox.isChecked()
        return

    def spicinessCheckBoxStateChanged(self):
        self.spiciness = self.spicinessCheckBox.isChecked()
        return

    def dynamicHeightCheckBoxStateChanged(self):
        self.dynamicHeight = self.dynamicHeightCheckBox.isChecked()
        return

    def latitudeDesiredBoxEditingFinished(self):
        self.latitudeDesired = self.latitudeDesiredBox.value()
        self.sCX0 = self.sCY * np.cos(self.latitudeDesired)
        return

    def longitudeDesiredBoxEditingFinished(self):
        self.longitudeDesired = self.longitudeDesiredBox.value()
        return

    def appendCheckBoxStateChanged(self):
        self.append = self.appendCheckBox.isChecked()
        return

    def verboseCheckBoxStateChanged(self):
        self.verbose = self.verboseCheckBox.isChecked()
        return

    def backToSettingsButtonClicked(self):
        self.progressLabel.setText("Waiting for Settings")     
        self.timeSeriesStackedWidget.setCurrentWidget(self.settingsPage)
        return

    def nextButtonClicked(self):
        # Check to make sure all params have been filled
        # self.timeSeriesStackedWidget.setCurrentWidget(self.pleaseWaitPage)
        # time.sleep(2)
        print "Next clicked!!"
        # self.progressLabel.setText("Interpolation in progress. Please wait...")
        # print "Label text is set"
        self.saveUsedSettings()
        self.prepareOutputFiles()
        self.mainLoop()
        self.timeSeriesStackedWidget.setCurrentWidget(self.calculationsPage)
        return

    def plotTempButtonClicked(self):
        self.plotContour(self.plotTemp, self.plotTime, self.plotDepth)
        return

    def plotSalinityButtonClicked(self):
        self.plotContour(self.plotSalinity, self.plotTime, self.plotDepth)
        return

    def plotSigmaTButtonClicked(self):
        self.plotContour(self.plotSigmaT, self.plotTime, self.plotDepth)
        return

    def plotSpicinessButtonClicked(self):
        self.plotContour(self.plotSpiciness, self.plotTime, self.plotDepth)
        return


    def updateNPress(self):
        self.nPress = 1 + int(self.maxInterpDepth / self.stepSize)
        # if self.nPress > 300:
        #     self.nPress = 300
        return

    ''' Pulls from self.startDate and self.endDate to get the 2 corresponding 
    Julian dates '''
    def getJulianStartAndEnd(self):
        yearDayStart = self.firstDayTuple.timetuple().tm_yday
        yearDayEnd = self.lastDayTuple.timetuple().tm_yday
        julStart = dateToJulian(self.startDate[2],
                                self.startDate[1],
                                self.startDate[0],
                                yearDayStart)
        julEnd = dateToJulian(self.endDate[2],
                              self.endDate[1],
                              self.endDate[0],
                              yearDayEnd)
        return julStart, julEnd

    def sampleWindowBoxEditingFinished(self):
        self.sampleWindow = self.sampleWindowBox.value()
        return

    def updateProgress(self, iterDayNum, julStart, julEnd):
        if iterDayNum != 0:
            self.progress = (iterDayNum / 
                float(ceil((julEnd - julStart) / float(self.dayStepSize))))
            self.calculationProgressBar.setValue(self.progress)
        return

    def initPlotArrays(self, julStart, julEnd):
        # print "ceiling days is ", ceil((julEnd - julStart) / float(self.dayStepSize))
        self.plotTemp = np.empty((ceil((julEnd - julStart) / 
                                    float(self.dayStepSize)), 
                                ceil(self.maxInterpDepth / 
                                    float(self.stepSize)) + 1))
        self.plotSalinity = np.empty((ceil((julEnd - julStart) / 
                                        float(self.dayStepSize)), 
                                    ceil(self.maxInterpDepth / 
                                        float(self.stepSize)) + 1))
        self.plotSigmaT = np.empty((ceil((julEnd - julStart) / 
                                        float(self.dayStepSize)), 
                                    ceil(self.maxInterpDepth / 
                                        float(self.stepSize)) + 1))
        self.plotSpiciness = np.empty((ceil((julEnd - julStart) / 
                                        float(self.dayStepSize)), 
                                    ceil(self.maxInterpDepth / 
                                        float(self.stepSize)) + 1))
        self.plotTime = np.arange(julStart, julEnd, self.dayStepSize)
        self.plotDepth = np.arange(0, self.maxInterpDepth + 1, self.stepSize)
        return

    ''' This is the main loop for the program. Once the user fills out 
    parameters and hits the Next button, we commence interpolation. The rest of 
    the program is called from mainLoop(). We use *_index.CSV files
    to determine which floats to use, and then pull the data from each float and
    manipulate it.'''
    def mainLoop(self):
        julStart, julEnd = self.getJulianStartAndEnd()
        self.initPlotArrays(julStart, julEnd)
        iterDayNum = 0
        for iDay in xrange(julStart, julEnd, self.dayStepSize):
            self.updateProgress(iterDayNum, julStart, julEnd)
            self.numFloats = 0
            self.xCoord = iDay
            julWindowStart = iDay - self.sampleWindow
            julWindowEnd = iDay + self.sampleWindow
            floats = []
            for cycleJulDate in xrange(julWindowStart, julWindowEnd):
                # Get a list of all the names of floats to be used
                floats, self.yearMonthDayPath0 = \
                    self.checkFloatsFromIndex(floats, cycleJulDate, self.path0)
            sTav = 0
            sTavSq = 0
            sTKnt = 0
            j = 0
            for singleFloat in floats:
                numRecs = getProfile(singleFloat, self.P, self.T, self.S)
                numRecs = self.sanityCheck(numRecs)
                self.floatUsable[j] = True
                sigT = getSigmaT(self.S[0], self.T[0])
                sTav, sTavSq, sTKnt = self.calcStats(
                    sigT, sTav, sTavSq, sTKnt, numRecs, j)
                j += 1
                if TESTING:
                    break
            iPres75 = self.formatResults(sTav, sTavSq, sTKnt)
            self.removeEmpties()
            self.computeDynHeight(iterDayNum)
            self.finishUp(iPres75)
            iterDayNum += 1
        self.cleanUp()
        return

    ''' Not in a utils file due to repeated use of class variables. ''' 
    def checkFloatsFromIndex(self, floats, cycleJulDate, path0):
        day, month, year = julianToDate(cycleJulDate)
        if self.verbose:
            print "Day, month, year are ", day, month, year
        yearMonthDayPath0 = (path0 + 
                             str(year) + '\\' + 
                             str(month).zfill(2) + '\\' + 
                             str(day).zfill(2) + '\\')
        inFileName = (yearMonthDayPath0 + 
                      str(year) + 
                      str(month).zfill(2) + 
                      str(day).zfill(2) + 
                      '_index.csv')
        print "inFileName is ", inFileName
        try:
            with open((inFileName),
                      'rb') as indexCSV:
                if self.verbose:
                    print "Opened index ", inFileName
                reader = csv.reader(indexCSV)
                # Skip the first entry since it's a header
                reader.next()
                for row in reader:
                    # Filer out the ******* lines that are in some index files 
                    if row[3] == "*******":
                        continue
                    lat = float(row[1])
                    lon = float(row[2])
                    if lon < 0:
                        lon += 360
                    press = float(row[3])
                    if(lat > self.firstLatitude and 
                       lat < self.secondLatitude and 
                       lon > self.firstLongitude and 
                       lon < self.secondLongitude and 
                       press > self.pressureCutOff):
                        # print "So it passed, press is ", press
                        floatNum = row[0]
                        floats.append(yearMonthDayPath0 + row[0])
                        self.passedFloat(floatNum, lat, lon)
                        self.numFloats += 1
        except (OSError, IOError), e:
            print "File was not found: ", inFileName
        return floats, yearMonthDayPath0

    def passedFloat(self, floatNum, lat, lon):
        self.Lat[self.numFloats] = lat
        self.Lon[self.numFloats] = lon
        dX = self.sCX0 * (float(lon) - self.longitudeDesired)
        dY = self.sCY * (float(lat) - self.latitudeDesired)
        rho = np.sqrt(dX * dX + dY * dY)
        if self.closestDist > rho:
            self.closestDist = rho
            self.closestFloatNum = floatNum
        # ToDo: A few lines left out here, looks completely unnecessary
        return

    ''' Creates a contour plot popup window from the input parameters:
    plotInput: 2D array with data at correct time/depth related positions
    plotTime: 1D array '''
    def plotContour(self, plotInput, plotTime, plotDepth):
        # print "Length of plotInput is ", len(plotInput), "plotInput[0] ", len(plotInput[0])
        # print "Length of plotTime is ", len(plotTime), "plotDepth is ", len(plotDepth)
        if len(plotTime) < 2:
            print "Error: plotting data. There are too few time entries"
            return
        origin = 'lower'
        CS = plt.contourf(plotTime, 
            plotDepth, 
            plotInput.T, 
            linewidths=(3,),
            cmap=plt.cm.bone)
        plt.clabel(CS, fmt='%2.1f', colors='w', fontsize=14)
        plt.gca().invert_yaxis()

        plt.show()
        return

    ''' Closes open files that were used for the TimeSeries app '''
    def cleanUp(self):
        self.tempCSV.close()
        self.salinityCSV.close()
        self.sigmaTCSV.close()
        self.spicinessCSV.close()
        self.dynamicHeightCSV.close()
        self.hgtCSV.close()
        # self.cLimCSV.close()
        # self.stratCSV.close()
        print "----------------------------------------------------------------"
        print "Interpolations are completed and stored."
        if SAVELOCALLY:
            lastRun(self.localSoftwarePath, -1)
        else:
            lastRun(self.drive, -1)
        return

    ''' Writes to final output file TS_Shgt.csv'''
    def finishUp(self, iPres75):
        if self.dynamicHeight:
            self.hgtCSV.write(str(self.xCoord) + ',' +
                              str(self.dynH[0]) + ',' +
                              str(self.weight) + '\n')
        stdDiff = self.St[iPres75 - 1] - self.St[0]
        if stdDiff < -0.002:
            stdDiff = 0.0
        # self.stratCSV.write(str(self.xCoord) + ',' + str(stdDiff) + ',' + 
        #                     str(self.St[0]) + ',' + str(self.St[iPres75]) + ',' + 
        #                     str(self.weight) + '\n')
        stdDiff = 0.001 * int(1000.0 * stdDiff + 0.5)
        if stdDiff < 0.001:
            stdDiff = 0.000
        St1 = 0.001 * int(0.5 + 1000.0 * self.St[0])
        dayOut, monthOut, yearOut = julianToDate(self.xCoord)
        finalMessage = ("| " + str(self.xCoord) + " " + str(dayOut) + "/" + 
                        str(monthOut) +  "/" + str(yearOut) + " | " + 
                        str(self.weight) + "  " + str(stdDiff) + "    " + 
                        str(St1))
        self.closestDist = 0.1 * int(10.0 * self.closestDist + 0.5)
        finalMessage += (" | " + str(self.closestDist) + "  " + 
                         self.closestFloatNum + " |")
        # ToDo: Optimize, open outside of loop
        # lstMsge = open((self.outPath + "lstmsge.csv"), 'w')
        # lstMsge.truncate()
        # lstMsge.write(finalMessage)
        # lstMsge.close()
        return

    ''' Calculates some statistical analyses on the data provided'''
    def calcStats(self, sigT, sTav, sTavSq, sTKnt, numRecs, iFloat):
        # Copied, not sure why Howard does this
        sigT = int((0.5 + 1000 * sigT)) / 1000
        if sigT > 5 and sigT < 30:
            sTav += sigT
            sTavSq += np.power(sigT, 2)
            sTKnt += 1
        if self.P[0] < 20:
            self.P[0] = 0
        for iPress in xrange(0, self.nPress):
            # Interpolate to pressureCalc
            pressureCalc = self.stepSize * (iPress + 1)
            if (pressureCalc < self.P[0]): 
                continue
            if (pressureCalc > self.P[numRecs]):
                return sTav, sTavSq, sTKnt
            eps = 1.0E-6
            for k in xrange(1, numRecs):
                if pressureCalc >= self.P[k - 1] and pressureCalc <= self.P[k]:
                    deltaP = self.P[k] - self.P[k - 1]
                    if deltaP < 100:
                        rho = (pressureCalc - self.P[k - 1]) / \
                            (self.P[k] - self.P[k - 1] + eps)
                        self.Te[iFloat, iPress] = self.T[k - 1] + \
                            rho * (self.T[k] - self.T[k - 1] + eps)
                        self.Sa[iFloat, iPress] = self.S[k - 1] + \
                            rho * (self.S[k] - self.S[k - 1] + eps)
                    break
            # Finished interpolation to standard pressures for Float iFloat
        # Finished interpolation to standard pressures for all floats
        return sTav, sTavSq, sTKnt

    ''' Calculates values and adds them to P, T, S, St, and Sp 
    Line 814 '''
    def formatResults(self, sTav, sTavSq, sTKnt):
        if sTKnt != 0:
            sTav = sTav / sTKnt
            sTavSq = sTavSq / sTKnt
            stdDev = np.sqrt(sTavSq - sTav * sTav)

        kDel = 0
        qSal = 0
        qTemp = 0
        qSpice = 0
        iPres75 = int(1.1 + 75 / self.stepSize)
        weightSumTMax = 0
        for iterPress in xrange(0, self.nPress):
            pressureCal = self.stepSize * (iterPress - 1)  
            weightSumS = 0
            weightSumT = 0
            salWeightSumT = 0
            tempWeightSumT = 0
            for iterFloat in xrange(0, self.numFloats):
                if self.floatUsable[iterFloat] == False:
                    continue
                avgLat = (self.Lat[iterFloat] + self.latitudeDesired) / 2
                sCX = self.sCY * np.cos(avgLat)
                dX = sCX * (self.longitudeDesired - self.Lon[iterFloat])
                dY = self.sCY * (self.latitudeDesired - self.Lat[iterFloat])
                rho = np.sqrt(dX * dX + dY * dY)
                z = rho / self.rho0

                if not (z > 5):
                    self.weight = np.exp(-z * z)
                    if (self.Te[iterFloat, iterPress] > -1.5 and
                            self.Te[iterFloat, iterPress] < 40):
                        weightSumT = weightSumT + self.weight
                        tempWeightSumT = tempWeightSumT + \
                            self.weight * self.Te[iterFloat, iterPress]

                    if (self.Sa[iterFloat, iterPress] > 20 and 
                            self.Sa[iterFloat, iterPress] < 40):
                        weightSumS = weightSumS + self.weight
                        salWeightSumT = salWeightSumT + \
                            self.weight * self.Sa[iterFloat, iterPress]

            if weightSumT > weightSumTMax:
                weightSumTMax = weightSumT
            if weightSumT > 0:
                qTemp = tempWeightSumT / weightSumT
            else:
                qTemp = tempWeightSumT
            if weightSumS > 0:
                qSal = salWeightSumT / weightSumS
            else:
                pass
            self.weight = (int((weightSumTMax * 1000) + 0.5)) / 1000.0
            qSigmaT = getSigmaT(qSal, qTemp)
            qSpice = getSpiciness(qTemp, qSal)
            sigma, svan = getSvanom(qSal, qTemp, 0)
            self.P[iterPress] = self.stepSize * (iterPress - 1)
            self.T[iterPress] = qTemp
            self.S[iterPress] = qSal
            self.St[iterPress] = qSigmaT
            self.Sp[iterPress] = qSpice
            kDel += 1
            self.arSva[kDel] = svan
        return iPres75

    ''' Remove entries in P, T, S, St, and Sp that are unusable ''' 
    def removeEmpties(self):
        if self.St[0] > self.St[1] and np.abs(self.St[0] - self.St[1]) > .05:
            self.St[0] = self.St[1]
            self.T[0] = self.T[1]
            self.S[0] = self.S[1]
            self.Sp[0] = self.Sp[1]
        return

    ''' Writes the values from P, T, S, St, and Sp to their respective files '''
    def computeDynHeight(self, iterDayNum):
        # Now compute DHgt relative to Pmax
        self.dynH[self.nPress - 1] = 0
        q = 5.6E-6  # q=0.5f/g at station Papa
        for k in xrange(1, self.nPress):
            iPress = self.nPress - k - 1
            self.dynH[iPress] = ((self.dynH[iPress + 1]) + q * 
                (self.arSva[iPress + 1] + self.arSva[iPress]) * self.stepSize)

        iterPressNum = 0
        for iPress in xrange(0, self.nPress):
            qTemp = 0.001 * int(0.5 + 1000 * self.T[iPress])
            qSal = 0.001 * int(0.5 + 1000 * self.S[iPress])
            qSigmaT = 0.001 * int(0.5 + 1000 * self.St[iPress])
            qSpice = 0.001 * int(0.5 + 1000 * self.Sp[iPress])
            qDynamicHeight = 0.001 * int(0.5 + 1000 * self.dynH[iPress])

            pressCount = self.stepSize * (iPress)

            if(not (qSigmaT < self.sigRefSigT[0]) and 
                not (qSigmaT > self.sigRefSigT[70])):
                for iCl in xrange(0, 69):
                    if (qSigmaT >= self.sigRefSigT[iCl] and 
                    qSigmaT <= self.sigRefSigT[iCl + 1]):
                        tep, sap = self.foundPair()
                        break
                            
            xCoAndpressCount = str(self.xCoord) + ',' + str(-pressCount) + ',' 
            if self.temp:
                self.tempCSV.write(xCoAndpressCount + str(qTemp) + '\n')
                self.plotTemp[iterDayNum][iterPressNum] = qTemp
            if self.salinity:
                self.salinityCSV.write(xCoAndpressCount + str(qSal) + '\n')
                self.plotSalinity[iterDayNum][iterPressNum] = qSal
            if self.sigmaT:
                self.sigmaTCSV.write(xCoAndpressCount + str(qSigmaT) + '\n')
                self.plotSigmaT[iterDayNum][iterPressNum] = qSigmaT
            if self.spiciness:
                self.spicinessCSV.write(xCoAndpressCount + str(qSpice) + '\n')
                self.plotSpiciness[iterDayNum][iterPressNum] = qSpice
            if self.dynamicHeight:
                self.dynamicHeightCSV.write(xCoAndpressCount + 
                    str(qDynamicHeight) + '\n')
            iterPressNum += 1
        return 

    ''' When a pair is found, performs calculations on that pair and writes to 
    Sigref.csv'''
    # Variable names conflict with standard due to direct port of code of 
    # Howard's HT Basic.
    def foundPair(self, qSigmaT, iCl):
        Ra = (qSigmaT - self.sigRefSigT[iCl]) / \
            (self.sigRefSigT[iCl + 1] - self.sigRefSigT[iCl])
        Ter = self.sigRefTemp[iCl] + Ra * \
            (self.sigRefTemp[iCl + 1] - self.sigRefTemp[iCl])
        Sar = self.sigRefSal[iCl] + Ra * \
            (self.sigRefSal[iCl + 1] - self.sigRefSal[iCl])
        Tep = qTemp - Ter
        Sap = qSal - Sar
        Tep = 0.001 * int(0.5 + 1000.0 * Tep)
        Sap = 0.001 * int(0.5 + 1000.0 * Sap)
        # self.cLimCSV.write(qSigmaT + ',' + qTemp + ',' + Tep + ',' + qSal + 
        #     ',' + Sap)
        return Tep, Sap

    # Used for the "Default Settings" check box option. Was not found useful
    # def defaultOptions(self):
    #     self.setEnabledParameters(False)
    #     return

    ''' Determines potential validity of data ''' 
    def sanityCheck(self, numRecs):
        numRecs = self.checkOutOfRange(numRecs)
        numRecs = checkPressureMonotonic(numRecs, self.P, self.T, self.S)
        return numRecs

    def checkOutOfRange(self, numRecs):
        # Loop, if out of range, remove and adjust all values to cover up
        for i in xrange(0, numRecs):
            if (self.P[i] < -0.5 or
                self.P[i] > 2200 or
                self.T[i] < -2.5 or
                self.T[i] > 30 or
                self.S[i] < 30 or
                self.S[i] > 39):
                numRecs = removeIndexFromPTS(i, numRecs, self.P, self.T, self.S)
        return numRecs

    def setEnabledParameters(self, isEnabled):
        self.dayStepSizeBox.setEnabled(isEnabled)
        self.sampleWindowBox.setEnabled(isEnabled)
        self.firstLatitudeBox.setEnabled(isEnabled)
        self.secondLatitudeBox.setEnabled(isEnabled)
        self.firstLongitudeBox.setEnabled(isEnabled)
        self.secondLongitudeBox.setEnabled(isEnabled)
        self.pressureCutOffBox.setEnabled(isEnabled)
        self.maxInterpDepthBox.setEnabled(isEnabled)
        self.stepSizeBox.setEnabled(isEnabled)
        self.tempCheckBox.setEnabled(isEnabled)
        self.salinityCheckBox.setEnabled(isEnabled)
        self.sigmaTCheckBox.setEnabled(isEnabled)
        self.spicinessCheckBox.setEnabled(isEnabled)
        self.dynamicHeightCheckBox.setEnabled(isEnabled)
        self.latitudeDesiredBox.setEnabled(isEnabled)
        self.longitudeDesiredBox.setEnabled(isEnabled)
        self.appendCheckBox.setEnabled(isEnabled)
        self.verboseCheckBox.setEnabled(isEnabled)
        self.startRangeDateEdit.setEnabled(isEnabled)
        self.endRangeDateEdit.setEnabled(isEnabled)
        if (not isEnabled):
            # self.loadDefaults() ToDo
            pass
        return

    def userDefinedSettings(self):
        self.setEnabledParameters(True)
        # If the user has NOT selected default settings
        # self.writeDefinedSettings()
        return


    ''' Reads in the settings from the last run of TimeSeries from 
    TimeSeriesSettings.cfg '''
    def loadOldSettings(self):
        cfg = configparser.SafeConfigParser()
        cfg.read('TimeSeriesSettings.cfg')
        startDateRaw = cfg.get("inputParams", "dayRangeStart")
        startDateList = [x.strip() for x in startDateRaw.split(',')]
        self.startRangeDateEdit.setDate(QDate(int(startDateList[0]), 
                                              int(startDateList[1]), 
                                              int(startDateList[2])))
        endDateRaw = cfg.get("inputParams", "dayRangeEnd")
        endDateList = [x.strip() for x in endDateRaw.split(',')]
        self.endRangeDateEdit.setDate(QDate(int(endDateList[0]), 
                                            int(endDateList[1]), 
                                            int(endDateList[2])))
        self.dayStepSizeBox.setValue(cfg.getint("inputParams", "dayStepSize"))
        self.sampleWindowBox.setValue(cfg.getint("inputParams", "timeWindow"))
        self.firstLatitudeBox.setValue(cfg.getint("inputParams", "firstLatitude"))
        self.secondLatitudeBox.setValue(cfg.getint("inputParams", "secondLatitude"))
        self.firstLongitudeBox.setValue(cfg.getint("inputParams", "firstLongitude"))
        self.secondLongitudeBox.setValue(cfg.getint("inputParams", "secondLongitude"))
        self.pressureCutOffBox.setValue(cfg.getint("inputParams", "pressureCutoff"))
        self.maxInterpDepthBox.setValue(cfg.getint("inputParams", "maxPressure"))
        self.stepSizeBox.setValue(cfg.getint("inputParams", "changeInPressure"))

        self.tempCheckBox.setChecked(cfg.getboolean("desiredResults", "dispTemp"))
        self.salinityCheckBox.setChecked(cfg.getboolean("desiredResults", "dispSalinity"))
        self.sigmaTCheckBox.setChecked(cfg.getboolean("desiredResults", "dispSigmaT"))
        self.spicinessCheckBox.setChecked(cfg.getboolean("desiredResults", "dispSpiciness"))
        self.dynamicHeightCheckBox.setChecked(cfg.getboolean("desiredResults", "dispDynamicHeight"))
        self.latitudeDesiredBox.setValue(cfg.getfloat("desiredResults", "latitudeDesired"))                
        self.longitudeDesiredBox.setValue(cfg.getfloat("desiredResults", "longitudeDesired"))        
        return


    ''' Writes the settings from this run to TimeSeriesSettings.cfg '''
    # From my understanding and research, writing to a *.cfg file needs to 
    # completely overwrite file. Because of this, every option will always
    # be written. This isn't TOO inefficient since the file is very short
    def saveUsedSettings(self):
        print "Saving settings"
        cfg = configparser.SafeConfigParser()
        date = self.startRangeDateEdit.date().getDate()
        cfg.add_section("inputParams")
        cfg.set("inputParams", "dayRangeStart", str(date[0]) + ',' + 
                                                str(date[1]) + ',' + 
                                                str(date[2]))
        date = self.endRangeDateEdit.date().getDate()
        cfg.set("inputParams", "dayRangeEnd", str(date[0]) + ',' + 
                                              str(date[1]) + ',' + 
                                              str(date[2]))
        cfg.set("inputParams", "dayStepSize", str(self.dayStepSizeBox.value()))
        cfg.set("inputParams", "timeWindow", str(self.sampleWindowBox.value()))
        cfg.set("inputParams", "firstLatitude", str(self.firstLatitudeBox.value()))
        cfg.set("inputParams", "secondLatitude", str(self.secondLatitudeBox.value()))
        cfg.set("inputParams", "firstLongitude", str(self.firstLongitudeBox.value()))
        cfg.set("inputParams", "secondLongitude", str(self.secondLongitudeBox.value()))
        cfg.set("inputParams", "pressureCutoff", str(self.pressureCutOffBox.value()))
        cfg.set("inputParams", "maxPressure", str(self.maxInterpDepthBox.value()))
        cfg.set("inputParams", "changeInPressure", str(self.stepSizeBox.value()))
        
        cfg.add_section("desiredResults")
        cfg.set("desiredResults", "dispTemp", str(self.tempCheckBox.isChecked()))
        cfg.set("desiredResults", "dispSalinity", str(self.salinityCheckBox.isChecked()))
        cfg.set("desiredResults", "dispSigmaT", str(self.sigmaTCheckBox.isChecked()))
        cfg.set("desiredResults", "dispSpiciness", str(self.spicinessCheckBox.isChecked()))
        cfg.set("desiredResults", "dispDynamicHeight", str(self.dynamicHeightCheckBox.isChecked()))
        cfg.set("desiredResults", "latitudeDesired", str(self.latitudeDesiredBox.value()))
        cfg.set("desiredResults", "longitudeDesired", str(self.longitudeDesiredBox.value()))

        with open('TimeSeriesSettings.cfg', 'wb') as cfgFile:
            cfg.write(cfgFile)
        return

    def prepInfoFile(self, outLon, outLat, outDepth):
        exists = False
        sameParams = False
        amount = 1
        infoName = ("TS_Info_" + str(outLon) + '_' + str(outLat) + '_' + 
                    str(outDepth) + "m.txt")
        tempName = ("TS_Temp_" + str(outLon) + '_' + str(outLat) + '_' + 
                    str(outDepth) + "m")
        info = open((self.outPath + infoName), 'a+')
        # Check if the info file is empty
        if os.stat(self.outPath + infoName).st_size == 0:
            self.fillEmptyInfoFile(info)
            exists = False
        # The info file has contents
        else: 
            exists = True
            for line in info:
                if tempName in line:
                    splitLine = [x.strip() for x in line.split('| ')]
                    print "Split line is ", splitLine
                    # If the input filename does not match
                    if not (int(splitLine[1]) == self.stepSize and 
                            int(splitLine[2]) == self.sampleWindow):
                        amount += 1
                        isIn = True
                        sameParams = False
                    # There is an exact match in parameters and filename
                    else:
                        sameParams = True
                        # Return immediately because the rest do not matter,
                        # we found the file we're appending to
                        return exists, sameParams, info, amount
        return exists, sameParams, info, amount

    def fillEmptyInfoFile(self, info):
        info.write("Time Series is a program that interpolates data at a") 
        info.write("specific location and time, " + '\n')
        info.write("based on the data available from ARGO floats within")
        info.write("the selected nearby area." + '\n' + '\n')
        info.write("FILEDATA" + '\n' + "Format:" + '\n')
        info.write("Name_Data Type_Latitude of Data_Longitude of data_Maximum ")
        info.write("interpolation depth.csv, Pressure Step Size, Day Range")
        info.write('\n')
        return

    def prepOutFileNames(self, outLon, outLat, outDepth):
        tempName = ("TS_Temp_" + str(outLon) + '_' + str(outLat) + '_' + 
                                   str(outDepth) + "m")
        salName = ("TS_Salinity_" + str(outLon) + '_' + str(outLat) + '_' + 
                                               str(outDepth) + "m")
        sigTName = ("TS_SigmaT_" + str(outLon) + '_' + str(outLat) + '_' + 
                                               str(outDepth) + "m")
        spiceName = ("TS_Spiciness_" + str(outLon) + '_' + str(outLat) + '_' + 
                                               str(outDepth) + "m")
        dynHName = ("TS_DynamicHeight_" + str(outLon) + '_' + str(outLat) + '_' + 
                                               str(outDepth) + "m")
        surfHeightName = ("TS_SurfaceHeight_" + str(outLon) + '_' + 
                            str(outLat) + '_' + str(outDepth) + "m")
        return tempName, salName, sigTName, spiceName, dynHName, surfHeightName

    def addOffsets(self, tempName, salName, sigTName, spiceName, dynHName, 
                   surfHeightName, offset):
        tempName += '_' + str(offset)
        salName += '_' + str(offset)
        sigTName += '_' + str(offset)
        spiceName += '_' + str(offset)
        dynHName += '_' + str(offset)
        surfHeightName += '_' + str(offset)
        return tempName, salName, sigTName, spiceName, dynHName, surfHeightName

    def writeToInfo(self, tempName, salName, sigTName, spiceName, dynHName, 
                    surfHeightName, infoFile):
        infoFile.write(tempName + " | " + str(self.stepSize) + ' | ' + 
            str(self.sampleWindow) + '\n')
        infoFile.write(salName + " | " + str(self.stepSize) + ' | ' + 
            str(self.sampleWindow) + '\n')
        infoFile.write(sigTName + " | " + str(self.stepSize) + ' | ' + 
            str(self.sampleWindow) + '\n')
        infoFile.write(spiceName + " | " + str(self.stepSize) + ' | ' + 
            str(self.sampleWindow) + '\n')
        infoFile.write(dynHName + " | " + str(self.stepSize) + ' | ' + 
            str(self.sampleWindow) + '\n')
        infoFile.write(surfHeightName + " | " + str(self.stepSize) + ' | ' + 
            str(self.sampleWindow) + '\n')
        return

    ''' Opens and saves the output destination files in the class scope ''' 
    def prepareOutputFiles(self):

        outLon = self.longitudeDesired
        outLat = self.latitudeDesired
        outDepth = self.maxInterpDepth
        
        exists, sameParams, infoFile, offset = \
                self.prepInfoFile(outLon, outLat, outDepth)
        doOffset = False
        appendToInfo = False
        if not exists:
            doOffset = False
            appendToInfo = True
        elif exists and sameParams:
            doOffset = False
            appendToInfo = False
        if exists and not sameParams:
            doOffset = True
            appendToInfo = True

        tempName, salName, sigTName, spiceName, dynHName, surfHeightName = \
            self.prepOutFileNames(outLon, outLat, outDepth)
        if doOffset:
            tempName, salName, sigTName, spiceName, dynHName, surfHeightName = \
                self.addOffsets(tempName, salName, sigTName, spiceName, 
                                dynHName, surfHeightName, offset)
        tempName += ".csv"
        salName += ".csv"
        sigTName += ".csv"
        spiceName += ".csv"
        dynHName += ".csv"
        surfHeightName += ".csv"
        if appendToInfo:
            self.writeToInfo(tempName, salName, sigTName, spiceName, dynHName, 
                             surfHeightName, infoFile)

        tempPath = self.outPath + tempName
        salinityPath = self.outPath + salName
        sigmaTPath = self.outPath +  sigTName
        spicinessPath = self.outPath + spiceName
        dynamicHeightPath = self.outPath + dynHName
        hgtPath = self.outPath + surfHeightName         # Surface Height = Shgt


        # cLimPath = self.sigPath + "Sigref.csv"
        # stratPath = self.outPath + "Strat.csv"
        self.outputFilesLabel.setText("Output Files:\n" +
                                      tempPath + '\n' + 
                                      salinityPath + '\n' + 
                                      sigmaTPath + '\n' + 
                                      spicinessPath + '\n' + 
                                      dynamicHeightPath + '\n' + 
                                      hgtPath + '\n') 
                                      # stratPath + '\n' +
                                      # cLimPath)
        # self.outputLocationLabel.setText("Output Location:\n" + 
        #                                  self.outPath)
        writeType = 'a'
        if not self.append:
            writeType = 'w'
        self.tempCSV = open((tempPath), writeType)
        self.salinityCSV = open((salinityPath), writeType)
        self.sigmaTCSV = open((sigmaTPath), writeType)
        self.spicinessCSV = open((spicinessPath), writeType)
        self.dynamicHeightCSV = open((dynamicHeightPath), writeType)
        self.hgtCSV = open((hgtPath), writeType)
        # self.cLimCSV = open((cLimPath), writeType)
        # self.stratCSV = open((stratPath), writeType)
        if SAVELOCALLY:
            self.filesInUse = open((self.outPath + "IOS_Files_In_Use.txt"), 'w+')
        return
