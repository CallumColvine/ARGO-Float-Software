''' 
MainApp.py
Callum Colvine - Research Assistant
Callum.Colvine@dfo-mpo.gc.ca

Following Pep 8 formatting with the following exceptions:
- There is no spacing between a docstring and a function

MainApp hosts the GUI window surrounding the TimeSeries Application. TimeSeries
has the functionality from Howard Freeland's TimeSeries written in HT Basic. 

TimeSeries reads ARGO data and outputs interpolated float data into several 
TS_*.csv files. 

'''

# This is to deal with path issues for the sake of project organizaiton
import sys
sys.path.append("..")

# Utility File imports
from utilities.TimeSeriesUtilities import formatToDateTime, \
    julianToDate, dateToJulian
# UI imports
from PySide import QtCore, QtGui
from PySide.QtCore import QSize, QDate
from PySide.QtGui import QMainWindow
from ui_Files.ui_MainApp import Ui_MainApp
# Calculations imports
import numpy as np
import os
import time
import csv
from datetime import datetime
import matplotlib.pyplot as plt
from math import ceil


# TESTING is used to define if the project should limit which files it takes in
# assuming the user does not have access to all ARGO files
TESTING = False
# SAVELOCALLY will allow the user to control if output is saved onto the network
# or if it is saved to a local directory. Useful for testing.
SAVELOCALLY = True

class MainApp(QMainWindow, Ui_MainApp):

    def __init__(self):
        super(MainApp, self).__init__()
        self.setupUi(self)
        self.setupSignals()
        self.initAllClassVariables()
        self.userDefinedSettings()
        return

    ''' This method holds references to all class variables. For anyone unclear,
    Python allows any method inside the class to access these variables. They 
    work similarly to "Global" variables, but in the class scope.'''
    def initAllClassVariables(self):
        # Declare input parameter variables
        self.defaultSettings = False
        self.xCoord = 0
        self.floatRadius = 0
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

        self.latitudeDesired = self.latitudeDesiredBox.value()
        self.longitudeDesired = self.longitudeDesiredBox.value()
        self.append = True
        self.verbose = False
        self.weight = 0

        self.tempCSV = None
        self.salinityCSV = None
        self.sigmaTCSV = None
        self.spicinessCSV = None
        self.dynamicHeightCSV = None
        self.hgtCSV = None
        self.cLimCSV = None

        self.yearMonthDayPath0 = ''
        self.Scy = 111.2  # ! Scy = km/degree of latitude
        self.Scx0 = self.Scy * np.cos(self.latitudeDesired)
        self.Rho0 = 300  # ! Decay scale of weighting function
        self.Pref = 1000  # ! Reference level for dynamic height calculations
        # Pdel never used
        # self.Pdel = np.empty((1000))
        self.arSva = np.zeros((1000))
        self.Dh = np.zeros((1000))
        # self.Dh[:] = 0

        if TESTING:
            self.curDay = 18
            self.curMonth = 05
            self.curYear = 2016
        else:
            self.curDay = int(time.strftime("%d"))
            self.curMonth = int(time.strftime("%m"))
            self.curYear = int(time.strftime("%Y"))
        todayInJul = self.todayInJulian()
        day, month, year = self.todayInDate()
        self.currentDayDateEdit.setDate(QDate(year, month, day))
        self.numFloats = 0

        self.Te = np.empty((500, 300))
        self.Sa = np.empty((500, 300))
        self.Lat = np.empty((500))
        self.Lon = np.empty((500))
        self.Accpt = np.empty((500))
        self.Fltnm = np.empty((500))
        self.Te[:] = 999.9
        self.Sa[:] = 999.9
        self.Lat[:] = 999.9
        self.Lon[:] = 999.9
        self.Accpt[:] = (-1)

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
        self.Closest = 9999.9
        self.Closest_fltnm = ''

        self.Rho_print = -1
        # 240 ! Initialize data matrices with nul values 999.9
        # Should I do the NumPy NULL or just stick with 999.9?
        # Double backslash because in Python you need to escape special
        # characters

        # Max is 100, indicating completed
        self.progress = 0
            
        self.drive = "P:\\"
        self.path0 = self.drive + "argo_mirror\\pacific_ocean\\"
        self.outPath = self.drive + "argo_out_TEST\\TimeSeries\\"
        self.sigPath = self.drive + "projects\\Sigma_Climate\\"
        if SAVELOCALLY:
            self.localSoftwarePath = "C:\\Users\\ColvineC\\IOS_DFO\\ARGO-Float-Software\\"
            self.outPath = self.localSoftwarePath + "argo_out_TEST\\TimeSeries\\"
            self.sigPath = self.localSoftwarePath + "projects\\Sigma_Climate\\"
            self.lastRun(self.localSoftwarePath, 2)
        else:
            self.lastRun(self.drive, 2)
        # Made this file empty here when he writes 3x71 empties
        # if problems: change this maybe
        csvF = open((self.sigPath + "Mp26_i.csv"), 'w')
        csvF.close()
        # self.timeSeriesStackedWidget.setCurrentWidget(self.settingsPage)
        return

    ''' Methods called by the user interacting with the GUI '''
    def setupSignals(self):
        self.exitButton.clicked.connect(self.close)
        self.timeSeriesButton.clicked.connect(self.timeSeriesButtonClicked)
        # TimeSeries stuff
        self.setupInputParameterSignals()
        self.backToProgramListButton.clicked.connect(
            self.backToProgramListButtonClicked)
        self.nextButton.clicked.connect(self.nextButtonClicked)
        return

    def setupInputParameterSignals(self):
        self.defaultSettingsCheckBox.stateChanged.connect(
            self.defaultSettingsCheckBoxStateChanged)
        self.floatRadiusBox.editingFinished.connect(
            self.floatRadiusBoxEditingFinished)
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

    def defaultSettingsCheckBoxStateChanged(self, boxInput):
        if boxInput == 2:
            self.defaultOptions()
        else:
            self.userDefinedSettings()
        return

    def floatRadiusBoxEditingFinished(self):
        self.floatRadius = self.floatRadiusBox.value()
        return

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
        self.pressureCutOff = self.pressureCutOffBox.value()
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
        self.Scx0 = self.Scy * np.cos(self.latitudeDesired)
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
        self.timeSeriesStackedWidget.setCurrentWidget(self.settingsPage)
        return

    def backToProgramListButtonClicked(self):
        self.listAllPages.setCurrentWidget(self.programListPage)
        return

    def timeSeriesButtonClicked(self):
        self.listAllPages.setCurrentWidget(self.timeSeriesPage)
        return

    def nextButtonClicked(self):
        # Check to make sure all params have been filled
        self.timeSeriesStackedWidget.setCurrentWidget(self.pleaseWaitPage)
        # time.sleep(2)
        self.prepareOutputFiles()
        self.commenceInterpolation()
        self.timeSeriesStackedWidget.setCurrentWidget(self.calculationsPage)
        return

    def updateNPress(self):
        self.nPress = 1 + int(self.maxInterpDepth / self.stepSize)
        if self.nPress > 300:
            self.nPress = 300
        return

    ''' ToDo: verify integrity of data'''
    def sanityCheck(self):

        return

    ''' Heavy calculations originated from Howard's code. Will not be converted
    to naming standards '''
    def getSigmaT(self, S, T):
        # 12640! T  = Temperature in C
        # 12650! S  = salinity in pss-78
        R4 = 4.8314E-4
        Dr350 = 28.106331
        salSqRoot = np.sqrt(S)
        R1 = ((((6.536332E-9 * T - 1.120083E-6) * T + 1.001685E-4)
               * T - 9.09529E-3) * T + 6.793952E-2) * T - 28.263737
        R2 = (((5.3875E-9 * T - 8.2467E-7) * T + 7.6438E-5)
              * T - 4.0899E-3) * T + 8.24493E-1
        R3 = (-1.6546E-6 * T + 1.0227E-4) * T - 5.72466E-3
        Sig = (R4 * S + R3 * salSqRoot + R2) * S + R1
        Sigma = Sig + Dr350
        Sigma_t = ((Sigma + 1000.) / .999975) - 1000.
        return Sigma_t

    ''' Called on float numbers that have been found in the index file. Get 
    profile will determine the order in which and fill up the P, T, S arrays '''
    def getProfile(self,
                   floatPath,
                   rFlag):
        rFlag = -1
        order = []
        numChannels = 0
        numRecords = 0
        fileTempQc = 0
        filePSalQc = 0
        numRecs = self.extractDataFromFloatFile(floatPath, order)
        if not self.checkIfReturn(numRecs, fileTempQc, filePSalQc):
            print "in getProfile, returning prematurely due to bad values"
            return
        return numRecs

    ''' Used by getProfile() to determine validity of float data '''
    def extractDataFromFloatFile(self, path, order):
        # print "Opening path ", path
        floatFile = open(path, 'r')
        passedEndOfHeader = False
        pFound = False
        tFound = False
        sFound = False
        numRecs = 0
        i = 0
        for line in floatFile:
            if not passedEndOfHeader:
                if line.find("NUMBER OF RECORDS") != -1:
                    numRecs = self.getNumRecs(line)
                order, pFound, tFound, sFound = \
                    self.findOrder(line, order, pFound, tFound, sFound)
                if not (self.verifyUsability(line)):
                    return
            else:
                i += 1
                self.appendInfo(line, order, i)
            if line.find("END OF HEADER") != -1:
                passedEndOfHeader = True
        return numRecs

    def getNumRecs(self, line):
        split = line.split()
        # print "numRecs is ", int(split[4])
        return int(split[4])

    ''' Adds data pulled form individual float files to the P, T, S arrays '''
    def appendInfo(self, line, order, i):
        split = line.split()
        # print "Split is ", split
        offset = 5
        position = 0
        # if order != ['P', 'T', 'S']:
        #     print "----- Order is different -----"
        for data in order:
            twoPos = float(split[2 + (position * 5)])
            zeroPos = float(split[0 + (position * 5)])
            if data == 'P':
                if np.abs(twoPos) < 9000:
                    self.P[i] = twoPos
                else:
                    self.P[i] = zeroPos
            elif data == 'T':
                if np.abs(twoPos) < 9000:
                    self.T[i] = twoPos
                else:
                    self.T[i] = zeroPos
            elif data == 'S':
                if np.abs(twoPos) < 9000:
                    self.S[i] = twoPos
                else:
                    self.S[i] = zeroPos
            position += 1

        return

    ''' Makes sure float has minimum requirements for data '''
    def verifyUsability(self, line):
        if line.find("NUMBER OF RECORDS") != -1:
            if line.split()[2] < 15:
                return False
        if line.find("NUMBER OF CHANNELS") != -1:
            if line.split()[2] < 8:
                return False
        return True

    ''' Used to determine the order of data held inside each float file'''
    # TODO: Improve runtime, reduce calls to this function
    def findOrder(self, line, order, pFound, tFound, sFound):
        resP = -1
        resT = -1
        resS = -1
        # print line
        if not pFound:
            resP = line.find("PRES_ADJUSTED ")
        if not tFound:
            resT = line.find("TEMP_ADJUSTED ")
        if not sFound:
            resS = line.find("PSAL_ADJUSTED ")
        if resP != -1:
            pFound = True
            order.append('P')
        elif resT != -1:
            tFound = True
            order.append('T')
        elif resS != -1:
            sFound = True
            order.append('S')
        return order, pFound, tFound, sFound

    ''' Quits dataset of float if it does not have all the required components '''
    def checkIfReturn(self, numRecords, fileTempQc, filePSalQc):
        if numRecords < 15:
            return False
        if fileTempQc == 'C':
            return False
        if fileTempQc == 'D':
            return False
        if fileTempQc == 'E':
            return False
        if fileTempQc == 'F':
            return False
        if filePSalQc == 'C':
            return False
        if filePSalQc == 'D':
            return False
        if filePSalQc == 'E':
            return False
        if filePSalQc == 'F':
            return False
        return True
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

    def updateProgress(self, iterDayNum, julStart, julEnd):
        if iterDayNum != 0:
            self.progress = (iterDayNum / 
                float(ceil((julEnd - julStart) / float(self.dayStepSize))))
            self.calculationProgressBar.setValue(self.progress)
        return

    def initPlotArrays(self, julStart, julEnd):
        # print "ceiling days is ", ceil((julEnd - julStart) / float(self.dayStepSize))
        plotTemp = np.empty((ceil((julEnd - julStart) / float(self.dayStepSize)), 
            ceil(self.maxInterpDepth / float(self.stepSize)) + 1))
        plotTime = np.arange(julStart, julEnd, self.dayStepSize)
        plotDepth = np.arange(0, self.maxInterpDepth + 1, self.stepSize)
        # plotDepth = plotDepth0[::-1]
        return plotTemp, plotTime, plotDepth

    ''' This is the main loop for the program. Once the user fills out 
    parameters and hits the Next button, we commence interpolation. The rest of 
    the program is called from commenceInterpolation(). We use *_index.CSV files
    to determine which floats to use, and then pull the data from each float and
    manipulate it.'''
    def commenceInterpolation(self):
        julStart, julEnd = self.getJulianStartAndEnd()
        plotTemp, plotTime, plotDepth = self.initPlotArrays(julStart, julEnd)
        # julStart = 0
        # julEnd = 0
        iterDayNum = 0
        for Dc in xrange(julStart, julEnd, self.dayStepSize):
            self.updateProgress(iterDayNum, julStart, julEnd)
            self.numFloats = 0
            self.xCoord = Dc
            julWindowStart = Dc - self.sampleWindow
            julWindowEnd = Dc + self.sampleWindow
            floats = []
            for cycleJulDate in xrange(julWindowStart, julWindowEnd):
                self.checkFloatsFromIndex(floats, cycleJulDate)
            # print "The number of floats inside the box", self.numFloats
            sTav = 0
            sTavSq = 0
            sTKnt = 0
            rFlag = -1
            j = 0
            for singleFloat in floats:
                numRecs = self.getProfile(singleFloat, rFlag)
                self.sanityCheck()
                self.Accpt[j] = True
                sigT = self.getSigmaT(self.S[0], self.T[0])
                sTav, sTavSq, sTKnt = self.calcStats(
                    sigT, sTav, sTavSq, sTKnt, numRecs, j)
                j += 1
                if TESTING:
                    break

            Ipr75 = self.formatResults(sTav, sTavSq, sTKnt)
            self.removeEmpties()
            self.computeDHgt(plotTemp, plotTime, plotDepth, iterDayNum)
            self.finishUp(Ipr75)
            iterDayNum += 1
        # self.plotData()
        self.plotContour(plotTemp, plotTime, plotDepth)
        self.cleanUp()
        return


    def plotContour(self, plotTemp, plotTime, plotDepth):
        if len(plotTime) == 1:
            print "Error: plotting temperatures. There is only 1 time entry"
            return

        origin = 'lower'
        # print "plotTemp dimensions", len(plotTemp), len(plotTemp[0]), plotTemp
        # print "plotTime dimensions", len(plotTime), plotTime
        # print "plotDepth dimensions", len(plotDepth), plotDepth
        CS = plt.contourf(plotTime, 
            plotDepth, 
            plotTemp.T, 
            # [2, 3, 4, 5, 6, 7, 8], 
            linewidths=(3,),
            cmap=plt.cm.bone)
        plt.clabel(CS, fmt='%2.1f', colors='w', fontsize=14)
        plt.gca().invert_yaxis()

        plt.show()
        return

    def plotData(self):
        plt.plot(self.plotT)
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
        self.cLimCSV.close()
        self.stratCSV.close()
        print "----------------------------------------------------------------"
        print "Interpolations are completed and stored."
        if SAVELOCALLY:
            self.lastRun(self.localSoftwarePath, -1)
        else:
            self.lastRun(self.drive, -1)
        return

    ''' Writes to final output files and write to status file that the program 
    was successful '''
    def finishUp(self, Ipr75):
        if self.dynamicHeight:
            self.hgtCSV.write(str(self.xCoord) + ',' +
                              str(self.Dh[0]) + ',' +
                              str(self.weight))
        Stdiff = self.St[Ipr75] - self.St[1]
        if Stdiff < -.002:
            Stdiff = 0.
        self.stratCSV.write(str(self.xCoord) + ',' + str(self.St[1]) + ','
                            + str(self.St[Ipr75]) + ',' + str(self.weight))

        Stdiff = .001 * int(1000. * Stdiff + .5)
        if Stdiff < .001:
            Stdiff = 0.000
        St1 = .001 * int(.5 + 1000. * self.St[1])
        Dayc, Monc, Yearc = julianToDate(self.xCoord)
        R = "| " + str(self.xCoord) + " " + str(Dayc) + "/" + str(Monc) + \
            "/" + str(Yearc) + " | " + str(self.weight) + "  " + str(Stdiff) + \
            "    " + str(St1)
        self.Closest = .1 * int(10. * self.Closest + .5)
        R += R + " | " + str(self.Closest) + "  " + \
            self.Closest_fltnm + " |"
        # ToDo: Optimize, open outside of loop
        lstMsge = open((self.outPath + "lstmsge.csv"), 'w+')
        lstMsge.truncate()
        lstMsge.write(R)
        lstMsge.close()
        return

    ''' Calculates some statistical analyses on the data provided'''
    def calcStats(self, sigT, sTav, sTavSq, sTKnt, numRecs, Iflt):
        # Copied, not sure why Howard does this
        sigT = int((0.5 + 1000 * sigT)) / 1000
        if sigT > 5 and sigT < 30:
            sTav += sigT
            sTavSq += np.power(sigT, 2)
            sTKnt += 1
        if self.P[0] < 20:
            self.P[0] = 0
        for Ipr in range(0, self.nPress):
            # Interpolate to pressureCalc
            pressureCalc = self.stepSize * (Ipr + 1)
            if (pressureCalc < self.P[0]): 
                continue
            if (pressureCalc > self.P[numRecs]):
                return sTav, sTavSq, sTKnt
            Eps = 1.0E-6
            for K in range(1, numRecs):
                if pressureCalc >= self.P[K - 1] and pressureCalc <= self.P[K]:
                    Del_p = self.P[K] - self.P[K - 1]
                    if Del_p < 100:
                        Rho = (pressureCalc - self.P[K - 1]) / \
                            (self.P[K] - self.P[K - 1] + Eps)
                        self.Te[Iflt, Ipr] = self.T[K - 1] + \
                            Rho * (self.T[K] - self.T[K - 1] + Eps)
                        self.Sa[Iflt, Ipr] = self.S[K - 1] + \
                            Rho * (self.S[K] - self.S[K - 1] + Eps)
                    break
            # Finished interpolation to standard pressures for Float Iflt
        # Finished interpolation to standard pressures for all floats
        return sTav, sTavSq, sTKnt

    ''' Calculates values and adds them to P, T, S, St, and Sp 
    Line 814 '''
    def formatResults(self, sTav, sTavSq, sTKnt):
        try:
            sTav = sTav / sTKnt
            sTavSq = sTavSq / sTKnt
            stdDev = np.sqrt(sTavSq - sTav * sTav)
        except ZeroDivisionError, e:
            print "Error: No floats in this day range set match parameters"
        Kdel = 0
        qSal = 0
        qTemp = 0
        qSpice = 0
        Ipr75 = int(1.1 + 75 / self.stepSize)
        weightSumTMax = 0
        for iterPress in xrange(0, self.nPress):
            pressureCal = self.stepSize * (iterPress - 1)  
            weightSumS = 0
            weightSumT = 0
            salWeightSumT = 0
            tempWeightSumT = 0
            for iterFloat in xrange(0, self.numFloats):
                if self.Accpt[iterFloat] == False:
                    continue
                Latav = (self.Lat[iterFloat] + self.latitudeDesired) / 2
                Scx = self.Scy * np.cos(Latav)
                Dx = Scx * (self.longitudeDesired - self.Lon[iterFloat])
                Dy = self.Scy * (self.latitudeDesired - self.Lat[iterFloat])
                Rho = np.sqrt(Dx * Dx + Dy * Dy)
                Z = Rho / self.Rho0

                if not (Z > 5):
                    self.weight = np.exp(-Z * Z)
                    if (self.Te[iterFloat, iterPress] > -1.5 and
                            self.Te[iterFloat, iterPress] < 40):
                        weightSumT = weightSumT + self.weight
                        tempWeightSumT = tempWeightSumT + \
                            self.weight * self.Te[iterFloat, iterPress]

                    if self.Sa[iterFloat, iterPress] > 20 and \
                            self.Sa[iterFloat, iterPress] < 40:
                        weightSumS = weightSumS + self.weight
                        salWeightSumT = salWeightSumT + \
                            self.weight * self.Sa[iterFloat, iterPress]

            if weightSumT > weightSumTMax:
                weightSumTMax = weightSumT
            if weightSumT > 0:
                qTemp = tempWeightSumT / weightSumT
            else:
                qTemp = tempWeightSumT
                print "float not usable to calculate qTemp"
            if weightSumS > 0:
                qSal = salWeightSumT / weightSumS
            else:
                print "float not usable to calculate qSal"

            qSigmaT = self.getSigmaT(qSal, qTemp)
            qSpice = self.getSpiciness(qTemp, qSal)
            Sigma, Svan = self.getSvanom(qSal, qTemp, 0)
            self.P[iterPress] = self.stepSize * (iterPress - 1)
            self.T[iterPress] = qTemp
            self.S[iterPress] = qSal
            self.St[iterPress] = qSigmaT
            self.Sp[iterPress] = qSpice
            Kdel += 1
            self.arSva[Kdel] = Svan
        return Ipr75

    ''' Remove entries in P, T, S, St, and Sp that are unusable ''' 
    def removeEmpties(self):
        if self.St[0] > self.St[1] and np.abs(self.St[0] - self.St[1]) > .05:
            self.St[0] = self.St[1]
            self.T[0] = self.T[1]
            self.S[0] = self.S[1]
            self.Sp[0] = self.Sp[1]
        return

    ''' Writes the values from P, T, S, St, and Sp to their respective files '''
    def computeDHgt(self, plotTemp, plotTime, plotDepth, iterDayNum):
        # Now compute DHgt relative to Pmax
        self.Dh[self.nPress - 1] = 0
        Q = 5.6E-6  # Q=0.5f/g at station Papa
        for K in xrange(1, self.nPress):
            Ipr = self.nPress - K - 1
            self.Dh[Ipr] = ((self.Dh[Ipr + 1]) + Q * 
                (self.arSva[Ipr + 1] + self.arSva[Ipr]) * self.stepSize)

        iterPressNum = 0
        for Ipr in xrange(0, self.nPress):
            qTemp = 0.001 * int(0.5 + 1000 * self.T[Ipr])
            qSal = 0.001 * int(0.5 + 1000 * self.S[Ipr])
            qSigmaT = 0.001 * int(0.5 + 1000 * self.St[Ipr])
            qSpice = 0.001 * int(0.5 + 1000 * self.Sp[Ipr])
            qDynamicHeight = 0.001 * int(0.5 + 1000 * self.Dh[Ipr])

            Pc = self.stepSize * (Ipr)

            Sigref_opt = -1
            if(not (qSigmaT < self.sigRefSigT[0]) and 
                not (qSigmaT > self.sigRefSigT[70])):
                for Icl in range(0, 69):
                    if (qSigmaT >= self.sigRefSigT[Icl] and 
                    qSigmaT <= self.sigRefSigT[Icl + 1]):
                        tep, sap = self.foundPair()
                        break
                            
            xCoAndPc = str(self.xCoord) + ',' + str(-Pc) + ',' 
            if self.temp:
                self.tempCSV.write(xCoAndPc + str(qTemp) + '\n')
                plotTemp[iterDayNum][iterPressNum] = qTemp
            if self.salinity:
                self.salinityCSV.write(xCoAndPc + str(qSal) + '\n')
            if self.sigmaT:
                self.sigmaTCSV.write(xCoAndPc + str(qSigmaT) + '\n')
            if self.spiciness:
                self.spicinessCSV.write(xCoAndPc + str(qSpice) + '\n')
            if self.dynamicHeight:
                self.dynamicHeightCSV.write(xCoAndPc + str(qDynamicHeight) + '\n')
            iterPressNum += 1
        return 

    def foundPair(self, qSigmaT, Icl):
        Ra = (qSigmaT - self.sigRefSigT[Icl]) / \
            (self.sigRefSigT[Icl + 1] - self.sigRefSigT[Icl])
        Ter = self.sigRefTemp[Icl] + Ra * \
            (self.sigRefTemp[Icl + 1] - self.sigRefTemp[Icl])
        Sar = self.sigRefSal[Icl] + Ra * \
            (self.sigRefSal[Icl + 1] - self.sigRefSal[Icl])
        Tep = qTemp - Ter
        Sap = qSal - Sar
        Tep = .001 * int(.5 + 1000. * Tep)
        Sap = .001 * int(.5 + 1000. * Sap)
        self.cLimCSV.write(qSigmaT + ',' + qTemp + ',' + Tep + ',' + qSal + 
            ',' + Sap)
        return Tep, Sap

    def getSvanom(self, S, T, P0):
        # Compute the density anomaly, sigma, in kg/m^3
        # Density anomaly is identical with sigma-t without pressure terms
        # P0 = Pressure in decibars
        # T  = Temperature in deg C
        # S  = salinity in pss-78
        R3500 = 1028.106331
        R4 = 4.8314E-4
        Dr350 = 28.106331
        P = P0 / 10
        salSqRoot = np.sqrt(S)
        R1 = ((((6.536332E-9 * T - 1.120083E-6) * T + 1.001685E-4)
               * T - 9.09529E-3) * T + 6.793952E-2) * T - 28.263737
        R2 = (((5.3875E-9 * T - 8.2467E-7) * T + 7.6438E-5)
              * T - 4.0899E-3) * T + 8.24493E-1
        R3 = (-1.6546E-6 * T + 1.0227E-4) * T - 5.72466E-3
        Sig = (R4 * S + R3 * salSqRoot + R2) * S + R1
        V350p = 1 / R3500
        Sva = -Sig * V350p / (R3500 + Sig)
        Sigma = Sig + Dr350
        # Scale specific volume anomaly to normally reported units
        Svan = Sva * 1.0E+8
        if P == 0:
            return Sigma, Svan
        E = (9.1697E-10 * T + 2.0816E-8) * T - 9.9348E-7
        Bw = (5.2787E-8 * T - 6.12293E-6) * T + 3.47718E-5
        B = Bw + E * S
        D = 1.91075E-4
        C = (-1.6078E-6 * T - 1.0981E-5) * T + 2.2838E-3
        Aw = ((-5.77905E-7 * T + 1.16092E-4) * T + 1.43713E-3) * T - .1194975
        A = (D * salSqRoot + C) * S + Aw
        B1 = (-5.3009E-4 * T + 1.6483E-2) * T + 7.944E-2
        A1 = ((-6.167E-5 * T + 1.09987E-2) * T - .603459) * T + 54.6746
        Kw = (((-5.155288E-5 * T + 1.360477E-2) * T - 2.327105)
              * T + 148.4206) * T - 1930.06
        K0 = (B1 * salSqRoot + A1) * S + Kw
        Dk = (B * P + A) * P + K0
        K35 = (5.03217E-5 * P + 3.359406) * P + 21582.27
        Gam = P / K35
        Pk = 1.0 - Gam
        Sva = Sva * Pk + \
            (V350p + Sva) * P * Dk / (K35 * (K35 + Dk))
        Svan = Sva * 1.0E+8
        V350p = V350p * Pk
        # Density anomaly computed relative to 1000 kg/m^3
        # DR350 = density anomaly at 35 pss, 0 deg C and 0 decibars
        # dr35p = density anomaly at 35 pss, 0 deg C and pressure = p0 decibars
        # Dvan  = Density anomaly variations involving spec vol anom
        Dr35p = Gam / V350p
        Dvan = Sva / (V350p * (V350p + Sva))
        Sigma = Dr350 + Dr35p - Dvan
        return Sigma, Svan
    ''' Returns the spiciness of the water. '''
    def getSpiciness(self, temp, salt):
        # Hardcoded by Howard
        B = np.matrix([[0.0, .77442, -.00585, .000984, -.000206],
                       [.051665, .002034, -.0002745, -.0000085, .0000136],
                       [-6.64783E-3, -2.4681E-4, -1.428E-5, 3.337E-5, 7.894E-6],
                       [-5.4023E-5, 7.326E-6, 7.0036E-6, -3.0412E-6, -1.0853E-6],
                       [3.949E-7, -3.029E-8, -3.8209E-7, 1.0012E-7, 4.7133E-8],
                       [-6.36E-10, -1.309E-9, 6.048E-9, -1.1409E-9, -6.676E-10]])
        spice = 0
        sp = salt - 35
        theta = temp
        # Reversed I and J here, are arrays backwards in HP Basic?
        for I in range(0, 5):
            for J in range(0, 4):
                spice = spice + B[I, J] * (np.power(theta, I)) * (np.power(sp, J))
        return spice

    def sanityCheck(self):
        # ToDo
        return

    def checkFloatsFromIndex(self, floats, cycleJulDate):
        day, month, year = julianToDate(cycleJulDate)
        if self.verbose:
            print "Day, month, year are ", day, month, year
        self.yearMonthDayPath0 = (self.path0 
            + str(year) + '\\' 
            + str(month).zfill(2) + '\\' 
            + str(day).zfill(2) + '\\')
        inFileName = (self.yearMonthDayPath0 + str(year)
                      + str(month).zfill(2)
                      + str(day).zfill(2)
                      + '_index.csv')
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
                        floatNum = row[0]
                        floats.append(self.yearMonthDayPath0 + row[0])
                        self.passedFloat(floatNum, lat, lon)
                        self.numFloats += 1
        except (OSError, IOError), e:
            print "File was not found: ", inFileName
        return floats

    def passedFloat(self, floatNum, lat, lon):
        self.Lat[self.numFloats] = lat
        self.Lon[self.numFloats] = lon
        Dx = self.Scx0 * (float(lon) - self.longitudeDesired)
        Dy = self.Scy * (float(lat) - self.latitudeDesired)
        Rho = np.sqrt(Dx * Dx + Dy * Dy)
        if self.Closest > Rho:
            self.Closest = Rho
            self.Closest_fltnm = floatNum
        # ToDo: A few lines left out here, looks completely unnecessary
        return

    ''' TODO 
    Pulls settings from an already existing file. INCOMPLETE'''
    def setSettings():
        # Now using a .cfg file to load settings
        # settingsF = open((self.drive + "argo_programs\\TSlast.txt"), 'r+')
        return

    def defaultOptions(self):
        self.setEnabledParameters(False)
        return

    def setEnabledParameters(self, isEnabled):
        self.floatRadiusBox.setEnabled(isEnabled)
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
        self.writeDefinedSettings()
        return

    def writeDefinedSettings(self):
        # Save settings to .cfg file
        return

    ''' Opens and saves the output destination files in the class scope ''' 
    def prepareOutputFiles(self):
        tempPath = self.outPath + "TS_Te.csv"
        salinityPath = self.outPath + "TS_Sa.csv"
        sigmaTPath = self.outPath + "TS_St.csv"
        spicinessPath = self.outPath + "TS_Sp.csv"
        dynamicHeightPath = self.outPath + "TS_Dh.csv"
        hgtPath = self.outPath + "TS_Shgt.csv"     # Surface Height = Shgt
        cLimPath = self.sigPath + "Sigref.csv"
        stratPath = self.outPath + "Strat.csv"
        self.outputFilesLabel.setText("Output Files:\n" +
                                      tempPath + '\n' + 
                                      salinityPath + '\n' + 
                                      spicinessPath + '\n' + 
                                      dynamicHeightPath + '\n' + 
                                      hgtPath + '\n' + 
                                      stratPath + '\n' +
                                      cLimPath)
        self.outputLocationLabel.setText("Output Location:\n" + 
                                         self.outPath)
        writeType = 'a'
        if not self.append:
            writeType = 'w'
        self.tempCSV = open((tempPath), writeType)
        self.salinityCSV = open((salinityPath), writeType)
        self.sigmaTCSV = open((sigmaTPath), writeType)
        self.spicinessCSV = open((spicinessPath), writeType)
        self.dynamicHeightCSV = open((dynamicHeightPath), writeType)
        self.hgtCSV = open((hgtPath), writeType)
        self.cLimCSV = open((cLimPath), writeType)
        self.stratCSV = open((stratPath), writeType)
        if SAVELOCALLY:
            self.filesInUse = open((self.outPath + "IOS_Files_In_Use.txt"), 'w+')
        return

    def todayInDate(self):
        return self.curDay, self.curMonth, self.curYear

    def todayInJulian(self):
        day = int(time.strftime("%d"))
        month = int(time.strftime("%m"))
        year = int(time.strftime("%Y"))
        dayOfYear = int(time.strftime("%j"))
        todayInJul = dateToJulian(day, month, year, dayOfYear)
        return todayInJul

    ''' Lastrun (I think) defines if the program was a success or not. 
    Starts with 2 being written, then writes -1 after "Normal End"-ing'''

    def lastRun(self, curDrive, status):
        # Overwrite the previous "lastrun.txt" file
        lastRunF = open((curDrive + r"projects\lastrun.txt"), 'w+')
        # Typecast status to a string because write() uses strings
        lastRunF.write(str(status))
        lastRunF.close()
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
