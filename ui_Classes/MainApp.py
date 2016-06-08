# This is to deal with path issues for the sake of project organizaiton
import sys
sys.path.append("..")

# UI Imports
from PySide import QtCore, QtGui
from PySide.QtCore import QSize, QDate
from PySide.QtGui import QMainWindow
from ui_Files.ui_MainApp import Ui_MainApp
# Calculations imports
import numpy as np
import os
import datetime
import time
import csv


class MainApp(QMainWindow, Ui_MainApp):

    def __init__(self):
        super(MainApp, self).__init__()
        # self.ui = Ui_MainWindow()
        self.setupUi(self)
        self.setupSignals()
        # Declare input parameter variables
        self.defaultSettings = False
        self.xCoord = 0
        self.floatRadius = 0
        self.startDate = self.startRangeDateEdit.date().getDate()
        self.endDate = self.endRangeDateEdit.date().getDate()
        self.firstDayTuple = self.formatToDateTime(self.startDate[0],
                                                   self.startDate[1],
                                                   self.startDate[2])
        self.lastDayTuple = self.formatToDateTime(self.endDate[0],
                                                  self.endDate[1], 
                                                  self.endDate[2])
        self.dayStepSize = self.dayStepSizeBox.value()
        self.sampleWindow = 0
        self.firstLatitude = self.firstLatitudeBox.value()
        self.secondLatitude = self.secondLatitudeBox.value()
        self.firstLongitude = self.firstLongitudeBox.value()
        self.secondLongitude = self.secondLongitudeBox.value()
        self.pressureCutOff = self.pressureCutOffBox.value()
        self.maxInterpDepth = self.maxInterpDepthBox.value()
        self.stepSize = self.stepSizeBox.value()
        self.nPress = 1 + int(0.001 + self.maxInterpDepth / self.stepSize)
        if self.nPress > 300:
            self.nPress = 300

        self.temp = self.tempCheckBox.isChecked()
        self.salinity = self.salinityCheckBox.isChecked()
        self.sigmaT = self.sigmaTCheckBox.isChecked()
        self.spiciness = self.spicinessCheckBox.isChecked()
        self.dynamicHeight = self.dynamicHeightCheckBox.isChecked()

        self.latitudeDesired = self.latitudeDesiredBox.value()
        self.longitudeDesired = self.longitudeDesiredBox.value()
        self.append = True
        self.verbose = False
        self.outPath = None
        self.drive = ''
        self.path0 = ''
        self.pathTemp = ''
        self.sigPath = "projects\\Sigma_Climate"
        
        self.weight = 0
        self.tempCSV = None
        self.salinityCSV = None
        self.sigmaTCSV = None
        self.spicinessCSV = None
        self.dynamicHeightCSV = None
        self.hgtCSV = None
        self.cLimCSV = None

        self.yearMonthDayPath0 = ''
        self.Scy=111.2  # ! Scy = km/degree of latitude
        self.Scx0 = self.Scy*np.cos(self.latitudeDesired)
        self.Rho0=300  # ! Decay scale of weighting function
        self.Pref=1000. # ! Reference level for dynamic height calculations
        # Actual code, for testing easier though
        # self.curDay = int(time.strftime("%d"))
        # self.curMonth = int(time.strftime("%m"))
        # self.curYear = int(time.strftime("%Y"))
        self.curDay = 18
        self.curMonth = 05
        self.curYear = 2016
        # Junk to make testing easier ^
        self.numFloats = 0
        self.userDefinedSettings()
        self.initVars()
        return

    def setupSignals(self):
        self.exitButton.clicked.connect(self.close)
        self.timeSeriesButton.clicked.connect(self.timeSeriesButtonClicked)
        # Declare input parameter variables
        self.setupInputParameterSignals()
        # GUI related to TimeSeries
        self.backToProgramListButton.clicked.connect(
            self.backToProgramListButtonClicked)
        self.nextButton.clicked.connect(self.nextButtonClicked)
        return

    def setupInputParameterSignals(self):
        self.defaultSettingsCheckBox.stateChanged.connect(
            self.defaultSettingsCheckBoxStateChanged)
        self.floatRadiusBox.editingFinished.connect(
            self.floatRadiusBoxEditingFinished)
        # self.dayRangeBox.editingFinished.connect(
        #         self.dayRangeBoxEditingFinished)
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


    def formatToDateTime(self, year, month, day):
        stringVersion = str(day) + '.' + str(month) + '.' + str(year)
        stringFormat = '%d.%m.%Y'
        dateTuple = datetime.datetime.strptime(stringVersion,
                                               stringFormat)
        return dateTuple

    # Start of parameter boxes
    def startDateChanged(self, date):
        self.startDate = date.getDate()
        self.firstDayTuple = self.formatToDateTime(self.startDate[0],
                                                   self.startDate[1],
                                                   self.startDate[2])
        return

    def endDateChanged(self, date):
        self.endDate = date.getDate()
        self.lastDayTuple = self.formatToDateTime(self.endDate[0],
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
        self.dayStepSizeBox = self.dayStepSizeBox.value()
        return

    def sampleWindowBoxEditingFinished(self):
        self.sampleWindowBox = self.sampleWindowBox.value()
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
        return

    def stepSizeBoxEditingFinished(self):
        self.stepSize = self.stepSizeBox.value()
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
    # End of parameter boxes

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
        self.yearMonthDayPath0 = self.path0 \
            + str(self.curYear) + '\\' \
            + str(self.curMonth).zfill(2) + '\\' \
            + str(self.curDay).zfill(2) + '\\'
        self.timeSeriesStackedWidget.setCurrentWidget(self.calculationsPage)
        self.prepareOutputFiles()
        self.commenceInterpolation()
        return

    def sanityCheck(self):

        return

    def getSigmaT(self, S, T):
        # 12620 SUB Sigma_t(S,T,Sigma_t)
        # 12630!
        # 12640! T  = Temperature in C
        # 12650! S  = salinity in pss-78
        # 12660!
        R4 = 4.8314E-4
        Dr350 = 28.106331
        Sr = np.sqrt(S)
        R1 = ((((6.536332E-9 * T - 1.120083E-6) * T + 1.001685E-4)
               * T - 9.09529E-3) * T + 6.793952E-2) * T - 28.263737
        R2 = (((5.3875E-9 * T - 8.2467E-7) * T + 7.6438E-5)
              * T - 4.0899E-3) * T + 8.24493E-1
        R3 = (-1.6546E-6 * T + 1.0227E-4) * T - 5.72466E-3
        Sig = (R4 * S + R3 * Sr + R2) * S + R1
        Sigma = Sig + Dr350
        Sigma_t = ((Sigma + 1000.) / .999975) - 1000.
        return Sigma_t

    # Line 1402-1623
    def getProfile(self,
                   floatNum,
                   pressureMatrix,
                   tempMatrix,
                   salinityMatrix,
                   rFlag):
        # Setup path (skip)
        rFlag = -1
        order = []
        path = self.yearMonthDayPath0 + floatNum
        numChannels = 0
        numRecords = 0
        fileTempQc = 0
        filePSalQc = 0
        numRecs = self.extractDataFromFloatFile(path, order)
        # He opens this path?
        # Line1 -> 16210 close file, stop logging, lastRun
        if not self.checkIfReturn(numRecs, fileTempQc, filePSalQc):
            print "in getProfile, returning prematurely due to bad values"
            return
        # self.channelOrder()
        return numRecs

    def extractDataFromFloatFile(self, path, order):
        print "Opening path ", path
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
        print "Order is ", order
        return numRecs 

    def getNumRecs(self, line):
        split = line.split()
        # print "numRecs is ", int(split[4]) 
        return int(split[4])         

    def appendInfo(self, line, order, i):
        split = line.split()
        # print "Split is ", split
        offset = 5
        position = 0
        for data in order:
            twoPos = float(split[2 + (position * 5)])
            zeroPos = float(split[0 + (position * 5)])
            if data == 'P':
                if np.abs(twoPos) < 9000:
                    # print "twoPos is ", twoPos
                    # print "i is ", i
                    self.P[i] = twoPos
                else:
                    # print "zero pso is ", zeroPos
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

    def verifyUsability(self, line):
        if line.find("NUMBER OF RECORDS") != -1:
            if line.split()[2] < 15:
                return False
        if line.find("NUMBER OF CHANNELS") != -1:
            if line.split()[2] < 8:
                return False
        return True

    # TODO: Improve runtime, reduce calls to this function
    def findOrder(self, line, order, pFound, tFound, sFound):
        resP = -1
        resT = -1
        resS = -1
        print line
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

    def commenceInterpolation(self):
        yearDayStart = self.firstDayTuple.timetuple().tm_yday
        yearDayEnd = self.lastDayTuple.timetuple().tm_yday
        julStart = self.dateToJulian(self.startDate[2],
                                     self.startDate[1],
                                     self.startDate[0],
                                     yearDayStart)
        julEnd = self.dateToJulian(self.endDate[2],
                                   self.endDate[1],
                                   self.endDate[0],
                                   yearDayEnd)
        # print "Starting date ", self.startDate[2], " Ending date ", self.endDate[2]
        # print "start Date is ", julStart, " end Date is ", julEnd
        for Dc in xrange(julStart, julEnd, self.dayStepSize):
            # print some stuff for user
            # print "INSIDE LOOP"
            self.xCoord = Dc
            sTav = 0
            sTavSq = 0
            sTKnt = 0
            floats = self.checkFloatsFromIndex(self.yearMonthDayPath0)
            # print "floats are ", floats
            rFlag = -1
            i = 0
            for singleFloat in floats:
                # print "there's a float for all the floats!"
                numRecs = self.getProfile(singleFloat, self.P, self.T, self.S, rFlag)
                self.sanityCheck()
                sigT = self.getSigmaT(self.S[1], self.T[1])
                sTav, sTavSq, sTKnt = self.randomCalcs(sigT, sTav, sTavSq, sTKnt, numRecs, i)
                i += 1
                break               # ToDo: This is for testing
            # secondRun has double for loop
            Ipr75 = self.secondRun(sTav, sTavSq, sTKnt)
            self.thirdRun()
            self.computeDHgt()
            self.finishUp(Ipr75)
        self.cleanUp()
        return

    def cleanUp(self):
        # IF Te_opt>.5 THEN ASSIGN @Pte TO *
        # IF Sa_opt>.5 THEN ASSIGN @Psa TO *
        # IF St_opt>.5 THEN ASSIGN @Pst TO *
        # IF Sp_opt>.5 THEN ASSIGN @Psp TO *
        # IF Dh_opt>.5 THEN ASSIGN @Pdh TO *
        # IF Dh_opt>.5 THEN ASSIGN @Psu TO *
        # ASSIGN @Pstrat TO *
        self.tempCSV.close()
        self.salinityCSV.close()
        self.sigmaTCSV.close()
        self.spicinessCSV.close()
        self.dynamicHeightCSV.close()
        self.hgtCSV.close()
        self.cLimCSV.close()
        self.stratCSV.close()

        print "-----------------------------------------------------------------------"
        # 
        # DISP
        print "Interpolations are completed and stored."
        # !
        # CALL Logit(2,"TimeSeries.prg",Drive$)
        # 
        # ON KEY 5 LABEL "Select  Programs" GOTO 10320
        # ON KEY 6 LABEL "Quit" GOTO 10370
        # 
        # GOTO 10300
        # 
        # GCLEAR
        # CLEAR SCREEN
        # 
        # LOAD "d:\projects\Argo\Selector.prg",1
        # 
        # PRINT "Normal End"
        self.Lastrun(self.drive,-1)
        # 
        return

    def finishUp(self, Ipr75):

        # S$=Xcoord$&","&VAL$(Dh(1))&","&Wgt$
        # IF Dh_opt>.5 THEN OUTPUT @Psu;S$
        if self.dynamicHeight:
            self.dynamicHeightCSV.write(str(self.xCoord) + ',' + 
                str(self.Dh[1]) + ',' + str(self.weight))
        # !
        Stdiff=self.St[Ipr75]-self.St[1]
        if Stdiff<-.002:
            Stdiff=0.
        self.stratCSV.write(str(self.xCoord) + ',' + str(self.St[1]) + ',' 
            + str(self.St[Ipr75]) + ',' + str(self.weight))
        # !
        Stdiff=.001*int(1000.*Stdiff+.5)
        # IF Stdiff<.9999999 THEN Stdiff$="0"&Stdiff$
        if Stdiff<.001:
            Stdiff=0.000
        St1 = .001*int(.5+1000.*self.St[1])
        Dayc, Monc, Yearc = self.julianToDate(self.xCoord)
        R = "| " + str(self.xCoord) + " " + str(Dayc) + "/" + str(Monc) + "/" + \
            str(Yearc) + " | " + str(self.weight) + "  " + str(Stdiff) + "    " + \
            str(St1)
        self.Closest=.1*int(10.*self.Closest+.5)
        # Pp=POS(Closest_fltnm$,".")
        R += R + " | " + str(self.Closest) + "  " + \
            self.Closest_fltnm + " |"
        print R

        # ToDo: Optimize, open outside of loop
        lstMsge = open((self.outPath + "lstmsge.csv"), 'w+')
        lstMsge.truncate()
        lstMsge.write(R)
        lstMsge.close()        
        # IF Verb_opt>.5 THEN PRINT "Done ";Dc1;Dc;Dc2
        # NEXT Dc
        return

    def randomCalcs(self, sigT, sTav, sTavSq, sTKnt, numRecs, Iflt):
        # Copied, not sure why Howard does this
        sigT = ((.5 + 1000. * sigT)) / 1000
        # print "sigT is ", sigT
        if sigT > 5 and sigT < 30:
            sTav += sigT
            sTavSq += np.power(sigT, 2)
            sTKnt += 1
        # print "P[1] is ", self.P[1]
        if self.P[1] < 20:
            self.P[1] = 0
        # Re-intitializes Te[x][y] and Sa to 999.9
        # print "nPress is ", self.nPress
        for Ipr in range(1, self.nPress):
            pressureCalc = self.stepSize * (Ipr - 1.)  # Interpolate to pressureCalc
            # print "pressureCalc is ", pressureCalc, " P[1] is ", self.P[1]
            if (pressureCalc < self.P[1]):
                break
            # print "numRecs is ", numRecs
            # print "self.P[numRecs] is ", self.P[numRecs]
            # print "P is ", self.P
            if (pressureCalc > self.P[numRecs]):
                break
            # ! Somewhere in P(Ipr) there is data at pressure pressureCalc
            Eps = 1.0E-6
            # print "numRecs is ", numRecs
            # print "nPress is ", self.nPress
            for K in range(2, numRecs):
                if pressureCalc >= self.P[K - 1] and pressureCalc <= self.P[K]:
                    Del_p = self.P[K] - self.P[K - 1]
                    # print "----- Del_p is ", Del_p 
                    if Del_p < 100:
                        Rho = (pressureCalc - self.P[K - 1]) / \
                            (self.P[K] - self.P[K - 1] + Eps)
                        self.Te[Iflt, Ipr] = self.T[K - 1] + \
                            Rho * (self.T[K] - self.T[K - 1] + Eps)
                        # print "Setting Te to ", self.T[K - 1] + \
                        #     Rho * (self.T[K] - self.T[K - 1] + Eps), " at ", Iflt, Ipr 
                        self.Sa[Iflt, Ipr] = self.S[K - 1] + \
                            Rho * (self.S[K] - self.S[K - 1] + Eps)
                    break
            # ! Finished interpolation to standard pressures for Float Iflt
        # Finished interpolation to standard pressures for all floats
        return sTav, sTavSq, sTKnt

    def secondRun(self, sTav, sTavSq, sTKnt):
        sTav = sTav / sTKnt
        sTavSq = sTavSq / sTKnt
        stdDev = np.power((sTavSq - sTav * sTav), 2)
        # if Verb_opt > .5:
            # print "Mean & std dev of St = ", sTav, stdDev
        # ! Interpolate to station
        # MAT Pdel=(-9999.9)
        # MAT Sva=(-9999.9)
        # MAT Dh=(-9999.9)
        Kdel = 0
        Qsa = 0
        Qte = 0
        Qsp = 0
        Ipr75 = int(1.1 + 75. / self.stepSize)
        weightSumTMax = 0.
        for iterPress in range(1, self.nPress):
            pressureCal = self.stepSize * (iterPress - 1)  # What does this line do
            # 8330 REM IF Pc>Wdep THEN GOTO 5000
            weightSumS = 0.
            weightSumT = 0.
            salWeightSumT = 0.
            tempWeightSumT = 0.
            # print "numFloats is ", self.numFloats
            for iterFloat in range(1, self.numFloats):
                # IF Accpt(Iflt)<0. THEN GOTO 8610
                # Latav = latitude average?
                Latav = (self.Lat[iterFloat] + self.latitudeDesired) / 2
                Scx = self.Scy * np.cos(Latav)
                Dx = Scx * (self.longitudeDesired - self.Lon[iterFloat])
                # print "Dx is ", Dx
                Dy = self.Scy * (self.latitudeDesired - self.Lat[iterFloat])
                # print "Dy is ", Dy
                Rho = np.sqrt(Dx * Dx + Dy * Dy)
                # print "Rho is ", Rho
                Z = Rho / self.Rho0
                # print "Z is ", Z
                # IF Z>5. THEN GOTO 8610
                if not (Z > 5):
                    # print "THERE IS AN OCCURANCE OF Z NOT > 5"
                    self.weight = np.exp(-Z * Z)
                    # !
                    # print "Te's are ", self.Te[iterFloat, iterPress], self.Te[iterFloat, iterPress]
                    # print "inside if at ", iterFloat, iterPress
                    if (self.Te[iterFloat, iterPress] > -1.5 and \
                            self.Te[iterFloat, iterPress] < 40.):
                        weightSumT = weightSumT + weight
                        tempWeightSumT = tempWeightSumT + \
                            weight * self.Te[iterFloat, iterPress]

                    if self.Sa[iterFloat, iterPress] > 20. and \
                            self.Sa[iterFloat, iterPress] < 40.:
                        weightSumS = weightSumS + weight
                        salWeightSumT = salWeightSumT + \
                            weight * self.Sa[iterFloat, iterPress]
            if weightSumT > weightSumTMax:
                weightSumTMax = weightSumT
            if weightSumT > 0:
                Qte = tempWeightSumT / weightSumT
            else:
                print "no floats usable to calculate Qte"
            if weightSumS > 0:
                Qsa = salWeightSumT / weightSumS
            else:
                print "no floats usable to calculate Qsa"

            Qst = self.getSigmaT(Qsa, Qte)
            Spiciness = self.getSpiciness(Qsp, Qte, Qsa)
            Sigma, Svan = self.getSvanom(Qsa, Qte, 0)
            self.P[iterPress] = self.stepSize * (iterPress - 1.)
            # print "things are ", self.stepSize * (iterPress - 1.)
            self.T[iterPress] = Qte
            self.S[iterPress] = Qsa
            self.St[iterPress] = Qst
            self.Sp[iterPress] = Qsp
            Kdel = Kdel + 1
            self.Pdel[Kdel] = pressureCal
            self.Sva[Kdel] = Svan
        return Ipr75

    def thirdRun(self):
        # 8880 Some kind of formatting?
        # ----- ^ ALL UNNECESSARY? -----

        if self.St[1]>self.St[2] and np.abs(self.St[1]-self.St[2])>.05:
            self.St[1]=self.St[2]
            self.T[1]=self.T[2]
            self.S[1]=self.S[2]
            self.Sp[1]=self.Sp[2]

        # # ! Now compute DHgt relative to Pmax
        # Dh(Npress)=0.
        # Q=5.6E-6 ! Q=0.5f/g at station Papa
        return


    def computeDHgt(self):
        
        # Now compute DHgt relative to Pmax
        # Dh(Npress)=0.
        Q=5.6E-6 # ! Q=0.5f/g at station Papa
        for K in range(2, self.nPress):
            # Ipr=Npress+1-K
            # Dh(Ipr)=Dh(Ipr+1)+Q*(Sva(Ipr+1)+Sva(Ipr))*Deltap
            # NEXT K
            pass
        
        for Ipr in range(1, self.nPress):
            Qte=.001*int(.5+1000.*self.T[Ipr])
            Qsa=.001*int(.5+1000.*self.S[Ipr])
            Qst=.001*int(.5+1000.*self.St[Ipr])
            Qsp=.001*int(.5+1000.*self.Sp[Ipr])
            Qdh=.001*int(.5+1000.*self.Dh[Ipr])
            Pc=self.stepSize*(Ipr-1.)
            
            Sigref_opt=-1
            # print "Qst is ", Qst
            # IF Qst<Sigref_st(1) THEN GOTO 9420
            # IF Qst>Sigref_st(71) THEN GOTO 9420

            # Updated indexes from HP references 1 -> 0, 71 -> 70
            if not (Qst < self.Sigref_st[0]) and not (Qst > self.Sigref_st[70]): 
                # for Icl=1 TO 70
                #     IF Qst>=Sigref_st(Icl) AND Qst<=Sigref_st(Icl+1) THEN GOTO 9300
                #     NEXT Icl
                for Icl in range(0, 69):
                    if Qst >= self.Sigref_st[Icl] and Qst <= self.Sigref_st[Icl+1]:
                        tep, sap = self.foundPair()
                        break
            
            # S$=Xcoord$&","&VAL$(-Pc)&","&VAL$(Qte)
            xCoAndPc = str(self.xCoord) + ',' + str(-Pc) + ','
            if self.temp:
                self.tempCSV.write(xCoAndPc + str(Qte))
            if self.salinity:
                self.salinityCSV.write(xCoAndPc + str(Qsa))
            if self.sigmaT:
                self.sigmaTCSV.write(xCoAndPc + str(Qst))
            if self.spiciness:
                self.spicinessCSV.write(xCoAndPc + str(Qsp))
            if self.dynamicHeight:
                self.dynamicHeightCSV.write(xCoAndPc + str(Qdh))
        return

    def foundPair(self, Qst, Icl):
        # found a pair
        Ra=(Qst-self.Sigref_st[Icl])/(self.Sigref_st[Icl+1]-self.Sigref_st[Icl])
        Ter=self.Sigref_te[Icl]+Ra*(self.Sigref_te[Icl+1]-self.Sigref_te[Icl])
        Sar=self.Sigref_sa[Icl]+Ra*(self.Sigref_sa[Icl+1]-self.Sigref_sa[Icl])
        Tep=Qte-Ter
        Sap=Qsa-Sar
        Tep=.001*int(.5+1000.*Tep)
        Sap=.001*int(.5+1000.*Sap)
        # S$=Xcoord$&","&VAL$(Qst)&","&VAL$(Qte)&","&VAL$(Tep)&","&VAL$(Qsa)&","&VAL$(Sap)
        # OUTPUT @Psclim;S$
        print "Qst: ", Qst, " Qte: ", Qte, " Tep: ", Tep, " Qsa: ", " Sap: ", Sap
        return Tep, Sap

    def getSvanom(self, S, T, P0):
        #! Compute the density anomaly, sigma, in kg/m^3
        #! Density anomaly is identical with sigma-t without pressure terms
        #!
        #! P0 = Pressure in decibars
        #! T  = Temperature in deg C
        #! S  = salinity in pss-78
        #!
        R3500 = 1028.106331
        R4 = 4.8314E-4
        Dr350 = 28.106331
        P = P0 / 10
        Sr = np.sqrt(S)
        R1 = ((((6.536332E-9 * T - 1.120083E-6) * T + 1.001685E-4)
               * T - 9.09529E-3) * T + 6.793952E-2) * T - 28.263737
        R2 = (((5.3875E-9 * T - 8.2467E-7) * T + 7.6438E-5)
              * T - 4.0899E-3) * T + 8.24493E-1
        R3 = (-1.6546E-6 * T + 1.0227E-4) * T - 5.72466E-3
        Sig = (R4 * S + R3 * Sr + R2) * S + R1
        V350p = 1 / R3500
        Sva = -Sig * V350p / (R3500 + Sig)
        Sigma = Sig + Dr350
        #! Scale specific volume anomaly to normally reported units
        Svan = Sva * 1.0E+8
        if P == 0:
                    # THEN GOTO 13710
            return 0, 0 # ToDo: Might not be correct
        E = (9.1697E-10 * T + 2.0816E-8) * T - 9.9348E-7
        Bw = (5.2787E-8 * T - 6.12293E-6) * T + 3.47718E-5
        B = Bw + E * S
        D = 1.91075E-4
        C = (-1.6078E-6 * T - 1.0981E-5) * T + 2.2838E-3
        Aw = ((-5.77905E-7 * T + 1.16092E-4) * T + 1.43713E-3) * T - .1194975
        A = (D * Sr + C) * S + Aw
        B1 = (-5.3009E-4 * T + 1.6483E-2) * T + 7.944E-2
        A1 = ((-6.167E-5 * T + 1.09987E-2) * T - .603459) * T + 54.6746
        Kw = (((-5.155288E-5 * T + 1.360477E-2) * T - 2.327105)
              * T + 148.4206) * T - 1930.06
        K0 = (B1 * Sr + A1) * S + Kw
        Dk = (B * P + A) * P + K0
        K35 = (5.03217E-5 * P + 3.359406) * P + 21582.27
        Gam = P / K35
        Pk = 1.0 - Gam
        Sva = Sva * Pk + (V350p + Sva) * P * Dk / (K35 * (K35 + Dk))
        Svan = Sva * 1.0E+8
        V350p = V350p * Pk
        # Density anomaly computed relative to 1000 kg/m^3
        # DR350 = density anomaly at 35 pss, 0 deg C and 0 decibars
        # dr35p = density anomaly at 35 pss, 0 deg C and pressure = p0 decibars
        # Dvan  = Density anomaly variations involving spec vol anom
        Dr35p = Gam / V350p
        Dvan = Sva / (V350p * (V350p + Sva))
        Sigma = Dr350 + Dr35p - Dvan
        # SUBEND
        return Sigma, Dvan

    def getSpiciness(self, spice, temp, salt):
        # Hardcoded by Howard
        B = np.matrix([[0.0, .77442, -.00585, .000984, -.000206],
                     [.051665, .002034, -.0002745, -.0000085, .0000136],
                     [-6.64783E-3, -2.4681E-4, -1.428E-5, 3.337E-5, 7.894E-6],
                     [-5.4023E-5, 7.326E-6, 7.0036E-6, -3.0412E-6, -1.0853E-6],
                     [3.949E-7, -3.029E-8, -3.8209E-7, 1.0012E-7, 4.7133E-8],
                     [-6.36E-10, -1.309E-9, 6.048E-9, -1.1409E-9, -6.676E-10]])
        spice = 0
        Sp = salt - 35
        Theta = temp
        # Reversed I and J here, are arrays backwards in HP Basic?
        for I in range(0, 5):
            Ii = I - 1
            for J in range(0, 4):
                Jj = J - 1
                # In Python, ^ is XOR 8 ^ 3 => 1000 ^ 0011 => 1011 = 11
                spice = spice + B[I, J] * (Theta ^ I) * (Sp ^ J)
        return spice

    def sanityCheck(self):
        # ToDo
        return

    def checkFloatsFromIndex(self, folder):
        with open((folder + str(self.curYear)
                   + str(self.curMonth).zfill(2)
                   + str(self.curDay).zfill(2)
                   + '_index.csv'),
                  'rb') as indexCSV:
            print "Opened index"
            reader = csv.reader(indexCSV)
            buoysToUse = []
            print self.firstLatitude, self.secondLatitude, self.firstLongitude, self.secondLongitude
            # [1:] means skip the first. Skip it because it's a header
            reader.next()
            for row in reader:
                if float(row[3]) < self.pressureCutOff:
                    print "Long ", float(row[2])
                    print "Lat ", float(row[1])
                
                # ToDo: Does not account for users with lat2 < lat1
                if (float(row[1]) > self.firstLatitude
                        and float(row[1]) < self.secondLatitude
                        and float(row[2]) > self.firstLongitude
                        and float(row[2]) < self.secondLongitude
                        and float(row[3]) > self.pressureCutOff):
                    # ToDo: FOR TESTING
                    buoysToUse.append("20160518_2901481.IOS")
                    # ToDo: Reactivate for real
                    # buoysToUse.append(row[0]) 
                    self.passedFloat(row)
                    self.numFloats += 1
                    break # ToDo: Only active for testing
                    # Save it
                # print row[0]
        print "Bouys to use is ", buoysToUse
        return buoysToUse

    def passedFloat(self, row):
        self.Lat[self.numFloats] = float(row[1])
        self.Lon[self.numFloats] = float(row[2])
        Dx=self.Scx0*(float(row[2])- self.longitudeDesired)
        Dy=self.Scy*(float(row[1])- self.latitudeDesired)
        Rho=np.sqrt(Dx*Dx+Dy*Dy)
        if self.Closest > Rho:
            self.Closest = Rho
            self.Closest_fltnm = row[0]
        # ToDo: A lot left out here, maybe unecessary? 
        return

    # Line 0-45
    def initVars(self):
        # 10 OPTION BASE 1
        # Te/Sa are temperature and salinity?
        self.Te = np.empty((500, 300))
        self.Sa = np.empty((500, 300))
        self.Lat = np.empty((500))
        self.Lon = np.empty((500))
        self.Accpt = np.empty((500))
        self.Fltnm = np.empty((500))

        self.Pdel = np.empty((1000))
        self.Sva = np.empty((1000))
        self.Dh = np.empty((1000))

        Sa_av = np.empty(12)
        Sa_knt = np.empty(12)

        self.Sigref_st = np.empty((71))
        self.Sigref_te = np.empty((71))
        self.Sigref_sa = np.empty((71))
        # 110! Dimensions of data to be read from files
        self.P = np.zeros((2000))
        self.T = np.empty((2000))
        self.S = np.empty((2000))
        self.St = np.empty((2000))
        self.Sp = np.empty((2000))
        # 140 ! ALPHA PEN 1
        # 150 Tabrow=11
        self.Closest = 9999.9
        # 240 ! Initialize data matrices with nul values 999.9
        # Should I do the NumPy NULL or just stick with 999.9?
        self.Te[:] = 999.9
        self.Sa[:] = 999.9
        self.Lat[:] = 999.9
        self.Lon[:] = 999.9
        self.Accpt[:] = (-1)
        Sa_av[:] = 0
        Sa_knt[:] = 0
        # Dout backslash because in Python you need to escape special
        # chearacters
        self.drive = "E:\\"
        # 350 CALL Lastrun(self.drive$,2)
        self.lastRun(self.drive, 2)
        self.path0 = self.drive + "argo_mirror\\pacific_ocean\\"
        self.outPath = self.drive + "argo_out_TEST\\TimeSeries\\"
        self.sigPath = self.drive + "projects\\Sigma_Climate\\"
        # Made this file empty here when he writes 3x71 empties
        # if problems: change this maybe
        csvF = open((self.sigPath + "Mp26_i.csv"), 'w')
        csvF.close()
        return

    ''' TODO 
    Pulls settings from an already existing file. INCOMPLETE'''
    def setSettings():
        # 470 ! Read previous settings
        # 480 ASSIGN @Pin TO Drive$&"argo_programs\TSlast.txt";FORMAT ON
        # ToDo: Change to 'a'
        settingsF = open((self.drive + "argo_programs\\TSlast.txt"), 'r+')
        for x in xrange(1, 20):
            pass
        # 490 ON END @Pin GOTO 660
        # 500 FOR L=1 TO 20
        # 510 ENTER @Pin;S$
        # 520 ! PRINT L;S$
        # 530 Pc=POS(S$,",")
        # 540 IF POS(S$,"G_opt")>.1 THEN G_opt$=TRIM$(S$[Pc+1])
        # 550 IF POS(S$,"c1dc")>.1 THEN Dc12$=TRIM$(S$[Pc+1])
        # 560 IF POS(S$,"Dt")>.1 THEN Dt$=TRIM$(S$[Pc+1])
        # 570 IF POS(S$,"windo")>.1 THEN Twin$=TRIM$(S$[Pc+1])
        # 580 IF POS(S$,"at1&L")>.1 THEN Lat$=TRIM$(S$[Pc+1])
        # 590 IF POS(S$,"on1&L")>.1 THEN Lon$=TRIM$(S$[Pc+1])
        # 600 IF POS(S$,"cutof")>.1 THEN Pcut$=TRIM$(S$[Pc+1])
        # 610 IF POS(S$,"max&Del")>.1 THEN Pmax$=TRIM$(S$[Pc+1])
        # 620 IF POS(S$,"ocatio")>.1 THEN Location$=TRIM$(S$[Pc+1])
        # 630 IF POS(S$,"ameter")>.1 THEN Param$=TRIM$(S$[Pc+1])
        # 640 NEXT L
        # 650 PRINT "Param$ = ";Param$
        # 660 ASSIGN @Pin TO *
        # 670 !
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

    # Line 65-285
    def userDefinedSettings(self):
        self.setEnabledParameters(True)
        # If the user has NOT selected default settings
        todayInJul = self.todayInJulian()
        print todayInJul
        day, month, year = self.todayInDate()
        print day, month, year
        self.currentDayDateEdit.setDate(QDate(year, month, day))
        self.writeDefinedSettings()
        return

    def writeDefinedSettings(self):
        # 2900 !
        # 2910 ! Write these settings for next time
        # 2920 PURGE Drive$&"argo_programs\TSlast.txt"
        # 2930 CREATE Drive$&"argo_programs\TSlast.txt",1
        # 2940 ASSIGN @Pout TO Drive$&"argo_programs\TSlast.txt";FORMAT ON
        # 2950 OUTPUT @Pout;"G_opt,"&TRIM$(VAL$(G_opt))
        # 2960 OUTPUT @Pout;"Dc1dc2,"&TRIM$(VAL$(Dc1))&","&TRIM$(VAL$(Dc2))
        # 2970 OUTPUT @Pout;"Dt,"&TRIM$(VAL$(Dt))
        # 2980 OUTPUT @Pout;"Twindow,"&TRIM$(VAL$(Twin))
        # 2990 OUTPUT @Pout;"Lat1&Lat2,"&TRIM$(VAL$(Lat1))&","&TRIM$(VAL$(Lat2))
        # 3000 OUTPUT @Pout;"Lon1&Lon2,"&TRIM$(VAL$(Lon1))&","&TRIM$(VAL$(Lon2))
        # 3010 OUTPUT @Pout;"Pcutoff,"&TRIM$(VAL$(Pcutoff))
        # 3020 OUTPUT @Pout;"Pmax&Deltap,"&TRIM$(VAL$(Pmax))&","&TRIM$(VAL$(Deltap))
        # 3030 Param$=TRIM$(VAL$(Te_opt))&","&TRIM$(VAL$(Sa_opt))&","&TRIM$(VAL$(St_opt))
        # 3040 Param$=Param$&","&TRIM$(VAL$(Sp_opt))&","&TRIM$(VAL$(Dh_opt))
        # 3050 OUTPUT @Pout;"Parameters,"&Param$
        # 3060 OUTPUT @Pout;"Location,"&TRIM$(VAL$(Latc))&","&TRIM$(VAL$(Lonc))
        # 3070 ASSIGN @Pout TO *
        # 3080 !
        return

    def prepareOutputFiles(self):

        self.tempPath = self.outPath + "TS_Te.csv"
        self.salinityPath = self.outPath + "TS_Sa.csv"
        self.sigmaTPath = self.outPath + "TS_St.csv"
        self.spicinessPath = self.outPath + "TS_Sp.csv"
        self.dynamicHeightPath = self.outPath + "TS_Dh.csv"
        self.hgtPath = self.outPath + "TS_Shgt.csv"
        self.cLimPath = self.sigPath + "Sigref.csv" 
        self.stratPath = self.outPath + "Strat.csv"
        self.outputFilesLabel.setText("Output Files:\n"
                                      + self.tempPath + '\n'
                                      + self.salinityPath + '\n'
                                      + self.spicinessPath + '\n'
                                      + self.dynamicHeightPath + '\n'
                                      + self.hgtPath + '\n'
                                      + self.cLimPath + '\n'
                                      + self.stratPath)
        writeType = 'a'
        if not self.append:
            writeType = 'w'
        self.tempCSV = open((self.tempPath), writeType)
        self.salinityCSV = open((self.salinityPath), writeType)
        self.sigmaTCSV = open((self.sigmaTPath), writeType)
        self.spicinessCSV = open((self.spicinessPath), writeType)
        self.dynamicHeightCSV = open((self.dynamicHeightPath), writeType)
        self.hgtCSV = open((self.hgtPath), writeType)
        self.cLimCSV = open((self.cLimPath), writeType)
        self.stratCSV = open((self.stratPath), writeType)
        # What the heck is hgt? Mercury...t?

        # File on sigma-levels
        if (self.append):
            self.prepareFilesAndAppend()
        else:
            self.prepareFilesAndClear()
            # Maybe not necessary?
        # Display "Starting interpolation value 1 to value 2"
        return

    def prepareFilesAndAppend(self):
        # 3850 !
        # ALL OF THIS IS JUST TO SET THE FILES UP FOR APPENDING -_______-
        # 5840 !
        # 5850 ! End of yes-append option
        # 5860 END IF
        # 5870 !
        return

    def prepareFilesAndClear(self):
        # 3230 !
        # 3240 IF Append_opt<.5 THEN
        # 3250 ! No Append
        # 3260 ON ERROR GOTO 3280
        # 3270 PURGE Sclim$
        # 3280 OFF ERROR
        # 3290 CREATE Sclim$,1
        # 3300 ASSIGN @Psclim TO Sclim$;FORMAT ON
        # 3310 !
        # 3320 IF Te_opt>.5 THEN
        # 3330 ON ERROR GOTO 3350
        # 3340 PURGE Ste$
        # 3350 OFF ERROR
        # 3360 CREATE Ste$,1
        # 3370 ASSIGN @Pte TO Ste$;FORMAT ON
        # 3380 END IF
        # 3390 !
        # 3400 IF Sa_opt>.5 THEN
        # 3410 ON ERROR GOTO 3430
        # 3420 PURGE Ssa$
        # 3430 OFF ERROR
        # 3440 CREATE Ssa$,1
        # 3450 ASSIGN @Psa TO Ssa$;FORMAT ON
        # 3460 END IF
        # 3470 !
        # 3480 IF St_opt>.5 THEN
        # 3490 ON ERROR GOTO 3510
        # 3500 PURGE Sst$
        # 3510 OFF ERROR
        # 3520 CREATE Sst$,1
        # 3530 ASSIGN @Pst TO Sst$;FORMAT ON
        # 3540 END IF
        # 3550 !
        # 3560 IF Sp_opt>.5 THEN
        # 3570 ON ERROR GOTO 3590
        # 3580 PURGE Ssp$
        # 3590 OFF ERROR
        # 3600 CREATE Ssp$,1
        # 3610 ASSIGN @Psp TO Ssp$;FORMAT ON
        # 3620 END IF
        # 3630 !
        # 3640 IF Dh_opt>.5 THEN
        # 3650 ON ERROR GOTO 3670
        # 3660 PURGE Sdh$
        # 3670 OFF ERROR
        # 3680 CREATE Sdh$,1
        # 3690 ASSIGN @Pdh TO Sdh$;FORMAT ON
        # 3700 !
        # 3710 ON ERROR GOTO 3730
        # 3720 PURGE Shgt$
        # 3730 OFF ERROR
        # 3740 CREATE Shgt$,1
        # 3750 ASSIGN @Psu TO Shgt$;FORMAT ON
        # 3760 END IF
        # 3770 !
        # 3780 ON ERROR GOTO 3800
        # 3790 PURGE Pathout$&"Strat.csv"
        # 3800 OFF ERROR
        # 3810 CREATE Pathout$&"Strat.csv",1
        # 3820 ASSIGN @Pstrat TO Pathout$&"Strat.csv";FORMAT ON
        # 3830 ! End of no-append option
        # Pte = self.pathTemp
        self.pathTemp = self.pathOut + "temp.csv"
        # Psclim = self.pathScLim
        self.pathScLim = self.pathOut + ""

        return

    def todayInDate(self):
        return self.curDay, self.curMonth, self.curYear

    # Line 1373-1400
    def todayInJulian(self):
        day = int(time.strftime("%d"))
        month = int(time.strftime("%m"))
        year = int(time.strftime("%Y"))
        dayOfYear = int(time.strftime("%j"))
        todayInJul = self.dateToJulian(day, month, year, dayOfYear)
        return todayInJul

    # Should be working !!! NOT tested
    def julianToDate(self, julDate):
        yearAmounts = 0.25
        curYear = 0
        daysIn = 0
        # 7306 is the Julian Date for 2020
        for i in xrange(0, 7306):
            if i <= julDate and (i + int(yearAmounts) + 365) > julDate:
                daysIn = i + int(yearAmounts)
                break
            curYear += 1
            yearAmounts += 0.25
            # Move i to the next julian year
            i += (yearAmounts + 365)
        year = curYear + 2001
        result = datetime.datetime(year, 1, 1) + datetime.timedelta(daysIn)
        month = result.month
        day = result.day
        # print "----- Date is ", type(result) 
        # print "Year is ", year
        # Take out the year <-- Highest subdivision
        # Take out the month
        # Take out the day
        return day, month, year

    def dateToJulian(self, day, month, year, dayOfYear):
        milleniaRemainder = year % 1000
        newYear = (milleniaRemainder - 1) * 365
        # Account for leap years
        newYear += ((year - 2001) / 4)
        newYear += dayOfYear
        return newYear

    # Line 1625-1633
    ''' Lastrun (I think) defines if the program was a success or not. 
    Starts with 2 being written, then writes -1 after "Normal End"-ing'''

    def lastRun(self, curDrive, status):
        # Overwrite the previous "lastrun.txt" file
        lastRunF = open((curDrive + r"projects\argo\status\lastrun.txt"), 'w')
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
