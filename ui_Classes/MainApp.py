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
from datetime import datetime
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
        self.floatRadius = 0
        self.dayRange = 0
        self.dayStepSize = 0
        self.sampleWindow = 0
        self.firstLatitude = self.firstLatitudeBox.value()
        self.secondLatitude = self.secondLatitudeBox.value()
        self.firstLongitude = self.firstLongitudeBox.value()
        self.secondLongitude = self.secondLongitudeBox.value()
        self.pressureCutOff = self.pressureCutOffBox.value()
        self.maxInterpDepth = self.maxInterpDepthBox.value()
        self.stepSize = self.stepSizeBox.value()
        self.nPress = 1 + int(0.001 + self.maxInterpDepth / self.stepSize)
        self.temp = False
        self.salinity = False
        self.sigmaT = False
        self.spiciness = False
        self.dynamicHeight = False
        self.latitudeDesired = 0
        self.longitudeDesired = 0
        self.append = True
        self.verbose = False
        self.outPath = None
        self.drive = ''
        self.path0 = ''
        self.yearMonthDayPath0 = ''
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
        return

    # Start of parameter boxes
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
        R4 =4.8314E-4
        Dr350 = 28.106331
        Sr = np.sqrt(S)
        R1 = ((((6.536332E-9 * T - 1.120083E-6) * T + 1.001685E-4) * T - 9.09529E-3) * T + 6.793952E-2) * T - 28.263737
        R2 = (((5.3875E-9 * T - 8.2467E-7) * T + 7.6438E-5) * T - 4.0899E-3) * T + 8.24493E-1
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
        self.extractDataFromFloatFile(path, order)          
        # He opens this path?
        numChannels = 0
        numRecords = 0
        fileTempQc = 0
        filePSalQc = 0
        # Line1 -> 16210 close file, stop logging, lastRun
        if not self.checkIfReturn(numRecords, fileTempQc, filePSalQc):
            print "in getProfile, returning prematurely due to bad values"
            return
        self.channelOrder()
        return

    def extractDataFromFloatFile(self, path, order):
        print "Opening path ", path
        floatFile = open(path, 'r')
        passedEndOfHeader = False
        pFound = False
        tFound = False
        sFound = False
        i = 0
        for line in floatFile:
            if not passedEndOfHeader:
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
        return 

    def appendInfo(self, line, order, i):
        split = line.split()
        print "Split is ", split
        offset = 5
        position = 0
        for data in order:
            if data == 'P':
                if abs(float(split[2 + (position * 5)])) < 9000:
                    self.P[i] = float(split[2 + (position * 5)])
                else:
                    self.P[i] = float(split[0 + (position * 5)])
            elif data == 'T':
                if abs(float(split[2 + (position * 5)])) < 9000:
                    self.T[i] = float(split[2 + (position * 5)])
                else:
                    self.T[i] = float(split[0 + (position * 5)])
            elif data == 'S':
                if abs(float(split[2 + (position * 5)])) < 9000:                
                    self.S[i] = float(split[2 + (position * 5)])
                else:
                    self.S[i] = float(split[0 + (position * 5)])
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

    def commenceInterpolation(self,):
        # print some stuff for user
        sTav = 0
        floats = self.checkFloatsFromIndex(self.yearMonthDayPath0)
        rFlag = -1
        for singleFloat in floats:
            self.getProfile(singleFloat, self.P, self.T, self.S, rFlag)
            self.sanityCheck()
            sigT = self.getSigmaT(self.S, self.T)
            self.randomCalcs(sigT)
        return

    def randomCalcs(self, sigT):
        # Copied, not sure why Howard does this
        sigT = .001 * int(.5 + 1000. * sigT)
        print "sigT is ", sigT
        if sigT > 5 and sigT < 30:
            sTav += sigT
            sTavSq += np.power(sigT, 2)
            sTKnt += 1
        print "P[1] is ", p[1]
        if P[1] < 20:
            p[1] = 0
        # Re-intitializes Te[x][y] and Sa to 999.9
        print "nPress is ", self.nPress
        for iterPress in xrange(1, self.nPress):
            pc = self.stepSize * (iterPress - 1) # What does this line do
            # 8330 REM IF Pc>Wdep THEN GOTO 5000
            Wgtsums=0.
            Wgtsumt=0.
            Swgtsum=0.
            Twgtsum=0.
            for iterFloat in range(1, numFloats):
                # IF Accpt(Iflt)<0. THEN GOTO 8610
                # Latav = latitude average?
                Latav = .5 * (Lat(iterFloat) + Latc)
                Scx=Scy*COS(Latav)
                Dx=Scx*(Lonc-Lon(iterFloat))
                Dy=Scy*(Latc-Lat(iterFloat))
                Rho=SQR(Dx*Dx+Dy*Dy)
                Z=Rho/Rho0
                print "Z is ", Z
                # IF Z>5. THEN GOTO 8610
                if Z > 5:
                    break
                Wgt=np.exp(-Z*Z)
                # !
                if Te(iterFloat,iterPress)>-1.5 and Te(iterFloat,iterPress)<40.:
                    Wgtsumt=Wgtsumt+Wgt
                    Twgtsum=Twgtsum+Wgt*Te(iterFloat,iterPress)

                if Sa(iterFloat,iterPress)>20. and Sa(iterFloat,iterPress)<40.:
                    Wgtsums=Wgtsums+Wgt
                    Swgtsum=Swgtsum+Wgt*Sa(iterFloat,iterPress)

            if Wgtsumt>Wgtsumt_max: 
                Wgtsumt_max=Wgtsumt
            Qte=Twgtsum/Wgtsumt
            Qsa=Swgtsum/Wgtsums
            
            Qst = self.getSigmaT(Qsa,Qte)
            self.Spiciness(Qsp,Qte,Qsa)
            CALL Svanom(Qsa,Qte,0.,Sigma,Svan)
            
            P(Ipr)=Deltap*(Ipr-1.)
            T(Ipr)=Qte
            S(Ipr)=Qsa
            St(Ipr)=Qst
            Sp(Ipr)=Qsp
            Kdel=Kdel+1
            Pdel(Kdel)=Pc
            Sva(Kdel)=Svan
            
        return

    def spiciness(self, spice, temp, salt):
        # Hardcoded by Howard
        B = np.array([0.0,.77442,-.00585,.000984,-.000206],
                [.051665,.002034,-.0002745,-.0000085,.0000136],
                [-6.64783E-3,-2.4681E-4,-1.428E-5,3.337E-5,7.894E-6],
                [-5.4023E-5,7.326E-6,7.0036E-6,-3.0412E-6,-1.0853E-6],
                [3.949E-7,-3.029E-8,-3.8209E-7,1.0012E-7,4.7133E-8],
                [-6.36E-10,-1.309E-9,6.048E-9,-1.1409E-9,-6.676E-10])
        Spice=0.
        Sp=Salt-35.
        Theta=Temp
        # Reversed I and J here, are arrays backwards in HP Basic?
        for I=1 in range(0, 5)
            Ii=I-1
            for J=1 in range(0, 4)
                Jj=J-1
                Spice = Spice + B(I, J)*(Theta^Ii)*(Sp^Jj)
                NEXT J
        NEXT I
        
        SUBEND
        return

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
                print "Interp depth ", float(row[3])
                # ToDo: Does not account for users with lat2 < lat1
                if (    float(row[1]) < self.firstLatitude 
                    and float(row[1]) > self.secondLatitude 
                    and float(row[2]) < self.firstLongitude
                    and float(row[2]) > self.secondLongitude
                    and float(row[3]) < self.pressureCutOff):
                        self.numFloats += 1
                        buoysToUse.append(row[0])
                        self.Lat[numFloats] = float(row[1])
                        self.Lon[numFloats] = float(row[2])
                    # Save it
                # print row[0]
        print "Bouys to use is ", buoysToUse
        return buoysToUse

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

        Pdel = np.empty((1000))
        Sva = np.empty((1000))
        Dh = np.empty((1000))
        # 40 DIM S$[100],T$[100],U$[500],Path0$[80],self.outPath$[80],Path$[80],R$[600],Closest_fltnm$[100]
        # Ignore all these declarations if they're strings, since I can do them at use-time
        Sa_av = np.empty(12)
        Sa_knt = np.empty(12)
        # 70 DIM Fil$[100],Filte$[100],Filsa$[100],Filst$[100],Out$[200],App$[200]
        # 80 DIM Ste$[500],Ssa$[500],Sst$[500],Ssp$[500],Sdh$[500],Shgt$[500],Pa$[500],Sclim$[500]
        # 90 DIM Location$[100],Param$[50],Lstmsg$[50]
        Sigref_st = np.empty((71))
        Sigref_te = np.empty((71))
        Sigref_sa = np.empty((71))
        # 110! Dimensions of data to be read from files
        self.P = np.empty((2000))
        self.T = np.empty((2000))
        self.S = np.empty((2000))
        St = np.empty((2000))
        Sp = np.empty((2000))
        # 140 ! ALPHA PEN 1
        # 150 Tabrow=11
        Scy = 111.2     # Scy = km/degree of latitude
        Rho0 = 300      # Decay scale of weighting function
        Pref = 1000     # Reference level for dynamic height calculations
        Closest = 9999.9
        # 240 ! Initialize data matrices with nul values 999.9
        # Should I do the NumPy NULL or just stick with 999.9?
        self.Te[:] = 999.9
        self.Sa[:] = 999.9    
        Lat[:] = 999.9
        Lon[:] = 999.9
        Accpt[:] = (-1)
        Sa_av[:] = 0
        Sa_knt[:] = 0
        # Dout backslash because in Python you need to escape special chearacters
        self.drive = "E:\\"
        # 350 CALL Lastrun(self.drive$,2)
        self.lastRun(self.drive, 2)
        self.path0 = self.drive + "argo_mirror\\pacific_ocean\\"
        self.outPath = self.drive + "argo_out\\"
        Sigpath = self.drive + "projects\\Sigma_Climate\\"
        # Made this file empty here when he writes 3x71 empties
        # if problems: change this maybe
        csvF = open((Sigpath + "Mp26_i.csv"), 'w')
        csvF.close()
        return

    ''' TODO 
    Pulls settings from an already existing file. INCOMPLETE'''
    def setSettings():
        # 470 ! Read previous settings
        # 480 ASSIGN @Pin TO Drive$&"argo_programs\TSlast.txt";FORMAT ON
        # ToDo: Change to 'a'
        settingsF = open((self.drive + "argo_programs\\TSlast.txt"), 'r+')
        for x in xrange(1,20):
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
        # self.dayRangeBox.setEnabled(isEnabled)
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
        if (not isEnabled):
            # self.loadDefaults()
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

        tempPath = self.outPath + "TS_Te.csv"
        salinityPath = self.outPath + "TS_Sa.csv"
        sigmaTPath = self.outPath + "TS_St.csv"
        spicinessPath = self.outPath + "TS_Sp.csv"
        dynamicHeightPath = self.outPath + "TS_Dh.csv"
        self.outputFilesLabel.setText("Output Files:\n" 
                + tempPath + '\n' 
                + salinityPath + '\n' 
                + spicinessPath + '\n' 
                + dynamicHeightPath)
        writeType = 'a'
        if not self.append:
            writeType = 'w'
        tempCSV = open((tempPath), writeType)
        salinityCSV = open((salinityPath), writeType)
        sigmaTCSV = open((sigmaTPath), writeType)
        spicinessCSV = open((spicinessPath), writeType)
        dynamicHeightCSV = open((dynamicHeightPath), writeType)
        # What the heck is hgt? Mercury...t?
        # hgtCSV = open((self.outPath + "TS_Shgt.csv"), writeType)
        # File on sigma-levels
        climCSV = open((self.outPath + "Sigref.csv"), writeType)    
        if (self.append):
            self.prepareFilesAndAppend()
        else:
            self.prepareFilesAndClear()
            # Maybe not necessary? 
        # Display "Starting interpolation value 1 to value 2" 
        tempCSV.close()
        salinityCSV.close()
        sigmaTCSV.close()
        spicinessCSV.close()
        dynamicHeightCSV.close()
        # hgtCSV.close()
        climCSV.close()
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
        yearAmounts = 0
        curYear = 0
        # 7306 is the Julian Date for 2020
        for i in xrange(0, 7306):
            if i <= julDate and (i + int(yearAmounts) + 365) > julDate:
                break
            curYear += 1
            yearAmounts += 0.25
            # Move i to the next julian year
            i += (yearAmounts + 365)
        year = curYear + 2001
        print "Year is ", year
        # Take out the year <-- Highest subdivision
        # Take out the month
        # Take out the day
        return

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

