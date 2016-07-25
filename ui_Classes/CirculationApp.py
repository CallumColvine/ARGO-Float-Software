''' TimeSeriesApp.py
Callum Colvine - Research Assistant
Callum.Colvine@dfo-mpo.gc.ca

Following Pep 8 formatting with the following exceptions:
- There is no spacing between a docstring and a function

'''

# UI Imports
from PySide import QtCore, QtGui
from PySide.QtGui import QWidget
# Plotting Imports
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
# Calculations imports
import numpy as np
from datetime import date, timedelta
import csv
# Util file Imports
from utilities.ARGO_Utilities import formatToDateTime, dateToJulian, \
    julianToDate, getProfile, despike, fillgaps, getSigmaT, getSvanom, \
    convertLatLonToNegative
# Formatting Imports
import re

from ui_Files.ui_circulationapp import Ui_CirculationApp



class CirculationApp(QWidget, Ui_CirculationApp):
    def __init__(self, parent):
        super(CirculationApp, self).__init__()
        self.setupUi(self)
        return

    def experimentSelected(self):
        self.initAllClassVariables()
        self.loadDefaultSettings()
        self.setupSignals()
        self.readEndOfFiles()
        return

    def initAllClassVariables(self):
        self.plotCentre = None
        self.rangeCentre = None
        self.firstLatitude = self.firstLatitudeBox.value()
        self.secondLatitude = self.secondLatitudeBox.value()
        self.firstLongitude = self.firstLongitudeBox.value()
        self.secondLongitude = self.secondLongitudeBox.value()
        self.sampleWindow = self.sampleWindowBox.value()
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

        self.Te = np.empty((500, 800))
        self.Sa = np.empty((500, 800))
        self.Lat = np.empty((800))
        self.Lon = np.empty((800))
        self.floats = []    # 800 long

        self.Dh = np.empty((800))
        self.Dhf = np.empty((800))
        self.Dhk = np.empty((800))
        self.Dhsave = np.empty((800))
        self.accept = [True] * 800    # 800 long

        self.S = ""         # 9000 length string
        self.T = ""         # 100 length string
        self.path0 = ""
        self.outPath = ""
        self.path = ""
        self.pathHg = ""
        self.gif = ""
        self.hlt = [[""]]
        
        self.dHPa = np.empty((100))
        self.vPa = np.empty((100))
        self.xKMa = np.empty((100))
        self.dList = np.empty((900))
        self.dListX = np.empty((900))
        self.dListY = np.empty((900))
        
        self.line1 = ""
        self.line2 = ""
        self.line3 = ""
        self.line4 = ""
        self.line5 = ""
        self.line6 = ""

        self.fil = ""
        self.filTe = ""
        self.filSa = ""
        self.filSt = ""

        self.sTe = ""
        self.sSa = ""
        self.sSt = ""
        self.sSp = ""
        self.sDh = ""
        self.sShgt = ""
        self.Pa = ""
        self.Qdt = ""

        self.listFil = ""
        self.param = ""

        self.P = np.empty((2000))
        self.T = np.empty((2000))
        self.S = np.empty((2000))

        self.Nx = 65
        self.Ny = 91

        self.Dx = 78.63
        self.Dy = 111.2 / 3.0

        # This seems dumb and non Pythonic. Redo
        self.E = np.empty((20, self.Ny, self.Nx))
        self.evals = np.empty((20))
        self.Dhfit = np.empty((self.Ny, self.Nx))
        # ,Deld(193),Xp(4),Yp(4),Maxlat(20),Maxval(20)
        # 150 INTEGER Depi(361,361)
        # 160 DIM Cov(20,20),Tr(20)
        # 170 MAT Cov=(0.)
        # 180 DIM Halat(133),Halon(193),Ha(193,133)
        # ...
        # ...
        # ...
        self.argoPath = "C:\Users\\ColvineC\\IOS_DFO\\ARGO-Float-Software\\"
        self.outPath = self.argoPath + "argo_out_TEST\\Circulation\\"
        self.drive = "P:\\"
        self.path0 = self.drive + "argo_mirror\\pacific_ocean\\"

        self.Scy=111.2 #! Scy = km/degree of latitude
        self.Scx=78.3
        self.Rho0=150 #! Decay scale of weighting function
        self.Pref=1000 #! Reference level for dynamic height calculations

        self.stepSize = self.stepSizeBox.value()
        self.nPress = 0
        self.updateNPress()
        return

    def readEndOfFiles(self):
        # Read eofs using non-HTBasic code
        self.E.fill(np.nan)
        modes = open((self.outPath + "hjf_modes.31"), 'r')
        n = 0
        k = 0
        modes.next()
        for line in modes:
            n += 1 
            if n == 4311:
                k += 1
                n = 0
                # This skips the intermediary lines between sectons of 4310
                continue
            # Default parameters for split is whitespace
            split = line.split()
            # print "it is ", split
            lon = float(split[1])
            lat = float(split[2])
            value = float(split[3])

            i = int(0.002 + 1.0 * (lon - 180.0))
            j = int(0.002 + 3.0 * (lat - 30.0))        
            # print "Settings to val", value
            self.E[k, j, i] = value
            # WHY IS THIS HERE?!?!
            # if (lat >= 51.7) and (lon <= 200):
            #     y = 51.7 + .033 * (lon - 180.0) * (lon - 180.0)
            #     if lat > y:
            #         if k == 0:
            #             print "K IS 0 and setting no nan"
            #             print "j problem is ", j, "i problem is ", i
            #             print "lat problem is ", lat, "lon problem is ", lon
            #         # print "lat is ", lat, " y is ", y
            #         # print "on iters k", k, "j", j, "i", i
            #         self.E[k, j, i] = np.nan
            #         # print ""
        return

    def setupSignals(self):
        self.setupInputParameterSignals()
        self.setupButtonPressSignals()
        return

    def setupButtonPressSignals(self):
        self.nextButton.clicked.connect(self.nextButtonClicked)
        return

    def setupInputParameterSignals(self):
        self.firstLatitudeBox.editingFinished.connect(
            self.firstLatitudeBoxEditingFinished)
        self.secondLatitudeBox.editingFinished.connect(
            self.secondLatitudeBoxEditingFinished)
        self.firstLongitudeBox.editingFinished.connect(
            self.firstLongitudeBoxEditingFinished)
        self.secondLongitudeBox.editingFinished.connect(
            self.secondLongitudeBoxEditingFinished)
        self.sampleWindowBox.editingFinished.connect(
            self.sampleWindowBoxEditingFinished)
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

    def sampleWindowBoxEditingFinished(self):
        self.sampleWindow = self.sampleWindowBox.value()
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

    def updateNPress(self):
        self.nPress = 1 + int(self.maxInterpDepth / self.stepSize)
        # if self.nPress > 300:
        #     self.nPress = 300
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

    def nextButtonClicked(self):
        self.saveChosenSettings()
        self.circulationStackedWidget.setCurrentWidget(self.calculatingPage)
        self.mainLoop()
        return

    def saveChosenSettings(self):
        # Saves the user picked settings to a *.cfg file
        return

    def loadDefaultSettings(self):
        # Read default/previous from a *.cfg file
        return

    def mainLoop(self):
        # self.numFloats = 0
        doSort = True
        julStart, julEnd = self.getJulianStartAndEnd()
        # Non-inclusive, should I do (julEnd + 1) ?
        # floats = []
        self.numFloats = 0
        for i in xrange(julStart, julEnd):
            yearMonthDayPath0 = self.checkFloatsFromIndex(i, self.path0)
        iFlt = 0
        if self.numFloats == 0:
            print "There are no floats within the provided range"
            return
        for flt in self.floats:
            numRecs = getProfile(flt, self.P, self.T, self.S)
            despike(self.T, numRecs)
            despike(self.S, numRecs)
            # fillgaps(self.T)
            # fillgaps(self.S)
            skip = self.checkSkip()
            if skip:
                continue
            self.dataForcing()
            self.storeData(numRecs, iFlt)
            iFlt += 1
        # Specific volume and Dynamic Heights
        # Replaced Deltap with self.stepSize
        self.specificVolAndDynH()
        self.meanAndStdDev()
        if doSort:
            self.removeDuplicates()
        dynHeightBar, dynHeightVar = self.recomputeMeanAndVariation()
        self.subtractMeanFromStored(dynHeightBar)
        self.bestFitMode(dynHeightVar)
        self.mapCirculations(dynHeightBar)
        divSl = self.findDivingStreamline()
        xPLL, xPLH, yPLL, yPLH = self.collectPlotData(divSl)
        self.contour(xPLL, xPLH, yPLL, yPLH)
        self.saveOutput()
        return

    def saveOutput(self):

        outFilePath = self.outPath + "Dh" + str(self.plotCentre) + "_0000.csv" 
        outFile = open(outFilePath, 'w')
        for I in xrange(0, self.Nx):
            Lonc = I - 180.
            for J in xrange(0, self.Ny):
                if not np.isnan(self.Dhfit[J, I]):
                    Latc = 30. + J / 3.
                    S = (str(Lonc) + "," + str(Latc) + "," + 
                         str(self.Dhfit[J, I]) + '\n')
                    outFile.write(S)
        return

    def collectPlotData(self, Divsl):
        Xp = np.empty((4))
        Yp = np.empty((4))
        XpListLow = []
        XpListHigh = []
        YpListLow = []
        YpListHigh = []

        for Idiv in xrange(0, 64):
            X0 = self.Dx * (Idiv - 1.)
            for Jdiv in xrange(0, 90):
                Y0 = self.Dy * (Jdiv - 1.)
                X1 = self.Dhfit[Jdiv, Idiv]
                X2 = self.Dhfit[Jdiv+1, Idiv]
                X3 = self.Dhfit[Jdiv+1, Idiv+1]
                X4 = self.Dhfit[Jdiv, Idiv+1]
                Xmax = max(X1,X2,X3,X4)
                Xmin = min(X1,X2,X3,X4)
                if(Xmax > 90. or Xmax < Divsl or Xmin > Divsl):
                    continue
                # Contour passes through the box
                Ip = 0

                Xmax = max(X1, X2)
                Xmin = min(X1, X2)
                if(Divsl >= Xmin and Divsl <= Xmax):
                    # print "Adding 1 to Xp and Yp"
                    Xp[Ip] = X0
                    Yp[Ip] = Y0 + self.Dy * (Divsl - X1) / (X2 - X1)
                    Ip = Ip + 1
                Xmax = max(X2,X3)
                Xmin = min(X2,X3)
                if(Divsl>=Xmin and Divsl<=Xmax):
                    # print "Adding 2 to Xp and Yp"
                    Xp[Ip] = X0 + self.Dx * (Divsl - X2) / (X3 - X2)
                    Yp[Ip] = Y0 + self.Dy
                    Ip=Ip+1
                Xmax = max(X3,X4)
                Xmin = min(X3,X4)
                if(Divsl >= Xmin and Divsl <= Xmax):
                    # print "Adding 3 to Xp and Yp"
                    Xp[Ip] = X0 + self.Dx
                    Yp[Ip] = Y0 + self.Dy * (Divsl - X4) / (X3 - X4)
                    Ip = Ip + 1
                Xmax = max(X1,X4)
                Xmin = min(X1,X4)
                if(Divsl >= Xmin and Divsl <= Xmax):
                    # print "Adding 4 to Xp and Yp"
                    Xp[Ip] = X0 + self.Dx * (Divsl - X1) / (X4 - X1)
                    Yp[Ip] = Y0
                    Ip = Ip + 1
                XpListLow.append(Xp[0])
                XpListHigh.append(Xp[1])
                YpListLow.append(Yp[0])                
                YpListHigh.append(Yp[1])
                # MOVE Xp(1),Yp(1)
                # DRAW Xp(2),Yp(2)
        # print "XpList is ", XpList
        # print "YpList is ", YpList
        return XpListLow, XpListHigh, YpListLow, YpListHigh 

    def findDivingStreamline(self):
        # Now search for dividing streamline
        Divsl = 0.
        Knt = 0
        # ToDo: Was 58, changed to 59
        for Latd in xrange(40, 59):
            J = int(1.01 + 3. * (Latd - 30.)) - 1
            for Im in xrange(0, 60):
                I = 64 - Im
                if np.isnan(self.Dhfit[J, I]): 
                    continue
                Dh1 = self.Dhfit[J, I]
                break
            # ! Have Dh1
            Divsl = Divsl + Dh1
            Knt = Knt + 1
        Divsl = Divsl / Knt
        return Divsl

    def prepareCallContour(self):
        # ! CLEAR SCREEN
        contourInterval = 0.05
        # CALL Contour(Dhfit(*),Dx,Dy,Nx,Ny,Conin)
        self.contour()
        # Xl = Dx * (self.Nx - 1.0)
        # Yl = Dy * (self.Ny - 1.0)
        return

    # self.Dhfit, Dx, Dy, self.Nx, self.Ny
    def contour(self, xPLL, xPLH, yPLL, yPLH):
        # setup polyconic basemap
        # by specifying lat/lon corners and central point.
        # area_thresh=1000 means don't plot coastline features less
        # than 1000 km^2 in area.
        latLeft, lonLeft = convertLatLonToNegative(self.firstLatitude, 
                                                   self.firstLongitude)
        latRight, lonRight = convertLatLonToNegative(self.secondLatitude,
                                                     self.secondLongitude)
        # print "lonLeft", lonLeft, "latLeft", latLeft, "latRight", latRight, "lonRight", lonRight
        m = Basemap(llcrnrlon=lonLeft,llcrnrlat=latLeft,urcrnrlon=lonRight, \
                    urcrnrlat=latRight,\
                    resolution='l',area_thresh=1000.,projection='poly',\
                    lat_0=50,lon_0=-140)
        # m = Basemap(llcrnrlon=-160,llcrnrlat=43,urcrnrlon=-100,urcrnrlat=57,\
        #             resolution='l',area_thresh=1000.,projection='poly',\
        #             lat_0=50,lon_0=-140)
        m.drawcoastlines()
        m.fillcontinents(color='darksage',lake_color='royalblue')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))
        m.drawmapboundary(fill_color='royalblue')
        plt.title("Circulation Data on the Southern Alaskan Coast")
        # for x in xrange(0, len(xPLL) - 1):
        #     plt.plot([xPLL[x], xPLH[x]], [yPLL[x], yPLH[x]])
         
        plt.show()
        return

    def mapCirculations(self, dynHeightBar):
        # Now add mean and modes
        self.Dhfit.fill(dynHeightBar)
        for i in xrange(0, self.Nx):
            for j in xrange(0, self.Ny):
                if not np.isnan(self.E[0, j, i]):
                    for k in xrange(0, self.totalModes): 
                        self.Dhfit[j, i] += (self.evals[k] * self.E[k, j, i])
                        # if i == 0 and j < 20:
                        #     # print "self.evals[k] is ", self.evals[k], "self.E[k, j, i] is ", self.E[k, j, i]
                        #     # print "dynHeightBar is ", dynHeightBar
                        #     pass
                else:
                    self.Dhfit[j, i] = np.nan
                    # print "Dhfit nan set on j", j, " i ", i 
        return

    def bestFitMode(self, Dhvar):
        Frvarlast = 1.
        for Mq in xrange(0, self.totalModes):
            if Mq > 3.5:
                M = Mq
            if Mq < 3.5:
                M = 4 - Mq
            M = Mq
            # ! Find best fit between Mode M and the data
            Nent = 0
            for Iflt in xrange(0, self.numFloats):
                if self.accept[Iflt]:
                    Ef = self.interpolate(self.E, M, self.Nx, self.Ny, 
                                          self.Lat[Iflt], self.Lon[Iflt])
                    # print "Ef is ", Ef
                    self.Dhf[Iflt] = Ef
                    Nent = Nent + 1
            Q1 = 0.
            Q2 = 0.
            Ntemp = 0
            for Iflt in xrange(0, self.numFloats):
                if not np.isnan(self.Dhf[Iflt]) and self.accept[Iflt]:
                    Q1 += self.Dh[Iflt] * self.Dhf[Iflt]
                    # print "self.Dh[Iflt]", self.Dh[Iflt], "self.Dhf[Iflt]", self.Dhf[Iflt]
                    Q2 += self.Dhf[Iflt] * self.Dhf[Iflt]
                    Ntemp = Ntemp + 1
            print "Q1 is ", Q1, " Q2 is ", Q2
            self.evals[M] = Q1 / Q2
            Rvar = 0.
            for Iflt in xrange(0, self.numFloats): 
                if not np.isnan(self.Dhf[Iflt]) and self.accept[Iflt]:
                    self.Dh[Iflt] = (self.Dh[Iflt] - 
                                     self.evals[M] * self.Dhf[Iflt])
                    Rvar = Rvar + self.Dh[Iflt] * self.Dh[Iflt]
            Rvar = Rvar / Ntemp
            Frvar = Rvar / Dhvar
            # IF Verb_opt>.5 THEN PRINT "After fitting mode #";M;" residual variance = ";Frvar
            Contvar = Frvarlast - Frvar
            Frvarlast = Frvar
            # OUTPUT @Pvar;VAL$(M)&","&VAL$(Frvar)&","&VAL$(Contvar)
        # ASSIGN @Pvar TO *
        return

    def interpolate(self, E, M, Nx, Ny, Lati, Loni):
        Eval = np.nan
        Y = 1. + 3. * (Lati - 30.)
        X = 1. + 1. * (Loni - 180.)
        Ix = int(X + .0001)
        Iy = int(Y + .0001)
        Rx = X - Ix
        Ry = Y - Iy
        
         # Outside bounds, so ignore
        if Ix < 1 or Ix >= Nx:
            print "returning early X"
            return Eval
        if Iy < 1 or Iy >= Ny:
            print "returning early Y"
            return Eval
         # Integer lat/lon, so easy
        if np.abs(Rx) < .01 and np.abs(Ry) < .01:
            Eval = self.E[M, Iy, Ix]
            print "Returning 3rd if Eval", Eval
            return Eval       
        I1 = Ix - 1
        I2 = Ix + 2
        J1 = Iy - 1
        J2 = Iy + 2
        
        sumWeight = 0.
        sumWeightE = 0.
        A2 = .25
        for I in xrange(I1, I2):
            if I < 1 or I > Nx:
                continue
            for J in xrange(J1, J2):
                if J < 1 or J > Ny or np.isnan(self.E[M, J, I]):
                    continue
                Rho2 = (I - X) * (I - X) + (J - Y) * (J - Y)
                Wgt = np.square(-Rho2 / A2)
                sumWeight += Wgt
                sumWeightE += self.E[M, J, I] * Wgt
                # print "self.E[M, J, I]", self.E[M, J, I], "Wgt", Wgt
        if sumWeight > .2:
            Eval = sumWeightE / sumWeight
        #     print "Returning 4th Eval", Eval
        # print "Returning end Eval"
        return Eval

    def subtractMeanFromStored(self, dynHeightBar):
        # MAT Dhsave=Dh
        for Iflt in xrange(0, self.numFloats):
            if self.accept[Iflt]:
                self.Dh[Iflt] = self.Dh[Iflt] - dynHeightBar
        # MAT Dhk=Dh
        return

    def recomputeMeanAndVariation(self):
        Dhbar = 0.
        Dhvar = 0.
        Mfloats = 0
        for iFlt in xrange(0,  self.numFloats):
            print "numFloats looped"
            if self.accept[iFlt]:
                print "Accept the iFlt one"
                Dhbar = Dhbar + self.Dh[iFlt]
                Dhvar = Dhvar + self.Dh[iFlt] * self.Dh[iFlt]
                Mfloats = Mfloats + 1
        Dhbar = Dhbar / Mfloats
        Dhvar = Dhvar / Mfloats
        Dhvar = Dhvar - Dhbar * Dhbar
        # PRINT "Mean Dh = ";Dhbar;" Variance in Dh = ";Dhvar;Mfloats
        return Dhbar, Dhvar

    def removeDuplicates(self):
        # ! Now have a list of good float observations, bad ones have Accpt(i)=-1
        print "Checking for duplicates"
        # ! Sort entries so that there is only one entry per float
        # ! First see if there are any duplicates
        for iFlt in xrange(0, self.numFloats-1):
            if not self.accept[iFlt]: 
                continue
            flt = self.floats[iFlt]
            splitName = re.split('[\W_]+', flt)
            # Gets 2nd last element from split string. Last element is "IOS"
            endOfName = splitName[-2]
            # print "Split name is ", splitName
            for jFlt in xrange(iFlt + 1, self.numFloats):
                if(not self.accept[jFlt]): 
                    continue
                # duplicate found
                if(self.floats[jFlt] == endOfName):
                    self.foundDuplicate(endOfName, iFlt)
        return

    def foundDuplicate(self, endOfName, dFlt):
        # ! at least one duplicate is present for float dFlt
        Dhbar = 0.
        Latav = 0.
        Lonav = 0.
        Nav = 0
        for jFlt in xrange(dFlt, self.numFloats):
            if (self.floats[jFlt] == endOfName): 
                continue
            if (self.Dh[jFlt] < 0 or self.Dh[jFlt]>5):
                self.accept[jFlt] = False
            if (not self.accept[jFlt]): 
                continue
            Dhbar = Dhbar + self.Dh[jFlt]
            Latav = Latav + self.Lat[jFlt]
            Lonav = Lonav + self.Lon[jFlt]
            Nav = Nav + 1
        if Nav > .5:
            self.Dh[dFlt] = Dhbar / Nav
            self.Lat[dFlt] = Latav / Nav
            self.Lon[dFlt] = Lonav / Nav
        # if Verb_opt>.5 THEN PRINT "Eliminating duplicate entries of "&Q$;Nav;" of them."
        for jFlt in xrange(dFlt+1, self.numFloats):
            if (self.floats[jFlt] == endOfName): 
                self.accept[jFlt] = False
        return

    def meanAndStdDev(self):
        # ! Now compute mean and stnd dev and remove mean from Dh(*)
        print "Computing mean and standard deviation"
        for i in xrange(0, 10):
            Dhbar = 0.
            Dhvar = 0.
            Mfloats = 0
            Nrej = 0
            for iFlt in xrange(0, self.numFloats):
                if self.accept[iFlt]:
                    Dhbar = Dhbar + self.Dh[iFlt]
                    Dhvar = Dhvar + self.Dh[iFlt] * self.Dh[iFlt]
                    Mfloats = Mfloats + 1
            Dhbar = Dhbar / Mfloats
            Dhvar = Dhvar / Mfloats
            Dhvar = Dhvar - Dhbar * Dhbar

            Stdev = np.sqrt(Dhvar)
            Devmax = 0
            for iFlt in xrange(1, self.numFloats):
                if self.accept[iFlt]:
                    Devn = np.fabs((self.Dh[iFlt] - Dhbar) / Stdev)
                    if (Devn > Devmax):
                        Devmax = Devn
                    if (Devn > 3):
                        self.accept[iFlt] = False
                        # IF Verb_opt > .5 THEN PRINT "Rejecting float # ";Iflt
                        Nrej = Nrej + 1

            # IF Verb_opt>.5 THEN PRINT "Maximum deviation = ";Devmax;" standard deviations."
            if(Nrej < .5): 
                break
        return

    def specificVolAndDynH(self):
        Ipref = 1 + self.Pref / self.stepSize
        for iFlt in xrange(0, self.numFloats):
            for iPress in xrange(0, 500):
                pressCount = self.stepSize * (iPress - 1.)
                qTe = self.Te[iFlt, iPress]
                qSa = self.Sa[iFlt, iPress]
                if (qSa > 900):
                    print "Skipping from qSa"
                    continue
                if (qTe > 900):
                    print "Skipping from qTe"
                    continue
                qSt = getSigmaT(qSa, qTe)
                sigma, Svan = getSvanom(qSa, qTe, 0)
                # print "Svan is ", Svan
                self.Sa[iFlt,iPress] = Svan
                
            # ! Compute dynamic height at surface P = Pmap relative to 1000
            Dhi = 0
            Q = 0.5 * 1.0E-5
            Ipr0 = int(0.01 + self.dynHeightsAtP / self.stepSize)
            Ipr1 = Ipr0 + 1

            for iPress in xrange(Ipr1, Ipref):
                Dhi += (Q * (self.Sa[iFlt, iPress] + 
                            self.Sa[iFlt, iPress-1]) * 
                        self.stepSize)
            self.Dh[iFlt] = Dhi
        return

    def storeData(self, numRecs, numFloat):
        for iPr in xrange(0, self.nPress):
            pressCount = self.stepSize * iPr
            if pressCount < self.P[0]: 
                continue
            if pressCount > self.P[numRecs]:
                return
            for i in xrange(1, numRecs):
                if pressCount >= self.P[i-1]:
                    rho = 0
                    if(np.abs(self.P[i] - self.P[i-1]) > 1.0E-2): 
                        rho = ((pressCount - self.P[i-1]) / 
                               (self.P[i] - self.P[i-1]))

                    self.Te[numFloat, iPr] = 999.9
                    if(np.abs(self.T[i-1]) < 50. and np.abs(self.T[i]) < 50.):
                        self.Te[numFloat, iPr] = \
                            self.T[i-1] + rho * (self.T[i] - self.T[i-1])
                    
                    self.Sa[numFloat, iPr] = 999.9
                    if (np.abs(self.S[i-1]) < 50. and np.abs(self.S[i]) < 50.): 
                        self.Sa[numFloat, iPr] = \
                            self.S[i-1] + rho * (self.S[i] - self.S[i-1])
        #     5320 ! Finished interpolation to standard pressures for Float Iflt
        # 5340 ! Finished interpolation to standard pressures for all floats
        return

    def dataForcing(self):
        # 4950 ! If the sample is "near-surface" then force it to be surface
        # 4960 IF P(1)<20. THEN P(1)=0.
        # 4970 !
        # 4980 ! If the sample is almost deep enough, force it to be deep enough
        # 4990 IF P(Nrecs)<Pmax THEN
        # 5000 Dtdp=(T(Nrecs)-T(Nrecs-1))/(P(Nrecs)-P(Nrecs-1))
        # 5010 Dsdp=(S(Nrecs)-S(Nrecs-1))/(P(Nrecs)-P(Nrecs-1))
        # 5020 Text=T(Nrecs)+Dtdp*(Pmax-P(Nrecs))
        # 5030 Sext=S(Nrecs)+Dsdp*(Pmax-P(Nrecs))
        # 5040 S(Nrecs)=Sext
        # 5050 T(Nrecs)=Text
        # 5060 P(Nrecs)=Pmax
        return

    def checkSkip(self):
        # 4830 IF Rflag<-.1 THEN GOTO 5330
        # 4840 Ifirst=0
        # 4850 Ifirst=Ifirst+1
        # 4860 IF Ifirst>20 THEN GOTO 5330
        # 4870 IF P(Ifirst)<.01 THEN GOTO 4850
        # 4880 IF P(Ifirst)>9000. THEN GOTO 4850
        # 4890 IF P(Ifirst)>18. THEN GOTO 5330
        # 4910 IF Rflag<0 THEN GOTO 5330
        return False

    def getJulianStartAndEnd(self):
        date = self.centreOfPlotDateEdit.date().getDate()
        dateTimeObj = formatToDateTime(date[0], date[1], date[2])
        dayOfYear = dateTimeObj.timetuple().tm_yday 
        julDate = dateToJulian(date[2], date[1], date[0], dayOfYear)
        self.plotCentre = julDate
        julStart = julDate - self.sampleWindow
        julEnd = julDate + self.sampleWindow
        print "julStart: ", julStart, " julEnd: ", julEnd
        return julStart, julEnd

    ''' Not in a utils file due to repeated use of class variables. ''' 
    def checkFloatsFromIndex(self, cycleJulDate, path0):
        day, month, year = julianToDate(cycleJulDate)
        yearMonthDayPath0 = (path0 + 
                             str(year) + '\\' + 
                             str(month).zfill(2) + '\\' + 
                             str(day).zfill(2) + '\\')
        inFileName = (yearMonthDayPath0 + 
                      str(year) + 
                      str(month).zfill(2) + 
                      str(day).zfill(2) + 
                      '_index.csv')
        # print "inFileName is ", inFileName
        try:
            with open((inFileName),
                      'rb') as indexCSV:
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
                        self.floats.append(yearMonthDayPath0 + row[0])
                        self.passedFloat(floatNum, lat, lon)
                        self.numFloats += 1
        except (OSError, IOError), e:
            print "File was not found: ", inFileName
        return yearMonthDayPath0

    def passedFloat(self, floatNum, lat, lon):
        self.Lat[self.numFloats] = lat
        self.Lon[self.numFloats] = lon
        # dX = self.sCX0 * (float(lon) - self.longitudeDesired)
        # dY = self.sCY * (float(lat) - self.latitudeDesired)
        # rho = np.sqrt(dX * dX + dY * dY)
        # if self.closestDist > rho:
        #     self.closestDist = rho
        #     self.closestFloatNum = floatNum
        # ToDo: A few lines left out here, looks completely unnecessary
        return
