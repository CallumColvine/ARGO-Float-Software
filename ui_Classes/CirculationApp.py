''' 
CirculationApp.py
Callum Colvine - Research Assistant
Callum.Colvine@dfo-mpo.gc.ca (until August 25th)
CallumColvine@gmail.com

CirculationApp builds on the functionality of Howard Freeland's CirculationApp 
which was written in HT Basic.

CirculationApp reads ARGO data from *_index.csv and *.IOS files to calculate the
dynamic height between longitudes: 180-224 and latitudes: 30-60. CirculationApp
then outputs this data into Dh****_0000.cav files and forms a contour plot with
the used float locations on it. 

Potential upgrades to CirculationApp:
- dhpa, xKMa investigate
- dList, x, y are part of plotting. Investigate
- O(n^3) loop in program, would be nice to rework
- Change mainLoop for loop calls to be inside functions (slightly faster)


Following Pep 8 formatting with the following exceptions:
- There is no spacing between a docstring and a function
- Using x[i][j] notation instead of x[i, j] since that is 
the format I was taught when initially learning in Java. 
This might (hopefully) make the code more readable 

'''

# UI Imports
from PySide import QtCore, QtGui
from PySide.QtGui import QWidget
from ui_Files.ui_circulationapp import Ui_CirculationApp
# Plotting Imports
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
# Calculations imports
from math import ceil
import numpy as np
from datetime import date, timedelta
import csv
# Util file Imports
from utilities.ARGO_Utilities import formatToDateTime, dateToJulian, \
    julianToDate, getProfile, despike, fillgaps, getSigmaT, getSvanom, \
    convertLatLonToNegative
# Formatting Imports
import re


class CirculationApp(QWidget, Ui_CirculationApp):
    def __init__(self, parent):
        super(CirculationApp, self).__init__()
        self.CLASS_INITTED = False
        self.signalsSetup = False
        self.setupUi(self)
        return

    ''' This is called when the user clicks the Circulation App button on the 
    homepage. This saves memory by not initializing lists unless the user 
    chooses this experiment. '''
    def experimentSelected(self):
        self.initAllClassVariables()
        self.loadDefaultSettings()
        if not self.signalsSetup:
            self.setupSignals()
        self.signalsSetup = True
        self.CLASS_INITTED = True
        print "Class init true"
        return

    ''' This method holds references to all class variables. For anyone unclear,
    Python allows any method inside the class to access these variables. They 
    work similarly to "Global" variables, but in the class scope.'''
    def initAllClassVariables(self):
        print "Initializing class variables"
        self.verbose = True
        # The first/second lat/lon are fixed due to modes file used
        self.firstLatitude = 30
        self.secondLatitude = 60
        self.firstLongitude = 180
        self.secondLongitude = 244
        self.sampleWindow = self.sampleWindowBox.value()
        self.pressureCutOff = self.pressureCutOffBox.value()
        self.maxInterpDepth = self.maxInterpDepthBox.value()
        self.pressureStepSize = self.stepSizeBox.value()
        self.dynHeightsAtP = self.dynHeightAtPBox.value()
        self.relativeToPref = self.relativeToPrefBox.value()
        self.totalModes = self.totalModesBox.value()

        self.Te = np.empty((800, 800))
        self.Sa = np.empty((800, 800))
        self.Lat = np.zeros((800))
        self.Lon = np.zeros((800))

        self.Dh = np.empty((800))
        self.Dhf = np.empty((800))
        self.Dhk = np.empty((800))
        self.Dhsave = np.empty((800))
        self.accept = [True] * 800    # 800 long

        self.S = ""         # 9000 length string
        self.T = ""         # 100 length string
        self.path0 = ""
        self.outPath = ""
        
        # self.dHPa = np.empty((100))
        # self.xKMa = np.empty((100))
        # self.dList = np.empty((900))
        # self.dListX = np.empty((900))
        # self.dListY = np.empty((900))
        
        # self.line1 = ""
        # self.line2 = ""
        # self.line3 = ""
        # self.line4 = ""
        # self.line5 = ""
        # self.line6 = ""

        self.P = np.empty((2000))
        self.T = np.empty((2000))
        self.S = np.empty((2000))

        self.nX = 65
        self.nY = 91

        self.dX = 78.63
        self.dY = 111.2 / 3.0

        # This seems non Pythonic. ToDo
        self.E = np.empty((20, self.nY, self.nX))
        self.evals = np.empty((20))
        self.dHFit = np.empty((self.nY, self.nX))

        # self.argoPath = "C:\Users\ColvineC\IOS_DFO\ARGO-Float-Software\\"
        # self.outPath = self.argoPath + "argo_out_TEST\\Circulation\\"
        # self.drive = "P:\\"
        # self.path0 = self.drive + "argo_mirror\pacific_ocean\\"

        print "Init my path in Circ"
        self.argoPath = ''
        self.path0 = ''
        self.outPath = ''

        # Used for plotting
        self.scy = 111.2        #! Scy = km/degree of latitude
        self.scx = 78.3
        self.stepSize = self.stepSizeBox.value()
        self.nPress = 0
        self.updateNPress()
        self.circulationStackedWidget.setCurrentWidget(self.settingsPage)
        return

    ''' Called initially, this function reads the necessary values from Howard's 
    modes file.'''
    def readEVals(self):
        self.E.fill(np.nan)
        modes = open((self.outPath + "hjf_modes.31"), 'r')
        n = 0
        k = 0
        modes.next()
        iMax = 0
        for line in modes:
            n += 1 
            if n == 4311:
                k += 1
                n = 0
                # This skips the intermediary lines between sectons of 4310
                continue
            # Default parameters for split is whitespace
            split = line.split()
            lon = float(split[1])
            lat = float(split[2])
            value = float(split[3])

            i = int(0.002 + 1.0 * (lon - 180.0))
            j = int(0.002 + 3.0 * (lat - 30.0))        
            if i > iMax:
                iMax = i
            self.E[k][j][i] = value
            # Why are these lines here in Howard's commented code? 
            # They create errors in E-values
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
        print "i max is", iMax
        return

    ''' Methods called by the user interacting with the GUI '''
    def setupSignals(self):
        self.setupInputParameterSignals()
        self.setupButtonPressSignals()
        return

    def setupButtonPressSignals(self):
        self.nextButton.clicked.connect(self.nextButtonClicked)
        return

    def setupInputParameterSignals(self):
        self.sampleWindowBox.editingFinished.connect(
            self.sampleWindowBoxEditingFinished)
        self.pressureCutOffBox.editingFinished.connect(
            self.pressureCutOffBoxEditingFinished)
        self.maxInterpDepthBox.editingFinished.connect(
            self.maxInterpDepthBoxEditingFinished)
        self.stepSizeBox.editingFinished.connect(
            self.stepSizeBoxEditingFinished)
        self.dynHeightAtPBox.editingFinished.connect(
            self.dynHeightAtPBoxEditingFinished)
        self.relativeToPrefBox.editingFinished.connect(
            self.relativeToPrefBoxEditingFinished)
        self.totalModesBox.editingFinished.connect(
            self.totalModesBoxEditingFinished)
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
        return

    def maxInterpDepthBoxEditingFinished(self):
        self.maxInterpDepth = self.maxInterpDepthBox.value()
        self.updateNPress()
        return

    def stepSizeBoxEditingFinished(self):
        self.stepSize = self.stepSizeBox.value()
        self.updateNPress()
        return

    def relativeToPrefBoxEditingFinished(self):
        self.relativeToPref = self.relativeToPrefBox.value()
        return

    def totalModesBoxEditingFinished(self):
        self.totalModes = self.totalModesBox.value()
        return

    def dynHeightAtPBoxEditingFinished(self):
        self.dynHeightAtP = self.dynHeightAtPBox.value()
        return

    def nextButtonClicked(self):
        if not self.CLASS_INITTED:
            print "Please wait for class variables to initialize"
            return
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

    ''' Every function call for the normal working of the program is made here.
    mainLoop is called at the start of the program and it does everything from
    collecting data, to plotting results.'''
    def mainLoop(self):
        floats = []    # 800 long
        # del floats[:]
        doSort = True
        julStart, julEnd = self.getJulianStartAndEnd()
        self.numFloats = 0
        for i in xrange(julStart, julEnd):
            yearMonthDayPath0 = self.checkFloatsFromIndex(i, self.path0, floats)
        if self.numFloats == 0:
            print ("There are no floats within the provided range. ",
                   "Terminating program")
            return
        if self.verbose:
            print "numFloats is ", self.numFloats
            print "len of floats is ", len(floats)
        iFloat = 0
        print "Processing data..."
        for flt in floats:
            numRecs = getProfile(flt, self.P, self.T, self.S)
            despike(self.T, numRecs)
            despike(self.S, numRecs)
            skip = self.checkSkip()
            if skip:
                continue
            self.dataForcing(numRecs)
            self.storeData(numRecs, iFloat)
            iFloat += 1
        print "Performing various calculations..."
        # Specific volume and Dynamic Heights
        # Replaced Deltap with self.stepSize
        self.specificVolAndDynH()
        self.meanAndStdDev()
        if doSort:
            self.removeDuplicates(floats) 
        dynHeightBar, dynHeightVar = self.recomputeMeanAndVariation()
        self.subtractMeanFromStored(dynHeightBar)
        self.bestFitMode(dynHeightVar)
        self.mapCirculations(dynHeightBar)
        divSl = self.findDivingStreamline()
        xPLL, xPLH, yPLL, yPLH = self.collectPlotData(divSl)
        plotLats, plotLons, plotData = self.saveOutput(dynHeightBar)
        self.contour(xPLL, xPLH, yPLL, yPLH, plotLats, plotLons, plotData)
        return

    def determineOutputFilename(self):
        # self.outPath
        outFileName = ("Dh_" + 
                       str(self.plotCentre) + "_" + 
                       str(self.sampleWindow) + "_" + 
                       str(self.totalModes) + "_" +  
                       str(self.maxInterpDepth))
        infoFile = open((self.outPath + "CirculationInfo.txt"), 'a+')
        amount = self.checkForFileName(outFileName, infoFile)
        if amount > 1:
            outFileName += "_" + str(amount)
        outFileName += ".csv"
        infoFile.write(outFileName + " | " + 
                       str(self.pressureCutOff) + " | " + 
                       str(self.stepSize) + " | " + 
                       str(self.dynHeightsAtP) + " | " + 
                       str(self.relativeToPref) + " | " + "\n")
        return self.outPath + outFileName

    def checkForFileName(self, outFileName, infoFile):
        amount = 1
        for line in infoFile:
            if outFileName in line:
                amount += 1
        return amount    

    ''' Writes results to .csv file and creates storage arrays for data to be 
    plotted.'''
    def saveOutput(self, dynHeightBar):
        plotLats = []

        minLat = 90
        maxLat = -90
        minLon = 360
        maxLon = -180

        outFileName = self.determineOutputFilename()
        outFile = open(outFileName, 'w')

        for i in xrange(0, self.nX):
            lonC = i - 180.
            for j in xrange(0, self.nY):
                if not np.isnan(self.dHFit[j][i]):
                    latC = 30. + j / 3.
                    outString = (str(lonC) + "," + str(latC) + "," + 
                                 str(self.dHFit[j][i]) + '\n')
                    outFile.write(outString)
                    # Recording the latitude saved vals to use for 
                    if len(plotLats) == 0:
                        plotLats.append(latC)
                    elif plotLats[-1] < latC:
                        plotLats.append(latC)
                    if latC < minLat:
                        minLat = latC
                    if latC > maxLat:
                        maxLat = latC
                    if lonC < minLon:
                        minLon = lonC
                    if lonC > maxLon:
                        maxLon = lonC
        outFile.write("-116,60.0," + str(dynHeightBar) + '\n')
        latDiff = ceil(abs(maxLat - minLat) * 3)
        lonDiff = int(abs(maxLon - minLon))
        plotLats = np.asarray(plotLats)
        plotLons = np.arange(minLon, maxLon + 1)

        plotData = self.savePlotData(len(plotLons), len(plotLats))        
        plotLats, plotLons = self.convertTo2D(plotLats, plotLons)

        outFile.close()
        return plotLats, plotLons, plotData

    ''' Pulls formats results in a way that there are no empty/wasted lines '''
    def savePlotData(self, dimX, dimY):
        jPlot = 0
        iPlot = 0
        iNext = False
        plotData = np.zeros((dimX, dimY))
        for i in xrange(0, self.nX - 1):
            iNext = True
            lonC = i - 180.
            for j in xrange(0, self.nY):
                if not np.isnan(self.dHFit[j][i]):
                    latC = 30. + j / 3.
                    plotData[i][j] = self.dHFit[j][i]
                    jPlot += 1
                    if iNext:
                        iPlot += 1
                        iNext = False
            jPlot = 0
        return plotData

    ''' Converts 1D latitude/longitude arrays into 2D arrays with longitude 
    iterating across the 2nd array, and latitude iterating across the 1st array
    '''
    def convertTo2D(self, plotLats, plotLons):
        plotLatsFix = np.zeros((len(plotLats), len(plotLons)))
        plotLonsFix = np.zeros((len(plotLats), len(plotLons)))
        for i in xrange(0, len(plotLats)):
            for j in xrange(0, len(plotLons)):
                plotLatsFix[i][j] = plotLats[i]

        for i in xrange(0, len(plotLats)):
            for j in xrange(0, len(plotLons)):
                plotLonsFix[i][j] = plotLons[j]
        return plotLatsFix, plotLonsFix

    def collectPlotData(self, divSl):
        xP = np.empty((4))
        yP = np.empty((4))
        xPListLow = []
        xPListHigh = []
        yPListLow = []
        yPListHigh = []

        for iDiv in xrange(0, 64):
            x0 = self.dX * (iDiv - 1.)
            for jDiv in xrange(0, 90):
                y0 = self.dY * (jDiv - 1.)
                x1 = self.dHFit[jDiv][iDiv]
                x2 = self.dHFit[jDiv+1][iDiv]
                x3 = self.dHFit[jDiv+1][iDiv+1]
                x4 = self.dHFit[jDiv][iDiv+1]
                xMax = max(x1,x2,x3,x4)
                xMin = min(x1,x2,x3,x4)
                if(xMax > 90. or xMax < divSl or xMin > divSl):
                    continue
                # Contour passes through the box
                iP = 0

                xMax = max(x1, x2)
                xMin = min(x1, x2)
                if(divSl >= xMin and divSl <= xMax):
                    xP[iP] = x0
                    yP[iP] = y0 + self.dY * (divSl - x1) / (x2 - x1)
                    iP = iP + 1
                xMax = max(x2,x3)
                xMin = min(x2,x3)
                if(divSl>=xMin and divSl<=xMax):
                    xP[iP] = x0 + self.dX * (divSl - x2) / (x3 - x2)
                    yP[iP] = y0 + self.dY
                    iP=iP+1
                xMax = max(x3,x4)
                xMin = min(x3,x4)
                if(divSl >= xMin and divSl <= xMax):
                    xP[iP] = x0 + self.dX
                    yP[iP] = y0 + self.dY * (divSl - x4) / (x3 - x4)
                    iP = iP + 1
                xMax = max(x1,x4)
                xMin = min(x1,x4)
                if(divSl >= xMin and divSl <= xMax):
                    xP[iP] = x0 + self.dX * (divSl - x1) / (x4 - x1)
                    yP[iP] = y0
                    iP = iP + 1
                xPListLow.append(xP[0])
                xPListHigh.append(xP[1])
                yPListLow.append(yP[0])                
                yPListHigh.append(yP[1])
        return xPListLow, xPListHigh, yPListLow, yPListHigh 

    ''' Calculates the dividing streamline '''
    def findDivingStreamline(self):
        # Now search for dividing streamline
        divSl = 0.
        knt = 0
        # Was 58, changed to 59. Change is valid
        for dLat in xrange(40, 59):
            j = int(1.01 + 3. * (dLat - 30.))
            for iM in xrange(0, 60):
                i = 64 - iM
                if np.isnan(self.dHFit[j][i]): 
                    continue
                dynHeight1 = self.dHFit[j][i]
                break
            # Have dynHeight1
            divSl = divSl + dynHeight1
            knt = knt + 1
        divSl = divSl / knt
        return divSl

    ''' Calls functions to create the Basemap plot, the contour lines on top, 
    and the red dots representing ARGO floats. '''
    def contour(self, xPLL, xPLH, yPLL, yPLH, plotLats, plotLons, plotData):
        m = self.plotBasemap()
        self.plotFloatDots(m)
        self.plotContour(m, plotLats, plotLons, plotData)
        plt.show()
        return

    ''' Adds contour lines to the basemap object'''
    def plotContour(self, m, plotLats, plotLons, plotData):
        mX, mY = m(plotLons, plotLats)
        try:
            cs = m.contour(mX, mY, plotData.T)
        except Exception, e:
            print "ERROR: Invalid parameters for contour plotting"
            raise e
        return

    ''' Adds dots for the positions of each float that dat was collected form'''
    def plotFloatDots(self, m):        
        latToDel = []
        lonToDel = []
        # Remove all values that do not fit within existing lat/lon possibilites
        for i in xrange(0, self.Lat.size):
            if self.Lat[i] == 0 or self.Lat[i] > 360:
                latToDel.append(i)                
        for i in xrange(0, self.Lon.size):
            if self.Lon[i] == 0 or self.Lon[i] > 360:
                lonToDel.append(i)
        self.Lat = np.delete(self.Lat, latToDel)
        self.Lon = np.delete(self.Lon, lonToDel)

        mapLon, mapLat = m(self.Lon, self.Lat) 

        if self.verbose:
            print "NOR Lon Max ", np.amax(self.Lon), " Lon Min ", np.amin(self.Lon)
            print "NOR Lat Max ", np.amax(self.Lat), " Lat Min ", np.amin(self.Lat)

            print "MAP Lon Max ", np.amax(mapLon), " Lon Min ", np.amin(mapLon)
            print "MAP Lat Max ", np.amax(mapLat), " Lat Min ", np.amin(mapLat)

            axes = plt.gca()
            print "X axis intervals ", axes.xaxis.get_view_interval()
            print "Y axis intervals ", axes.yaxis.get_view_interval()

        plt.plot(mapLon, mapLat, 'ro')
        return 

    ''' Creates the basemap plot object. Many of the parameters are editable,
    so if the user is interested look at different projections to start. '''
    def plotBasemap(self):
        # setup polyconic basemap
        # by specifying lat/lon corners and central point.
        # area_thresh=1000 means don't plot coastline features less
        # than 1000 km^2 in area.
        latLeft = self.firstLatitude
        lonLeft = self.firstLongitude
        latRight = self.secondLatitude
        lonRight = self.secondLongitude

        m = Basemap(llcrnrlon=lonLeft,llcrnrlat=latLeft,urcrnrlon=lonRight, \
                    urcrnrlat=latRight, resolution='l',area_thresh=1000., \
                    projection='poly',lat_0=50,lon_0=-140)

        m.drawcoastlines()
        m.fillcontinents(color='darksage',lake_color='royalblue')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))
        m.drawmapboundary(fill_color='royalblue')
        plt.title("Circulation Data on the Southern Alaskan Coast")

        if self.verbose:
            print ("BASEMAP PLOT lonLeft/latLeft/lonRight/latRight ", 
                lonLeft, latLeft, lonRight, latRight)
        return m

    ''' Saves final dHFit results. '''
    def mapCirculations(self, dynHeightBar):
        # Now add mean and modes
        self.dHFit.fill(dynHeightBar)
        for i in xrange(0, self.nX):
            for j in xrange(0, self.nY):
                if not np.isnan(self.E[0][j][i]):
                    for k in xrange(0, self.totalModes): 
                        self.dHFit[j][i] += (self.evals[k] * self.E[k][j][i])
                else:
                    self.dHFit[j][i] = np.nan
        return

    ''' Various calculations. Edits values of self.Dhf, self.Dh, self.evals. '''
    def bestFitMode(self, dynHeightVar):
        frVarLast = 1.
        for qMode in xrange(0, self.totalModes):
            if qMode > 3.5:
                mode = qMode
            if qMode < 3.5:
                mode = 4 - qMode
            mode = qMode
            # Find best fit between Mode M and the data
            numEnt = 0
            for iFloat in xrange(0, self.numFloats):
                if self.accept[iFloat]:
                    eF = self.interpolate(self.E, mode, self.Lat[iFloat], self.Lon[iFloat])
                    self.Dhf[iFloat] = eF
                    numEnt = numEnt + 1
            q1 = 0.
            q2 = 0.
            nTemp = 0
            for iFloat in xrange(0, self.numFloats):
                if not np.isnan(self.Dhf[iFloat]) and self.accept[iFloat]:
                    q1 += self.Dh[iFloat] * self.Dhf[iFloat]
                    q2 += self.Dhf[iFloat] * self.Dhf[iFloat]
                    nTemp = nTemp + 1
            self.evals[mode] = q1 / q2
            rVar = 0.
            for iFloat in xrange(0, self.numFloats): 
                if not np.isnan(self.Dhf[iFloat]) and self.accept[iFloat]:
                    self.Dh[iFloat] = (self.Dh[iFloat] - 
                                     self.evals[mode] * self.Dhf[iFloat])
                    rVar = rVar + self.Dh[iFloat] * self.Dh[iFloat]
            rVar = rVar / nTemp
            frVar = rVar / dynHeightVar
            Contvar = frVarLast - frVar
            frVarLast = frVar
        return

    ''' Calculates the e-value. '''
    def interpolate(self, E, mode, Lati, Loni):
        eVal = np.nan
        y = 1. + 3. * (Lati - 30.)
        x = 1. + 1. * (Loni - 180.)
        iX = int(x + .0001)
        iY = int(y + .0001)
        rX = x - iX
        rY = y - iY
        
        # Outside bounds, so ignore
        if iX < 1 or iX >= self.nX:
            print "returning early X"
            return eVal
        if iY < 1 or iY >= self.nY:
            print "returning early Y"
            return eVal
        # Integer lat/lon, so easy
        if np.abs(rX) < .01 and np.abs(rY) < .01:
            eVal = self.E[mode][iY][iX]
            print "Returning 3rd if eVal", eVal
            return eVal       
        i1 = iX - 1
        i2 = iX + 1
        j1 = iY - 1
        j2 = iY + 1
        
        sumWeight = 0.
        sumWeightE = 0.
        a2 = .25
        for i in xrange(i1, i2):
            if i < 1 or i > self.nX:
                continue
            for j in xrange(j1, j2):
                if j < 1 or j > self.nY or np.isnan(self.E[mode][j][i]):
                    continue
                rho2 = (i - x) * (i - x) + (j - y) * (j - y)
                weight = np.square(-rho2 / a2)
                sumWeight += weight
                sumWeightE += self.E[mode][j][i] * weight
        if sumWeight > .2:
            eVal = sumWeightE / sumWeight
        return eVal

    def subtractMeanFromStored(self, dynHeightBar):
        for iFloat in xrange(0, self.numFloats):
            if self.accept[iFloat]:
                self.Dh[iFloat] = self.Dh[iFloat] - dynHeightBar
        return

    def recomputeMeanAndVariation(self):
        dynHeightBar = 0
        dynHeightVar = 0
        mFloats = 0
        for iFloat in xrange(0,  self.numFloats):
            if self.accept[iFloat]:
                dynHeightBar = dynHeightBar + self.Dh[iFloat]
                dynHeightVar = dynHeightVar + self.Dh[iFloat] * self.Dh[iFloat]
                mFloats = mFloats + 1
        dynHeightBar = dynHeightBar / mFloats
        dynHeightVar = dynHeightVar / mFloats
        dynHeightVar = dynHeightVar - dynHeightBar * dynHeightBar
        if self.verbose:
            print "Mean dynamic height ", dynHeightBar
            print "Variance in dynamic Height ", dynHeightVar
        return dynHeightBar, dynHeightVar

    ''' The goal of this functuon is to remove the listings for floats that 
    are very close together. This would help avoid map clutter. This function
    is still in testing due to its non-crucial nature. '''
    def removeDuplicates(self, floats):
        # Now have a list of good float observations, bad ones have Accpt(i)=-1
        if self.verbose:
            print "Checking for duplicates"
        # Sort entries so that there is only one entry per float
        # First see if there are any duplicates
        for iFloat in xrange(0, self.numFloats - 1):
            if not self.accept[iFloat]: 
                continue
            flt = floats[iFloat]
            splitName = re.split('[\W_]+', flt)
            # Gets 2nd last element from split string. Last element is "IOS"
            endOfName = splitName[-2]
            for jFlt in xrange(iFloat + 1, self.numFloats):
                if(not self.accept[jFlt]): 
                    continue
                # duplicate found
                if(floats[jFlt] == endOfName):
                    self.foundDuplicate(endOfName, iFloat, floats)
        return

    ''' This function is called when very similar floats are being used. '''
    def foundDuplicate(self, endOfName, dFlt, floats):
        # At least one duplicate is present for float dFlt
        dynHeightBar = 0
        avgLat = 0
        avgLon = 0
        Nav = 0
        for jFlt in xrange(dFlt, self.numFloats):
            if (floats[jFlt] == endOfName): 
                continue
            if (self.Dh[jFlt] < 0 or self.Dh[jFlt] > 5):
                self.accept[jFlt] = False
            if (not self.accept[jFlt]): 
                continue
            dynHeightBar = dynHeightBar + self.Dh[jFlt]
            avgLat = avgLat + self.Lat[jFlt]
            avgLon = avgLon + self.Lon[jFlt]
            Nav = Nav + 1
        if Nav > 0.5:
            self.Dh[dFlt] = dynHeightBar / Nav
            self.Lat[dFlt] = avgLat / Nav
            self.Lon[dFlt] = avgLon / Nav
        for jFlt in xrange(dFlt+1, self.numFloats):
            if (floats[jFlt] == endOfName): 
                self.accept[jFlt] = False
        return


    def meanAndStdDev(self):
        # Now compute mean and stnd dev and remove mean from dynHeight
        print "Computing mean and standard deviation"
        for i in xrange(0, 10):
            dynHeightBar = 0.
            dynHeightVar = 0.
            mFloats = 0
            nRej = 0
            for iFloat in xrange(0, self.numFloats):
                if self.accept[iFloat]:
                    dynHeightBar += self.Dh[iFloat]
                    dynHeightVar += self.Dh[iFloat] * self.Dh[iFloat]
                    mFloats = mFloats + 1
            dynHeightBar = dynHeightBar / mFloats
            dynHeightVar = dynHeightVar / mFloats
            dynHeightVar = dynHeightVar - dynHeightBar * dynHeightBar

            stDev = np.sqrt(dynHeightVar)
            devMax = 0
            for iFloat in xrange(1, self.numFloats):
                if self.accept[iFloat]:
                    devN = np.fabs((self.Dh[iFloat] - dynHeightBar) / stDev)
                    if (devN > devMax):
                        devMax = devN
                    if (devN > 3):
                        self.accept[iFloat] = False
                        if self.verbose:
                            "Rejecting float #", iFloat
                        nRej = nRej + 1
            if self.verbose:
                "Max deviation: ", devMax
            if(nRej < 0.5): 
                break
        return

    ''' Calculate specific volume and then dynamic height '''
    def specificVolAndDynH(self):
        iPressRef = 1 + self.relativeToPref / self.stepSize
        for iFloat in xrange(0, self.numFloats):
            for iPress in xrange(0, 500):
                # pressCount = self.stepSize * (iPress - 1.0)
                qTe = self.Te[iFloat][iPress]
                qSa = self.Sa[iFloat][iPress]
                if (qSa > 900):
                    print "Skipping from qSa"
                    continue
                if (qTe > 900):
                    print "Skipping from qTe"
                    continue
                qSt = getSigmaT(qSa, qTe)
                sigma, sVan = getSvanom(qSa, qTe, 0)
                self.Sa[iFloat][iPress] = sVan
            # Compute dynamic height at surface P = Pmap relative to 1000
            dHi = 0
            q = 0.5 * 1.0E-5
            iPress0 = int(0.01 + self.dynHeightsAtP / self.stepSize)
            iPress1 = iPress0 + 1

            for iPress in xrange(iPress1, iPressRef):
                dHi += (q * (self.Sa[iFloat][iPress] + 
                             self.Sa[iFloat][iPress-1]) * 
                        self.stepSize)
            self.Dh[iFloat] = dHi
        return

    ''' Save results into self.P, self.T, and self.S '''
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
                        self.Te[numFloat, iPr] = (
                            self.T[i-1] + rho * (self.T[i] - self.T[i-1]))
                    
                    self.Sa[numFloat, iPr] = 999.9
                    if (np.abs(self.S[i-1]) < 50. and np.abs(self.S[i]) < 50.): 
                        self.Sa[numFloat, iPr] = (
                            self.S[i-1] + rho * (self.S[i] - self.S[i-1]))
                # Finished interpolation to standard pressures for Float iFloat
            # Finished interpolation to standard pressures for all floats
        return

    def dataForcing(self, numRecs):
        # See if numRecs -1 overall works, since it's used as an index
        # 4950 ! If the sample is "near-surface" then force it to be surface
        if self.P[0] < 20:
            self.P[0] = 0

        # 4980 ! If the sample is almost deep enough, force it to be deep enough
        if self.P[numRecs] < self.maxInterpDepth:
            dTdP = (self.T[numRecs] - self.T[numRecs - 1] /
                    self.P[numRecs] - self.P[numRecs - 1])
            dSdP = (self.S[numRecs] - self.S[numRecs - 1] /
                    self.P[numRecs] - self.P[numRecs - 1])
            text = (self.T[numRecs] + dTdP * 
                   (self.maxInterpDepth - self.P[numRecs]))
            sext = (self.S[numRecs]  + dSdP *
                   (self.maxInterpDepth - self.P[numRecs]))            
            self.S[numRecs] = sext
            self.T[numRecs] = text
            self.P[numRecs] = self.maxInterpDepth
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
    def checkFloatsFromIndex(self, cycleJulDate, path0, floats):
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
        if self.verbose:
            print "inFileName is ", inFileName
        try:
            with open((inFileName), 'rb') as indexCSV:
                reader = csv.reader(indexCSV)
                # Skip the first entry since it's a header
                reader.next()
                for row in reader:
                    # Filer out the '*******' lines that are in some index files 
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
                        floats.append(yearMonthDayPath0 + row[0])
                        self.passedFloat(lat, lon)
                        self.numFloats += 1
        except (OSError, IOError), e:
            print "File was not found: ", inFileName
        return yearMonthDayPath0

    def passedFloat(self, lat, lon):
        self.Lat[self.numFloats] = lat
        self.Lon[self.numFloats] = lon
        # ToDo: A few lines left out here, looks completely unnecessary
        return
