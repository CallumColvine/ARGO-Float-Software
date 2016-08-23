import datetime
import numpy as np
import time

def formatToDateTime(year, month, day):
    stringVersion = str(day) + '.' + str(month) + '.' + str(year)
    stringFormat = '%d.%m.%Y'
    dateTuple = datetime.datetime.strptime(stringVersion,
                                           stringFormat)
    return dateTuple

def julianToDate(julDate):
    yearAmounts = 0
    curYear = 0
    daysIn = 0
    # 7306 is the Julian Date for 2020
    # 14612 is the Howard Julian date for 2040
    for i in xrange(0, 14612, (365 + int(yearAmounts))):
        if (i <= julDate) and ((i + 365 + int(yearAmounts)) > julDate):
            daysIn = julDate - i
            break
        curYear += 1
        yearAmounts += 0.25
        if yearAmounts > 1:
            yearAmounts -= 1
    year = curYear + 2001
    result = datetime.datetime(year, 1, 1) + datetime.timedelta(daysIn)
    month = result.month
    day = result.day
    return day, month, year


''' Converts a year and the day of that year to Howard's Julian date version.
NOTE: Does NOT use day, or month to calculate the Julian date. These are passed 
in for debugging purposes. It only uses year, and dayOfYear. '''
def dateToJulian(day, month, year, dayOfYear):
    milleniaRemainder = year % 1000
    newYear = int((milleniaRemainder - 1) * 365)
    # Account for leap years
    newYear += ((year - 2001) / 4)
    newYear += dayOfYear
    return newYear


''' Lastrun defines if the program was a success or not. 
Starts with 2 being written, then writes -1 after "Normal End"-ing'''
def lastRun(curDrive, status):
    # Overwrite the previous "lastrun.txt" file
    lastRunF = open((curDrive + r"projects\lastrun.txt"), 'w+')
    # Typecast status to a string because write() uses strings
    lastRunF.write(str(status))
    lastRunF.close()
    return

def todayInJulian():
    day = int(time.strftime("%d"))
    month = int(time.strftime("%m"))
    year = int(time.strftime("%Y"))
    dayOfYear = int(time.strftime("%j"))
    todayInJul = dateToJulian(day, month, year, dayOfYear)
    return todayInJul


def todayInDate():
    curDay = int(time.strftime("%d"))
    curMonth = int(time.strftime("%m"))
    curYear = int(time.strftime("%Y"))
    return curDay, curMonth, curYear



''' Heavy calculations originated from Howard's code. Will not be converted
to naming standards '''
def getSigmaT(S, T):
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

# Variable names conflict with standard due to direct port of code of 
# Howard's HT Basic.
def getSvanom(S, T, P0):
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
def getSpiciness(temp, salt):
    # Hardcoded by Howard
    b = np.matrix([[0.0, .77442, -.00585, .000984, -.000206],
                   [.051665, .002034, -.0002745, -.0000085, .0000136],
                   [-6.64783E-3, -2.4681E-4, -1.428E-5, 3.337E-5, 7.894E-6],
                   [-5.4023E-5, 7.326E-6, 7.0036E-6, -3.0412E-6, -1.0853E-6],
                   [3.949E-7, -3.029E-8, -3.8209E-7, 1.0012E-7, 4.7133E-8],
                   [-6.36E-10, -1.309E-9, 6.048E-9, -1.1409E-9, -6.676E-10]])
    spice = 0
    sp = salt - 35
    theta = temp
    # Reversed I and J here, are arrays backwards in HP Basic?
    for i in range(0, 5):
        for j in range(0, 4):
            spice = spice + b[i, j] * (np.power(theta, i)) * (np.power(sp, j))
    return spice


''' Called on float numbers that have been found in the index file. Get 
profile will determine the order in which and fill up the P, T, S arrays '''
def getProfile(floatPath, P, T, S):
    order = []
    numChannels = 0
    numRecords = 0
    fileTempQc = 0
    filePSalQc = 0
    numRecs = extractDataFromFloatFile(floatPath, order, P, T, S)
    if not checkIfReturn(numRecs, fileTempQc, filePSalQc):
        print "Returning profile prematurely due to bad values"
        return numRecs
    return numRecs

''' Quits dataset of float if it does not have all the required components '''
def checkIfReturn(numRecords, fileTempQc, filePSalQc):
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


''' Used by getProfile() to determine validity of float data '''
def extractDataFromFloatFile(path, order, P, T, S):
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
                numRecs = getNumRecs(line)
                # if numRecs > 500:
                #     print "Num recs > 500 in file: ", path
            order, pFound, tFound, sFound = \
                findOrder(line, order, pFound, tFound, sFound)
            if not (verifyUsability(line)):
                return
        else:
            i += 1
            appendInfo(line, order, i, P, T, S)
        if line.find("END OF HEADER") != -1:
            passedEndOfHeader = True
    return numRecs

def getNumRecs(line):
    split = line.split()
    # print "numRecs is ", int(split[4])
    return int(split[4])


''' Used to determine the order of data held inside each float file'''
# TODO: Improve runtime, reduce calls to this function
def findOrder(line, order, pFound, tFound, sFound):
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


''' Makes sure float has minimum requirements for data '''
def verifyUsability(line):
    if line.find("NUMBER OF RECORDS") != -1:
        if line.split()[2] < 15:
            return False
    if line.find("NUMBER OF CHANNELS") != -1:
        if line.split()[2] < 8:
            return False
    return True


''' Adds data pulled form individual float files to the P, T, S arrays '''
def appendInfo(line, order, i, P, T, S):
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
                P[i] = twoPos
            else:
                P[i] = zeroPos
        elif data == 'T':
            if np.abs(twoPos) < 9000:
                T[i] = twoPos
            else:
                T[i] = zeroPos
        elif data == 'S':
            if np.abs(twoPos) < 9000:
                S[i] = twoPos
            else:
                S[i] = zeroPos
        position += 1

    return



# ToDo: could be moved to the Utils file
''' Checks if the pressure array is monotonic'''
def checkPressureMonotonic(numRecs, P, T, S):
    for i in xrange(1, (numRecs - 2)):
        if ((P[i] > P[i - 1] and P[i] > P[i + 1]) or 
            (P[i] < P[i - 1] and P[i] < P[i + 1])):
            numRecs = removeIndexFromPTS(i, numRecs, P, T, S)
    return numRecs


def removeIndexFromPTS(i, numRecs, P, T, S):
    np.delete(P, i)
    np.delete(T, i)
    np.delete(S, i)
    numRecs -= 1
    return numRecs


def despike(T, N):
    Nm1 = N - 1
    for i in xrange(1, Nm1):
        Dt1 = T[i] - T[i - 1]
        Dt2 = T[i + 1] - T[i]
        Pr = Dt1 * Dt2
        if(Pr < 0. and np.abs(Dt1) > 3. and np.abs(Dt2) > 3.): 
            T[i]= .5 * (T[i + 1] + T[i - 1])
    return

def convertLatLonToNegative(latIn, lonIn):
    if latIn > 90:
        latIn -= 90
        latIn = -(latIn)
    if lonIn > 180:
        lonIn = lonIn % 180
        lonIn -= 180
        # lonIn = -(lonIn)
    return latIn, lonIn

# Maybe not necessary?
def fillgaps(T, Nrecs1):
    # 31630 Nrecs=Nrecs1
    # 31640 !
    # 31650 Ifirst=0
    # 31660 Line_ifirst: Ifirst=Ifirst+1
    # 31670 IF T(Ifirst)>900. THEN GOTO Line_ifirst
    # 31680 !
    # 31690 ! Ifirst is first good value
    # 31700 Nrecs=Nrecs+1
    # 31710 Line_nrecs: Nrecs=Nrecs-1
    # 31720 IF Nrecs<2 THEN Line_nogaps
    # 31730 IF T(Nrecs)>900. THEN GOTO Line_nrecs
    # 31740 !
    # 31750 ! Nrecs is the last good value
    # 31760 !
    # 31770 ! Now search for a gap
    # 31780 I1=Ifirst
    # 31790 !
    # 31800 Line_i1: I1=I1+1
    # 31810 IF I1>=Nrecs THEN GOTO Line_nogaps
    # 31820 IF T(I1)<900. THEN GOTO Line_i1
    # 31830 !
    # 31840 ! Have a gap, 1st bad value is at I1
    # 31850 I2=I1
    # 31860 Line_i2: I2=I2+1
    # 31870 IF T(I2)>900. THEN GOTO Line_i2
    # 31880 !
    # 31890 ! 1st good value after I1 is I2
    # 31900 !
    # 31910 I2=I2-1
    # 31920 Wid=I2-I1+1
    # 31930 Tdif=T(I2+1)-T(I1-1)
    # 31940 FOR I=I1 TO I2
    # 31950 T(I)=T(I-1)+Tdif/Wid
    # 31960 NEXT I
    # 31970 I1=I2
    # 31980 GOTO Line_i1
    # 31990 !
    # 32000 Line_nogaps: ! No gaps left
    return