'''
TimeSeries.py
 
TimeSeries is a program that calculates information at a given point based on
data available from other Argo floats nearby that transmitted recently.
 
TimeSeries uses NumPy external math libraries accessed by using the Python 3.5
distribution Anaconda.
 
Callum Colvine
Callum.Colvine@dfo-mpo.gc.ca
May 10th 2016
'''
 
# Line 45-65
# ToDo: This function. I don't have the file so I don't kno how to read it
# 370 Path0$=Drive$&"argo_mirror\pacific_ocean\"
# 380 Pathout$=Drive$&"argo_out\"
# 390 Sigpath$=Drive$&"projects\Sigma_Climate\"

 
import numpy as np
import sys 

# from os import unlink
import os
from datetime import time 

# Line 0-45
def initVars():
    # 10 OPTION BASE 1
    # 20 DIM Te(500,300),Sa(500,300),Lat(500),Lon(500),Accpt(500),Fltnm$(500)[50]
    Te = np.empty((500, 300))
    Sa = np.empty((500, 300))
    Lat = np.empty((500))
    Lon = np.empty((500))
    Accpt = np.empty((500))
    # ----- Is A$[50] a string of length 50? -----
    Fltnm = np.empty((500))
    # ----- Not sure ^^^ -----    
    Pdel = np.empty((1000))
    Sva = np.empty((1000))
    Dh = np.empty((1000))
    # 40 DIM S$[100],T$[100],U$[500],Path0$[80],Pathout$[80],Path$[80],R$[600],Closest_fltnm$[100]
    # Ignore all these declarations if they're strings, since I can do them at use-time
    Sa_av = np.empty(12)
    Sa_knt = np.empty(12)
    # 70 DIM Fil$[100],Filte$[100],Filsa$[100],Filst$[100],Out$[200],App$[200]
    # 80 DIM Ste$[500],Ssa$[500],Sst$[500],Ssp$[500],Sdh$[500],Shgt$[500],Pa$[500],Sclim$[500]
    # 90 DIM Location$[100],Param$[50],Lstmsg$[50]
    Sigref_st = np.empty((71))
    Sigref_te = np.empty((71))
    Sigref_sa = np.empty((71))
    Sigpath = np.empty((71))
    # 110! Dimensions of data to be read from files
    P = np.empty((2000))
    T = np.empty((2000))
    S = np.empty((2000))
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
    Te[:] = 999.9
    Sa[:] = 999.9    
    Lat[:] = 999.9
    Lon[:] = 999.9
    Accpt[:] = (-1)
    Sa_av[:] = 0
    Sa_knt[:] = 0
    # Dout backslash because in Python you need to escape special chearacters
    Drive = "E:\\"
    # 350 CALL Lastrun(Drive$,2)
    lastRun(Drive, 2)
    Path0 = Drive + "argo_mirror\\pacific_ocean\\"
    Pathout = Drive + "argo_out\\"
    Sigpath = Drive + "projects\\Sigma_Climate\\"
    # Made this file empty here when he writes 3x71 empties
    # if problems: change this maybe
    csvF = open((Sigpath + "Mp26_i.csv"), 'w')
    csvF.close()
    return 
''' Pulls settings from an already existing file. INCOMPLETE'''
def setSettings():
    # 470 ! Read previous settings
    # 480 ASSIGN @Pin TO Drive$&"argo_programs\TSlast.txt";FORMAT ON
    settingsF = open((Drive + "argo_programs\\TSlast.txt"), 'r+')
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

# Line 65-285
def userDefinedSettings():
    
    if valid:
        Rho_print = input("List floats within radius (km) = ... ")
    # 1160 CLEAR SCREEN
    clear()
    # 1170! Return Point
    # 1180 CALL today(Td)
    today(1)
    # 1190 Q$=Dc12$
    # 1200 Pc=POS(Q$,",")
    # 1210 Dc2=VAL(Q$[Pc+1])
    # 1220 Dc1=Dc2-30.
    # 1230 PRINT "Suggested last end day number "&Q$[Pc+1]
    # 1240 PRINT "Today is day number ";Td
    # 1250 PRINT "Earliest recommended start day is 190"
    # 1260 !
    # 1270 ASSIGN @Plst TO Drive$&"argo_out\lstmsge.csv";FORMAT ON
    # 1280 ENTER @Plst;R$
    # 1290 ASSIGN @Plst TO *
    # 1300 Last_final_day=VAL(R$[3,6])
    # 1310 PRINT "Final day last time program was run = ";Last_final_day
    # 1320!
    # 1330 Td5=5.0*INT(Td/5.0)
    # 1340 IF Last_final_day>Td5 THEN Td5=Last_final_day
    # 1350 IF Def_opt<.1 THEN
    # 1360 OUTPUT KBD;VAL$(Dc1)&","&VAL$(Td5)&"ÿ<ÿ<ÿ<ÿ<";
    # 1370 INPUT "Range of days (day numbers of 2001) = ",Dc1,Dc2
    # 1380 ELSE
    # 1390 Dc2=Td5
    # 1400 END IF
    # 1410!
    # 1420 IF Def_opt<.1 THEN
    # 1430 OUTPUT KBD;Dt$&"ÿ<ÿ<ÿ<ÿ<";
    # 1440 INPUT "Day step = Dt = (centre estimates on successive Dt days)...",Dt
    # 1450 ELSE
    # 1460 Dt=VAL(Dt$)
    # 1470 END IF
    # 1480!
    # 1490 IF Def_opt<.1 THEN
    # 1500 OUTPUT KBD;Twin$&"ÿ<ÿ<ÿ<ÿ<";
    # 1510 INPUT "Sample window i.e. include data centre date +/- (10 is usual)Twin days...",Twin
    # 1520 ELSE
    # 1530 Twin=VAL(Twin$)
    # 1540 END IF
    # 1550!
    # 1560 CLEAR SCREEN
    # 1570!
    # 1580 CALL Julday(1,Dc1,Dayf,Monf,Yearf)
    # 1590 CALL Julday(1,Dc2,Dayl,Monl,Yearl)
    # 1600 S$=VAL$(Dc1)&" = ("&VAL$(Dayf)&"/"&VAL$(Monf)&"/"&VAL$(Yearf)&") to "&VAL$(Dc2)&" = ("&VAL$(Dayl)&"/"&VAL$(Monl)&"/"&VAL$(Yearl)&")"
    # 1610 PRINT "Start and end dates are:- "&S$
    # 1620 PRINT "Julian day centre +/- window of "&VAL$(Twin)
    # 1630 PRINT "Output data at intervals of "&VAL$(Dt)&" days."
    # 1640!
    # 1650 IF Def_opt<.1 THEN
    # 1660 OUTPUT KBD;Lat$&"ÿ<ÿ<ÿ<ÿ<ÿ<";
    # 1670 INPUT "Accept data between Lat1 and Lat2....",Lat1,Lat2
    # 1680 ELSE
    # 1690 Pc=POS(Lat$,",")
    # 1700 Lat1=VAL(Lat$[1,Pc-1])
    # 1710 Lat2=VAL(Lat$[Pc+1])
    # 1720 END IF
    # 1730 !
    # 1740 IF ABS(Lat1)>90. OR ABS(Lat2)>90. THEN GOTO 1170
    # 1750 !
    # 1760 IF Def_opt<.1 THEN
    # 1770 OUTPUT KBD;Lon$&"ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<";
    # 1780 INPUT "Accept data between Lon1 and Lon2....",Lon1,Lon2
    # 1790 ELSE
    # 1800 Pc=POS(Lon$,",")
    # 1810 Lon1=VAL(Lon$[1,Pc-1])
    # 1820 Lon2=VAL(Lon$[Pc+1])
    # 1830 END IF
    # 1840 !
    # 1850 IF ABS(Lon1)>360.001 OR ABS(Lon2)>360.001 THEN GOTO 1170
    # 1860 IF Lon1<-.000001 THEN Lon1=Lon1+360
    # 1870 IF Lon2<-.000001 THEN Lon2=Lon2+360
    # 1880!
    # 1890 IF Lat1>Lat2 THEN
    # 1900 Latq=Lat1
    # 1910 Lat1=Lat2
    # 1920 Lat2=Latq
    # 1930 END IF
    # 1940!
    # 1950 IF Lon1>Lon2 THEN
    # 1960 Lonq=Lon1
    # 1970 Lon1=Lon2
    # 1980 Lon2=Lonq
    # 1990 END IF
    # 2000!
    # 2010 PRINT "Latitude range for accepting data = "&VAL$(Lat1)&"°N - "&VAL$(Lat2)&"°N"
    # 2020 PRINT "Longitude range for accepting data = ";VAL$(Lon1)&"°E, "&VAL$(Lon2)&"°E = ("&VAL$(360.-Lon1)&"°W, "&VAL$(360-Lon2)&"°W)"
    # 2030!
    # 2040 IF Def_opt<.1 THEN
    # 2050 OUTPUT KBD;Pcut$&"ÿ<ÿ<ÿ<ÿ<";
    # 2060 INPUT "Reject floats with maximum pressure less then Pcutoff = ",Pcutoff
    # 2070 ELSE
    # 2080 Pcutoff=VAL(Pcut$)
    # 2090 END IF
    # 2100 PRINT "Reject floats that don't sample deeper than "&VAL$(Pcutoff)&" dbar"
    # 2110!
    # 2120 IF Def_opt<.1 THEN
    # 2130 OUTPUT KBD;Pmax$&"ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<";
    # 2140 INPUT "Interpolate to Pmax with steps Delta-P...Pmax,DeltaP = ",Pmax,Deltap
    # 2150 ELSE
    # 2160 Pc=POS(Pmax$,",")
    # 2170 Pmax=VAL(Pmax$[1,Pc-1])
    # 2180 Deltap=VAL(Pmax$[Pc+1])
    # 2190 END IF
    # 2200 !
    # 2210 PRINT "Interpolate vertically to Pmax = "&VAL$(Pmax)&"dbar with steps of DeltaP = "&VAL$(Deltap)&" dbar."
    # 2220 !
    # 2230 IF Pmax<0. OR Deltap<0. THEN GOTO 1170
    # 2240 Npress=1+INT(.001+Pmax/Deltap)
    # 2250 !
    # 2260 IF Def_opt<.1 THEN
    # 2270 OUTPUT KBD;Param$&"ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<";
    # 2280 INPUT "Enter 1 (else 0) to output Temp, Salt, Sig-t, Spiciness, Dyn Hgt......",Te_opt,Sa_opt,St_opt,Sp_opt,Dh_opt
    # 2290 ELSE
    # 2300 Pc1=POS(Param$,",")
    # 2310 Param$[Pc1,Pc1]=";"
    # 2320 Pc2=POS(Param$,",")
    # 2330 Param$[Pc2,Pc2]=";"
    # 2340 Pc3=POS(Param$,",")
    # 2350 Param$[Pc3,Pc3]=";"
    # 2360 Pc4=POS(Param$,",")
    # 2370 Param$[Pc4,Pc4]=";"
    # 2380 Te_opt=VAL(Param$[1,Pc1-1])
    # 2390 Sa_opt=VAL(Param$[Pc1+1,Pc2-1])
    # 2400 St_opt=VAL(Param$[Pc2+1,Pc3-1])
    # 2410 Sp_opt=VAL(Param$[Pc3+1,Pc4-1])
    # 2420 Dh_opt=VAL(Param$[Pc4+1])
    # 2430 END IF
    # 2440 !
    # 2450 Out$="Output files are:- "
    # 2460 IF Te_opt>.5 THEN Out$=Out$&"Temp"
    # 2470 IF Sa_opt>.5 THEN Out$=Out$&", Salt"
    # 2480 IF St_opt>.5 THEN Out$=Out$&", Sigt"
    # 2490 IF Sp_opt>.5 THEN Out$=Out$&", Spice"
    # 2500 IF Dh_opt>.5 THEN Out$=Out$&", Dyn. Hgt."
    # 2510 PRINT Out$
    # 2520 !
    # 2530 IF Def_opt<.1 THEN
    # 2540 OUTPUT KBD;Location$&"ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<ÿ<";
    # 2550 INPUT "Enter latitude and longitude of the station to interpolate to.....",Latc,Lonc
    # 2560 ELSE
    # 2570 Pc=POS(Location$,",")
    # 2580 Latc=VAL(Location$[1,Pc-1])
    # 2590 Lonc=VAL(Location$[Pc+1])
    # 2600 END IF
    # 2610 !
    # 2620 PRINT "Will interpolate to the station:- "&TRIM$(VAL$(Latc))&"°N, "&TRIM$(VAL$(Lonc))&"°E ("&TRIM$(VAL$(360.-Lonc))&"°W)"
    # 2630 !
    # 2640 DEG
    # 2650 Scx0=Scy*COS(Latc)
    # 2660 !
    # 2670 IF Def_opt<.1 THEN
    # 2680 !IF ABS(Dc1-Dc2)<151. THEN OUTPUT KBD;"1ÿ<";
    # 2690 !IF ABS(Dc1-Dc2)>151. THEN OUTPUT KBD;"0ÿ<";
    # 2700 OUTPUT KBD;"1ÿ<";
    # 2710 INPUT "Enter 1 to append to existing files, otherwise 0....",Append_opt
    # 2720 ELSE
    # 2730 Append_opt=1
    # 2740 END IF
    # 2750 !
    # 2760 IF Def_opt<.1 THEN
    # 2770 OUTPUT KBD;"0ÿ<";
    # 2780 INPUT "Enter 1 for verbose option....",Verb_opt
    # 2790 ELSE
    # 2800 Verb_opt=0
    # 2810 END IF
    # 2820 !
    # 2830 IF Def_opt<.1 THEN
    # 2840 OUTPUT KBD;"1ÿ<";
    # 2850 INPUT "Enter 1 if all entries are OK, 0 to return....",Ret_opt
    # 2860 ELSE
    # 2870 Ret_opt=1
    # 2880 END IF
    # 2890 IF ABS(Ret_opt-1.)>.1 THEN GOTO 1170
    # 2900 !
    return 

# Line 285-308
def saveUserSetting():
 
    return
# Line 308-321
def defineOutputFiles():
 
    return 
# Line 321-382
def purgeAndAssign():
 
    return 
# Line 382-583
def prepareCSV():
 
    return 
# Line 583-752
def interpolations():
              
    return  
# Line 752-814
def interpolatePressures():
              
    return  
# Line 814-905
def interpolateStation():

    return  
# Line 905-1010
def computeDHgtRelativeToPmax():
              
    return 
# Line 1010-1040
def interpolationsComplete():
 
    return  
# Line 1042-1259
# (Opt,Jd,Dom,Mon,Year)
def julDateToDate():
 
    return  
# Line 1042-1259
# (Opt,Jd,Dom,Mon,Year)
def dateToJulDate():
 
    return  
# Line 1262-1276
def sigmaT():

    return  
# Line 1278-1286
def despike():

    return  
# 1288-1312
def spiciness():
 
    return  
# Line 1314-1371
def svanom():

    return  
# Line 1373-1400
def today(Td):
    print(time())
    return 
# Line 1402-1623
def getProfile():
 
    return  

def clear():
    os.system('cls' if os.name=='nt' else 'clear')
    return

# Line 1625-1633
''' Lastrun (I think) defines if the program was a success or not. 
Starts with 2 being written, then writes -1 after "Normal End"-ing'''
def lastRun(curDrive, status):
    # Overwrite the previous "lastrun.txt" file
    lastRunF = open((curDrive + r"projects\argo\status\lastrun.txt"), 'w')
    # Typecast status to a string because write() uses strings
    lastRunF.write(str(status))
    lastRunF.close()
    return 
# Line 1635-1687
def sanityCheck():
 
    return 
# Line 1689-1748
def log():
              
    return  
def main():
    print ("Hello world")     
    initVars()   
    # setSettings()
    userDefinedSettings()

 
if (__name__ == "__main__"):
    main()