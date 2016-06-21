import datetime


def formatToDateTime(year, month, day):
    stringVersion = str(day) + '.' + str(month) + '.' + str(year)
    stringFormat = '%d.%m.%Y'
    dateTuple = datetime.datetime.strptime(stringVersion,
                                           stringFormat)
    return dateTuple


# Should be working !!! NOT tested
def julianToDate(julDate):
    # print "INPUT julian date is ", julDate
    yearAmounts = 0
    curYear = 0
    daysIn = 0
    # 7306 is the Julian Date for 2020
    for i in xrange(0, 7306, (365 + int(yearAmounts))):
        if (i <= julDate) and ((i + 365 + int(yearAmounts)) > julDate):
            daysIn = julDate - i
            break
        curYear += 1
        yearAmounts += 0.25
        if yearAmounts > 1:
            yearAmounts -= 1
        # Move i to the next julian year
        # i += (yearAmounts + 365)
    year = curYear + 2001
    # print "daysIn are ", daysIn
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
    print "Returning the year: ", newYear, " from ", day, month, year
    # ----- Sandy -----
    # newYear = int((year - 2001) * 365.25) # <- Sandy's solution
    # ----- Howard -----
    # if (int(year == 2004) or int(year == 2008) or int(year == 2012) or 
    #     int(year == 2016) or int(year == 2020)):
    #     leap = 1
    # else: 
    #     leap = 0
    # if int(month) > 1:  newYear += 31
    # if int(month) > 2:  newYear += 28 + leap
    # if int(month) > 3:  newYear += 31
    # if int(month) > 4:  newYear += 30
    # if int(month) > 5:  newYear += 31
    # if int(month) > 6:  newYear += 30
    # if int(month) > 7:  newYear += 31
    # if int(month) > 8:  newYear += 31
    # if int(month) > 9:  newYear += 30
    # if int(month) > 10:  newYear += 31
    # if int(month) > 11:  newYear += 30
    # newYear += int(day)    
    # ----- Callum ----- 
    return newYear




