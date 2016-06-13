import datetime


def formatToDateTime(year, month, day):
    stringVersion = str(day) + '.' + str(month) + '.' + str(year)
    stringFormat = '%d.%m.%Y'
    dateTuple = datetime.datetime.strptime(stringVersion,
                                           stringFormat)
    return dateTuple


# Should be working !!! NOT tested
def julianToDate(julDate):
    print "INPUT julian date is ", julDate
    yearAmounts = 0.25
    curYear = 0
    daysIn = 0
    # 7306 is the Julian Date for 2020
    for i in xrange(0, 7306, int(yearAmounts + 365)):
        if (i <= julDate) and ((i + int(yearAmounts) + 365) > julDate):
            daysIn = julDate - i
            break
        curYear += 1
        yearAmounts += 0.25
        # Move i to the next julian year
        # i += (yearAmounts + 365)
    year = curYear + 2001
    print "daysIn are ", daysIn
    result = datetime.datetime(year, 1, 1) + datetime.timedelta(daysIn)
    month = result.month
    day = result.day
    return day, month, year

def dateToJulian(day, month, year, dayOfYear):
    milleniaRemainder = year % 1000
    newYear = (milleniaRemainder - 1) * 365
    # Account for leap years
    newYear += ((year - 2001) / 4)
    newYear += dayOfYear
    return newYear

