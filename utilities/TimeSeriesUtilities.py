import datetime


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
    for i in xrange(0, 7306, (365 + int(yearAmounts))):
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




