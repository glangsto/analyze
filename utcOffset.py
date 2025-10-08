# function to return the offset between utc and local time zone
# HISTORY
# 25Oct06 GIL initial version

from datetime import datetime
import tzlocal
import pytz

def strAst( utc):
    """
    Convert a datetime Utc value to the AST format time string
    """
    utcYYmmmDD = utc.strftime("%y%b%dT%H%M%S.%f")
    # only keep milliseconds
    nutc = len(utcYYmmmDD)
    utcYYmmmDD = utcYYmmmDD[0:(nutc-3)]
    # replace the space with a T
#    utcYYmmmDD = utcYYmmmDD.replace(' ', 'T', 1)
    return utcYYmmmDD
# end of strAst()

def getTimeZoneName():
    """
    returns the string name of your local time zone
    something like "America/New York" etc
    """
    # get local time zone
    now = datetime.now()
    local_now = now.astimezone()
    timeZoneName = tzlocal.get_localzone_name()
    return timeZoneName
# end of getTimeZoneName()

def getUtcOffsetTimeZone(date_str, timeZoneName):
    """
    Calculates the UTC offset for a given date string and timezone name.

    Args:
        date_str (str): The date in 'YYYY-MM-DD' format.
        timeZoneName (str): The name of the local timezone (e.g., 'America/New_York').

    Returns:
        float: The UTC offset in hours for the specified date and timezone.
    """
    # Create a naive datetime object from the date string
    naive_date = datetime.strptime(date_str, '%Y-%m-%d')

    # Get the timezone object
    tz = pytz.timezone(timeZoneName)

    # Localize the naive datetime to the specified timezone
    aware_datetime_local = tz.localize(naive_date)

    # Get the UTC offset as a timedelta object
    offset_timedelta = aware_datetime_local.utcoffset()

    # Convert the timedelta to hours
    offset_hours = offset_timedelta.total_seconds() / 3600

    return offset_hours

def getUtcOffset( dateStr):
    """
    getUtcOffset() gets the UTC offset in hours between your timezone and UTC for an input date
    where
    dateStr = date for which the offset is needed. 
    return utcOffsetHours.    Updates due to daylight savings time
    """
    timeZoneName = getTimeZoneName()
    
    offset = getUtcOffsetTimeZone(dateStr, timeZoneName)
    print(f"UTC offset for {dateStr} in {timeZoneName}: {offset} hours")
    return offset

def main():
    """
    Test function
    """
    # Example usage:
    dateStr = '2024-03-10'  # A date before DST starts in many regions
    timeZoneName = 'America/New_York'
    offset_march = getUtcOffsetTimeZone(dateStr, timeZoneName)
    print(f"UTC offset for {dateStr} in {timeZoneName}: {offset_march} hours")

    # now see what today's offset is
    now = datetime.now()
    dateStrNow = now.strftime("%Y-%m-%d")
    offset_now = getUtcOffset( dateStrNow)

    utcYYmmmDD = strAst( now)
    print("UTC in AST format: %s" % (utcYYmmmDD))
    
    # finally check a summer date
    dateStr_dst = '2024-06-15' # A date during DST
    offset_june = getUtcOffsetTimeZone(dateStr_dst, timeZoneName)
    print(f"UTC offset for {dateStr_dst} in {timeZoneName}: {offset_june} hours")
    return

if __name__ == "__main__":
    main()

