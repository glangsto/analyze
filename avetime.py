# python code to average times and arrays of times
# Glen Langston
# HISTORY 
# 180106 GIL first version

from datetime import datetime
import numpy as np

def ave2datetimes( datetime1, datetime2):
    """
    Return the average fo two datetime values
    """
    t3 = np.datetime64(datetime1)
    t4 = np.datetime64(datetime2)
    dt = t4-t3
#    print dt
    outdatetime = t3 + (dt/2.)

    return outdatetime

def avedatetimes( datetimes):
    """
    Return the average fo two datetime values
    """
    count = 0
    diff = np.timedelta64( '0','s')

    for adate in datetimes:
        if count == 0:
            t1 = np.datetime64( adate)
#            print 'First: ', t1
            count = 1
        else:
            t2 = np.datetime64( adate)
            diff = diff + (t2 - t1)
            count = count + 1

    outtime = t1 + (diff/float(count))

    return outtime

if __name__=='__main__':
    time1 = "2017-12-30T16:33:35.177875"
    time2 = "2017-12-30T16:37:35.372924"
    time3 = "2017-12-30T16:48:35.372924"
    print 'Inputs:  ', time1,time2

    outtime = ave2datetimes( time1, time2)
    print 'Average: ',outtime

    # test second version
    times = [ time1, time2, time3]
    print 'Inputs:  ', times
    outtime2 = avedatetimes( times)
    print 'Average: ', outtime2

    

