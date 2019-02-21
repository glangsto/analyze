# Function to compute the mean, average deviation, standard deviation
# variance, skew and kurtosis
# grid with a convolving function
# HISTORY
# 17JAN28 GIL finish initial version

import numpy as np

def moment(  xs, ys):
    """
    moment returns the ave, average deviation, standard deviation, variance, skew and curtosis) of an input array
    """

#void moment( data, int n, float *ave, float *adev, float *sdev, float *var, float *skew, float *curt)
#Given an array of data[1..n], this routine returns its mean ave, average deviation adev, standard deviation sdev, variance var, skewness skew, and kurtosis curt.
    n = len( ys)

    ave = 0
    adev = 0
    sdev = 0
    var = 0
    skew = 0
    curt = 0
    sum = 0
    if (n <= 2):
        print "Array length n must be at least 2 in moment"
        return (ave,adev,sdev,var,skew,curt)

# First pass to get the mean.
    for iii in range(n):
        sum += ys[iii]
    ave=sum/n;
    ep = 0

# Second pass to get the first (absolute), second, third, and fourth moments of the deviation from the mean.
    for iii in range(n):
        s = ys[iii]-ave
        adev = adev + np.fabs(s)
        ep = ep + s
        p = s * s
        var = var + p
        p = p * s
        skew = skew + p
        p = p * s
        curt = curt + p

#Put the pieces together according to the conventional definitions.
    adev = adev / n
    var  = (var-(ep*ep)/n)/(n-1); 
    sdev = np.sqrt(var);
    if var > 0:
        skew = skew / (n * var * sdev)

# Corrected two-pass formula.
    curt = curt / ((n*var*var)-3.0)
    return (ave,adev,sdev,var,skew,curt)

def test_main():
    """
    Test the Moment code
    """
    xs = np.linspace(-10, 10, 20)
    ys = np.linspace(-10, 10, 20)
    ys[5] = 0.    
    ys[15] = 0.

    # test the code
    (ave,adev,sdev,var,skew,curt) = moment(xs, ys)
    # display the inputs

    print 'ys[0-4]  : ',ys[0:5]
    print 'ys[5-9]  : ',ys[5:10]
    print 'ys[10-14]: ',ys[10:15]
    print 'ys[15-19]: ',ys[15:20]

    # display the outputs
    print 'ave,adev,sdev: ', ave,adev,sdev
    print 'var,skew,curt: ', var,skew,curt

if __name__ == "__main__":
    """
    If run as a standalone program, run the test
    """
    test_main()

# original C code is below:
#void moment( data, int n, float *ave, float *adev, float *sdev, float *var, float *skew, float *curt)
#Given an array of data[1..n], this routine returns its mean ave, average deviation adev, standard deviation sdev, variance var, skewness skew, and kurtosis curt.
#{
#void nrerror(char error_text[]); int j;
#float ep=0.0,s,p;
#if (n <= 1) nrerror("n must be at least 2 in moment");
#    j=1
#      s=0.0;
#for (j=1;j<=n;j++) s += data[j]; *ave=s/n; *adev=(*var)=(*skew)=(*curt)=0.0; for (j=1;j<=n;j++) {
#*adev += fabs(s=data[j]-(*ave)); ep += s;
#*var += (p=s*s);
#*skew += (p *= s);
#*curt += (p *= s); }
#*adev /= n; *var=(*var-ep*ep/n)/(n-1); *sdev=sqrt(*var);
#if (*var) {
#*skew /= (n*(*var)*(*sdev));
#First pass to get the mean.
#$Second pass to get the first (absolute), sec- ond, third, and fourth moments of the deviation from the mean.
#$Corrected two-pass formula.
#Put the pieces together according to the con-
#ventional definitions.
#*curt=(*curt)/(n*(*var)*(*var))-3.0;
#} else nrerror("No skew/kurtosis when variance = 0 (in moment)");


