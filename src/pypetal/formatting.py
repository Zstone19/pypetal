import numpy as np
import astropy.units as u



def str2unit(x):

    if (x == 'Arbitrary Units'):
        unit = u.dimensionless_unscaled
    else:
        unit = u.Unit(x)

    return unit

def err2str(val, up_err, lo_err, dec=2):
    return "%(val)^{+%(up_err)}_{-%(lo_err)}" % {'val': round(val, dec), 'up_err': round(up_err, dec), 'lo_err': round(lo_err, dec)}



def write_data(arr, fname, header=None):

    arr = np.array(arr, dtype=object)

    ndim = len( arr.shape )
    assert ndim <= 2


    if ndim == 2:
        cols = len(arr)
        rows = len(arr[0])    

        with open(fname, "w") as file:

            if header is not None:
                file.write(header + "\n")

            for i in range(rows):

                string = "{},".format(arr[0][i])
                for j in range(1, cols-1):
                    string += "{},".format(arr[j][i])

                string += "{}\n".format(arr[-1][i])

                file.write(string)

    elif ndim == 1:
        rows = len(arr)

        with open(fname, "w") as file:

            if header is not None:
                file.write(header + "\n")

            for i in range(rows):
                file.write( "{}\n".format(arr[i]) )        

    return
            
    