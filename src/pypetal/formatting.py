import numpy as np
import astropy.units as u
import re



def str2unit(x):

    if (x == 'Arbitrary Units'):
        unit = u.dimensionless_unscaled
    else:
        unit = u.Unit(x)

    return unit


#From https://stackoverflow.com/questions/5807952/removing-trailing-zeros-in-python/5808661#5808661
def number_shaver(ch,
                  regx = re.compile('(?<![\d.])0*(?:'
                                    '(\d+)\.?|\.(0)'
                                    '|(\.\d+?)|(\d+\.\d+?)'
                                    ')0*(?![\d.])')  ,
                  repl = lambda mat: mat.group(mat.lastindex)
                                     if mat.lastindex!=3
                                     else '0' + mat.group(3) ):
    return regx.sub(repl,ch)

def err2str(val, up_err, lo_err, dec=2):
    val = number_shaver( str(round(val, dec)) )
    up_err = number_shaver( str(round(up_err, dec)) )
    lo_err = number_shaver( str(round(lo_err, dec)) )
    
    return val + '^{+' + up_err + '}_{-' + lo_err + '}'


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
            
    