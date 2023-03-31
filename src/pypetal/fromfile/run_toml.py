import toml
import numpy as np

import os
import glob


def str2none(x):
    if x == 'None':
        return None
    else:
        return x

#########################################################################################
################################# ASSIST FUNCTIONS ######################################
#########################################################################################

def make_toml(output_dir, fnames, 
              run_arr, param_arr,
              line_names=None, filename=None):

    tot_keys = ['drw_rej', 'detrend', 'pyccf', 'pyzdcf', 'pyroa', 'javelin', 'weighting']
    toml_dict = {}


    #Inputs
    toml_dict['inputs'] = {}
    toml_dict['inputs']['output_dir'] = output_dir
    toml_dict['inputs']['filenames'] = fnames

    if line_names is not None:
        toml_dict['inputs']['line_names'] = line_names





    toml_dict['params'] = {}

    #General parameters
    if param_arr[0] != {}:
        toml_dict['params']['general'] = param_arr[0]


    #Module parameters
    for i in range(len(run_arr)):
        if run_arr[i]:
            toml_dict['params'][tot_keys[i]] = param_arr[i+1]

    #Write to file
    if filename is not None:
        toml.dump(toml_dict, open(filename, 'w+'))

    return toml_dict





def get_toml_modules(filename):

    tot_keys = ['drw_rej', 'detrend', 'pyccf', 'pyzdcf', 'pyroa', 'javelin', 'weighting']
    run_arr = []

    toml_dat = toml.load(filename)['params']

    for key in tot_keys:
        if key in toml_dat.keys():
            if 'run' in toml_dat[key]:
                run_val = toml_dat[key]['run']
            else:
                run_val = False
        else:
            run_val = False


        run_arr.append(run_val)

    return run_arr




def get_toml_params(filename, run_arr):
    
    tot_keys = ['drw_rej', 'detrend', 'pyccf', 'pyzdcf', 'pyroa', 'javelin', 'weighting']
    param_arr = []


    toml_dat = toml.load(filename)['params']

    #General parameters
    if 'general' in toml_dat.keys():
        param_arr.append(toml_dat['general'])
    else:
        param_arr.append({})


    #Module parameters
    for i in range(len(run_arr)):

        if run_arr[i]:
            param_arr.append(toml_dat[tot_keys[i]])
        else:
            param_arr.append({})

    return param_arr








def get_toml_inputs(filename):

    toml_dat = toml.load(filename)

    assert 'inputs' in toml_dat.keys(), 'No inputs found in toml file'
    assert 'output_dir' in toml_dat['inputs'].keys(), 'No output directory found in toml file'        
    assert 'filenames' in toml_dat['inputs'].keys(), 'No filenames found in toml file'


    #################################
    #Output directory
    
    output_dir = os.path.abspath(toml_dat['inputs']['output_dir']) + r'/'


    #################################
    #Filenames

    #If input is a list
    if isinstance( toml_dat['inputs']['filenames'], list ):
        for i in range(len(toml_dat['inputs']['filenames'])):
            assert os.path.isfile(toml_dat['inputs']['filenames'][i]), 'File not found: ' + toml_dat['inputs']['filenames'][i]

        filenames = toml_dat['inputs']['filenames']



    #If input is a string
    if isinstance( toml_dat['inputs']['filenames'], str ):
        exists = os.path.exists(toml_dat['inputs']['filenames'])
        isdir = os.path.isdir(toml_dat['inputs']['filenames'])


        #Assume all files in directory are inputs
        if exists and isdir:
            filenames = glob.glob(toml_dat['inputs']['filenames']+'*')

        #Assume input is a glob pattern
        elif exists and (not isdir):
            filenames = glob.glob(toml_dat['inputs']['filenames'])
            assert len(filenames) > 0, 'No files found matching glob pattern: ' + toml_dat['inputs']['filenames']


        for i in range(len(filenames)):
            assert os.path.isfile(filenames[i]), 'File not found: ' + filenames[i]




    #################################
    #Line names

    if str2none(toml_dat['inputs']['line_names']) is None:
        line_names = []
        for i in range(len(filenames)):
            line_names.append('Line {}'.format(i+1))

    else:
        assert isinstance( toml_dat['inputs']['line_names'], list ), 'Line names must be a list'
        line_names = toml_dat['inputs']['line_names']

    return output_dir, filenames, line_names



#########################################################################################
################################ TO RUN SINGLE OBJECTS ##################################
#########################################################################################

#Before javelin
def run_from_toml1(filename):
    import pypetal.pipeline as pl

    output_dir, filenames, line_names = get_toml_inputs(filename)
    run_arr = get_toml_modules(filename)
    param_arr = get_toml_params(filename, run_arr)

    res = pl.run_pipeline(output_dir, filenames, line_names,
                          run_drw_rej=run_arr[0], drw_rej_params=param_arr[1],
                          run_detrend=run_arr[1], detrend_params=param_arr[2],
                          run_pyccf=run_arr[2], pyccf_params=param_arr[3],
                          run_pyzdcf=run_arr[3], pyzdcf_params=param_arr[4],
                          run_pyroa=run_arr[4], pyroa_params=param_arr[5],
                          **param_arr[0])

    return res


#After javelin
def run_from_toml2(filename):
    import pypetal.pipeline as pl

    output_dir, _, line_names = get_toml_inputs(filename)
    run_arr = get_toml_modules(filename)
    param_arr = get_toml_params(filename, run_arr)

    if run_arr[-1]:
        res = pl.run_weighting(output_dir, line_names,
                               run_pyccf=run_arr[2], pyccf_params=param_arr[3],
                               run_pyroa=run_arr[4], pyroa_params=param_arr[5],
                               run_javelin=run_arr[5], javelin_params=param_arr[6],
                               weighting_params=param_arr[7],
                               **param_arr[0])

        return res

    else:
        return
    
#########################################################################################
################################ TO RUN MULTIPLE OBJECTS ################################
#########################################################################################

def run_all1(filenames):
    for f in filenames:
        _ = run_from_toml1(f)
        
    return

def run_all3(filenames):
    for f in filenames:
        _ = run_from_toml2(f)
        
    return
