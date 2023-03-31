import os
import numpy as np

from run_toml import make_toml
from pypetal.utils.petalio import write_data


def write_bash(toml_fname, env1, env2):
      
    """Need to use virtualenv. This won't work otherwise."""
        
    txt = """source {}/bin/activate
python run_toml.py "pre" {}
deactivate

source {}/bin/activate
python run_toml.py "jav" {}
deactivate

source {}/bin/activate
python run_toml.py "post" {}
deactivate
""".format( env1, toml_fname, env2, toml_fname, env1, toml_fname ) 
        
        
    
    with open('run.sh', 'w+') as f:
        f.write(txt)    
    
    return



def bash_from_python(output_dir, arg2, 
                    toml_fname, env1, env2,
                    line_names=None,
                    run_drw_rej=False, drw_rej_params={},
                    run_detrend=False, detrend_params={},
                    run_pyccf=False, pyccf_params={},
                    run_pyzdcf=False, pyzdcf_params={},
                    run_pyroa=False, pyroa_params={},
                    run_javelin=False, javelin_params={},
                    run_weighting=False, weighting_params={},
                    **kwargs):

    if arg2 is None:
        raise Exception('Please provide a list of light curve filenames or the light curves themselves')


    if not isinstance(arg2[0], str):
        os.makedirs( output_dir + 'input_lcs/', exist_ok=True )
        fnames = []

        for i in range( len(arg2) ):

            if i == 0:
                name = 'continuum'
            else:
                name = 'line{}'.format(i+1)

            write_data( arg2[i], output_dir + 'input_lcs/' + name + '.dat' )
            fnames.append( output_dir + 'input_lcs/' + name + '.dat' )

        fnames = np.array(fnames)
        kwargs['file_fmt'] = 'csv'
    else:
        fnames = arg2



    run_arr = [run_drw_rej, run_detrend, run_pyccf, run_pyzdcf, run_pyroa, run_javelin, run_weighting]
    param_arr = [kwargs, drw_rej_params, detrend_params, pyccf_params, pyzdcf_params, pyroa_params, javelin_params, weighting_params]

    make_toml(output_dir, fnames, run_arr, param_arr, line_names, filename=toml_fname)
    write_bash(toml_fname, env1, env2)

    return
