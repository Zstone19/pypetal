import petl.pipeline as pl
import glob

main_dir = '../examples/dat/javelin_'
filenames = [main_dir + 'continuum.dat', main_dir + 'yelm.dat', main_dir + 'zing.dat']

output_dir = '.tmp/'
line_names = ['continuum', 'yelm', 'zing']


params = {
    'reject_data': [True, False, True]
}

res = pl.run_pipeline(filenames, output_dir, line_names,
                      run_drw_rej=True,
                      drw_rej_params=params)


main_directories = glob.glob('.tmp/*/')
subdirs = ['.tmp/continuum/', '.tmp/yelm/', '.tmp/zing/']


#Make sure each line has a subdirectory
if '.tmp/continuum/' not in main_directories:
    raise Exception('continuum directory not found')

if '.tmp/yelm/' not in main_directories:
    raise Exception('yelm directory not found')

if '.tmp/zing/' not in main_directories:
    raise Exception('zing directory not found')




#Make sure each line has a "drw_rej" subdirectory
for fdir in sub_dirs:
    subdirs = glob.glob(fdir + '*')
    if fdir + 'drw_rej/' not in subdirs:
        raise Exception('drw_rej directory not found in ' + fdir)


    

#Make sure each "drw_rej" subdirectory has a chain,drw_fit, and mask file
for i, fdir in enumerate(sub_dirs):
    rej_dir = fdir + 'drw_rej/'
    files = glob.glob(rej_dir + '*')
    
    if rej_dir + line_names[i] + '_chain.dat' not in files:
        raise Exception('chain file not found in ' + rej_dir)
    
    if rej_dir + line_names[i] + '_drw_fit.dat' not in files:
        raise Exception('drw_fit file not found in ' + rej_dir)
    
    if rej_dir + line_names[i] + '_mask.dat' not in files:
        raise Exception('mask file not found in ' + rej_dir)
    
    if rej_dir + line_names[i] + '_drw_fit.pdf' not in files:
        raise Exception('figure not found in ' + rej_dir)
    
    
    
    
#Make sure light curves get saved
if '.tmp/light_curves/' not in main_directories:
    raise Exception('light curve directory not found')

lc_files = glob.glob('.tmp/light_curves/*')
if '.tmp/light_curves/continuum.dat' not in lc_files:
    raise Exception('continuum light curve not found')

if '.tmp/light_curves/yelm.dat' not in lc_files:
    raise Exception('yelm light curve not found')

if '.tmp/light_curves/zing.dat' not in lc_files:
    raise Exception('zing light curve not found')




#Make sure "reject_data" feature works
if '.tmp/continuum_data.dat' not in glob.glob('.tmp/*'):
    raise Exception('continuum_data.dat not found')

if '.tmp/yelm_data.dat' not in glob.glob('.tmp/*'):
    raise Exception('yelm_data.dat not found')
