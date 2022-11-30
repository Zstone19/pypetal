import pypetal.pipeline as pl
import glob
import shutil


main_dir = '../examples/dat/pyccf_'
filenames = [main_dir + 'lc1.dat', main_dir + 'lc2.dat']

output_dir = '.tmp/'
line_names = ['continuum', 'line1']

res = pl.run_pipeline(filenames, output_dir, line_names,
                      run_pyccf=True)


main_directories = glob.glob('.tmp/*/')

line_dirs = glob.glob( '.tmp/line1/*' )
#Make sure pyccf directory exists for the line
if './tmp/line1/pyccf/' not in line_dirs:
    raise Exception('pyccf directory not found')




files = glob.glob( './tmp/line1/pyccf/*' )

#Make sure auxiliary files exist
if '.tmp/line1/pyccf/line1_ccf_dists.dat' not in files:
    raise Exception('ccf_dists file not found')

if '.tmp/line1/pyccf/line1_ccf.dat' not in files:
    raise Exception('ccf file not found')

if '.tmp/line1/pyccf/line1_ccf.pdf' not in files:
    raise Exception('ccf figure not found')




shutil.rmtree('.tmp')