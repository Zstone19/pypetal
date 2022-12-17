import numpy as np
import matplotlib.pyplot as plt
from linmix import LinMix

from pypetal.petalio import write_data


import matplotlib as mpl

mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.direction'] = 'in'

mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = 'in'

mpl.rcParams["figure.autolayout"] = False

mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.format'] = 'pdf'



def detrend(x, y, yerr, K=2, parallelize=False, 
            time_unit='d', lc_unit='',
            output_dir=None, verbose=False, plot=False):
    
    
    if lc_unit == '':
        flux_txt = 'Flux'
    elif lc_unit == 'mag':
        flux_txt = 'Magnitude'
    else:
        flux_txt = 'Flux [{}]'.format(lc_unit)
    
    
    lm = LinMix(x, y, ysig=yerr, K=K, parallelize=parallelize)
    lm.run_mcmc(silent=True)
    
    b_chains = []
    m_chains = []
    sigsqr_chains = []
    
    for i in range(len(lm.chain)):
        b_chains.append( lm.chain[i][0] )
        m_chains.append( lm.chain[i][1] )
        sigsqr_chains.append( lm.chain[i][1] )

        
    b_med = np.median(b_chains)
    b_err_hi = np.percentile(b_chains, 84) - b_med
    b_err_lo = b_med - np.percentile(b_chains, 16)

    m_med = np.median(m_chains)
    m_err_hi = np.percentile(m_chains, 84) - m_med
    m_err_lo = b_med - np.percentile(m_chains, 16)

    sigsqr_med = np.median(sigsqr_chains)
    sigsqr_err_hi = np.percentile(sigsqr_chains, 84) - sigsqr_med
    sigsqr_err_lo = sigsqr_med - np.percentile(sigsqr_chains, 16)


    if verbose:
        print('m = {:.3f} + {:.3f} - {:.3f}'.format(m_med, m_err_hi, m_err_lo) )
        print('b = {:.3f} + {:.3f} - {:.3f}'.format(b_med, b_err_hi, b_err_lo) )
        print('sigsqr = {:.3f} + {:.3f} - {:.3f}'.format(sigsqr_med, sigsqr_err_hi, sigsqr_err_lo) )   
        
        
        
    y0_vals = []
    for i in range(len(lm.chain)):
        m_val = lm.chain[i][1]
        b_val = lm.chain[i][0]
        x0 = np.min(x)
        
        y0_vals.append( m_val*x0 + b_val )


    lo_ind = np.argmin(y0_vals)
    m_lo = lm.chain[lo_ind][1]
    b_lo = lm.chain[lo_ind][0]

    hi_ind = np.argmax(y0_vals)
    m_hi = lm.chain[hi_ind][1]
    b_hi = lm.chain[hi_ind][0]
        

        
    xsim = np.linspace(x.min(), x.max(), 100)        
    ysim_avg = m_med*xsim + b_med
    ysim_lo = m_lo*xsim + b_lo
    ysim_hi = m_hi*xsim + b_hi

        
        
    fig, ax = plt.subplots()
        
    _, _, bars = ax.errorbar(x, y, yerr, fmt='.k')
    [bar.set_alpha(.3) for bar in bars]

        
    ax.plot(xsim, ysim_avg, 'r')
    ax.fill_between(xsim, ysim_lo, ysim_hi, color='r', alpha=.3)

    if lc_unit == 'mag':
        ax.invert_yaxis()
    
    ax.set_ylabel(flux_txt)
    ax.set_xlabel('Time [{}]'.format(time_unit))
    
    ax.tick_params('both', which='major', length=8)
    ax.tick_params('both', which='minor', length=3)
    
    if output_dir is not None:    
        plt.savefig(output_dir + 'detrend.pdf', dpi=200, bbox_inches='tight')
    
    if plot:    
        plt.show()
        
    plt.cla()
    plt.clf()
    plt.close()
    
    
    y_dt = y - (m_med*x + b_med)
    
    return y_dt, [m_err_lo, m_med, m_err_hi], \
        [b_err_lo, b_med, b_err_hi], \
        [sigsqr_err_lo, sigsqr_med, sigsqr_err_hi]     









def detrend_tot(output_dir, cont_fname, line_fnames, line_names, general_kwargs):

    verbose = general_kwargs['verbose']
    plot = general_kwargs['plot']
    threads = general_kwargs['threads']
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']
    
    fnames_tot = np.hstack([ [cont_fname], line_fnames])
    
    if threads > 1:
        parallelize = True
    else:
        parallelize = False    
        
        
    txt = """
    Running detrending
    -------------------
    parallelize: {}
    K = 2
    """.format(parallelize)
    
    if verbose:
        print(txt)
            
    
    m_vals = []
    b_vals = []
    sigsqr_vals = []
    
    for i, fname in enumerate(fnames_tot):        
        x, y, yerr = np.loadtxt( fname, unpack=True, usecols=[0,1,2], delimiter=',' )
        
        y_dt, m_res, b_res, sigsqr_res = detrend(x, y, yerr, K=2, parallelize=parallelize, 
                       time_unit=time_unit, lc_unit=lc_unit[i],
                       output_dir=output_dir + line_names[i] + r'/',
                       verbose=verbose, plot=plot)
                
        write_data( output_dir + 'processed_lcs/' + line_names[i] + '_data.dat' )
        
        m_vals.append(m_res)
        b_vals.append(b_res)
        sigsqr_vals.append(sigsqr_res)
    
        
    output = {
        'm': m_vals,
        'b': b_vals,
        'sigsqr': sigsqr_vals
    }
        
    return output
        
    
    
    