import numpy as np

from pypetal.pyroa.plotting import (plot_fits, plot_histograms,
                                    pyroa_corner_plot, pyroa_trace_plot)
from pypetal.pyroa.utils import run_pyroa
from pypetal.utils import defaults


def pyroa_tot(cont_fname, line_fnames, line_names, output_dir,
              general_kwargs, kwargs):

   #--------------------------------------------------
    #Read general kwargs

    verbose = general_kwargs['verbose']
    plot = general_kwargs['plot']
    time_unit = general_kwargs['time_unit']
    lc_unit = general_kwargs['lc_unit']
    lag_bounds = general_kwargs['lag_bounds']

    #--------------------------------------------------
    #Read kwargs

    nchain, nburn, init_tau, subtract_mean, div_mean, \
            add_var, delay_dist, psi_types, together, \
                objname = defaults.set_pyroa( kwargs, len(line_names) )

    if verbose:

        txt_str = """
Running PyROA
----------------
nburn: {}
nchain: {}
init_tau: {}
subtract_mean: {}
div_mean: {}
add_var: {}
delay_dist: {}
psi_types: {}
together: {}
objname: {}
----------------
        """.format(nburn, nchain, init_tau, subtract_mean, div_mean, add_var,
                     delay_dist, psi_types, together, objname)

        print(txt_str)

    tot_fnames = np.hstack( [ [cont_fname], line_fnames ] )

    lc_dir = output_dir + 'pyroa_lcs/'
    if not together:
        line_dir = [ output_dir + x + '/pyroa/' for x in line_names[1:]]
    else:
        line_dir = output_dir + 'pyroa/'



    res = run_pyroa( tot_fnames, lc_dir, line_dir, line_names,
                           nburn, nchain, lag_bounds, init_tau,
                           together=together, subtract_mean=subtract_mean,
                           div_mean=div_mean, add_var=add_var,
                           delay_dist=delay_dist, psi_types=psi_types,
                           objname=objname, verbose=verbose)

    lc_fnames = [ lc_dir + objname + '_' + x + '.dat' for x in line_names ]

    if together:
        pyroa_trace_plot( res.samples, line_names, add_var=add_var,
                                delay_dist=delay_dist, nburn=nburn,
                                fname = output_dir + 'pyroa/trace_plot.pdf',
                                show=plot)

        plot_histograms( res.samples, line_names, nburn=nburn,
                               add_var=add_var, delay_dist=delay_dist,
                               fname= output_dir + 'pyroa/histogram_plot.pdf',
                               show=plot)

        pyroa_corner_plot( res.samples, line_names, nburn=nburn,
                                 add_var=add_var, delay_dist=delay_dist,
                                 split=False,
                                 fname = output_dir + 'pyroa/corner_plot.pdf',
                                 show=plot)

        plot_fits( lc_fnames, line_names, res.samples, res.models,
                         nburn=nburn, add_var=add_var, delay_dist=delay_dist,
                         psi_types=psi_types, delimiter=None,
                         time_unit=time_unit, lc_unit=lc_unit,
                         output_fname = output_dir + 'pyroa/fits_plot.pdf',
                         show=plot)


    else:

        for i, res_i in enumerate(res):
            pyroa_trace_plot( res_i.samples, line_names, add_var=add_var[i],
                                    delay_dist=delay_dist[i], nburn=nburn,
                                    fname = output_dir + line_names[i+1] + '/pyroa/trace_plot.pdf',
                                    show=plot)

            plot_histograms( res_i.samples, line_names, nburn=nburn,
                                   add_var=add_var[i], delay_dist=delay_dist[i],
                                   fname= output_dir + line_names[i+1] + '/pyroa/histogram_plot.pdf',
                                   show=plot)

            pyroa_corner_plot( res_i.samples, line_names, nburn=nburn,
                                     add_var=add_var[i], delay_dist=delay_dist[i],
                                     split=False,
                                     fname = output_dir + line_names[i+1] + '/pyroa/corner_plot.pdf',
                                     show=plot)

            plot_fits( [lc_fnames[0], lc_fnames[i+1]], line_names, res_i.samples, res_i.models,
                             nburn=nburn, add_var=add_var[i], delay_dist=delay_dist[i],
                             psi_types=psi_types[i], delimiter=None,
                             time_unit=time_unit, lc_unit=lc_unit,
                             output_fname = output_dir + line_names[i+1] + '/pyroa/fits_plot.pdf',
                             show=plot)


    return res
