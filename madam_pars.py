#!/usr/bin/env python3

def madam_pars(nside=1024, fsample=1000):
    pars = {}

    pars['info'] = 2
    #pars['nthreads'] = 1
    pars['nsubchunk'] = 0
    #pars['isubchunk'] = 1
    #pars['time_unit'] =
    pars['base_first'] = 1
    pars['nshort'] = 10
    pars['nside_map'] = nside
    pars['nside_cross'] = nside // 2
    pars['nside_submap'] = nside // 32
    #pars['good_baseline_fraction'] =
    #pars['concatenate_messages'] =
    pars['allreduce'] = True
    #pars['reassign_submaps'] =
    #pars['pix_mode_map'] =
    #pars['pix_mode_cross'] =
    #pars['pixlim_map'] =
    #pars['pixlim_cross'] =
    #pars['incomplete_matrices'] =
    #pars['allow_decoupling'] =
    #pars['kfirst'] = False
    pars['basis_func'] = 'polynomial'
    pars['basis_order'] = 0
    #pars['iter_min'] =
    pars['iter_max'] = 10
    #pars['cglimit'] =
    pars['fsample'] = fsample
    #pars['mode_detweight'] =
    #pars['flag_by_horn'] =
    #pars['write_cut'] =
    #pars['checknan'] =
    #pars['sync_output'] =
    #pars['skip_existing'] =
    pars['temperature_only'] = False
    #pars['force_pol'] = False
    pars['noise_weights_from_psd'] = False
    #pars['radiometers'] =
    #pars['psdlen'] =
    #pars['psd_down'] =
    #pars['kfilter'] = False
    #pars['diagfilter'] = 0.0
    #pars['precond_width_min'] =
    #pars['precond_width_max'] =
    #pars['use_fprecond'] =
    #pars['use_cgprecond'] =
    #pars['rm_monopole'] = True
    #pars['path_output'] = '/home/klee_ext/kmlee/hpc_data/madam_test/'
    pars['path_output'] = './' 
    pars['file_root'] = 'madam_test'

    pars['write_map'] = True
    pars['write_binmap'] = True
    pars['write_hits'] = True
    pars['write_matrix'] = False#True
    pars['write_wcov'] = False#True
    pars['write_base'] = True
    pars['write_mask'] = True
    pars['write_leakmatrix'] = False
    #pars['write_tod'] = True

    #pars['unit_tod'] =
    #pars['file_gap_out'] =
    #pars['file_mc'] =
    #pars['file_inmask'] =
    #pars['file_spectrum'] =
    #pars['file_gap'] =
    #pars['binary_output'] =
    #pars['nwrite_binary'] =
    #pars['file_covmat'] =
    #pars['detset'] =
    #pars['detset_nopol'] =
    #pars['survey'] = ['hm1:{} - {}'.format(0, nsample / 2),]
    #pars['bin_subsets'] = True
    #pars['mcmode'] =

    return pars


if __name__=="__main__":
    par = madam_pars()
    print (par)
