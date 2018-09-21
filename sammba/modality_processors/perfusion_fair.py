# -*- coding: utf-8 -*-
"""
Created on Sunday 4th March 2018

@author: Nachiket Nadkarni
"""


import pandas as pd
import nibabel as nib
import numpy as np
# from scipy.optimize import least_squares as ls
import lmfit
from multiprocessing import cpu_count, Pool
from functools import partial
import tqdm


def _perf_fair_read_ptbl(nii_fname):
    """
    Extract perfusion FAIR acquisition parameters (TIs, order of selective and
    non-selective inversions) from the *_ptbl.txt file of a *.nii.gz file. Note
    that both are created simultaneously when sammba.io_conversions.dcm_to_nii
    is used for DICOM to NIfTI-1 conversion. If a different converter was used
    that does not export this data, or does it differently, this function will
    not work.

    Parameters
    ----------
    nii_fname : str
        Path to the perfusion *.nii.gz file converted by
        sammba.io_conversions.dcm_to_nii, NOT to the *_ptbl.txt file itself.
        The *_ptbl.txt file name is determined automatically from that of the
        *.nii.gz.

    Returns
    -------
    Dictionary of parameters needed for perfusion fitting:
    ti_list:
        simple list of all TIs
    long_ti_list:
        full list of TIs, the same length and in the same order as the
        acquisitions (the final dimension of the NIfTI-1 file)
    fc_list:
        full list of inversion type descriptions ("Selective" or
        "Non-selective"), the same length and in the same order as the
        acquisitions (the final dimension of the NIfTI-1 file). fc is a DICOM
        acronym for frame comment, the DICOM field used for storing the
        inversion type descriptions
    picker_sel:
        list of positions of selective inversions in the acquisition
    picker_nonsel:
        list of positions of non-selective inversions in the acquisition
        
    Notes
    -----    
    Assuming the acqusition is as currently specified, for a given TI, there is
    always one selective inversion and one non-selective inversion. The
    perf_fair_fitter function assumes this. If this is not the case then it
    will fail, and so will the whole perfusion processing procedure. So the
    length of ti must be the same as picker_sel and picker_nonsel. The long_ti
    and fc outputs are diagnostic in case something goes wrong.
    """
    
    # the *_ptbl.txt file will have the same name as the perfusion file itself,
    # but with _ptbl.txt at the end rather than just .nii.gz. arguably this
    # function should directly take the text file name as input, but it is
    # difficult to envisage a scenario where this would be needed. the
    # *_ptbl.txt file has a specific format created by
    # sammba.io_conversions.dcm_to_nii. if someone supplies their own file
    # name, it would mean they changed the *_ptbl.txt name (why?) or they have
    # their own formatted text file which this function would not be able to
    # process
    
    bfx = nii_fname.split('_')[-1].split('.')[0] # bfx = bruker folder number
    ptbl = pd.read_table(nii_fname.replace(bfx + '.nii.gz', bfx + '_ptbl.txt'))
    ti_list = ptbl.TI[ptbl.slice == 1][ptbl.FC == 'Selective Inversion'].tolist()
    long_ti_list = ptbl.TI[ptbl.slice == 1].tolist()
    fc_list = ptbl.FC[ptbl.slice == 1].tolist()
    picker_sel = [n for n,x in enumerate(fc_list) if x == 'Selective Inversion']
    picker_nonsel = [n for n,x in enumerate(fc_list) if x == 'Non-selective Inversion']
    return {'TI':ti_list, 'long_TI':long_ti_list, 'FC':fc_list,
            'picker_sel':picker_sel, 'picker_nonsel':picker_nonsel}


# the perfusion fluid-attenuated inversion-recovery function
# the jacobian should be callable but I do not know enough maths to create it
# pars = parameters
def _fair_t1_func(pars, s0, ti):
    
    # return pars[0] + np.absolute(pars[1] * (1 - 2 * np.exp(-ti / pars[2]))) - s0
    bias, m0, t1 = pars['bias'], pars['M0'], pars['T1']
    m = bias + np.absolute(m0 * (1.0 - 2.0 * np.exp(-ti / t1)))
    return m - s0

def _fair_t1_fit(s0, ti, t1_guess):
    """
    Fit the perfusion FAIR equation:
        
        signal = bias + abs(M0 * (1 - 2 * (exp(-TI/T1))))
        
    using scipy.optimize.least_squares (default Trust Region Reflective
    algorithm). Bias, M0 and T1 are all estimated with a zero lower bound.
    Starting values for bias (zero) and M0 (mean of input signals) are
    calculated automatically; T1 has to be supplied, with the mean of TIs
    being a good guess.
    
    Parameters
    ----------
    s0 : numpy array of int or float
        The aquired signals.
        
    ti : list of int or float
        The inversion times. Must have the same length as s0.
        
    t1_guess : int or float
        An initial starting value for T1. The mean of TIs is a good guess.

    Returns
    -------
    class scipy.optimize.OptimizeResult    
    
    Notes
    ----- 
    Levenberg-Marquardt (method='lm'), the usual literature algorithm, cannot
    be used here as the implementation in scipy does not accept bounds. From
    experience with the R MINPACK implementation (minpack.lm function nlsLM),
    these bounds are necessary to assure good fitting. The lmfit python package
    does implement bounds for MINPACK Levenberg-Marquardt (method=leastsq;
    default), though the results appear to be worse than for TRR
    (method=least_squares), whose results bizarrely look closer to those of R
    minpack.lm nlsLM (to confirm), and are of course identical to
    scipy.optimize.least_squares with method='trf' (default).
    """

    # return ls(_fair_t1_func, np.array([0, np.mean(s0), t1_guess]),
    #           bounds=([0, 0, 0], np.inf), args=(s0, np.array(ti)))
    p = lmfit.Parameters()
    p.add('bias', value=0, min=0.0)
    p.add('M0', value=np.mean(s0), min=0.0)
    p.add('T1', value=t1_guess, min=0.0)
    minner = lmfit.Minimizer(_fair_t1_func, params=p,
                             fcn_args=(s0, np.asarray(ti)),
                             nan_policy='propagate')
    return minner.minimize()


def _perf_fair_fit(s0, t1_blood, ti, t1_guess, picker_sel, picker_nonsel,
                   lambda_blood=0.9, multiplier=6000000, outtype='simple'):
    """
    Wrapper to execute fair_t1_fit on a real signal vector containing both
    selective and non-selective inversions. Also calculates rCBF and absolute
    CBF from the resulting selective and non-selective T1 values.
    
    Parameters
    ----------
    s0 : numpy array of int or float
        The aquired signals.
        
    t1_blood : int or float
        T1 of blood in ms at the acquisition field strength. Empirical
        constant; an example value is 2800 at 11.7T
        (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3752414/table/T1,
        Table 1 of Blood longitudinal (T1) and transverse (T2) relaxation time
        constants at 11.7 Tesla, MAGMA. 2012 Jun; 25(3): 245–249,
        Ai-Ling Lin, Qin Qin, Xia Zhao, and Timothy Q. Duong).
        
    ti : list of int or float
        The inversion times. Must have the same length as s0.
        
    t1_guess : int or float
        An initial starting value for T1. The mean of TIs is a good guess.
    
    picker_sel : list or numpy array of int
        Vector indicating positions of selectively inverted signals in s0. The
        lengths of picker_sel and picker_sel must sum to the length of s0.
        
    picker_nonsel : list or numpy array of int
        Vector indicating positions of non-selectively inverted signals in s0.
        The lengths of picker_sel and picker_sel must sum to the length of s0.

    lambda_blood : float, optional
        The assumed blood tissue partition coefficient of water in ml per g.
        Empirical constant; usually 0.9.
        
    multiplier : int or float, optional
        The absolute CBF result is initially produced in units of ml per g per
        ms. The multiplier converts to desired units. Usually it is 6000000 to
        convert to ml per 100g per min.
        
    outtype : {'simple', 'complicated'}, optional
        Changes return type. If 'simple', return a list of parameter, rCBF and
        absolute CBF values. If 'complicated', return the two classes of
        scipy.optimize.OptimizeResult plus rCBF and CBF values.
    
    Returns
    -------
    If outtype='simple':
        List of 14 floats- selective inversion bias, M0 and T1, non-selective
        inversion bias, M0 and T1, all paired with their standard errors,
        plus rCBF and CBF. Failed fits produce zeroes.
    If outtype='complicated':
        List of 4- two scipy.optimize.OptimizeResult classes (selective fit
        then non-selective) plus floats for rCBF then CBF.
    """
    
    # there are no standard errors as these are not in the current output of
    # scipy.optimize.least_squares. the lmfit library can produce them, but
    # during testing they were usually incalculable, and lmfit was twice as
    # slow (or maybe I was not using it properly?) it is unlikely that the
    # errors are of much use any way so will stay with
    # scipy.optimize.least_squares.
    
    if outtype not in ['simple', 'complicated']:
        raise ValueError('unrecognized outtype')
    
    s0_sel = np.array(s0)[np.array(picker_sel)]
    s0_nonsel = np.array(s0)[np.array(picker_nonsel)]

    r_sel = _fair_t1_fit(s0_sel, ti, t1_guess)
    r_nonsel = _fair_t1_fit(s0_nonsel, ti, t1_guess)
    
    if r_sel.success and r_nonsel.success:
        
        # t1_sel = r_sel.x[2]
        # t1_nonsel = r_nonsel.x[2]
        t1_sel = r_sel.params['T1'].value
        t1_nonsel = r_nonsel.params['T1'].value
        
        rCBF = 100 * (t1_nonsel - t1_sel) / t1_nonsel
        CBF = multiplier * lambda_blood * (
              (t1_nonsel / t1_blood) * ((1 / t1_sel) - (1 / t1_nonsel)))
     
        if outtype == 'simple':
            # r = [x for x in r_sel.x] + [x for x in r_nonsel.x] + [rCBF, CBF]
            r = [r_sel.params['bias'].value, r_sel.params['bias'].stderr,
                 r_sel.params['M0'].value, r_sel.params['M0'].stderr,
                 r_sel.params['T1'].value, r_sel.params['T1'].stderr,
                 r_nonsel.params['bias'].value, r_nonsel.params['bias'].stderr,
                 r_nonsel.params['M0'].value, r_nonsel.params['M0'].stderr,
                 r_nonsel.params['T1'].value, r_nonsel.params['T1'].stderr,
                 rCBF, CBF]
            return [x if x is not None else 0 for x in r]
        else:
            return [r_sel, r_nonsel, rCBF, CBF]
    
    else:
        
        if outtype == 'simple':
            # return np.repeat(0, 8)
            return np.repeat(0, 14)
        else:
            return ['r_sel fit ' + r_sel.success,
                    'r_nonsel fit ' + r_nonsel.success,
                    'rCBF nan', 'CBF nan']


def _perf_fair_fit_mp(all_s0, t1_blood, ti, t1_guess, picker_sel,
                      picker_nonsel, ncpu=cpu_count() - 1, **kwargs):
    """
    Wrapper to execute perf_fair_fitter (outtype='simple') in parallel on
    multiple signal vectors, tracking execution with a progress bar.

    Parameters
    ----------
    all_s0 : 2D numpy array of int or float
        The aquired signals. The first dimension is the number of s0 vectors,
        the second dimension the length of all of them. So for a real image,
        usually len(D1) >> len(D2).

    t1_blood : int or float
        T1 of blood in ms at the acquisition field strength. Empirical
        constant; an example value is 2800 at 11.7T
        (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3752414/table/T1,
        Table 1 of Blood longitudinal (T1) and transverse (T2) relaxation time
        constants at 11.7 Tesla, MAGMA. 2012 Jun; 25(3): 245–249,
        Ai-Ling Lin, Qin Qin, Xia Zhao, and Timothy Q. Duong).

    ti : list of int or float
        The inversion times. Must have the same length as s0.

    t1_guess : int or float
        An initial starting value for T1. The mean of TIs is a good guess.

    picker_sel : list or numpy array of int
        Vector indicating positions of selectively inverted signals in s0. The
        lengths of picker_sel and picker_sel must sum to the length of s0.

    picker_nonsel : list or numpy array of int
        Vector indicating positions of non-selectively inverted signals in s0.
        The lengths of picker_sel and picker_sel must sum to the length of s0.

    ncpu : int, optional
        Number of processes to launch in parallel. Defaults to using all but
        one of the available CPUs.

    kwargs : dict, optional
        Additional keyword arguments passed to perf_fair_fitter.

    Returns
    -------
    List of all_s0.shape[0] lists, each of 14 floats (see _perf_fair_fit).
    """
    
    # inspiration provided by:
    # https://stackoverflow.com/q/5442910
    # https://stackoverflow.com/a/45276885
    
    # cannot use a with statement in python 2.7; may cause problems
    pool = Pool(processes=ncpu)
    return list(tqdm.tqdm(pool.imap(partial(_perf_fair_fit, t1_blood=t1_blood,
                                            ti=ti, t1_guess=t1_guess,
                                            picker_sel=picker_sel,
                                            picker_nonsel=picker_nonsel,
                                            **kwargs),
                                    all_s0), total = len(all_s0)))


def perf_fair_nii_proc(nii_in_fname, t1_blood, ti, t1_guess, picker_sel,
                       picker_nonsel, nii_out_fname=None, **kwargs):
    """
    Wrapper to execute perf_fair_fitter_mp on a NIfTI-1 image.

    Parameters
    ----------
    nii_in_fname : str
        Input file path.

    t1_blood : int or float
        T1 of blood in ms at the acquisition field strength. Empirical
        constant; an example value is 2800 at 11.7T
        (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3752414/table/T1,
        Table 1 of Blood longitudinal (T1) and transverse (T2) relaxation time
        constants at 11.7 Tesla, MAGMA. 2012 Jun; 25(3): 245–249,
        Ai-Ling Lin, Qin Qin, Xia Zhao, and Timothy Q. Duong).

    ti : list of int or float
        The inversion times. Must have the same length as s0.

    t1_guess : int or float
        An initial starting value for T1. The mean of TIs is a good guess.

    picker_sel : list or numpy array of int
        Vector indicating positions of selectively inverted signals in s0. The
        lengths of picker_sel and picker_sel must sum to the length of s0.

    picker_nonsel : list or numpy array of int
        Vector indicating positions of non-selectively inverted signals in s0.
        The lengths of picker_sel and picker_sel must sum to the length of s0.

    nii_out_fname : str
        Output file path. If none, will be the same as the input file path,
        but suffixed with _proc.

    kwargs : dict, optional
        Additional keyword arguments passed to _perf_fair_fit_mp.

    Returns
    -------
    NIfTI-1 file saved to nii_out_fname. There are eight images in the time
    dimension, for:
        selective inversion bias, M0 and T1, non-selective inversion bias, M0
        and T1, rCBF (relative cerebral blood flow) and CBF (absolute).
    Failed fits produce zero-valued voxels.
    """

    nii_in = nib.load(nii_in_fname)
    in_mat = nii_in.get_data()
    all_s0 = in_mat.reshape((np.product(in_mat.shape[:-1]), in_mat.shape[-1]))
    r = _perf_fair_fit_mp(all_s0, t1_blood, ti, t1_guess, picker_sel,
                          picker_nonsel, **kwargs)
    r = np.array(r)
    img = nib.Nifti1Image(np.reshape(r, in_mat.shape[:-1] + (r.shape[1],)),
                          nii_in.get_affine())
    
    if nii_out_fname is None:
        nii_out_fname = nii_in_fname.replace('.nii.gz', '_proc.nii.gz')
    return img.to_filename(nii_out_fname)


def perf_fair_niiptbl_proc(nii_in_fname, t1_blood, **kwargs):
    """
    Wrapper to execute perf_fair_nii_proc supplied with acquisition parameters
    automatically extracted using perf_fair_read_ptbl. So TI and inversion type
    are automatically determined for all signals, as is a T1 starting guess.
    
    Parameters
    ----------
    nii_in_fname : str
        Input file path.
        
    t1_blood : int or float
        T1 of blood in ms at the acquisition field strength. Empirical
        constant; an example value is 2800 at 11.7T
        (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3752414/table/T1,
        Table 1 of Blood longitudinal (T1) and transverse (T2) relaxation time
        constants at 11.7 Tesla, MAGMA. 2012 Jun; 25(3): 245–249,
        Ai-Ling Lin, Qin Qin, Xia Zhao, and Timothy Q. Duong).

    kwargs : dict, optional
        Additional keyword arguments passed to perf_fair_nii_proc.
    
    Returns
    -------
    NIfTI-1 file saved to nii_out_fname. There are eight images in the time
    dimension, for:
        selective inversion bias, M0 and T1, non-selective inversion bias, M0
        and T1, rCBF (relative cerebral blood flow) and CBF (absolute).
    Failed fits produce zero-valued voxels.
    """

    ptbl_dict = _perf_fair_read_ptbl(nii_in_fname)
    ti = ptbl_dict['TI']
    t1_guess = np.mean(ptbl_dict['TI'])
    picker_sel = ptbl_dict['picker_sel']
    picker_nonsel = ptbl_dict['picker_nonsel']
    perf_fair_nii_proc(nii_in_fname, t1_blood, ti, t1_guess, picker_sel,
                       picker_nonsel, **kwargs)
