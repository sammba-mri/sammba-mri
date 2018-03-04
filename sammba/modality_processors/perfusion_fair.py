# -*- coding: utf-8 -*-
"""
Created on Sunday 4th March 2018

@author: Nachiket Nadkarni
"""


import pandas as pd
import nibabel as nib
import numpy as np
from scipy.optimize import least_squares as ls
from multiprocessing import Pool
from functools import partial
import tqdm


def perf_fair_read_ptbl(perf_fair_fname):
    """
    Extract important information about the perfusion FAIR acquisition
    parameters from the *_ptbl.txt file created by
    sammba.io_conversions.dcm_to_nii.

    Parameters
    ----------
    filename : str
        Path to the perfusion *.nii.gz file converted by
        sammba.io_conversions.dcm_to_nii, NOT to the *_ptbl.txt file itself.
        The *_ptbl.txt file name is determined automatically fro; that of the
        *.nii.gz.

    Returns
    -------
    Dictionary of parameters useful to perfusion fitting functions.    
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
    
    bfx = perf_fair_fname.split('_')[-1].split('.')[0]
    ptbl = pd.read_table(perf_fair_fname.replace(bfx + '.nii.gz', bfx + '_ptbl.txt'))
    ti = ptbl.TI[ptbl.slice == 1][ptbl.FC == 'Selective Inversion'].tolist()
    long_ti = ptbl.TI[ptbl.slice == 1].tolist()
    fc = ptbl.FC[ptbl.slice == 1].tolist()
    picker_sel = [n for n,x in enumerate(fc) if x == 'Selective Inversion']
    picker_nonsel = [n for n,x in enumerate(fc) if x == 'Non-selective Inversion']
    return {'TI':ti, 'long_TI':long_ti, 'FC':fc,
            'picker_sel':picker_sel, 'picker_nonsel':picker_nonsel}


#the perfusion fluid-attenuated inversion-recovery function
#the jacobian should be callable but I do not know enough maths to create it
def _fair_t1_func(x, ti, s0):
    return x[0] + np.absolute(x[1] * (1.0 - 2.0 * np.exp(-ti / x[2]))) - s0


def fair_t1_fit(s0, ti, t1_guess):
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
    
    return ls(_fair_t1_func, np.array([0, np.mean(s0), t1_guess]),
              bounds=([0, 0, 0], np.inf), args=(np.array(ti), s0))


def perf_fair_fitter(s0, t1_blood, lambda_blood, ti, multiplier, t1_guess,
                     picker_sel, picker_nonsel, outtype='simple'):
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
        constant; an example value is 2800 at 11.7T (Duong ref).
        
    lambda_blood : int or float
        The assumed blood tissue partition coefficient of water in ml per g.
        Empirical constant; usually 0.9.
        
    ti : list of int or float
        The inversion times. Must have the same length as s0.
        
    multiplier : int or float
        The absolute CBF result is initially produced in units of ml per g per
        ms. The multiplier converts to desired units. Usually it is 6000000 to
        convert to ml per 100g per min.    
        
    t1_guess : int or float
        An initial starting value for T1. The mean of TIs is a good guess.
    
    picker_sel : list or numpy array of int
        Vector indicating positions of selectively inverted signals in s0. The
        lengths of picker_sel and picker_sel must sum to the length of s0.
        
    picker_nonsel : list or numpy array of int
        Vector indicating positions of non-selectively inverted signals in s0.
        The lengths of picker_sel and picker_sel must sum to the length of s0.
        
    outtype : {'simple', 'complicated'}, optional
        Changes return type. If 'simple', return a list of parameter, rCBF and
        absolute CBF values. If 'complicated', return the two classes of
        scipy.optimize.OptimizeResult plus rCBF and CBF values.
    
    Returns
    -------
    If outtype='simple':
        List of 8 floats- selective inversion bias, M0 and T1, non-selective
        inversion bias, M0 and T1, rCBF and CBF. Failed fits produce zeroes.
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
    
    if outtype not in ['simple', 'complicated']: return 'unrecognized outtype'
    
    s0_sel = np.array(s0)[np.array(picker_sel)]
    s0_nonsel = np.array(s0)[np.array(picker_nonsel)]

    r_sel = fair_t1_fit(s0_sel, ti, t1_guess)
    r_nonsel = fair_t1_fit(s0_nonsel, ti, t1_guess)
    
    if r_sel.success and r_nonsel.success:
        
        t1_sel = r_sel.x[2]
        t1_nonsel = r_nonsel.x[2]
        
        rCBF = 100 * (t1_nonsel - t1_sel) / t1_nonsel
        CBF = multiplier * lambda_blood * (
              (t1_nonsel / t1_blood) * ((1 / t1_sel) - (1 / t1_nonsel)))
     
        if outtype == 'simple':
            r = [x for x in r_sel.x] + [x for x in r_nonsel.x] + [rCBF, CBF]
            return [x if x is not None else 0 for x in r]
        else: return [r_sel, r_nonsel, rCBF, CBF]
    
    else:
        
        if outtype == 'simple': return np.repeat(0, 8)
        else:
            return ['r_sel fit ' + r_sel.success,
                    'r_nonsel fit ' + r_nonsel.success,
                    'rCBF nan', 'CBF nan']


def perf_fair_fitter_mp(all_s0, t1_blood, lambda_blood, ti, multiplier,
                        t1_guess, picker_sel, picker_nonsel, ncpu):
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
        constant; an example value is 2800 at 11.7T (Duong ref).
        
    lambda_blood : int or float
        The assumed blood tissue partition coefficient of water in ml per g.
        Empirical constant; usually 0.9.
        
    ti : list of int or float
        The inversion times. Must have the same length as s0.
        
    multiplier : int or float
        The absolute CBF result is initially produced in units of ml per g per
        ms. The multiplier converts to desired units. Usually it is 6000000 to
        convert to ml per 100g per min.    
        
    t1_guess : int or float
        An initial starting value for T1. The mean of TIs is a good guess.
    
    picker_sel : list or numpy array of int
        Vector indicating positions of selectively inverted signals in s0. The
        lengths of picker_sel and picker_sel must sum to the length of s0.
        
    picker_nonsel : list or numpy array of int
        Vector indicating positions of non-selectively inverted signals in s0.
        The lengths of picker_sel and picker_sel must sum to the length of s0.
        
    ncpu : int
        Number of processes to launch in parallel.    
    
    Returns
    -------
    List of all_s0.shape[0] lists, each of 8 floats (see perf_fair_fitter).
    """
    
    # inspiration provided by:
    # https://stackoverflow.com/q/5442910
    # https://stackoverflow.com/a/45276885
    
    # if __name__ == '__main__': this might be dangerous
    if 2 > 1: # to keep the indent in case above line is reinstated
        # cannot use a with statement in python 2.7; may cause problems
        pool = Pool(processes=ncpu)
        return list(tqdm.tqdm(pool.imap(partial(perf_fair_fitter,
                                                t1_blood=t1_blood,
                                                lambda_blood=lambda_blood,
                                                ti=ti,
                                                multiplier=multiplier,
                                                t1_guess=t1_guess,
                                                picker_sel=picker_sel,
                                                picker_nonsel=picker_nonsel),
                                        all_s0), total = len(all_s0)))


def perf_fair_nii_proc(nii_in_fname, nii_out_fname, t1_blood, lambda_blood, ti,
                       multiplier, t1_guess, picker_sel, picker_nonsel, ncpu):
    """
    Wrapper to execute perf_fair_fitter_mp on a NIfTI-1 image.
    
    Parameters
    ----------
    nii_in_fname : str
        Input file path.
        
    nii_out_fname : str
        Output file path.
        
    t1_blood : int or float
        T1 of blood in ms at the acquisition field strength. Empirical
        constant; an example value is 2800 at 11.7T (Duong ref).
        
    lambda_blood : int or float
        The assumed blood tissue partition coefficient of water in ml per g.
        Empirical constant; usually 0.9.
        
    ti : list of int or float
        The inversion times. Must have the same length as s0.
        
    multiplier : int or float
        The absolute CBF result is initially produced in units of ml per g per
        ms. The multiplier converts to desired units. Usually it is 6000000 to
        convert to ml per 100g per min.    
        
    t1_guess : int or float
        An initial starting value for T1. The mean of TIs is a good guess.
    
    picker_sel : list or numpy array of int
        Vector indicating positions of selectively inverted signals in s0. The
        lengths of picker_sel and picker_sel must sum to the length of s0.
        
    picker_nonsel : list or numpy array of int
        Vector indicating positions of non-selectively inverted signals in s0.
        The lengths of picker_sel and picker_sel must sum to the length of s0.
        
    ncpu : int
        Number of processes to launch in parallel.    
    
    Returns
    -------
    NIfTI-1 file saved to nii_out_fname.
    """

    nii_in = nib.load(nii_in_fname)
    in_mat = nii_in.get_data()
    all_s0 = in_mat.reshape((np.product(in_mat.shape[:-1]), in_mat.shape[-1]))
    r = perf_fair_fitter_mp(all_s0, t1_blood, lambda_blood, ti, multiplier,
                            t1_guess, picker_sel, picker_nonsel, ncpu)
    r = np.array(r)
    img = nib.Nifti1Image(np.reshape(r, in_mat.shape[:-1] + (r.shape[1],)),
                          nii_in.get_affine())
    img.to_filename(nii_out_fname)


def perf_fair_niiptbl_proc(nii_in_fname, nii_out_fname, t1_blood, lambda_blood,
                           multiplier, ncpu):
    """
    Wrapper to execute perf_fair_nii_proc supplied with acquisition parameters
    automatically extracted using perf_fair_read_ptbl. So TI and inversion type
    are automatically determined for all signals, as is a T1 starting guess.
    
    Parameters
    ----------
    nii_in_fname : str
        Input file path.
        
    nii_out_fname : str
        Output file path.
        
    t1_blood : int or float
        T1 of blood in ms at the acquisition field strength. Empirical
        constant; an example value is 2800 at 11.7T (Duong ref).
        
    lambda_blood : int or float
        The assumed blood tissue partition coefficient of water in ml per g.
        Empirical constant; usually 0.9.
        
    multiplier : int or float
        The absolute CBF result is initially produced in units of ml per g per
        ms. The multiplier converts to desired units. Usually it is 6000000 to
        convert to ml per 100g per min.
        
    ncpu : int
        Number of processes to launch in parallel.    
    
    Returns
    -------
    NIfTI-1 file saved to nii_out_fname.
    """

    ptbl_dict = perf_fair_read_ptbl(nii_in_fname)
    ti=ptbl_dict['TI']
    picker_sel=ptbl_dict['picker_sel']
    picker_nonsel=ptbl_dict['picker_nonsel']
    t1_guess=np.mean(ptbl_dict['TI'])
    perf_fair_nii_proc(nii_in_fname, nii_out_fname, t1_blood, lambda_blood, ti,
                       multiplier, t1_guess, picker_sel, picker_nonsel, ncpu)

