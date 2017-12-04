# Author: Nachiket Nadkarni, 2017
# License: CeCILL-B

import numpy as np
from sammba.externals.nmrglue import varian
import nibabel
from .utils import _rotate_affine


def fid_to_nii(fid_directory, save_filename, fftzpsize):
    """ Converts FID to NIFTI.

    Parameters
    ----------
    fid_directory : str
        Path to the FID directory.

    save_filename : str
        Path save the extracteed NIFTI image.

    fftzpsize : int

    Returns
    -------
    save_filename : str
        Path to the created NIFTI image.
    """
    dic, data = varian.read(fid_directory, as_2d=1)
    pss = np.array([float(value) for value in dic['procpar']['pss']['values']])
    pssorder = pss.argsort()
    pelist = np.array(
        [int(value) for value in dic['procpar']['pelist']['values']])
    peorder = np.argsort(pelist - min(pelist))

    nro = dic['np'] / 2  # number of readout lines
    npe = pelist.size  # number of phase encodes. in the procpar I think this
    # is nv so (better?) alternative is maybe
    # int(dic['procpar']['nv']['values'][0])
    ETL = int(dic['procpar']['etl']['values'][0])  # echo train length
    nslices = pss.size  # number of slices
    nsegments = npe / ETL
    nblocks = dic['nblocks']
    ntraces = dic['ntraces']

    rawarray = np.empty((fftzpsize, fftzpsize, nslices, nblocks))

    nrozp = (fftzpsize - nro) / 2
    npezp = (fftzpsize - npe) / 2

    for block in range(nblocks):
        for slicen in pssorder:
            # block and slice are not tiled or repeated as we only want
            # one of each
            linepicker = (np.repeat(range(nsegments), ETL) * nslices * ETL) + \
                (block * ntraces) + (slicen * ETL) + \
                (np.tile(range(ETL), nsegments))
            kspace = np.pad(data[linepicker][peorder], (nrozp, npezp),
                            'constant', constant_values=(0, 0))
            rawarray[:, :, slicen, block] = np.absolute(
                np.fft.fftshift(np.fft.fft2(kspace)))

    rawarray = np.average(rawarray, axis=3)

    # the x, y, z, ro, pe and slice relationships are prob all mixed here
    dx = float(dic['procpar']['lro']['values'][0]) * 10 / fftzpsize  # in cm not mm!
    dy = float(dic['procpar']['lpe']['values'][0]) * 10 / fftzpsize  # in cm not mm!
    dz = float(dic['procpar']['thk']['values'][0])

    affine = np.matrix(
        [[dx,  0,  0, 0],
         [0, dz,  0, 0],
         [0,  0, dy, 0],
         [0,  0,  0, 1]])

    affine = affine * _rotate_affine(270, 'x') * _rotate_affine(180, 'y') * \
        _rotate_affine(180, 'z')

    header = nibabel.Nifti1Header()
    header.set_xyzt_units('mm', 'msec')
    header.set_data_dtype(np.int16)
    img = nibabel.Nifti1Image(rawarray, affine, header=header)
    img.to_filename(save_filename)
    return save_filename
