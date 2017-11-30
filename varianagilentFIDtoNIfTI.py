import numpy as np
import sys
import nmrglue
import nibabel

def rotate_affine(angle, axis):
    """
    rotate an affine matrix by the given angle in degrees about the given axis
    x, y or z used for orientation correction
    https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.28_z-y.E2.80.99-x.E2.80.B3_intrinsic.29_.E2.86.92_Rotation_matrix
    """
    a = angle * np.pi / 180
    s = np.sin(a)
    c = np.cos(a)
    if axis == 'x':
        return np.matrix([[ 1,  0,  0,  0],
                          [ 0,  c, -s,  0],
                          [ 0,  s,  c,  0],
                          [ 0,  0,  0,  1]])
    if axis == 'y':
        return np.matrix([[ c,  0,  s,  0],
                          [ 0,  1,  0,  0],
                          [-s,  0,  c,  0],
                          [ 0,  0,  0,  1]])
    if axis == 'z':
        return np.matrix([[ c, -s,  0,  0],
                          [ s,  c,  0,  0],
                          [ 0,  0,  1,  0],
                          [ 0,  0,  0,  1]])


#%%
fiddir = sys.argv[1]
savefilename = sys.argv[2]
fftzpsize = int(sys.argv[3])

#%%
dic, data = nmrglue.varian.read(fiddir, as_2d = 1)

#%%

pss = np.array([float(value) for value in dic['procpar']['pss']['values']])
pssorder = pss.argsort()

pelist = np.array([int(value) for value in dic['procpar']['pelist']['values']])
peorder = np.argsort(pelist - min(pelist))

#%%

nro = dic['np']/2 #number of readout lines
npe = pelist.size #number of phase encodes. in the procpar I think this is nv
#so (better?) alternative is maybe int(dic['procpar']['nv']['values'][0])
ETL = int(dic['procpar']['etl']['values'][0]) #echo train length
nslices = pss.size #number of slices
nsegments = npe / ETL
nblocks = dic['nblocks']
ntraces = dic['ntraces']

rawarray = np.empty((fftzpsize, fftzpsize, nslices, nblocks))

nrozp = (fftzpsize - nro) / 2
npezp = (fftzpsize - npe) / 2

#%%
for block in range(nblocks):
    for slicen in pssorder:
        #block and slice are not tiled or repeated as we only want one of each
        linepicker = (np.repeat(range(nsegments), ETL) * nslices * ETL) + (
                     block * ntraces) + (
                     slicen * ETL) + (
                     np.tile(range(ETL), nsegments))
        kspace = np.pad(data[linepicker][peorder], (nrozp, npezp),
                        'constant', constant_values = (0, 0))
        rawarray[:, :, slicen, block] = np.absolute(np.fft.fftshift(np.fft.fft2(kspace)))

rawarray = np.average(rawarray, axis = 3)
    
#%%
#the x, y, z, ro, pe and slice relationships are prob all fucked here

dx = float(dic['procpar']['lro']['values'][0]) * 10 / fftzpsize #in cm not mm!
dy = float(dic['procpar']['lpe']['values'][0]) * 10 / fftzpsize #in cm not mm!
dz = float(dic['procpar']['thk']['values'][0])

affine = np.matrix(
    [[dx,  0,  0, 0],
     [ 0, dz,  0, 0],
     [ 0,  0, dy, 0],
     [ 0,  0,  0, 1]])

affine = affine * rotate_affine(270, 'x') * rotate_affine(180, 'y') * rotate_affine(180, 'z')

header = nibabel.Nifti1Header()
header.set_xyzt_units('mm', 'msec')
header.set_data_dtype(np.int16)
img = nibabel.Nifti1Image(rawarray, affine, header = header)
img.to_filename(savefilename)

