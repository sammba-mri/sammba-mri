# -*- coding: utf-8 -*-

#convert Bruker Paravision enhanced (multiframe) DICOM files to NIfTI-1 format
#recursively searches a given folder for EnIm*.dcm files and saves *.nii to
#a(nother) given folder. names the .nii based on several tags found in the .dcm
#for FAIR_EPI also extracts TIs and saves as a .txt with same name as the .nii;
#will eventually do the same for some other sequences (such as T1/T2 mapping)
#only tested on PV6 .dcm and a limited number of sequences. little or no
#error-checking. depends on dcmdump, see http://dcmtk.org/dcmtk.php.en

#%%

#numpy is used for general image matrix and affine transform management
#os handles file/folder paths and file deletion
#subprocess interfaces to dcmdump
#itertools and pandas generate and manage parameter tables used to sort image
#frames from the dicom file into a "correct" order
#nibabel creates the NIfTI-1 file
#sys grabs shell flags to the python script

import numpy as np
import os
import subprocess
import itertools
import pandas
import nibabel
import sys

def rotate_affine(angle, axis):
    """
    rotate an affine matrix by the given angle in degrees about the given axis
    x, y or z used for orientation correction (I disagree with Paravision about
    the definition of the superior-inferior and anterior-posterior axes)
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
    
def EnDCM_to_NII(dcmdump_path, EnDCM, save_directory, SIAPfix, valstart, splt1, splt2):
    """
    the actual converter
    assumes int16; this will be changed in the future.
    dcmdump_path is the path to the compiled dcmdump (part of OFFIS dcmtk),
    which does the parsing and (later on) extraction of the image matrix.
    EnDCM is the DICOM file (typically named EnIm1.dcm).
    save_directory is what it says it is
    SIAPfix is a yes or no to swapping the superior-inferior and
    anterior-posterior axes to be how I think they should be. In rodents,
    Paravision sets dorsal-ventral as AP and rostral-caudal as SI. I think they
    should be the other way round https://en.wikipedia.org/wiki/Anatomical_terms_of_location#Main_terms
    valstart and splt1 are parameters used to determine how to parse dcmdump
    output, which differs according to the os. for unix these are respectively
    15 and \\, for windows 17 and \\\\. splt2 splits the EnDCM path to help
    determine the experiment folder/number, which is always 3 levels higher.
    of course, file/folder separators differ according to os, with unix being /
    and windows \\.
    """

    #regularise/standardise input paths
    input_paths = [dcmdump_path, EnDCM, save_directory]
    input_paths = [os.path.abspath(input_path) for input_path in input_paths]
    dcmdump_path, EnDCM, save_directory = input_paths

    #parse EnDCM using dcmdump
    fields = subprocess.check_output([dcmdump_path, EnDCM, '+L',
        '+P', '0008,0021',   # SeriesDate
        '+P', '0008,0031',   # SeriesTime
        '+P', '0010,0020',   # PatientID
        '+P', '0018,0023',   # MRAcquisitionType
        '+P', '0018,0050',   # SliceThickness
        '+P', '0018,0080',   # RepetitionTime
        '+P', '0018,1030',   # ProtocolName
        '+P', '0018,5100',   # PatientPosition
        '+P', '0018,9005',   # PulseSequenceName
        '+P', '0018,9079',   # InversionTimes
        '+P', '0020,0032',   # ImagePositionPatient
        '+P', '0020,0037',   # ImageOrientationPatient
        '+P', '0020,9057',   # InStackPositionNumber
        '+P', '0020,9128',   # TemporalPositionIndex
        '+P', '0020,9158',   # FrameComments
        '+P', '0028,0008',   # NumberOfFrames
        '+P', '0028,0010',   # Rows
        '+P', '0028,0011',   # Columns
        '+P', '0028,0030'])  # PixelSpacing
    
    #might be useful in the future
    #(0028,0100)  # BitsAllocated
    #(0028,0101)  # BitsStored
    #(0028,0102)  # HighBit
    #(0028,0103)  # PixelRepresentation

    TR = []  # is not always in DICOM
    protocol = []  # cos of DKI
    bruker_sequence = []  # cos of DKI
    TIs = []
    frame_comments = []
    repno = []
    IPPs = []
    ISPnums = []
    
    for line in fields.splitlines():
    
        #dcmdump output is in bytes, easier to parse strings
        #I am guessing that a consequence of the DICOM spec and/or dcmdump mean 
        #that the field value ALWAYS starts at the same position (valstart)
        line = str(line)
        valsfieldsep = '#  '
        splitline = line[valstart:].split(valsfieldsep)
        if len(splitline) != 2:
            print('assumed delimiter "' + valsfieldsep + '" not present or ' +
                  'present more than once in ' + line)
        vals, field = splitline
        vals = vals.rstrip()
        
        #remove [] that dcmdump delimits strings with
        if vals[0] == '[':
            vals = vals[1:-1]
        
        vals = vals.split(splt1)
        
        #turn single value lists into single values
        if type(vals) == list:
            if len(vals) == 1:
                vals = vals[0]
    
        #assign output to variable names
        if 'SeriesDate' in field:  #                                  0008,0021
            acqdate = vals
        if 'SeriesTime' in field:  #                                  0008,0031
            acqtime= vals
        if 'PatientID' in field:  #                                   0010,0020
            patID = vals
        if 'MRAcquisitionType' in field:  #                           0018,0023
            acqdims = vals
        if 'SliceThickness' in field:  #                              0018,0050
            thk = float(vals)
        if 'RepetitionTime' in field:  #                              0018,0080
            TR = float(vals)
        if 'ProtocolName' in field:  #                                0018,1030
            protocol = vals
        if 'PatientPosition' in field:  #                             0018,5100
            patpos = vals
        if 'PulseSequenceName' in field:  #                           0018,9005
            bruker_sequence = vals        
        if 'InversionTimes' in field:  #                              0018,9079
            TIs.append(float(vals))
        if 'FrameComments' in field:  #                               0018,9158
            frame_comments.append(vals)
        if 'ImagePositionPatient' in field:  #                        0020,0032
            IPPs.append([float(val) for val in vals])
        if 'ImageOrientationPatient' in field:  #                     0020,0037
            cosines = [float(val) for val in vals]
        if 'InStackPositionNumber' in field:  #                       0020,9057
            ISPnums.append(int(vals))
        if 'TemporalPositionIndex' in field:  #                       0020,9128
            repno.append(int(vals))
        if 'NumberOfFrames' in field:  #                              0028,0008
            frames = int(vals)
        if 'Rows' in field:  #                                        0028,0010
            rows = int(vals)
        if 'Columns' in field:  #                                     0028,0011
            cols = int(vals)
        if 'PixelSpacing' in field:  #                                0028,0030
            pixspac = [float(val) for val in vals]
    
    #InStackPositionNumber (ISPnum) is not necessarily needed. sorted unique 
    #slice positions (used for calculating interslice gap and the affine matrix,
    #plus doing any potential data reordering) can be determined from the slice 
    #position information itself. however, that sounds like too much hard work
    
    uIPPs = []  # u for unique
    slices = max(ISPnums)
    #I think this handles the (unlikely) possibility of InStackPositionNumber
    #not being ordered
    for ISPnum in range(1, slices + 1):
        sliceindex = ISPnums.index(ISPnum)
        uIPPs.append(IPPs[sliceindex])
    uIPPs = np.array(uIPPs)  # find an append method to avoid this line
    
    #https://en.wikipedia.org/wiki/Euclidean_distance
    #not used yet
    IPPsdiffs = uIPPs[:-1] - uIPPs[1:]
    eucliddists = (IPPsdiffs[:,0]**2 + IPPsdiffs[:,1]**2 + IPPsdiffs[:,2]**2)**0.5
    if slices > 1:
        slicegap = sum(eucliddists) / len(eucliddists)
    
    #Documentation
    #C.7.6.2.1.1
    #10.7.1.3 (no letter beforehand!!)
    #common.py of dicom2nifiti by Arne Brys, icometrix
    #http://dicom2nifti.readthedocs.io
    
    if len(uIPPs) == 1:  # single slice
        step = [0, 0, -1]
    else:
        step = (uIPPs[0] - uIPPs[-1]) / (1 - len(uIPPs))
    #if, for some evil reason, the slices were not stored in order, I think the
    #above line would fail. the solution would be to use InStackPositionNumber,
    #or being even more careful, calculating which are the two outer slices
    
    affine = np.matrix(
    [[-cosines[0]*pixspac[1], -cosines[3]*pixspac[0], -step[0], -uIPPs[0,0]],
     [-cosines[1]*pixspac[1], -cosines[4]*pixspac[0], -step[1], -uIPPs[0,1]],
     [ cosines[2]*pixspac[1],  cosines[5]*pixspac[0],  step[2],  uIPPs[0,2]],
     [                     0,                      0,        0,           1]])
    
    affineident = np.matrix(
    [[1, 0, 0, 0],
     [0, 1, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 1]])
	
    if patpos == 'HFS':
        ppaff = rotate_affine(180, 'y')
    #not sure this is correct: a reflection may be necessary too
    if patpos == 'FFS':
        ppaff = affineident
    if patpos == 'HFP':
        ppaff = rotate_affine(180, 'x')
    
    if SIAPfix == 'yes':
        affine = rotate_affine(270,'x') * rotate_affine(180,'y') * rotate_affine(180,'z') * ppaff * affine
    else:
        if SIAPfix != 'no':
            print('SIAPfix is neither yes nor no, assuming no')
    
	#next a table of certain frame parameters is generated that can be sorted
	#by pandas. once correctly ordered, the vector describing the original
	#positions within the DICOM image matrix is later used to extract and
	#reorder the frames correctly
	
    #there must be a more efficient way
    if len(TIs) != len(ISPnums):
        TIs = list(itertools.repeat(0, len(ISPnums)))
    if len(frame_comments) != len(ISPnums):
        frame_comments = list(itertools.repeat(0, len(ISPnums)))
    if len(repno) != len(ISPnums):
        repno = list(itertools.repeat(0, len(ISPnums)))
    
    paramtable = pandas.DataFrame(
        data = {'TI': TIs, 'FC': frame_comments, 'slice': ISPnums, 'repno': repno})
    paramtable['FC'] = paramtable['FC'].astype('category')
    
    if bruker_sequence == 'Bruker:FAIR_EPI':
        paramtable['FC'].cat.reorder_categories(
            ['Selective Inversion', 'Non-selective Inversion'], inplace=True)
    
    paramtable = paramtable.sort_values(list(['repno', 'TI', 'FC', 'slice']))
    
    #use check_output rather than call to avoid stdout filling up the terminal
    rawfile = os.path.join(save_directory, 'EnIm1.dcm.0.raw')
    if os.path.exists(rawfile):
        os.remove(rawfile)
    sink = subprocess.check_output([dcmdump_path, EnDCM, '+W', save_directory])
    rawarray = np.fromfile(rawfile, dtype = np.int16)
    os.remove(rawfile)
    
	#rawarray is actually just a vector, need to reshape into a run of frames
    rawarray = np.reshape(rawarray, (frames, rows, cols))
    rawarray = rawarray[paramtable.index.values] # reorder frames using table
	#make into 4D file. some acquisitions (such as simultaneous T1 and T2
	#mapping) are >4D, and the NIFTI-1 standard can handle up to 7?, but I do
	#not know any visualization software that can, so stick to 4D with a
	#sensible use of temporal order for supplementary dimensions
    rawarray = np.reshape(rawarray, (int(frames / slices), slices, rows, cols))
    rawarray = np.transpose(rawarray, (3,2,1,0))
    
    header = nibabel.Nifti1Header()
    header.set_xyzt_units('mm', 'msec')
    header.set_data_dtype(np.int16)
    if type(TR) != list:
        header['pixdim'][4] = TR    
    img = nibabel.Nifti1Image(rawarray, affine, header=header)

    bf = 'bf' + EnDCM.split(splt2)[-5]  # Paravision experiment folder number
    if SIAPfix == 'yes':
        SIAPfixres = 'fixedSIAP'
    else:
        SIAPfixres = 'origSIAP'
    #the joins deal with a problem of DKI recon where one of the patID and
    #protocol seems to be a list
    NIIname = (''.join(patID) + '__' +
               ''.join(acqdate) + '__' +
               ''.join(acqtime) + '__' +
               ''.join(protocol) + '__' +
               ''.join(SIAPfixres) + '__' +
               ''.join(bf) + '.nii.gz')

    img.to_filename(os.path.join(save_directory, NIIname))

    #need to find a better way to extract TIs
    if bruker_sequence == 'Bruker:FAIR_EPI': 
        f = open(os.path.join(save_directory, NIIname[:-7] + '_TIs.txt'), 'w')        
        f.write(str(paramtable.TI[paramtable.slice == 1][paramtable.FC == 'Selective Inversion'].tolist())[1:-1].replace(' ',''))
        f.close()

#%%

dcmdump_path = sys.argv[1]
sessdir = sys.argv[2]
save_directory = sys.argv[3]
SIAPfix = sys.argv[4]

#%%

if os.name == 'nt':
    valstart = 17
    splt1 = '\\\\'
    splt2 = '\\'
if os.name == 'posix':
    valstart = 15
    splt1 = '\\'
    splt2 = '/'

for root, dirs, files in os.walk(sessdir):
    for file in files:
        if file.startswith('EnIm'):
            if file.endswith('.dcm'):
                EnDCM_to_NII(dcmdump_path,
                             os.path.join(root, file),
                             save_directory,
                             SIAPfix,
                             valstart,
                             splt1,
                             splt2)
