# -*- coding: utf-8 -*-

#functions to convert Bruker Paravision enhanced multiframe DICOM files into
#NIfTI-1 format.

#depends on the presence of a compiled version of dcmdump (part of OFFIS dcmtk;
#http://dcmtk.org/dcmtk.php.en and http://support.dcmtk.org/docs/dcmdump.html)
#which does the initial parsing: extraction of the header/metadata as text and
#the image data as a raw vector

#so this is effectively a python wrapper to dcmdump, processing and passing
#its output to nibabel and a text file

#useful documents include DICOM spec C.7.6.2.1.1 and 10.7.1.3
#common.py of dicom2nifiti by Arne Brys, icometrix saved me when it comes to
#specifiying the affine for nibabel! see http://dicom2nifti.readthedocs.io.
#he effectively did the hard part, interpreting nibabel's DICOM tutorial for me

#exports to a .txt text file some meta-data useful to the processing of certain
#protocols such as perfusion, T1/T2 mapping and diffusion

#only tested on PV6 .dcm files and a limited number of sequences. little or no
#error-checking. there are a lot of circumstances where this converter will
#fail or be sub-optimal

#library uses:
#numpy: affine and image matrix manipulation
#os: standardise file/folder paths, file deletion
#subprocess: interface dcmdump
#itertools and pandas: image frame and parameter sorting
#nibabel: NIfTI-1 format


#%%

import numpy as np
import os
import subprocess
import itertools
import pandas
import nibabel

#%%

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
    dcmdump_path is the path to the compiled dcmdump
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
#%%

    #regularise/standardise input paths
    input_paths = [dcmdump_path, EnDCM, save_directory]
    input_paths = [os.path.abspath(input_path) for input_path in input_paths]
    dcmdump_path, EnDCM, save_directory = input_paths

#%%

    #read EnDCM header/metadata using dcmdump
    #+L ensures long tags are fully printed
    dcmdump_output = subprocess.check_output([dcmdump_path, EnDCM, '+L'])

#%%

    #fields that are likely to be empty or multiple are awkward to dynamically
    #declare, so do it here rather than in the parser loop below
    TR = []  # is not always in DICOM
    protocol = []  # cos of DKI
    bruker_sequence = []  # cos of DKI
    diffdir = []
    TIs = []
    bval = []
    bvalXX = []
    bvalXY = []
    bvalXZ = []
    bvalYY = []
    bvalYZ = []
    bvalZZ = []
    IPPs = []
    ISPnums = []
    repno = []
    DIVs = []
    frame_comments = []
    
    #via subprocess, dcmdump produces numerous lines of string output.
    #for some parameters such as SeriesDate, there should hopefully be only one
    #line. for others, such as InStackPositionNumber, there may be several.
    #three example lines are shown below:
    
    #'(0008,0021) DA [20160219]                               #   8, 1 SeriesDate'
    #'(0018,9079) FD 35                                       #   8, 1 InversionTimes'
    #'(0020,0032) DS [-9\\-7.417465618\\1.68366613]             #  26, 3 ImagePositionPatient'
    
    #fields interpreted as strings are delimited with square brackets
    #(incidentally, several numeric fields get interpreted as strings, I do not
    #know why, and I assume it does not matter). in the loop below, each line
    #is parsed to identify what field it represents and extract the value

#%%

    for line in dcmdump_output.splitlines():
    
        #pline = processed line. some lines have leading whitespace
        pline = line.lstrip()
        
        if len(pline) >= 11:
        
            if  pline[0] == '(' and pline[5] == ',' and pline[10] == ')':

                tag = pline[0:11]
                val = pline[:pline.rfind('#')].rstrip()[valstart:]
                
                #remove [] that dcmdump delimits strings with
                if val[0] == '[':
                    val = val[1:-1]
                
                vals = val.split(splt1)
                
                #turn single value lists into single values
                if type(vals) == list:
                    if len(vals) == 1:
                        vals = vals[0]
            
                #assign output to variable names
                if tag == '(0008,0021)':  #Series Date
                    acqdate = vals
                if tag == '(0008,0031)':  #Series Time
                    acqtime = vals
                if tag == '(0010,0020)':  #Patient ID​
                    patID = vals
                if tag == '(0018,0023)':  #MR Acquisition Type
                    acqdims = vals
                if tag == '(0018,0050)':  #Slice Thickness​
                    thk = float(vals)
                if tag == '(0018,0080)':  #Repetition Time
                    TR = float(vals)
                if tag == '(0018,1030)':  #Protocol Name
                    protocol = vals
                if tag == '(0018,5100)':  #Patient Position
                    patpos = vals
                if tag == '(0018,9005)':  #Pulse Sequence Name
                    bruker_sequence = vals
                if tag == '(0018,9075)':  #Diffusion Directionality​
                    diffdir.append(vals)
                if tag == '(0018,9079)':  #Inversion Times
                    TIs.append(float(vals))
                if tag == '(0018,9087)':  #Diffusion b-value
                    bval.append(float(vals))
                if tag == '(0018,9602)':  #Diffusion b-value XX
                    bvalXX.append(float(vals))
                if tag == '(0018,9603)':  #Diffusion b-value XY​
                    bvalXY.append(float(vals))
                if tag == '(0018,9604)':  #Diffusion b-value XZ
                    bvalXZ.append(float(vals))
                if tag == '(0018,9605)':  #Diffusion b-value YY
                    bvalYY.append(float(vals))
                if tag == '(0018,9606)':  #Diffusion b-value YZ
                    bvalYZ.append(float(vals))
                if tag == '(0018,9607)':  #Diffusion b-value Z
                    bvalZZ.append(float(vals))
                if tag == '(0020,0032)':  #Image Position (Patient)
                    IPPs.append([float(i) for i in vals])
                if tag == '(0020,0037)':  #Image Orientation (Patient)
                    cosines = [float(i) for i in vals]
                if tag == '(0020,9057)':  #In-Stack Position Number
                    ISPnums.append(int(vals))
                if tag == '(0020,9128)':  #Temporal Position Index​
                    repno.append(int(vals))
                if tag == '(0020,9157)':  #Dimension Index Values
                    DIVs.append([int(i) for i in vals])
                if tag == '(0020,9158)':  #Frame Comments
                    frame_comments.append(vals)                
                if tag == '(0028,0008)':  #Number of Frames
                    frames = int(vals)
                if tag == '(0028,0010)':  #Rows
                    rows = int(vals)
                if tag == '(0028,0011)':  #Columns
                    cols = int(vals)
                if tag == '(0028,0030)':  #Pixel Spacing
                    pixspac = [float(i) for i in vals]

                #might be useful in the future
                #(0028,0100)  # Bits Allocated
                #(0028,0101)  # Bits Stored
                #(0028,0102)  # High Bit
                #(0028,0103)  # Pixel Representation

#%%    
    
    #if they have length zero (or just unequal to ISPnums, though no idea how
    #that could be possible), populate vectors that will be included in ptbl
    #(parameter table). there must be a more efficient way
    if len(diffdir) != len(ISPnums):
        diffdir = list(itertools.repeat('NA', len(ISPnums)))
    if len(TIs) != len(ISPnums):
        TIs = list(itertools.repeat('NA', len(ISPnums)))
    if len(bval) != len(ISPnums):
        bval = list(itertools.repeat('NA', len(ISPnums)))
    if len(bvalXX) != len(ISPnums):
        bvalXX = list(itertools.repeat('NA', len(ISPnums)))
    if len(bvalXY) != len(ISPnums):
        bvalXY = list(itertools.repeat('NA', len(ISPnums)))
    if len(bvalXZ) != len(ISPnums):
        bvalXZ = list(itertools.repeat('NA', len(ISPnums)))
    if len(bvalYY) != len(ISPnums):
        bvalYY = list(itertools.repeat('NA', len(ISPnums)))
    if len(bvalYZ) != len(ISPnums):
        bvalYZ = list(itertools.repeat('NA', len(ISPnums)))
    if len(bvalZZ) != len(ISPnums):
        bvalZZ = list(itertools.repeat('NA', len(ISPnums)))
    if len(repno) != len(ISPnums):
        repno = list(itertools.repeat('NA', len(ISPnums)))
    if len(frame_comments) != len(ISPnums):
        frame_comments = list(itertools.repeat('NA', len(ISPnums)))
    
    #ptbl = parameter table
    ptbl = pandas.DataFrame(data = {'repno': repno, 'slice': ISPnums,
        'diffdir': diffdir, 'TI': TIs, 'bval': bval, 'bvalXX': bvalXX,
        'bvalXY': bvalXY, 'bvalXZ': bvalXZ, 'bvalYY': bvalYY, 'bvalYZ': bvalYZ,
        'bvalZZ': bvalZZ, 'slicepos': IPPs, 'FC': frame_comments})

    DIVsdf = pandas.DataFrame({'DIVs': DIVs})
    DIVsdf = pandas.DataFrame([x[0] for x in DIVsdf.values])
    DIVsdf = DIVsdf.rename(columns = lambda x: 'DIV' + str(x))

    dfs = [DIVsdf, ptbl]
    ptbl = pandas.concat(dfs, axis = 1)


#%%

    slices = max(ISPnums) #maybe a bit dangerous
    
    #use check_output rather than call to avoid stdout filling up the terminal
    rawfile = os.path.join(save_directory, 'EnIm1.dcm.0.raw')
    if os.path.exists(rawfile):
        os.remove(rawfile)
    sink = subprocess.check_output([dcmdump_path, EnDCM, '+W', save_directory])
    rawarray = np.fromfile(rawfile, dtype = np.int16)
    os.remove(rawfile)
    
    #rawarray is actually just a vector, need to reshape into a run of frames
    rawarray = np.reshape(rawarray, (frames, rows, cols))
    rawarray = np.reshape(rawarray, (int(frames / slices), slices, rows, cols))
    rawarray = np.transpose(rawarray, (3,2,1,0))


#%%

    #fsp = first slice position, lsp = last slice position
    fsp = np.array(ptbl.slicepos[ptbl.slice == 1].tolist()[0])
    lsp = np.array(ptbl.slicepos[ptbl.slice == slices].tolist()[0])
    
    #https://en.wikipedia.org/wiki/Euclidean_distance
    #not used yet
    flspdiff = fsp - lsp
    eucliddist = (flspdiff[0]**2 + flspdiff[1]**2 + flspdiff[2]**2)**0.5
    slicegap = eucliddist / (slices -1)
   
    if slices == 1:  # single slice
        step = [0, 0, -1]
    else:
        step = (fsp - lsp) / (1 - slices)
    
    affine = np.matrix(
    [[-cosines[0]*pixspac[1], -cosines[3]*pixspac[0], -step[0], -fsp[0]],
     [-cosines[1]*pixspac[1], -cosines[4]*pixspac[0], -step[1], -fsp[1]],
     [ cosines[2]*pixspac[1],  cosines[5]*pixspac[0],  step[2],  fsp[2]],
     [                     0,                      0,        0,      1]])
    
    affineident = np.eye(4, dtype=int)
    
    #not sure if any of these patpos specs is really correct
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

#%%

    header = nibabel.Nifti1Header()
    header.set_xyzt_units('mm', 'msec')
    header.set_data_dtype(np.int16)
    if type(TR) != list:
        header['pixdim'][4] = TR    
    img = nibabel.Nifti1Image(rawarray, affine, header = header)

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

#%%

    ptbl.to_csv(path_or_buf = os.path.join(save_directory, NIIname[:-7] +
                                           '_ptbl.txt'),
                sep = '\t', index = 0)
