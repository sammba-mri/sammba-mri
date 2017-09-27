==================
sammba-MRI scripts
==================
*Caution: The following scripts are meant to be rewritten in python
as examples after including the recurrent functionalities in the API.
They are provided here as a temporary tool and will no more be maintained.*
    The scripts are in `the dev branch of sammba-MRI
    <https://github.com/sammba-mri/sammba-mri/tree/dev>`_.


Resting state pipeline
----------------------

    The scripts described here are best controlled by a master_projectname.bash script which defines several input parameters (including lots that may not be used though this is not a problem), customised to the needs of the input data and processing requirements. See individual project pages for examples of these (for instance `MLA master
    <https://github.com/sammba-mri/sammba-mri/blob/dev/sammba/projects/MLA/master_MLA.bash>`_)

    The following files/parameters need to be decided/known/prepared/available for full function (a fuller explanation for some of the stranger ones follows):

1. pipeline: script file describing the order and number of processing stages
2. MRIsessionslist: text file listing all the raw data directories
3. analysisdir: directory where all the results will be stored
4. brain: NIfTI of a template brain
5. atlas: NIfTI of a brain atlas
6. mask: NIfTI of a brain mask
7. head: NIfTI of a template head
8. headweight: NIfTI of a mask that covers the brain plus some amount of head
9. basetype: whether to register to the head or brain (defined by above files)
10. dofolderoverwrite: yes/no to overwriting pre-existing directories
11. tmpdir: directory where temporary files are stored
12. registerfunctional: yes/no to coarse registration of non-anatomical data 
13. subpipeline: script that does any per-slice registration
14. Urad: parameter used in bias correction
15. brainvol: approximate brain volume for the species in ml
16. scale: ratio of the above to a human brain

Perfusion-related, need to be supplied even if not used as scripts might not
always function unless they exist (just supply nonsense values if necessary):

17. T1blood: in ms
18. lambda: for the species
19. multiplier: depends on output units
20. T1guess: in ms
21. mccores: maximum number of processors to use

The intended functioning of the pipeline scripts starts with MRIsessionslist.
Each path in it should be to a Bruker ParaVision session folder. This is looped
through: for each folder, a new folder with the same name is created in the
analysisdir. The session folder is searched for DICOM files, which are converted
to NIfTI and stored in the new folder just created in analysisdir. Based on the
NIfTI file names (which are determined using information found in the DICOM),
anatomical, perfusion and BOLD-weighted acquisitions are identified, and given
simpler names including acquisition order, for example anat_n0.nii.gz,
perfFAIREPI_n0.nii.gz, perfFAIREPI_n1.nii.gz, perfFAIREPI_n2.nii.gz and
rs_n0.nii.gz. Unidentified files are left unmodified. A text file called
namechanges.txt is created to record the name changes made. The classic fMRI
processing pipeline can then be, and is, launched. The anatomical scan is
aligned to a template, with options to do this in various ways using the brain
and/or head with weights. Any existing atlas is also transformed to the
individual by application of inverted alignment parameters. Perfusion maps are
calculated for any perfusion data. These perfusion acquisitions and results are
then aligned to the anatomical on a per-slice basis (only a fine correction- in
theory they are already well-aligned given the animal should not have moved,
this is mostly for correction of EPI distortion) before application of the
anatomical-to-template results, such that the perfusion data is now also aligned
to the template. The same spatial correction is made for BOLD-weighted data.

If the analysis was already run in a particular analysisdir, the folders will
already exist, and are not overwritten unless you explicitly demand it. A
deliberate bug is that even this will not actually work, since the way the
script that loops through MRIsessionslist is written, it cannot ever overwrite
the folder. At best it will create a new one with the suffix yyymmdd_hhmmss. You
should just delete old ones manually before running anything.

Certain stages (especially per-slice registration) generate lots of temporary
files. In theory we could just hardcode their storage in /tmp as is normal.
However, the server we use has the OS on a fast SSD with little capacity, so we
added this option so that an alternative location could be chosen (this should
ideally be a fast system drive).

Occasionally we have had to analyse external and non-Bruker data where only
NIfTI files with incorrect orientation information in the header were available.
An option was added to allow an initial coarse registration of functional data
in these cases. The results are not great since functional data is often
distorted, partial coverage and with poor anatomical contrast, but it can work
well enough.

The subpipeline is actually just the name of the subscript that does generalised
per-slice registration. There is only and there will only ever be one, so this
option is not really necessary.

Exactly why there is only control or Urad and not other parameters in
bias-correction is lost in the mists of time. Scale is needed as some paramaters
in 3dAllineate default to values in mm that are suitable for humans but not
animals. To be honest this could be calculated automatically by dividing
brainvol by the size of an average human brain. This change will be made
eventually.

List of pipeline scripts
------------------------

* pipeline_wrapper.bash

Loops through MRIsessionslist, creating a new folder in analysisdir with the
same name as the ParaVision session folder containing the raw data. Launches a
pipeline for each one, saving results in the new folder created in analysisdir.

* pipeline*.bash

List of things to do (which usually means scripts to launch). This is usually
DICOM to NIfTI conversion, MRI protocol identification and simplification of
NIfTI files' names, anatomical to template registration, perfusion calculations,
registration of perfusion data, and finally registration and initial processing
of BOLD-weighted data.

* PVEnDCMtoNIfTI.py

DICOM to NIfTI conversion for Bruker ParaVision enhanced multiframe DICOM files.

* IDsequencetypes.bash

Based on the NIfTI name, guesses what kind of acquisition a file is and renames
it to something simpler, including its acquisition order compared to other files
of the same protocol.

* fixobliquity.bash

Restores obliquity information to the NIfTI header; AFNI often strips this out.

* anattotemplate.bash

Register anatomical image to a template.

* perfFAIREPI.bash

Wrapper to an R script that does the perfusion calculation.

* perfFAIREPItoNIfTI.R

Wrapper to an R source that contains functions for perfusion calculation.

* variousRfunctions.R

Functions that, amongst other things, do the perfusion calculation.

* pipeline_reg.bash

Per-slice registration of a multislice image.

* perfFAIREPI_spatnorm.bash

Register perfusion images to the template using the anatomical registration results.

* tpattern.R

Generate slice timing file for a BOLD-weighted acquisition.

* rs.bash

Slice timing correction and spatial normalization of BOLD-weighted data.

* videomaker.bash

Make videos to help rapidly check spatial registration quality.


Template creation
-----------------
Any given template project needs a master_projectname.bash (a script which
custom organises template generation according to the needs of the project) and
a convert_projectname.bash which does the initial data conversion (also
customised to the needs of the project). See individual project pages for these.

Scripts explained here: MRIT2_extrcen.bash, MRIT3_shr.bash, MRIT4_aff.bash,
MRIT5_Qw.bash, MRIT6_Qw.bash, MRIT7_origproc.bash and MRIT8_TBM.bash. 

These scripts assume that a directory for making the template has already been
created which 1) contains several subdirectories with unique names, one for each
image acquisition that will contribute to the template, and 2) inside each
subdirectory is a single NIfTI-1 format file of the actual contributing image,
named anat.nii.gz. No other files or subdirectories should be present.

The directory and data within it are usually made by a conversion script
specific to the project; data sources differ so much from project to project
that making a universal script would be a lot of hard effort. Assuming it worked
correctly, each anat.nii.gz will be the same for the following (reasonable)
criteria: field-of-view, matrix, approximate orientation (large amounts of tilt
are tolerated, but S-I/D-V, A-P/R-C and L-R axes all need to be the same), plus
the head being reasonably centered and dominant in the image. The template
scripts will likely fail if the conversion script does not achieve this.

Initially, registration is of extracted brains. Once these are reasonably
aligned, whole heads are registered, weighted by masks that, if parameters are
chosen well, include some scalp. The amount of scalp is hopefully enough to help
in differentiating the brain-scalp boundary without including so much head
tissue that it starts to contaminate the registration with the highly variable
head tissue.

It is also useful to read the powerpoint presentation of the mouse lemur
template which visualises the procedure.

List of template scripts
------------------------

* MRIT2_extrcen.bash

This script carries out an initial coarse registration using brain centre of
mass (CoM). It loops through all the subdirectories, processing the anat.nii.gz
in each one. These are bias-corrected then brain-extracted (aided by an
approximate guessed brain volume), with the NIfTI image centre (as defined in
the header) then being set to the CoM of the extracted brain. The head is then
shifted within the image to place the CoM at the image center. This is
effectively a translation-only registration of the anat.nii.gz images to each
other's brain's (as defined by the brain extractor) CoMs. Images are
additionally concatenated at all stages to produce videos that allow rapid
passage through acquisitions to check for quality. The intermediate and final
results are named after their processing steps: Un = 3dUnifize (bias
correction), Bm = brain mask (produced by RATS), Be = brain extracted (using the
brain mask), CC = translated to brain CoM (using CC is a historical accident or
I have forgotten what its meaning is). Individual results are saved in the
respective subdirectories. Group results (means, videos, emptytemplate) are
saved in the main directory.


* MRIT3_shr.bash

Rigid-body registration of CoM brains produced by MRIT2_extrcen.bash, and
application of this registration to CoM heads. This registration requires a
target template. MRIT2_extrcen.bash produces a potential one called
UnBmBeCC_mean.nii.gz (mean of all bias-corrected, brain-extracted, mass-centered
images). Other possibilities include an externally-sourced image or, more
biased, a nicely-aligned individual. Loops through subdirectories and records
videos of all stages as for MRIT2_extrcen.bash. If the script needs to be run
more than once, changing the n option prevents overwriting of files from
previous runs. Al = 3dAllineate, Aa = apply 3dAllineate, shr = shift rotate (the
flag in 3dAllineate to do rigid-body alignment). The Aa images are more for
show, they are not really used. The count mask is useful for looking at brain
extraction efficiency and differences in brain size, and can be used by
MRIT4_aff.bash.


* MRIT4_aff.bash

As for MRIT3_shr.bash but 1) affine registration, 2) inputs are not brains but
aligned heads from MRIT3_shr.bash weighted by a mask (more on the mask later),
3) each resulting registration matrix is concatenated to the corresponding
MRIT3_shr.bash result then directly applied to the CoM brain and head produced
by MRIT2_extrcen.bash, reducing reslice errors in the final result, and 4) there
are two n options: the first is to identify the MRIT3_shr.bash run, the second
to give a number to the affine results to make them easier to tell apart from
the previous linear registrations. It is important to note that there is not a
unique weight mask for each head. Instead, since heads are already roughly
aligned, a single generous weight mask can be used for all of them (typically
the count mask produced by the previous MRIT3_shr.bash), which should also cover
some scalp for most individuals. If brains vary massively in size this may not
work well, but from experience it is fine. In any case, if size is that
variable, the extraction (which assumes brain are around a certain size) would
not have worked well either, messing up the entire strategy. Doing instead an
affine transform of extracted brains (or whole heads weighted by an
individual-specific mask presumably based on brain extraction results) is worse
for two reasons. Firstly, it is less robust if extraction is sometimes poor. The
use of a mask calculated from the results of the whole group effectively
provides safety in numbers. Secondly, it tends to make them too large, at least
with 3dAllineate on small mammal data in the context of how it is being
processed here. This is acceptable for mutual registration but gives the final
template an incorrect size. An alternative is to skip MRIT4_aff.bash and go
straight to MRIT5_Qw.bash, which seems robust enough to handle the rougher
inputs from MRIT3_shr.bash, and does not inflate average brain size.


* MRIT5_Qw.bash

The first non-linear registration. Like MRIT3_shr.bash and MRIT4_aff.bash it
requires a target template. It also requires a weight mask, ideally one that
extends beyond the brain, incorporating some surrounding tissue to help better
define the brain head boundary. The input source images can be initially
transformed prior to registration in three ways: ident (meaning none), an
.aff12.1D, or a warp. Such possibilities can be exploited to avoid building up
reslice errors. If the input transform is a warp, it is concatenated to IDENT
initially; I forget why, I think it is to avoid some weird bug. Source images
should already be quite well-aligned (or quite well-aligned after the initial
transform) to the template. Registration is from a given level to another,
though any level below a patch size of 25 will not be done, see 3dQwarp help for
further detail.


* MRIT6_Qw.bash

Simpler version of MRIT5_Qw.bash with no potential for .aff12.1D input warp, and
going to a minpatch rather than a maxlev.


* MRIT7_origproc.bash

Apply non-linear registration results to uncorrected images, and transform any
image in template space to individual spaces. The first result is just out of
curiosity, the second is useful for individualising atlases, which can then be
used for morphological measurements in each individual, so morphometry of input
data.


* MRIT8_TBM.bash

Calculate and extract per-voxel bulk, shear and vorticity measures from the
non-linear registration results.

