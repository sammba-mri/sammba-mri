=========== TEMPLATE ===========

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


MRIT2_extrcen.bash

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


MRIT3_shr.bash

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


MRIT4_aff.bash

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


MRIT5_Qw.bash

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


MRIT6_Qw.bash

Simpler version of MRIT5_Qw.bash with no potential for .aff12.1D input warp, and
going to a minpatch rather than a maxlev.


MRIT7_origproc.bash

Apply non-linear registration results to uncorrected images, and transform any
image in template space to individual spaces. The first result is just out of
curiosity, the second is useful for individualising atlases, which can then be
used for morphological measurements in each individual, so morphometry of input
data.


MRIT8_TBM.bash

Calculate and extract per-voxel bulk, shear and vorticity measures from the
non-linear registration results.