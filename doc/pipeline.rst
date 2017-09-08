=========== PIPELINE ===========

The scripts described here are best controlled by a master_projectname.bash
script which defines several input parameters (including lots that may not be
used though this is not a problem), customised to the needs of the input data
and processing requirements. See individual project pages for examples of these.

The following files/parameters need to be decided/known/prepared/available for
full function (a fuller explanation for some of the stranger ones follows):
1) pipeline: script file describing the order and number of processing stages
2) MRIsessionslist: text file listing all the raw data directories
3) analysisdir: directory where all the results will be stored
4) brain: NIfTI of a template brain
5) atlas: NIfTI of a brain atlas
6) mask: NIfTI of a brain mask
7) head: NIfTI of a template head
8) headweight: NIfTI of a mask that covers the brain plus some amount of head
9) basetype: whether to register to the head or brain (defined by above files)
10) dofolderoverwrite: yes/no to overwriting pre-existing directories
11) tmpdir: directory where temporary files are stored
12) registerfunctional: yes/no to coarse registration of non-anatomical data 
13) subpipeline: script that does any per-slice registration
14) Urad: parameter used in bias correction
15) brainvol: approximate brain volume for the species in ml
16) scale: ratio of the above to a human brain

Perfusion-related, need to be supplied even if not used as scripts might not
always function unless they exist (just supply nonsense values if necessary):
17) T1blood: in ms
18) lambda: for the species
19) multiplier: depends on output units
20) T1guess: in ms
21) mccores: maximum number of processors to use

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


pipeline_wrapper.bash

Loops through MRIsessionslist, creating a new folder in analysisdir with the
same name as the ParaVision session folder containing the raw data. Launches a
pipeline for each one, saving results in the new folder created in analysisdir.


pipeline*.bash

List of things to do (which usually means scripts to launch). This is usually
DICOM to NIfTI conversion, MRI protocol identification and simplification of
NIfTI files' names, anatomical to template registration, perfusion calculations,
registration of perfusion data, and finally registration and initial processing
of BOLD-weighted data.


PVEnDCMtoNIfTI.py

DICOM to NIfTI conversion for Bruker ParaVision enhanced multiframe DICOM files.
Also outputs a text file of useful metadata.

IDsequencetypes.bash

Based on the NIfTI name, guesses what kind of acquisition a file is and renames
it to something simpler, including its acquisition order compared to other files
of the same protocol.


fixobliquity.bash

Restores obliquity information to the NIfTI header; AFNI often strips this out.


anattotemplate.bash

Register anatomical image to a template.


perfFAIREPI.bash

Wrapper to an R script that does the perfusion calculation.


perfFAIREPItoNIfTI.R

Wrapper to an R source that contains functions for perfusion calculation.


variousRfunctions.R

Functions that, amongst other things, do the perfusion calculation.


pipeline_reg.bash

Per-slice registration of a multislice image.


perfFAIREPI_spatnorm.bash

Register perfusion images to the template using the anatomical registration
results.


tpattern.R

Generate slice timing file for a BOLD-weighted acquisition.


rs.bash

Slice timing correction and spatial normalization of BOLD-weighted data.


videomaker.bash

Make videos to help rapidly check spatial registration quality.
