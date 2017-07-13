=========== PIPELINE ===========

PVEnDCMtoNIfTI.py

DICOM to NIfTI conversion for Bruker Paravision multiframe DICOMs

pipeline_wrapper.bash

loop a pipeline script across several sessions given in a text file

pipeline*.bash

list of things to do to converted data

IDsequencetypes.bash

based on the NIfTI name, guesses what kind of acquisition a file is and renames
it to something simple

fixobliquity.bash

restores obliquity information to the NIfTI header; AFNI often strips this out

anattotemplate.bash

register anatomical image to a template

perfFAIREPI.bash

wrapper to an R script that does the perfusion calculation

perfFAIREPItoNIfTI.r

wrapper to an R source that contains functions for perfusion calculation

perfFAIREPI+T1_T2map_RARE_fitters.r

functions to do perfusion calculation

pipeline_reg.bash

per-slice registration of a multislice image

perfFAIREPI_spatnorm.bash

register perfusion images to template

tpattern.r

generate slice timing fole for a functional acquisition

rs.bash

slice timing correction and spatial normalization of functional MRI acquisition

videomaker.bash

make videos to help rapidly check spatial registration quality
