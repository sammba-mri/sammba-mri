=========== TEMPLATE ===========

MLAmaster.bash

commands to make the template

MRIT1_Mcbconv.bash

convert data to NIfTI

MRIT2_extrcen

bias correction, brain extraction, placing image origin at centre of mass (CoM),
registration to CoM

MRIT3_shr

rigid body registration

MRIT4_aff

affine registration

MRIT5_Qw

non-linear registration from an intial level to a maximum level

MRIT6_Qw

non-linear registration from an intial level to a minimum patch

MRIT7_origproc

apply non-linear registration results to uncorrected images

MRIT8_TBM

calculate and extract per-voxel bulk, shear and vorticity measures from the
non-linear registration results
