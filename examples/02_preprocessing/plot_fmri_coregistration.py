"""
Functional and anatomical coregistration
========================================
Standard functional preprocessing and registration of functional image to the
anatomical.
"""

##############################################################################
# Retrieve data
# -------------
from sammba import data_fetchers

retest = data_fetchers.fetch_zurich_test_retest(subjects=[0],
                                                correct_headers=True)

##############################################################################
# retest contains paths to images and data description
anat_filename = retest.anat[0]
func_filename = retest.func[0]
print(func_filename)

##############################################################################
# We use the `Coregistrator`, which coregisters the anatomical to
# a given modality
from sammba.registration import Coregistrator

coregistrator = Coregistrator(output_dir='animal_1366', brain_volume=400,
                              caching=True)
print(coregistrator)

##############################################################################
# `Coregistrator` comes with a parameter `mask_clipping_fraction=.2` which
# sometimes needs to be changed to get a good brain mask. You can check how
# this parameter impacts the brain segmentation
from sammba.segmentation import brain_mask_report

print(brain_extraction_report(anat_filename,
                              clipping_fractions=[.1, .2, .9, None]))

##############################################################################
# Anatomical to functional registration
# -------------------------------------
coregistrator.fit_anat(anat_filename)
coregistrator.fit_modality(func_filename, 'func',
                           slice_timing=False,
                           prior_rigid_body_registration=True)

###############################################################################
# Check out the results
# ---------------------
from nilearn import plotting, image

display = plotting.plot_epi(image.mean_img(coregistrator.undistorted_func_),
                            title='coreg anat edges on top of mean coreg EPI')
display.add_edges(coregistrator.anat_in_func_space_)
plotting.show()