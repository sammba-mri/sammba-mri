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
# We encapsulate them in an object called `FMRISession`
from sammba.registration import FMRISession

animal_session = FMRISession(anat=anat_filename, func=func_filename,
                             animal_id='1366')

##############################################################################
# Define the write directory
# --------------------------
# Note that this directory needs to be empty, to ovoid overwrting conflicts.
import os

write_dir = os.path.join(os.getcwd(), 'zurich_fmri')
if not os.path.exists(write_dir):
    os.makedirs(write_dir)

##############################################################################
# Functional to anatomical registration
# -------------------------------------
from sammba.registration import coregister_fmri_session

coregister_fmri_session(animal_session, 1., write_dir,
                        slice_timing=True,
                        prior_rigid_body_registration=True)

###############################################################################
# Check out the results
# ---------------------
from nilearn import plotting, image

display = plotting.plot_epi(image.mean_img(animal_session.coreg_func_),
                            title='coreg anat edges on top of mean coreg EPI')
display.add_edges(animal_session.coreg_anat_)
plotting.show()
