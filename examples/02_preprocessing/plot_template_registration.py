"""
Registration to template
============================
Registration of anatomical and functional scans to a given template.

"""
##############################################################################
# Retrieve the fMRI data
# ----------------------
from sammba import data_fetchers

retest = data_fetchers.fetch_zurich_test_retest(subjects=[0],
                                                correct_headers=True)

##############################################################################
# retest contains paths to images and data description
anat_filename = retest.anat[0]
func_filename = retest.func[0]
print(func_filename)

##############################################################################
# We encapsulate them through an object called `FMRISession`
from sammba.registration import FMRISession

animal_session = FMRISession(anat=anat_filename, func=func_filename,
                             animal_id='1366')

##############################################################################
# Choose the template
# -------------------
# We choose the downsampled version of DORR average, to avoid memory issues
dorr = data_fetchers.fetch_atlas_dorr_2008(downsample='100')

##############################################################################
# dorr contains paths to DORR atlas and average T2 image
template_filename = dorr.t2
print(template_filename)

##############################################################################
# Define the writing directory
# ----------------------------
import os

write_dir = os.path.join(os.getcwd(), 'zurich_template')
if not os.path.exists(write_dir):
    os.makedirs(write_dir)

##############################################################################
# Perform the registration
# ------------------------
# We purposely choose a low maxlev, to speed-up computations
from sammba.registration import fmri_sessions_to_template

animal_session.register_to_template(t_r=1., brain_volume=400,
                                    head_template_filename=template_filename,
                                    write_dir=write_dir,
                                    slice_timing=True,
                                    prior_rigid_body_registration=True,
                                    maxlev=1, caching=True)

##############################################################################
# The paths to the registered images are added to `animal_session`
print('Registrered anatomical image is ' + animal_session.registered_anat_)
print('Registrered functional image is ' + animal_session.registered_func_)

##############################################################################
# Visalize results
# ----------------
# We plot the edges of the template on top of the registered anatomical image
from nilearn import plotting, image

display = plotting.plot_anat(animal_session.registered_anat_, dim=-1.95,
                             title='template edges on top of registered anat')
display.add_edges(template_filename)

##############################################################################
# and on top of the average registered functional
from nilearn import image

mean_registered_func_filename = image.mean_img(animal_session.registered_func_)
display = plotting.plot_epi(mean_registered_func_filename,
                            title='template edges on top of registered func')
display.add_edges(template_filename)
plotting.show()
