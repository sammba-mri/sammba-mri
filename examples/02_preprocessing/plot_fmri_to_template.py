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
# Define the output directory
import os

output_dir = os.path.abspath('zurich_template')

##############################################################################
# We encapsulate them through an object called `FMRISession`
from sammba.registration import FuncSession

animal_session = FuncSession(anat=anat_filename, func=func_filename,
                             animal_id='1366', brain_volume=400, t_r=1.,
                             output_dir=output_dir)

##############################################################################
# Choose the template
# -------------------
# We choose the downsampled version of DORR average, to avoid memory issues
dorr = data_fetchers.fetch_atlas_dorr_2008(downsample='100')

##############################################################################
# dorr contains paths to DORR atlas and average T2 image
head_template_filename = dorr.t2
print(head_template_filename)

##############################################################################
# Prepare the registration
# ------------------------
# First we need to coregister the functional and anatomical to the same space
animal_session.coregister(prior_rigid_body_registration=True)

##############################################################################
# Then we compute a brain template
from nilearn import image

brain_template_filename = os.path.join(output_dir, 'dorr_brain.nii.gz')
brain_template_img = image.math_img(
    'img1 * (img2>0)', img1=head_template_filename, img2=dorr.maps)
brain_template_img.to_filename(brain_template_filename)

##############################################################################
# Register to the template
# ------------------------
# We purposely choose to zero `maxlev` of the nonlinear registration step,
# to speed-up computations
animal_session.register_to_template(
    head_template_filename=head_template_filename,
    brain_template_filename=brain_template_filename,
    func_voxel_size=(.2, .2, .2), maxlev=0)

##############################################################################
# The paths to the registered images are added to `animal_session`
print('Registrered anatomical image is ' + animal_session.registered_anat_)
print('Registrered functional image is ' + animal_session.registered_func_)

##############################################################################
# Visalize results
# ----------------
# We plot the edges of the template on top of the registered anatomical image
from nilearn import plotting

display = plotting.plot_anat(animal_session.registered_anat_, dim=-1.95,
                             title='template edges on top of registered anat')
display.add_edges(head_template_filename)

##############################################################################
# and on top of the average registered functional
from nilearn import image

mean_registered_func_filename = image.mean_img(animal_session.registered_func_)
display = plotting.plot_epi(mean_registered_func_filename,
                            title='template edges on top of registered func')
display.add_edges(head_template_filename)
plotting.show()
