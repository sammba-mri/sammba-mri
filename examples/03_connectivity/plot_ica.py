"""
ICA on mouse
============
Independent components analysis on 5 mice.

"""
##############################################################################
# Retrieve the fMRI data
# ----------------------
from sammba import data_fetchers

retest = data_fetchers.fetch_zurich_test_retest(subjects=range(5),
                                                correct_headers=True)

##############################################################################
# Define the writing directory
# ----------------------------
import os

output_dir = os.path.abspath('zurich_ica')

##############################################################################
# Link data per animal
# --------------------
# We encapsulate data for each mouse through an object called `FuncSession`
from sammba.registration import FuncSession

sessions = []
for (animal_id, anat_file, func_file) in zip(['0', '1', '2', '3', '4'],
                                             retest.anat, retest.func):
    session = FuncSession(anat=anat_file, func=func_file,
                          animal_id=animal_id, t_r=1., brain_volume=400,
                          output_dir=output_dir)
    sessions.append(session)

##############################################################################
# Load the head template
# ----------------------
dorr = data_fetchers.fetch_atlas_dorr_2008(downsample='100')
head_template_filename = dorr.t2

##############################################################################
# Compute the brain template
# --------------------------
import os
from nilearn import image

brain_template_filename = os.path.join(output_dir, 'dorr_brain.nii.gz')
brain_template_img = image.math_img(
    'img1 * (img2>0)', img1=head_template_filename, img2=dorr.maps)
brain_template_img.to_filename(brain_template_filename)

##############################################################################
# Coregister anat to func
# -----------------------
for session in sessions:
    session.coregister(prior_rigid_body_registration=True, caching=True)

##############################################################################
# Register to the template
# ------------------------
for session in sessions:
    session.register_to_template(
        head_template_filename=head_template_filename,
        brain_template_filename=brain_template_filename,
        func_voxel_size=(.2, .2, .2), maxlev=0, caching=True)

registered_funcs = [session.registered_func_ for session in sessions]

##############################################################################
# Run ICA
# -------
from nilearn.decomposition import CanICA

canica = CanICA(n_components=30, smoothing_fwhm=.3, n_jobs=-1)
canica.fit(registered_funcs)

##############################################################################
# Retrieve the independent components in brain space.
components_img = canica.masker_.inverse_transform(canica.components_)

##############################################################################
# Visualize the components
# ------------------------
# We can plot the outline of all components on one figure.
from nilearn import plotting

plotting.plot_prob_atlas(components_img,
                         anat_img=head_template_filename,
                         title='All ICA components')
