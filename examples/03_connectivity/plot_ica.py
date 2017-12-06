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
# We encapsulate data for each mouse through an object called `FMRISession`
from sammba.registration import FMRISession

sessions = [FMRISession(anat=anat_filename, func=func_filename)
            for anat_filename, func_filename in zip(retest.anat, retest.func)]

##############################################################################
# Load the template
# -------------------
dorr = data_fetchers.fetch_atlas_dorr_2008(downsample='100')
template_filename = dorr.t2

##############################################################################
# Define the writing directory
# ----------------------------
import os

write_dir = os.path.join(os.getcwd(), 'zurich_ica')
if not os.path.exists(write_dir):
    os.makedirs(write_dir)

##############################################################################
# Register to the template
# ------------------------
from sammba.registration import fmri_sessions_to_template

fmri_sessions_to_template(sessions, head_template_filename=dorr.t2,
                          t_r=1., write_dir=write_dir,
                          brain_volume=400,
                          prior_rigid_body_registration=True, caching=True)

registered_funcs = [sess.registered_func_ for sess in sessions]

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
                         anat_img=dorr.t2,
                         title='All ICA components')
