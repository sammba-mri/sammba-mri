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
# Load the template
# -------------------
dorr = data_fetchers.fetch_atlas_dorr_2008(downsample='100')
template = dorr.t2

##############################################################################
# Load the brain mask of the template
# -----------------------------------
import os

dorr_masks = data_fetchers.fetch_masks_dorr_2008(downsample='100')
template_brain_mask = os.path.join('ica', 'dorr_100_mask.nii.gz')
if not os.path.isfile(template_brain_mask):
    dorr_masks.brain.to_filename(template_brain_mask)

##############################################################################
# Register to the template
# ------------------------
from sammba.registration import TemplateRegistrator

registrator = TemplateRegistrator(brain_volume=400, caching=True,
                                  template=template,
                                  template_brain_mask=template_brain_mask)

registered_funcs = []
for anat, func in zip(retest.anat, retest.func):
    animal_id = os.path.basename(os.path.dirname(anat))
    registrator.output_dir = os.path.join('ica', animal_id)
    registrator.fit_anat(anat)
    registrator.fit_modality(func, 'func', t_r=1., voxel_size=(.3, .3, .3),
                             prior_rigid_body_registration=True)
    registered_funcs.append(registrator.registered_func_)

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
                         bg_img=registrator.template_brain_,
                         display_mode='z',
                         title='All ICA components')
