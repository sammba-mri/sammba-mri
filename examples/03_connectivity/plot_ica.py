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
# Load the template and its brain mask
# ------------------------------------
dorr = data_fetchers.fetch_atlas_dorr_2008(downsample='100')
dorr_masks = data_fetchers.fetch_masks_dorr_2008(downsample='100')
print('Path to template is {} and to the brain mask is {}'.format(dorr.t2,
      dorr_masks.brain))

##############################################################################
# Register to the template
# ------------------------
import os
from sammba.registration import TemplateRegistrator

registrator = TemplateRegistrator(brain_volume=400, caching=True,
                                  template=dorr.t2, use_rats_tool=False,
                                  template_brain_mask=dorr_masks.brain)

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
