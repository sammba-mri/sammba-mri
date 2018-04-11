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

write_dir = os.path.join(os.getcwd(), 'zurich_ica')
write_dir = os.path.join('/tmp', 'zurich_ica')
if not os.path.exists(write_dir):
    os.makedirs(write_dir)

##############################################################################
# Load the template
# -------------------
dorr = data_fetchers.fetch_atlas_dorr_2008(downsample='100')
template_filename = dorr.t2

##############################################################################
# Extract the brain template
# --------------------------
import os
from sammba.registration.base import _apply_mask

# XXX apply_mask_to
dorr_masks = data_fetchers.fetch_masks_dorr_2008(downsample='100')
brain_extracted_template = '/home/bougacha/zurich_ica/Dorr_2008_average_100um_brain.nii.gz'
if not os.path.isfile(brain_extracted_template):
    dorr_masks.brain.to_filename('/tmp/brain_mask_100u.nii.gz')
    brain_extracted_template = _apply_mask(template_filename,
                                 '/tmp/brain_mask_100u.nii.gz',
                                 write_dir)

##############################################################################
# Register to the template
# ------------------------
from sammba.registration import template_registrator

registered_funcs = []
for anat, func in zip(retest.anat, retest.func):
    output_dir = os.path.join(write_dir,
                              os.path.basename(os.path.dirname(anat)))
    registrator = template_registrator.TemplateRegistrator(
        anat=anat,
        template=template_filename,
        brain_extracted_template=brain_extracted_template,
        caching=False, output_dir=output_dir,
        brain_volume=400)
    registrator.segment()
    registrator.normalize()
    normalized_func = registrator.normalize_modality(func, 'func', t_r=1.,
                                                     voxel_size=(.3, .3, .3),
                                                 prior_rigid_body_registration=True)
    registered_funcs.append(normalized_func)

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
                         anat_img='/home/Promane/2014-ROMANE/5_Experimental-Plan-Experiments-Results/mouse/XNAT_CSD_MRI_MOUSE/20170925_goodposttemplate/Qw4_meanhead200.nii.gz',
                         title='All ICA components')
