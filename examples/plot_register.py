"""
Registration to common space
============================

Here we show how to create a template from multiple anatomical scans and
register all of them to it.
Initially, registration is of extracted brains. Once these are reasonably
aligned, whole heads are registered, weighted by masks that, if parameters
are chosen well, include some scalp. The amount of scalp is hopefully enough
to help in differentiating the brain-scalp boundary without including so much
head tissue that it starts to contaminate the registration with the highly
variable head tissue.
"""
##############################################################################
# Retrieve the data
# -----------------
from sammba import data_fetchers

retest = data_fetchers.fetch_zurich_test_retest(subjects=range(6))

##############################################################################
# Define the writing directory
# ----------------------------
import os

write_dir = os.path.join(os.getcwd(), 'zurich_registration')
if not os.path.exists(write_dir):
    os.makedirs(write_dir)

##############################################################################
# Correct orientation
# -------------------
# This data has wrong sforms and qforms in the headers, so we correct them.
from sammba.preprocessors import correct_affines

anat_filenames = []
for raw_anat_finame in retest.anat:
    mouse_id = os.path.basename(os.path.dirname(raw_anat_finame))
    anat_filename = os.path.join(write_dir, mouse_id + '_anat.nii.gz')
    correct_affines(in_file=raw_anat_finame,
                    axes_to_permute=[(1, 2)],
                    axes_to_flip=[0],
                    out_file=anat_filename)
    anat_filenames.append(anat_filename)

##############################################################################
# Register: Rigid-body
# --------------------
from sammba.preprocessors import register_to_common

rigid_body = register_to_common(anat_filenames, write_dir,
                                registration_kind='rigid-body',
                                caching=True)
stop
##############################################################################
# We set caching to True, so that this step computations are not restarted.
##############################################################################
# Visalize results
# ----------------
# We plot the edges of one individual anat on top of the average image
from nilearn import plotting, image

average_img = image.mean_img(rigid_body.registered_anats)
display = plotting.plot_anat(average_img, dim=-1.7, title='rigid-body')
display.add_edges(rigid_body.registered_anats[0])


##############################################################################
# Register: Affine
# ----------------
# Besides translation and rotation, we allow zoom and shear in the transforms.
affine = register_to_common(anat_filenames, write_dir,
                            registration_kind='affine',
                            caching=True)

average_img = image.mean_img(affine.registered_anats)
display = plotting.plot_anat(average_img, dim=-1.7, title='affine')
display.add_edges(affine.registered_anats[0])

##############################################################################
# Register: Nonlinear 
# -------------------
# We also allow nonlinear warps.
nonlinear = register_to_common(anat_filenames, write_dir,
                               registration_kind='nonlinear',
                               nonlinear_iters=2,
                               caching=True)

average_img = image.mean_img(nonlinear.registered_anats)
display = plotting.plot_anat(average_img, dim=-1.6, title='nonlinear')
display.add_edges(nonlinear.registered_anats[0])
plotting.show()
