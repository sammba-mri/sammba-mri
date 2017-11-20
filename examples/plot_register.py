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

retest = data_fetchers.fetch_zurich_test_retest(subjects=range(3),
                                                correct_headers=True)

##############################################################################
# Define the writing directory
# ----------------------------
import os

write_dir = os.path.join(os.getcwd(), 'zurich')
if not os.path.exists(write_dir):
    os.makedirs(write_dir)

##############################################################################
# Register: Affine registration
# -----------------------------
from sammba.preprocessors import register_to_common

affine = register_to_common(retest.anat, write_dir, caching=True, verbose=1)

##############################################################################
# We set caching to True, so that this step computations are not restarted.

##############################################################################
# Visalize results
# ----------------
# We plot the edges of one individual anat on top of the average image
from nilearn import plotting, image

average_img = image.mean_img(affine.registered)
display = plotting.plot_anat(average_img, dim=-1.7, title='affine')
display.add_edges(affine.affine_anats[0])


##############################################################################
# Register: Rigid-body
# --------------------
# Compare to rigid-body registration
rigid = register_to_common(retest.anat, write_dir, registration_kind='rigid')

average_img = image.mean_img(rigid.registered_anats)
display = plotting.plot_anat(average_img, dim=-1.6, title='rigid-body')
display.add_edges(rigid.registered_anats[0])
plotting.show()
