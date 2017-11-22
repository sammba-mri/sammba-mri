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

retest = data_fetchers.fetch_zurich_test_retest(subjects=[0, 1, 2],
                                                correct_headers=True)

##############################################################################
# retest contains paths to images and data description
print(retest.anat)

##############################################################################
# Define the writing directory
# ----------------------------
import os

write_dir = os.path.join(os.getcwd(), 'zurich_common')
if not os.path.exists(write_dir):
    os.makedirs(write_dir)

##############################################################################
# Affine registration
# -------------------
from sammba.registration import anats_to_common

affine_register = anats_to_common(retest.anat, write_dir, caching=True)

##############################################################################
# We set caching to True, so that this step computations are not restarted.
# The paths to the registered images can be accessed easilly
registered_anats = affine_register.registered
print(registered_anats)

##############################################################################
# Visalize results
# ----------------
# We plot the edges of one individual anat on top of the average image
from nilearn import plotting, image

average_img = image.mean_img(registered_anats)
display = plotting.plot_anat(average_img, dim=-1.8, title='affine register')
display.add_edges(registered_anats[0])
plotting.show()

##############################################################################
# Visualize pipeline steps
# -------------------------
from sammba.registration import create_pipeline_graph

graph_file = os.path.join(write_dir, 'affine_registration_graph')
create_pipeline_graph('anats_to_common_affine', graph_file)
