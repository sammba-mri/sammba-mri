"""
fMRI registration
=================

Here we show how to create a template from multiple anatomical scans.
Initially, registration is of extracted brains. Once these are reasonably
aligned, whole heads are registered, weighted by masks that, if parameters
are chosen well, include some scalp. The amount of scalp is hopefully enough
to help in differentiating the brain-scalp boundary without including so much
head tissue that it starts to contaminate the registration with the highly
variable head tissue.
"""


##############################################################################
# Retrieve data
# -------------
from sammba import data_fetchers

retest = data_fetchers.fetch_zurich_test_retest(subjects=[0],
                                                correct_headers=True)

##############################################################################
# Define the write directory
# --------------------------
# Note that this directory needs to be empty, to ovoid overwrting conflicts.
import os

write_dir = os.path.join(os.getcwd(), 'results_1366')
if not os.path.exists(write_dir):
    os.makedirs(write_dir)

##############################################################################
# Functional to anatomical registration
# -------------------------------------
from sammba import preprocessors

registered_func_filename, resampled_anat_filename = \
    preprocessors.register_func_to_anat(retest.func[0],
                                        retest.anat[0],
                                        tr='1',
                                        write_dir=write_dir)

###############################################################################
# Check out the results
# ---------------------
from nilearn import plotting, image

display = plotting.plot_epi(image.mean_img(registered_func_filename),
                            title='anat edges on top of registered mean EPI')
display.add_edges(resampled_anat_filename)
plotting.show()
