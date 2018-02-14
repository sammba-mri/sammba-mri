"""
Functional and anatomical coregistration
========================================

Standard functional preprocessing and registration of functional image to the
anatomical.

"""

##############################################################################
# Retrieve data
# -------------
from sammba import data_fetchers

retest = data_fetchers.fetch_zurich_test_retest(subjects=[0],
                                                correct_headers=True)

##############################################################################
# retest contains paths to images and data description
anat_filename = retest.anat[0]
func_filename = retest.func[0]
print(func_filename)

##############################################################################
# We encapsulate them in an object called `FMRISession`
from sammba.registration import FuncSession

animal_session = FuncSession(anat=anat_filename, func=func_filename,
                             animal_id='1366', brain_volume=400, t_r=1.,
                             output_dir='zurich_fmri')

##############################################################################
# Functional to anatomical registration
# -------------------------------------
animal_session.coregister(prior_rigid_body_registration=True)

###############################################################################
# Check out the results
# ---------------------
from nilearn import plotting, image

display = plotting.plot_epi(image.mean_img(animal_session.coreg_func_),
                            title='coreg anat edges on top of mean coreg EPI')
display.add_edges(animal_session.coreg_anat_)
plotting.show()
