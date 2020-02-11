================================
List of commands from sammba-MRI
================================

Sammba-MRI is an application programming interface (API). Its modules are
presented underneath.

.. contents:: **List of modules**
   :local:


.. _io_conversions_ref:

:mod:`sammba.io_conversions`: From scanner (DICOM files from Bruker) to NIFTI
=============================================================================

.. automodule:: sammba.io_conversions
   :no-members:
   :no-inherited-members:

**Functions**:

.. currentmodule:: sammba.io_conversions

.. autosummary::
   :toctree: generated/
   :template: function.rst

   dcm_to_nii
   recursive_dcm_to_nii


.. _data_fetchers_ref:

:mod:`sammba.data_fetchers`: To fetch specific datasets and atlas
=================================================================

.. automodule:: sammba.data_fetchers
   :no-members:
   :no-inherited-members:

**Functions**:

.. currentmodule:: sammba.data_fetchers

.. autosummary::
   :toctree: generated/
   :template: function.rst

   fetch_zurich_test_retest
   fetch_zurich_anesthesiant
   fetch_atlas_dorr_2008
   fetch_atlas_waxholm_rat_2014
   fetch_masks_dorr_2008
   fetch_atlas_lemur_mircen_2019
   fetch_lemur_mircen_2019_t2

.. _registration_ref:

:mod:`sammba.registration`: Registration utilities
==================================================

.. automodule:: sammba.registration
   :no-members:
   :no-inherited-members:

**Functions**:

.. currentmodule:: sammba.registration

.. autosummary::
   :toctree: generated/
   :template: function.rst

   anats_to_common
   anats_to_template
   fmri_sessions_to_template
   coregister_fmri_session


**Classes**:

.. currentmodule:: sammba.registration

.. autosummary::
   :toctree: generated/
   :template: class.rst

   FMRISession
   TemplateRegistrator
   Coregistrator

.. _segmentation_ref:

:mod:`sammba.segmentation`: Segmentation utilities
==================================================

.. automodule:: sammba.segmentation
   :no-members:
   :no-inherited-members:

**Functions**:

.. currentmodule:: sammba.segmentation

.. autosummary::
   :toctree: generated/
   :template: function.rst

   brain_extraction_report


.. _graphs_ref:

:mod:`sammba.graphs`: Pipeline graphs
==================================================

.. automodule:: sammba.graphs
   :no-members:
   :no-inherited-members:

**Functions**:

.. currentmodule:: sammba.graphs

.. autosummary::
   :toctree: generated/
   :template: function.rst

   create_pipeline_graph


External tools wrapped in python
================================

.. toctree::
   :maxdepth: 2

   interfaces/generated/sammba.segmentation.interfaces.rst

