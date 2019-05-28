MIRCEN 2019 atlas


Notes
-----
MRI template and 80-region atlas for the mouse lemur primate Microcebus
murinus. Generated from 34 animals aged 15-60 months old scanned at 7T using
a T2-weighted sequence, resolution 115 × 115 × 230 µm3.
Aided by brain identification using the IIBI RATS skull-stripper,
images were registered and averaged with AFNI through linear then non-linear
stages to produce a final template. This was upsampled to 91 µm3 isotropic for
hand-segmentation of structures.

At the time of writing (September 2017) this an imperfect draft release.
The template is not in an optimal orientation, and the atlas is still rough.
Improved versions will be uploaded over the next few months.

Content
-------
    :"t2": str. path to mnc file containing the T2 weighted average.
    :"maps": str. path to mnc file containing regions.
    :"labels": numpy recarray containing the ids and names of each region
    :"description": description about the atlas.

Licence
-------
CeCill v2