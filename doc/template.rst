=================
Template creation
=================

.. contents:: **Contents**
    :local:
    :depth: 1


Small mammals templates 
=======================

The use of standard anatomical templates facilitates data analysis and 
comparison across animals of a given population. Templates can be created by
averaging multiple anatomical scans after a proper alignment. The registration
precision is critical for the template quality.


sammba-MRI strategy
===================

Initially, registration is of extracted brains. Once these are reasonably aligned, whole heads are registered, weighted by masks that, if parameters are chosen well, include some scalp. The amount of scalp is hopefully enough to help in differentiating the brain-scalp boundary without including so much head tissue that it starts to contaminate the registration with the highly variable head tissue.
