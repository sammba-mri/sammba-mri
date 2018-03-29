import os
from sklearn.utils import deprecated
from nilearn._utils.compat import _basestring


@deprecated("Class 'FMRISession' has been refactored and renamed to "
            "'FuncSession'. It will be removed in future release. ")
class FMRISession(object):
    """
    Encapsulation for fMRI data, relative to preprocessing.

    Parameters
    ----------
    func : str
        Path to the functional 4D image

    anat : str
        Path to anatomical image

    animal_id : str
        Animal id
    """

    def __init__(self, func=None, anat=None, animal_id=None):
        self.func = func
        self.anat = anat
        self.animal_id = animal_id

    def _set_items(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

    def _check_inputs(self):
        if self.func is not None:
            if not os.path.isfile(self.func):
                raise IOError('func must be an existing image file,'
                              'you gave {0}'.format(self.func))

        if self.anat is not None:
            if not os.path.isfile(self.anat):
                raise IOError('anat must be an existing image file,'
                              'you gave {0}'.format(self.anat))

        if not isinstance(self.animal_id, _basestring):
            raise ValueError('animal_id must be a string, you provided '
                             '{0}'.format(self.animal_id))

    def _set_output_dir_(self, output_dir):
        setattr(self, 'output_dir_', output_dir)
        if hasattr(self, 'output_dir_'):
            if not os.path.isdir(self.output_dir_):
                os.makedirs(self.output_dir_)
