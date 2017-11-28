import os


class Mammal(object):
    """
    Encapsulation for subject data, relative to preprocessing.

    Parameters
    ----------
    func : str
        Path to the functional 4D image

    anat : str
        Path to anatomical image

    mammal_id : str, optional
        Mammal id
    """

    def __init__(self, func=None, anat=None, mammal_id="mam001"):
        self.func = func
        self.anat = anat
        self.mammal_id = mammal_id

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