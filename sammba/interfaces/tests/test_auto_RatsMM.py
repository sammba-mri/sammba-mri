from __future__ import unicode_literals
from sammba.interfaces import RatsMM


def test_RatsMM_inputs():
    input_map = dict(args=dict(argstr='%s',
    ),
    environ=dict(nohash=True,
    usedefault=True,
    ),
    ignore_exception=dict(nohash=True,
    usedefault=True,
    ),
    in_file=dict(argstr='%s',
    mandatory=True,
    position=0,
    ),
    intensity_threshold=dict(argstr='-t %s',
    ),
    out_file=dict(argstr='%s',
    position=1,
    ),
    terminal_output=dict(deprecated='1.0.0',
    nohash=True,
    ),
    volume_threshold=dict(argstr='-v %s',
    ),
    )
    inputs = RatsMM.input_spec()

    for key, metadata in list(input_map.items()):
        if key != 'terminal_output':
            for metakey, value in list(metadata.items()):
                assert getattr(inputs.traits()[key], metakey) == value


def test_RatsMM_outputs():
    output_map = dict(out_file=dict(),
    )
    outputs = RatsMM.output_spec()

    for key, metadata in list(output_map.items()):
        for metakey, value in list(metadata.items()):
            assert getattr(outputs.traits()[key], metakey) == value
