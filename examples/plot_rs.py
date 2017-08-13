from nipype.interfaces import afni, fsl
from nipype.caching import Memory

memory = Memory('/tmp')
tshift = memory.cache(afni.TShift)
basename = '/tmp/test_func/rs'
out_tshift = tshift(in_file=basename + '.nii.gz',
                    out_file=basename + '_Ts.nii.gz',
#                    tpattern=basename + '_tpattern.txt'
                    )

clip_level = memory.cache(afni.ClipLevel)
out_clip_level = clip_level(in_file=out_tshift.outputs.out_file)

threshold = memory.cache(fsl.Threshold)
out_threshold = threshold(in_file=out_tshift.outputs.out_file,
                          thresh=out_clip_level.outputs.clip_val,
                          out_file=basename + '_TsTm.nii.gz')

volreg = memory.cache(afni.Volreg)
out_volreg = volreg(in_file=out_threshold.outputs.out_file,
                    oned_file=basename + '_TsTmVr.1Dfile.1D',
                    oned_matrix_save=basename + '_TsTmVr.aff12.1D',
                    out_file=basename + '_TsTmVr.nii.gz')

allineate = memory.cache(afni.Allineate)
out_allineate = allineate(in_file=out_tshift.outputs.out_file,
                          master=out_tshift.outputs.out_file,
                          in_matrix=out_volreg.outputs.oned_matrix_save,
                          out_file=basename + '_TsAv.nii.gz')

copy_geom = memory.cache(fsl.CopyGeom)
out_copy_geom = copy_geom(in_file=out_volreg.outputs.out_file,
                          dest_file=out_allineate.outputs.out_file)

tstat = memory.cache(afni.TStat)
out_tstat = tstat(in_file=out_allineate.outputs.out_file,
                  args='-mean',
                  out_file=basename + '_TsAvAv.nii.gz')

#N3BiasFieldCorrection 3 "$p"_TsAvAv.nii.gz "$p"_TsAvAvN3.nii.gz #N4 fails for some reason. Not tried 3dUnifize yet.
#1d_tool.py -infile "$p"_TsTmVr.1Dfile.1D -derivative -censor_prev_TR -collapse_cols euclidean_norm -moderate_mask -1.2 1.2 -show_censor_count -write_censor "$p"_censor.1D -write_CENSORTR "$p"_CENSORTR.1D -overwrite
