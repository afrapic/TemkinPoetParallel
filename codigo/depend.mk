# DO NOT DELETE THIS LINE - used by make depend
1Dinteg.o: sturmians_base.mod sturmians_data.mod
1Dinteg_f0.o: problemdata_module.mod sturmians_base.mod sturmians_data.mod
2Dinteg.o: size_module.mod sturmians_base.mod sturmians_data.mod
2Dinteg_f0.o: problemdata_module.mod size_module.mod sturmians_base.mod sturmians_data.mod
BasisSet.o: debug_module.mod problemdata_module.mod size_module.mod sturmians_base.mod sturmians_data.mod
baseLag_module.o: complex_module.mod gaulag_module.mod gauleg_module.mod inf_nan_detection.mod pi_module.mod sturmians_base.mod sturmians_data.mod
calcula_base.o: baselag_module.mod sturmians_base.mod sturmians_data.mod
cross_section.o: array_module.mod gauleg_module.mod parablock.mod pi_module.mod problemdata_module.mod size_module.mod sturmians_data.mod
dump_matrix.o: debug_module.mod mapping_module.mod parablock.mod procdata.mod
lin_pzgesv.o: debug_module.mod matrixdat.mod parablock.mod procdata.mod
mapping_module.o: matrixdat.mod procdata.mod
matrices.o: pi_module.mod problemdata_module.mod size_module.mod sturmians_base.mod sturmians_data.mod
problemdata_module.o: debug_module.mod file_module.mod sturmians_base.mod sturmians_data.mod
quadrature.o: gauleg_module.mod problemdata_module.mod sturmians_base.mod sturmians_data.mod
radial_sturm.o: baselag_module.mod sturmians_base.mod sturmians_data.mod
scatt_flux.o: array_module.mod complex_module.mod parablock.mod problemdata_module.mod size_module.mod sturmians_base.mod sturmians_data.mod
scatt_func.o: array_module.mod parablock.mod pi_module.mod problemdata_module.mod size_module.mod sturmians_base.mod sturmians_data.mod
size_module.o: sturmians_base.mod sturmians_data.mod
sturmians_base.o: sturmians_data.mod
sturmians_calc.o: complex_module.mod inf_nan_detection.mod sturmians_base.mod sturmians_data.mod
sturmians_data.o: file_module.mod
sturmians_spline.o: array_module.mod problemdata_module.mod sturmians_base.mod sturmians_data.mod
array_module.mod: .//array_module.o
baselag_module.mod: .//baseLag_module.o
crashblock.mod: .//crashblock.o
debug_module.mod: .//debug_module.o
mapping_module.mod: .//mapping_module.o
matrixdat.mod: .//matrixdat.o
parablock.mod: .//parablock.o
problemdata_module.mod: .//problemdata_module.o
procdata.mod: .//procdata.o
size_module.mod: .//size_module.o
sturmians_base.mod: .//sturmians_base.o
sturmians_calc.mod: .//sturmians_calc.o
sturmians_data.mod: .//sturmians_data.o
