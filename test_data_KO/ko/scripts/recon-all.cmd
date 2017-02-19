\n\n#---------------------------------
# New invocation of recon-all Wed May  4 14:27:33 PDT 2016 
\n mri_convert /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/raw/ko/T1w.nii.gz /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/orig/001.mgz \n
#--------------------------------------------
#@# T2/FLAIR Input Wed May  4 14:27:48 PDT 2016
\n mri_convert --no_scale 1 /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/raw/ko/T2w.nii.gz /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/orig/T2raw.mgz \n
#--------------------------------------------
#@# MotionCor Wed May  4 14:28:25 PDT 2016
\n cp /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/orig/001.mgz /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/rawavg.mgz \n
\n mri_convert /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/rawavg.mgz /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/orig.mgz --conform \n
\n mri_add_xform_to_header -c /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/transforms/talairach.xfm /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/orig.mgz /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/orig.mgz \n
#--------------------------------------------
#@# Talairach Wed May  4 14:28:40 PDT 2016
\n mri_nu_correct.mni --no-rescale --i orig.mgz --o orig_nu.mgz --n 1 --proto-iters 1000 --distance 50 \n
\n talairach_avi --i orig_nu.mgz --xfm transforms/talairach.auto.xfm \n
talairach_avi log file is transforms/talairach_avi.log...
\n cp transforms/talairach.auto.xfm transforms/talairach.xfm \n
#--------------------------------------------
#@# Talairach Failure Detection Wed May  4 14:36:40 PDT 2016
\n talairach_afd -T 0.005 -xfm transforms/talairach.xfm \n
\n awk -f /Applications/freesurfer/bin/extract_talairach_avi_QA.awk /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/transforms/talairach_avi.log \n
\n tal_QC_AZS /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/transforms/talairach_avi.log \n
#--------------------------------------------
#@# Nu Intensity Correction Wed May  4 14:36:42 PDT 2016
\n mri_nu_correct.mni --i orig.mgz --o nu.mgz --uchar transforms/talairach.xfm --n 2 \n
\n mri_add_xform_to_header -c /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/transforms/talairach.xfm nu.mgz nu.mgz \n
#--------------------------------------------
#@# Intensity Normalization Wed May  4 14:43:42 PDT 2016
\n mri_normalize -g 1 -mprage nu.mgz T1.mgz \n
#--------------------------------------------
#@# Skull Stripping Wed May  4 14:45:43 PDT 2016
\n mri_em_register -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mri_em_register.skull.dat -skull nu.mgz /Applications/freesurfer/average/RB_all_withskull_2016-03-21.gca transforms/talairach_with_skull.lta \n
\n mri_watershed -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mri_watershed.dat -T1 -brain_atlas /Applications/freesurfer/average/RB_all_withskull_2016-03-21.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz \n
\n cp brainmask.auto.mgz brainmask.mgz \n
\n\n#---------------------------------
# New invocation of recon-all Wed May  4 15:22:25 PDT 2016 
#--------------------------------------------
#@# Sphere lh Wed May  4 15:22:30 PDT 2016
\n mris_sphere -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mris_sphere.lh.dat -seed 1234 ../surf/lh.inflated ../surf/lh.sphere \n
\n\n#---------------------------------
# New invocation of recon-all Wed May  4 15:23:02 PDT 2016 
#-------------------------------------
#@# EM Registration Wed May  4 15:23:05 PDT 2016
\n mri_em_register -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mri_em_register.dat -uns 3 -mask brainmask.mgz nu.mgz /Applications/freesurfer/average/RB_all_2016-03-21.gca transforms/talairach.lta \n
#--------------------------------------
#@# CA Normalize Wed May  4 15:33:04 PDT 2016
\n mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /Applications/freesurfer/average/RB_all_2016-03-21.gca transforms/talairach.lta norm.mgz \n
#--------------------------------------
#@# CA Reg Wed May  4 15:34:22 PDT 2016
\n mri_ca_register -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mri_ca_register.dat -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /Applications/freesurfer/average/RB_all_2016-03-21.gca transforms/talairach.m3z \n
#--------------------------------------
#@# Remove Neck Wed May  4 17:32:29 PDT 2016
\n mri_remove_neck -radius 25 nu.mgz transforms/talairach.m3z /Applications/freesurfer/average/RB_all_2016-03-21.gca nu_noneck.mgz \n
#--------------------------------------
#@# SubCort Seg Wed May  4 17:33:33 PDT 2016
\n mri_ca_label -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mri_ca_label.dat -relabel_unlikely 9 .3 -prior 0.5 -align norm.mgz transforms/talairach.m3z /Applications/freesurfer/average/RB_all_2016-03-21.gca aseg.auto_noCCseg.mgz \n
\n mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/transforms/cc_up.lta ko \n
#--------------------------------------
#@# Merge ASeg Wed May  4 18:03:44 PDT 2016
\n cp aseg.auto.mgz aseg.presurf.mgz \n
#--------------------------------------------
#@# Intensity Normalization2 Wed May  4 18:03:44 PDT 2016
\n mri_normalize -mprage -aseg aseg.presurf.mgz -mask brainmask.mgz norm.mgz brain.mgz \n
#--------------------------------------------
#@# Mask BFS Wed May  4 18:06:59 PDT 2016
\n mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz \n
#--------------------------------------------
#@# WM Segmentation Wed May  4 18:07:02 PDT 2016
\n mri_segment -mprage brain.mgz wm.seg.mgz \n
\n mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.presurf.mgz wm.asegedit.mgz \n
\n mri_pretess wm.asegedit.mgz wm norm.mgz wm.mgz \n
#--------------------------------------------
#@# Fill Wed May  4 18:08:46 PDT 2016
\n mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz \n
#--------------------------------------------
#@# Tessellate lh Wed May  4 18:09:21 PDT 2016
\n mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz \n
\n mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix \n
\n rm -f ../mri/filled-pretess255.mgz \n
\n mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix \n
#--------------------------------------------
#@# Tessellate rh Wed May  4 18:09:34 PDT 2016
\n mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz \n
\n mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix \n
\n rm -f ../mri/filled-pretess127.mgz \n
\n mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix \n
#--------------------------------------------
#@# Smooth1 lh Wed May  4 18:09:47 PDT 2016
\n mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix \n
#--------------------------------------------
#@# Smooth1 rh Wed May  4 18:09:58 PDT 2016
\n mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix \n
#--------------------------------------------
#@# Inflation1 lh Wed May  4 18:10:09 PDT 2016
\n mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix \n
#--------------------------------------------
#@# Inflation1 rh Wed May  4 18:10:39 PDT 2016
\n mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix \n
#--------------------------------------------
#@# QSphere lh Wed May  4 18:11:10 PDT 2016
\n mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix \n
#--------------------------------------------
#@# QSphere rh Wed May  4 18:13:57 PDT 2016
\n mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix \n
#--------------------------------------------
#@# Fix Topology Copy lh Wed May  4 18:16:43 PDT 2016
\n cp ../surf/lh.orig.nofix ../surf/lh.orig \n
\n cp ../surf/lh.inflated.nofix ../surf/lh.inflated \n
#--------------------------------------------
#@# Fix Topology Copy rh Wed May  4 18:16:44 PDT 2016
\n cp ../surf/rh.orig.nofix ../surf/rh.orig \n
\n cp ../surf/rh.inflated.nofix ../surf/rh.inflated \n
#@# Fix Topology lh Wed May  4 18:16:45 PDT 2016
\n mris_fix_topology -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mris_fix_topology.lh.dat -mgz -sphere qsphere.nofix -ga -seed 1234 ko lh \n
#@# Fix Topology rh Wed May  4 18:34:18 PDT 2016
\n mris_fix_topology -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mris_fix_topology.rh.dat -mgz -sphere qsphere.nofix -ga -seed 1234 ko rh \n
\n mris_euler_number ../surf/lh.orig \n
\n mris_euler_number ../surf/rh.orig \n
\n mris_remove_intersection ../surf/lh.orig ../surf/lh.orig \n
\n rm ../surf/lh.inflated \n
\n mris_remove_intersection ../surf/rh.orig ../surf/rh.orig \n
\n rm ../surf/rh.inflated \n
#--------------------------------------------
#@# Make White Surf lh Wed May  4 18:50:55 PDT 2016
\n mris_make_surfaces -aseg ../mri/aseg.presurf -noaparc -whiteonly -mgz -T1 brain.finalsurfs ko lh \n
#--------------------------------------------
#@# Make White Surf rh Wed May  4 18:55:29 PDT 2016
\n mris_make_surfaces -aseg ../mri/aseg.presurf -noaparc -whiteonly -mgz -T1 brain.finalsurfs ko rh \n
#--------------------------------------------
#@# Smooth2 lh Wed May  4 18:59:58 PDT 2016
\n mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white ../surf/lh.smoothwm \n
#--------------------------------------------
#@# Smooth2 rh Wed May  4 19:00:08 PDT 2016
\n mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white ../surf/rh.smoothwm \n
#--------------------------------------------
#@# Inflation2 lh Wed May  4 19:00:18 PDT 2016
\n mris_inflate -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mris_inflate.lh.dat ../surf/lh.smoothwm ../surf/lh.inflated \n
#--------------------------------------------
#@# Inflation2 rh Wed May  4 19:00:48 PDT 2016
\n mris_inflate -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mris_inflate.rh.dat ../surf/rh.smoothwm ../surf/rh.inflated \n
#--------------------------------------------
#@# Curv .H and .K lh Wed May  4 19:01:18 PDT 2016
\n mris_curvature -w lh.white \n
\n mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 lh.inflated \n
#--------------------------------------------
#@# Curv .H and .K rh Wed May  4 19:02:15 PDT 2016
\n mris_curvature -w rh.white \n
\n mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 rh.inflated \n
\n#-----------------------------------------
#@# Curvature Stats lh Wed May  4 19:03:13 PDT 2016
\n mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm ko lh curv sulc \n
\n#-----------------------------------------
#@# Curvature Stats rh Wed May  4 19:03:22 PDT 2016
\n mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm ko rh curv sulc \n
#--------------------------------------------
#@# Sphere lh Wed May  4 19:03:31 PDT 2016
\n mris_sphere -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mris_sphere.lh.dat -seed 1234 ../surf/lh.inflated ../surf/lh.sphere \n
#--------------------------------------------
#@# Sphere rh Wed May  4 19:46:23 PDT 2016
\n mris_sphere -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mris_sphere.rh.dat -seed 1234 ../surf/rh.inflated ../surf/rh.sphere \n
#--------------------------------------------
#@# Surf Reg lh Wed May  4 20:05:22 PDT 2016
\n mris_register -curv -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mris_register.lh.dat ../surf/lh.sphere /Applications/freesurfer/average/lh.curvature.buckner40.2016-03-20.tif ../surf/lh.sphere.reg \n
#--------------------------------------------
#@# Surf Reg rh Wed May  4 20:58:13 PDT 2016
\n mris_register -curv -rusage /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/touch/rusage.mris_register.rh.dat ../surf/rh.sphere /Applications/freesurfer/average/rh.curvature.buckner40.2016-03-20.tif ../surf/rh.sphere.reg \n
#--------------------------------------------
#@# Jacobian white lh Wed May  4 21:38:05 PDT 2016
\n mris_jacobian ../surf/lh.white ../surf/lh.sphere.reg ../surf/lh.jacobian_white \n
#--------------------------------------------
#@# Jacobian white rh Wed May  4 21:38:08 PDT 2016
\n mris_jacobian ../surf/rh.white ../surf/rh.sphere.reg ../surf/rh.jacobian_white \n
#--------------------------------------------
#@# AvgCurv lh Wed May  4 21:38:10 PDT 2016
\n mrisp_paint -a 5 /Applications/freesurfer/average/lh.curvature.buckner40.2016-03-20.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv \n
#--------------------------------------------
#@# AvgCurv rh Wed May  4 21:38:13 PDT 2016
\n mrisp_paint -a 5 /Applications/freesurfer/average/rh.curvature.buckner40.2016-03-20.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv \n
#-----------------------------------------
#@# Cortical Parc lh Wed May  4 21:38:15 PDT 2016
\n mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 ko lh ../surf/lh.sphere.reg /Applications/freesurfer/average/lh.DKatlas.2016-03-20.gcs ../label/lh.aparc.annot \n
#-----------------------------------------
#@# Cortical Parc rh Wed May  4 21:38:31 PDT 2016
\n mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 ko rh ../surf/rh.sphere.reg /Applications/freesurfer/average/rh.DKatlas.2016-03-20.gcs ../label/rh.aparc.annot \n
#--------------------------------------------
#@# Make Pial Surf lh Wed May  4 21:38:47 PDT 2016
\n mris_make_surfaces -orig_white white -orig_pial white -aseg ../mri/aseg.presurf -nowhite -mgz -T1 brain.finalsurfs ko lh \n
#--------------------------------------------
#@# Make Pial Surf rh Wed May  4 21:46:52 PDT 2016
\n mris_make_surfaces -orig_white white -orig_pial white -aseg ../mri/aseg.presurf -nowhite -mgz -T1 brain.finalsurfs ko rh \n
#--------------------------------------------
#@# Refine Pial Surfs w/ T2/FLAIR Wed May  4 21:54:59 PDT 2016
\n bbregister --s ko --mov /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/orig/T2raw.mgz --lta /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/transforms/T2raw.auto.lta --init-coreg --T2 \n
\n cp /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/transforms/T2raw.auto.lta /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/transforms/T2raw.lta \n
\n mri_convert -odt float -at /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/transforms/T2raw.lta -rl /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/orig.mgz /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/orig/T2raw.mgz /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/T2.prenorm.mgz \n
\n mri_normalize -sigma 0.5 -nonmax_suppress 0 -min_dist 1 -aseg /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/aseg.presurf.mgz -surface /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/surf/rh.white identity.nofile -surface /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/surf/lh.white identity.nofile /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/T2.prenorm.mgz /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/T2.norm.mgz \n
\n mri_mask /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/T2.norm.mgz /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/brainmask.mgz /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/mri/T2.mgz \n
\n cp -v /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/surf/lh.pial /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/surf/lh.woT2.pial \n
\n mris_make_surfaces -orig_white white -orig_pial woT2.pial -aseg ../mri/aseg.presurf -nowhite -mgz -T1 brain.finalsurfs -T2 ../mri/T2 -nsigma_above 2 -nsigma_below 5 ko lh \n
\n cp -v /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/surf/rh.pial /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/surf/rh.woT2.pial \n
\n mris_make_surfaces -orig_white white -orig_pial woT2.pial -aseg ../mri/aseg.presurf -nowhite -mgz -T1 brain.finalsurfs -T2 ../mri/T2 -nsigma_above 2 -nsigma_below 5 ko rh \n
#--------------------------------------------
#@# Surf Volume lh Wed May  4 23:06:31 PDT 2016
#--------------------------------------------
#@# Surf Volume rh Wed May  4 23:06:40 PDT 2016
#--------------------------------------------
#@# Cortical ribbon mask Wed May  4 23:06:48 PDT 2016
\n mris_volmask --aseg_name aseg.presurf --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon ko \n
#-----------------------------------------
#@# Parcellation Stats lh Wed May  4 23:18:53 PDT 2016
\n mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab ko lh white \n
\n mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.pial.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab ko lh pial \n
#-----------------------------------------
#@# Parcellation Stats rh Wed May  4 23:20:09 PDT 2016
\n mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab ko rh white \n
\n mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.pial.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab ko rh pial \n
#-----------------------------------------
#@# Cortical Parc 2 lh Wed May  4 23:21:26 PDT 2016
\n mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 ko lh ../surf/lh.sphere.reg /Applications/freesurfer/average/lh.CDatlas.2016-03-20.gcs ../label/lh.aparc.a2009s.annot \n
#-----------------------------------------
#@# Cortical Parc 2 rh Wed May  4 23:21:49 PDT 2016
\n mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 ko rh ../surf/rh.sphere.reg /Applications/freesurfer/average/rh.CDatlas.2016-03-20.gcs ../label/rh.aparc.a2009s.annot \n
#-----------------------------------------
#@# Parcellation Stats 2 lh Wed May  4 23:22:11 PDT 2016
\n mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab ko lh white \n
#-----------------------------------------
#@# Parcellation Stats 2 rh Wed May  4 23:22:53 PDT 2016
\n mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab ko rh white \n
#-----------------------------------------
#@# Cortical Parc 3 lh Wed May  4 23:23:38 PDT 2016
\n mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 ko lh ../surf/lh.sphere.reg /Applications/freesurfer/average/lh.DKTatlas.2016-03-20.gcs ../label/lh.aparc.DKTatlas.annot \n
#-----------------------------------------
#@# Cortical Parc 3 rh Wed May  4 23:23:56 PDT 2016
\n mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 ko rh ../surf/rh.sphere.reg /Applications/freesurfer/average/rh.DKTatlas.2016-03-20.gcs ../label/rh.aparc.DKTatlas.annot \n
#-----------------------------------------
#@# Parcellation Stats 3 lh Wed May  4 23:24:15 PDT 2016
\n mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.DKTatlas.stats -b -a ../label/lh.aparc.DKTatlas.annot -c ../label/aparc.annot.DKTatlas.ctab ko lh white \n
#-----------------------------------------
#@# Parcellation Stats 3 rh Wed May  4 23:24:54 PDT 2016
\n mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.DKTatlas.stats -b -a ../label/rh.aparc.DKTatlas.annot -c ../label/aparc.annot.DKTatlas.ctab ko rh white \n
#-----------------------------------------
#@# WM/GM Contrast lh Wed May  4 23:25:30 PDT 2016
\n pctsurfcon --s ko --lh-only \n
#-----------------------------------------
#@# WM/GM Contrast rh Wed May  4 23:25:48 PDT 2016
\n pctsurfcon --s ko --rh-only \n
#-----------------------------------------
#@# Relabel Hypointensities Wed May  4 23:26:05 PDT 2016
\n mri_relabel_hypointensities aseg.presurf.mgz ../surf aseg.presurf.hypos.mgz \n
#-----------------------------------------
#@# AParc-to-ASeg aparc Wed May  4 23:26:25 PDT 2016
\n mri_aparc2aseg --s ko --volmask --aseg aseg.presurf.hypos \n
#-----------------------------------------
#@# AParc-to-ASeg a2009s Wed May  4 23:28:06 PDT 2016
\n mri_aparc2aseg --s ko --volmask --aseg aseg.presurf.hypos --annot aparc.a2009s \n
#-----------------------------------------
#@# AParc-to-ASeg DKTatlas Wed May  4 23:29:52 PDT 2016
\n mri_aparc2aseg --s ko --volmask --aseg aseg.presurf.hypos --annot aparc.DKTatlas \n
#-----------------------------------------
#@# APas-to-ASeg Wed May  4 23:31:46 PDT 2016
\n apas2aseg --i aparc+aseg.mgz --o aseg.mgz \n
#--------------------------------------------
#@# ASeg Stats Wed May  4 23:31:53 PDT 2016
\n mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --empty --brainmask mri/brainmask.mgz --brain-vol-from-seg --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol --totalgray --euler --ctab /Applications/freesurfer/ASegStatsLUT.txt --subject ko \n
#-----------------------------------------
#@# WMParc Wed May  4 23:35:07 PDT 2016
\n mri_aparc2aseg --s ko --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz \n
\n mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --subject ko --surf-wm-vol --ctab /Applications/freesurfer/WMParcStatsLUT.txt --etiv \n
#--------------------------------------------
#@# BA_exvivo Labels lh Wed May  4 23:42:50 PDT 2016
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA1_exvivo.label --trgsubject ko --trglabel ./lh.BA1_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA2_exvivo.label --trgsubject ko --trglabel ./lh.BA2_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA3a_exvivo.label --trgsubject ko --trglabel ./lh.BA3a_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA3b_exvivo.label --trgsubject ko --trglabel ./lh.BA3b_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA4a_exvivo.label --trgsubject ko --trglabel ./lh.BA4a_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA4p_exvivo.label --trgsubject ko --trglabel ./lh.BA4p_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA6_exvivo.label --trgsubject ko --trglabel ./lh.BA6_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA44_exvivo.label --trgsubject ko --trglabel ./lh.BA44_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA45_exvivo.label --trgsubject ko --trglabel ./lh.BA45_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.V1_exvivo.label --trgsubject ko --trglabel ./lh.V1_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.V2_exvivo.label --trgsubject ko --trglabel ./lh.V2_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.MT_exvivo.label --trgsubject ko --trglabel ./lh.MT_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.entorhinal_exvivo.label --trgsubject ko --trglabel ./lh.entorhinal_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.perirhinal_exvivo.label --trgsubject ko --trglabel ./lh.perirhinal_exvivo.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA1_exvivo.thresh.label --trgsubject ko --trglabel ./lh.BA1_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA2_exvivo.thresh.label --trgsubject ko --trglabel ./lh.BA2_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA3a_exvivo.thresh.label --trgsubject ko --trglabel ./lh.BA3a_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA3b_exvivo.thresh.label --trgsubject ko --trglabel ./lh.BA3b_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA4a_exvivo.thresh.label --trgsubject ko --trglabel ./lh.BA4a_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA4p_exvivo.thresh.label --trgsubject ko --trglabel ./lh.BA4p_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA6_exvivo.thresh.label --trgsubject ko --trglabel ./lh.BA6_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA44_exvivo.thresh.label --trgsubject ko --trglabel ./lh.BA44_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.BA45_exvivo.thresh.label --trgsubject ko --trglabel ./lh.BA45_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.V1_exvivo.thresh.label --trgsubject ko --trglabel ./lh.V1_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.V2_exvivo.thresh.label --trgsubject ko --trglabel ./lh.V2_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.MT_exvivo.thresh.label --trgsubject ko --trglabel ./lh.MT_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.entorhinal_exvivo.thresh.label --trgsubject ko --trglabel ./lh.entorhinal_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/lh.perirhinal_exvivo.thresh.label --trgsubject ko --trglabel ./lh.perirhinal_exvivo.thresh.label --hemi lh --regmethod surface \n
\n mris_label2annot --s ko --hemi lh --ctab /Applications/freesurfer/average/colortable_BA.txt --l lh.BA1_exvivo.label --l lh.BA2_exvivo.label --l lh.BA3a_exvivo.label --l lh.BA3b_exvivo.label --l lh.BA4a_exvivo.label --l lh.BA4p_exvivo.label --l lh.BA6_exvivo.label --l lh.BA44_exvivo.label --l lh.BA45_exvivo.label --l lh.V1_exvivo.label --l lh.V2_exvivo.label --l lh.MT_exvivo.label --l lh.entorhinal_exvivo.label --l lh.perirhinal_exvivo.label --a BA_exvivo --maxstatwinner --noverbose \n
\n mris_label2annot --s ko --hemi lh --ctab /Applications/freesurfer/average/colortable_BA.txt --l lh.BA1_exvivo.thresh.label --l lh.BA2_exvivo.thresh.label --l lh.BA3a_exvivo.thresh.label --l lh.BA3b_exvivo.thresh.label --l lh.BA4a_exvivo.thresh.label --l lh.BA4p_exvivo.thresh.label --l lh.BA6_exvivo.thresh.label --l lh.BA44_exvivo.thresh.label --l lh.BA45_exvivo.thresh.label --l lh.V1_exvivo.thresh.label --l lh.V2_exvivo.thresh.label --l lh.MT_exvivo.thresh.label --l lh.entorhinal_exvivo.thresh.label --l lh.perirhinal_exvivo.thresh.label --a BA_exvivo.thresh --maxstatwinner --noverbose \n
\n mris_anatomical_stats -th3 -mgz -f ../stats/lh.BA_exvivo.stats -b -a ./lh.BA_exvivo.annot -c ./BA_exvivo.ctab ko lh white \n
\n mris_anatomical_stats -th3 -mgz -f ../stats/lh.BA_exvivo.thresh.stats -b -a ./lh.BA_exvivo.thresh.annot -c ./BA_exvivo.thresh.ctab ko lh white \n
#--------------------------------------------
#@# BA_exvivo Labels rh Wed May  4 23:47:20 PDT 2016
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA1_exvivo.label --trgsubject ko --trglabel ./rh.BA1_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA2_exvivo.label --trgsubject ko --trglabel ./rh.BA2_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA3a_exvivo.label --trgsubject ko --trglabel ./rh.BA3a_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA3b_exvivo.label --trgsubject ko --trglabel ./rh.BA3b_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA4a_exvivo.label --trgsubject ko --trglabel ./rh.BA4a_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA4p_exvivo.label --trgsubject ko --trglabel ./rh.BA4p_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA6_exvivo.label --trgsubject ko --trglabel ./rh.BA6_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA44_exvivo.label --trgsubject ko --trglabel ./rh.BA44_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA45_exvivo.label --trgsubject ko --trglabel ./rh.BA45_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.V1_exvivo.label --trgsubject ko --trglabel ./rh.V1_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.V2_exvivo.label --trgsubject ko --trglabel ./rh.V2_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.MT_exvivo.label --trgsubject ko --trglabel ./rh.MT_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.entorhinal_exvivo.label --trgsubject ko --trglabel ./rh.entorhinal_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.perirhinal_exvivo.label --trgsubject ko --trglabel ./rh.perirhinal_exvivo.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA1_exvivo.thresh.label --trgsubject ko --trglabel ./rh.BA1_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA2_exvivo.thresh.label --trgsubject ko --trglabel ./rh.BA2_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA3a_exvivo.thresh.label --trgsubject ko --trglabel ./rh.BA3a_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA3b_exvivo.thresh.label --trgsubject ko --trglabel ./rh.BA3b_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA4a_exvivo.thresh.label --trgsubject ko --trglabel ./rh.BA4a_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA4p_exvivo.thresh.label --trgsubject ko --trglabel ./rh.BA4p_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA6_exvivo.thresh.label --trgsubject ko --trglabel ./rh.BA6_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA44_exvivo.thresh.label --trgsubject ko --trglabel ./rh.BA44_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.BA45_exvivo.thresh.label --trgsubject ko --trglabel ./rh.BA45_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.V1_exvivo.thresh.label --trgsubject ko --trglabel ./rh.V1_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.V2_exvivo.thresh.label --trgsubject ko --trglabel ./rh.V2_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.MT_exvivo.thresh.label --trgsubject ko --trglabel ./rh.MT_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.entorhinal_exvivo.thresh.label --trgsubject ko --trglabel ./rh.entorhinal_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mri_label2label --srcsubject fsaverage --srclabel /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/fsaverage/label/rh.perirhinal_exvivo.thresh.label --trgsubject ko --trglabel ./rh.perirhinal_exvivo.thresh.label --hemi rh --regmethod surface \n
\n mris_label2annot --s ko --hemi rh --ctab /Applications/freesurfer/average/colortable_BA.txt --l rh.BA1_exvivo.label --l rh.BA2_exvivo.label --l rh.BA3a_exvivo.label --l rh.BA3b_exvivo.label --l rh.BA4a_exvivo.label --l rh.BA4p_exvivo.label --l rh.BA6_exvivo.label --l rh.BA44_exvivo.label --l rh.BA45_exvivo.label --l rh.V1_exvivo.label --l rh.V2_exvivo.label --l rh.MT_exvivo.label --l rh.entorhinal_exvivo.label --l rh.perirhinal_exvivo.label --a BA_exvivo --maxstatwinner --noverbose \n
\n mris_label2annot --s ko --hemi rh --ctab /Applications/freesurfer/average/colortable_BA.txt --l rh.BA1_exvivo.thresh.label --l rh.BA2_exvivo.thresh.label --l rh.BA3a_exvivo.thresh.label --l rh.BA3b_exvivo.thresh.label --l rh.BA4a_exvivo.thresh.label --l rh.BA4p_exvivo.thresh.label --l rh.BA6_exvivo.thresh.label --l rh.BA44_exvivo.thresh.label --l rh.BA45_exvivo.thresh.label --l rh.V1_exvivo.thresh.label --l rh.V2_exvivo.thresh.label --l rh.MT_exvivo.thresh.label --l rh.entorhinal_exvivo.thresh.label --l rh.perirhinal_exvivo.thresh.label --a BA_exvivo.thresh --maxstatwinner --noverbose \n
\n mris_anatomical_stats -th3 -mgz -f ../stats/rh.BA_exvivo.stats -b -a ./rh.BA_exvivo.annot -c ./BA_exvivo.ctab ko rh white \n
\n mris_anatomical_stats -th3 -mgz -f ../stats/rh.BA_exvivo.thresh.stats -b -a ./rh.BA_exvivo.thresh.annot -c ./BA_exvivo.thresh.ctab ko rh white \n
#--------------------------------------------
#@# Hippocampal Subfields processing (T1 + T2 volume) left Wed May  4 23:51:50 PDT 2016
\n /Applications/freesurfer/bin/segmentSF_T1T2.sh /Applications/freesurfer/MCRv80 /Applications/freesurfer ko /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir /Volumes/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/raw/ko/T2w.nii.gz left hippocampus \n
See log file: /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/scripts/hippocampal-subfields-T1T2.log
#--------------------------------------------
#@# Hippocampal Subfields processing (T1 + T2 volume) right Thu May  5 00:43:32 PDT 2016
\n /Applications/freesurfer/bin/segmentSF_T1T2.sh /Applications/freesurfer/MCRv80 /Applications/freesurfer ko /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir /Volumes/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/raw/ko/T2w.nii.gz right hippocampus \n
See log file: /Volumes/group/iang/users/camachoMC/freesurfer_playground/freesurfer6.0/subjs_dir/ko/scripts/hippocampal-subfields-T1T2.log
