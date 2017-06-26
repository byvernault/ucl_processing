python Spider_Mprage2Gad_Reg_v1_0_0.py \
-d ~/data/jobsdir/test_mprage2gad/ \
-a IMPLAN-x-thok-x-thok_20151203-x-11-x-Mprage2Gad_Reg_v1 \
--exe perform_mprage2gad_registration.py \
--host $XE \
--gad xnat:/project/IMPLAN/subject/thok/experiment/thok_20151203/scan/10/resource/NIFTI \
--mprage xnat:/project/IMPLAN/subject/thok/experiment/thok_20151203/scan/5/resource/NIFTI \
--labels xnat:/project/IMPLAN/subject/thok/experiment/thok_20151203/assessor/IMPLAN-x-thok-x-thok_20151203-x-5-x-GIF_Parcellation_v3/resource/LABELS \
--brain xnat:/project/IMPLAN/subject/thok/experiment/thok_20151203/assessor/IMPLAN-x-thok-x-thok_20151203-x-5-x-GIF_Parcellation_v3/resource/BRAIN