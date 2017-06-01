python Spider_Compute_ADC_Verdict_v1_0_0.py \
-a INNOVATE-x-INN-047-FCO-x-INN-047-FCO_WALK__20160630-x-Compute_ADC_Verdict_v1 \
-d /Users/byvernault/data/jobsdir/test_adc_verdict \
--dcm_file xnat:/project/INNOVATE/subject/INN-047-FCO/experiment/INN-047-FCO_WALK__20160630/scan/801/resource/DICOM \
--acquisitions xnat:/project/INNOVATE/subject/INN-047-FCO/experiment/INN-047-FCO_WALK__20160630/assessor/INNOVATE-x-INN-047-FCO-x-INN-047-FCO_WALK__20160630-x-Registration_Verdict_v1/out/resource/ACQ1,xnat:/project/INNOVATE/subject/INN-047-FCO/experiment/INN-047-FCO_WALK__20160630/assessor/INNOVATE-x-INN-047-FCO-x-INN-047-FCO_WALK__20160630-x-Registration_Verdict_v1/out/resource/ACQ2 \
--matlab_src /Users/byvernault/home-local/code/matlab/VERDICT/ADC_MAP/v1/ \
--camino /Users/byvernault/home-local/softwares/caminoLaura \
--scheme_filename /Users/byvernault/home-local/code/matlab/VERDICT/Noptimised/NOptimisedADC_IN.scheme \
--host https://prostate-xnat.cs.ucl.ac.uk