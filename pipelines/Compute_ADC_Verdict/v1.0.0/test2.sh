python Spider_Compute_ADC_Verdict_v1_0_0.py \
-a INNOVATE-x-INN-001-KRE-x-INN-001-KRE_20160101-x-Compute_ADC_Verdict_v1 \
-d /Users/byvernault/data/jobsdir/test_adc_verdict \
--dcm_file xnat:/project/INNOVATE/subject/INN-001-KRE/experiment/INN-001-KRE_20160101/scan/801/resource/DICOM \
--acquisitions xnat:/project/INNOVATE/subject/INN-001-KRE/experiment/INN-001-KRE_20160101/assessor/INNOVATE-x-INN-001-KRE-x-INN-001-KRE_20160101-x-Registration_Verdict_v1/out/resource/ACQ1 \
--matlab_src /Users/byvernault/home-local/code/matlab/VERDICT/ADC_MAP/v1/ \
--camino /Users/byvernault/home-local/softwares/caminoLaura \
--scheme_filename /Users/byvernault/home-local/code/matlab/VERDICT/Noptimised/NOptimisedADC_IN.scheme \
--host https://prostate-xnat.cs.ucl.ac.uk