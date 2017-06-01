python Spider_GIF_Parcellation_v3_0_0.py \
-a DIAN-x-09NXRZ-x-09NXRZ_v00_mr-x-6-x-GIF_Parcellation_v2 \
-d /Users/byvernault/data/run_niftipipe/prion_GIF \
--gif perform_gif_propagation.py \
--dbt /Users/byvernault/data/db/db.xml \
--omp 4 \
--working_dir /tmp/prion_GIF/tmp_files \
--t1 xnat:/project/DIAN/subject/09NXRZ/experiment/09NXRZ_v00_mr/scan/6/resource/NIFTI \
--host $XC
