# loop through each file in the first column of infiles.txt and run cibersort 
# on each file using the CIBERSORT docker image
#	- infiles.txt is a list of files to run cibersort on
#	- outdir is the directory to save the output files to	
#       example command to run CIBERSORT: 
# 	docker run -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction:/src/data -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/cibersort_outputs:/src/outdir cibersortx/fractions --username <username> --token <token> --refsample subtype_data.txt --phenoclasses phenoclasses.txt --mixture inhouse_data.txt --perm 100 --QN FALSE --rmbatchBmode TRUE

# first try on one file
cd /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction
start_time=$SECONDS
outdir_base=/home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs_20231211_different_inhouse_samples/cibersort_outputs
# loop through each row in the dirname column of infiles.csv and concatenate the string to outdir_base
while IFS=, read -r filename_subpop filename_inhouse dirname
do
    outdir=$outdir_base/$dirname
    echo $outdir
    # run cibersort on the file
    docker run -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction:/src/data \
           -v $outdir:/src/outdir \
           cibersortx/fractions \
           --username k.ryan45@nuigalway.ie \
           --token 2ea740c2ffa1a6c10741bc7799dcbacf \
           --mixture $filename_inhouse \
           --refsample $filename_subpop \
           --phenoclasses phenoclasses_caf.txt \
           --QN FALSE \
           --perm 500 \
           --G.min 300 \
           --G.max 500 
done < 2023-12-11_filenames_limma_remove_batch_effect_batchcorrection_combined_inhouse_subpop.csv

  #--rmbatchBmode TRUE \


#docker run -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction:/src/data \
#	   -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/cibersort_outputs:/src/outdir \
#           cibersortx/fractions \
#           --username k.ryan45@nuigalway.ie \
#           --token 2ea740c2ffa1a6c10741bc7799dcbacf \
#           --refsample 2023-12-08_limma_remove_batch_effect_remove_study_remove_patient_.txt \
#           --phenoclasses phenoclasses_caf.txt \
#           --QN FALSE
#           #--mixture 2023-12-08_inhouse_data_for_cibersort.txt \
#
#cp ../cibersort_outputs/CIBERSORTx_phenoclasses_caf.CIBERSORTx_2023-12-08_limma_remove_batch_effect_remove_study_remove_patient_.bm.K999.txt .
##CIBERSORTx_phenoclasses_caf.CIBERSORTx_2023-12-08_limma_remove_batch_effect_remove_study_remove_patient_.bm.K999.txt
#
#docker run -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction:/src/data \
#           -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/cibersort_outputs/cibersort_remove_study_remove_patient_outputs:/src/outdir \
#           cibersortx/fractions \
#           --username k.ryan45@nuigalway.ie \
#           --token 2ea740c2ffa1a6c10741bc7799dcbacf \
#           --mixture 2023-12-08_inhouse_data_for_cibersort.txt \
#	   --sigmatrix CIBERSORTx_phenoclasses_caf.CIBERSORTx_2023-12-08_limma_remove_batch_effect_remove_study_remove_patient_.bm.K999.txt \
#	   --rmbatchBmode TRUE \
#           --perm 500 \
#           --G.min 300 \
#           --G.max 500

#docker run -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/limma_batch_correction:/src/data \
#           -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/nov2023-match-patient-origin/outputs/cibersort_outputs/cibersort_remove_patient_only_outputs:/src/outdir \
#           cibersortx/fractions \
#           --username k.ryan45@nuigalway.ie \
#           --token 2ea740c2ffa1a6c10741bc7799dcbacf \
#           --mixture 2023-12-08_inhouse_data_for_cibersort.txt \
#           --rmbatchBmode TRUE \
#           --refsample 2023-12-08_limma_remove_batch_effect_remove_patient_.txt \
#           --phenoclasses phenoclasses_caf.txt \
#	   --QN FALSE \
#           --perm 500 \
#           --G.min 300 \
#           --G.max 500	
