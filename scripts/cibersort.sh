#!/bin/bash

docker run -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/intermediate_files/cibersort:/src/data -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/intermediate_files/cibersort/outputs:/src/outdir cibersortx/fractions --username k.ryan45@nuigalway.ie --token b7f03b943ade9b4146dc2126b4ac9d19 --refsample caf_tpm_for_signature_engs_version_batch_corrected_2022-10-03.txt --mixture caf_tpm_mixture_engs_version_batch_corrected_2022-10-03.txt --rmbatchBmode FALSE --phenoclasses phenoclasses_caf_2022-10-03.txt --perm 500 --verbose TRUE --G.min 300 --G.max 500
