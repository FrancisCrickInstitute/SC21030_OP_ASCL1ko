#!/bin/sh

module purge
module load Anaconda3/2020.02
source activate scvelo-0.2.2

Rscript -e "rmarkdown::render('/camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/analysis/0-2.Generate_index_for_Alevin.Rmd')"

conda deactivate