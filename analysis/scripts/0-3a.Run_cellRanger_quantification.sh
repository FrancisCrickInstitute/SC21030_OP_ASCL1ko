#!/bin/sh

module purge
module load CellRanger/5.0.0-bcl2fastq-2.20.0


sample_id=$1
fastq_loc=$2
mkdir -p /camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/derived_data/CellRanger_quantification
cd /camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/derived_data/CellRanger_quantification
transcriptome="/camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/CellRanger_Homo_sapiens.GRCh37.75" 

instruction="cellranger count --id=$sample_id \
	--transcriptome=$transcriptome \
	--fastqs=$fastq_loc \
	--sample=$sample_id --project=SC21030"

echo $instruction

cellranger count --id=$sample_id \
	--transcriptome=$transcriptome \
	--fastqs=$fastq_loc \
	--sample=$sample_id --project=SC21030


