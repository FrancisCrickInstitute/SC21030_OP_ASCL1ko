#!/bin/sh

module purge
module load CellRanger/5.0.0-bcl2fastq-2.20.0

cd /camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/

instruction="cellranger mkref \
    --genome=CellRanger_Homo_sapiens.GRCh37.75 \
    --fasta=/camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa \
    --genes=/camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/Homo_sapiens.GRCh37.75_filtered.gtf \
    --nthreads 16"

echo $instruction

cellranger mkref \
    --genome=CellRanger_Homo_sapiens.GRCh37.75 \
    --fasta=/camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa \
    --genes=/camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/Homo_sapiens.GRCh37.75_filtered.gtf \
    --nthreads 16

