#!/bin/sh

module purge
module load Anaconda3/2020.02
source activate salmon

#cd ../data/raw_data/genome_info/human/

#grep ">" Homo_sapiens.GRCh38.dna_sm.primary_assembly_CreERT.fa | cut -d ">" -f 2 | cut -d " " -f 1 > Homo_sapiens.GRCh38.dna_sm.primary_assembly_CreERT.chrnames.txt

cat /camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.expanded.fa /camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa > /camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/combined.fa

salmon index \
    --transcripts /camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/combined.fa \
    --index /camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.expanded \
    --gencode \
    -p 12 \
    --decoys /camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.chrnames.txt
    
    
rm /camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/combined.fa

conda deactivate