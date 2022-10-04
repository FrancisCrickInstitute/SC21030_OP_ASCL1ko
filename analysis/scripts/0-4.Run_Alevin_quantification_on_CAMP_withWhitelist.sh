#!/bin/sh

module purge
module load Anaconda3/2020.02

source activate salmon

read1_info="$1"
read2_info="$2"
sample_name="$3"
cellranger_dir="$4"


genome_index="/camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.expanded"
tgMap_loc="/camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/raw_data/genome_info/human/Homo_sapiens.GRCh37.75.annotation.expanded.tx2gene.tsv"
alevin_out="/camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/data/derived_data/alevin_quantification_whitelist/${sample_name}"
cellranger_out="/camp/stp/babs/working/ghanata/projects/guillemotf/oana.paun/ASCL1_KO_cortical_neural_differentiation_day24_with_Notch_inhibition/${cellranger_dir}/${sample_name}/outs/filtered_feature_bc_matrix/"
cellranger_zbarcode="${cellranger_out}/barcodes.tsv.gz"
cellranger_barcode_temp="${cellranger_out}/barcodes-temp.tsv"
cellranger_barcode="${cellranger_out}/barcodes.tsv"

# Unzip barcode file from Cellranger output
gzip -d $cellranger_zbarcode
cp $cellranger_barcode $cellranger_barcode_temp

# Gunzip barcode file from Cellranger output
gzip $cellranger_barcode

# Remove trailing "-1"
sed -i "s/-1//g" $cellranger_barcode_temp

# Note ~10 threads are recommended by Alevin documentation
instruction="salmon alevin -l ISR -i $genome_index \
    -1 $read1_info \
    -2 $read2_info \
    -o $alevin_out \
    -p 12 \
    --tgMap $tgMap_loc \
    --chromiumV3 --dumpFeatures \
    --whitelist $cellranger_barcode_temp"

echo $instruction

salmon alevin -l ISR -i $genome_index \
    -1 $read1_info \
    -2 $read2_info \
    -o $alevin_out \
    -p 12 \
    --tgMap $tgMap_loc \
    --chromiumV3 --dumpFeatures \
    --whitelist $cellranger_barcode_temp
    
# Remove barcode file
rm $cellranger_barcode_temp