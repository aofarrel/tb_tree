#/bin/bash
#_FILES = /more_storage/lily/tb_tree/parse_vcf/sample_files/
bcftools query -l /more_storage/lily/tb_tree/temp_loc_for_megavcf/merged.vcf.bgz > /more_storage/lily/tb_tree/temp_loc_for_megavcf/samples.txt
split -l 100 -d /more_storage/lily/tb_tree/temp_loc_for_megavcf/samples.txt /more_storage/lily/tb_tree/parse_vcf/sample_files/split_

for file in /more_storage/lily/tb_tree/parse_vcf/sample_files/*;
do 
    echo $file;
done