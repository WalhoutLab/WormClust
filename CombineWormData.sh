#!/bin/bash
# https://libsleipnir.bitbucket.io/index.html#se
source /project/csbio/shivani_umass/.bashrc
export MATLABPATH=/project/csbio/shivani_umass/Worms/Final_scripts_070920/

## set up directories
#cd /project/csbio/shivani_umass/Worms
#mkdir FinalData_070920
#mkdir FinalData_070920/all_data
#mkdir FinalData_070920/10_conditions
## Step 1. Download data
# total 400 datasets
cd /project/csbio/shivani_umass/Worms/FinalData_070920/all_data
wget ftp://caltech.wormbase.org/pub/wormbase/spell_download/datasets/WB*csv
tar -xpvf AllDatasetsDownload.tgz

## Step 2.Filter gene expression data
# only keep gene expression data that includes at least 100 genes and at least 10 conditions
# copy these files to /project/csbio/shivani_umass/Worms/FinalData_070920/10_conditions
# these criteria can be modified
nice matlab -nodisplay -nodesktop -nosplash -r "check_10_csv_file(10,100);exit" </dev/null> /dev/null

## Step 3. convert CSV file to PCL file
cd /project/csbio/shivani_umass/Worms/FinalData_070920/10_conditions
nice matlab -nodisplay -nodesktop -nosplash -r "csv2pcl;exit" </dev/null> /dev/null

## Step 4. Normalize expression data of C elegans individually using z scores
cd /project/csbio/shivani_umass/Worms/FinalData_070920/10_conditions
for pclfile in `ls WBPaper*ce*.pcl`
do
        if [ ! -f "z_normalized_$pclfile.dab" ]; then 
                Normalizer -t pcl -z -i $pclfile -o z_normalized_$pclfile      
        fi
done

## Step 5 Keep only microarray and RNASeq datasets. Remove any other kind of datasets
rm z_normalized*ms*pcl
rm z_normalized*tr.pcl

## Step 6. Combine the z-normalised datasets 
Combiner -t pcl  -o /project/csbio/shivani_umass/Worms/FinalData_070920/10_conditions/combined_z_normalised.pcl   /project/csbio/shivani_umass/Worms/FinalData_070920/10_conditions/z_normalized_WBPaper000*.pcl

# Step 7. Impute the missing data with 10 nearest neighbours and remove genes with less than 70% of missing data (optional step)
#KNNImputer -i combined_z_normalised.pcl -o imputed_combined_z_normalised.pcl

# Step 8. Calculate pearson and spearman correlations between all genes
Distancer -d pearson -z off  -g /project/csbio/shivani_umass/Worms/MetabolicGenes120120.txt -i combined_z_normalised.pcl -o pearson_combined_z_normalised.dab
Distancer -d pearsig -z off -g /project/csbio/shivani_umass/Worms/MetabolicGenes120120.txt -i combined_z_normalised.pcl -o pearsig_combined_z_normalised.dab
Dat2Dab -i pearson_combined_z_normalised.dab -o pearson_combined_z_normalized.dat
Dat2Dab -i pearsig_combined_z_normalised.dab -o pearsig_combined_z_normalized.dat

Distancer -d spearman -z off  -g /project/csbio/shivani_umass/Worms/MetabolicGenes120120.txt -i combined_z_normalised.pcl -o Spearman_combined_z_normalised.dab
Dat2Dab -i Spearman_combined_z_normalised.dab -o Spearman_combined_z_normalized.dat

#step 9. Extract information for only metabolic genes
##python MetabolicRelationshipsOnly.py pearson_imputed_combined_z_normalized.dat pearsoncorr_metabolic.dat
#python MetabolicRelationshipsOnly.py pearsig_imputed_combined_z_normalized.dat pearsig_metabolic.dat
#python MetabolicRelationshipsOnly.py
