# M1-get_assnp_local
## ==================================================================================
## Function: get_assnp_local
###  Aim: To get allele-specific SNP using local file
###  INPUT 1: index SNP list (a .csv file), with the 1st col the rsID and the 2nd col Population
###  INPUT 2: snp information table (a .csv file), contains column names with rsID and population
###  OUTPUT: potential allele-specific effect (ase) snps
## ==================================================================================

# source function
# source("./src/scripts/M1-get_assnp.R")


# wrapper function: get_assnp_byBam ---------------------------------------------------------------------------

# get ase snp just by .bam files, not using .vcf files; no information of peak, and cnv added
get_assnp_byBam = function(index_snp_file,
                               snp_info_file = NA,
                               sample_name = "",
                               bam_dir,
                               merge_replicates = F) {

        get_assnp(index_snp_file = index_snp_file,
                  snp_info_file = snp_info_file,
                  sample_name = sample_name,
                  bam_dir = bam_dir,
                  genotype_by_sample = F,
                  merge_replicates = merge_replicates)
}


# wrapper function: get_assnp_sample ------------------------------------------------------------------------------

# get ase snp just by .bam files and .vcf files; with information of peak, and cnv added; ase are corrected by cnv
get_assnp_bySample = function(index_snp_file,
                            snp_info_file = NA,
                            sample_name = "",
                            sample_dir = "./data/samples/DDBJ_A549/",
                            merge_replicates = F) {

        bam_dir = paste0(sample_dir, "/bam_files/")
        peak_dir = paste0(sample_dir, "/peak_files/")
        vcf_dir = paste0(sample_dir, "/vcf_files/")
        vcf_file_for_cnv = list.files(vcf_dir, "[Ww][Gg][Ss].*vcf$", full.names = T)

        get_assnp(index_snp_file = index_snp_file,
                  snp_info_file = snp_info_file,
                  sample_name = sample_name,
                  bam_dir = bam_dir,
                  peak_dir = peak_dir,
                  vcf_dir = vcf_dir,
                  vcf_file_for_cnv = vcf_file_for_cnv,
                  genotype_by_sample = T,
                  merge_replicates = merge_replicates)
}
