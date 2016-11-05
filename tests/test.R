for (file in list.files("./R/", full.names = T)) {
        source(file)
}

assnp_dir = .libPaths()
index_snp = paste0(assnp_dir, "/AlleleSNP/data/input_snps/LUC_Index_SNPs_20160607_short.csv")
# get_assnp_singleBam(index_snp_file = "../00-AlleleSNP/data/input_snps/LUC_Index_SNPs_20160607_short.csv",
                    # bam_dir = "../00-AlleleSNP/data/samples/DDBJ_A549/bam_files/", sample_name = "A549_singleBam")

# get_assnp_sample(index_snp_file = "../00-AlleleSNP/data/input_snps/LUC_Index_SNPs_20160607_short.csv",
#                  sample_name = "A549_bySample", sample_dir = "../00-AlleleSNP/data/samples/DDBJ_A549/")

# get_assnp_encodeDGF(index_snp_file = "../00-AlleleSNP/data/input_snps/LUC_Index_SNPs_20160607_short.csv",
                    # cell_sel = "A549")

# get_assnp_encodeDGF(index_snp_file = "../00-AlleleSNP/data/input_snps/LUC_Index_SNPs_20160607_short.csv",
#                     cell_sel = "A549", use_encode_cnv = F,
#                     vcf_file_for_cnv = "../00-AlleleSNP/data/samples/DDBJ_A549/vcf_files/A549_snv.GATK.formatcor.vcf")
