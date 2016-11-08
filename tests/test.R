# library(devtools)
#
# install_github("foreverycc/AlleleSNP_Package")
# for (file in list.files("./R/", pattern = ".*R", full.names = T)) {
#         print (file)
#         source(file)
# }

library(AlleleSNP)
assnp_dir = .libPaths()
index_snp_file = paste0(assnp_dir, "/AlleleSNP/extdata/input_snps/input_snp_example2.csv")
read.csv(index_snp_file, header = F)

sample_dir = paste0(assnp_dir, "/AlleleSNP/extdata/sample/A549")

# get_assnp_singleBam(index_snp_file = "../00-AlleleSNP/data/input_snps/LUC_Index_SNPs_20160607_short.csv",
                    # bam_dir = "../00-AlleleSNP/data/samples/DDBJ_A549/bam_files/", sample_name = "A549_singleBam")

get_assnp_sample(index_snp_file = "../00-AlleleSNP/data/input_snps/LUC_Index_SNPs_20160607_short.csv",
                 sample_name = "A549_bySample", sample_dir = "../00-AlleleSNP/data/samples/DDBJ_A549/")

# get_assnp_sample(index_snp_file = index_snp_file,
#                  sample_name = "A549", sample_dir = sample_dir)

# get_assnp_encodeDGF(index_snp_file = "../00-AlleleSNP/data/input_snps/LUC_Index_SNPs_20160607_short.csv",
                    # cell_sel = "A549")

# get_assnp_encodeDGF(index_snp_file = "../00-AlleleSNP/data/input_snps/LUC_Index_SNPs_20160607_short.csv",
#                     cell_sel = "A549", use_encode_cnv = F,
#                     vcf_file_for_cnv = "../00-AlleleSNP/data/samples/DDBJ_A549/vcf_files/A549_snv.GATK.formatcor.vcf")
