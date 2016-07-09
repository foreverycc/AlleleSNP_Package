# F2.5-get_encodeCnv_info
## ==================================================================================
## Function: get_encodeCnv_info_main
###  Aim: To get the cnv of input SNPs by vcf/encode cnv data files
###  INPUT: input snp file (.csv), vcf/encode cnv data
###  OUTPUT: input snp file with information of cnv
## ==================================================================================

# helper functions ------------------------------------------------------------------------------------------------

# source function: gen_output_file_cnvInfo
# source("./src/scripts/F2.4-get_cnv_info.R")

# 1. source function: read_inputSNP_file
# source ("./src/scripts/F2.1-get_alleleDist_info.R")

# 2. generate snp-cnv table
# Aim: to generate snp-cnv table: rows - snps; cols - encode samples;
gen_snp_cnv_table = function(snp_info_df, snp_info_gr) {

        cat("Add ENCODE cnv information...\n")
        # encode_cnv_gr_list = readRDS(cnv_grList_file)
        data("encode_cnv_gr_list")
        snp_cnv_table = data.frame(matrix(nrow = length(snp_info_gr), ncol = length(encode_cnv_gr_list)))
        colnames(snp_cnv_table) = names(encode_cnv_gr_list)

        for (i in 1:length(encode_cnv_gr_list)) {
                # ----- test ------ #
                # print (i)
                # ----------------- #
                overlaps_df_i = as.data.frame(findOverlaps(snp_info_gr, encode_cnv_gr_list[[i]]))
                cnv_snp_info_i = encode_cnv_gr_list[[i]][overlaps_df_i$subjectHits]$name
                snp_cnv_table[overlaps_df_i$queryHits, i] = cnv_snp_info_i
        }

        # merge data frame
        snp_cnv_df = data.frame(snp_info_df, snp_cnv_table)

        return(snp_cnv_df)
}


# main function ---------------------------------------------------------------------------------------------------

get_encodeCnv_info_main = function(snp_info_file, output_dir = "./", sample_name = "", output_file = NA) {

# Test #
# snp_info_file = "./data/haploreg_files/LUC_Index+LD_SNPs_20160607.csv"
# snp_info_alleleDist_file = "./tests/test6/LUC_Index_SNPs_20160607_short_ENCODE_DNase_assnp/LUC_Index_SNPs_20160607_short_riskPop_0.5_ENCODE_DNase_alleleDist.csv"
# snp_info_alleleDist_df = read.csv(snp_info_alleleDist_file)

        # 1. read snp info file
        snp_info_list = read_inputSNP_file(snp_info_file)

        # 2. generate snp-cnv table (row: snps, col: all encode cnv samples)
        snp_info_addCnv_df = gen_snp_cnv_table(snp_info_df = snp_info_list$snp_info_df,
                                               snp_info_gr = snp_info_list$snp_info_gr)
        # head(snp_cnv_df)

        # 3. write down snp-cnv table
        if (is.na(output_file)) {
                output_file = gen_output_file_cnvInfo(snp_info_file = snp_info_file,
                                                      output_dir = output_dir,
                                                      sample_name = sample_name)
        }
        cat("    output file name:", output_file, '\n')

        if (output_file != F) {
                write.csv0(snp_info_addCnv_df, output_file)
        }
        cat("cnv information added ... \n")

        return(list(snp_info_addCnv_df = snp_info_addCnv_df, output_file = output_file))
}


# Test #
# get_encodeCnv_info_main("./data/haploreg_files/LUC_Index+LD_SNPs_20160607.csv", sample_name = "ENCODE")
