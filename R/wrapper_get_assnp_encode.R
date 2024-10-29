# M3-get_assnp_encode
## ==================================================================================
## Function: get_ase_snp
###  Aim: To get allele-specific SNP
###  INPUT 1: index SNP list (a .csv file), with the 1st col the rsID and the 2nd col Population
###  INPUT 2: snp information table (a .csv file), contains column names with rsID and population
###  OUTPUT: potential allele-specific effect (ase) snps
## ==================================================================================

# source function
# source("./src/scripts/M1-get_assnp.R")

# wrapper function
get_assnp_encodeDGF = function(index_snp_file,
                               snp_info_file = NA,
                               cell_sel = "",
                               sample_name = "ENCODE_DGF",
                               server = "http",
                               bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDgf/",
                               use_encode_cnv = T,
                               ...) {

        get_assnp(index_snp_file = index_snp_file,
                  snp_info_file = snp_info_file,
                  sample_name = sample_name,
                  cell_sel = cell_sel,
                  server = server,
                  bam_dir = bam_dir,
                  genotype_by_sample = F,
                  use_encode_cnv = use_encode_cnv,
                  ...)
}


# wrapper function
get_assnp_encodeDNase = function(index_snp_file,
                                 snp_info_file = NA,
                                 cell_sel = "",
                                 sample_name = "ENCODE_DNase",
                                 server = "http",
                                 bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/",
                                 use_encode_cnv = T,
                                 ...) {

        get_assnp(index_snp_file = index_snp_file,
                  snp_info_file = snp_info_file,
                  sample_name = sample_name,
                  cell_sel = cell_sel,
                  server = server,
                  bam_dir = bam_dir,
                  genotype_by_sample = F,
                  use_encode_cnv = use_encode_cnv,
                  ...)
}


# wrapper function
get_assnp_encodeHaibTFBS = function(index_snp_file,
                                    snp_info_file = NA,
                                    cell_sel = "",
                                    sample_name = "ENCODE_HaibTFBS",
                                    server = "http",
                                    bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/",
                                    use_encode_cnv = T,
                                    ...) {

        get_assnp(index_snp_file = index_snp_file,
                  snp_info_file = snp_info_file,
                  sample_name = sample_name,
                  cell_sel = cell_sel,
                  server = server,
                  bam_dir = bam_dir,
                  genotype_by_sample = F,
                  use_encode_cnv = use_encode_cnv,
                  ...)
}


# wrapper function
get_assnp_encodeSydhTFBS = function(index_snp_file,
                                    snp_info_file = NA,
                                    cell_sel = "",
                                    sample_name = "ENCODE_SydhTFBS",
                                    server = "http",
                                    bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/",
                                    use_encode_cnv = T,
                                    ...) {

        get_assnp(index_snp_file = index_snp_file,
                  snp_info_file = snp_info_file,
                  sample_name = sample_name,
                  cell_sel = cell_sel,
                  server = server,
                  bam_dir = bam_dir,
                  genotype_by_sample = F,
                  use_encode_cnv = use_encode_cnv,
                  ...)
}


# wrapper function
get_assnp_encodeUChicagoTFBS = function(index_snp_file,
                                        snp_info_file = NA,
                                        cell_sel = "",
                                        sample_name = "ENCODE_UChicagoTFBS",
                                        server = "http",
                                        bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUchicagoTfbs/",
                                        use_encode_cnv = T,
                                        ...) {

        get_assnp(index_snp_file = index_snp_file,
                  snp_info_file = snp_info_file,
                  sample_name = sample_name,
                  cell_sel = cell_sel,
                  server = server,
                  bam_dir = bam_dir,
                  genotype_by_sample = F,
                  use_encode_cnv = use_encode_cnv,
                  ...)
}


# wrapper function
get_assnp_encodeBroadHistone = function(index_snp_file,
                                        snp_info_file = NA,
                                        cell_sel = "",
                                        sample_name = "ENCODE_BroadHistone",
                                        server = "http",
                                        bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/",
                                        use_encode_cnv = T,
                                        ...) {

        get_assnp(index_snp_file = index_snp_file,
                  snp_info_file = snp_info_file,
                  sample_name = sample_name,
                  cell_sel = cell_sel,
                  server = server,
                  bam_dir = bam_dir,
                  genotype_by_sample = F,
                  use_encode_cnv = use_encode_cnv,
                  ...)
}


# wrapper function
get_assnp_encodeUWHistone = function(index_snp_file,
                                     snp_info_file = NA,
                                     cell_sel = "",
                                     sample_name = "ENCODE_UWHistone",
                                     server = "http",
                                     bam_dir = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/",
                                     use_encode_cnv = T,
                                     ...) {

        get_assnp(index_snp_file = index_snp_file,
                  snp_info_file = snp_info_file,
                  sample_name = sample_name,
                  cell_sel = cell_sel,
                  server = server,
                  bam_dir = bam_dir,
                  genotype_by_sample = F,
                  use_encode_cnv = use_encode_cnv,
                  ...)
}

