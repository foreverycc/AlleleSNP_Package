# define the "indexSNP" class
setClass("indexSNP", 
         representation(rsID = "character", 
                        ldSNPs = "vector",
                        ldInfo = "vector")
)

setGeneric("ldRelation", function(snp1, snp2) standardGeneric("ldRelation"))
setMethod("ldRelation", "indexSNP", function(snp1, snp2) {
        if (snp1@rsID %in% snp2@ldSNPs) {
                pos = which(snp2@ldSNPs == snp1@rsID )                
                return (ifelse(snp2@ldInfo[pos] > 0, 1, -1))
        } else if (any (snp1@ldSNPs %in% snp2@ldSNPs)) {
                pos1 = which(snp1@ldSNPs %in% snp2@ldSNPs)[1]
                inter_snp = snp1@ldSNPs[pos1]
                pos2 = which(snp2@ldSNPs == inter_snp)
                return (ifelse(snp1@ldInfo[pos1] * snp2@ldInfo[pos2] > 0, 1, -1))
        } else return (NA)
})