with_meCpG = read.table("barcodes_with_meCpGs.txt")
with_or_without_meCpG = read.table("barcodes_with_or_without_meCpGs.txt")

AC1 = read.table("mapping/AC1.ins.gz")
AC2 = read.table("mapping/AC2.ins.gz")
AG1 = read.table("mapping/AG1.ins.gz")
AG2 = read.table("mapping/AG2.ins.gz")
GT1 = read.table("mapping/GT1.ins.gz")
GT2 = read.table("mapping/GT2.ins.gz")
TC1 = read.table("mapping/TC1.ins.gz")
TC2 = read.table("mapping/TC2.ins.gz")
TG1 = read.table("mapping/TG1.ins.gz")
TG2 = read.table("mapping/TG2.ins.gz")

AC = rbind(AC1, AC2)
AG = rbind(AG1, AG2)
GT = rbind(GT1, GT2)
TC = rbind(TC1, TC2)
TG = rbind(TG1, TG2)

AC_meCpG = subset(AC, V1 %in% with_meCpG$V2[with_meCpG$V1 == "AC"])
AC_all = subset(AC, V1 %in% with_or_without_meCpG$V2[with_meCpG$V1 == "AC"])
AG_meCpG = subset(AG, V1 %in% with_meCpG$V2[with_meCpG$V1 == "AG"])
AG_all = subset(AG, V1 %in% with_or_without_meCpG$V2[with_meCpG$V1 == "AG"])
GT_meCpG = subset(GT, V1 %in% with_meCpG$V2[with_meCpG$V1 == "GT"])
GT_all = subset(GT, V1 %in% with_or_without_meCpG$V2[with_meCpG$V1 == "GT"])
TC_meCpG = subset(TC, V1 %in% with_meCpG$V2[with_meCpG$V1 == "TC"])
TC_all = subset(TC, V1 %in% with_or_without_meCpG$V2[with_meCpG$V1 == "TC"])
TG_meCpG = subset(TG, V1 %in% with_meCpG$V2[with_meCpG$V1 == "TG"])
TG_all = subset(TG, V1 %in% with_or_without_meCpG$V2[with_meCpG$V1 == "TG"])

mappings_with_meCpG = rbind(
   AC_meCpG, AG_meCpG, GT_meCpG, TC_meCpG, TG_meCpG)
mappings_with_or_without_meCpG = rbind(
   AC_all, AG_all, GT_all, TC_all, TG_all)
mappings_without_meCpG = subset(mappings_with_or_without_meCpG,
   !(V1 %in% mappings_with_meCpG$V1))

write.table(mappings_with_meCpG, file="mappings_with_meCpG.txt",
   sep="\t", quote=FALSE, row.names=FALSE)
write.table(mappings_without_meCpG, file="mappings_without_meCpG.txt",
   sep="\t", quote=FALSE, row.names=FALSE)
