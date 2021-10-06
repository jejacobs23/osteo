##library(survminer)
library(survival)

args = commandArgs(trailingOnly=TRUE)
working_dir <- args[1]
Pname <- args[2]
Sname <- args[3]
outfile_name <- args[4]

COM_DIR <- "/home/exacloud/lustre1/jjacobs"
WORKING_DIR <- paste(COM_DIR, "/data/osteo/Reactome_Analyses/", working_dir, sep = "")

pID <- c()
p_val <- c()
aberrant <- c()
non_aberrant <- c()
P <- read.table(paste(WORKING_DIR, Pname, sep = ""), header = T, sep = "\t")

print(P$Reaction_ID)

for (i in P$Reaction_ID) {
    pID <- append(pID, i)
    S <- read.table(paste(WORKING_DIR, "/Survival_DataFrame_", i, Sname, sep = ""), header = T)

##    ggsurvplot(
##      fit = survfit(Surv(death, status) ~ group, data = S),
##      pval = TRUE,
##      title = pID,
##      pval.method = TRUE,
##      xlab = "Months",
##      ylab = "Overall Survival Probability")

    print(i)
    group1 <- length(S$group[S$group == "aberrant"])
    group2 <- length(S$group[S$group == "non-aberrant"])
    aberrant <- append(aberrant, group1)
    non_aberrant <- append(non_aberrant, group2)
    if (group1 == 0) {
        p <- "Only_One_Group"
    } else if (group2 == 0) {
        p <- "Only_One_Group"
    } else {
        F <- survdiff(Surv(death, status) ~ group, data = S)
        p <- pchisq(F$chisq, df=1, lower.tail=FALSE)
    }
    print(p)
    p_val <- append(p_val, p)
}
fdr <- p.adjust(p_val, "BH")
df <- data.frame(pID, aberrant, non_aberrant, p_val, fdr)
write.table(df, file = paste(WORKING_DIR, "/Pathway_p-values_", outfile_name, sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
