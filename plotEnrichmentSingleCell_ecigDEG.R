
# Read the data

setwd("/genomics/chiragProjects/bgiStereoseq/analysisGwas/3months_snRNAseq/singleCellDEG")
A <- read.table("count", header = FALSE, stringsAsFactors = FALSE)
colnames(A) <- c("File", "A", "B", "C", "D")  # Rename columns for clarity

# Initialize vectors for p-value and odds ratio
p_values <- numeric(nrow(A))
odds_ratios <- numeric(nrow(A))

# Perform Fisher's test for each row
for (i in 1:nrow(A)) {
  contingency_table <- matrix(c(A$A[i], A$B[i], A$C[i], A$D[i]), nrow = 2)
  test_result <- fisher.test(contingency_table)
  p_values[i] <- test_result$p.value
  odds_ratios[i] <- test_result$estimate
}

# Add results to the data frame
A$p_value <- p_values
A$odds_ratio <- odds_ratios
# Adjust p-values for multiple testing (e.g., Benjamini-Hochberg)
A$adjusted_p_value <- p.adjust(A$p_value, method = "BH")
write.table(A, file = "fisher_test_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)


cleaned_names <- gsub("_DEGcombined\\.csv", "", A$File)
Aorder <- A[order(A$adjusted_p_value), ]
pdf (file="FigDEG.pdf")
par (mar=c(5,9,3,4))
# Create the barplot with cleaned names
barplot(-log10(rev(Aorder$adjusted_p_value)), horiz = TRUE, xlim = c(0, 20),names.arg = rev(cleaned_names), las = 1)
box()
abline(v = -log10(0.05), lty = 2)

dev.off()

S = read.table("overalpdeg_ASD_DEGcombined", header=FALSE)
S1=table(S$V8)
S1=sort(S1, decreasing = TRUE)
pdf (file="FigSCell.pdf")
barplot(S1, las=2, ylim=c(0,250) )
dev.off()





