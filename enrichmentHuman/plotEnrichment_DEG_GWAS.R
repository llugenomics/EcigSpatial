

setwd("/genomics/chiragProjects/bgiStereoseq/analysisGwas/3months_snRNAseq/gwas")
inFile="outputOverlapDEG_GWAS"
outFile=paste(inFile,"_pvalue",sep="")

df = read.table(inFile, header=FALSE)
df=df[df$V2>0 & df$V3 >0 & df$V4>0 & df$V5>0,]

# Create an empty data frame to store the results
df_out <- data.frame(Variant = character(nrow(df)),
                     Count1 = numeric(nrow(df)),
                     Count2 = numeric(nrow(df)),
                     Count3 = numeric(nrow(df)),
                     Count4 = numeric(nrow(df)),
                     P_Value_Raw = numeric(nrow(df)),
                     P_Value_Corrected = numeric(nrow(df)),
                     Odds_Ratio = numeric(nrow(df)))
# Store raw p-values
p_values_raw <- numeric(nrow(df))

for (i in 1:nrow(df)) {
  values <- as.numeric(df[i, -1])  # Exclude the first column (variant names) and convert to numeric
  result <- fisher.test(matrix(values, nrow = 2))
  p_values_raw[i] <- result$p.value
  odds_ratio <- result$estimate
  # Store the results in the df_out data frame
  df_out$Variant[i] <- df$V1[i]
  df_out$Count1[i] <- values[1]
  df_out$Count2[i] <- values[2]
  df_out$Count3[i] <- values[3]
  df_out$Count4[i] <- values[4]
  df_out$P_Value_Raw[i] <- p_values_raw[i]
  df_out$P_Value_Corrected[i] <- p_values_raw[i]  # Initialize corrected p-values with raw p-values
  df_out$Odds_Ratio[i] <- odds_ratio
}

# Correct p-values using Bonferroni correction
corrected_p_values <- p.adjust(p_values_raw, method ="fdr")
df_out$P_Value_Corrected <- corrected_p_values
write.table(df_out, file=outFile, sep="\t", row.names=FALSE, quote=FALSE)



# Figure to show fold change and p-value for GWAS enrichment 

outFile=paste("Fig_", inFile,".pdf",sep="")
pdf(file=outFile, width=14)

A = read.table (inFile, header=TRUE)
A1=A[order(A$Odds_Ratio),]

# Define the color range
color_range <- colorRampPalette(c("blue", "red"))(20)  # You can adjust the number of colors
colors <- color_range[cut(A1$P_Value_Corrected, breaks = 20)]

layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
par(mar = c(12, 4, 4, 2) + 0.1)
ylab="Odds ratio (Overlap GWAS with eCig DEG)"
ylim=c(0,8)
barplot(A1$Odds_Ratio, ylim = ylim, las = 1, col = colors, names.arg = A1$Variant, las=2, ylab=ylab, cex.names=0.9 )
abline(h=1, lty=2)
# Add color legend
par(mar = c(5, 0, 5, 2) + 0.1)
image(1:10, 1:10, matrix(1:100, nrow = 10), col = rev(color_range), axes = FALSE, xlab = "", ylab = "")
mtext(formatC(max(A1$P_Value_Corrected), format = "e", digits = 2), side = 3, line = -28)
mtext(formatC(min(A1$P_Value_Corrected), format = "e", digits = 2), side = 3, line = -4)

dev.off()




