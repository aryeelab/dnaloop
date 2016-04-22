#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("ChIA-PET set directory must be specified.\n", call.=FALSE)
}

library(ggplot2)
library(scales)
library(foreach)
library(gridExtra)
library(reshape2)
library(dplyr)

#dir <- "/PHShome/ma695/work/projects/chiapet_design/output/dnaloop-0.5.14/esc_010/"
dir <- args[1]

if (length(args)==2) {
    pdf <- args[2]
} else {
    pdf <- file.path(dir, "qc-report.pdf")    
}

message("Processing: ", dir)
message("Saving QC report to: ", pdf)

samples <- dir(dir, pattern = ".loop_counts.bedpe")
samples <- sub(".loop_counts.bedpe", "", samples)

# Creates a dataframe of summary statistics from the individual sample log output files
readstats <- foreach(sample = samples, .combine="rbind") %do% {
    sfilename <- file.path(dir, paste("log/", sample, ".read_stats.txt", sep=""))
    rs <- read.table(sfilename, header=FALSE, stringsAsFactors = FALSE)
    rs <- cbind(sample=sample, rs)
    colnames(rs) <- c("sample", "metric", "count")
    rs
}
metrics <- c("Total_PETs", "PETs_with_linker", "Mapped_PETs_q30", "Mapped_unique_PETs_q30", "Mapped_unique_intrachromosal_PETs_q30_5kb")
readstats <- subset(readstats, metric %in% metrics)

# Get loop lengths and counts
loop_pets <- foreach(sample = samples, .combine="rbind") %do% {
    sfilename <- file.path(dir, paste(sample , ".loop_counts.bedpe",sep=""))
    x <- read.table(sfilename, stringsAsFactors = FALSE)
    intra <- x[,1]==x[,4]
    x <- x[intra,]
    loop_length <- rep(x[, 5] - x[, 3], x[, 8])
    data.frame(sample = sample, loop_length = pmax(0, loop_length))
}
head(loop_pets)

# Add long-range loop counts to the read stats dataframe.
grouped <- group_by(loop_pets, sample)
summ <- summarise(grouped, long_range=sum(loop_length>=5000))
df <- data.frame(sample=as.character(summ$sample), metric="Anchor_mapped_PETs_5kb", count=summ$long_range)
readstats <- rbind(readstats, df)
metrics <- c(metrics, "Anchor_mapped_PETs_5kb")
readstats$metric <- factor(readstats$metric, levels=metrics)

# Plot read stats
p <- ggplot(readstats, aes(x = sample, y = count, fill=metric))  + 
    geom_bar(stat='identity', position=position_dodge()) + theme_bw() + 
    ggtitle("Counts of PET Types") + xlab("") + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
#p_metrics <- p +  scale_y_log10(labels=comma, breaks=10^(1:12))
p_metrics <- p + scale_y_continuous(labels=comma)

# Plot PET anchor separation distribution
p_hist <- ggplot(loop_pets, aes(loop_length)) + geom_histogram() + scale_x_log10(labels=comma, breaks=10^(3:9)) + facet_wrap(~sample, ncol=1) + theme_bw()

# Output graphics and table to PDF file
pdf(pdf, height=6, width=11, onefile = TRUE)
p_metrics
plot.new()
grid.table(format(acast(readstats, metric~sample, sum), big.mark=","))
tab <- acast(readstats, metric~sample, sum)
tab_percent <- 100*sweep(tab, 2, tab["Total_PETs",], FUN="/")
plot.new()
grid.table(format(tab_percent, digits=2, nsmall=2))
p_hist
dev.off()