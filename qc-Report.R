library(ggplot2)
library(foreach)
library(gridExtra)
library(reshape2)

dir <- "/Users/lareauc/Desktop/esc_jurkat/"

samples <- dir(dir, pattern = ".loop_counts.bedpe")
samples <- sub(".loop_counts.bedpe", "", samples)

# Creates a list of a numeric matrix for each sample
# Matrix has 2 columns... distance between anchors and number of counts 
dist_counts <- foreach(sample = samples) %do% {
    sfilename <- paste(dir, sample , ".loop_counts.bedpe",sep="")
    awkcmd <- paste("awk '$1 == $4 {print $5-$3 \" \" $8}' " , sfilename, sep = "")
    awkout <- system(awkcmd, intern=TRUE)
    d <- matrix(as.numeric(matrix(t(as.data.frame(strsplit(awkout, " "))))),ncol=2)
    colnames(d) <- c("dist", "counts")
    d
}

# Create binned matrix
binned <- sapply(seq(1,length(dist_counts)), function(j){
    temp <- dist_counts[[j]]
    df <- as.data.frame(temp)
    vals <- c(-Inf, 0, 1000, 2000, 3000, 4000, 5000, Inf)
    sapply(seq(1, length(vals)-1), function(i){
        with(df, sum(df[dist >= vals[i] & dist < vals[i+1], "counts"]))
    })
})
colnames(binned) <- samples
rownames(binned) <- c("Self", "0-1kb", "1-2kb", "2-3kb", "3-4kb", "4-5kb", ">5kb")

# Make ggplot object
m1a <- melt(binned)
colnames(m1a) <- c("LoopType", "Sample", "Counts")
p1a <- ggplot(m1a, aes(x = Sample, y = Counts, fill=LoopType))  + 
    geom_bar(stat='identity', position=position_dodge()) + theme_bw() + 
    ggtitle("PETs with Distance between Anchors") + xlab("") + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 



# Creates a dataframe of the summary statistics from the log output files
readstatsraw <- foreach(sample = samples) %do% {
    sfilename <- paste(dir, "samples/", sample, "/read_stats.txt", sep="")
    rs <- read.table(sfilename, header=TRUE)
    colnames(rs) <- c("description", sample)
    rs
}
readstats <- Reduce(function(...) merge(..., all=T), readstatsraw)
rd <- readstats[,-1]
rownames(rd) <- readstats[,1]
readstats <- t(rd)

# Make data frame for second plot
rs2 <- cbind(readstats[,c(15,5,8,14,13)],binned[7,])
colnames(rs2) <- c("Total", "withLinker", "Mapped", "Unique", "Intrachromosomal", ">5kb")
rs2

# Make Second Plot
m2a <- melt(rs2)
colnames(m2a) <- c("Sample", "Type", "Counts")
p2a <- ggplot(m2a, aes(x = Sample, y = Counts, fill=Type))  + 
    geom_bar(stat='identity', position=position_dodge()) + theme_bw() + 
    ggtitle("Counts of PET Types") + xlab("") + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 

# Output graphics and table to PDF file
pdf("qc-Report.pdf", height=8.5, width=11, onefile = TRUE)

#First plot and table
p2a
plot.new()
grid.table(rs2)

#Second plot and table
p1a
Total <- colSums(binned)
b<- rbind(binned, Total)
plot.new()
grid.table(t(b))

dev.off()