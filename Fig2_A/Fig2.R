library(RColorBrewer)
library(ComplexHeatmap)
library(ggpubr)
library(RColorBrewer)
library(circlize)

# Load data.
mut <- read.table('1_Mutations_per_Mb.txt', header = TRUE)
drivers <- t(as.matrix(read.table('2_Drivers_table.txt', header = TRUE, sep = '\t', fill = TRUE)))
pam <- t(as.matrix(read.table('PAM50_sample_JBH.txt', header = TRUE, fill = TRUE)[1:2]))
columns <- sort(unique(c(as.character(mut$NEW_timepoint), drivers[1, ], pam[1, ])))

# Load PAM50 and create order index.
pam[pam == 'N006'] <- NA
all.pam <- matrix(NA, nrow = nrow(pam) - 1, ncol = length(columns))
all.pam[, columns %in% pam[1, ]] <- matrix(pam[-1, ], nrow = 1)
rownames(all.pam) <- 'PAM50'
colnames(all.pam) <- columns

# Make patient separators.
pat <- substr(columns, 1, 4)
starts <- sapply(lapply(lapply(unique(pat) , '==', pat), which), min)
ends <- sapply(lapply(lapply(unique(pat) , '==', pat), which), max)

# Compute order(all.pam[, starts]) and manually modify the 2 exceptions.
idx.pat <- c(3, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36)
idx.sam <- unlist(Map(seq, starts, ends))

#PDF_FirstPart

all.mut <- rep(NA, nrow(mut))
all.mut[columns %in% mut$NEW_timepoint] <- mut$MutMb
names(all.mut) <- columns
color <- c("#FFCD00B2", 'chartreuse3',"#CC0C00B2", "black")[as.factor(substr(columns, 6, 7))]

pdf('./1_Plot1_mutMb_VISITS_Colours.pdf', height = 6, width = 16)
barplot(all.mut[idx.sam], ylab = 'Mutations per Mb', col = color, las = 2, cex.names = 0.7, 
        space=c(0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,0,1,0,1,0,1,0,0,0,1,
                0,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,1))
legend("topright", 
       legend = c("V1", "V2", "V3", "VR"), 
       fill = c("#FFCD00B2", 'chartreuse3',"#CC0C00B2", "black"))
dev.off()

# Mid panel.
# Start with drivers.
drivers[drivers == 'INFRAME;INFRAME'] <- 'INFRAME'
drivers[drivers == 'MISSENSE;MISSENSE'] <- 'MISSENSE'
drivers[drivers == ''] <- NA
idx.dri <- order(rowSums(!is.na(drivers[-1, ])))
all.dri <- matrix(NA, nrow = nrow(drivers) - 1, ncol = length(columns))
all.dri[, columns %in% drivers[1, ]] <- drivers[-1, ][idx.dri, ]
rownames(all.dri) <- rownames(drivers)[-1][idx.dri]
colnames(all.dri) <- columns

# Merge drivers and PAM50 data.
all.dat <- rbind(all.dri[nrow(all.dri):1, ], all.pam)

pdf('./2_Plot1_drivers_pam50.pdf', height = 4.6, width = 16)
hml <- Heatmap(all.dat[, starts[1]:ends[1]], cluster_rows = FALSE, cluster_columns = FALSE, split = rep(1:2, times = c(nrow(all.dri), nrow(all.pam))),
               col = c('INFRAME' = 'darkgrey', 'MISSENSE' = 'forestgreen', 'TRUNC' = 'black',
                       'BasalLike' = 'red', 'Her2Enriched' = 'pink', 'LuminalA' = 'darkblue', 'LuminalB' = 'blue'), na_col = 'beige',
               rect_gp = gpar(col = 'white'), show_row_names = FALSE, column_names_gp = gpar(fontsize = 7), show_heatmap_legend = FALSE)
for (i in 2:(length(starts) - 1)) {
    hml <- hml + Heatmap(all.dat[, starts[i]:ends[i]], cluster_rows = FALSE, cluster_columns = FALSE, split = rep(1:2, times = c(nrow(all.dri), nrow(all.pam))),
                         col = c('INFRAME' = 'darkgrey', 'MISSENSE' = 'forestgreen', 'TRUNC' = 'black',
                                 'BasalLike' = 'red', 'Her2Enriched' = 'pink', 'LuminalA' = 'darkblue', 'LuminalB' = 'blue'), na_col = 'beige',
                         rect_gp = gpar(col = 'white'), show_row_names = FALSE, column_names_gp = gpar(fontsize = 7), show_heatmap_legend = FALSE)
}
hml <- hml + Heatmap(all.dat[, starts[length(starts)]:ends[length(starts)]], cluster_rows = FALSE, cluster_columns = FALSE, split = rep(1:2, times = c(nrow(all.dri), nrow(all.pam))),
                     col = c('INFRAME' = 'darkgrey', 'MISSENSE' = 'forestgreen', 'TRUNC' = 'black',
                             'BasalLike' = 'red', 'Her2Enriched' = 'pink', 'LuminalA' = 'darkblue', 'LuminalB' = 'blue'), na_col = 'beige',
                     rect_gp = gpar(col = 'white'), row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7), show_heatmap_legend = FALSE)
hml
dev.off()

#Purity panel.

DF <- read.csv("4_Purity.txt", sep= ";", row.names = 1)
mat <- as.matrix(DF)
col_fun = colorRamp2(c(0, 0.5, 1), c("plum1", "purple1", "purple4"))
hm = Heatmap(mat, cluster_rows = FALSE, col = col_fun, row_names_side = "left")
pdf('./3_Purity_EDIT_Filtered.pdf', height = 12, width = 12)
ht = grid.grabExpr(draw(hm))
grid.newpage()
pushViewport(viewport(angle = 90))
grid.draw(ht)
popViewport()
dev.off()

# Bottom panel.
# % signature (merged).

sig <- read.table('5_Signature.txt')
sig.idx <- c(1, 2, 3, 6, 9, 10, 13, 15, 30)

sig.mer <- data.frame('Aging' = sig[, 1])
sig.mer <- cbind(sig.mer, data.frame('APOBEC' = rowSums(sig[, c(2, 13)])))
sig.mer <- cbind(sig.mer, data.frame('BRCA' = sig[, 3]))
sig.mer <- cbind(sig.mer, data.frame('Defective DNA repair' = rowSums(sig[, c(6, 15)])))
sig.mer <- cbind(sig.mer, data.frame('Polymerase n' = sig[, 9]))
sig.mer <- cbind(sig.mer, data.frame('POLE mutations' = sig[, 10]))
sig.mer <- cbind(sig.mer, data.frame('Unknown breast' = sig[, 30]))
sig.mer <- cbind(sig.mer, data.frame('Others' = rowSums(sig[, -c(sig.idx, ncol(sig))])))
all.mer <- matrix(0, nrow = ncol(sig.mer), ncol = length(columns))
all.mer[, columns %in% sub('_PT', '', rownames(sig.mer))] <- t(sig.mer)
all.mer <- rbind(all.mer, 0)
rownames(all.mer) <- c(colnames(sig.mer), 'empty')
colnames(all.mer) <- columns
all.mer['empty', apply(all.mer == 0, 2, all)] <- 1

pdf('figures_JBH/Plot1_signatures_merged.pdf', height = 3, width = 16)
par(mar = c(5, 3, 2, 11) + 0.1, xpd = TRUE)
barplot(all.mer[, idx.sam], col = c("#C0C0C0", '#BEBADA', '#FF7F00', '#80B1D3', '#FDB462',
                                   '#B3DE69', '#FCCDE5', '#BC80BD'), yaxt = 'n', las = 2, cex.names = 0.5,
       space=c(0,0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,1,0,0,1,0,1,0,1,0,0,0,1,
                0,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,1))
legend(x = 93, y = 1.1,
       legend = c('Clock-like, age related (1, 5)', 'APOBEC (2, 13)', 'BRCA (3)', 'Defective DNA repair (6, 15)', 'Polymerase n (9)', 'POLE mutations (10)', 'Unknown breast (30)', 'Others'),
       fill = c("#C0C0C0", '#BEBADA', '#FF7F00', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD'), bty = 'n')
par(mar = c(5, 4, 4, 2) + 0.1)
dev.off()