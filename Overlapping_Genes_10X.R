library(dplyr)
library(tidyverse)

AF <- load_af("/Users/florin/r_project/data/matrices/Tenx/Real/AF/SCD-WP-g001-filtered")

KB <- load_kb("/Users/florin/r_project/data/matrices/Tenx/Real/KB/SCD-WP-g001-filtered")

SS <- load_star("/Users/florin/r_project/data/matrices/Tenx/Real/SS1")

CR <- load_star("/Users/florin/r_project/data/matrices/Tenx/Real/CR")


AF <- data.frame(as.matrix(AF))

KB <- data.frame(as.matrix(KB))

SS <- data.frame(as.matrix(SS))

CR <- data.frame(as.matrix(CR))


AF <- AF %>% mutate(count=rowSums(.!=0))

KB <- KB %>% mutate(count=rowSums(.!=0))

SS <- CR %>% mutate(count=rowSums(.!=0))

CR <- CR %>% mutate(count=rowSums(.!=0))


AF <- filter(AF, count > 0)

KB <- filter(KB, count > 0)

SS <- filter(SS, count > 0)

CR <- filter(CR, count > 0)


common_genes_af_kb <- intersect(row.names(AF), row.names(KB))

common_genes_af_ss <- intersect(row.names(AF), row.names(SS))

common_genes_af_cr <- intersect(row.names(AF), row.names(CR))


# Check which rows differ, to find out which Genes are found differently
# These values appear in m_kb but do not appear in m_af
setdiff(row.names(maf_non_zero), row.names(mkb_non_zero))


# By reversing the dataframes we will get different answers
# These values appear in m_af but do not appear in m_kb
setdiff(row.names(mkb_non_zero), row.names(maf_non_zero))




# Get common genes from the two dataframes and assign it to a new dataframe
common_genes <- intersect(row.names(AF), row.names(KB))
common_genes <- intersect(row.names(SS), common_genes)
common_genes <- intersect(row.names(CR), common_genes)

common_genes_af <- AF[common_genes,]
common_genes_kb <- KB[common_genes,]
common_genes_ss <- SS[common_genes,]
common_genes_cr <- CR[common_genes,]


common_genes_af <- common_genes_af[8804:8804]
common_genes_kb <- common_genes_kb[9215:9215]
common_genes_ss <- common_genes_ss[9423:9423]
common_genes_cr <- common_genes_cr[9423:9423]

common_genes_all <- data.frame(common_genes_af$count, common_genes_kb$count, common_genes_ss$count, common_genes_cr$count)


cor(common_genes_all, method = "spearman")


cor(common_genes_all, method = "pearson")

#cor(common_genes_af1$count, common_genes_kb1$count, method = 'spearman')
#cor(common_genes_af1$count, common_genes_kb1$count, method="pearson")

#cor(common_genes_af2$count, common_genes_kb2$count, method = 'spearman')
#cor(common_genes_af2$count, common_genes_kb2$count, method="pearson")

#cor(common_genes_af3$count, common_genes_kb3$count, method = 'spearman')
#cor(common_genes_af3$count, common_genes_kb3$count, method="pearson")

# Plot 
#plot(common_genes_af1$count, common_genes_kb1$count, main="Correlation between common genes af - kb", xlab="common_genes_af1", ylab="common_genes_kb1")
#plot(common_genes_af2$count, common_genes_kb2$count, main="Correlation between common genes af - kb", xlab="common_genes_af2", ylab="common_genes_kb2")
#plot(common_genes_af3$count, common_genes_kb3$count, main="Correlation between common genes af - kb", xlab="common_genes_af3", ylab="common_genes_kb3")

# Merged DF
merged_af_kb <- merge(common_genes_af, common_genes_kb, by='row.names', all=TRUE)
merged_af_ss <- merge(common_genes_af, common_genes_ss, by='row.names', all=TRUE)
merged_af_cr <- merge(common_genes_af, common_genes_cr, by='row.names', all=TRUE)
merged_kb_ss <- merge(common_genes_kb, common_genes_ss, by='row.names', all=TRUE)
merged_kb_cr <- merge(common_genes_kb, common_genes_cr, by='row.names', all=TRUE)
merged_ss_cr <- merge(common_genes_ss, common_genes_cr, by='row.names', all=TRUE)

# ALEVIN FRY AND KALLISTO BUSTOOLS CORRELATION GRAPHS
# Plot the Merge 
ggplot(data = merged_af_kb, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Alevin-Fry Counts', y='Kallisto Bustools Counts', title='Alevin-Fry and Kallisto Bustools Genes Count Correlation\nSCD-WP-g001 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_af_ss, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Alevin-Fry Counts', y='STARsolo Counts', title='Alevin-Fry and STARsolo Genes Count Correlation\nSCD-WP-g001 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_af_cr, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Alevin-Fry Counts', y='Cell Ranger Counts', title='Alevin-Fry and Cell Ranger Genes Count Correlation\nSCD-WP-g001 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_kb_ss, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Kallisto Bustools Counts', y='STARsolo Counts', title='Kallisto Bustools and STARsolo Genes Count Correlation\nSCD-WP-g001 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_kb_cr, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Kallisto Bustools Counts', y='Cell Ranger Counts', title='Kallisto Bustools and Cell Ranger Genes Count Correlation\nSCD-WP-g001 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_ss_cr, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='STARsolo Counts', y='Cell Ranger Counts', title='STARsolo and Cell Ranger Genes Count Correlation\nSCD-WP-g001 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))
