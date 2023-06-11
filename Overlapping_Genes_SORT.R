library(dplyr)
library(tidyverse)

AF1 <- load_af("/Users/florin/r_project/data/matrices/Sort/Real/AF/SCD-TEST-s701-filtered")
AF2 <- load_af("/Users/florin/r_project/data/matrices/Sort/Real/AF/SCD-TEST-s702-filtered")
AF3 <- load_af("/Users/florin/r_project/data/matrices/Sort/Real/AF/SCD-TEST-s703-filtered")

KB1 <- load_kb("/Users/florin/r_project/data/matrices/Sort/Real/KB/SCD-TEST-s701-filtered")
KB2 <- load_kb("/Users/florin/r_project/data/matrices/Sort/Real/KB/SCD-TEST-s702-filtered")
KB3 <- load_kb("/Users/florin/r_project/data/matrices/Sort/Real/KB/SCD-TEST-s703-filtered")

SS1 <- load_star("/Users/florin/r_project/data/matrices/Sort/Real/SS/SCD-TEST-s701-filtered")
SS2 <- load_star("/Users/florin/r_project/data/matrices/Sort/Real/SS/SCD-TEST-s702-filtered")
SS3 <- load_star("/Users/florin/r_project/data/matrices/Sort/Real/SS/SCD-TEST-s703-filtered")


maf1 <- data.frame(as.matrix(AF1))
maf2 <- data.frame(as.matrix(AF2))
maf3 <- data.frame(as.matrix(AF3))

mkb1 <- data.frame(as.matrix(KB1))
mkb2 <- data.frame(as.matrix(KB2))
mkb3 <- data.frame(as.matrix(KB3))

mss1 <- data.frame(as.matrix(SS1))
mss2 <- data.frame(as.matrix(SS2))
mss3 <- data.frame(as.matrix(SS3))


maf1_non_zero <-  maf1 %>% mutate(count=rowSums(.!=0))
maf2_non_zero <-  maf2 %>% mutate(count=rowSums(.!=0))
maf3_non_zero <-  maf3 %>% mutate(count=rowSums(.!=0))

mkb1_non_zero <- mkb1 %>% mutate(count=rowSums(.!=0))
mkb2_non_zero <- mkb2 %>% mutate(count=rowSums(.!=0))
mkb3_non_zero <- mkb3 %>% mutate(count=rowSums(.!=0))

mss1_non_zero <- mss1 %>% mutate(count=rowSums(.!=0))
mss2_non_zero <- mss2 %>% mutate(count=rowSums(.!=0))
mss3_non_zero <- mss3 %>% mutate(count=rowSums(.!=0))


maf1_non_zero <- filter(maf1_non_zero, count > 0)
maf2_non_zero <- filter(maf2_non_zero, count > 0)
maf3_non_zero <- filter(maf3_non_zero, count > 0)

mkb1_non_zero <- filter(mkb1_non_zero, count > 0)
mkb2_non_zero <- filter(mkb2_non_zero, count > 0)
mkb3_non_zero <- filter(mkb3_non_zero, count > 0)

mss1_non_zero <- filter(mss1_non_zero, count > 0)
mss2_non_zero <- filter(mss2_non_zero, count > 0)
mss3_non_zero <- filter(mss3_non_zero, count > 0)


common_genes_af_kb1 <- intersect(row.names(maf1_non_zero), row.names(mkb1_non_zero))
common_genes_af_kb2 <- intersect(row.names(maf2_non_zero), row.names(mkb2_non_zero))
common_genes_af_kb3 <- intersect(row.names(maf3_non_zero), row.names(mkb3_non_zero))

common_genes_af_ss1 <- intersect(row.names(maf1_non_zero), row.names(mss1_non_zero))
common_genes_af_ss2 <- intersect(row.names(maf2_non_zero), row.names(mss2_non_zero))
common_genes_af_ss3 <- intersect(row.names(maf3_non_zero), row.names(mss3_non_zero))

common_genes_all1 <- intersect(common_genes_af_ss1, common_genes_af_kb1)
common_genes_all2 <- intersect(common_genes_af_ss2, common_genes_af_kb2)
common_genes_all3 <- intersect(common_genes_af_ss3, common_genes_af_kb3)


# Check which rows differ, to find out which Genes are found differently
# These values appear in m_kb but do not appear in m_af
setdiff(row.names(maf1_non_zero), row.names(mkb1_non_zero))
setdiff(row.names(maf1_non_zero), row.names(mkb1_non_zero))
setdiff(row.names(maf1_non_zero), row.names(mkb1_non_zero))

# By reversing the dataframes we will get different answers
# These values appear in m_af but do not appear in m_kb
setdiff(row.names(mkb1_non_zero), row.names(maf1_non_zero))
setdiff(row.names(mkb1_non_zero), row.names(maf1_non_zero))
setdiff(row.names(mkb1_non_zero), row.names(maf1_non_zero))




# Get common genes from the two dataframes and assign it to a new dataframe
common_genes1 <- intersect(row.names(maf1_non_zero), row.names(mkb1_non_zero))
common_genes2 <- intersect(row.names(maf2_non_zero), row.names(mkb2_non_zero))
common_genes3 <- intersect(row.names(maf3_non_zero), row.names(mkb3_non_zero))

common_genes1 <- intersect(row.names(mss1_non_zero), common_genes1)
common_genes2 <- intersect(row.names(mss2_non_zero), common_genes2)
common_genes3 <- intersect(row.names(mss3_non_zero), common_genes3)


common_genes_af1 <- maf1_non_zero[common_genes1,]
common_genes_kb1 <- mkb1_non_zero[common_genes1,]
common_genes_ss1 <- mss1_non_zero[common_genes1,]

common_genes_af2 <- maf2_non_zero[common_genes2,]
common_genes_kb2 <- mkb2_non_zero[common_genes2,]
common_genes_ss2 <- mss2_non_zero[common_genes2,]

common_genes_af3 <- maf3_non_zero[common_genes3,]
common_genes_kb3 <- mkb3_non_zero[common_genes3,]
common_genes_ss3 <- mss3_non_zero[common_genes3,]


common_genes_af1 <-common_genes_af1[291:291]
common_genes_kb1 <-common_genes_kb1[294:294]
common_genes_ss1 <-common_genes_ss1[304:304]

common_genes_af2 <-common_genes_af2[338:338]
common_genes_kb2 <-common_genes_kb2[345:345]
common_genes_ss2 <-common_genes_ss2[346:294]

common_genes_af3 <-common_genes_af3[333:333]
common_genes_kb3 <-common_genes_kb3[335:335]
common_genes_ss3 <-common_genes_ss3[340:340]


common_genes_1 <- data.frame(common_genes_kb1$count, common_genes_af1$count,common_genes_ss1$count)
common_genes_2 <- data.frame(common_genes_kb2$count, common_genes_af2$count,common_genes_ss2$count)
common_genes_3 <- data.frame(common_genes_kb3$count, common_genes_af3$count,common_genes_ss3$count)

cor(common_genes_1, method = "spearman")
cor(common_genes_2, method = "spearman")
cor(common_genes_3, method = "spearman")

cor(common_genes_1, method = "pearson")
cor(common_genes_2, method = "pearson")
cor(common_genes_3, method = "pearson")

#cor(common_genes_af1$count, common_genes_kb1$count, method = 'spearman')
#cor(common_genes_af1$count, common_genes_kb1$count, method="pearson")

#cor(common_genes_af2$count, common_genes_kb2$count, method = 'spearman')
#cor(common_genes_af2$count, common_genes_kb2$count, method="pearson")

#cor(common_genes_af3$count, common_genes_kb3$count, method = 'spearman')
#cor(common_genes_af3$count, common_genes_kb3$count, method="pearson")

# Plot 
plot(common_genes_af1$count, common_genes_kb1$count, main="Correlation between common genes af - kb", xlab="common_genes_af1", ylab="common_genes_kb1")
plot(common_genes_af2$count, common_genes_kb2$count, main="Correlation between common genes af - kb", xlab="common_genes_af2", ylab="common_genes_kb2")
plot(common_genes_af3$count, common_genes_kb3$count, main="Correlation between common genes af - kb", xlab="common_genes_af3", ylab="common_genes_kb3")

# Merged DF
merged_af_kb1 <- merge(common_genes_kb1, common_genes_af1, by='row.names', all=TRUE)
merged_af_kb2 <- merge(common_genes_kb2, common_genes_af2, by='row.names', all=TRUE)
merged_af_kb3 <- merge(common_genes_kb3, common_genes_af3, by='row.names', all=TRUE)

merged_af_ss1 <- merge(common_genes_ss1, common_genes_af1, by='row.names', all=TRUE)
merged_af_ss2 <- merge(common_genes_ss2, common_genes_af2, by='row.names', all=TRUE)
merged_af_ss3 <- merge(common_genes_ss3, common_genes_af3, by='row.names', all=TRUE)

merged_kb_ss1 <- merge(common_genes_kb1, common_genes_ss1, by='row.names', all=TRUE)
merged_kb_ss2 <- merge(common_genes_kb2, common_genes_ss2, by='row.names', all=TRUE)
merged_kb_ss3 <- merge(common_genes_kb3, common_genes_ss3, by='row.names', all=TRUE)

# ALEVIN FRY AND KALLISTO BUSTOOLS CORRELATION GRAPHS
# Plot the Merge 
ggplot(data = merged_af_kb1, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Kallisto Bustools Counts', y='Alevin-Fry Counts', title='Alevin-Fry and Kallisto Bustools Genes Count Correlation\nSCD701 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_af_kb2, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Kallisto Bustools Counts', y='Alevin-Fry Counts', title='Alevin-Fry and Kallisto Bustools Genes Count Correlation\nSCD702 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_af_kb3, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Kallisto Bustools Counts', y='Alevin-Fry Counts', title='Alevin-Fry and Kallisto Bustools Genes Count Correlation\nSCD703 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))


# ALEVIN FRY AND STARSOLO CORRELATION GRAPHS
# Plot the Merge 
ggplot(data = merged_af_ss1, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='STARsolo Counts', y='Alevin-Fry Counts', title='Alevin-Fry and STARsolo Genes Count Correlation\nSCD701 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_af_ss2, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='STARsolo Counts', y='Alevin-Fry Counts', title='Alevin-Fry and STARsolo Genes Count Correlation\nSCD702 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_af_ss3, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='STARsolo Counts', y='Alevin-Fry Counts', title='Alevin-Fry and STARsolo Genes Count Correlation\nSCD703 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))


#STARSOLO AND KALLISTO BUSTOOLS CORRELATION GRAPHS
# Plot the Merge 
ggplot(data = merged_kb_ss1, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Kallisto Bustools Counts', y='STARsolo Counts', title='STARsolo and Kallisto Bustools Genes Count Correlation\nSCD701 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_kb_ss2, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Kallisto Bustools Counts', y='STARsolo Counts', title='STARsolo and Kallisto Bustools Genes Count Correlation\nSCD702 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))

# Plot the Merge 
ggplot(data = merged_kb_ss3, mapping = aes(x=count.x, y=count.y)) + 
  geom_point(alpha = 0.1, color ='blue') +
  geom_smooth(method='lm') +
  labs(x='Kallisto Bustools Counts', y='STARsolo Counts', title='STARsolo and Kallisto Bustools Genes Count Correlation\nSCD703 Dataset') + 
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold'))
