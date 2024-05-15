library('foreign')
library('tidyr')
library('readxl')
library('ggplot2')
library('stats')
library('dplyr')
library('factoextra')
library('stringr')
library('umap')
library('patchwork')
library('sva')
library('tibble')
library('data.table')

#read in normalized data
normalized1 <- read.csv("Data/Data_Round1/twf-53299_belgium_proteomics_somalogic_somascan_results_20211124_normalized.csv") %>%
  filter(RowCheck == 'PASS')

normalized2 <- read.csv("data_round2/Data_Round2/twf-53299_belgium_proteomics_somalogic_somascan_results_20230322_anml_normalized.csv") %>%
  filter(RowCheck == 'PASS')

#read in annotation data
annot1 <- read.csv('Data/belgium_somalogic_annotation_20211124.csv')
annot2 <- read.csv('data_round2/Data_Round2/converted_files/somalogic_aptamer_anno.csv')

#read in sample info
sample_info1 <- read.csv('Data/belgium_somalogic_sample_info_20211124.csv')
sample_info2 <- read.csv('data_round2/Data_Round2/converted_files/somalogic_sample_info.csv')

#select annotated proteins
norm_values1 <- normalized1[, c('SampleId', annot1$SeqIndex)]
norm_values2 <- normalized2[, c('SampleId', annot2$SeqIndex)]

#batch removal
phase <- data.frame(c(rep(1, nrow(norm_values1)), rep(2, nrow(norm_values2))))
colnames(phase) <- 'phase'

common_cols <- data.frame(names = intersect(colnames(norm_values1), colnames(norm_values2)))%>%
  filter(substr(names,1,3) == 'seq')

comb <- rbind(norm_values1[, common_cols$names], norm_values2[, common_cols$names])

removal <- ComBat(t(comb), phase$phase)
removal <- t(removal)

#IMPORTANT: ONE PROTEIN FROM PHASE 2 IS REMOVED FOR THE FURTHER PROCESSING
norm_values_removal1 <- data.frame(removal[1:nrow(norm_values1), ])
norm_values_removal1$SampleId <- normalized1$SampleId
norm_values_removal2 <- data.frame(removal[(nrow(norm_values1)+1):nrow(removal), ])
norm_values_removal2$SampleId <- normalized2$SampleId


#somascan data contains visit_id as well as SampleId, we can use this to link the normalized data to the metadata for phase 1
somadata_phase_1 <- read.spss("Data/SOMAscan_data_phase1.sav", to.data.frame = TRUE)
metadata1 <- read_xlsx("Data/Sample_metadata.xlsx")
somadata_phase_1$visit_id <- paste(0, somadata_phase_1$visit_id, sep = "")
somadata_ids1 <- somadata_phase_1[ ,c('visit_id', 'SampleId')]

#use MATCHED file to link sampleids from second phase with visit_ids (which are the same in phase 1 and 2)
matched2 <- read.csv("./data_round2/Data_Round2/converted_files/MATCHEDround1_somalogic_customid_20230322.csv") %>%
  mutate(visit_id = substr(SubjectID, 4, 11))
metadata2 <- read_xlsx("data_round2/Data_Round2/Sample_metadata_phase2.xlsx")


link_ids_round2 <- merge(somadata_phase_1[ ,c('visit_id','SampleId')], matched2[ ,c('visit_id','SampleId')], by = 'visit_id') %>%
  mutate(SampleId = SampleId.y) %>%
  select(-SampleId.x, -SampleId.y)


#normalized protein intensities with visit_id and SampleId
norm_values_removal1 <- merge(somadata_ids1, norm_values_removal1, by = 'SampleId')
norm_values_removal2 <- merge(link_ids_round2, norm_values_removal2, by = 'SampleId')

#normalized values merged with metadata
norm_meta1 <- merge(norm_values_removal1, metadata1, by = 'visit_id')
sum(norm_meta1$Gender2 == 'Male')
sum(norm_meta1$Gender2 == 'Female')
norm_meta1$Gender <- factor(norm_meta1$Gender, levels = c('F', 'M'), labels = c('Female', 'Male'))

norm_meta2 <- merge(norm_values_removal2, metadata2, by = 'visit_id')
sum(norm_meta2$Sex == 'male')
sum(norm_meta2$Sex == 'female')
norm_meta2$Gender <- factor(norm_meta2$Sex, levels = c('female', 'male'), labels = c('Female', 'Male'))
norm_meta2$Age_R2_years = norm_meta2$`Age R2 (year)`


#-----------------------------------------------------------------------------------------------------------
#plotting intensities phase 1 and phase 2 separately
intensitydata1 <- norm_meta1 %>% select(starts_with(('seq')))
intensitydata1_batch_removed <- norm_meta1 %>% 
  filter(!str_sub(VisitB_Date, 1, 7) %in% c("2003-04", "2003-05", "2003-06", "2003-07")) %>%
  select(starts_with('seq'))
intensitydata2 <- norm_meta2 %>% select(starts_with(('seq')))

#transpose data frame
intensities1.T <- transpose(intensitydata1)
rownames(intensities1.T) <- colnames(intensitydata1)
colnames(intensities1.T) <- rownames(intensitydata1)

w.plot1 <- melt(intensities1.T)  %>%
  mutate(value = log2(value))

ggplot(aes(x=value, color = variable), data=w.plot1) +
  geom_density() + 
  theme(legend.position = "none") +
  labs(x = 'log2(intensities)', y = 'Density', title = 'Distribution of intensities for phase 1')

#transpose data frame
intensities1_batch_removed.T <- transpose(intensitydata1_batch_removed)
rownames(intensities1_batch_removed.T) <- colnames(intensitydata1_batch_removed)
colnames(intensities1_batch_removed.T) <- rownames(intensitydata1_batch_removed)

w.plot1_batch_removed <- melt(intensities1_batch_removed.T)  %>%
  mutate(value = log2(value))

ggplot(aes(x=value, color = variable), data=w.plot1_batch_removed) +
  geom_density() + 
  theme(legend.position = "none") +
  labs(x = 'log2(intensities)', y = 'Density', title = 'Distribution of intensities for phase 1 (batches removed)')


#transpose data frame
intensities2.T <- transpose(intensitydata2)
rownames(intensities2.T) <- colnames(intensitydata2)
colnames(intensities2.T) <- rownames(intensitydata2)

w.plot2 <- melt(intensities2.T)  %>%
  mutate(value = log2(value))

p <- ggplot(aes(x=value, color = variable), data=w.plot2) +
  geom_density() + 
  theme(legend.position = "none") +
  labs(x = 'log2(intensities)', y = 'Density', title = 'Distribution of intensities for phase 2')

density_data <- ggplot_build(p)$data[[1]]
density_data[density_data$density == max(density_data$density), ]


intensities2_removed_person.T <- transpose(intensitydata2[-186 ,])
rownames(intensities2_removed_person.T) <- colnames(intensitydata2[-186 ,])
colnames(intensities2_removed_person.T) <- rownames(intensitydata2[-186 ,])

w.plot2_removed_person <- melt(intensities2_removed_person.T)  %>%
  mutate(value = log2(value))

ggplot(aes(x=value, color = variable), data=w.plot2_removed_person) +
  geom_density() + 
  theme(legend.position = "none") +
  labs(x = 'log2(intensities)', y = 'Density', title = 'Distribution of intensities for phase 2 (without outlying person)')

#------------------------------------------------------------------------------------------------------------
#plotting intensities for both phases combined
w.plot2$variable <- as.factor(as.numeric(w.plot2$variable) + max(as.numeric(w.plot1$variable)))
w.plot1$phase <- "phase_1"
w.plot2$phase <- "phase_2"
comb <- rbind(w.plot1, w.plot2)

ggplot(data = comb, aes(x = value, group = variable)) +
  geom_density(aes(color = phase)) +
  labs(x = 'log2(intensities)', y = 'Density', title = 'Distribution of intensities for both phases')


w.plot2$variable <- as.factor(as.numeric(w.plot2$variable) + max(as.numeric(w.plot1_batch_removed$variable)))
w.plot1_batch_removed$phase <- "phase_1"
w.plot2$phase <- "phase_2"
comb2 <- rbind(w.plot1_batch_removed, w.plot2)

ggplot(data = comb2, aes(x = value, group = variable)) +
  geom_density(aes(color = phase)) +
  labs(x = 'log2(intensities)', y = 'Density', title = 'Distribution of intensities for both phases (batch effect from phase 1 removed)')
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#PCA for phase 1
pca_data1 <- norm_meta1 %>%
  select(c(VisitB_Date, starts_with('seq'))) %>%
  na.omit(pca_data1) %>%
  mutate(time = as.factor(substr(VisitB_Date, start = 1, stop= 7))) %>%
  select(-VisitB_Date)

pca1 <- prcomp(pca_data1 %>% select(-time), center = TRUE, scale = TRUE)
summary1 <- summary(pca1)
variances1 <- summary1$importance[2,1:2]

plotdata1_1 <- as.data.frame(bind_cols(pca1$x[,1:2], pca_data1$time)) %>%
  mutate(time = ...3) %>% select(-...3)
ggplot(plotdata1_1, aes(x = PC1, y = PC2, color = time)) +
  labs(x = paste0('PC1. Explained variance: ', round(variances1[1],4)*100,'%'), 
       y = paste0('PC2. Explained variance: ', round(variances1[2],4)*100, '%')) +
  ggtitle("First round - check for batch effects") +
  geom_point()+
  scale_x_continuous(limits = c(-225, 50)) +
  scale_y_continuous(limits = c(-40, 130))

plotdata1_2 <- aggregate(cbind(PC1, PC2) ~ time, data = plotdata1_1, FUN = mean)
ggplot(plotdata1_2, aes(x = PC1, y = PC2, color = time)) +
  labs(x = paste0('PC1. Explained variance: ', round(variances1[1],4)*100,'%'), 
       y = paste0('PC2. Explained variance: ', round(variances1[2],4)*100, '%')) +
  ggtitle("First round - check for batch effects") +
  geom_point() + geom_label(aes(label = time)) +
  scale_x_continuous(limits = c(-225, 50)) +
  scale_y_continuous(limits = c(-40, 130)) +
  theme(legend.position = 'none')

ggplot() +
  geom_point(data = plotdata1_1, aes(x = PC1, y = PC2, color = time)) +
  geom_text(data = plotdata1_2, aes(x = PC1, y = PC2, label = time), vjust = -0.5) +
  scale_x_continuous(limits = c(-225, 50)) +
  scale_y_continuous(limits = c(-40, 130)) +
  labs(x = paste0('PC1. Explained variance: ', round(variances1[1],4)*100,'%'), 
       y = paste0('PC2. Explained variance: ', round(variances1[2],4)*100, '%')) +
  ggtitle("First round - check for batch effects")

#----------------------------------------------------------------------------------------------------
#PCA for phase 2
pca_data2 <- norm_meta2 %>%
  select(c(`Examination date R2`, starts_with('seq'))) %>%
  na.omit(pca_data2) %>%
  mutate(time = as.factor(substr(`Examination date R2`, start = 1, stop= 7))) %>%
  select(-`Examination date R2`)

pca2 <- prcomp(pca_data2 %>% select(-time), center = TRUE, scale = TRUE)
summary2 <- summary(pca2)
variances2 <- summary2$importance[2,1:2]

plotdata2_1 <- as.data.frame(bind_cols(pca2$x[,1:2], pca_data2$time)) %>%
  mutate(time = ...3) %>% select(-...3)

ggplot(plotdata2_1, aes(x = PC1, y = PC2, color = time)) +
  labs(x = paste0('PC1. Explained variance: ', round(variances2[1],4)*100,'%'), 
       y = paste0('PC2. Explained variance: ', round(variances2[2],4)*100, '%')) +
  ggtitle("Second round - check for batch effects") +
  geom_point() +
  scale_x_continuous(limits = c(-50, 150)) +
  scale_y_continuous(limits = c(-100, 50))

plotdata2_2 <- aggregate(cbind(PC1, PC2) ~ time, data = plotdata2_1, FUN = mean)

ggplot(plotdata2_2, aes(x = PC1, y = PC2, color = time)) +
  labs(x = paste0('PC1. Explained variance: ', round(variances2[1],4)*100,'%'), 
       y = paste0('PC2. Explained variance: ', round(variances2[2],4)*100, '%')) +
  ggtitle("Second round - check for batch effects") +
  geom_point() + geom_label(aes(label = time))+
  scale_x_continuous(limits = c(-50, 150)) +
  scale_y_continuous(limits = c(-100, 50)) +
  theme(legend.position = 'none')

ggplot() +
  geom_point(data = plotdata2_1, aes(x = PC1, y = PC2, color = time)) +
  geom_text(data = plotdata2_2, aes(x = PC1, y = PC2, label = time), vjust = -0.5) +
  scale_x_continuous(limits = c(-50, 150)) +
  scale_y_continuous(limits = c(-100, 50)) +
  labs(x = paste0('PC1. Explained variance: ', round(variances2[1],4)*100,'%'), 
       y = paste0('PC2. Explained variance: ', round(variances2[2],4)*100, '%')) +
  ggtitle("Second round - check for batch effects")


#phase1+phase2 plots-----------------------------------------------------------------------------------------------------------
pcadata1_and_2 <- bind_rows(pca_data1, pca_data2)

pca1_and_2 <- prcomp(pcadata1_and_2 %>% select(-time), center = TRUE, scale = TRUE) #leave out last protein in dataframe because there are no measurements for this one in phase1
summary1_and_2 <- summary(pca1_and_2)
variances1_and_2 <- summary1_and_2$importance[2,1:2]

plotdata1_and_2 <- as.data.frame(bind_cols(pca1_and_2$x[,1:2], pcadata1_and_2$time)) %>%
  mutate(time = ...3) %>% select(-...3) %>% mutate(month = substr(time, 6,7))
plotdata1_and_2$phase <- factor(as.numeric(substr(plotdata1_and_2$time, 1, 4)) > 2004, levels = c(FALSE, TRUE), labels = c('Phase1', 'Phase2'))

ggplot(plotdata1_and_2, aes(x = PC1, y = PC2, color = phase)) +
  labs(x = paste0('PC1. Explained variance: ', round(variances1_and_2[1],4)*100,'%'), 
       y = paste0('PC2. Explained variance: ', round(variances1_and_2[2],4)*100, '%')) +
  ggtitle("Both rounds - check for batch effects") +
  geom_point() +
  scale_x_continuous(limits = c(-200, 50)) +
  scale_y_continuous(limits = c(-50, 125))


plotdata_means <- aggregate(cbind(PC1, PC2) ~ month, data = plotdata1_and_2, FUN = mean) %>%
  mutate(month = factor(month, levels = month, labels= c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')))

ggplot(plotdata_means, aes(x = PC1, y = PC2, color = month)) +
  labs(x = paste0('PC1. Explained variance: ', round(variances1_and_2[1],4)*100,'%'), 
       y = paste0('PC2. Explained variance: ', round(variances1_and_2[2],4)*100, '%')) +
  ggtitle("Both rounds - check for batch effects") +
  geom_point() + geom_label(aes(label = month), show.legend = F)+
  scale_x_continuous(limits = c(-200, 50)) +
  scale_y_continuous(limits = c(-50, 125))

plotdata_means_sep <- aggregate(cbind(PC1, PC2) ~ time, data = plotdata1_and_2, FUN = mean)
ggplot() +
  geom_point(data = plotdata1_and_2, aes(x = PC1, y = PC2, color = phase)) +
  geom_text(data = plotdata_means_sep, aes(x = PC1, y = PC2, label = time), vjust = -0.5) +
  scale_x_continuous(limits = c(-200, 50)) +
  scale_y_continuous(limits = c(-50, 125)) +
  labs(x = paste0('PC1. Explained variance: ', round(variances1_and_2[1],4)*100,'%'), 
       y = paste0('PC2. Explained variance: ', round(variances1_and_2[2],4)*100, '%')) +
  ggtitle("Both rounds - check for batch effects")
