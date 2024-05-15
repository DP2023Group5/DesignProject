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
norm_values1 <- merge(somadata_ids1, norm_values1, by = 'SampleId')
norm_values2 <- merge(link_ids_round2, norm_values2, by = 'SampleId')

#normalized values merged with metadata
norm_meta1 <- merge(norm_values1, metadata1, by = 'visit_id')
sum(norm_meta1$Gender2 == 'Male')
sum(norm_meta1$Gender2 == 'Female')
norm_meta1$Gender <- factor(norm_meta1$Gender, levels = c('F', 'M'), labels = c('Female', 'Male'))

norm_meta2 <- merge(norm_values2, metadata2, by = 'visit_id')
sum(norm_meta2$Sex == 'male')
sum(norm_meta2$Sex == 'female')
norm_meta2$Gender <- factor(norm_meta2$Sex, levels = c('female', 'male'), labels = c('Female', 'Male'))
norm_meta2$Age_R2_years = norm_meta2$`Age R2 (year)`



ageplot1 <-ggplot(norm_meta1, aes(x = Age_years, fill = Gender)) +
  geom_histogram(binwidth  =3, alpha = 0.4, position = 'identity', color = 'black') +
  ggtitle('Age distribution')

ageplot2 <-ggplot(norm_meta2, aes(x = Age_R2_years, fill = Gender)) +
  geom_histogram(binwidth  =3, alpha = 0.4, position = 'identity', color = 'black') +
  ggtitle('Age distribution')
# 
# 
# nicplot <- ggplot(norm_meta, aes(x = Nicotin_use_any_nicotin)) + geom_bar(aes(fill = Gender)) + 
#   ggtitle('Nicotin usage') + xlab('')
# 
# norm_meta$BMI_class <- factor(norm_meta$BMI_class, 
#                               levels = c('BMI<20 kg/m2', 'BMI 20-25 kg/m2', 'BMI 25-27 kg/m2', 'BMI 27-30 kg/m2', 'BMI 30-35 kg/m2', 'BMI >35 kg/m2'))
# BMIplot<- ggplot(norm_meta, aes(x = BMI_class)) + geom_bar(aes(fill = Gender)) + ggtitle('BMI')
# 
# IMTdata <- norm_meta %>%
#   select(c(IMTCLMAX, IMTCRMAX, IMTFLMAX, IMTFRMAX, Gender)) %>%
#   mutate(IMTCLMAX = ifelse(is.na(IMTCLMAX), 0, IMTCLMAX),
#          IMTCRMAX = ifelse(is.na(IMTCRMAX), 0, IMTCRMAX),
#          IMTFLMAX = ifelse(is.na(IMTFLMAX), 0, IMTFLMAX),
#          IMTFRMAX = ifelse(is.na(IMTFRMAX), 0, IMTFRMAX)) %>%
#   mutate(maxIMT = pmax(IMTCLMAX, IMTCRMAX, IMTFLMAX, IMTFRMAX))
# 
# IMTplot <- ggplot(IMTdata, aes(x = maxIMT, fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black') + 
#   ggtitle('Maximum IMT') + xlab('MaxIMT (mm)')
# 
# logIMT <- IMTdata %>%
#   mutate(maxIMT = log2(maxIMT))
# logIMTplot <- ggplot(logIMT, aes(x = maxIMT, fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black') +
#   ggtitle('Log Maximum IMT') +  xlab('Log(MaxIMT)')
# 
# alcohol <- ggplot(norm_meta, aes(x = Alcohol2)) + geom_bar(aes(fill = Gender)) + ggtitle('Alcohol usage') +
#   xlab('')
# 
# norm_meta$Glycemic_state <- factor(norm_meta$Glycemic_state, 
#                                    levels = c('Euglycemic', 'IFG >= 100 mg/dl', 'Known type 2 diabetes', 'Discovered type 2 diabetes'))
# glycemia <- ggplot(norm_meta, aes(x = factor(Glycemic_state))) + geom_bar(aes(fill = Gender)) +
#   xlab(label ='') + ggtitle(label = 'Glycemic state')
# 
# 
# blood1 <- ggplot(norm_meta, aes(x = Systolic_BP_mmHg, fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# blood2 <- ggplot(norm_meta, aes(x = Diastolic_BP_mmHg, fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# blood3 <- ggplot(norm_meta, aes(x = Pulse_pressure_mmHg, fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# 
# bpplot <- blood1+ ggtitle(label = 'Blood pressure') + blood2 + blood3 + plot_layout(nrow = 3) 
# 
# lipid1 <- ggplot(norm_meta, aes(x = Total_cholesterol_mgdl, fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# lipid2 <- ggplot(norm_meta, aes(x = HDLcholesterol_mgdl, fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# lipid3 <- ggplot(norm_meta, aes(x = LDLcholesterol_mgdl, fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# lipid4 <- ggplot(norm_meta, aes(x = log2(Triglycerides_mgdl), fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# 
# lipidplot <- lipid1 + ggtitle('Lipid levels') + lipid2 + lipid3 + lipid4
# 
# stress1 <- ggplot(norm_meta, aes(x = log2(Highsensitive_CRP_mgl), fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# stress2 <- ggplot(norm_meta, aes(x = log2(Interleukin6_pgml), fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# stress3 <- ggplot(norm_meta, aes(x = Whitebloodcellcount10e3microl, fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# stress4 <- ggplot(norm_meta, aes(x = Fibrinogen_mgdl, fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# stress5 <- ggplot(norm_meta, aes(x = log2(Homocystein_micromoll), fill = Gender)) + 
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# 
# stress <- stress1 + ggtitle('Inflammatory and oxidative stress markers') + stress2 + stress3 + stress4 + stress5 + plot_layout(nrow = 3, ncol = 2)
# 
# renal1 <- ggplot(norm_meta, aes(x = Creatinin_mgdl, fill = Gender)) +
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# renal2 <- ggplot(norm_meta, aes(x = Uric_acid_mgdl, fill = Gender)) +
#   geom_histogram(bins = 30, alpha = 0.4, position = 'identity', color = 'black')
# 
# renal <- renal1 +ggtitle('Renal function') + renal2 
# 
# 
# #------------------------------------------------------------------------------------------------------------------
# #density plot of intensities 
# #(proteins still need to be factorized and plots need to be filled by proteins,
# #but this won't work because computationally very expensive?)
# # intensity_plotdata <- norm_values %>%
# #   pivot_longer(cols = starts_with('seq'), names_to = 'proteins', values_to = 'intensities') %>%
# #   select(c('proteins', 'intensities')) %>%
# #   mutate(proteins = factor(proteins))
# # 
# # intensity_plot <- ggplot(intensity_plotdata, aes(x = intensities, color = proteins)) + geom_density()
# # intensity_plot
# # 
# # logintensity_plotdata <- intensity_plotdata %>%
# #   mutate(intensities = log2(intensities))
# # logintensity_plot <- ggplot(logintensity_plotdata, aes(x = intensities, color = proteins)) + geom_density()
# # logintensity_plot
# 

#PCA for phase 1
pca_data1 <- norm_meta1 %>%
  select(c(VisitB_Date, starts_with('seq'))) %>%
  na.omit(pca_data1) %>%
  mutate(time = as.factor(substr(VisitB_Date, start = 1, stop= 7))) %>%
  select(-VisitB_Date)

pca1 <- prcomp(pca_data1[ , 1:7523], center = TRUE, scale = TRUE)
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

pca2 <- prcomp(pca_data2[ , 1:7524], center = TRUE, scale = TRUE)
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

pca1_and_2 <- prcomp(pcadata1_and_2[ , 1:7523], center = TRUE, scale = TRUE) #leave out last protein in dataframe because there are no measurements for this one in phase1
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
