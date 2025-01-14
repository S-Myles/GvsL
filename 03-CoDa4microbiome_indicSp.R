# Packages 
library(phyloseq)
library(coda4microbiome)
library(tidyverse)

load("input_objects.RData")
#load("results.RData")

# These indicator psecies analyses are to be done only on 1 size-fraction datasets because
# RDAs showed fundametal impact on the data. 

# Only investigate sample attributes showing a significant association to community organization 
# from rd-RDA anova analysis

# Year should not be a factor to investigate. So only Season, Depth, or Station if significant.

# Depth should be binned to a binary as photic (1 and 20m) or aphotic (250 m)


# So for 
# 12S : none (only year significant with 1 size fraction, station when all samples present)
# 16S : Depth & Season
# 18S : Depth & Season & Station
# COI : Depth & Season 


####
# Data Prep ------------  16S
####
S16_df <- otu_table(S16_filt_data) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")
metadata <-  read_tsv("./data/metadata.txt")
metadata <- merge(S16_df, metadata, "sampleid") # Merge by sampleid

# Create sample atteibute binary vectors from metadata
S16_y_season <- metadata[,"Season"] %>% 
  as.factor()
S16_y_depth <- ifelse(metadata[,"Depth"] %in% c(1, 20), "photic", "aphotic") %>% 
  as.factor()

S16_df <- column_to_rownames(S16_df, var = "sampleid")



#####
# Cross-sectional (binary) glm penalized net
#####
# Running the model
# Function recognizes the binary outcome so, implements a penalized logistic regression
# First plot generated monitors AUC drop off as feature space is reduced iteratively
# Dashed verticals show highest AUC score and cutoff selection for dimensionality reduction
# This cutoff value is decided based on the lambda argument (default is 1 standard dev)

# 16S
#S16_coda_glmnet_seasonal<-coda_glmnet(x=S16_df,y=S16_y_season, nfolds = 3)
#save.image("results.RData")
S16_coda_glmnet_depth<-coda_glmnet(x=S16_df,y=S16_y_depth, nfolds = 3)
save.image("results.RData")



####
# Model evaluation
###

# taxa kept for model
S16_coda_glmnet_seasonal$taxa.num
# their names
S16_coda_glmnet_seasonal$taxa.name
# Their log contrast coefficients
(coef <- S16_coda_glmnet_seasonal$`log-contrast coefficients`)
# The sum of which is = 0, because this model is built as a log-contrast
# log-contrast preserve the (subcompositional) invariance principle

# Extracting and ordering the taxa positively associated with the prediction
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
S16_coda_glmnet_seasonal$taxa.name[positives[op]]

# The negatives
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
S16_coda_glmnet_seasonal$taxa.name[negatives[on]]

# Same thing Vizualized on a plot
S16_coda_glmnet_seasonal$`signature plot`

# Total signature per sample (prediction)
head(S16_coda_glmnet_seasonal$`predictions`)

# Classification accuracy scores
S16_coda_glmnet_seasonal$`apparent AUC`
S16_coda_glmnet_seasonal$`mean cv-AUC`
S16_coda_glmnet_seasonal$`sd cv-AUC`

# Distributions of all samples signatures (predictions)
S16_coda_glmnet_seasonal$`predictions plot`




# PERMUTATIONAL TEST OF SIGNIFICANCE
# Performs permutational cross-validation tests to evaluate the model ito the null (y permuted)
#S16_coda_PERM_glmnet_seasonal<-coda_glmnet_null(x=S16_df, y=S16_y_season, niter=10) # for final evaluation, increase niter!

# 
#summary(S16_coda_PERM_glmnet_seasonal$"accuracy")
#S16_coda_PERM_glmnet_seasonal$"confidence interval"





####
# Data Prep ------------   18S
####

S18_df <- otu_table(S18_filt_data) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")
metadata <-  read_tsv("./data/metadata.txt")
metadata <- merge(S18_df, metadata, "sampleid") # Merge by sampleid

# Create sample atteibute binary vectors from metadata
S18_y_season <- metadata[,"Season"] %>% 
  as.factor()
S18_y_depth <- ifelse(metadata[,"Depth"] %in% c(1, 20), "photic", "aphotic") %>% 
  as.factor()

S18_df <- column_to_rownames(S18_df, var = "sampleid")



#####
# Cross-sectional (binary) glm penalized net
#####
# Running the model
# Function recognizes the binary outcome so, implements a penalized logistic regression
# First plot generated monitors AUC drop off as feature space is reduced iteratively
# Dashed verticals show highest AUC score and cutoff selection for dimensionality reduction
# This cutoff value is decided based on the lambda argument (default is 1 standard dev)

# 18S
S18_coda_glmnet_seasonal<-coda_glmnet(x=S18_df,y=S18_y_season, nfolds = 3)
save.image("results.RData")
S18_coda_glmnet_depth<-coda_glmnet(x=S18_df,y=S18_y_depth, nfolds = 3)
save.image("results.RData")


####
# Model evaluation
###

# taxa kept for model
S18_coda_glmnet_seasonal$taxa.num
# their names
S18_coda_glmnet_seasonal$taxa.name
# Their log contrast coefficients
(coef <- S18_coda_glmnet_seasonal$`log-contrast coefficients`)
# The sum of which is = 0, because this model is built as a log-contrast
# log-contrast preserve the (subcompositional) invariance principle

# Extracting and ordering the taxa positively associated with the prediction
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
S18_coda_glmnet_seasonal$taxa.name[positives[op]]

# The negatives
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
S18_coda_glmnet_seasonal$taxa.name[negatives[on]]

# Same thing Vizualized on a plot
S18_coda_glmnet_seasonal$`signature plot`

# Total signature per sample (prediction)
head(S18_coda_glmnet_seasonal$`predictions`)

# Classification accuracy scores
S18_coda_glmnet_seasonal$`apparent AUC`
S18_coda_glmnet_seasonal$`mean cv-AUC`
S18_coda_glmnet_seasonal$`sd cv-AUC`

# Distributions of all samples signatures (predictions)
S18_coda_glmnet_seasonal$`predictions plot`




# PERMUTATIONAL TEST OF SIGNIFICANCE
# Performs permutational cross-validation tests to evaluate the model ito the null (y permuted)
#S18_coda_PERM_glmnet_seasonal<-coda_glmnet_null(x=S18_df, y=S18_y_season, niter=10) # for final evaluation, increase niter!

# 
#summary(S18_coda_PERM_glmnet_seasonal$"accuracy")
#S18_coda_PERM_glmnet_seasonal$"confidence interval"






####
# Data Prep ------------  COI
####
COI_df <- otu_table(COI_filt_data) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")
metadata <-  read_tsv("./data/metadata.txt")
metadata <- merge(COI_df, metadata, "sampleid") # Merge by sampleid

# Create sample atteibute binary vectors from metadata
COI_y_season <- metadata[,"Season"] %>% 
  as.factor()
COI_y_depth <- ifelse(metadata[,"Depth"] %in% c(1, 20), "photic", "aphotic") %>% 
  as.factor()

COI_df <- column_to_rownames(COI_df, var = "sampleid")



#####
# Cross-sectional (binary) glm penalized net
#####
# Running the model
# Function recognizes the binary outcome so, implements a penalized logistic regression
# First plot generated monitors AUC drop off as feature space is reduced iteratively
# Dashed verticals show highest AUC score and cutoff selection for dimensionality reduction
# This cutoff value is decided based on the lambda argument (default is 1 standard dev)

# COI
COI_coda_glmnet_seasonal<-coda_glmnet(x=COI_df,y=COI_y_season, nfolds = 3)
save.image("results.RData")
COI_coda_glmnet_depth<-coda_glmnet(x=COI_df,y=COI_y_depth, nfolds = 3)
save.image("results.RData")


####
# Model evaluation
###

# taxa kept for model
COI_coda_glmnet_seasonal$taxa.num
# their names
COI_coda_glmnet_seasonal$taxa.name
# Their log contrast coefficients
(coef <- COI_coda_glmnet_seasonal$`log-contrast coefficients`)
# The sum of which is = 0, because this model is built as a log-contrast
# log-contrast preserve the (subcompositional) invariance principle

# Extracting and ordering the taxa positively associated with the prediction
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
COI_coda_glmnet_seasonal$taxa.name[positives[op]]

# The negatives
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
COI_coda_glmnet_seasonal$taxa.name[negatives[on]]

# Same thing Vizualized on a plot
COI_coda_glmnet_seasonal$`signature plot`

# Total signature per sample (prediction)
head(COI_coda_glmnet_seasonal$`predictions`)

# Classification accuracy scores
COI_coda_glmnet_seasonal$`apparent AUC`
COI_coda_glmnet_seasonal$`mean cv-AUC`
COI_coda_glmnet_seasonal$`sd cv-AUC`

# Distributions of all samples signatures (predictions)
COI_coda_glmnet_seasonal$`predictions plot`




# PERMUTATIONAL TEST OF SIGNIFICANCE
# Performs permutational cross-validation tests to evaluate the model ito the null (y permuted)
#COI_coda_glmnet_seasonal<-coda_glmnet_null(x=COI_df, y=COI_y_season, niter=10) # for final evaluation, increase niter!

# 
#summary(COI_coda_glmnet_seasonal$"accuracy")
#COI_coda_glmnet_seasonal$"confidence interval"


save.image("results.RData")