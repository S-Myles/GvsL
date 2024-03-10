# Packages 
library(phyloseq)
library(coda4microbiome)
library(tidyverse)

####
# Data Prep
####

# 12S
S12_df <- otu_table(S12_physeq_filt) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")

metadata <-  read_tsv("./data/metadata.txt")
metadata <- merge(metadata, S12_df, "sampleid")[,1:6]  # This is to make sure it's aligned
# Create binary vectors from metadata
S12_y_site <- metadata[,3] %>% 
  as.factor()
S12_y_season <- metadata[,5] %>% 
  as.factor()

S12_df <- column_to_rownames(S12_df, var = "sampleid")


# 16S
S16_df <- otu_table(S16_physeq_filt) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")

metadata <-  read_tsv("./data/metadata.txt")
metadata <- merge(metadata, S16_df, "sampleid")[,1:6]  # This is to make sure it's aligned
# Create binary vectors from metadata
S16_y_site <- metadata[,3] %>% 
  as.factor()
S16_y_season <- metadata[,5] %>% 
  as.factor()

S16_df <- column_to_rownames(S16_df, var = "sampleid")


# 18S
S18_df <- otu_table(S18_physeq_filt) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")

metadata <-  read_tsv("./data/metadata.txt")
metadata <- merge(metadata, S18_df, "sampleid")[,1:6]  # This is to make sure it's aligned
# Create binary vectors from metadata
S18_y_site <- metadata[,3] %>% 
  as.factor()
S18_y_season <- metadata[,5] %>% 
  as.factor()

S18_df <- column_to_rownames(S18_df, var = "sampleid")


# COI
COI_df <- otu_table(COI_physeq_filt) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sampleid")

metadata <-  read_tsv("./data/metadata.txt")
metadata <- merge(metadata, COI_df, "sampleid")[,1:6]  # This is to make sure it's aligned
# Create binary vectors from metadata
COI_y_site <- metadata[,3] %>% 
  as.factor()
COI_y_season <- metadata[,5] %>% 
  as.factor()

COI_df <- column_to_rownames(COI_df, var = "sampleid")



#####
# Log-ratio exploratory analysis
#####

# 12S
# Compute all log ratios (from x) and their association to the binary outcome y based on the glm function results
S12_seasonal_logratios<-explore_logratios(x=S12_df,y=S12_y_season, measure = "glm") # can change glm to be AUC instead
S12_site_logratios<-explore_logratios(x=S12_df,y=S12_y_site, measure = "glm")

# The correlation matrix given shows the strength of pairwise associations to outcome
# And the features are ordered in terms of their summed associations (across all pairwise ratios)
# So, this object stores that order
S12_seasonal_logratios$`order of importance`
# find out taxonomy
S12_seasonal_logratios$`name of most important variables`


# 16S
# Compute all log ratios (from x) and their association to the binary outcome y based on the glm function results
S16_seasonal_logratios<-explore_logratios(x=S16_df,y=S16_y_season, measure = "glm") 
S16_site_logratios<-explore_logratios(x=S16_df,y=S16_y_site, measure = "glm") 


# 18S
# Compute all log ratios (from x) and their association to the binary outcome y based on the glm function results
S18_seasonal_logratios<-explore_logratios(x=S18_df,y=S18_y_season, measure = "glm") 
S18_site_logratios<-explore_logratios(x=S18_df,y=S18_y_site, measure = "glm") 

# COI
# Compute all log ratios (from x) and their association to the binary outcome y based on the glm function results
COI_seasonal_logratios<-explore_logratios(x=COI_df,y=COI_y_season, measure = "glm") 
COI_site_logratios<-explore_logratios(x=COI_df,y=COI_y_site, measure = "glm") 
################################################################################

#####
# TUTORIAL Cross-sectional studies
#####
# Running the model
# Function recognizes the binary outcome so, implements a penalized logistic regression
# First plot generated monitors AUC drop off as feature space is reduced iteratively
# Dashed verticals show highest AUC score and cutoff selection for dimensionality reduction
# This cutoff value is decided based on the lambda argument (default is 1 standard dev)

# 12S
S12_coda_glmnet_seasonal<-coda_glmnet(x=S12_df,y=S12_y_season, nfolds = 3)
S12_coda_glmnet_site<-coda_glmnet(x=S12_df,y=S12_y_site, nfolds = 3)

# With unfiltered data
S16_coda_glmnet_seasonal<-coda_glmnet(x=S16_df,y=S16_y_season)
S16_coda_glmnet_site<-coda_glmnet(x=S16_df,y=S16_y_site, nfolds = 3)


# With unfiltered data
S18_coda_glmnet_seasonal<-coda_glmnet(x=S18_df,y=S18_y_season)
S18_coda_glmnet_site<-coda_glmnet(x=S18_df,y=S18_y_site, nfolds = 3)

# With unfiltered data
COI_coda_glmnet_seasonal<-coda_glmnet(x=COI_df,y=COI_y_season)
COI_coda_glmnet_site<-coda_glmnet(x=COI_df,y=COI_y_site)


####
# Model evaluation
###

# taxa kept for model
coda_glmnet_seasonal$taxa.num
# their names
coda_glmnet_seasonal$taxa.name
# Their log contrast coefficients
(coef <- coda_glmnet_seasonal$`log-contrast coefficients`)
# The sum of which is = 0, because this model is built as a log-contrast
# log-contrast preserve the (subcompositional) invariance principle

# Extracting and ordering the taxa positively associated with the prediction
positives<-which(coef>=0)
op<-order(coef[positives], decreasing = TRUE)
coda_glmnet_seasonal$taxa.name[positives[op]]

# The negatives
negatives<-which(coef<0)
on<-order(abs(coef[coef<0]), decreasing = TRUE)
coda_glmnet_seasonal$taxa.name[negatives[on]]

# Same thing Vizualized on a plot
coda_glmnet_seasonal$`signature plot`

# Total signature per sample (prediction)
head(coda_glmnet_seasonal$`predictions`)

# Classification accuracy scores
coda_glmnet_seasonal$`apparent AUC`
coda_glmnet_seasonal$`mean cv-AUC`
coda_glmnet_seasonal$`sd cv-AUC`

# Distributions of all samples signatures (predictions)
coda_glmnet_seasonal$`predictions plot`






# PERMUTATIONAL TEST OF SIGNIFICANCE
# Performs permutational cross-validation tests to evaluate the model ito the null (y permuted)
null_acc_seasonal<-coda_glmnet_null(x=ASV_table,y=y_season, niter=10) # for final evaluation, increase niter!

# 
summary(null_acc$"accuracy")
null_acc$"confidence interval"
