#Note ==========================================================================
#It remains a question whether this script can be applied to multiple treatments


##R version 3.5, platform x86_64-redhat-linux-gnu
options(stringsAsFactors = F)

source("/home/rstudio/rdata/PSM_testing/SexImm/Code/GIPW_function_omega.R")
source("/home/rstudio/rdata/PSM_testing/SexImm/Code/utils_functions.R")
#Input==========================================================================
#Two parts
#Phenotype+Molecular data 
#NA is not acceptable
####

#Define key parameters =========================================================
topic <- 'Microbiome'
analysis="gender" 
#Create analysis folders =======================================================
##setup new folder, deposit result in this folder
setwd('/home/rstudio/rdata/PSM_testing/SexImm/')
folder <- "PSM_Output"
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- folder
setwd(scripts.dir)


#Read data =====================================================================
#No duplicates
#data and feature must be well aligned
relative_abundance <- readr::read_rds('../Data/curatedMetagenomicData_relative_abundance.2021.10.19.value.merged.rds')
relative_abundance
sum(duplicated(colnames(relative_abundance))) #check any duplicates

meta.info <- readr::read_csv('../Data/Microbiome.Metadata.csv')
meta.info[1:5,1:10]
sum(duplicated(meta.info$sample_id)) ##some duplicates, watching out!!!!

meta.info.subset <- meta.info %>%
  filter(!is.na(age)) %>%
  select(sample_id,body_site,age,gender) %>%
  unique()

sum(colnames(relative_abundance) %in% meta.info.subset$sample_id)


sum.FeatureAll <- data.frame()

data <- meta.info.subset
data <- as.data.frame(data)
#data <- data[!data$age_at_initial_pathologic_diagnosis=='#N/A',]

#data$age_at_initial_pathologic_diagnosis <- as.numeric(data$age_at_initial_pathologic_diagnosis)
#data <- unique(data) 
rownames(data) <- data[,1]
data <- data[,-1]

# Feature
Feature <- t(relative_abundance)

valid.sample <- intersect(rownames(data),rownames(Feature))
data <- data[valid.sample,]
Feature <- Feature[valid.sample,]

if(analysis=="gender"){
  # convert female and male to numeric 1,0 to suppress the warning message in lm
  data$gender <- ifelse(data$gender=="female",1,0)
  colnames(data)[which(colnames(data)=="gender")] <- "Z"
}

# convert to dummy - for catorgrical variable
dummy.feature <- setdiff(colnames(data),c("Z","age"))#which column needs to be in dummy
if(length(dummy.feature)>0){
  data.dum <- dummy.data.frame(data, names=dummy.feature)
  dummy.list <- attr(data.dum,"dummies")
  rm.col <- c()
  for (i in 1:length(dummy.list))
  {
    rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
  }
  data.dum <- data.dum[,-rm.col]
  data.dum$X0 <- rep(1, nrow(data.dum))
  #form <- as.formula("Z~.") # should exclude X0
  exclude.col <- match(c("Z","X0"), colnames(data.dum))
  colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
  form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
}else{
  data.dum <- data
  data.dum$X0 <- rep(1, nrow(data.dum))
  #form <- as.formula("Z~.") # should exclude X0
  exclude.col <- match(c("Z","X0"), colnames(data.dum))
  colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
  form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
}


# perform calculation

Feature.pri <- rm.zero.col(Feature)
Feature.pri <- apply(Feature.pri,2,function(x){log2(x+1)})
folder <- paste(topic,"_",analysis,sep="")
if (!file.exists(folder)) { dir.create(folder) }

Feature.result <- weight.test(data.dum, form, Feature.pri, is.continuous=TRUE,
                              weight=ifelse(analysis=="gender","MW","ATT"),
                              mirror.plot=FALSE, topic, data.type= "Feature", 
                              outdir=paste(scripts.dir, "/",topic,"_",analysis,sep=""),
                              perm=FALSE)

sum.Immune <- summarize.p(Feature.pri, Feature.result, print=TRUE)
summarize.p(Feature.pri, Feature.result, print=TRUE, cutoff=0.05)
write.summary(sum.Immune, topic, analysis,"Immune")
write.result(Feature.result, topic, analysis,"Immune")
save(Feature.result, file=paste(topic,"_", analysis,"_result.RData",sep=""))
if(length(which(Feature.result$pvalue < 0.05)) > 0){
  sum.Immune <- data.frame(sum.Immune)
  sum.Immune$class <- rep(topic,times=nrow(sum.Immune))
  if(nrow(sum.FeatureAll) == 0){
    sum.FeatureAll <- sum.Immune
  }else{
    sum.FeatureAll <- rbind(sum.FeatureAll,sum.Immune)
  }
  
}
write.table(sum.FeatureAll,file="feature difference.across.topic.typesAll.txt",quote = F,sep="\t",row.names = F)

