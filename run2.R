source("KM.R")
setwd("/srv/shiny-server/precog")
# setwd("/Volumes/Thibaud/BMIR/precog/precog")
precog <- read.csv("data/precog.csv", header = T, sep = "\t")[c(-1,-2),]

# Handling command line arguments
args <- commandArgs(TRUE)
cancer = args[1]
computePlots = if (!is.na(args[2])) TRUE else (FALSE)

# TODO: handle cancers with subtypes
# TODO: `Breast cancer` issue with accession data E-TABM-158 that cannot be downloaded

if (cancer == 'Bladder_cancer'){
  cancer = 'Bladder Cancer'
  subset = precog[which(precog$Cancer == "Bladder cancer"),]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    cat("------------------------ \n")
    cat("cancer: ", cancer, "\n")
    cat("accession: ", accession, "\n")
    cat("platform: ", platform, "\n")
    
    # print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
    # KM plot computation
    if (computePlots & accession == 'GSE1827'){
      cat('Start computing KM plots... \n')
      KMdata <- importKMdata(cancer, accession, platform)
      maxIndexIter = 5
      generateKMPlot(KMdata, TRUE, 65, 71)
    }
  }
} # accession GSE13507 is very slow...

# FIXME
if (cancer == 'Brain_cancer_astrocytoma'){
  cancer = 'Brain Cancer Astrocytoma'
  subset0 = precog[which(precog$Cancer == "Brain cancer"),] # TODO: check
  subset = subset0[which(precog$Subtype == 'Astrocytoma'), ]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    cat("------------------------ \n")
    cat("cancer: ", cancer, "\n")
    cat("accession: ", accession, "\n")
    cat("platform: ", platform, "\n")
    
    print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
  }
  
}
if (cancer == 'Brain_cancer_glioblastoma'){
  cancer = 'Brain Cancer Astrocytoma'
  subset0 = precog[which(precog$Cancer == "Brain cancer"),] # TODO: check
  subet = subset0[which(precog$Subtype == 'Glioblastoma'), ]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    cat("------------------------ \n")
    cat("cancer: ", cancer, "\n")
    cat("accession: ", accession, "\n")
    cat("platform: ", platform, "\n")
    
    print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
  }
  
}
if (cancer == 'Brain_cancer_glioma'){
  cancer = 'Brain Cancer Glioma'
  subset0 = precog[which(precog$Cancer == "Brain cancer"),] # TODO: check
  subet = subset0[which(precog$Subtype == 'Glioma'), ]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    cat("------------------------ \n")
    cat("cancer: ", cancer, "\n")
    cat("accession: ", accession, "\n")
    cat("platform: ", platform, "\n")
    
    print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
  }
  
}
if (cancer == 'Brain_cancer_medulloblastoma'){
  cancer = 'Brain Cancer Medulloblastoma'
  subset0 = precog[which(precog$Cancer == "Breast cancer"),] # TODO: check
  subet = subset0[which(precog$Subtype == 'Medulloblastoma'), ]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    cat("------------------------ \n")
    cat("cancer: ", cancer, "\n")
    cat("accession: ", accession, "\n")
    cat("platform: ", platform, "\n")
    
    print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
  }
  
}
if (cancer == 'Brain_cancer_meningioma'){
  cancer = 'Brain Cancer Meningioma'
  subset0 = precog[which(precog$Cancer == "Brain cancer"),] # TODO: check
  subet = subset0[which(precog$Subtype == 'Meningioma'), ]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    cat("------------------------ \n")
    cat("cancer: ", cancer, "\n")
    cat("accession: ", accession, "\n")
    cat("platform: ", platform, "\n")
    
    print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
  }
  
}
if (cancer == 'Brain_cancer_neuroblastoma'){
  cancer = 'Brain Cancer Neuroblastoma'
  subset0 = precog[which(precog$Cancer == "Brain cancer"),] # TODO: check
  subet = subset0[which(precog$Subtype == 'Neuroblastoma'), ]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    cat("------------------------ \n")
    cat("cancer: ", cancer, "\n")
    cat("accession: ", accession, "\n")
    cat("platform: ", platform, "\n")
    
    print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
  }
  
}

if (cancer == 'Breast_cancer'){
  cancer = 'Breast Cancer'
  subset = precog[which(precog$Cancer == "Breast cancer"),]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    # cat("------------------------ \n")
    # cat("cancer: ", cancer, "\n")
    #  cat("accession: ", accession, "\n")
    # cat("platform: ", platform, "\n")
    
    # print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
    
    if (computePlots & accession == 'GSE29174'){
      cat('cancer: ', cancer, '\n')
      cat('accession: ', accession, '\n')
      cat('Start computing KM plots... \n')
      KMdata <- importKMdata(cancer, accession, platform)
      # maxIndexIter = 5
      generateKMPlot(KMdata, TRUE, 1, 23287)
    }
  }
  
}

if (cancer == 'Colon_cancer'){
  cancer = 'Colon Cancer'
  subset = precog[which(precog$Cancer == "Colon cancer"),]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    cat("------------------------ \n")
    cat("cancer: ", cancer, "\n")
    cat("accession: ", accession, "\n")
    cat("platform: ", platform, "\n")
    
   #  print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
    if (computePlots & accession == 'GSE16125'){
      cat('cancer: ', cancer, '\n')
      cat('accession : ', accession, '\n')
      cat('Start computing KM plots... \n')
      KMdata <- importKMdata(cancer, accession, platform)
      # maxIndexIter = 5
      generateKMPlot(KMdata, TRUE, 1, 23287) 
    }
  }
  
} # OK

if (cancer == 'Gastric_cancer'){
  cancer = 'Gastric Cancer'
  subset = precog[which(precog$Cancer == "Gastric cancer"),]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    cat("------------------------ \n")
    cat("cancer: ", cancer, "\n")
    cat("accession: ", accession, "\n")
    cat("platform: ", platform, "\n")
    
    print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
  }
  
} 

if (cancer == 'Ovarian_cancer'){
  cancer = 'Ovarian Cancer'
  subset = precog[which(precog$Cancer == "Ovarian cancer"),]
  accessions = as.character(subset$Accession)
  platforms =  as.character(subset$Platform)
  
  for (i in 1:length(accessions)){
    accession = accessions[i]
    platform = platforms[i]
    cat("------------------------ \n")
    cat("cancer: ", cancer, "\n")
    cat("accession: ", accession, "\n")
    cat("platform: ", platform, "\n")
    
    print(importKMdata(cancer, accession, platform, return = FALSE))
    cat("------------------------ \n")
  }
  
} # OK but one accession file cannot be downloade.
