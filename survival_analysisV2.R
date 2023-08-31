options(warn=-1)
suppressMessages(library(survival))
library(survminer)
library(sjmisc)
suppressMessages(library(ggfortify))
suppressMessages(library(stringr))
library(preprocessCore)

# VM use only 
# setwd('/srv/shiny-server/precog')

# ----- Function to import data for a given study -----
importData <- function(cancer, accession, platform, return = TRUE, RDS_format = FALSE){
  
  id <- paste(str_replace_all(cancer, ' ', '.'), accession, platform, 'data.RData',sep = '.')
  cancer <- str_replace(cancer, 'Cancer', 'cancer')
  
  ##### IMPORT GENE EXPRESSION DATA #####
  cat('Start importing all files... \n')
  # precog main table 
  precog <- read.csv("precogdata/PRECOG_V2/original_precog_V2.csv", header = T, sep = ",")[, -1]
  precog$Accession <- as.character(precog$Accession)
  precog$Platform <- as.character(precog$Platform)
  outcome <- as.character(precog$Outcome[which(precog$Accession == accession & precog$Platform == platform)])
  subset.file <- read.csv('precogdata/SUBSETFILE.txt', sep = '\t')
  subset.file$X.is.na..Array.. <- as.character(subset.file$X.is.na..Array..)
  scc = as.character(subset.file$SCC[which(subset.file$TCGA_LUSC == accession)])
  cat('Sub cancer class:', scc, '\n')
  ## meta z (not preprocessed yet)
  metaz0 <- read.csv("precogdata/PRECOG_V2/metaz_V2.csv", header = T ,sep = ",")
  metaz0$Gene <- as.character(metaz0$Gene)
  metaz0$Adrenocortical_cancer <- as.character(metaz0$Adrenocortical_cancer)
  metaz0$Name <- as.character(metaz0$Name)
  cat("✓ PRECOG master table and Metaz Score table succesfully imported!  \n")
  
  # Import matrix and platform data 
  # probe <- read.table("PRECOG_probe_level/Bladder_cancer.GSE5287.HGU133A_EntrezCDF.OS.unisurv_cns.txt", sep = "\t", fill = T, header =T)
  
  if (grepl("EntrezCDF", platform)){
    # Handling multiple files
    if (grepl("AB", platform)){
      platform2 <- str_replace(platform, "B", "")
    }else {platform2 <- platform}
  }else if (grepl('GPL96_97', platform) || grepl('GPL96/GPL97', platform) || grepl('GPL96/97', platform)){
    platform <- 'GPL96_GPL97'
    platform2 <- 'GPL96'
  }else{
    platform2 <- platform
  }
  # Handle different Hematopoietic sub cancers. 
  if (grepl('Hematopoietic', cancer)){
    if (grepl('B-ALL', cancer)){
      id <- paste(str_replace(cancer, 'Hematopoietic cancer', 'ALL.'), accession, platform2, 'data.RData',sep = '.')
      id <- str_replace_all(id, ' ', '')
    }else if (grepl('Burkitt lymphoma', cancer)){
      # TODO: add sub cancer check
      id <- paste(str_replace(cancer, 'Hematopoietic cancer', 'BL_DLBCL.'), accession, platform2, 'data.RData',sep = '.')
      id <- str_replace_all(id, ' ', '')
      id <- str_replace(id, 'Burkittlymphoma', 'Burkitts.lymphoma')
    }else if('Burkitts lymphoma' %in% scc){
      id <- paste(str_replace(cancer, 'Hematopoietic cancer', 'BL_DLBCL.'), accession, platform, 'data.RData',sep = '.')
      id <- str_replace_all(id, ' ', '')
    }else{
      id <- paste(str_remove(cancer, 'Hematopoietic cancer'), accession, platform2, 'data.RData',sep = '.')
      id <- str_replace_all(id, ' ', '.')
      id <- substr(id, 2, nchar(id))
    }
  }else if (cancer == 'Melanoma' && 'Primary' %in% scc){
    id <- paste(str_replace_all(cancer, ' ', '.'), "Primary",accession, platform2, 'data.RData',sep = '.')
  }else if (cancer == 'Ewing sarcoma' && 'Primary' %in% scc){
    id <- paste(str_replace_all(cancer, ' ', '.'), "Primary",accession, platform2, 'data.RData',sep = '.')
  }else if (cancer == 'Liver cancer' && 'Primary' %in% scc){
    id <- paste(str_replace_all(cancer, ' ', '.'), "Primary",accession, platform2, 'data.RData',sep = '.')
  }else if (grepl('Hypopharyngeal', cancer)){
    id <- paste("Hypopharyngeal.cancer",accession, platform2, 'data.RData',sep = '.')
  }else if (grepl('Oesophageal', cancer)){
    id <- paste("Oesophageal.cancer",accession, platform2, 'data.RData',sep = '.')
  }
  else if (grepl('Oral', cancer) && grepl('SCC', scc)){
    id <- paste("Oral.cancer.OralSCC",accession, platform2, 'data.RData',sep = '.')
  }
  else{
    # if (grepl(cancer, 'Mantle cell lymphoma')){
    #   
    # }
    id <- paste(str_replace_all(cancer, ' ', '.'), accession, platform2, 'data.RData',sep = '.')
  }
  # id <- paste(str_replace_all(cancer, ' ', '.'), accession, platform2, 'data.RData',sep = '.')
  cat('path: ', paste('precogdata/Data/', id, sep = ''), '\n')
  if (!RDS_format){
    matrix <- readRDS(paste('precogdata/Data/', id, sep = ''))
  }
  else{
    load(paste('precogdata/Data/', id, sep = ''))
    assign('matrix', d)
    rm(d)
  }
  colnames(matrix)[c(1,2)] <- c("Gene", "Description")
  matrix$Description <- as.character(matrix$Description)
  cat("✓ Gene expression data succesfully imported! \n")
  
  ##### IMPORT PLATFORM FILE #####
  path.platform <- paste('precogdata/Platforms', '/', platform2, '_annot.txt', sep ='')
  cat('path platform: ', path.platform, '\n')
  platform_data <- read.csv(path.platform, header =T, sep="\t")
  cat("✓ Platform data succesfully imported! \n")
  ##### IMPORT ANNOTATION FILE #####
  
  cancerid = tolower(paste(tolower(str_replace_all(cancer, ' ', '.')), accession, platform, sep = "."))
  cancerid <- str_replace_all(cancerid, '\\.\\.', '\\.')
  # print(cancerid)
  cat('annot path: ', paste("precogdata/filtered_clinical_annotations/",cancerid, ".tsv \n", sep = ""))
  annotation <- read.table(paste("precogdata/filtered_clinical_annotations/",cancerid, ".tsv", sep = ""), 
                           header =T)
  
  cat("✓ Clinical annotation data succesfully imported! \n")
  
  # IMPORT THE INFO DT THAT WILL BE UPDATED
  # km_info = data.frame(Gene = metaz0$Gene)
  km_info <- read.csv("precogdata/km_info.csv", header = T)[,-1]
  col = tolower(str_replace_all(paste(cancer, accession, platform, sep ="."), " ", "."))
  if (!col %in% colnames(km_info)){
    cat(col, "NOT PROCESSED YET \n")
    km_info[[col]] = -1
  } else {cat(col, "ALREADY PROCESSED \n")}
  
  # isSuccess$km_info = T
  
  # Print output
  # cat("Import Data Completion \n")
  # cat(rep("-", 11), "\n")
  # for (i in 1:6){
  #   cat(names(isSuccess)[i],
  #       rep("", 12 - nchar(names(isSuccess)[i])),
  #       ": " ,
  #       isSuccess[[i]],
  #       "\n")
  # }
  
  
  # ----- Preprocess -----
  cat('Start preprocessing matrix data... \n')
  # match patients for all datasets
  if (dim(annotation)[1] > dim(matrix)[2]){
    patients = colnames(matrix)[which(colnames(matrix) %in% annotation$Array)]
  }
  else{
    patients = annotation$Array[which(annotation$Array %in% colnames(matrix))]
  }
  cat("No of patients = ", length(patients), '(Sub cancer class:', scc ,')', "\n")
  annot = matrix[, c(1:2)]
  matrix <- matrix[, which(colnames(matrix) %in% patients)]
  matrix <- matrix[, order(colnames(matrix))]
  # annot <- annot[, order(colnames(matrix))]
  matrix <- cbind(annot, matrix)
  annotation <- annotation[which(annotation$Array %in% patients), ]
  annotation <- annotation[order(annotation$Array), ]
  
  data = matrix(as.numeric(unlist(matrix)),nrow=nrow(matrix))[, -c(1:2)]
  
  # View(data)
  # print(sapply(matrix, class))
  maxex=max(data, na.rm=T)
  # print(maxex)
  if (maxex > 100) {  # Check if in log2 space
    data=log2(data)
    # matrix[, -c(1:2)]=ifelse(is.finite(matrix[, -c(1:2)]),matrix[, -c(1:2)],NA)
  }
  # print(class(matrix))
  data=normalize.quantiles(data,copy=T)
  
  sds=apply(data, 1,sd,na.rm=T)
  mns= apply(data,1,mean,na.rm=T)
  mds= apply(data,1,median,na.rm=T)
  
  # Filter out genes with no variation
  lo=which(sds<0.00001)
  if (length(lo)>0) {
    # annot=annot[-lo,]
    data = data[-lo, ]
    matrix = matrix[-lo, ]
    sds=sds[-lo]
    mns=mns[-lo]
    mds=mds[-lo]
  }
  
  # Filter out rows/cols with more than 80% NAs
  xy=dim(data)
  rownas=rowSums(is.na(data))
  colnas=colSums(is.na(data))
  rrow=which(rownas>=0.8*xy[2])
  rcol=which(colnas>=0.8*xy[1])
  if (length(rrow>0)) {
    # annot=annot[-rrow,]
    data = data[-rrow, ]
    matrix <- matrix[-rrow, ]
    sds = sds[-rrow]
    mns = mns[-rrow]
    mds = mds[-rrow]
  }
  if (length(rcol)>0) {
    data =  data[, -rcol]
    annotation = annotation[-rcol, ]
    # infon=infon[-rcol,]
  }
  
  xy=dim(data)
  
  # Standardize genes
  data = (data-mns)/sds
  
  cat("✓ Preprocessed done! \n")
  # end of preprocess 
  # ----- Return output -----
  res = list()
  res$precog = precog
  res$outcome = outcome
  res$metaz0 = metaz0
  res$matrix = matrix
  res$preprocessed.data <- data
  res$platform = platform_data
  res$annotation = annotation
  res$km_info = km_info
  res$colname <- col
  
  if (return){
    return (res)
  }
  else{
    return ("Nothing returned :)")
  }
}

# ----- Function to generate KM plot -----
generateKMPlots <- function(data, updateKM = TRUE, start, end, spe.gene = NULL){
  # Import parameters 
  precog = data$precog
  outcome = data$outcome
  metaz0 = data$metaz0
  matrix = as.data.frame(data$matrix) # matrix data file
  clinical.annot = data$annotation[, c(1:3)] # clinical annotations
  clinical.annot[, 3] <- as.numeric(as.character(clinical.annot[, 3]))
  platform_data = data$platform 
  km_info = as.data.frame(data$km_info) # helper file for tracking KM status
  datan <- data$preprocessed.data # gene X data without annotations (subset of `matrix`)
  col <- data$colname # specific id for this study
  
  # Fit survival object
  # bruteforce to handle DSS status in annotation files
  if ("DSS_Time" %in% colnames(clinical.annot)){
    # print('step here')
    for (i in 1:length(colnames(clinical.annot))){
      val <- colnames(clinical.annot)[i]
      colnames(clinical.annot)[i] <- str_replace(val, "DSS", "OS")
    }
    # outcome = 'DSS'
  }
  if ("DFS_Time" %in% colnames(clinical.annot)){
    # print('step here')
    for (i in 1:length(colnames(clinical.annot))){
      val <- colnames(clinical.annot)[i]
      colnames(clinical.annot)[i] <- str_replace(val, "DFS", "OS")
    }
    # outcome = 'DSS'
  }
  Y = Surv(clinical.annot$OS_Time, clinical.annot$OS_Status) 
  # Loop over rows of gene X data
  for (index in start:min(dim(matrix)[1], end)){
    probe = as.character(matrix[index, 1])
    cat('Probe sequence: ', probe, '(', index, 'out of', end, ') \n')
    gene = as.character(sub("\\ -.*", "", matrix[index, 2]))
    # gene = str_replace_all(gene, '/', '-')
    if (gene == "" || 
        gene == 'NA' || 
        grepl('OK', gene) ||
        grepl('/', probe) ||
        is.na(gene) || 
        gene == '---' || 
        str_contains(gene, 'or') || 
        gene == '#NAME?' ||
        !(gene %in% km_info$gene)){next}
    
    # When a specific gene is given as an argument
    if (!is.null(spe.gene)){
      if (spe.gene != gene){next}
    }
    
    cat('Gene: ', gene, '\n')
    
    n_col = dim(datan)[2]
    # We compute the median using the matrix file :
    median = median(as.numeric(datan[index,c(1:n_col)]), na.rm = T)
    # Number of patients using the annotation file :
    n_patients = min(length(clinical.annot[,1]), dim(matrix)[2] - 2)
    # Copy of annotation in order to work at a gene scale
    # so we can handle patients with too many missing values for this specific gene
    # clinical.annot <- annotation[, c(1,2,3)]
    # If there are missing values for this patient on this gene, we remove the patient :
    list_of_patient_tobe_removed = c()
    if (length(list_of_patient_tobe_removed) == (dim(matrix)[2] - 2)){
      cat('but NA values for all patients... \n')
      next}
    # Add median split variable : HIGH = 1 and LOW = 0
    clinical.annot$split = "Low"
    for (i in c(1:n_patients)){
      # print(i)
      if (datan[index, i] > median &
          !is.na(datan[index, i])){
        clinical.annot$split[i] = "High"}
    }
    # Remove all these patients for there are missing values with this gene :
    if (length(list_of_patient_tobe_removed) != 0){
      clinical.annot <- clinical.annot[-list_of_patient_tobe_removed,]
    }
    
    # # bruteforce to handle DSS status in annotation files
    # if ("DSS_Time" %in% colnames(clinical.annot)){
    #   # print('step here')
    #   for (i in 1:length(colnames(clinical.annot))){
    #     val <- colnames(clinical.annot)[i]
    #     colnames(clinical.annot)[i] <- str_replace(val, "DSS", "OS")
    #   }
    #   # outcome = 'DSS'
    # }
    # if ("DFS_Time" %in% colnames(clinical.annot)){
    #   # print('step here')
    #   for (i in 1:length(colnames(clinical.annot))){
    #     val <- colnames(clinical.annot)[i]
    #     colnames(clinical.annot)[i] <- str_replace(val, "DSS", "OS")
    #   }
    #   # outcome = 'DSS'
    # }
    
    if (!('High' %in% clinical.annot$split)){
      legend = c('Low')
    }
    else{
      legend = c('High', 'Low')
    }
    
    ######################################
    ##### KM Plot after median split #####
    ######################################
    title = paste('GeneID: ', matrix$Gene[index], '- Gene Symbol:',
                  matrix$Description[[index]], sep = " ")
    
    # Y = Surv(clinical.annot$OS_Time, clinical.annot$OS_Status) 
    fit <- survfit(Surv(OS_Time, OS_Status) ~ split, 
                   data = clinical.annot)
    
    # ----- Cox fit -----
    x = datan[index, ]
    cox.res <- coxph(Y ~ x)
    screg <- summary(cox.res)
    HR <- round(coef(screg)[,2], 3)
    HR.conf.L <- round(screg$conf.int[3], 2)
    HR.conf.H <- round(screg$conf.int[4], 2)
    # print(screg$coef)
    # HR = NA
    # HR.conf.L = NA
    # HR.conf.H = NA
    # ------ Compute KM plot ------
    gene = str_replace_all(gene, '/', '-')
    file_name = tolower(paste(col, gene, probe, outcome, 'png',sep = "."))
    # cancer_name <- str_replace_all(cancer, " ", "_")
    folder <- paste(col, tolower(outcome), 'KMplots', sep = '.')
    path <- paste("precogdata/KMplots", folder, file_name, sep ="/")
    cat('EXPORT PATH: ', path, '\n')
    png(file=path, width=500, height=500, res = 65)
    km <- ggsurvplot(fit, 
                     clinical.annot,
                     # main = 'Survival',
                     xlab = 'Time (months)',
                     # conf.int = T,
                     # conf.int.alpha = 0.15,
                     censor = T, 
                     censor.shape = 124,
                     censor.size = 6, 
                     pval = T, 
                     pval.coord = c(.1, .075),
                     pval.size = 5,
                     title  = title,
                     font.main = c(16, "bold", "darkblue"),
                     font.tickslab = c(15, "bold", "black"),
                     
                     surv.scale = 'percent', 
                     legend.title = 'Median split',
                     legend.labs = legend, 
                     font.legend = c(12, "bold", "black"),
                     legend = c(.9, .9),
                     font.x = c(15, "bold", "black"),
                     font.y = c(15, "bold", "black"), 
                     risk.table = T, 
                     risk.table.y.text.col = T, 
                     risk.table.x.text = F, 
                     risk.table.fontsize = 6,
                     risk.table.col = 'strata',
                     tables.theme = theme_cleantable(),
                     # risk.table.fontsize.x = c(10, 15, 10),
                     
                     # y.text = F,   
                     # ylab = "",  
                     # xlab = "",
                     # font.tickslab = c(20, "bold"),
                     
                     surv.plot.height = .8, 
                     risk.table.height = 0.2,
                     # palette = c("#E7B800", "#2E9FDF"),
                     ggtheme = theme_light())
    # customised the plot
    km$plot <- km$plot +
      ggplot2::annotate(
        geom = "text",
        x = 0,
        y = 0,
        vjust = 0,
        hjust = 0,
        label = paste("HR = ", HR,
                      ' (95% CI = ', HR.conf.L, '-', HR.conf.H, ')',
                      sep = ''),
        size = 5
      ) 
    
    print(km, main = title)
    
    
    
    dev.off()
    
    # } # end of distinction with gene not available in the study.
    
    # Update km info DT
    # km_info <- read.csv("data/km_info.csv", header = T)[,-1]
    # col = tolower(str_replace_all(paste(cancer, accession, platform, sep ="."), " ", "."))
    
    # Should we write in the info datatable ?
    
    km_info[which(km_info$gene == gene),][col] <- 1
    cat('---------------------------------------- \n')
    
  } # Loop end
  # ----- Export the updated info datatable -----
  if (updateKM){
    write.csv(km_info,"precogdata/km_info.csv")
  }
}

# ----- generate Zscore -----
generateCoxFitNEW <- function(data, start = 1 ,
                              end,
                              gene.list = NULL){
  
  precog = data$precog
  metaz0 = data$metaz0
  matrix = as.data.frame(data$matrix) # matrix data file
  clinical.annot = data$annotation # clinical annotations
  platform_data = data$platform 
  km_info = as.data.frame(data$km_info) # helper file for tracking KM status
  datan <- data$preprocessed.data # matrix data preprocessed and without annotations (subset of `matrix`)
  col <- data$colname
  
  
  
  if ("DSS_Time" %in% colnames(clinical.annot)){
    # print('step here')
    for (i in 1:length(colnames(clinical.annot))){
      val <- colnames(clinical.annot)[i]
      colnames(clinical.annot)[i] <- str_replace(val, "DSS", "OS")
    }
    # outcome = 'DSS'
  }
  if ("DFS_Time" %in% colnames(clinical.annot)){
    # print('step here')
    for (i in 1:length(colnames(clinical.annot))){
      val <- colnames(clinical.annot)[i]
      colnames(clinical.annot)[i] <- str_replace(val, "DFS", "OS")
    }
    # outcome = 'DSS'
  }
  
  # clinical.annot$OS_Time <- as.numeric(clinical.annot$OS_Time)
  
  # Fitting cox regression
  Y = Surv(clinical.annot$OS_Time, clinical.annot$OS_Status) 
  zscore_matrix <- cbind(metaz0[, 1], NA, NA, NA, NA, NA)
  
  for (index in start:min(dim(matrix)[1], end)){
    probe = as.character(matrix[index, 1])
    cat('Probe sequence: ', probe, '(', index, 'out of', end, ') \n') # for debugging purposes
    gene = as.character(sub("\\ -.*", "", matrix[index, 2]))
    if (gene == "" || 
        gene == 'NA' || 
        grepl('OK', gene) ||
        is.na(gene) || 
        gene == '---' || 
        str_contains(gene, 'or') || 
        gene == '#NAME?' ||
        !(gene %in% km_info$gene)){next}
    
    
    cat('Gene: ', gene, '\n') # for debugging purposes
    cox_res <- coxph(Y ~ datan[index, ])
    screg <- summary(cox_res)
    zscore <- screg$coef[4]
    # cat('zscore = ', zscore, '\n')
    gene_row = as.numeric(rownames(metaz0[which(metaz0$Gene == gene), ]))
    for (i in 2:6){
      if (is.na(zscore_matrix[gene_row, i])){
        zscore_matrix[gene_row, i] = zscore
        break
      } else{next}
    }
  }
  
  # HR <- as.numeric(screg$coef[2])
  # HR <- round(coef(screg)[,2], 3)
  # HR.conf.L <- round(screg$conf.int[3], 2)
  # HR.conf.H <- round(screg$conf.int[4], 2)
  
  # Output
  # if (gene_row_number_in_matrix == -1){
  #   return(NA)
  # }
  # else{
  #   res = list("zscore"= screg$coef[4], "HR" = HR, 
  #              "HR.low" = HR.conf.L, "HR.95.high" = HR.conf.H,
  #              "screg" = screg)
  #   
  #   return(res)
  # }
  
  zscore_matrix <- as.data.frame(zscore_matrix)
  colnames(zscore_matrix)[1] <- 'Gene' 
  
  return(zscore_matrix)
  
}

# ----- test -----
# cancer = 'Head and neck cancer'
# accession = 'GSE686'
# platform = 'GPL503'
# # 
# data <- importData(cancer, accession, platform, TRUE, TRUE)

# generateKMPlot(data, FALSE, 1, 23287, 'AKT3')
# generateKMPlot(data, FALSE, 1, 1)
# s = Sys.time()
# generateKMPlots(data, TRUE, 1, 2)
# print(difftime(Sys.time(), s, units = 'secs'))

# s = Sys.time()
# output <- generateCoxFitNEW(data, 1, 23287)
# output2 <- sapply(output, function(x){as.numeric(as.character(x))})
# print(difftime(Sys.time(), s, units = 'secs'))

