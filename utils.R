library('stringr')
# ----- Function to turn text into link -----
ToLink <- function(txt,link) {
  paste0('<a href=',link,">",txt,'</a>')
}
# ----- Function to reduce the length of a string -----
reduce_char_lenght <- function(char, max_size){
  
  char = str_remove_all(char, "\\)")
  char = str_remove_all(char, "\\(")
  char = str_remove_all(char, "\\+")
  
  
  if(nchar(char) > max_size){
    

    # char = str_remove_all(char, "\\-")
    # char = str_remove_all(char, "1")
    if (nchar(char) == max_size +1){
      to_remove = substr(char, max_size, nchar(char))
      # print(to_remove)
    }
    else{to_remove = substr(char, max_size+1, nchar(char))}
    to_remove = substr(char, max_size+1, nchar(char))
    return(str_replace(char, to_remove, "..."))
  }
  
  else(
    return(char)
  )
}
# ----- Function to add one row in a dataframe -----
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  return(existingDF)
}

# ----- Function to get the sample size of a study -----
getSampleSize <- function(cancer, accession, platform){
  cancerid <- paste(tolower(str_replace_all(cancer, ' ', '.')), accession, platform, sep = ".")
  clinical.annot <- read.table(paste("data/filtered.infofiles/",cancerid, ".tsv", sep = ""),header =T)
  return(dim(clinical.annot)[1])
}
# getSampleSize(cancer, accession, platform)
# ----- Function to compute column of metaz score for a given cancer -----
computeMetaZScore <- function(cancer){
  cancer <- str_replace(cancer, 'Cancer', 'cancer')
  cancer <- str_replace_all(cancer, '_', ' ')
  subtype <- gsub("^.*\\cancer","", cancer)
  subtype <- str_replace(subtype, " ", "")
  # Import PRECOG list of studies
  original_precog <- read.csv('data/original_precog.csv', header =T, sep = ',')[, -1]
  original_precog <- original_precog[which(original_precog$Cancer == cancer & original_precog$Subtype == subtype),]
  studies.id = as.character(original_precog$PRECOG.ID)
  # Import PRECOG individual zscores
  data <- read.csv('data/precog_ind_score.csv', header = T, sep = ',')
  metaZ.col = rep(0, dim(data)[1])
  norm.C = 0
  for (id in studies.id){
    cat('id: ', id, '\n')
    accession <- as.character(original_precog$Accession[which(original_precog$PRECOG.ID == id)])
    platform <- as.character(original_precog$Platform[which(original_precog$PRECOG.ID == id)])
    # precog.id <- as.character(original_precog$PRECOG.ID[which(original_precog$Accession == acc)])
    zscore.col <- data[, which(colnames(data) == paste(str_replace_all(cancer, ' ', '_'), id, sep = '_'))]
    ind <- ifelse(zscore.col != 0, 1, 0)
    
    sampleSize = getSampleSize(cancer, accession, platform)
    cat('Sample size:', sampleSize, '\n')
    # Stouffer's Z-score method
    w = sqrt(sampleSize)
    norm.C = norm.C + ind * w^2
    metaZ.col = metaZ.col + w * zscore.col 
  }
  return(cbind(as.matrix(data$Gene), as.matrix(metaZ.col/sqrt(norm.C))))
}
# metaz.col <- computeMetaZScore('Bladder cancer')
# View(metaz.col)
# original_metaz <- read.csv('data/metaz.csv', header =T, sep = ';')
# order(metaz.col)
# metaz.col <- metaz.col[match(original_metaz$Gene, metaz.col),]

# ----- Function to generate precog information table -----
generatePRECOGinfo <- function(precog){
  precog[, -1] <- lapply(precog[, -1], function(x) as.character(x))
  precog_info <- precog[,c(2,6)]
  precog_info$Accession_num <- 0
  precog$Probe.Level <- as.character(precog$Probe.Level)
  
  i = 1
  for (cancer in precog_info$Cancer){
    if (precog$Subtype[i] != " "){
      precog_info$Cancer[i] <- paste(cancer, precog$Subtype[i], sep = " " )
    }
    i = i + 1
  }
  i= 1
  for (cancer in precog_info$Cancer){
    precog_info$Accession_num[i] <- length(precog_info$Accession[precog_info$Cancer == cancer])
    i = i + 1
  }
  rownames(precog_info) <- c(1:nrow(precog))
  precog_info$acc_and_plat <- paste(precog$Accession, precog$Platform, sep = '::')
  precog_info$firstRow <- rep(1, nrow(precog))
  old <- precog_info$Cancer[1]
  offset = 1
  for (i in c(2:nrow(precog))){
    if (precog_info$Cancer[i] != old){
      precog_info$firstRow[i] <- i 
      old <- precog_info$Cancer[i]
      offset = 1
    }
    else{
      precog_info$firstRow[i] <- precog_info$firstRow[i-1]
      offset = offset + 1
    }
  }
  return(precog_info)
}


# ----- Function to import all data for a given cancer/accession/platform -----
importData <- function(cancer, accession, platform, return = TRUE){
  
  id <- paste(str_replace_all(cancer, ' ', '.'), accession, platform, 'data.RData',sep = '.')
  
  #######################
  ##### IMPORT DATA #####
  #######################
  cat('Start importing all files... \n')
  # precog main table 
  precog <- read.csv("data/precog.csv", header = T, sep = ",")[,-1]
  ## meta z (not preprocessed yet)
  metaz0 <- read.csv("data/metaz.csv", header = T ,sep = ";")
  metaz0$Gene <- as.character(metaz0$Gene)
  metaz0$Adrenocortical_cancer <- as.character(metaz0$Adrenocortical_cancer)
  metaz0$Name <- as.character(metaz0$Name)
  
  id <- paste(str_replace_all(cancer, ' ', '.'), accession, platform, 'data.RData',sep = '.')
  cat("✓ PRECOG master table and Metaz Score table succesfully imported!  \n")
  
  # Import matrix and platform data 
  # probe <- read.table("PRECOG_probe_level/Bladder_cancer.GSE5287.HGU133A_EntrezCDF.OS.unisurv_cns.txt", sep = "\t", fill = T, header =T)
  if (grepl("EntrezCDF", platform)){
    # Handling multiple files
    if (grepl("AB", platform)){
      platform2 <- str_replace(platform, "B", "")
    }
    else {platform2 <- platform}
    id <- paste(str_replace_all(cancer, ' ', '.'), accession, platform2, 'data.RData',sep = '.')
  }
  
  matrix <- readRDS(paste('precogdata/Data/', id, sep = ''))
  matrix$Name <- as.character(matrix$Name)
  colnames(matrix)[c(1,2)] <- c("Gene", "Description")
  cat("✓ Gene expression data succesfully imported! \n")
  
  #########################
  ##### PLATFORM FILE #####
  #########################
  path.platform <- paste('precogdata/Platforms', '/', platform, '_annot.txt', sep ='')
  platform_data <- read.csv(path.platform, header =T, sep="\t")
  ###########################
  ##### ANNOTATION FILE #####
  ###########################
  cancerid = tolower(paste(tolower(str_replace_all(cancer, ' ', '.')), accession, platform, sep = "."))
  # print(cancerid)
  annotation <- read.table(paste("precogdata/filtered.infofiles/",cancerid, ".tsv", sep = ""), 
                           header =T)
  # fill = T, header = T, quote ="")
  # annotation <- na.exclude(annotation)
  # annotation[,2] <- as.numeric(annotation[,2])
  # annotation[,1] <- as.character(annotation[,1])
  # cat('annot dim: ', dim(annotation), '\n')
  # isSuccess$platform = T
  # isSuccess$annotation = T
  cat("✓ Platform data succesfully imported! \n")
  
  # IMPORT THE INFO DT THAT WILL BE UPDATED
  # km_info = data.frame(Gene = metaz0$Gene)
  km_info <- read.csv("data/km_info.csv", header = T)[,-1]
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
  
  
  # ----- Prepprocess -----
  cat('Start preprocessing matrix data... \n')
  # ----- match patients for all datasets -----
  if (dim(annotation)[1] > dim(matrix)[2]){
    patients = colnames(matrix)[which(colnames(matrix) %in% annotation$Array)]
  }
  else{
    patients = annotation$Array[which(annotation$Array %in% colnames(matrix))]
  }
  cat("No of patients = ", length(patients), "\n")
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
  
  
  # if (length(which(is.na(data)))>0) {
  #   library(samr)
  #   xx=impute.knn(data)
  #   data=xx$data
  #   # cat(paste("# IMPUTEDMISSING: TRUE\n"),file="ALL.B-ALL.GSE18497.HGU133Plus2_EntrezCDF.OS.unisurv_cns.txt")
  # } else {
  #   # cat(paste("# IMPUTEDMISSING: FALSE\n"),file="ALL.B-ALL.GSE18497.HGU133Plus2_EntrezCDF.OS.unisurv_cns.txt")
  # }
  
  # names <- colnames(matrix)
  # datan <- as.data.frame(cbind(matrix[, 1:2], data))
  # colnames(datan) <- names
  cat("✓ Preprocessed done! \n")
  # ----- end of preprocess -----
  
  # ----- Return output -----
  res = list()
  res$precog = precog
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
# ----- Function to generate KM plot given cancer/accession/platform and gene index or gene name -----
generateKMPlot <- function(data, updateKM = TRUE, start_gene = 1 ,
                           maxIndexIter = dim(metaz0)[1], spe.gene = NULL){
  ####################
  # Note:
  # function depends on
  # - platform_data
  # - metaz0
  # - km_info
  # - matrix
  ####################
  precog = data$precog
  metaz0 = data$metaz0
  matrix = as.data.frame(data$matrix) # matrix data file
  clinical.annot = data$annotation[, c(1:3)] # clinical annotations
  clinical.annot[, 3] <- as.numeric(as.character(clinical.annot[, 3]))
  platform_data = data$platform 
  km_info = as.data.frame(data$km_info) # helper file for tracking KM status
  datan <- data$preprocessed.data # matrix data without annotations (subset of `matrix`)
  col <- data$colname
  #### LOOP START ####
  for (index in start_gene : maxIndexIter){
    # PREPROCESS AND RETRIEVE THE INFORMATION NEEDED
    # index = 3
    gene = metaz0[["Gene"]][index] # TODO: load list of genes instead
    if (!is.null(spe.gene)){
      if (spe.gene != gene){next}
    }
    gene_row_number_in_matrix = -1
    cat('GENE: ', gene, '-------------------------------', 'Number: ', index, '\n')
    
    # # If already processed :
    # if (km_info[km_info$gene == gene,][col]==1){
    #   cat("gene already processed \n")
    #   next}
    
    # ----- Find the gene row in `matrix` by using the platform data -----
    sequence = -1
    for (k in 1:dim(platform_data)[1]){
      if ((as.character(platform_data[k,3]) == gene) & (!is.na(as.character(platform_data[k,3])))){
        sequence = as.character(platform_data[k,1])
      }
    }
    
    for (k in 1:dim(matrix)[1]){
      if (matrix$Gene[k] == sequence){
        gene_row_number_in_matrix = k
      }
    }
    
    if (gene_row_number_in_matrix > 0){
      cat("Gene found in matrix file -------------> SUCCESS \n")
      # cat(gene_row_number_in_matrix, "\n")
      n_col = dim(datan)[2]
      # We compute the median using the matrix file :
      # median = median(as.numeric(matrix[gene_row_number_in_matrix,c(3:n_col)]), na.rm = T)
      median = median(as.numeric(datan[gene_row_number_in_matrix,c(1:n_col)]), na.rm = T)
      # Number of patients using the annotation file :
      n_patients = min(length(clinical.annot[,1]), dim(matrix)[2] - 2)
      # Copy of annotation in order to work at a gene scale
      # so we can handle patients with too many missing values for this specific gene
      # clinical.annot <- annotation[, c(1,2,3)]
      # If there are missing values for this patient on this gene, we remove the patient :
      list_of_patient_tobe_removed = c()
      # for (i in c(1:n_patients)){
      #   if (is.na(matrix[gene_row_number_in_matrix, i+2])){
      #     patient = colnames(matrix)[i+2]
      #     # print(patient)
      #     row_idx_in_annotations <- as.numeric(rownames(clinical.annot[clinical.annot$Array==patient,]))
      #     list_of_patient_tobe_removed = c(row_idx_in_annotations, list_of_patient_tobe_removed)
      #   }
      # }
      # cat('list of patients to be removed length: ', length(list_of_patient_tobe_removed), '\n')
      if (length(list_of_patient_tobe_removed) == (dim(matrix)[2] - 2)){
        cat('but NA values for all patients... \n')
        next}
      # Add median split variable : HIGH = 1 and LOW = 0
      clinical.annot$split = "Low"
      for (i in c(1:n_patients)){
        # print(i)
        if (datan[gene_row_number_in_matrix, i] > median &
            !is.na(datan[gene_row_number_in_matrix, i])){
          clinical.annot$split[i] = "High"}
      }
      # Remove all these patients for there are missing values with this gene :
      if (length(list_of_patient_tobe_removed) != 0){
        clinical.annot <- clinical.annot[-list_of_patient_tobe_removed,]
      }
      
      # print(list_of_patient_tobe_removed)
      
      # bruteforce to handle DSS status in annotation files
      if ("DSS_Time" %in% colnames(clinical.annot)){
        # print('step here')
        for (i in 1:length(colnames(clinical.annot))){
          val <- colnames(clinical.annot)[i]
          colnames(clinical.annot)[i] <- str_replace(val, "DSS", "OS")
        }
      }
      
      if (!('High' %in% clinical.annot$split)){
        legend = c('Low')
      }
      else{
        legend = c('High', 'Low')
      }
      
      ######################################
      ##### KM Plot after median split #####
      ######################################
      title = paste(matrix$Gene[gene_row_number_in_matrix], 
                    matrix$Description[[gene_row_number_in_matrix]], sep = " ")
      
      Y = Surv(clinical.annot$OS_Time, clinical.annot$OS_Status) 
      fit <- survfit(Surv(OS_Time, OS_Status) ~ split, 
                     data = clinical.annot)
      
      # ----- Cox fit -----
      # TODO: do it in a separate function
      x = datan[gene_row_number_in_matrix, ]
      cox.res <- coxph(Y ~ x)
      screg <- summary(cox.res)
      HR <- round(coef(screg)[,2], 3)
      HR.conf.L <- round(screg$conf.int[3], 2)
      HR.conf.H <- round(screg$conf.int[4], 2)
      # print(screg$coef)
      #----------
      
      gene = str_replace(gene, '/', '-')
      name_on_file = paste(gene, col, 'png',sep = ".")
      cancer_name <- str_replace_all(cancer, " ", "_")
      # old path 
      # path <- paste("www/KM_plots/", cancer_name ,"/", accession, "/",
      #               name_on_file, sep ="")
      path <- paste("precogdata/KM_plots/", name_on_file, sep ="")
      
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
    } # end of distinction with gene not available in the study.
    
    # Update km info DT
    # km_info <- read.csv("data/km_info.csv", header = T)[,-1]
    # col = tolower(str_replace_all(paste(cancer, accession, platform, sep ="."), " ", "."))
    
    # Should we write in the info datatable ?
    if (gene_row_number_in_matrix > 0 ){
      km_info[index,][col] = 1
      # cat('km[x, y] = ', km_info[index,][col], '\n')
    }
    else {
      km_info[index, ][col] == 0 
    }
    
    
  }
  #### LOOP END #####
  
  # Export the updated info datatable.
  if (updateKM){
    write.csv(km_info,"data/km_info.csv")
  }
  
  # return(km_fit_split)
  
}