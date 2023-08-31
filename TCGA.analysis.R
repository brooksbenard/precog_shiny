##################
# TCGA analysis #
#################
# ----- File meant to generate KM plots for TCGA individual studies -----

# setwd("/Volumes/Thibaud/BMIR/precog/precog")
# ----- on VM -----
setwd("/srv/shiny-server/precog/")

options(warn=-1)
suppressMessages(library(survival))
library(survminer)
suppressMessages(library(ggfortify))
suppressMessages(library(stringr))

# ------ Import files -----
TCGA.annot.file <- read.csv('precogdata/TCGA/TCGA-clinical-extraendpoints.txt', 
                       header =T, sep = '\t')
TCGA.annot.file$type <- as.character(TCGA.annot.file$type)
TCGA.annot.file$OS_Time <- (12/365) * TCGA.annot.file$OS_Time # months to days
# TODO: split up files


##########################
# TCGA cancer dictionary #
##########################
cancer.names <- list('ACC' = 'Adrenocorticol carcinoma',
                     'BLCA' = 'Bladder Urothelial Carcinoma',
                     'BRCA' = 'Breast Invasive Carcinoma',
                     'CESC' = 'Cervical Squamous Cell Carcinoma and Endocervical and Adenoarcinoma',
                     'CHOL' = 'Cholangiocarcinoma',
                     'COAD' = 'Colon Adenocarcinoma',
                     'DLBC' = 'Diffuse Large B Cell Lymphoma',
                     'ESCA' = 'Esophageal carcinoma',
                     'GBM' = 'Glioblastoma multiforme', 
                     'HNSC' = 'Head and Neck Squamous Cell Carcinoma', 
                     'KICH' = 'Kidney Chromophobe', 
                     'KIRC' = 'Kidney Renal Clear Cell Carcinoma', 
                     'KIRP' = 'Kidney Renal Papillary Cell Carcinoma', 
                     'LAML' = 'Acute Myeloid Leukemia', 
                     'LGG' = 'Brain Lower Grade Glioma', 
                     'LIHC' = 'Liver Hepatocellular Carcinoma',
                     'LUAD' = 'Lung Adenocarcinoma', 
                     'LUSC' = 'Lung Squamous Cell Carcinoma',
                     'MESO' = 'Mesothelioma',
                     'OV' = 'Ovarian Serous Cystadenocarcinoma',
                     'PAAD' = 'Pancreatic Adeno Carcimona', 
                     'PCPG' = 'Phechormocyotma and Paraganglioma',
                     'PRAD' = 'Prostate Adenocarcinoma', 
                     'READ' = 'Rectum Adenocarcinoma', 
                     'SARC' = 'Sarcoma', 
                     'SKCM' = 'Skim Cutaneous Melanoma',
                     'STAD' = 'Stomach adenocarcinoma',
                     # 'TCHA' = 'Tyroid Carcinoma', 
                     'TGCT' = 'Testicular Germ Cell Tumors',
                     'THCA' = 'Thyroid carcinoma',
                     'THYM' = 'Thymoma',
                     'UCEC' = 'Uterine Corpus Endometrial Carcinoma',
                     # 'UCF' = 'Uterine Corpus Endometrial Carcinoma', 
                     'UCS' = 'UCS',
                     'UVM' = 'Uveal Melanoma')


metaz <- read.csv('precogdata/TCGA/TCGA.metaz.csv', header =T)
gene.list = as.vector(metaz$Name)
cancer.list = names(cancer.names)[-c(1:10)]
gene.list = gene.list
# gene.list = c('AATK.AS1')

# -----------------------------
# ----- Generate KM plots -----
# -----------------------------
for (cancer.id in cancer.list){
  # ----- Data -----
  # cancer.id = 'ACC'
  if (cancer.id != 'KICH'){next}
  cancer = cancer.names[cancer.id]
  clinical.annot <- TCGA.annot.file[which(TCGA.annot.file$type == cancer.id), ]
  clinical.annot <- clinical.annot[, c('ID', 'type', 'OS_Status', 'OS_Time')]
  path <- paste('precogdata/TCGA/TCGA_', cancer.id, '_tpm.fullIDs.remapped.tsv', sep = '')
  data <- read.csv(path, header = T, sep = '')
  path2 <- str_replace('precogdata/TCGA/TCGA_tumors_prognostic_zscores/TCGA-XXX-prognostic-zscore-all.csv', 
                      'XXX', cancer.id)
  cxres <- read.csv(path2, header = T, sep = ',')
  
  cat('cancer: ', as.character(cancer), '\n')
  
  # ----- Loop start -----
  for (gene in gene.list){
    if (gene!= 'KRT76'){next}
    # index = 1
    # gene_row_number_in_matrix = -1
    # cat('GENE: ', gene, '-------------------------------', 'Number: ', index, '\n')
    # cat('gene: ', class(gene), '... \n')
    gene = str_replace_all(gene, '\\.', '-')
    cat('gene: ', gene, '... \n')
    gene.data <- data[which(data$Gene == gene), ]
    n_patients <- dim(clinical.annot)[1] #FIXME
    n_col = ncol(gene.data)
    list_of_patient_tobe_removed = c()
    gene_row_number_in_matrix <- as.numeric(rownames(gene.data))
    median = median(as.numeric(gene.data[1, 2:n_col]), na.rm = T)
    if (is.na(median)){next}
    
    # ----- use Cox regression results -----
    HR = round(cxres[which(cxres$gene == gene), ]$HR, 3)
    HR.conf.H = round(cxres[which(cxres$gene == gene), ]$CI_down, 2) 
    HR.conf.L = round(cxres[which(cxres$gene == gene), ]$CI_up, 2)

  
    # Add median split variable : HIGH = 1 and LOW = 0
    annotation_gene_scale <- clinical.annot
    annotation_gene_scale$ID <- as.character(annotation_gene_scale$ID)
    annotation_gene_scale$split <- "Low"
    for (i in c(1:n_patients)){
      patient.id <- annotation_gene_scale$ID[i]
      if (length(which(grepl(patient.id, colnames(gene.data)))) == 0){
        list_of_patient_tobe_removed = c(i, list_of_patient_tobe_removed)
        score = 0
      } else if (is.na(annotation_gene_scale$OS_Time[i])){
        list_of_patient_tobe_removed = c(i, list_of_patient_tobe_removed)
      } else{
          if (length(which(grepl(patient.id, colnames(gene.data)))) == 1){
            score <- as.numeric(gene.data[, which(grepl(patient.id, colnames(gene.data)))])
          } else{
            score <- as.numeric(gene.data[, which(grepl(patient.id, colnames(gene.data)))[1]])
            }
      }
      # score <- as.numeric(gene.data[, which(grepl(patient.id, colnames(gene.data)))])
      if (score > median){annotation_gene_scale$split[i] <- "High"}
    }
    # Remove all these patients for there are missing values with this gene :
    if (length(list_of_patient_tobe_removed) != 0){
      annotation_gene_scale <- annotation_gene_scale[-list_of_patient_tobe_removed, ]
    }
    
    if (!('High' %in% annotation_gene_scale$split)){
      # print('skip')
      next
    }
    
    
    
    #### Create suvival object ####
    title = paste(gene, 'in', cancer, sep = ' ')
    

  
    fit <- survfit(Surv(OS_Time, OS_Status) ~ split, 
                   data = annotation_gene_scale)
    
    name_on_file = paste(gene, 'png', sep = ".")
    cancer_name <- str_replace(cancer, " ", "_")
    path <- paste("www/precogdata/TCGA.KMplots/", cancer.id , "/",
                  name_on_file, sep ="")

    png(file=path, width=500, height=500, res = 65)
    km <- ggsurvplot(fit, 
                     annotation_gene_scale,
                     # main = 'Survival',
                     xlab = 'Time (months)',
                     conf.int = T,
                     conf.int.alpha = 0.15,
                     censor = T, 
                     pval = T, 
                     pval.coord = c(.1, .075),
                     pval.size = 5,
                     title  = title,
                     font.main = c(16, "bold", "darkblue"),
                     font.tickslab = c(15, "bold", "black"),
                     
                     surv.scale = 'percent', 
                     legend.title = 'Median split',
                     legend.labs = c("High", "Low"), 
                     font.legend = c(12, "bold", "black"),
                     legend = c(.9, .9),
                     font.x = c(15, "bold", "black"),
                     font.y = c(15, "bold", "black"), 
                     risk.table = T, 
                     risk.table.y.text.col = T, 
                     risk.table.fontsize = 6,
                     risk.table.col = 'strata',
                     # risk.table.fontsize.x = c(10, 15, 10),
                     
                     # y.text = F,   
                     # ylab = "",  
                     # xlab = "",
                     # font.tickslab = c(20, "bold"),
                     
                     surv.plot.height = .8, 
                     risk.table.height = 0.2,
                     ggtheme = theme_light())
    # customized the plot
    if (length(HR) == 0 || !is.na(HR)){
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
    }
    
    print(km, main = title)
    dev.off()
    
    cat('Done! \n')
    
  }
  cat('-------------- \n')
}

# ----- Loop end -----




