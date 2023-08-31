# setwd('/Volumes/Thibaud/BMIR/precog/precog/precogdata/TCGA.KMplots')
# On VM
setwd('/srv/shiny-server/precog/precogdata/TCGA.KMplots')

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

for (i in names(cancer.names)){
  command = paste('sudo mkdir ', i, '_DSS', sep = '')
  # print(command)
  system(command)
}