setwd('/srv/shiny-server/precog')
source('code/utils.R')
precog <- read.csv('precogdata/original_precog.csv', header = T, sep = ',')[, -1]
# cancer = 'Head and neck cancer'
# subtype = ''
# accession = 'GSE78060'
# platform = 'GPL571'

precog <- read.csv('precogdata/PRECOG_V2/original_precog_V2.csv', header = T, ',')[, -1]
for (i in 1:dim(precog)[1]){
  cancer = as.character(precog$Cancer[i])
  if (cancer == 'Head and neck cancer'){
    # cancer = 'Head.and.neck.cancer'
    accession = as.character(precog$Accession[i])
    platform = as.character(precog$Platform[i])
    subtype = ''
    print(cancer)
    print(accession)
    print(platform)
    
    study_id <- paste(cancer, subtype, accession, platform, sep = '.')
    study_id <- str_replace_all(study_id, ' ', '.')
    study_id <- str_replace_all(study_id, '\\.\\.', '.')
    
    gene_expression_matrix = paste('precogdata/Data/', study_id, '.data.RData', sep ='')
    outcome = 'OS'
    analysis_type = 'Univariate COX Continuous'
    output_file = ''
    n_samples_used = getSampleSize(cancer, accession, platform)
    
    # ----- step 1: Generate file header -----
    filename <- paste(study_id, outcome, 'unisurv_cns.txt',sep = '.')
    probe_file <- file(paste('www/PRECOG_probe_level/', filename, sep =''))
    lines <- c(paste('# DISEASE: ', cancer, ''), 
               paste("# SUBSET: ", subtype, ''),
               paste('# SERIES: ', accession, ''),
               paste('# PLATFORM: ', platform, ''),
               paste('# DataMatrix: ', gene_expression_matrix, ''),
               "# Quantile Normalized: YES",
               paste("# OUTCOME: ", outcome, ''),
               paste('# Analysis Type: ', analysis_type, ''),
               '# COVARIATES: NONE', 
               paste('# OUTPUT FILE: ', filename, ''),
               paste('# NSAMPLES USED: ', n_samples_used, ''),
               '----------------------------------------------------------------------------------'
    )
    
    writeLines(lines, probe_file)
    close(probe_file)
    
    # ----- step 2: Fit Cox regression at probe level -----
    source('code/survival_analysisV2.R')
    # source('code/survival_analyisV2.R')
    
    data <- importData(cancer, accession, platform, TRUE, TRUE)
    s = Sys.time()
    output <- generateCoxFit_ProbeLevel(data, 1)
    print(difftime(Sys.time(), s, units = 'secs'))
    
    # add it to file
    write.table(output, file= paste('www/PRECOG_probe_level/', filename, sep =''), append=TRUE)

  }
}

# study_id <- paste(cancer, subtype, accession, platform, sep = '.')
# study_id <- str_replace_all(study_id, ' ', '.')
# study_id <- str_replace_all(study_id, '\\.\\.', '.')
# 
# gene_expression_matrix = paste('precogdata/Data/', study_id, '.data.RData', sep ='')
# outcome = 'OS'
# analysis_type = 'Univariate COX Continuous'
# output_file = ''
# n_samples_used = getSampleSize(cancer, accession, platform)
# 
# # ----- step 1: Generate file header -----
# filename <- paste(study_id, outcome, 'unisurv_cns.txt',sep = '.')
# probe_file <- file(paste('www/PRECOG_probe_level/', filename, sep =''))
# lines <- c(paste('# DISEASE: ', cancer, ''), 
#            paste("# SUBSET: ", subtype, ''),
#            paste('# SERIES: ', accession, ''),
#            paste('# PLATFORM: ', platform, ''),
#            paste('# DataMatrix: ', gene_expression_matrix, ''),
#            "# Quantile Normalized: YES",
#            paste("# OUTCOME: ", outcome, ''),
#            paste('# Analysis Type: ', analysis_type, ''),
#            '# COVARIATES: NONE', 
#            paste('# OUTPUT FILE: ', filename, ''),
#            paste('# NSAMPLES USED: ', n_samples_used, ''),
#            '----------------------------------------------------------------------------------'
#            )
# 
# writeLines(lines, probe_file)
# close(probe_file)
# 
# # ----- step 2: Fit Cox regression at probe level -----
# source('code/survival_analysisV2.R')
# # source('code/survival_analyisV2.R')
# 
# data <- importData(cancer, accession, platform, TRUE, TRUE)
# s = Sys.time()
# output <- generateCoxFit_ProbeLevel(data, 1)
# print(difftime(Sys.time(), s, units = 'secs'))
# 
# # add it to file
# write.table(output, file= paste('www/PRECOG_probe_level/', filename, sep =''), append=TRUE)


