setwd('/srv/shiny-server/precog')
source('code/survival_analysisV2.R')
# args <- commandArgs(TRUE)
# cancer = as.character(args[1])
# accession = as.character(args[2])
# platform = as.character(args[3])
# 
# if (is.na(cancer) || is.na(accession) || is.na(platform)){
#   print('----- Usage: Rscript computeKMplots -cancer -accession -platform')
#   stop('Wrong or missing arguments')
# }

# print(cancer)
# print(accession)
# print(platform)


# step 1: import data
precog <- read.csv('precogdata/PRECOG_V2/original_precog_V2.csv', header = T, ',')[, -1]
for (i in 1:dim(precog)[1]){
  cancer = as.character(precog$Cancer[i])
  if (cancer == 'Head and neck cancer'){
    cancer = 'Head.and.neck.cancer'
    accession = as.character(precog$Accession[i])
    platform = as.character(precog$Platform[i])
    print(cancer)
    print(accession)
    print(platform)
    cat('Importing data before computing plots.... \n')
    data <- importData(cancer, accession, platform, TRUE, TRUE)
    cat('Data import completed \n ')
    cat('Sanity checks passed! \n')
    cat('Start computing KMplots...')
    s = Sys.time()
    generateKMPlots(data, TRUE, 1, 23287)
    print(difftime(Sys.time(), s, units = 'secs'))
    cat('Plots are computed! \n')
  }
}
# cat('Importing data before computing plots.... \n')
# data <- importData(cancer, accession, platform, TRUE, TRUE)
# cat('Data import completed \n ')
# cat('Sanity checks passed! \n')
# cat('Start computing KMplots...')
# s = Sys.time()
# generateKMPlots(data, TRUE, 1, 23287)
# print(difftime(Sys.time(), s, units = 'secs'))
# cat('Plots are computed! \n')
