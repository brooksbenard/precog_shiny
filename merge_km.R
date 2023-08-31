library(stringr)
# ----- Option 1: merge `km_info` files -----
# km <- read.csv('data/km_info.csv', header = T, sep = ',')[, -1]
# km_s <- read.csv('precogdata/km_info.csv', header = T, sep = ',')[, -1]
# new_km <- merge(km_s, km, all.x = T)
# new_km[is.na(new_km)] = -1
# write.csv(new_km, 'data/km_info.csv')  

# ----- Option 2: go through plots name -----
# file to update KM_info file
# setwd('/Volumes/Thibaud/BMIR/precog/precog')
# on VM only
setwd('/srv/shiny-server/precog')

precog <- read.csv("precogdata/original_precog.csv", header = T, sep = ",")[, -1]
precog$Subtype <- as.character(precog$Subtype)
km <- read.csv('precogdata/km_info.csv', header = T, sep = ',')[, -1]

cancer = 'Head and neck cancer'
subtype = 'Oesophageal cancer'

subset = precog[which((precog$Cancer == cancer) & (precog$Subtype == subtype)),]
accessions = as.character(subset$Accession)
platforms =  as.character(subset$Platform)
outcomes = as.character(subset$Outcome)
cancer <- str_replace_all(cancer, ' ', '.')
for (i in 1:length(accessions)){
  accession = accessions[i]
  platform = platforms[i]
  platform <- str_replace_all(platform, '/', '_')
  outcome = outcomes[i]
  cat("------------------------ \n")
  cat("cancer: ", cancer, "\n")
  cat("subtype: ", subtype, "\n")
  cat("accession: ", accession, "\n")
  cat("platform: ", platform, "\n")
  folder = paste(tolower(paste(cancer, subtype, accession, platform, outcome, sep = '.')),
                 'KMplots/', sep = '.')
  folder <- str_replace_all(folder, '\\.\\.', '\\.')
  folder <- str_replace_all(folder, ' ', '\\.')
  id = str_replace_all(tolower(paste(cancer, subtype, accession, platform)), ' ', '\\.')
  id = str_replace_all(id, '\\.\\.', '\\.')
  plots_list <- system(paste('ls precogdata/KMplots/',
                             folder,
                             sep = ''), intern = T)
  cat('Study has', length(plots_list), 'KM plots \n')
  # Update `km_info` table
  cat('Let us update the `km_info` table... \n')
  if (length(plots_list) == 0){
    cat('No KM plot yet for this study! \n')
    next
  }
  if (!(id %in% colnames(km))){
    cat('Study with id:', id, 'not in file yet \n')
    km[[id]] = -1
  }
  s = Sys.time()
  for (i in 1:nrow(km)){
    gene = tolower(km[i, 1])
    if (TRUE %in% grepl(gene, plots_list)){
      km[i, ][id] = 1
    }
  }
  print(difftime(Sys.time(), s, units = 'secs'))
  cat('...Complete! \n')
}
cat('-------------------- \n')
cat('Exporting `km_info` table \n')
write.csv(km, 'precogdata/km_info.csv')
cat('-------------------- \n')
