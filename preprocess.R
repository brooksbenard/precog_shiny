source("utils.R")

##########################################################
#################### PRECOG MAIN DATA ####################
##########################################################

precog <- read.csv("data/precog.csv", header = T, sep = "\t")[c(-1,-2),]
precog$Pubmed.ID <- as.character(precog$Pubmed.ID)
precog$Accession <- as.character(precog$Accession)
precog$Platform <- as.character(precog$Platform)
precog$Annotation <- as.character(precog$Annotation)
precog$Subtype <- as.character(precog$Subtype)
precog$Cancer <- as.character(precog$Cancer)

precog_info <- precog[,c(2,6)]
precog_info$Accession_num <- 0
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

write.csv(precog_info, "data/precog_info.csv")

###################################### METAZ Data ######################################

metaz0 <- read.csv("data/metaz.csv", header = T ,sep = ";")
metaz0$Gene <- as.character(metaz0$Gene)
metaz0$Adrenocortical_cancer <- as.character(metaz0$Adrenocortical_cancer)
metaz0$Name <- as.character(metaz0$Name)
gene_list <- metaz0$Gene

# Mapping idx to (gene, cancer)
# It is used for the customized KM display when there is a click on the MetaZ DT.

# e <- new.env()
# e$my_key <- 10
# ls(e)
# 
# idx_GC <- vector(mode = "list", length = 23287)
# names(idx_GC) <-c(1:23287)
# col_idx = 4
# row_idx = 1
# for (i in c(1:23287)){
#   gene = as.character(metaz0$Gene[row_idx])
#   cancer = as.character(colnames(metaz0)[col_idx])
#   # print(coord)
#   coord = c(row_idx, col_idx - 1)
#   idx_GC[[coord]] = list(gene, cancer)
#   if (row_idx == 23287){
#     col_idx = col_idx + 1
#     row_idx = 1
#   }
#   else {row_idx = row_idx+1}
# }
# save(idx_GC, file = "data/idx2GC.rdata")


# Modify the values of the DT to handle the display of the KM plots.
# style_1 = " type=\"button\" style=\"background-color:IndianRed;color:black;width:300px;height:80px\" \"> "
# ref_link = "<a id=\"button\" href=\"#\" class=\"action-button\" onclick=\"Shiny.onInputChange(&quot;select_button&quot;,  this.id)"
# col_idx = 4
# row_idx = 1
# for (i in c(1:23287)){
#   link = str_replace(ref_link, "button", paste("button_", as.character(i), sep = ""))
#   link = paste(link, style_1, sep = "")
#   new_val = paste(link, as.character(metaz0[row_idx,col_idx]), "</a>", sep ="")
#   metaz0[row_idx, col_idx] = as.character(new_val)
#   if (row_idx == 23287){
#     col_idx = col_idx + 1
#     row_idx = 1
#   }
#   else {row_idx = row_idx+1}
# }

# HANDLE NA 

# Mapping a cancer name to its list of accessions : 


# Gene Symbol hyperlink
i = 1
for (a in metaz0$Gene){
  
  link = paste("https://www.ncbi.nlm.nih.gov/gene/?term=",a,"[gene]+AND+%22human%22[orgn]", sep = "")
  metaz0$Gene[i] <- as.character(ToLink(a, link))
  i = i+1

}

# Handle `Name` length
i=1
for(name in  metaz0$Name){
  metaz0$Name[i] <- reduce_char_lenght(metaz0$Name[i],15)
  i = i+1
} 


write.table(metaz0, "data/metaz0.txt", sep = "\t")
write.csv(metaz0, "data/metaz0.csv")

