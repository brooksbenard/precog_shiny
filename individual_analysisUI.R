##############################################################
###### UI for Individual Analysis Panels (PRECOG & TCGA) #####
##############################################################
# ----- Import TCGA individual studies result -----
TCGA.metaz <- read.csv('precogdata/TCGA/TCGA.metaz.csv', header = T, sep =',')
TCGA.metaz <- TCGA.metaz[, -1]
TCGA.metaz$Name <- as.character(TCGA.metaz$Name)
# ----- Import PRECOG individual studies result -----
# precog.indiv <- read.csv('precogdata/precog_ind_score.csv', header = T) # V1
precog.indiv <- read.csv('precogdata/PRECOG_V2/precog_ind_score_V2.csv', header = T) # V2
precog.indiv <- precog.indiv[, -c(1, 3)] # get rid of gene name
colnames(TCGA.metaz)[1] <- 'Name'

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

TCGA.cancer.id.list = colnames(TCGA.metaz)[-c(1)]
PRECOG.cancer.list = colnames(precog.indiv)[-c(1:2)]

# ------ UI -----
individual_analysisItem <- tabItem(
  
  tabName = "individual_analysis",
  
  fluidPage(
    
    h1("Individual Datasets Analysis"),
    
    tabsetPanel(type = "pills",
                id = "ind_analysis",
                # ----- PRECOG individual -----
                tabPanel(title = "PRECOG",value = "panel1" ,
                         h1("Visualize Z-scores across individual PRECOG datasets"),
                         box(
                           collapsible = T,
                           status = 'info', 
                           width = 30,
                           htmlOutput("precog.ind.analysis1"),
                           br(),
                           box(textOutput("precog.ind.analysis2"), width = 30, height = 90,  background = "red"),
                           
                           br()
                         ),
                         
                         
                         br(),
                         fluidRow(
                           box(
                             title = 'Z-score to P-value',
                             status = 'info', 
                             collapsible = T,
                             collapsed = T, 
                             width = 4,
                             # height = 10, 
                             htmlOutput("zscore.conversion.precog")
                             # br()
                           )
                         ),
                         # DTOutput("precog.zscore")
                         dataTableOutput('precog.zscore')
                ),
                       
                # ----- TCGA individual -----
                tabPanel("TCGA", value = "panel2" ,
              
                         h1("Visualize survival Z-scores across TCGA RNA-Seq datasets"),
                     
                         # fluidRow(                         
                         box(
                           collapsible = T,
                           status = 'info', 
                           width = 30,
                           htmlOutput("TCGAanalysis1"),
                           br(),
                           box(textOutput("TCGAanalysis2"), width = 30, height = 90,  background = "red"),

                           br()
                         ),
                         
                         
                         br(),
                         fluidRow(
                         box(
                           title = 'Z-score to P-value',
                           status = 'info', 
                           collapsible = T,
                           collapsed = T, 
                           width = 4,
                           # height = 10, 
                           htmlOutput("zscore.conversion.tcga")
                           # br()
                         )
                         ),
                         # DTOutput("TCGA.metaz")
                         dataTableOutput('TCGA.metaz')
                )
                ################
    )
    
  )
  
)
# Callback function (for formatting datatable headers)
headerCallback <- c(
  "function(thead, data, start, end, display){",
  "  var $ths = $(thead).find('th');",
  "  $ths.css({'height' : '220px', 'width': '21px', 'vertical-align': 'bottom', 'white-space': 'nowrap'});",
  "  var betterCells = [];",
  "  $ths.each(function(){",
  "    var cell = $(this);",
  "    var newDiv = $('<div>', {height: '90px'});",
  "    var newInnerDiv = $('<div>', {text: cell.text()});",
  "    newDiv.css({margin: 'auto'});",
  "    newInnerDiv.css({",
  "      'transform-origin': 'top right',",
  "      'transform': 'translate(-25px, 80px) rotate(222deg)',",
  "      'writing-mode': 'tb-rl',",
  "      'white-space': 'nowrap'",
  "    });",
  "    newDiv.append(newInnerDiv);",
  "    betterCells.push(newDiv);",
  "  });",
  "  $ths.each(function(i){",
  "    $(this).html(betterCells[i]);",
  "  });",
  "}")
headerCallback.TCGA <- c(
  "function(thead, data, start, end, display){",
  "  var $ths = $(thead).find('th');",
  "  $ths.css({'height' : '90px', 'width': '21px', 'vertical-align': 'bottom', 'white-space': 'nowrap'});",
  "  var betterCells = [];",
  "  $ths.each(function(){",
  "    var cell = $(this);",
  "    var newDiv = $('<div>', {height: '90px'});",
  "    var newInnerDiv = $('<div>', {text: cell.text()});",
  "    newDiv.css({margin: 'auto'});",
  "    newInnerDiv.css({",
  "      'transform-origin': 'top right',",
  "      'transform': 'translate(-25px, 80px) rotate(222deg)',",
  "      'writing-mode': 'tb-rl',",
  "      'white-space': 'nowrap'",
  "    });",
  "    newDiv.append(newInnerDiv);",
  "    betterCells.push(newDiv);",
  "  });",
  "  $ths.each(function(i){",
  "    $(this).html(betterCells[i]);",
  "  });",
  "}")
# ----- Server outputs ------ 
# ----- PRECOG -----
analysis1.text.precog <- renderText("Survival z-scores for each individual dataset in PRECOG 
are shown. Column headers include a numerical PRECOG ID which crossreferences the individual 
dataset (including citation information) in the datasets table. In the alternative meta-Z analysis, 
they are collapsed by cancer/cancer subtype as described in <a href=http://www.nature.com/doifinder/10.1038/nm.3909>Gentles/Newman et al</a>. 
FDRs for individual studies can be found in the Downloads.
Download the tab-delimited (PCL) file, or view the table below.
Filter the table based on search terms (case insensitive) in the search field below. If you wish to search 
multiple terms, separate each term with a pipe | (do not add any spaces). 
For example, searching with  <font color=\"#FF0000\"><b> foxm1|klrb1 </b></font>  should yield 2 results")

analysis2.text.precog <- renderText("Click on the number within any given cell to display the corresponding Kaplan Meier plot(s) for the 
                   corresponding gene and dataset(s). 
                    KM plots were generated using a median split.")

zscore.output.precog <- renderUI(HTML('<img src="zscore.svg", height="50px, width = 50px" style="vertical-align:middle"/>'))
ncol = dim(precog.indiv)[2]
brks <- quantile(precog.indiv[,c(2:ncol)], probs = seq(.05, .95, .05), na.rm = TRUE)
# clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>% {paste0("rgb(255,", ., ",", ., ")")}
clrs_a <- round(seq(255, 40, length.out = (length(brks) + 1) / 2) , 0) %>% {paste0("rgb(", 295 - . ,",", 295 - ., ",", 255, ")")}
clrs_b <- round(seq(255, 40, length.out = (length(brks) + 1) / 2) , 0) %>% {paste0("rgb(", 255 ,",", ., ",", ., ")")}
clrs <- c(clrs_a, clrs_b)

precog.zscore.Output <- renderDataTable(
  datatable(
    precog.indiv,
    style = 'default', 
    selection = list(mode="single", target="cell", selected = NULL), 
    escape = F,
    rownames = F,
    filter = list(positiomn = "top", clear = T, plain = T),
    options = list(scrollX = T,
                   headerCallback = JS(headerCallback),
                   search = list(regex = TRUE), 
                   pageLength = 20 , 
                   lengthMenu =  list(c(1, 10, 20, 50, 100, 300, 500), c('1', '10', '20', '50', '100', '300', '500')),
                   autoWidth = T, 
                   columnDefs = list(list(width = '30px', targets = "_all")
                   )
    ),
  ) %>% formatStyle(names(precog.indiv[,c(2:ncol)]), backgroundColor = styleInterval(brks, clrs))
  %>% formatRound(columns = colnames(precog.indiv)[-c(1)], digits=2)
)

# ----- TCGA -----
analysis1.text.TCGA <-renderText("Survival z-scores are collapsed by cancer/cancer subtype as described in
<a href=http://www.nature.com/doifinder/10.1038/nm.3909>Gentles/Newman et al.</a>.
You can download the tab-delimited (PCL) file, or view the table below. The survival scores 
for individual studies are also <a href=#shiny-tab-home data-toggle=tab data-value=home>available</a>.
False discovery rates corresponding to the meta-Z scores are available here, c
alculated by the method of Storey and Tibshirani (2003).Filter the table by gene symbol based on search terms (case insensitive) in the search field below. 
If you wish to search multiple terms, separate each term with a pipe | (do not add any spaces). 
For example, searching with <font color=\"#FF0000\"><b> foxm1|klrb1 </b></font> should yield 2 
results.")

analysis2.text.TCGA <- renderText("Click on the number within any given cell to display the corresponding Kaplan Meier plot(s) for the 
                   corresponding gene and dataset(s). 
                    KM plots were generated using a median split.")



zscore.output.tcga <- renderUI(HTML('<img src="zscore.svg", height="50px, width = 50px" style="vertical-align:middle"/>'))
# zscore.output <- renderUI(HTML('<img src="zscore.svg", style = "width:50px;height:50px;"/>'))
# margin:0px 0px
# style="vertical-align:middle;margin:0px 40px"/>
brks_tcga <- quantile(TCGA.metaz[,c(2:36)], probs = seq(.05, .95, .05), na.rm = TRUE)
# clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>% {paste0("rgb(255,", ., ",", ., ")")}
clrs_a_tcga <- round(seq(255, 40, length.out = (length(brks_tcga) + 1) / 2) , 0) %>% {paste0("rgb(", 295 - . ,",", 295 - ., ",", 255, ")")}
clrs_b_tcga <- round(seq(255, 40, length.out = (length(brks_tcga) + 1) / 2) , 0) %>% {paste0("rgb(", 255 ,",", ., ",", ., ")")}
clrs_tcga <- c(clrs_a_tcga, clrs_b_tcga)

TCGA.metazOutput <- renderDataTable(
  datatable(
    TCGA.metaz,
    style = 'default', 
    selection = list(mode="single", target="cell", selected = NULL), 
    escape = F,
    rownames = F,
    filter = list(positiomn = "top", clear = T, plain = T),
    options = list(scrollX = T,
                   headerCallback = JS(headerCallback.TCGA),
                   search = list(regex = TRUE), 
                   pageLength = 20 , 
                   lengthMenu =  list(c(1, 10, 20, 50, 100, 300, 500), c('1', '10', '20', '50', '100', '300', '500')),
                   autoWidth = T, 
                   columnDefs = list(list(width = '30px', targets = "_all")
                   )
    ),
  ) %>% formatStyle(names(TCGA.metaz[,c(2:36)]), backgroundColor = styleInterval(brks_tcga, clrs_tcga))
    %>% formatRound(columns = colnames(TCGA.metaz)[-c(1)], digits=2)
)