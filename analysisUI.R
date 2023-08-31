# library(dplyr)
# ----- PRECOG -----
# ----- IMPORT AND PREPROCESS THE DATA 
# metaz <- read.csv("precogdata/metaz0.csv", header = T)
metaz <- read.csv("precogdata/PRECOG_V2/metaz0_V2.csv", header = T)
metaz <- metaz[,-c(1,3)]
colnames(metaz) <- str_replace_all(colnames(metaz), "_"," ")
colnames(metaz)[1] <- "Gene Symbol"
precog_info <- read.csv("precogdata/PRECOG_V2/precog_info_V2.csv", header = T)
precog_info$Cancer <- as.character(precog_info$Cancer)
trim.trailing <- function (x) sub("\\s+$", "", x)
precog_info$Cancer <- trim.trailing(precog_info$Cancer)
precog_info$Accession <- as.character(precog_info$Accession)
# # ----- get rid of $Name 
# metaz <- metaz[, -2]
# load('data/precog.RData')
# Save the cancer list
cancer_list <- colnames(metaz)

# Gene list 

# Adjust the names of the columns so we can have small box in the DT. 
# FIXME: find a better way to dit
# colnames(metaz)[4] <- "Adrenoc ortical cancer"
# colnames(metaz)[6] <- "Brain cancer Astro cytoma"
# colnames(metaz)[7] <- "Brain cancer Gliobl astoma"
# colnames(metaz)[9] <- "Brain cancer Medullo blastoma"
# colnames(metaz)[10] <- "Brain cancer Menin gioma"

# TEXTS
analysis1_text <-"Survival z-scores are collapsed by cancer/cancer subtype as described in
<a href=http://www.nature.com/doifinder/10.1038/nm.3909>Gentles/Newman et al.</a>.
You can download the tab-delimited (PCL) file, or view the table below. The survival scores 
for individual studies are also <a href=#shiny-tab-home data-toggle=tab data-value=home>available</a>.
False discovery rates corresponding to the meta-Z scores are available here, c
alculated by the method of Storey and Tibshirani (2003).Filter the table by gene symbol based on search terms (case insensitive) in the search field below. 
If you wish to search multiple terms, separate each term with a pipe | (do not add any spaces). 
For example, searching with <font color=\"#FF0000\"><b> foxm1|klrb1 </b></font> should yield 2 
results."

analysis2_text <- "Click on the number within any given cell to display the corresponding Kaplan Meier plot(s) for the 
                   corresponding gene and dataset(s). 
                    KM plots were generated using a median split."


# SHINY COMPONENTS
## UI
analysisItem <- tabItem(
  tabName = "metaz",
  fluidPage(
    includeCSS('code/css/analysis.css'),
    # theme = shinytheme("cosmo"),
    h1("Visualize meta-Z analysis across PRECOG", id = 'analysis_title'),
    box(
    collapsible = T,
    status = 'info', 
    # background = "light-blue",
    width = 30,
    htmlOutput("analysis1"), 
    br(), 
    box(textOutput("analysis2"), width = 30, height = 90,  background = "red"),
    
    # br(), 
    h5('Columns can be sorted by clicking on the cancer name in the table header.'),
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
#         tags$head(tags$style("table.dataTable thead .sorting {
#     /* background-image:url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABMAAAATCAQAAADYWf5HAAAAkElEQVQoz7XQMQ5AQBCF4dWQSJxC5wwax1Cq1e7BAdxD5SL+Tq/QCM1oNiJidwox0355mXnG/DrEtIQ6azioNZQxI0ykPhTQIwhCR+BmBYtlK7kLJYwWCcJA9M4qdrZrd8pPjZWPtOqdRQy320YSV17OatFC4euts6z39GYMKRPCTKY9UnPQ6P+GtMRfGtPnBCiqhAeJPmkqAAAAAElFTkSuQmCC) */
# }")),
        htmlOutput("zscore.conversion")
        # br()
      ),

      # box(
      #   title = 'Number of rows in table',
      #   status = 'info', 
      #   collapsible = T,
      #   collapsed = T, 
      #   width = 4,
      #   selectInput('numrows1', 'Select number of rows to display', 
      #               c(1,2,10,11,12,13,14,15,16,17,18,19,20), 
      #               selected = 10, 
      #               multiple = FALSE
      #               # selectize = TRUE, width = NULL, size = NULL
      #               )
      #   # br()
      # ),
    
    ),
    

    dataTableOutput("metaz"),
# style = "overflow-y: scroll;"

    )

  )
  

headerCallback <- c(
  "function(thead, data, start, end, display){",
  "  var $ths = $(thead).find('th');",
  "  $ths.css({'height' : '250px', 'width': '21px', 'vertical-align': 'bottom', 'white-space': 'nowrap'});",
  "  var betterCells = [];",
  "  $ths.each(function(){",
  "    var cell = $(this);",
  "    var newDiv = $('<div>', {height: '190px'});",
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

# ------ Server outputs 
p = dim(metaz)[2]
analysis1 <- renderText(analysis1_text)
analysis2 <- renderText(analysis2_text)

zscore.output <- renderUI(HTML('<img src="zscore.svg", height="50px, width = 50px" style="vertical-align:middle"/>'))

# brks <- quantile(metaz[,c(2:p)], probs = seq(.05, .95, .05), na.rm = TRUE)
# clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>% {paste0("rgb(255,", ., ",", ., ")")}

m_brks <- quantile(metaz[,c(2:ncol(metaz))], probs = seq(.01, .99, .01), na.rm = TRUE)
m_clrs_a <- round(seq(255, 40, length.out = (length(m_brks) + 1) / 2) , 0) %>% {paste0("rgb(", 295 - . ,",", 295 - ., ",", 255, ")")}
m_clrs_b <- round(seq(255, 40, length.out = (length(m_brks) + 1) / 2) , 0) %>% {paste0("rgb(", 255 ,",", ., ",", ., ")")}
m_clrs <- c(m_clrs_a, m_clrs_b)

metazOutput <- renderDataTable(
                        datatable(
                        metaz,
                        style = 'default', 
                        selection = list(mode="single", target="cell", selected = NULL),
                        escape = F,
                        rownames = F,
                        filter = list(positiomn = "top", clear = T, plain = T),
                        options = list(scrollX = T,
                                       scrollY = T, 
                                       headerCallback = JS(headerCallback), 
                                       search = list(regex = TRUE), 
                                       pageLength = 20, 
                                       lengthMenu =  list(c(1, 10, 20, 50, 100, 300, 500), c('1', '10', '20', '50', '100', '300', '500')),
                                       autoWidth = F, 
                                       columnDefs = list(list(className = 'dt-center', width = '30px', targets = "_all")
                                                          )
                                      ),
                        ) %>% formatStyle(names(metaz[,c(2:p)]), backgroundColor = styleInterval(m_brks, m_clrs))
                          %>% formatRound(columns = colnames(metaz)[-c(1)], digits=2)
                        )
                



# ----- iPRECOG (Signature matrix) -----

# ui
ianalysis_text1 <- renderText('The iPRECOG signature matrix, termed LM22, identifies the genes that are specific to different immune cell types. Expression levels of these genes are presented visually below and are also available in a <a href = \'precogdata/iprecog/iPRECOG_sigmatrix.pcl\'>PCL </a> file. More information about LM22 is available in the CIBERSORT publication (Nature Methods 2015).')
ianalysisItem <- tabItem(
  tabName = "signature",
  fluidPage(
    # includeCSS('code/css/analysis.css'),
    h1('iPRECOG Signature Matrix - LM22', id = 'ianalysis_title'),
    box(
      collapsible = T,
      status = 'info', 
      # background = "light-blue",
      width = 30,
      htmlOutput("ianalysis1"), 
      # br(), 
      # box(textOutput("analysis2"), width = 30, height = 90,  background = "red"),
      br(), 
    ),
    br(), 
    dataTableOutput("signature_matrix"),
  )
)
# server
sig_matrix <- read.csv('precogdata/iprecog/signature_matrix.csv', header =T, sep = ',')[, -1]
ncol_sig = ncol(sig_matrix)
# sig_brks_2 <- quantile(sig_matrix[,c(2:ncol_sig)], probs = seq(.01, .99, .01), na.rm = TRUE)
# sig_clrs_2 <- round(seq(255, 40, length.out = length(sig_brks_2) + 1), 0) %>% {paste0("rgb(", 295 - . ,",", 40, ",", ., ")")}
s_brks <- quantile(sig_matrix[,c(2:ncol(sig_matrix))], probs = seq(.01, .99, .01), na.rm = TRUE)
s_clrs_a <- round(seq(255, 40, length.out = (length(s_brks) + 1) / 2) , 0) %>% {paste0("rgb(", 295 - . ,",", 295 - ., ",", 255, ")")}
s_clrs_b <- round(seq(255, 40, length.out = (length(s_brks) + 1) / 2) , 0) %>% {paste0("rgb(", 255 ,",", ., ",", ., ")")}
s_clrs <- c(s_clrs_a, s_clrs_b)
iprecog_sig_matrix_Output <- renderDataTable(
  datatable(
    sig_matrix,
    style = 'default', 
    selection = list(mode="single", target="cell", selected = NULL),
    escape = F,
    rownames = F,
    filter = list(positiomn = "top", clear = T, plain = T),
    options = list(scrollX = T,
                   scrollY = T,
                   headerCallback = JS(headerCallback),
                   search = list(regex = TRUE), 
                   pageLength = 20 , 
                   lengthMenu =  list(c(1, 10, 20, 50, 100, 300, 500), c('1', '10', '20', '50', '100', '300', '500')),
                   autoWidth = F, 
                   columnDefs = list(list(className = 'dt-center', width = '30px', targets = "_all")
                   )
    ),
  ) %>% formatStyle(names(sig_matrix[,c(2:ncol_sig)]), backgroundColor = styleInterval(s_brks, s_clrs))
  # ) %>% formatStyle(names(sig_matrix, 
  #                   width='25%', backgroundColor = styleInterval(seq(-1, 1, 0.25), rev(brewer.pal(n = 10, name = "RdBu")))))
  %>% formatRound(columns = colnames(sig_matrix)[-c(1)], digits=2)
)

# ----- iPRECOG (MetaZ table) -----
# ui
iMetazItem <- tabItem(
  tabName = "iMetaz",
  fluidPage(
    # includeCSS('code/css/analysis.css'),
    h1('iPRECOG meta-Z matrix'),
    box(
      collapsible = T,
      status = 'info', 
      # background = "light-blue",
      width = 30,
      htmlOutput("iMetaz1"), 
      br(), 
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
        htmlOutput("izscore.conversion")
        # br()
      )
    ),
    dataTableOutput("iMetaz_table"),
  )
)
# server
iprecog_metaz_text <- renderText('Inferred fractions of immune cell types can be associated with survival outcomes in a meta-Z analogous to the individual gene meta-Zs found in PRECOG. Z-score associations of these CIBERSORT-computed fractions are presented visually below and are also available in an <a href = \'precogdata/iprecog/iPRECOG_metaZ.xlsx\'> Excel file </a>.')
metaz_table <- read.csv('precogdata/iprecog/iPRECOG_metaZNEW.csv', header =T, sep = ',')[, -1]
metaz_table$Cancer.type <- as.character(metaz_table$Cancer.type)
sig_brks <- quantile(metaz_table[,c(4:ncol(metaz_table))], probs = seq(.01, .99, .01), na.rm = TRUE)
sig_clrs_a <- round(seq(255, 40, length.out = (length(sig_brks) + 1) / 2) , 0) %>% {paste0("rgb(", 295 - . ,",", 295 - ., ",", 255, ")")}
sig_clrs_b <- round(seq(255, 40, length.out = (length(sig_brks) + 1) / 2) , 0) %>% {paste0("rgb(", 255 ,",", ., ",", ., ")")}
sig_clrs <- c(sig_clrs_a, sig_clrs_b)
iprecog_metaz_table_Output <- renderDataTable(
  datatable(
    metaz_table,
    style = 'default', 
    selection = list(mode="single", target="cell", selected = NULL),
    # filter = 'top',
    escape = F,
    rownames = F,
    # filter = list(position = "top", clear = T, plain = T),
    options = list(scrollX = T,
                   scrollY = T,
                   headerCallback = JS(headerCallback),
                   search = list(regex = TRUE), 
                   pageLength = 20 , 
                   lengthMenu =  list(c(1, 10, 20, 50, 100, 300, 500), c('1', '10', '20', '50', '100', '300', '500')),
                   autoWidth = TRUE, 
                   columnDefs = list(
                                     list(className = 'dt-center', width = '30px', targets = c(1:24)),
                                     list(width = '200px', targets = 0))
    ),
  ) %>% formatStyle(names(metaz_table[,c(4:ncol(metaz_table))]), backgroundColor = styleInterval(sig_brks, sig_clrs))
  # ) %>% formatStyle(names(sig_matrix, 
  #                   width='25%', backgroundColor = styleInterval(seq(-1, 1, 0.25), rev(brewer.pal(n = 10, name = "RdBu")))))
  %>% formatRound(columns = colnames(metaz_table)[-c(1:3)], digits=2)
)

# ----- iPRECOG (immune fraction table) -----
# ui
iPRECOG_immune_fraction_Item <- tabItem(
  tabName = "immune-fraction",
  fluidPage(
    # includeCSS('code/css/analysis.css'),
    h1('iPRECOG cancer-specific immune fractions'),
    box(
      collapsible = T,
      status = 'info', 
      # background = "light-blue",
      width = 30,
      htmlOutput("immune_fraction_text"), 
      br(), 
    ),
    br(),
    dataTableOutput("immune_fraction_table"),
  )
)
# server
immune_fraction_text_output <- renderText('Relative RNA fractions of 22 leukocyte subsets in human cancers. Results from applying CIBERSORT to 5,782 tumor specimens, selected from 68 PRECOG cancer data sets profiled on Affymetrix HGU133 platforms. For clarity, values are presented as the mean fraction for each cell population and cancer type. In accordance with our sensitivity and specificity analysis, a CIBERSORT p-value threshold of 0.005 was employed to filter out insignificant deconvolution results. Results also available in <a href = \'precogdata/iprecog/iPRECOG_cancer-immunefractions.xlsx\'> Excel </a> format.')
immune_fraction_dat <- read.csv('precogdata/iprecog/iPRECOG_cancer-immunefractionsNEW.csv', header = T, sep =',')[, -1]
# metaz_table$Cancer.type <- as.character(metaz_table$Cancer.type)
i_brks <- quantile(immune_fraction_dat[,c(5:ncol(immune_fraction_dat) - 2)], probs = seq(.01, .99, .01), na.rm = TRUE)
# i_clrs <- round(seq(255, 40, length.out = length(i_brks) + 1), 0) %>% {paste0("rgb(", . ,",", 40, ",", 295 - ., ")")}
# i_brks <- quantile(immune_fraction_dat[,c(5:ncol(immune_fraction_dat))], probs = seq(.01, .99, .01), na.rm = TRUE)
i_clrs_a <- round(seq(255, 40, length.out = (length(i_brks) + 1) / 2) , 0) %>% {paste0("rgb(", 295 - . ,",", 295 - ., ",", 255, ")")}
i_clrs_b <- round(seq(255, 40, length.out = (length(i_brks) + 1) / 2) , 0) %>% {paste0("rgb(", 255 ,",", ., ",", ., ")")}
i_clrs <- c(i_clrs_a, i_clrs_b)
iprecog_immune_fraction_table_Output <- renderDataTable(
  datatable(
    immune_fraction_dat,
    style = 'default', 
    selection = list(mode="single", target="cell", selected = NULL),
    # filter = 'top',
    escape = F,
    rownames = F,
    # filter = list(position = "top", clear = T, plain = T),
    options = list(scrollX = T,
                   scrollY = T, 
                   headerCallback = JS(headerCallback),
                   search = list(regex = TRUE), 
                   pageLength = 20 , 
                   lengthMenu =  list(c(1, 10, 20, 50, 100, 300, 500), c('1', '10', '20', '50', '100', '300', '500')),
                   autoWidth = TRUE, 
                   columnDefs = list(
                     list(className = 'dt-center', width = '30px', targets = c(3:26)),
                     list(width = '100px', targets = 0))
    ),
  ) %>% formatStyle(names(immune_fraction_dat[,c(5:ncol(immune_fraction_dat) - 2)]), backgroundColor = styleInterval(i_brks, i_clrs))
  # ) %>% formatStyle(names(sig_matrix, 
  #                   width='25%', backgroundColor = styleInterval(seq(-1, 1, 0.25), rev(brewer.pal(n = 10, name = "RdBu")))))
  %>% formatRound(columns = colnames(immune_fraction_dat)[-c(1, 2, 3, 27, 28)], digits=2)
)
