# ----- Rscript for dashboard panel: data -----
# ----- Import precog overview data -----
precog <- read.csv("precogdata/PRECOG_V2/precog_V2.csv", header = T, sep = ",")[, -c(1)]
original_precog <- read.csv('precogdata/PRECOG_V2/original_precog_V2.csv', header = T, sep = ',')[, -1]
# ----- UI items -----
# texts
intro_data <- "Download all annotation tables as a"
text_data <- "Please note that for studies performed on multiple array platforms, each platform 
              is analyzed separately and then combined into the final meta-Z score. Hence each 
              platform within a study has a separate line below."

# datatable 
dataItem <- tabItem(
  
  tabName = "data",
  fluidPage(
    includeCSS('code/css/data.css'),
    fluidRow(
      uiOutput("text_data"),
      br(), 
      downloadButton("downloadAnnotations", label = "Download Annotations Tables")
    ),
    br(), 
    column(
      width = 12,
      dataTableOutput('precog_data')
    )
  )
)

# ----- Creating Data outputs -----
data_output <- DT::renderDataTable(
                               datatable(precog,
                               escape = F,
                               rownames = F, 
                               filter = 'top', 
                               options = list(
                                 autoWidth = TRUE,
                                 scrollX = T, 
                                 pageLength = 20,
                                 columnDefs = list(list(width = '150px', targets = c(1)),
                                                   list(width = '85px', targets = c(8:9)))
                               )
                               ) %>% formatStyle(
                                 colnames(precog),
                                 target = 'row',
                                 backgroundColor = styleEqual(c(53:60), rep('lightyellow', 8))
                               )
)

annotations_download <- downloadHandler(
                        filename = "AnnotationsTable.zip",
                        content = function(file){
                          fs = c("www/precogdata/filtered_clinical_annotations")
                          zip(zipfile = file, files = fs)
                        }
                        )