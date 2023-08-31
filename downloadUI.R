# ----- UI for Download Panel -----
downloadItem <- tabItem(tabName = 'downloads',
                        value = 'Download',
                        h4('Please find below old PRECOG (and TCGA) data. New studies may have been added to current tables.'),
                        fluidPage(
                          fluidRow(
                            column(
                              width = 4,
                              h1(HTML('<b> File </b>'))
                            ),
                            column(
                              width = 4,
                              h1(HTML('<b> Description </b>'))
                            ),
                            column(
                              width = 4, 
                              h1(HTML('<b> Size </b>'))
                            )
                          ),
                          
                          fluidRow(
                            column(
                              width = 4,
                              downloadButton('dwld1', 'Probe-level results'),
                              # tags$style(type='text/css', "#button { vertical-align- middle; }")
                            ),
                            column(
                              width = 4,
                              h4('Probe-level results (including FDRs) for individual studies/subsets in PRECOG,
prior to collapsing to gene symbols. See header of files for details')
                            ),
                            column(
                              width = 4, 
                              h4('1.28 GB')
                            )
                          ),
                          
                          fluidRow(
                            column(
                              width = 4,
                              downloadButton('dwld2', 'MetaZ', class = 'butt2')
                            ),
                            column(
                              width = 4,
                              h4('Matrix of Zscore in PRECOG')
                            ),
                            column(
                              width = 4, 
                              h4('7.6 MB')
                            )
                          ), 
                          
                          fluidRow(
                            column(
                              width = 4,
                              downloadButton('dwld3', 'Annotation Files', class = 'butt2')
                            ),
                            column(
                              width = 4,
                              h4('Overall survival times and outcomes for studies (0=ALIVE; DEAD=1 for censoring status)')
                            ),
                            column(
                              width = 4, 
                              h4('482 KB')
                            )
                          ), 
                          
                          fluidRow(
                            column(
                              width = 4,
                              downloadButton('dwld4', 'Individual Z', class = 'butt2')
                            ),
                            column(
                              width = 4,
                              h4('Matrix of gene-level Z-scores for individual studies/subset of studies in PRECOG')
                            ),
                            column(
                              width = 4, 
                              h4('23.8 MB')
                            )
                          ),
            
                          
                          fluidRow(
                            column(
                              width = 4,
                              downloadButton('dwld5', 'TCGA', class = 'butt2')
                            ),
                            column(
                              width = 4,
                              h4('Z-score matrix for TCGA RNA-seq datasets ordered by gene the same as the main PRECOG matrix')
                            ),
                            column(
                              width = 4, 
                              h4('12.2 MB')
                            )
                          )
                        )
)

dwld1_button <- downloadHandler(
  filename = "probe_level_results.zip",
  content = function(file){
    fs = c("www/PRECOG_probe_level")
    zip(zipfile = file, files = fs)
  }
)

dwld2_button <- downloadHandler(
  filename = "precog_metazscore.csv",
  content = function(file){
    write.csv('www/metaz0.csv',  file)
  }
)

dwld3_button <- downloadHandler(
  filename = "AnnotationsTable.zip",
  content = function(file){
    fs = c("www/precogdata/filtered.infofiles")
    zip(zipfile = file, files = fs)
  }
)

dwld4_button <- downloadHandler(
  filename = "precog_individual_zscore.csv",
  content = function(file){
    write.csv("www/precogdata/precog_ind_score.csv", file)
  }
)

dwld5_button <- downloadHandler(
  filename = "TCGA_metaz.csv",
  content = function(file){
    write.csv("www/precogdata/TCGA/TCGA_metaz.csv", file)
  }
)





