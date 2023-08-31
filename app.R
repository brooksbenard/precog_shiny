#################
##### app.R #####
#################

# ----- Server settings -----
# options(shiny.port = 8080)
# options(shiny.host = "127.0.0.1")

suppressPackageStartupMessages({
    library(shiny)
    library(shinydashboard)
    library(shinydashboardPlus)
    # library(shinythemes)
    # library(shinyjs)
    library(shinyBS)
    library(DT)
    library(stringr)
    library(XML)
})


# load("data/idx2GC.rdata") #FIXME

# ----- R scripts ------
source('code/utils.R')
source('code/dashboardHeader.R')
source('code/dashboardSidebar.R')
source("code/homeUI.R")
source('code/dataUI.R')
source('code/analysisUI.R')
source('code/individual_analysisUI.R')
source('code/contactUI.R')
source('code/downloadUI.R')

# info about all studies used in precog
accession.main.list <- as.character(original_precog$Accession)
platform.main.list <- as.character(original_precog$Platform)
study.id.list <- str_replace(str_replace_all(tolower(paste(original_precog$Cancer, 
                                  original_precog$Subtype,
                                  original_precog$Accession, 
                                  original_precog$Platform,
                                  sep ='.')), ' ', '\\.'), '\\.\\.', '\\.')
# ----- ui -----
ui <- dashboardPagePlus(
    # shinyjs::useShinyjs(),
    skin = 'black',
    # ----- HEADER -----
    header,
    # ----- SIDEBAR -----
    Sidebar, 
    # ----- BODY -----
    dashboardBody(
        # IPRECOG_home,
        tabItems(
            # ----- see code/home.R -----
            home, 
            iPRECOG_home,
            # About tab content
            tabItem(tabName = "about",
                    value = 'about', 
                    h2("About")
            ),
            # ----- see code/download.R -----
            downloadItem,
            
            # ----- see code/data.R -----
            dataItem,
            # ----- see code/analysisUI.R -----
            analysisItem, 
            ianalysisItem, 
            iMetazItem,
            iPRECOG_immune_fraction_Item,
            # ----- see code/individual_analysisUI.R
            individual_analysisItem,
            
            contactItem
            
        ),
        

        
        # Modal window for terms of use
        # bsModal(id = "use_modal",
        #         title ="Terms of Use",
        #         trigger = "use_link",
        #         textOutput("terms_of_use"))
    ),

    # ----- footer ----- #TODO: get rid of it
    # footer = dashboardFooter(
    #             left_text = span(img(src = "stanford.svg", height = 25), "Quantitative Science Unit - Gentle's lab"),
    #             right_text = actionLink("use_link", "Terms of Use")
    #          )
) 

server <- function(input, output, session) {
    
    # ----- Buttons for `Home` panel -----
    observeEvent(input$button1, {
        updateTabItems(session, "tabs", "data")
    })
    observeEvent(input$explore, {
        updateTabItems(session, "tabs", "metaz")
    })
    observeEvent(input$button2, {
        updateTabItems(session, "tabs", "metaz")
    })
    observeEvent(input$button3, {
        updateTabItems(session, "tabs", "individual_analysis")
        updateTabsetPanel(session, "ind_analysis",
                          selected = "panel1"
        )
    })
    observeEvent(input$button4, {
        updateTabItems(session, "tabs", "individual_analysis")
        updateTabsetPanel(session, "ind_analysis",
                          selected = "panel2"
        )
    })
    

    # observeEvent(input$Ibutton4, {
    #     updateTabItems(session, "tabs", "individual_analysis")
    #     updateTabsetPanel(session, "ind_analysis",
    #                       selected = "panel2"
    #     )
    # })
    
    #######################
    ####### OUTPUTS #######
    #######################
    
    # ----- Header (see code/dashboardHeader.R) -----
    output$login.user <- renderUI(HTML(session$userData$auth0_info$nickname))
    observeEvent(input$precog_home_page, {
        updateTabItems(session, "tabs", "home")
    })
    observeEvent(input$iprecog_home_page, {
        updateTabItems(session, "tabs", "iPRECOG_home")
    })
    output$cibersort_link <- renderUI(HTML('<a id = cybersort_link_header href = https://cibersort.stanford.edu> CIBERSORT </a>'))

    
    
    # ----- panel: Home (see code/home.R) -----
    output$home_text <- renderText({home_text1})
    
    # ----- panel: Home (iPRECOG) (see code/home.R) -----
    output$ihome_text <- renderText({ihome_text1})
    # Buttons
    observeEvent(input$ibutton1, {
        updateTabItems(session, "tabs", "signature")
    })
    observeEvent(input$iexplore, {
        updateTabItems(session, "tabs", "iMetaz")
    })
    observeEvent(input$ibutton2, {
        updateTabItems(session, "tabs", "iMetaz")
    })
    observeEvent(input$ibutton3, {
        updateTabItems(session, "tabs", "immune-fraction")
        # updateTabsetPanel(session, "ind_analysis",
        #                   selected = "panel1"
        # )
    })
    
    # ----- panel: Data. See code/dataUI.R -----
    output$downloadAnnotations <- annotations_download #FIXME
    output$text_data <- renderUI(text_data)
    output$precog_data <- data_output
    
    # FIXME
    # ----- Analysis - see `code/analysisUI.R` ----- 
    output$analysis1 <- analysis1
    output$analysis2 <- analysis2
    output$zscore.conversion <- zscore.output
    output$metaz <- metazOutput
    # output$metaz <- renderDataTable(
    #     datatable(
    #         metaz,
    #         style = 'default', 
    #         selection = list(mode="single", target="cell", selected = NULL),
    #         escape = F,
    #         rownames = F,
    #         filter = list(positiomn = "top", clear = T, plain = T),
    #         options = list(scrollX = T,
    #                        scrollY = T, 
    #                        headerCallback = JS(headerCallback), 
    #                        search = list(regex = TRUE), 
    #                        pageLength = input$numrows1, 
    #                        autoWidth = F, 
    #                        columnDefs = list(list(className = 'dt-center', width = '30px', targets = "_all")
    #                        )
    #         ),
    #     ) %>% formatStyle(names(metaz[,c(2:p)]), backgroundColor = styleInterval(brks, clrs))
    #     %>% formatRound(columns = colnames(metaz)[-c(1)], digits=2)
    # )
    
    
    # ----- Signature matrix (iPRECOG) -----
    output$ianalysis1 <- ianalysis_text1
    output$signature_matrix <- iprecog_sig_matrix_Output
    # ----- Metaz (iPRECOG) -----
    output$iMetaz1 <- iprecog_metaz_text
    output$izscore.conversion <- zscore.output
    output$iMetaz_table <- iprecog_metaz_table_Output
    # ----- Immune fractions (iPRECOG) ------
    output$immune_fraction_text <- immune_fraction_text_output
    output$immune_fraction_table <- iprecog_immune_fraction_table_Output
    
    # ----- Download tab -----
    output$dwld1 <- dwld1_button
    output$dwld2 <- dwld2_button
    output$dwld3 <- dwld3_button
    output$dwld4 <- dwld4_button
    output$dwld5 <- dwld5_button
    
    ######################################################
    ####### Table Cell selection in ANALYSIS panel #######
    ######################################################
    
    # ----- Modal for metaZscore table -----
    # Output 0
    htmlOutput("modal_intro")
    
    # AnyCellSelected <- reactive({
    #     RowSelected <- input$metaz_cells_selected[0,0]
    #     ColSelected <- input$metaz_cells_selected[1,1]
    # })
    
    selection <- reactiveValues()
    observe({
        
        if (length(input$metaz_cells_selected) != 0){
            selection$x <- input$metaz_cells_selected[1,1]
            selection$y <- input$metaz_cells_selected[1,2]
        }
        if (is.null(input$metaz_cells_selected)){
            selection$x <- -1
            selection$y <- -1
        }
        

    })

    # modal function
    myModal <- function(generate = F){
        
        x <- isolate(selection$x)
        y <- isolate(selection$y)
        gene_url = as.character(metaz$`Gene Symbol`[x]) # gene name with url 
        gene <- xpathSApply(htmlParse(gene_url), "//text()", xmlValue)
        cancer = as.character(cancer_list[y+1])
        num_accessions <- unique(precog_info$Accession_num[precog_info$Cancer == cancer])
        # list_of_accessions = str_replace_all(as.vector(precog_info$Accession[precog_info$Cancer == cancer]),"-",".")
        list_of_accessions = as.vector(precog_info$Accession[precog_info$Cancer == cancer])
        list_of_platforms = str_replace_all(as.vector(original_precog$Platform[original_precog$Cancer == cancer]),"-",".")
        selection$gene <- gene
        selection$cancer <- cancer
        # FIXME: get rid of it 
        km_info <- read.csv("precogdata/km_info.csv", header = T)[,-1] # FIXME: change location 
        paths = list() # list of file name for each accession
        # Start distinction of different cases
        if (y == 1){modal_text_1 <- "Unweighted meta Z scores"}
        # else if (y == 1){modal_text_1 <- "Gene Name"}
        else if (y == 0){modal_text_1 <- "Gene Symbol"}
        else {
            i = 1
            # Loop over all the available accessions for a given cancer :
            for (accession in precog_info$Accession[which(precog_info$Cancer == cancer)]){
                accession = as.character(accession)
                # accession = str_replace_all(accession, '-', '_')
                # FIXME
                cancer_first_row = precog_info$firstRow[precog_info$Cancer == cancer][1]
                for (j in c(1:dim(precog)[1])){
                    if(grepl(accession, precog$Accession[j])){
                        row = j ;}}
                # Retrieve the Accession URL
                accession_url = precog$Accession[row]
                # accession = str_replace_all(accession,"-",".")
                author = precog$First.Author[row] # FIXME: not working so good
                platform <- sub(".*\\::", "", precog_info$acc_and_plat[cancer_first_row + i -1])
                outcome <- as.character(original_precog$Outcome[which(original_precog$Accession == accession & original_precog$Platform == platform)])
                # Define `col` as study ID (used a lot of times below)
                col = tolower(str_replace_all(paste(cancer, accession, str_replace_all(platform, '/', '_'), sep ="."), " ", "."))
                # Import clinical annotation data
                cat('path to annot file: ', paste('precogdata/filtered.infofiles/', col, '.tsv', sep =''), '\n')
                annotation <- read.table(paste('precogdata/filtered_clinical_annotations/', col, '.tsv', sep =''), header = T)
                n_patients = 0
                n_patients = length(annotation[,1])
                
                # Outputs per accession: 
                # ----- Output 1 - information about the study -----
                id_str = paste(str_replace_all(accession, '-', '_') ,"_", as.character(i), sep="")
                cat('ID = ', id_str, '\n')
                htmlOutput(id_str)
                accession_text = paste("<b>Accession</b>: ", accession, " ",sep="")
                author_text = paste("<b>Author</b>: ", author, " ",sep="")
                patient_text = paste(" <b>Number of patients (OS/DSS)</b>: ", as.character(n_patients), "<br/>", "<br/>" ,sep = "")
                ### Store content in an auxiliary variable 
                assign(paste(id_str, "_text", sep=""), 
                       paste(accession_text, author_text , patient_text, sep="<br/>"))
                
                # ----- Output 2 - message for studies with less than 30 patients -----
                # cat('debug: ', paste(accession, "_2_", as.character(i),sep=""), '\n')
                id_str_2 = paste(str_replace_all(accession, '-', '_'), "_2_",as.character(i),sep="")
                htmlOutput(id_str_2)
                if (n_patients < 30){
                    output_2_text =  paste("<div class=\"alert alert-error\"> <strong>",
                                           " KM plot not shown for studies with less than 30 patients - 
                                           minimum needed for median split survival analysis.", "</strong></div>",
                                            sep="")
                    assign(paste(id_str_2, "_text", sep=""),  output_2_text)
                }
                else {
                    output_3_text = ""
                    assign(paste(id_str_2, "_text", sep=""),  output_3_text)
                }
                # ----- Output 3 - message for KM plots status -----
                id_str_3 = paste(str_replace_all(accession, '-', '_'),"_3_",as.character(i),sep="")
                button_id <- paste(accession,"_but_",as.character(i),sep="")
                # Where we should look for the plots:
                folder <- paste(col, tolower(outcome), 'KMplots/',sep = '.')
                file_name_list = system(paste('ls precogdata/KMplots/', folder, sep = ''), intern = T)[which(grepl(tolower(gene), system(paste('ls precogdata/KMplots/', folder, sep = ''), intern = T)))]
                # We might have several plots for this gene 
                path_list <- paste("precogdata/KMplots/", folder, file_name_list, sep ="")
                paths[[col]] = path_list # we keep track of it

                if (((!col %in% tolower(colnames(km_info))) && n_patients >= 30 ) || ((km_info[[col]][km_info$gene==gene]==-1) && n_patients >= 30)) {
                    output_3_text <- paste("<div class=\"alert alert-warning\"> <strong>",
                                           "Plots are not available for this cancer/study/gene",
                                           "</strong></div>",
                                           # "<button id=", col, 
                                           # " type=\"button\" class=\"btn btn-default action-button\">Generate KM plot</button>",
                                           "<hr>", sep="")
                    assign(paste(id_str_3, "_text", sep=""),  output_3_text)
                    plot_id = paste(str_replace_all(accession, '-', '_'),"_plot_",as.character(i),sep="")
                    # assign(paste(plot_id, "_path_1", sep=""),  "")
                    for (k in 1:length(path_list)){
                        assign(paste(plot_id, "_path_", k ,sep=""),  "")
                    }
                }
                else if ((km_info[[col]][km_info$gene==gene]==1) && n_patients >= 30){
                    # In this case, we want to insert KM plot
                    output_3_text <- ''
                    assign(paste(id_str_3, "_text", sep=""),  output_3_text)
                    plot_id = paste(accession,"_plot_",as.character(i),sep="") 
                    uiOutput(plot_id)
                    cancer_name <- str_replace(cancer, " c", "_C")
                    
                    # ----- Case 1: display KM plots at gene level only (speed up modal opening)
                    # file_name = paste(gene, col, tolower(outcome), 'png', sep = '.')
                    # assign(paste(id_plot, "_path", sep=""),  str_replace('<img src="logo.png", height="500px, width = 500px", alt = "KM plot", style="vertical-align:middle;margin:0px 10px", loading="lazy"/>', "logo.png", path_list[1]))
                    # ----- Case 2: display KM plots at probe level
                    for (k in 1:length(path_list)){
                        assign(paste(plot_id, "_path_", k ,sep=""),  str_replace('<img src="logo.png" height="500px width = 500px" alt = "KM plot", style="vertical-align:middle;margin:0px 10px" loading="lazy"/>', "logo.png", 
                                                                                 path_list[k])
                               ) 
                    }
                }
                else if ((km_info[[col]][km_info$gene==gene]==0) && n_patients >= 30){
                   
                    output_3_text <- paste("<div class=\"alert alert-error\"> <strong>",
                                           "KM plot not available for this gene and study.",
                                           "</strong></div>","<hr>", sep="")
                    assign(paste(id_str_3, "_text", sep=""),  output_3_text)
                    
                    plot_id = paste(accession,"_plot_",as.character(i),sep="")
                    assign(paste(plot_id, "_path_1", sep=""),  "")
                }
                else{
                    assign(paste(id_str_3, "_text", sep=""),  "")
                    plot_id = paste(accession,"_plot_",as.character(i),sep="")
                    assign(paste(plot_id, "_path_1", sep=""),  "")
                }
                # Iteration 
                i = i+1
            }
        }
          
        modalDialog(
            
         #    if (generate)
         #        # Ask confirmation for generating a specific plot
         #        div(HTML("Are you sure to generate KM plot for <b>" , gene, 
         #                   "</b> in <b>", cancer,   "</b> for <b>", isolate(selection$accession) ,
         #                 "</b> (platform: ", isolate(selection$platform) , ")? <br>")),
         #    if (generate)
         #        fluidRow(
         #            column(6, align="center", offset = 3,
         #                   actionButton("generate", label = "OK"),
         # tags$style(type='text/css', "#button { vertical-align- middle; height- 50px; width- 100%; font-size- 30px;}")
         #                   )),
                
            if (y==1 && !generate)
            title = "KM Plot(s)",
            if (y>1 && !generate)
            title = div(HTML(paste("KM Plot(s) for gene <b>", gene, "</b> in <i>", cancer, "</i>","<hr>", sep=" "))),
            ## Output 0 (for all modal windows)
            if (y>1 && !generate)
            output$modal_intro <- renderUI(HTML(
                paste("<div class=\"alert alert-info\"> <strong>", as.character(num_accessions), 
                      " dataset(s) represented below", "</strong></div>", sep="")
            )),
            
            # Note: max 2 plots per accession
            if (y>1 && !generate)
            lapply(1:num_accessions, function(i) {
                
                acc = str_replace_all(list_of_accessions[i], '-', '_')
                # acc = list_of_accessions[i]
                col = tolower(str_replace_all(paste(cancer, acc, list_of_platforms[i], sep ="."), " ", "."))
                
                outputId <- paste(acc,"_",as.character(i),sep="")
                cat('OUTPUT_ID =', outputId, '\n')
                outputId_2 <- paste(acc,"_2_",as.character(i),sep="")
                outputId_3 <- paste(acc,"_3_",as.character(i),sep="")
                button_id <- paste(acc,"_but_",as.character(i),sep="")
                outputId_plot <- paste(acc,"_plot_",as.character(i),sep="")
                num_plots = length(paths[[col]])
                # cat('col: ', col, '\n')
                # cat(list_of_accessions[i], 'num of plots to display: ', num_plots, '\n')
                
                if (num_plots <= 1){
                    output[[outputId]] <- renderUI(tagList(
                        # Output 1
                        HTML(eval(parse(text = paste(outputId, "_text", sep="")))),
                        # Output 2
                        HTML(eval(parse(text = paste(outputId_2, "_text", sep="")))),
                        # Output 3
                        HTML(eval(parse(text = paste(outputId_3, "_text", sep="")))),
                        # Output 3 - plot
                        HTML(eval(parse(text = paste(outputId_plot, "_path_", 1, sep=""))),
                             '<p style="color:black"></p>', "<hr>"),
                        # HTML(eval(parse(text = paste(outputId_plot, "_path_", 2, sep=""))),
                        #      '<p style="color:black"></p>', "<hr>")
                        
                    ))
                } else{ #TODO: add more than 2 plots (if applicable) per accession
                    output[[outputId]] <- renderUI(tagList(
                        # Output 1
                        HTML(eval(parse(text = paste(outputId, "_text", sep="")))),
                        # Output 2
                        HTML(eval(parse(text = paste(outputId_2, "_text", sep="")))),
                        # Output 3
                        HTML(eval(parse(text = paste(outputId_3, "_text", sep="")))),
                        # Output 3 - plots
                        HTML(eval(parse(text = paste(outputId_plot, "_path_", 1, sep=""))),
                             '<p style="color:black"></p>', "<hr>"),
                        HTML(eval(parse(text = paste(outputId_plot, "_path_", 2, sep=""))),
                             '<p style="color:black"></p>', "<hr>")

                    ))}
            }),
                
            easyClose = TRUE, 
            footer = tagList(modalButton("Close"))
        )
    }
    
    # DISPLAY KM plots
    observeEvent(
        ignoreInit = T,         
        ignoreNULL = T, 
        suspended = F, 
        once = F,
        input$metaz_cells_selected,
        {
            suppressWarnings({
            if (dim(input$metaz_cells_selected) != 0){
                showModal(myModal())
                }
            })
    }
    )
    
    # # Handle generate KM plot button
    # observe({
    #     # Identify all of the buttons in the table.
    #     lapply(study.id.list,
    #            function(x){
    #                observeEvent(
    #                    input[[x]],
    #                    {
    #                      x.split <- str_split(x, '\\.')[[1]]
    #                      accession <- accession.main.list [match(x.split, tolower(accession.main.list))]
    #                      accession <- accession[!is.na(accession)]
    #                      platform <- platform.main.list[match(x.split, tolower(platform.main.list))]
    #                      platform <- platform[!is.na(platform)]
    #                      selection$accession <- accession
    #                      selection$platform <- platform
    #                      showModal(myModal(generate = TRUE))
    #                    }
    #                )
    #            })
    # })
    # observeEvent(input$generate, {
    #     cancer <- isolate(selection$cancer)
    #     accession <- isolate(selection$accession)
    #     platform <- isolate(selection$platform)
    #     gene <- isolate(selection$gene)
    #     withProgress(message = 'Making KM plot', value = 0, {
    #         accession <- isolate(selection$accession)
    #         platform <- isolate(selection$platform)
    #         gene <- isolate(selection$gene)
    #         incProgress(1/4, detail = paste('Sourcing all data...'))
    #         library('preprocessCore')
    #         data <- importData(cancer, accession, platform)
    #         incProgress(1/4, detail = 'Succesfully imported!')
    #         Sys.sleep(1)
    #         incProgress(1/4, detail = 'Computing KM plot')
    #         library('survival')
    #         library('survminer')
    #         generateKMPlotNEW(data, TRUE, 1, 23287, spe.gene = gene) 
    #         incProgress(1/4, detail = 'Done!')
    #         Sys.sleep(1)
    #     })
    #     removeModal()
    # })
    
    ###################################
    ####### INDIVIDUAL ANALYSIS #######
    ###################################
    # ----- see code/individual_analysisUI.R -----
    output$precog.ind.analysis1 <- analysis1.text.precog
    output$precog.ind.analysis2 <- analysis2.text.precog
    output$zscore.conversion.precog <- zscore.output.precog 
    output$precog.zscore <- precog.zscore.Output
    output$TCGAanalysis1 <- analysis1.text.TCGA
    output$TCGAanalysis2 <- analysis2.text.TCGA
    output$TCGA.metaz <- TCGA.metazOutput
    output$zscore.conversion.tcga <- zscore.output.tcga
    
    # ----- Modal for TCGA zscore table -----
    
    # Outputs for ALL KM modal window 
    # Output 1
    # htmlOutput("TCGA.modal.intro")
    
    selection <- reactiveValues()
    observe({
        
        if (length(input$TCGA.metaz_cells_selected) != 0){
            selection$x.tcga <- input$TCGA.metaz_cells_selected[1,1]
            selection$y.tcga <- input$TCGA.metaz_cells_selected[1,2]
        }
        if (is.null(input$TCGA.metaz_cells_selected)){
            selection$x.tcga <- -1
            selection$y.tcga <- -1
        }
        
    })
    
    # TCGA MODAL FOR KM PLOTS 
    myTCGAModal <- function(outcome = 'OS'){
        
        x <- isolate(selection$x.tcga)
        y <- isolate(selection$y.tcga)
        gene <- TCGA.metaz$Name[x]
        cancer.id = as.character(TCGA.cancer.id.list[y])
        cancer = cancer.names[cancer.id]


        if (y == 0){modal_text_1 <- 'Gene Name'}
        else {

                ##############################################################################
                # ----- insert plot -----
                # id_plot = paste(gene, '.', cancer,"_plot", sep="")
                    # uiOutput(id_plot)
                uiOutput('KM.plot')
                if (outcome == 'OS'){
                    link <- paste("precogdata/TCGA.KMplots/", cancer.id, '_OS/', gene, '.png',  sep ="") 
                } else if (outcome == 'DSS'){
                    link <- paste("precogdata/TCGA.KMplots/", cancer.id, '_DSS/', gene, '.png',  sep ="")
                }
                
                # id_plot = paste(cancer.id, gene,"_plot_", sep="")
                    # 
                    # ----- PARAMETERS FOR THE KM_PLOTS ------
                text <- str_replace('<img src="logo.png", height="500px, width = 500px"
                            style="vertical-align:middle;margin:0px 10px", loading="lazy"/>', "logo.png", link)
                    # 
                    # # cat("NEW TEXT: ", text, '\n')
                    # 
                # assign(paste(id_plot, "_path", sep=""),  text)
                path2plot = text
                
                # DSS.button <- "<button id= DSS type=\"button\" class=\"btn btn-default action-button\"> Display DSS </button>"
        }
        
        modalDialog(
            if (y == 0)
                title = modal_text_1,
            if (y ==1) 
                title = 'PRECOG meta Z-score',
            if (y ==2)
                title = 'TCGA metaz Z-score',
            if (y>2 & outcome == 'OS')
                title = div(HTML(paste("KM Plot for gene <b>", gene, "</b> in <i>", cancer, "</i> (OS outcome)","<hr>", sep=" "))),
            if (y>2 & outcome == 'DSS')
                title = div(HTML(paste("KM Plot for gene <b>", gene, "</b> in <i>", cancer, "</i> (DSS outcome)","<hr>", sep=" "))),
            if (y > 2) 
                output$KM.plot <- renderUI(HTML(path2plot, '<p style="color:black"></p>', "<hr>")),
            if (y >2 & outcome == 'OS')
                fluidRow(
                    column(6, align="center", offset = 3,
                           actionButton('DSS_button', 'Display DSS outcome'),
                           tags$style(type='text/css', "#button { vertical-align- middle; height- 50px; width- 100%; font-size- 30px;}")
                    )),
            if (y >2 & outcome == 'DSS')
                fluidRow(
                    column(6, align="center", offset = 3,
                           actionButton('OS_button', 'Display OS outcome'),
                           tags$style(type='text/css', "#button { vertical-align- middle; height- 50px; width- 100%; font-size- 30px;}")
                    )),

            easyClose = TRUE, 
            footer = tagList(modalButton("Close"))
        )
        
    }
    
    
    # DISPLAY TCGA KM plots 
    observeEvent(
        ignoreInit = T,         
        ignoreNULL = T, 
        suspended = F, 
        once = F,
        input$TCGA.metaz_cells_selected,
        {
            suppressWarnings({
            if (dim(input$TCGA.metaz_cells_selected) != 0){
                showModal(myTCGAModal())  
            }
            })
        }
    )
    
    observeEvent(input$DSS_button, {
        showModal(myTCGAModal(outcome = 'DSS'))
    })
    
    observeEvent(input$OS_button, {
        showModal(myTCGAModal(outcome = 'OS'))
    })
    
    # ----- Modal for PRECOG individual analysis table -----
    selection <- reactiveValues()
    observe({
        
        if (length(input$precog.zscore_cells_selected) != 0){
            selection$x <- input$precog.zscore_cells_selected[1,1]
            selection$y <- input$precog.zscore_cells_selected[1,2]
        }
        if (is.null(input$precog.zscore_cells_selected)){
            selection$x <- -1
            selection$y <- -1
        }
        
    })
    
    # PRECOG MODAL FOR KM PLOTS 
    myPRECOGModal <- function(failed=T){
        
        x <- isolate(selection$x)
        y <- isolate(selection$y)
        # gene <- as.character(precog.indiv$Gene[x])
        # cancer = as.character(PRECOG.cancer.list[max(0, y-1)]) #FIXME
        # # Retrieve all cancer information 
        # id = sub(".*cancer_", "", cancer)
        # cancer <- gsub('[[:digit:]]+', "", cancer)
        # #   x = gsub('[[:digit:]]+', '', colnames(ind)[col])
        # cat('CANCER NAME: ', cancer, '\n')
        # accession <- as.character(original_precog$Accession[which(original_precog$PRECOG.ID == id)])
        # platform <- as.character(original_precog$Platform[which(original_precog$PRECOG.ID == id)])
        # outcome <- as.character(original_precog$Outcome[which(original_precog$Accession == accession & original_precog$Platform == platform)])
        # km_info <- read.csv("precogdata/km_info.csv", header = T)[,-1]

        if (y == 0){modal_text_1 <- 'Gene Symbol'}
        else if (y == 1){modal_text_1 <- 'Mean Z-Score'}
        else {
            gene <- as.character(precog.indiv$Gene[x])
            cancer = as.character(PRECOG.cancer.list[max(0, y-1)]) #FIXME
            # Retrieve all cancer information 
            id = strsplit(cancer, '_')[[1]][length(strsplit(cancer, '_')[[1]])]
            # cat('id = ', id, '\n')
            cancer <- gsub('[[:digit:]]+', "", cancer)
            cancer <- substr(cancer, 1, nchar(cancer)-1)
            #   x = gsub('[[:digit:]]+', '', colnames(ind)[col])
            # cat('CANCER NAME: ', cancer, '\n')
            accession <- as.character(original_precog$Accession[which(original_precog$PRECOG.ID == id)])
            # cat('ACCESSION:', accession, '\n')
            platform <- as.character(original_precog$Platform[which(original_precog$PRECOG.ID == id)])
            outcome <- as.character(original_precog$Outcome[which(original_precog$Accession == accession & original_precog$Platform == platform)])
            km_info <- read.csv("precogdata/km_info.csv", header = T)[,-1]
            
            cancer <- str_replace_all(cancer, '_', '.')
            col = tolower(str_replace_all(paste(cancer, accession, platform, sep ="."), " ", "."))
            col <- str_replace_all(col, '\\.\\.', '\\.')
            # Message output
            htmlOutput('info_message')
            isAvailable = TRUE
            if (!col %in% tolower(colnames(km_info)) || km_info[[col]][km_info$gene==gene]==-1){
                info_message <- paste("<div class=\"alert alert-warning\"> <strong>",
                                      "Plots are not available for this cancer/study/gene",
                                      "</strong></div>",
                                      "<hr>", sep="")
                isAvailable <- FALSE
            }
            # KM Plot(s) output
            
            folder <- paste(col, tolower(outcome), 'KMplots/',sep = '.')
            pathToFile = paste('ls precogdata/KMplots/', folder, sep = '')
            # cat('PATH TO FILE: ', pathToFile, '\n')
            file_name_list = system(pathToFile, intern = T)[which(grepl(tolower(gene), system(paste('ls precogdata/KMplots/', folder, sep = ''), intern = T)))]
            # We might have several plots for this gene 
            path_list <- paste("precogdata/KMplots/", folder, file_name_list, sep ="")
            # cat('number of plots for this gene/study: ', length(path_list), '\n')
            num_plots <- length(path_list)
            uiOutput('KM_plot_1')
            path2plot_1 <- str_replace('<img src="logo.png" height="500px width = 500px" alt = "KM plot", style="vertical-align:middle;margin:0px 10px" loading="lazy"/>', "logo.png", 
                                       path_list[1])
            uiOutput('KM_plot_2')
            path2plot_2 <- str_replace('<img src="logo.png" height="500px width = 500px" alt = "KM plot", style="vertical-align:middle;margin:0px 10px" loading="lazy"/>', "logo.png", 
                                       path_list[2])
            uiOutput('KM_plot_3')
            path3plot_3 <- str_replace('<img src="logo.png" height="500px width = 500px" alt = "KM plot", style="vertical-align:middle;margin:0px 10px" loading="lazy"/>', "logo.png", 
                                       path_list[3])
            # we can add more paths if we want more plots to be displayed...
               }
        
        modalDialog(
            if (y == 0)
                title = modal_text_1,
            if (y == 1)
                title = modal_text_1, 
            if (y > 1)
                title = div(HTML(paste("KM Plot(s) for gene <b>", gene, 
                                       "</b> in <i>", str_replace_all(cancer, '.', ' '), 
                                       "</i>",' (', accession ,')',"<hr>", sep=" "))),
            if (y > 1 && !isAvailable)
                output$info_message <- renderUI(HTML(info_message)),
            if (y > 1 && isAvailable)
                output$KM_plot_1 <- renderUI(HTML(path2plot_1, '<p style="color:black"></p>', "<hr>")),
            if (y > 1 && num_plots > 1 && isAvailable)
                output$KM_plot_2 <- renderUI(HTML(path2plot_2, '<p style="color:black"></p>', "<hr>")),
            if (y > 1 && num_plots > 2 && isAvailable)
                output$KM_plot_3 <- renderUI(HTML(path3plot_3, '<p style="color:black"></p>', "<hr>")),
            
            easyClose = TRUE, 
            footer = tagList(modalButton("Close"))
        )
        
    }
    
    # DISPLAY PRECOG KM plots 
    observeEvent(
        ignoreInit = T,         
        ignoreNULL = T, 
        suspended = F, 
        once = F,
        input$precog.zscore_cells_selected,
        {
            suppressWarnings({
            if (dim(input$precog.zscore_cells_selected) != 0){showModal(myPRECOGModal())}
            })
        }
    )
    
    
    
    
    ####################### 
    ####### CONTACT ####### 
    ####################### 
    
    output$andrew <- andrew.output
    output$aaron <- aaron.output
    output$chih <- chih.output
    output$sylvia <- sylvia.output
    output$ash <- ash.output
    
    ############################
    ####### TERMS OF USE #######
    ############################
    
    observeEvent(input$use_link,
            {
                terms_of_use_text = "The Board of Trustees of the Leland Stanford Junior University 
                (“Stanford”) provides PRECOG website features and 
                services (“Service”) free of charge for non-commercial use only.
                Use of the Service by any commercial entity for any purpose, including research, is prohibited."
                output$terms_of_use <- renderText(terms_of_use_text)
            }
               )
}

shinyApp(ui, server)
#auth0::shinyAppAuth0(ui, server)
