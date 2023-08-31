# ----- PRECOG -----
# RScript for dashboard home page 
home_title <- "PRECOG - PREdiction of Clinical Outcomes from Genomic Profiles"
home_text1 <- "PRECOG is a system for querying associations between 
              genomic profiles and cancer outcomes. It enables researchers to query whether,
              for example, high expression of a gene is prognostic for shorter or longer patient survival."

learn_more <- "We assembled, curated, and integrated cancer gene expression and clinical outcome data from the 
              public domain into a new resource for PREdiction of Clinical Outcomes from Genomic profiles (or PRECOG). 
              PRECOG encompasses 166 cancer expression data sets, including overall survival data for ~18,000 patients 
              diagnosed with 39 distinct malignancies. For more information, please refer to our publication by Gentles/Newman 
              et al. in Nature Medicine (2015)."

des1 <- "Find further information on the data included in PRECOG and the annotation and analysis steps."
des2 <- "Explore associations between gene expression and overall survival collapsed by cancer type."
des3 <- "Examine z-score survival associations in the individual datasets, which are used to compute the meta-Z scores."

# UI 
home <-  tabItem(tabName = "home",
                 # h2("Home"),
                 fluidPage(
                   includeCSS('code/css/home.css'),
                   titlePanel(h2(home_title, id = "home_title")),
                   br(),
                   fluidRow(textOutput("home_text")),
                   br(),
                   fluidRow(
                     accordion(
                       inputId = 'acc1',
                       accordionItem(
                         
                         title = "Learn More",
                         color = NULL,
                         collapsed = TRUE,
                         # id = 'acc1',
                         tags$span(learn_more, id = 'learn_more'),
                         div(br()), 
                         div(
                           # actionButton("explore", "Explore PRECOG"),
                           # icon("wpexplorer"),
                           actionButton("explore", "Explore PRECOG"),
                           style="text-align: center;")
                         
                       )
                       
                     ),
                     # br()
                   ),
                   
                   
                   fluidRow(
                     column(div(h2('PRECOG Datasets') ,style="text-align: center;"),
                            width = 4,
                            des1,
                            tags$br(),
                            tags$br(), 
                            div(actionButton("button1", "View Details"), style="text-align: center;")
                        ),
                     
                     column(div(h2('Meta-Z Analysis') ,style="text-align: center;"),
                            width = 4,
                            des2,
                            tags$br(),
                            tags$br(), 
                            div(actionButton("button2", "View Details"), style="text-align: center;")
                     ), 
                     
                     column(div(h2('Individual Datasets Analysis'),style="text-align: center;"),
                            width = 4,
                            des3,
                            tags$br(),
                            tags$br(), 

                            div(actionButton(
                              style = "display:inline-block",
                              "button3", "View Details (PRECOG)"), style="float:right"),
                            div(
                              style = "display:inline-block",
                              actionButton("button4", "View Details (TCGA)"), style="float:right")
                     )
                     
                   )
                 )
)
# ----- iPRECOG -----
#  RScript for dashboard home page 
ihome_title <- "iPRECOG - immune PREdiction of Clinical Outcomes from Genomic Profiles"
ihome_text1 <- "iPRECOG provides an immune-oriented view and analysis of PRECOG. For more information, please refer to our publication by Gentles/Newman et al. in Nature Medicine (2015)."

ilearn_more <- "iPRECOG applies CIBERSORT to PRECOG datasets using immune cell specific gene expression profiles, deconvolving their contributions in each patient’s tumor sample. These inferred immune fractions can then be associated with patient outcome in different malignancies, similarly to how PRECOG relates raw gene expression to survival."

ides1 <- "The signature matrix specifies the immune cell type specific genes that underlie iPRECOG."
ides2 <- "Associations between inferred cell fractions and survival outcomes."
ides3 <- "Fractions of 22 leukocyte subsets estimated in PRECOG bulk tumor samples."

#  UI 
iPRECOG_home <- tabItem(tabName = "iPRECOG_home",
                        # h2("Home"),
                        fluidPage(
                          includeCSS('code/css/home.css'),
                          titlePanel(h2(ihome_title, id = "ihome_title")),
                          br(),
                          fluidRow(textOutput("ihome_text")),
                          br(),
                          fluidRow(
                            accordion(
                              inputId = 'Iacc1',
                              accordionItem(
                                
                                title = "Learn More",
                                color = NULL,
                                collapsed = TRUE,
                                # id = 'acc1',
                                tags$span(ilearn_more, id = 'ilearn_more'),
                                div(br()),
                                div(
                                  # actionButton("explore", "Explore PRECOG"),
                                  # icon("wpexplorer"),
                                  actionButton("Iexplore", "Explore iPRECOG"),
                                  style="text-align: center;")
                                
                              )
                              
                            ),
                            # br()
                          ),
                          
                          
                          fluidRow(
                            column(div(h2('iPRECOG signature matrix') ,style="text-align: center;"),
                                   width = 4,
                                   ides1,
                                   tags$br(),
                                   tags$br(),
                                   div(actionButton("ibutton1", "View Details"), style="text-align: center;")
                            ),
                            
                            column(div(h2('iPRECOG meta-Z matrix') ,style="text-align: center;"),
                                   width = 4,
                                   des2,
                                   tags$br(),
                                   tags$br(),
                                   div(actionButton("ibutton2", "View Details"), style="text-align: center;")
                            ),
                            
                            column(div(h2('iPRECOG cancer-specific immune fractions'),style="text-align: center;"),
                                   width = 4,
                                   des3,
                                   tags$br(),
                                   tags$br(),
                                   
                                   div(actionButton(
                                     style = "display:inline-block",
                                     "ibutton3", "View Details"), style="text-align: center")
                                   # div(
                                   #   style = "display:inline-block",
                                   #   actionButton("Ibutton4", "View Details (TCGA)"), style="float:right")
                            )
                            
                          )
                        )
)
