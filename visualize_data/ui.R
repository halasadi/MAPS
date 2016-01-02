library(shiny)

# Define UI for application that draws a histogram
shinyUI(navbarPage(title = 'Visualization of IBD segments',
                   theme = 'flatly.css',
                   inverse = FALSE,
                   
                   tabPanel('Home',
                            
                            sidebarLayout(
                              sidebarPanel(
                                
                                radioButtons("datatype", "Dataset:",
                                             c("Western European" = 'popres',
                                               "Eastern European" = 'ee',
                                               "POBI" = 'pobi'), selected = 'popres'),
                                
                                sliderInput("lower_thr",
                                            "Minimum threshold:",
                                            min = 2,
                                            max = 12,
                                            value = 6),
                                
                                radioButtons("isFinite", label = "Is the upper bound finite?", 
                                             choices = list("Finite" = 'finite', "Infinite" = 'infinite'),
                                             selected = 'infinite'),
                                
                                conditionalPanel(
                                  condition = "input.isFinite == 'finite'",
                                  sliderInput("upper_thr",
                                              "Maximum threshold (must be greater than min. thr.):",
                                              min = 3,
                                              max = 12,
                                              value = 4)
                                )
                                
                              ),
                              
                              # Show a plot of the generated distribution
                              mainPanel(
                                
                                tabsetPanel(
                                  tabPanel('PCA', plotOutput("pca_plot", width = '600px', height = '600px')),
                                  tabPanel('Log Similarity Matrix', plotOutput("sim_plot", width = '600px', height = '600px')),
                                  tabPanel('Similarity Graph', plotOutput("graph_plot", width = '600px', height = '600px')),
                                  tabPanel('MAPS - Migration', plotOutput("maps_plot_migration", width = '600px', height = '600px')),
                                  tabPanel('MAPS - Population Size', plotOutput("maps_plot_population", width = '600px', height = '600px'))
                                )
                              )
                            )
                            
                   ),
                   
                   tabPanel('About',
                            fluidRow(column(width = 10, offset = 1,
                                            includeMarkdown('include/about.md')
                                            ))
                            )
                   
))
