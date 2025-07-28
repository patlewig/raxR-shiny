if (!requireNamespace("raxR", quietly=TRUE)){
install.packages("devtools")
devtools::install_github("patlewig/raxR")
}


library(shiny)
library(ggplot2)
library(plotly)
library(rcdk)
library(tidyverse)
library(bslib)
library(raxR)


data_path <- system.file("extdata", "small_acute_processed.csv", package = "raxR")
data_path_1 <- system.file("extdata", "smi_acute.csv", package = "raxR")
toxicity_data <- read.csv(data_path)
toxicity_data <- toxicity_data %>% select(dtxsid,  LD50_LM)
source_analogues <- process_substances(data_path_1)
source_FP <- generate_fingerprints(source_analogues)
source_chem <- generate_physchem(generate_mol(source_analogues))

ui <- page_fluid(
  navset_tab(
  nav_panel("Dataset",
              h4('Overview of the acute toxicity potency values distribution (expressed as LD50_LM) and the fulldataset'),
  plotlyOutput(outputId = "hist"),
  dataTableOutput("dynamic")), 
nav_panel("Analogue ID", 
          h5('Finding analogues on the basis of circular chemical fingerprints'),
          textInput(inputId = "smiles", "Get analogues for a Target chemical", placeholder = "Enter SMILES..."),
          textInput(inputId = "dtxsid", "", placeholder = "Enter identifier.."),
          numericInput(inputId = "k", "no of analogues", value = 10, min =1, max = 15),
          dataTableOutput('static')),
nav_panel("Physicochemical Profile",
          h5('Exploring whether any source analogues exceed a 3SD threshold with respect to their physicochemical profile'),
          dataTableOutput('static_2'),
          plotlyOutput(outputId = 'box')),
nav_panel("Read-across prediction",
          h5('Making a similarity weighted activity calculation to predict the toxicity of a target chemical based on its source analogues'), 
          dataTableOutput('static_3')
          )
)
)

server <- function(input, output, session) {
  # Panel 1 - histogram and plot
  output$hist <- renderPlotly({toxicity_data %>% ggplot(aes(LD50_LM)) + geom_histogram()})
  output$dynamic <- renderDataTable(toxicity_data, options = list(pageLength = 5))
  
  # Panel 2 - Get Analogues
  
  # Reactive: generate target fingerprint table once
  target_fp_tbl <- reactive({
    req(input$smiles, input$dtxsid)
    generate_fp(input$smiles, input$dtxsid)
  })
  
  # Reactive: calculate analogues using the fingerprint
  analogues <- reactive({
    req(target_fp_tbl(), input$k)
    if (nrow(target_fp_tbl()) == 0 || is.null(target_fp_tbl()$fingerprint[[1]])) {
      return(NULL)
    }
  
  get_analogues(source_FP, target_fp_tbl()$fingerprint[[1]], input$k)
})
  
  output$static <- renderDataTable({
    req(target_fp_tbl())
    
    if (is.null(analogues())) {
      return(tibble(Error = "No valid analogues found."))
    }
    
    analogues()
  })
  
  static2_data <- reactive({req(input$smiles, input$dtxsid, analogues())
              
                test_chem <- generate_physchem(generate_mol_obj(input$smiles, input$dtxsid))
                all_pchem <- create_pchem(source_chem, analogues(), test_chem)                      
                long_pchem <- prep_df(all_pchem)                
                })
  
  

  
  output$static_2 <- renderDataTable({
    filter_sd(static2_data())
  })
  
  output$box <- renderPlotly({
    add_scale(static2_data()) %>% ggplot(aes(x = property, y = property_scaled)) + geom_boxplot(outlier.shape = NA) + geom_jitter()
  })

 df <- reactive({req(input$dtxsid, analogues())
                   with_data <- analogues() %>% left_join(toxicity_data, by = 'dtxsid') %>% drop_na("LD50_LM")
                     wtavg(id = input$dtxsid, with_data, outcome_col = 'LD50_LM')
 })
 
 output$static_3 <- renderDataTable({
   df()
 })

}
shinyApp(ui, server)
