#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Hmisc)
library(janitor)
library(Matrix)
library(viridis)
library(factoextra)

cell_types <- c("Acinar", "Alpha", "Beta", "Cycling Alpha", "Delta", "Ductal", "Endothelial", "Gamma + Epsilon", "Immune (Macrophages)", "MUC5B+ Ductal", "Quiescent Stellate")

color_vars <- c("Description of diabetes status", "Age (years)", "Sex", "BMI", "Ethnicity", "Cause of death", "HbA1C percentage", "Chemistry", "Program")

#set up color palette
all_palette <- colorRampPalette(c("#FFBE0B", "#FB5607", "#FF006E", "#8338EC", "#3A86FF"))


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  shinybrowser::detect(),
  tags$head(
    tags$style(HTML("@import url('https://fonts.googleapis.com/css2?family=Open+Sans:ital,wght@0,300..800;1,300..800&display=swap');
                            .tabbable > .nav > li > a {color:black; border-top: 5px solid transparent; border-left: 0px solid transparent; border-right: 0px solid transparent; border-bottom: 5px solid transparent}
                            .tabbable > .nav > li[class=active] > a {color:rgba(33, 145, 151, 1); border-top: 5px solid transparent; border-left: 0px solid transparent; border-right: 0px solid transparent; border-bottom: 5px solid rgba(33, 145, 151, 1)}
                            .tabbable > .nav > li > a:hover {color:rgba(33, 145, 151, 1); background-color:transparent}
                            .well {background-color: rgba(33, 145, 151, 0.8); color:black}
                            input[type='checkbox'] {accent-color:rgba(148, 201, 94, 1); color:black}
                            .selectize-dropdown-content {background-color: rgba(148, 201, 94, 1); color:black}
                            .selectize-dropdown-content .active {background-color:rgba(148, 201, 94, 1) !important; color:white !important}
                            .selectize-dropdown-content .selected {background-color:rgba(148, 201, 94, 1); color:rgba(33, 145, 151, 0.8)}
                            * {font-family:'google', 'Open Sans', sans-serif; font-weight: 600}
                            a {color:rgba(33, 145, 151, 1)}
}
  "))),

    # Application title
    titlePanel("PanKbase PCA of pseudobulk RNA data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("CellType",
                        "Cell Type:",
                        choices = cell_types),
            selectInput("Color",
                        "Color by:",
                        choices = color_vars),
            downloadButton("DownloadPCA", 
                           "Download PCA")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("PCAPlot")
        )
    ),
    sidebarLayout(
      sidebarPanel(
        selectInput("PC",
                    "Principal Component",
                    choices = setNames(1:10,paste0("PC", 1:10))),
        downloadButton("DownloadContribs",
                       "Download Contributions")),
      mainPanel(
        p(),
        p(),
        plotOutput("contribsPlot")
      )
    )
)

# Define server logic required to draw PCA plots
server <- function(input, output) {
  
    PCA_fxn <- function() {
      # Get pseudobulk data for this cell type and filter for only untreated samples
      fname <- paste0("code/PCA_results/", input$CellType, "_PCA_results.csv")
      pca_res <- read.csv(fname)
      
      # Make values friendly, combine categories, etc.
      # Change capitalization and naming of some metadata fields
      pca_res$PanKbase_Description.of.diabetes.status <- recode(pca_res$PanKbase_Description.of.diabetes.status,
                                                                 "non-diabetic" = "No diabetes",
                                                                 "type 1 diabetes" = "Type 1 diabetes",
                                                                 "type 2 diabetes" = "Type 2 diabetes",
                                                                 "cystic fibrosis diabetes" = "Cystic fibrosis diabetes",
                                                                 "gestational diabetes" = "Gestational diabetes",
                                                                 "diabetes unspecified" = "Diabetes unspecified",
                                                                 "steroid-induced diabetes" = "Steroid-induced diabetes",
                                                                 "monogenic diabetes" = "Monogenic diabetes")

      pca_res$PanKbase_Ethnicities <- recode(na_if(pca_res$PanKbase_Ethnicities, ""),
                                              "African American,Black" = "African American, Black",
                                              "Caucasian" = "White",
                                              .missing = "Unknown")

      pca_res$PanKbase_Sex <- recode(pca_res$PanKbase_Sex,
                                      female = "Female",
                                      male = "Male")

      pca_res$source <- recode(pca_res$source,
                                "IIDP,Prodo" = "IIDP")

      # Rename values that have different capitalization so they get grouped together and change blanks to unknown
      pca_res$PanKbase_Cause.of.Death <- recode(na_if(pca_res$PanKbase_Cause.of.Death, ""),
                                                 "Cerebrovascular/stroke" = "Cerebrovascular/Stroke",
                                                 "Head Trauma" = "Head trauma",
                                                 "ICH/stroke" = "Cerebrovascular/Stroke",
                                                 "Cerebral Edema (DKA)" = "Cerebral edema (DKA)",
                                                 "IHC" = "Cerebrovascular/Stroke",
                                                 .missing = "Unknown")

      # Rename columns to friendly names
      pca_res <- pca_res %>% rename("Program" = source,
                                      "Age (years)" = PanKbase_Age..years.,
                                      "C. Peptide (ng/ml)" = PanKbase_C.Peptide..ng.ml.,
                                      "AAB-GADA value (unit/ml)" = PanKbase_AAB.GADA.value..unit.ml.,
                                      "AAB-IA2 value (unit/ml)" = PanKbase_AAB.IA2.value..unit.ml.,
                                      "AAB-IAA value (unit/ml)" = PanKbase_AAB.IAA.value..unit.ml.,
                                      "AAB-ZNT8 value (unit/ml)" = PanKbase_AAB.ZNT8.value..unit.ml.,
                                      "HbA1C percentage" = PanKbase_HbA1C..percentage.,
                                      "Hospital stay (hours)" = PanKbase_Hospital.Stay..hours.,
                                      "Description of diabetes status" = PanKbase_Description.of.diabetes.status,
                                      "Cause of death" = PanKbase_Cause.of.Death,
                                      "Sex" = PanKbase_Sex,
                                      "BMI" = PanKbase_BMI,
                                      "Ethnicity" = PanKbase_Ethnicities,
                                      "Chemistry" = chemistry)

      # Change NA AAB status to "Unknown"
      pca_res <- pca_res %>% mutate(across(starts_with("PanKbase_AAB_"), ~ifelse(is.na(.x), "Unknown", .x)))

      # Now change variable names to be friendly
      pca_res <- pca_res %>% rename("AAB-GADA Positive" = PanKbase_AAB_GADA_POSITIVE,
                                      "AAB-IA2 Positive" = PanKbase_AAB_IA2_POSITIVE,
                                      "AAB-IAA Positive" = PanKbase_AAB_IAA_POSITIVE,
                                      "AAB-ZNT8 Positive" = PanKbase_AAB_ZNT8_POSITIVE)
      

      if (is.numeric(pca_res[,input$Color]) == TRUE) {
        p <- ggplot(pca_res, aes(x = PC1, y = PC2, color = !!sym(input$Color))) +
          geom_point(size=2) +
          scale_color_gradientn(colors = viridis(length(unique(pca_res[,input$Color])))) +
          theme_minimal() +
          labs(color = input$Color) +
          theme(text = element_text(size=16),
                plot.title = element_text(hjust = 0.5)) +
          ggtitle(paste0(input$CellType, " Cells: CPM-normalized Counts"))
        print(p)
      }
      
      else {
        # Set up color palettes
        collections <- unique(pca_res[,input$Color]) %>% sort()
        collection_pal <- all_palette(length(collections))
        p <- ggplot(pca_res, aes(x = PC1, y = PC2, color = !!sym(input$Color))) + 
          geom_point(size=2) + 
          theme_minimal() +
          scale_color_manual(values = collection_pal) +
          labs(color = input$Color) +
          theme(text = element_text(size=14),
                plot.title = element_text(hjust = 0.5)) +
          ggtitle(paste0(input$CellType, " Cells: CPM-normalized Counts"))
        print(p)
      }

    }

    output$PCAPlot <- renderPlot({
        PCA_fxn()
    })
    
    contribs_fxn <- function() {
      fname <- paste0("code/PCA_results/", input$CellType, "_all_PCA_results.rds")
      pca_res_all <- readRDS(fname)

      p <- fviz_contrib(pca_res_all, choice = "var", axes = as.numeric(input$PC), top = 10, fill = "#219197", color = "#219197", ggtheme = theme_classic()) +
        ggtitle(paste0("Contributions of genes to PC", input$PC)) +
        theme(axis.text = element_text(size=12),
              axis.title = element_text(size=14),
              plot.title = element_text(size=16, hjust = 0.5))
      print(p)
    }
    
    output$contribsPlot <- renderPlot({
      contribs_fxn()
    })
    
    output$DownloadPCA <- downloadHandler(
      filename = function() {paste0("PanKbase_GeneExpr_PCA_", Sys.time(), ".png")},
      content = function(file) {
        ggsave(file, plot = PCA_fxn(), width = 12, height = 8, units = "in", bg = "white")
      }
    )
    
    output$DownloadContribs <- downloadHandler(
      filename = function() {paste0("PanKbase_GeneExpr_PCA_Contributions_", Sys.time(), ".png")},
      content = function(file) {
        ggsave(file, plot = contribs_fxn(), width = 12, height = 8, units = "in", bg = "white")
      }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
