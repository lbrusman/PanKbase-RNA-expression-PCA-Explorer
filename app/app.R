#
# This is the Shiny web application that runs the PCA Explorer for PanKbase
# See it live at https://pankbase.org/donor-metadata.html
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

## Set up some global variables ================================================

# Define Cell types
cell_types <- c("Acinar", "Active Stellate", "Alpha", "Beta", "Cycling Alpha", "Delta", "Ductal", "Endothelial", "Gamma + Epsilon", "Immune", "MUC5B+ Ductal", "Quiescent Stellate")

# Define variables to color by
color_vars <- c("Description of diabetes status", "Age (years)", "Sex", "BMI", "Ethnicity", "Cause of death", "HbA1C percentage", "Treatment", "Chemistry", "Program")

# Set up universal color palette
all_palette <- colorRampPalette(c("#FFBE0B", "#FB5607", "#FF006E", "#8338EC", "#3A86FF"))


## Define UI for application that plots PCAs ===================================
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


    # To plot PCA -------------------------------------------------------------- 
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
        mainPanel(
           p(),
           p(),
           plotOutput("PCAPlot")
        )
    ),
    # To plot PCA contributions ------------------------------------------------
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

## Define server logic required to draw PCA plots ==============================

server <- function(input, output) {
    
    # Make PCA plot ------------------------------------------------------------
    PCA_fxn <- function() {
      # Make sure name for "Immune (Macrophages)" is changed to "Immune"
      if (input$CellType == "Immune") {
        cell_name <- "Immune (Macrophages)"
      }
      else {
        cell_name <- input$CellType
      }
      
      # Read in pseudobulk data for this cell type
      fname <- paste0("outputs/PCA_results/", cell_name, "_PCA_results.csv")
      pca_res <- read.csv(fname)
      
      # Make values friendly, combine categories, change capitalization, etc.
      pca_res$pan_kbase_description_of_diabetes_status <- recode(pca_res$pan_kbase_description_of_diabetes_status,
                                                                 "non-diabetic" = "No diabetes",
                                                                 "type 1 diabetes" = "Type 1 diabetes",
                                                                 "type 2 diabetes" = "Type 2 diabetes",
                                                                 "cystic fibrosis diabetes" = "Cystic fibrosis diabetes",
                                                                 "gestational diabetes" = "Gestational diabetes",
                                                                 "diabetes unspecified" = "Diabetes unspecified",
                                                                 "steroid-induced diabetes" = "Steroid-induced diabetes",
                                                                 "monogenic diabetes" = "Monogenic diabetes")

      pca_res$pan_kbase_ethnicities <- recode(na_if(pca_res$pan_kbase_ethnicities, ""),
                                              "African American,Black" = "African American, Black",
                                              "Caucasian" = "White",
                                              .missing = "Unknown")

      pca_res$pan_kbase_sex <- recode(pca_res$pan_kbase_sex,
                                      female = "Female",
                                      male = "Male")

      pca_res$source <- recode(pca_res$source,
                                "IIDP,Prodo" = "IIDP")

      pca_res$pan_kbase_cause_of_death <- recode(na_if(pca_res$pan_kbase_cause_of_death, ""),
                                                 "Cerebrovascular/stroke" = "Cerebrovascular/Stroke",
                                                 "Head Trauma" = "Head trauma",
                                                 "ICH/stroke" = "Cerebrovascular/Stroke",
                                                 "Cerebral Edema (DKA)" = "Cerebral edema (DKA)",
                                                 "IHC" = "Cerebrovascular/Stroke",
                                                 .missing = "Unknown")

      # Rename columns to friendly names
      pca_res <- pca_res %>% rename("Program" = source,
                                    "Age (years)" = pan_kbase_age_years,
                                    "HbA1C percentage" = pan_kbase_hb_a1c_percentage,
                                    "Hospital stay (hours)" = pan_kbase_hospital_stay_hours,
                                    "Description of diabetes status" = pan_kbase_description_of_diabetes_status,
                                    "Cause of death" = pan_kbase_cause_of_death,
                                    "Sex" = pan_kbase_sex,
                                    "BMI" = pan_kbase_bmi,
                                    "Ethnicity" = pan_kbase_ethnicities,
                                    "Chemistry" = chemistry,
                                    "Treatment" = treatments)

      # Draw plots to color by continuous variables
      if (is.numeric(pca_res[,input$Color]) == TRUE) {
        # Make PCA plot
        p <- ggplot(pca_res, aes(x = PC1, y = PC2, color = !!sym(input$Color))) +
          geom_point(size=2) +
          scale_color_gradientn(colors = viridis(length(unique(pca_res[,input$Color])))) +
          theme_minimal() +
          labs(color = input$Color) +
          theme(text = element_text(size=14),
                plot.title = element_text(hjust = 0.5)) +
          ggtitle(paste0(input$CellType, " Cells"))
        print(p)
      }
      
      # Draw plots to color by categorical variables
      else {
        # Set up color palettes
        collections <- unique(pca_res[,input$Color]) %>% sort()
        collection_pal <- all_palette(length(collections))
        # Make PCA plot
        p <- ggplot(pca_res, aes(x = PC1, y = PC2, color = !!sym(input$Color))) + 
          geom_point(size=2) + 
          theme_minimal() +
          scale_color_manual(values = collection_pal) +
          labs(color = input$Color) +
          theme(text = element_text(size=14),
                plot.title = element_text(hjust = 0.5)) +
          ggtitle(paste0(input$CellType, " Cells"))
        print(p)
      }

    }

    output$PCAPlot <- renderPlot({
        PCA_fxn()
    })
    
    # Make PCA contributions plot ----------------------------------------------
    contribs_fxn <- function() {
      # Make sure name for "Immune (Macrophages)" is changed to "Immune"
      if (input$CellType == "Immune") {
        cell_name <- "Immune (Macrophages)"
      }
      else {
        cell_name <- input$CellType
      }
      
      # Read in RDS file with PCA results
      fname <- paste0("outputs/PCA_results/", cell_name, "_all_PCA_results.rds")
      pca_res_all <- readRDS(fname)
      
      # Plot top contributing genes
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
    
    # Make downloaders ---------------------------------------------------------
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

## Run the application ---------------------------------------------------------
shinyApp(ui = ui, server = server)
