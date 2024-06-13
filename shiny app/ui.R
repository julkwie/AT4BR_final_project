# Loading required libraries
library(shiny)      # For creating the Shiny web application


# Define UI for the Shiny application
fluidPage(
  titlePanel("FASTA File Processing"),
  
  # Sidebar layout
  sidebarLayout(
    sidebarPanel(
      align = "center", # nice and centered :)
      
      # File input section
      h4("Choose the FASTA file", style = "font-weight: bold;"),
      fileInput("file", ""),
      div(style = "height: 20px;"), # to space out sections 
      
      h4("Choose what you want to do with this file:", style = "font-weight: bold;"),
      div(style = "height: 40px;"),
      
      # Action button: Clustal Omega dendrogram
      h5("Draw me a dendrogram with use of Clustal Omega alignment"),
      actionButton("clustal", "Clustal Omega dendrogram!", icon = icon("sitemap"), style="color: #fff; background-color: #808080"),
      div(style = "height: 40px;"),
      
      # Action button: General sequence information
      h5("Show me some general info about sequences inside"),
      actionButton("info", "Inform me in general!", icon = icon("gears"), style="color: #fff; background-color: #808080"),
      div(style = "height: 40px;"),
      
      # Action button: Table of organisms and sequences
      h5("Show me a table of all organisms and sequences inside"),
      actionButton("table", "What's inside?", icon = icon("bars"), style="color: #fff; background-color: #808080"),
      div(style = "height: 40px;"),
      
      # Action button: Detailed sequence information
      h5("Show me detailed information about sequences inside"),
      actionButton("info2", "More info about the details!", icon = icon("file-lines"), style="color: #fff; background-color: #808080"),
      div(style = "height: 40px;"),
      
    ),
    
    # Main panel for displaying dynamic outputs
    mainPanel(
      align = "center",
      uiOutput("dynamic_output")
    )
  )
)
