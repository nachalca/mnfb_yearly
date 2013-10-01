
library(shiny)
bird.tot <- read.csv('bird_yeartotal.csv', header=T)
  
shinyUI(bootstrapPage(
  headerPanel("Bird Count Trend"), 
  
  sidebarPanel( 
    selectInput(inputId = "specie",
                label = "Select Bird specie:",
                choices = levels( bird.tot$abbrev )
    )

  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Average", plotOutput("plot1")),
      tabPanel("Total", plotOutput("plot2")),
      tabPanel("Summary", tableOutput("summary"))
      
    )
  )
  
))
