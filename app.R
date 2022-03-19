library(shiny)
library(IntCal)

source('calibration help function.R')
##### load the prior data ######################################
# Define UI for application 
ui <- fluidPage(
  
  # Application title
  titlePanel("14C Age Calibration  of Radiocarbon Age"),
  # Sidebar with a slider input
  #sidebarLayout(
  sidebarLayout(
    # Show the input 
    sidebarPanel(
      textInput('rdata','Radiocarbon Age:',value='3000 '),
      textInput('error','error of Radiocarbon Age:',value='10 '),
      actionButton('calibrate','Calibration'),
      br(),
      br(),
      textOutput('idata'),
      br(),
      textOutput('thetamin'),
      textOutput('thetamax')
      
    ),
    # Show two tab for help and computation results
    mainPanel(
      tabsetPanel(
        id = "tabset",
        tabPanel("Computation", 
                 plotOutput("distPlot",width="642px",height = "442px"),
                 verbatimTextOutput('summary')),
        tabPanel("Help", 
                 br(),
                 br(),
                 "Under Construction!")
        
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  rdata<-eventReactive(input$calibrate,{
    c(as.numeric(input$rdata),as.numeric(input$error))
    })
  
  output$summary <- renderText({
    paste("Calibration Results:",rdata()[1],"+/-",rdata()[2])
    
  })
  output$thetamin <- renderText({
    paste("the theta min is :",Mu_To_Theta_Min(
      rdata = rdata()[1],
      error=rdata()[2],
      cc=ca.curve))
 })
  
  output$thetamax <- renderText({
    paste("the theta min is :",Mu_To_Theta_Max(
      rdata = rdata()[1],
      error=rdata()[2],
      cc=ca.curve))
  }) 
  
  output$idata <- renderText({
    paste("You Data is :",as.numeric(rdata()[1]))
    
  })
  #---------------------

  output$distPlot <- renderPlot({
    
    #---------------------

    theta.min <- Mu_To_Theta_Min(
      rdata = rdata()[1],
      error=rdata()[2],
      cc=ca.curve)
    theta.max <- Mu_To_Theta_Max(
      rdata = rdata()[1],
      error=rdata()[2],
      cc=ca.curve)
    #---------------------
    
    theta <- seq(theta.min,
                 theta.max,
                 length.out=(theta.max-theta.min)*2)

    post <- Theta.MuER(theta=theta,cc=ca.curve)
    post <- Grid_Post(post,rdata=rdata()[1],error=rdata()[2])
    
    pm <- round(sum(post$theta*post$prob),0) # mean
    psd <- round(sqrt(sum((post$theta-pm)^2*post$prob)),0) # sd
    
    eps <- 1e-5
    xpost <- post[post$prob>eps,]
    #---------------------
    PlotCurve(post=xpost,pm=pm,psd=psd,rdata=rdata()[1],error=rdata()[2])
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
