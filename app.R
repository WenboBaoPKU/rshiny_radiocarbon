library(shiny)
library(IntCal)

source('calibration help function.R')
##### load the prior data ######################################
# intcal<-ccurve(1)
# ca.curve<-as.data.frame(intcal)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("14C Age Cl"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    
    sidebarPanel(
      
      textInput('rdata','Radiocarbon Age:',value='3000 '),
      textInput('error','error of Radiocarbon Age:',value='10 '),
      actionButton('calibrate','Calibration'),
      textOutput('idata'),
      br(),
      br(),
      textOutput('thetamin')
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot",width="642px",height = "442px"),
      verbatimTextOutput('summary')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  rdata<-eventReactive(input$calibrate,{
    c(as.numeric(input$rdata),as.numeric(input$error))
    })
  
  output$summary<- renderText({
    paste("Calibration Results:",rdata()[1],"+/-",rdata()[2])
    
  })
  output$idata <- renderText({
    paste("You got :",2*as.numeric(rdata()[1]))
    
  })
  
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
