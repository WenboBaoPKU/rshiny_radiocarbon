library(shiny)
#library(IntCal)
#setwd("D:/bayesian/openbugs/radiocarbon calibration/shiny_radiocarbon/Radiocarbon_Shiny_Test")
source('calibration help function.R')
##### load the prior data ######################################
# Define UI for application 
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "darkly"),
  
  # Application title
  titlePanel("14C Age Calibration  of Radiocarbon Age"),
  br(),
  # Sidebar with a slider input
  #sidebarLayout(
  sidebarLayout(
    # Show the input 
    sidebarPanel(
      textInput('rdata','Radiocarbon Age:',value='4340 '),
      textInput('error','error of Radiocarbon Age:',value='30 '),
      actionButton('calibrate','Calibration',class = "btn-primary btn"),
      br(),
      br(),
      textOutput('idata'),
      br(),
      textOutput('thetamin'),
      textOutput('thetamax'),
      width = 3),
    # Show two tab for help and computation results
    mainPanel(
      tabsetPanel(
        id = "tabset",
        tabPanel("Computation", 
                 plotOutput("distPlot",width="642px",height = "442px"),
                 br(),
                 verbatimTextOutput('summary')),
        tabPanel("Help", 
                 br(),
                 br(),
                 "Under Construction!!!")
      )
    ,width = 5)
  )
)


# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  #rdata <- reactive(c(as.numeric(input$rdata),as.numeric(input$error)))
  rdata<-eventReactive(input$calibrate,{
    c(as.numeric(input$rdata),as.numeric(input$error))
  })
  #-------preparetion for computation -------
  theta.min <- eventReactive(input$calibrate,{
    Mu_To_Theta_Min(rdata = rdata()[1],error=rdata()[2],cc=ca.curve)
  })
  theta.max <- eventReactive(input$calibrate,{
    Mu_To_Theta_Max(rdata = rdata()[1],error=rdata()[2],cc=ca.curve)
  })
  #------ posterior density---------------
  theta <- eventReactive(input$calibrate,{
    seq(theta.min(),theta.max(),length.out=(theta.max()-theta.min())*2)
  })
  
  post1 <- eventReactive(input$calibrate,{
    Theta.MuER(theta=theta(),cc=ca.curve)
  })
  post <- eventReactive(input$calibrate,{
    Grid_Post(post=post1(),rdata=rdata()[1],error=rdata()[2])
  })
  #---- summary of posterior -----------------
  hpd <- eventReactive(input$calibrate,{
    HPD(calib = post()[,c(1,5)],rounded = 2,prob = 0.682)
  })
  
  pm <- eventReactive(input$calibrate,{
    round(sum(post()$theta*post()$prob),0)
  }) # mean
  psd <- eventReactive(input$calibrate,{
    round(sqrt(sum((post()$theta-pm())^2*post()$prob)),0)
  }) # sd
  #-----show the result to interface ---------

  observeEvent(input$calibrate,{
    updateTabsetPanel(session, "tabset", selected = paste0("Computation"))
  })
  
  
  output$summary <- renderText({
    
    hpdtext <- ''
    if(is.null(hpd())){
      hpdtext <- paste0(hpdtext,'There is no results!!!!!!!!!!!!!!!!!!!')
     
    }else{
      for(i in 1:dim(hpd())[1]){
        hpdtext <- paste0(hpdtext,round(hpd()[i,1],digits = 0),"(",round(hpd()[i,3],1),'%)',round(hpd()[i,2],0),'\n ')
      }
    }
    #---------------------
    paste('Calibration Results:\n',
          'Mean: ',pm(),'; SD: ',psd(),'\n',
          'HPD for 68.2% range: \n' , hpdtext,'\n ')
    
  })

  output$thetamin <- renderText({
    paste("the theta min is :",theta.min())
  })
  output$thetamax <- renderText({
    paste("the theta min is :",theta.max())
  }) 
  output$idata <- renderText({
    paste("You Data is :",as.numeric(rdata()[1]),'+/-',as.numeric(rdata()[2]))
  })


  output$distPlot <- renderPlot({
    PlotCurve(post=post(),pm=pm(),psd=psd(),rdata=rdata()[1],error=rdata()[2],hpd=hpd())
  })
}
# Run the application 
shinyApp(ui = ui, server = server)
