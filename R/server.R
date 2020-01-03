#' overlapper shiny server function
#'
#' @return a shiny server function
#' 
#' @import shiny ggplot2
#' @importFrom cowplot theme_cowplot
#' @importFrom plotly renderPlotly
#' @importFrom VennDiagram venn.diagram
#' @importFrom UpSetR upset fromList
#' @importFrom ggsci pal_nejm
#' @export
overlapper.server <- function(){
  
shinyServer(function(input, output, session) {

  theme_set(theme_cowplot())
  
  getWords <- function(x){
    x <- strsplit(gsub(","," ",gsub("\n"," ",x))," ",fixed=TRUE)[[1]]
    x[x!=""]
  }
  
  nbLists <- reactiveVal(2)
  
  buildListInput <- function(i){
    name <- input[[paste0("name_list",i)]]
    items <- input[[paste0("items_list",i)]]
    en <- input[[paste0("enabled_list",i)]]
    box(width=4, id=paste0("list",i),
        textInput(paste0("name_list",i), "", value=ifelse(is.null(name), paste("List",i), name)),
        textAreaInput(paste0("items_list",i), "Items", height="300px",
                      value=ifelse(is.null(items), "", items),
                      placeholder="Enter items separated by space, commas, or linebreaks"),
        checkboxInput(paste0("enabled_list",i), "Use list", value=(is.null(en) || en))
    )
  }
  
  output$list_inputs <- renderUI({
    lapply(seq_len(nbLists()), buildListInput)
  })
  
  observeEvent(input$addlist, nbLists(isolate(nbLists())+1))
  
  background <- reactive({
    x <- getWords(input$background)
    if(input$extend_bg) x <- union(x,unique(unlist(rawLists())))
    if(length(x)==0) return(NULL)
    x
  })
  
  rawLists <- reactive({
    i <- seq_len(nbLists())
    isActive <- sapply(paste0("enabled_list",i), FUN=function(x){ input[[x]] })
    if(any(sapply(isActive,is.null))) return(list())
    i <- i[isActive]
    ll <- lapply(paste0("items_list",i), FUN=function(x) input[[x]])
    nn <- sapply(paste0("name_list",i), FUN=function(x) input[[x]])
    w <- which(nn=="")
    if(length(w)>0) nn[w] <- paste("list",w)
    names(ll) <- nn
    ll <- lapply(ll, getWords)
    ll[sapply(ll,length)>0]
  })
  
  output$badgeText_lists <- renderUI({
    ll <- rawLists()
    if(length(ll)==0) return(list(tags$span("no list"), icon("exclamation-triangle")))
    tags$span(paste0(length(activeLists()),"/",nbLists()))
  })
  
  activeLists <- reactive({
    ll <- rawLists()
    names(ll) <- breakStrings(names(ll), minSizeForBreak = 15)
    if(!is.null(background())) ll <- lapply(ll, intersect, y=background())
    ll[sapply(ll,length)>0]
  })
  
  observeEvent(input$random, {
    u <- c(LETTERS, letters)
    updateTextAreaInput(session, "background", value=paste(u,collapse=","))
    for(i in seq_len(nbLists())){
      updateTextAreaInput(session, paste0("items_list",i), 
                          value=sample(u, sample(c(8,12,15),1)))
    }
  })
  
  output$venn <- renderPlot({
    n <- length(activeLists())
    if(n<2) return(NULL)
    cols <- pal_nejm()(n)
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    grid.draw(venn.diagram( activeLists(), imagetype="svg", filename = NULL, 
                            euler.d=n<=3, scaled=n<=3, alpha=0.35,
                            cex = 1.4, cat.cex=1.3, cat.fontface="bold", 
                            col=cols, fill=cols ))
  },height=reactive(input$h_venn))
  
  output$sigres <- renderUI({
    f <- function(x){
      if(x<0.0001) return(scales::scientific(x, 3))
      e <- ceiling(-log10(x))
      round(10^e*x,3-1)/10^e
    }
    n <- length(activeLists())
    if(n==1 || is.null(background())) return("")
    if(n==2){
      return(HTML(
        paste0("<ul><b>Hypergeometric test:</b><li>over-representation p~", 
              f(overlap.prob(activeLists()[[1]], activeLists()[[2]], background())),
              "</li><li>under-representation p~",
              f(overlap.prob(activeLists()[[1]], activeLists()[[2]], background(), lower = TRUE)),
              "</li></ul>")
        ))
              
    }
    HTML(paste0( "Probability of having an intersection of all sets at least ",
                "this great (permutation test):<br/>",
             f(multintersect.test(activeLists(), background())), "(+/-0.002)"))
  })
  
  multint <- reactive({
    ll <- activeLists()
    if(is.null(ll) || length(ll)<2) return(NULL)
    multintersect(ll,universe=background())
  })
  
  output$heatmap <- renderPlotly({
    m <- multint()
    if(is.null(m)) return(m)
    plot.multintersect(m, cluster=FALSE) %>% event_register("plotly_click")
  })
  
  output$dotplot <- renderPlot({
    m <- multint()
    if(is.null(m)) return(m)
    dotplot.multintersect(m, sizeRange=c(input$minSize,input$maxSize))    
  },height=reactive(input$h_dotplot))

  output$upset <- renderPlot({
    ll <- activeLists()
    if(is.null(ll) || length(ll)<2) return(NULL)
    upset(fromList(ll))
  }, height=reactive(input$h_upset))
  
  output$selection <- renderText({
    if(is.null(multint())) return(NULL)
    s <- event_data("plotly_click")
    if (length(s) == 0) {
      "Click on a cell in the heatmap to see the intersecting features"
    }
    else {
      multint()$of[s$x, s$y]
    }
  })
  
})

}


#' overlapper
#' 
#' Launches the overlapper shiny app.
#'
#' @param ... passed to `shiny::shinyApp()`
#'
#' @return a shiny app
#' @export
#' @import shiny
overlapper <- function(...){
  shiny::shinyApp(ui=overlapper.ui(), server=overlapper.server(), ...)
}