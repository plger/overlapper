#' UI for the overlapper app
#'
#' @return a shiny UI
#' @import shiny shinydashboard shinycssloaders
#' @importFrom plotly plotlyOutput
#' @export
overlapper.ui <- function(){

shinyUI( dashboardPage(

  dashboardHeader(title="overlapper"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Lists", tabName="tab_lists", badgeLabel=uiOutput("badgeText_lists"), badgeColor="blue"),
      menuItem("Background", tabName="tab_bg"),
      menuItem("Venn", tabName="tab_venn"),
      menuItem("Heatmap", tabName="tab_heatmap"),
      menuItem("Dotplot", tabName="tab_dotplot"),
      menuItem("UpSet", tabName="tab_upset"),
      menuItem("Help", tabName="tab_help")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem("tab_lists",
        tags$p(style="text-align: right;", 
               actionButton("addlist", "Add list"), 
               actionButton("random","Random"), 
               checkboxInput("icase", "Convert all to uppercase", value=FALSE)),
        uiOutput("list_inputs")
      ),
      tabItem("tab_bg", box(width=12, 
        helpText(paste( "The background/universe is critical for significance,", 
                        "and typically represents the items that could have",
                        "made it into any of your sets.")),
        textAreaInput("background", "Background:", height="450px",
                      placeholder="Enter items separated by space, commas, or linebreaks"),
        checkboxInput("extend_bg", "Extend background to include all list items", TRUE),
        helpText("If not checked, any list item absent from the background will be removed from the lists."))),
      tabItem("tab_venn", box(width=12, plotOutput("venn", height="auto"),
                              tags$hr(),
                              numericInput("h_venn", "Plot height", 300, min=200, max=1000, step=50),
                              withSpinner(htmlOutput("sigres")))), 
      tabItem("tab_heatmap", box(width=12, withSpinner(plotlyOutput("heatmap", height="auto")),
                                 tags$hr(),
                              tags$b("Selected intersection:"), textOutput("selection"))), 
      tabItem("tab_dotplot", box(width=12, plotOutput("dotplot", height="auto"),
                                 tags$hr(),
                                 column(4, numericInput("h_dotplot", "Plot height", 400, min=200, max=1000, step=50)),
                                 column(4, numericInput("minSize", "Min dot size", 0, min=0, max=20, step=2)),
                                 column(4, numericInput("maxSize", "Max dot size", 20, min=2, max=100, step=2))
                                 )), 
      tabItem("tab_upset", box(width=12, withSpinner(plotOutput("upset", height="auto")),
                               tags$hr(),
                               numericInput("h_upset", "Plot height", 500, min=200, max=1000, step=50))),
      tabItem("tab_help", box(width=12, title="Help", 
                              tags$p("The overlapper app is meant to browse and test for significance the overlap between lists (e.g. of genes). The usage is simple:"),
                              tags$ol(lapply(list("In the first tab ('Lists'), enter the content of each of your lists (as many as you like) and eventually give them names.",
                                           "If you're interested in enrichment and significance of overlaps, you need to specify the universe or background in the second ('Background') tab. The appropriate background depends on the question, but it's usually the set of everything that could potentially have ended up in your overlap (for instance, if comparing two differential expressions, it'd be the intersection of the tested genes in both experiments).",
                                           list("Then you can go to any of the other tabs to view a representation:",
                                           tags$ul(lapply(list("the Venn diagram displays the size of the sets and their various overlaps; the significance of the overlap of all sets will also be reported.",
                                                        "the heatmap is interactive; the numbers shown are the number of items, but if you hover on a cell you'll see the enrichment/p-value, and if you click on a cell you'll see its content (i.e. the items of the intersection) displayed below the plot.",
                                                        "the dotplot is similar but represents significance, enrichment and size of the overlap in a static plot (better for papers and such).",
                                                        "the upset plot is an alternative to the Venn diagram to represent sets and their overlaps, which is especially better when you have many sets."
                                                        ), tags$li))) 
                                           ), tags$li))),
                          box(width=12, tags$p(paste("overlapper version",packageVersion("overlapper"),"- Pierre-Luc Germain")),
                              tags$p("In case of bugs/problems, create an issue here:"), 
                              tags$a(href="https://github.com/plger/overlapper/issues","github.com/plger/overlapper"))
      )
    )
  )
))

}