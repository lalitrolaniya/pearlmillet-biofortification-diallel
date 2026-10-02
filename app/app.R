# =====================================================================
# Shiny app: Pooled diallel analysis (Griffing Method 2, Model 1)
# Companion app for Punia et al. (2026) Theor Appl Genet 139:265
# https://doi.org/10.1007/s00122-026-05378-4
#
# Runs locally with shiny::runApp("app"), from GitHub with
# shiny::runGitHub(), or in the browser via shinylive (GitHub Pages).
# =====================================================================

library(shiny)
source("griffing_engine.R")

ui <- fluidPage(
  titlePanel("Pooled diallel analysis - Griffing Method 2, Model 1"),
  sidebarLayout(
    sidebarPanel(width = 3,
      fileInput("file", "Entry x environment means (CSV)",
                accept = ".csv"),
      helpText("Columns: Genotype, P1, P2, Env, then one column per",
               "trait. One row per entry per environment",
               "(see data/input_format.md). For the pearl millet",
               "dataset of the paper, use Supplementary Table S1",
               "from the article page."),
      numericInput("reps", "Replications per environment", 3, min = 2),
      numericInput("errdf", "Pooled error df", 216, min = 1),
      fileInput("errfile",
                "Pooled error mean squares (optional CSV: Trait,MS)",
                accept = ".csv"),
      helpText("Needed only for significance of effects and of the",
               "interaction terms."),
      uiOutput("trait_ui"),
      downloadButton("dl3", "Table 3 CSV"),
      downloadButton("dl4", "Table 4 CSV"),
      downloadButton("dls", "SCA CSV")
    ),
    mainPanel(width = 9,
      tabsetPanel(
        tabPanel("Table 3: ANOVA & genetic parameters",
                 tableOutput("t3")),
        tabPanel("Table 4: GCA effects & stability",
                 tableOutput("t4")),
        tabPanel("SCA effects", tableOutput("ts")),
        tabPanel("About", br(),
          p("This app performs combining ability analysis for a",
            "half diallel (parents + F1s, no reciprocals) pooled",
            "over environments, following Griffing (1956) Method 2,",
            "Model 1. GCA is tested against GCA x E, SCA against",
            "SCA x E, and the interactions against the pooled error."),
          p(strong("Citation: "),
            "Punia M, Sharma LD, Gothwal DK, Kajla SL, Rolaniya LK,",
            "Sharma V, Jat RL (2026) Combining ability and gene",
            "action for grain yield and biofortification traits in",
            "pearl millet. Theoretical and Applied Genetics 139:265.",
            a("doi:10.1007/s00122-026-05378-4",
              href = "https://doi.org/10.1007/s00122-026-05378-4")),
          p(em("In memory of Dr. Monika Punia, who conceived and",
               "led this research.")))
      )
    )
  ),
  tags$footer(style = "margin-top:2em;color:#666;font-size:85%;",
    "Based on the doctoral research of the late Dr. Monika Punia. ",
    "Published in her memory, in fulfilment of her wish to see it in print.")
)

server <- function(input, output, session) {
  dat <- reactive({
    req(input$file)
    read.csv(input$file$datapath, check.names = FALSE)
  })
  err <- reactive({
    if (!is.null(input$errfile)) {
      e <- read.csv(input$errfile$datapath)
      setNames(e[[2]], e[[1]])
    } else setNames(numeric(0), character(0))
  })
  traits <- reactive(setdiff(names(dat()),
                             c("Genotype", "P1", "P2", "Env")))
  output$trait_ui <- renderUI(
    checkboxGroupInput("traits", "Traits", traits(),
                       selected = traits()))
  tabs <- reactive({
    req(input$traits)
    build_tables(dat(), input$traits, input$reps, err(), input$errdf)
  })
  output$t3 <- renderTable(tabs()$table3, rownames = TRUE, striped = TRUE)
  output$t4 <- renderTable(tabs()$table4, rownames = TRUE, striped = TRUE)
  output$ts <- renderTable(tabs()$sca, striped = TRUE)
  output$dl3 <- downloadHandler("Table3_combining_ability.csv",
    function(f) write.csv(tabs()$table3, f))
  output$dl4 <- downloadHandler("Table4_GCA_effects.csv",
    function(f) write.csv(tabs()$table4, f))
  output$dls <- downloadHandler("SCA_effects.csv",
    function(f) write.csv(tabs()$sca, f, row.names = FALSE))
}

shinyApp(ui, server)
