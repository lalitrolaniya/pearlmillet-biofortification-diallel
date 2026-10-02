# =====================================================================
# Shiny app: Pooled diallel analysis (Griffing Method 2, Model 1)
# Companion app for Punia et al. (2026) Theor Appl Genet 139:265
# https://doi.org/10.1007/s00122-026-05378-4
#
# Accepts REPLICATED plot-level data (recommended) or entry means.
# Any number of parents, environments, replications and traits.
# =====================================================================

library(shiny)
source("griffing_engine.R")

ui <- fluidPage(
  titlePanel("Pooled diallel analysis - Griffing Method 2, Model 1"),
  sidebarLayout(
    sidebarPanel(width = 3,
      fileInput("file", "Data file (CSV)", accept = ".csv"),
      helpText(strong("Replicated data (recommended): "),
               "columns Genotype, P1, P2, Env, Rep, then one column",
               "per trait; one row per plot. The app computes the",
               "pooled error itself: no other input is needed."),
      helpText(strong("Entry means: "), "same columns without Rep;",
               "one row per entry per environment. Then set the",
               "number of replications and, for significance tests,",
               "upload pooled error mean squares (CSV: Trait,MS)."),
      uiOutput("means_ui"),
      uiOutput("trait_ui"),
      downloadButton("dl3", "ANOVA CSV"),
      downloadButton("dl4", "GCA CSV"),
      downloadButton("dls", "SCA CSV")
    ),
    mainPanel(width = 9,
      verbatimTextOutput("design"),
      tabsetPanel(
        tabPanel("ANOVA & genetic parameters", tableOutput("t3")),
        tabPanel("GCA effects & stability", tableOutput("t4")),
        tabPanel("SCA effects", tableOutput("ts")),
        tabPanel("About", br(),
          p("Combining ability analysis for a half diallel",
            "(parents + F1s, no reciprocals) pooled over",
            "environments, following Griffing (1956) Method 2,",
            "Model 1. GCA is tested against GCA x E, SCA against",
            "SCA x E, and the interactions against the pooled error",
            "from the replication-level RCBD ANOVA."),
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
  raw <- reactive({
    req(input$file)
    read.csv(input$file$datapath, check.names = FALSE)
  })
  replicated <- reactive("Rep" %in% names(raw()))

  # Extra inputs only when entry means (no Rep column) are uploaded
  output$means_ui <- renderUI({
    req(input$file)
    if (replicated()) return(NULL)
    tagList(
      numericInput("reps", "Replications per environment", 3, min = 2),
      numericInput("errdf", "Pooled error df", NA, min = 1),
      fileInput("errfile", "Pooled error mean squares (CSV: Trait,MS)",
                accept = ".csv")
    )
  })

  prep <- reactive({
    d <- raw()
    req_cols <- c("Genotype", "P1", "P2", "Env")
    validate(need(all(req_cols %in% names(d)),
      paste("Missing required columns:",
            paste(setdiff(req_cols, names(d)), collapse = ", "))))
    if (replicated()) {
      prep_from_plot_data(d)
    } else {
      traits <- names(d)[!(names(d) %in% req_cols) & sapply(d, is.numeric)]
      ems <- if (!is.null(input$errfile)) {
        e <- read.csv(input$errfile$datapath)
        setNames(e[[2]], e[[1]])
      } else setNames(numeric(0), character(0))
      list(means = d, error_ms = ems,
           error_df = if (is.null(input$errdf) || is.na(input$errdf))
                        NA else input$errdf,
           r = if (is.null(input$reps)) 3 else input$reps,
           traits = traits)
    }
  })

  output$trait_ui <- renderUI({
    req(input$file)
    checkboxGroupInput("traits", "Traits", prep()$traits,
                       selected = prep()$traits)
  })

  output$design <- renderText({
    req(input$file)
    p  <- prep()
    np <- length(unique(c(p$means$P1, p$means$P2)))
    ne <- length(unique(p$means$Env))
    ng <- length(unique(p$means$Genotype))
    sprintf(paste0("Detected design: %d parents | %d entries | %d ",
                   "environment(s) | %d replication(s) | error df: %s | ",
                   "input: %s"),
            np, ng, ne, p$r,
            ifelse(is.na(p$error_df), "-", p$error_df),
            ifelse(replicated(), "replicated plot data (error MS computed)",
                   "entry means"))
  })

  tabs <- reactive({
    req(input$traits)
    p <- prep()
    build_tables(p$means, input$traits, p$r, p$error_ms, p$error_df)
  })
  output$t3 <- renderTable(tabs()$table3, rownames = TRUE, striped = TRUE)
  output$t4 <- renderTable(tabs()$table4, rownames = TRUE, striped = TRUE)
  output$ts <- renderTable(tabs()$sca, striped = TRUE)
  output$dl3 <- downloadHandler("Table_ANOVA_combining_ability.csv",
    function(f) write.csv(tabs()$table3, f))
  output$dl4 <- downloadHandler("Table_GCA_effects.csv",
    function(f) write.csv(tabs()$table4, f))
  output$dls <- downloadHandler("Table_SCA_effects.csv",
    function(f) write.csv(tabs()$sca, f, row.names = FALSE))
}

shinyApp(ui, server)
