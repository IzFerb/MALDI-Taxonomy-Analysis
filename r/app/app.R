# Gums Taxonomy Analysis — point-and-click interface
#
# Launch from project root:
#   source("START_APP.R")        # recommended
#   source("r/launch_app.R")

# Paths are set by launch_app.R via options; fallback if Run App is clicked directly
project_root <- getOption("gums.project_root")
r_code_dir   <- getOption("gums.r_code_dir")

if (is.null(project_root) || is.null(r_code_dir)) {
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  if (basename(cwd) == "app" && basename(dirname(cwd)) == "r") {
    r_code_dir   <- normalizePath(file.path(cwd, "..", "R"))
    project_root <- normalizePath(file.path(cwd, "../.."))
  } else if (dir.exists(file.path(cwd, "r", "R"))) {
    r_code_dir   <- normalizePath(file.path(cwd, "r", "R"))
    project_root <- cwd
  } else {
    stop(
      "Could not locate r/R/.\n",
      "Working directory: ", cwd, "\n",
      "Fix: setwd('/Users/isaac/Desktop/art_institute_code'); source('r/launch_app.R')"
    )
  }
}

pipeline_file <- file.path(r_code_dir, "run_pipeline.R")
if (!file.exists(pipeline_file)) {
  stop("Pipeline not found at: ", pipeline_file)
}

message("Using project root: ", project_root)
source(pipeline_file, local = FALSE)

library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("MALDI Gum Taxonomy Analysis"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Data"),
      helpText("Point to the folder containing sample subfolders (e.g. sample_data/) or the project root."),
      textInput("data_dir", "Data folder", value = project_root),
      textInput("output_dir", "Output folder", value = file.path(project_root, "r", "output")),
      hr(),
      h4("Key settings"),
      numericInput("mz_min", "m/z minimum", value = 800, min = 0),
      numericInput("mz_max", "m/z maximum", value = 1300, min = 0),
      numericInput("max_shift", "Max spectral shift (icoshift)", value = 34, min = 1),
      numericInput("peak_prominence", "Peak prominence factor", value = 10, min = 1),
      selectInput("alignment_target", "Alignment reference", choices = c("median", "average")),
      hr(),
      actionButton("run_btn", "Run analysis", class = "btn-primary btn-lg", width = "100%"),
      br(), br(),
      uiOutput("status_ui")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Overview", verbatimTextOutput("log_text")),
        tabPanel("PCA plots", uiOutput("plot_gallery")),
        tabPanel("Results table", DT::dataTableOutput("results_table"))
      )
    )
  )
)

server <- function(input, output, session) {
  results <- reactiveVal(NULL)
  log_text <- reactiveVal("Ready. Click 'Run analysis' to start.\n")

  observeEvent(input$run_btn, {
    log_text("Running analysis...\n")
    output$status_ui <- renderUI(tags$p(class = "text-info", "Analysis in progress..."))

    data_dir <- input$data_dir
    if (basename(data_dir) != "sample_data" && dir.exists(file.path(data_dir, "sample_data"))) {
      data_dir <- file.path(data_dir, "sample_data")
    }

    cfg <- list(
      data_dir = data_dir,
      output_dir = input$output_dir,
      mz_min = input$mz_min,
      mz_max = input$mz_max,
      max_shift = as.integer(input$max_shift),
      peak_prominence_factor = input$peak_prominence,
      alignment_target = input$alignment_target,
      reference_groups = c(4),
      acacia_group_id = 4,
      historical_group_id = 2,
      pigmented_group_id = 3,
      exclude_group_id = 2,
      save_alignment_plot = TRUE,
      save_snv_plot = TRUE
    )

    res <- tryCatch(
      run_maldi_pipeline(config = cfg, verbose = TRUE),
      error = function(e) e
    )

    if (inherits(res, "error")) {
      log_text(paste0("ERROR: ", res$message, "\n"))
      output$status_ui <- renderUI(tags$p(class = "text-danger", "Analysis failed. See log for details."))
      results(NULL)
      return()
    }

    results(res)
    log_text(paste0(
      "Analysis complete.\n",
      "Samples kept: ", nrow(res$feature_matrix_filtered), "\n",
      "Acacia samples: ", length(res$acacia_classes %||% character()), "\n",
      "Output: ", normalizePath(res$output_dir), "\n"
    ))
    output$status_ui <- renderUI(tags$p(class = "text-success", "Analysis finished successfully."))
  })

  output$log_text <- renderText(log_text())

  output$plot_gallery <- renderUI({
    res <- results()
    if (is.null(res)) return(p("No plots yet."))
    plot_list <- res$plots %||% list()
    extra <- list()
    if (!is.null(res$historical_plot)) extra <- c(extra, list(res$historical_plot))
    if (!is.null(res$pigmented_plot)) extra <- c(extra, list(res$pigmented_plot))
    plot_list <- c(plot_list, extra)
    if (length(plot_list) == 0L) return(p("No plots generated."))
    lapply(seq_along(plot_list), function(i) {
      plotOutput(paste0("pca_plot_", i), height = "500px")
    })
  })

  observe({
    res <- results()
    req(res)
    plot_list <- res$plots %||% list()
    extra <- list()
    if (!is.null(res$historical_plot)) extra <- c(extra, list(res$historical_plot))
    if (!is.null(res$pigmented_plot)) extra <- c(extra, list(res$pigmented_plot))
    plot_list <- c(plot_list, extra)
    for (i in seq_along(plot_list)) {
      local({
        idx <- i
        p <- plot_list[[idx]]
        output[[paste0("pca_plot_", idx)]] <- renderPlot({
          print(p)
        })
      })
    }
  })

  output$results_table <- DT::renderDataTable({
    res <- results()
    if (is.null(res) || is.null(res$acacia_classes)) {
      return(data.frame(Message = "Run the analysis to see sample classes."))
    }
    classes_file <- file.path(res$output_dir, "acacia_sample_classes.csv")
    if (file.exists(classes_file)) {
      return(DT::datatable(utils::read.csv(classes_file), rownames = FALSE))
    }
    DT::datatable(data.frame(Message = "No results file found."), rownames = FALSE)
  })
}

shinyApp(ui, server)
