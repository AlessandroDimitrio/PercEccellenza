#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(shinydashboard)
library(TCGAbiolinks)
library(survminer)
library(survival)
library(pROC)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(DT)
library(shinyjs)

# merged_data, stratified_data, etc. should be available in the global environment

merged_data <- readRDS("data/merged_data.rds")
stratified_data <- readRDS("data/stratified_data.rds")

# UI
ui <- dashboardPage(
  title = "Percorso di Eccellenza",
  skin = "purple",
  dashboardHeader(
    title = div(h3('Percorso di Eccellenza', style="margin: 0; margin-top: 10px; font-size: 12px; font-weight: bold;"), h3('Alessandro Dimitrio', style="margin: 0;font-size: 12px;"))
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduzione", tabName = "int", icon = icon("house")),
      menuItem("Differential Expression", tabName = "de", icon = icon("dna")),
      menuItem("Analisi AUC", tabName = "roc", icon = icon("chart-line")),
      menuItem("Regressione Binomiale", tabName = "reg", icon = icon("chart-bar"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
      tags$head(tags$link(rel="shortcut icon", href="favicon.ico"))
    ),
    tabItems(
      
      # Introduzione Tab
      tabItem(
        tabName = "int",
        fluidRow(
          box(
            width = 12,
            title = "Analisi differenziale dei miRNA",
            status = "primary",
            solidHeader = TRUE,
            img(src = "intro.png", alt = "Descrizione dell'immagine", class="image_intro"),
            p("Ciao Giovanni, spero che questo tool possa essere utile per analizzare insieme i dati, senza passare per R."),
          )
        )
      ),
      
      # ROC Analysis Tab
      tabItem(tabName = "roc",
              fluidRow(
                box(
                  width = 12,
                  title = "ROC Curve",
                  selectizeInput("mirna_select", "Scegli i miRNA:", 
                                 choices = unique(merged_data$miRNA_ID), 
                                 multiple = TRUE,
                                 options = list(maxItems = 9)),
                  # actionButton("update_roc", "Aggiorna Grafici ROC"),
                  plotOutput("roc_plots", height = "800px")
                )
              )
      ),
      
      # Differential Expression Tab
      tabItem(tabName = "de",
              fluidRow(
                box(
                  width = 12,
                  title = "Differential Expression Analysis (DESeq2)",
                  selectizeInput("design_elements", "Seleziona gli elementi da considerare nell'analisi:",
                                 choices = c("Muscle stratification" = "muscle_stratification", 
                                             "Gender" = "gender", 
                                             "Age" = "age",
                                             "Cancer stage" = "cancer_stage"),
                                 multiple = TRUE,
                                 selected = "muscle_stratification",
                                 options = list(placeholder = 'Seleziona uno o più elementi')),
                  selectInput("padj_select", "Scegli il cut-off di significatività (valore p_adjusted):",
                              choices = c("0.001", "0.025", "0.05", "0.75")),
                  checkboxInput("use_quartiles", "Per lo SMA, stratifico per quartili? (lascia non spuntato per i terzili)", value = TRUE),
                  actionButton("run_deseq", "Fai la differential expression"),
                  hidden(
                    div(id = "loading",
                        style = "color: black;",
                        p("Analisi in corso... Per favore, attendi...")
                    )
                  )
                )
              ),
              fluidRow(
                box(
                  width = 6,
                  title = "Volcano Plot",
                  plotOutput("volcano_plot")
                ),
                box(
                  width = 6,
                  fluidRow(
                    column(8, 
                           h4(
                             style = "margin-top: 0;",
                             tags$span("Risultati: "), 
                             tags$span(textOutput("results_count", inline = TRUE))
                           )
                    ),
                    column(4, downloadButton("download_results", "Scarica", class = "pull-right"))
                  ),
                  checkboxInput("only_sig", "Mostra solo i significativi", value = TRUE),
                  DTOutput("deseq_results_table")
                )
              )
      ),
      tabItem(tabName = "reg",
              fluidRow(
                box(
                  width = 12,
                  title = "Regressione Binomiale",
                  selectizeInput("reg_vars", "Seleziona le variabili indipendenti per la regressione:", 
                                 choices = c("age", "gender", "Muscle", "SAT", "VAT"), 
                                 multiple = TRUE),
                  selectizeInput("reg_mirna", "Includi i miRNA per la regressione:", 
                                 choices = unique(merged_data$miRNA_ID), 
                                 multiple = TRUE,
                                 options = list(maxItems = 5)),
                  checkboxInput("use_quartiles_reg", "Per lo SMA, stratifico per quartili? (lascia non spuntato per i terzili)", value = TRUE),
                  actionButton("run_regression", "Esegui Regressione"),
                  DTOutput("reg_results_table"),
                  verbatimTextOutput("reg_results_raw")
                )
              )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  output$download_results <- downloadHandler(
    filename = function() {
      paste("deseq_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(res_df(), file, row.names = FALSE)
    }
  )
  
  
  
  # ROC Analysis
  output$roc_plots <- renderPlot({
    req(input$mirna_select)
    
    par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
    
    for(mirna in input$mirna_select) {
      mirna_data <- merged_data %>%
        filter(miRNA_ID == mirna) %>%
        mutate(muscle_stratification = ifelse(Muscle > median(Muscle, na.rm = TRUE), "High", "Low"))
      
      roc_obj <- roc(mirna_data$muscle_stratification, mirna_data$expression)
      auc_value <- auc(roc_obj)
      
      plot(roc_obj, main = paste("ROC for", mirna), col = "blue", lwd = 2)
      abline(a = 0, b = 1, lty = 2, col = "red")
      text(0.5, 0.1, paste("AUC:", round(auc_value, 3)), 
           col = "black", cex = 1.5, font = 2)
    }
  })
  
  # Differential Expression Analysis
  observeEvent(input$run_deseq, {
    
    shinyjs::show("loading")
    
    design_formula <- paste("~", paste(input$design_elements, collapse = " + "))
    
    # Stratify data based on user selection
    stratify_data <- reactive({
      merged_data %>%
        group_by(gender) %>%
        mutate(
          muscle_stratification = case_when(
            Muscle > quantile(Muscle, ifelse(input$use_quartiles, 0.75, 2/3), na.rm = TRUE) ~ "High",
            Muscle < quantile(Muscle, ifelse(input$use_quartiles, 0.25, 1/3), na.rm = TRUE) ~ "Low",
            TRUE ~ NA_character_
          )
        ) %>%
        ungroup() %>%
        filter(!is.na(muscle_stratification))
    })
    
    # Prepare data for DESeq2
    deseq_data <- stratify_data() %>%
      select(submitter_id, miRNA_ID, expression, muscle_stratification, gender, age, cancer_stage) %>%
      group_by(miRNA_ID, submitter_id) %>%
      summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop")  %>%
      pivot_wider(names_from = submitter_id, values_from = expression) %>%
      column_to_rownames("miRNA_ID")
    
    deseq_metadata <- stratify_data() %>%
      select(submitter_id, muscle_stratification, gender, age,cancer_stage) %>%
      distinct() %>%
      column_to_rownames("submitter_id")
    
    common_samples <- intersect(colnames(deseq_data), rownames(deseq_metadata))
    deseq_data <- deseq_data[, common_samples]
    deseq_metadata <- deseq_metadata[common_samples, ]
    
    # Run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = deseq_data,
                                  colData = deseq_metadata,
                                  design = as.formula(design_formula))
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("muscle_stratification", "Low", "High"))
    res <- lfcShrink(dds, contrast = c("muscle_stratification", "Low", "High"), res=res, type = 'normal')
    
    # Generate Volcano Plot
    output$volcano_plot <- renderPlot({
      EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'padj',
                      title = 'Differential Expression Analysis',
                      pCutoff = as.numeric(input$padj_select),
                      FCcutoff = 0.0001,
                      subtitle = input$design_select)
    })
    
    res_df <- reactive({
      req(res)
      df <- as.data.frame(res) %>%
        rownames_to_column("miRNA") %>%
        arrange(padj)
      
      if(input$only_sig) {
        df <- df %>% filter(padj < as.numeric(input$padj_select))
      }
      df
    })
    
    
    output$results_count <- renderText({
      if(is.null(res_df())) {
        return("0")
      } else {
        return(nrow(res_df()))
      }
    })
    
    
    
    output$deseq_results_table <- renderDT({
      datatable(res_df(),
                options = list(pageLength = 10, scrollX = TRUE),
                filter = 'top',
                selection = 'single')
    })
    
    output$download_results <- downloadHandler(
      filename = function() {
        paste("deseq_results_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(res_df(), file, row.names = FALSE)
      }
    )
    
    shinyjs::hide("loading")
  })
  
  format_regression_results <- function(model) {
    summary_model <- summary(model)
    coef_table <- as.data.frame(summary_model$coefficients)
    coef_table$Variable <- rownames(coef_table)
    coef_table <- coef_table[, c("Variable", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
    names(coef_table) <- c("Variabile", "Coefficiente", "Errore Std", "Valore z", "p-value")
    coef_table$Coefficiente <- round(coef_table$Coefficiente, 4)
    coef_table$`Errore Std` <- round(coef_table$`Errore Std`, 4)
    coef_table$`Valore z` <- round(coef_table$`Valore z`, 4)
    coef_table$`p-value` <- format.pval(coef_table$`p-value`, digits = 4)
    return(coef_table)
  }
  
  
  
  # Regressione Binomiale
  observeEvent(input$run_regression, {
    req(input$reg_vars)
    
    # Prepara i dati per la regressione
    regression_data <- merged_data %>%
      group_by(gender) %>%
      mutate(
        muscle_stratification = case_when(
          Muscle > quantile(Muscle, ifelse(input$use_quartiles_reg, 0.75, 2/3), na.rm = TRUE) ~ "High",
          Muscle < quantile(Muscle, ifelse(input$use_quartiles_reg, 0.25, 1/3), na.rm = TRUE) ~ "Low",
          TRUE ~ NA_character_
        )
      ) %>%
      ungroup() %>%
      filter(!is.na(muscle_stratification)) %>%
      mutate(muscle_stratification = as.factor(muscle_stratification)) %>%
      select(muscle_stratification, all_of(input$reg_vars), miRNA_ID, expression)
    
    
    # Se sono stati selezionati dei miRNA, aggiungi le loro espressioni al dataset
    if (length(input$reg_mirna) > 0) {
      mirna_data <- regression_data %>%
        filter(miRNA_ID %in% input$reg_mirna) %>%
        pivot_wider(names_from = miRNA_ID, values_from = expression, names_prefix = "miRNA_")
      
      # Rinomina le colonne dei miRNA per renderle valide come nomi di variabili
      names(mirna_data) <- make.names(names(mirna_data))
      
      regression_data <- regression_data %>%
        select(-miRNA_ID, -expression) %>%
        distinct() %>%
        left_join(mirna_data, by = names(regression_data)[names(regression_data) != "miRNA_ID" & names(regression_data) != "expression"])
    } else {
      regression_data <- regression_data %>%
        select(-miRNA_ID, -expression) %>%
        distinct()
    }
    
    # Crea la formula per la regressione
    mirna_vars <- names(regression_data)[grep("^miRNA_", names(regression_data))]
    formula_vars <- c(input$reg_vars, mirna_vars)
    formula_str <- paste("muscle_stratification ~", paste(formula_vars, collapse = " + "))
    
    # Esegui la regressione
    tryCatch({
      model <- glm(as.formula(formula_str), data = regression_data, family = binomial())
      
      # Tabella formattata
      output$reg_results_table <- renderDT({
        datatable(format_regression_results(model),
                  options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE)
      })
      
      # Output grezzo
      output$reg_results_raw <- renderPrint({
        summary(model)
      })
    }, error = function(e) {
      output$reg_results_table <- renderDT(NULL)
      output$reg_results_raw <- renderPrint({
        cat("Si è verificato un errore durante l'esecuzione della regressione:\n")
        cat(as.character(e))
      })
    })
  })
}

# Run the app
shinyApp(ui, server)