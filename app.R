############################
#                          #
#     Kimberle Shen        #
#     November 2020        #
#  Immuneering challenge   #
#                          #
############################

library(shiny)
library(tidyverse)
library(plotly)
library(Biobase)
library(DT)
library(ggplot2)
library(shinythemes)

# get eSet data
eSet <- readRDS("generated_eSet.RDS")

# phenoData
ph <- phenoData(eSet)
ph <- as(ph, "data.frame")

# countsData
as <- assayData(eSet)
as_expr <- (as[["exprs"]])

# Display a first interactive table called "Patient metadata" which shows the sample metadata from the eSet.
# The user should be able to select specific subjects from this table based on filters like batch or day. 

# A second interactive table named "Expression data" displays the expression data filtered to show only selected subjects from the first table.
# The user should be able to select specific genes from this table. 
# The user can also apply filters to the data based on disease, batch or sample date. 

# A third section of the app shows an interactive boxplot of counts of the selected genes from table 2.
# The boxplot should also display the points as individual counts per subject, 
# and the ID of the subject should be displayed when hovering over the points. 

ui <- fluidPage(
    theme = shinytheme("cerulean"),
    h2("Patient metadata"),
    DT::dataTableOutput("metadataTable"),
    h2("Expression data"),
    DT::dataTableOutput("filteredTable"),
    h2("Box plot data"),
    plotlyOutput("boxplot")
)

server <- function(input, output, session) {
    
    # First table of patient metadata, with filter at top of table
    output$metadataTable = DT::renderDT({
        datatable(
            ph, 
            filter = "top"
        )
    })
    
    # Use selection from first table to create second table with counts
    metadataTable_selected <- reactive({
        patients <- input$metadataTable_rows_selected
        
        # get pheno data for selected patients
        one <- ph[sort(patients),]
        one <-tibble::rownames_to_column(one, "subject")
        
        # get expression data for selected patients 
        two <- as_expr[, sort(patients)]
        two <- t(two)
        two <- as.data.frame(two)
        two <-tibble::rownames_to_column(two, "subject")
        two <- as.data.frame(two)
        t <- reshape(two,idvar = "subject",
                     varying = list(names(two)[-1]),
                     direction = 'long',times = names(two)[-1],
                     timevar = 'gene',v.names = 'counts')
        rownames(t) <- NULL
        
        # merge pheno data and expression data
        m <- merge(one, t,by= "subject")
        m <- m[c("gene", "subject", "counts", "Disease", "Batch", "Sample_date")]
        m
    })
    
    # Second table of expression data, with filter at top
    output$filteredTable = DT::renderDT({
        validate(
            need(length(input$metadataTable_rows_selected) > 1, message = "Please select more rows from Patient metadata")
        )
        datatable(
            metadataTable_selected(), 
            rownames = FALSE,
            filter = "top"
        )
    })

    # Create boxplot
    output$boxplot <- renderPlotly({
        validate(
            need(length(input$filteredTable_rows_selected) > 0, message = "Please select more rows from Expression data")
        )
        p <- ggplot(data = metadataTable_selected()[sort(input$filteredTable_rows_selected), ], aes(x=gene, y=counts)) + geom_boxplot() + geom_point(aes(text=sprintf("subject: %s", subject)))
        ggplotly(p)
    })
}



# Run the application 
shinyApp(ui = ui, server = server)
