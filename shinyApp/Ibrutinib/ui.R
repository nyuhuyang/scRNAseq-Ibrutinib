source("util.R")
library(shinyhelper)
library(shinycssloaders)
sc1conf = readRDS("sc1conf.rds")
sc1def = readRDS("sc1def.rds")
sc1gene = readRDS("sc1gene.rds")


### Start server code
shinyUI(fluidPage(
### HTML formatting of error messages

tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))),
list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))),


### Page title
titlePanel("PBMC from Ibrutinib plus Nivolumab on a phase 1"),
navbarPage(
    NULL,
    ### Tab.a0: cellInfo on dimRed
    tabPanel(
        HTML("CellInfo"),
        h4("Cell information on reduced dimensions"),
        "In this tab, users can visualise cell information on a single or multiple low-dimensional representions.",
        br(),br(),
        fluidRow(
            column(
                3, h4("Dimension Reduction"),
                fluidRow(
                    column(
                        12, selectInput("sc1a0drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                        selected = sc1def$dimred[1]),
                        selectInput("sc1a0drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                selected = sc1def$dimred[2]))
                )
            ), # End of column (3 space)
            column(
                3, actionButton("sc1a0tog1", "Toggle to subset cells"),
                conditionalPanel(
                    condition = "input.sc1a0tog1 % 2 == 1",
                    selectInput("sc1a0sub1_1", "Cell information to subset:",
                                            choices = sc1conf[grp == TRUE]$UI,
                                            selected = sc1def$grp1),
                    uiOutput("sc1a0sub1.ui"),
                    #---------------------
                    actionButton("sc1a0tog2", "Toggle to further subset cells"),
                    conditionalPanel(
                        condition = "input.sc1a0tog2 % 2 == 1",
                        selectInput("sc1a0sub2_1", "Cell information to subset:",
                                                choices = sc1conf[grp == TRUE]$UI,
                                                selected = sc1def$grp2),
                        uiOutput("sc1a0sub2.ui"),
                        #---------------------
                        actionButton("sc1a0tog3", "Toggle to further subset cells"),
                        conditionalPanel(
                            condition = "input.sc1a0tog3 % 2 == 1",
                            selectInput("sc1a0sub3_1", "Cell information to subset:",
                                        choices = sc1conf[grp == TRUE]$UI,
                                        selected = sc1def$grp2),
                            uiOutput("sc1a0sub3.ui")
                        )
                    )
                )
            ), # End of column (3 space)
            column(
                6, actionButton("sc1a0tog4", "Toggle graphics controls"),
                conditionalPanel(
                    condition = "input.sc1a0tog4 % 2 == 1",
                    fluidRow(
                        column(
                            6, sliderInput("sc1a0siz", "Point size:",
                                                         min = 0, max = 4, value = 1.25, step = 0.25),
                            radioButtons("sc1a0psz", "Plot size:",
                                                     choices = c("Small", "Medium", "Large", "Extra Large"),
                                                     selected = "Large", inline = TRUE),
                            radioButtons("sc1a0fsz", "Font size:",
                                                     choices = c("Extra Small","Small", "Medium", "Large"),
                                                     selected = "Medium", inline = TRUE),
                            radioButtons("sc1a0asp", "Aspect ratio:",
                                                     choices = c("Square", "Fixed", "Free"),
                                                     selected = "Square", inline = TRUE),
                            checkboxInput("sc1a0txt", "Show axis text", value = FALSE)
                        ),
                        column(
                            6,
                            selectInput("sc1a0col1","Select a color palette",
                                                    choices=rownames(pal.info),
                                                    selected="default"),
                            radioButtons("sc1a0ord1", "Plot order:",
                                                     choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                     selected = "Original", inline = TRUE),
                            radioButtons("sc1a0lab1", "Labels:",
                                                     choices = c("No labels", "black text","black labels", "color labels"),
                                                     selected = "color labels", inline = FALSE),
                            checkboxInput("sc1a0leg", "Show Legend", value = FALSE),
                            radioButtons("sc1a0legpos", "Legend positions:",
                                         choices = c("top", "right", "bottom"),
                                         selected = "bottom", inline = TRUE)
                        )
                    )
                )
            ) # End of column (6 space)
        ),     # End of fluidRow (4 space)
        fluidRow(
            column(
                12, h4("Cell information"),
                fluidRow(
                    column(
                        4, selectInput("sc1a0inp1", "group by:",
                                                     choices = sc1conf$UI,
                                                     selected = sc1def$grp1) %>%
                            helper(type = "inline", size = "m", fade = TRUE,
                                         title = "Cell information to colour cells by",
                                         content = c("Select cell information to colour cells,
                                                                 equivalent to 'group.by' in Seurat DimPlot"))
                    ),
                    column(
                        4, selectInput("sc1a0inpsplt", "split by:",
                                                     choices = c("no split",sc1conf$ID[!is.na(sc1conf$fID)]),
                                                     selected = "no split") %>%
                            helper(type = "inline", size = "m", fade = TRUE,
                                         title = "Cell information to split cells by",
                                         content = c("Select cell information to split cells,
                                                                 equivalent to split.by in Seurat DimPlot",
                                                                 "- Default is no split"))
                    ),
                    column(
                        4,
                        radioButtons("sc1a0arrange", NULL,
                                     choices = c("auto","1row","1column"),
                                     selected = "auto", inline = TRUE)

                    ),
                ),
                fluidRow(column(12, withSpinner(uiOutput("sc1a0oup1.ui")))),
                downloadButton("sc1a0oup1.png", "Download png"),
                downloadButton("sc1a0oup1.jpeg", "Download jpeg"),
                downloadButton("meta_data.rds", "Download annotation and UMAP"),
                downloadButton("csr_gexpr.h5ad", "Download normalized gene expression(h5ad)"), br(),
                div(style="display:inline-block",
                        numericInput("sc1a0oup1.h", "png / jpeg height:", width = "138px",
                                                 min = 4, max = 50, value = 13, step = 0.5)),
                div(style="display:inline-block",
                        numericInput("sc1a0oup1.w", "png / jpeg width:", width = "138px",
                                                 min = 4, max = 50, value = 13, step = 0.5)), br(),
                actionButton("sc1a0tog9", "Toggle to show cell numbers "),
                conditionalPanel(
                    condition = "input.sc1a0tog9 % 2 == 1",
                    h4("Cell numbers"),
                    dataTableOutput("sc1a0.dt")
                )
            ), # End of column (12 space)
        )        # End of fluidRow (4 space)
    ),         # End of tab (2 space)
    ### Tab.a1: GeneExpr on dimRed
    tabPanel(
        HTML("GeneExpr"),
        h4("GeneExpr on reduced dimensions"),
        "In this tab, users can visualise Gene expression on a single or multiple low-dimensional representions.",
        br(),br(),
        fluidRow(
            column(
                3, h4("Dimension Reduction"),
                fluidRow(
                    column(
                        12, selectInput("sc1a1drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                        selected = sc1def$dimred[1]),
                        selectInput("sc1a1drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                    selected = sc1def$dimred[2]))
                )
            ), # End of column (3 space)
            column(
                3, actionButton("sc1a1tog1", "Toggle to subset cells"),
                conditionalPanel(
                    condition = "input.sc1a1tog1 % 2 == 1",
                    selectInput("sc1a1sub1_1", "Cell information to subset:",
                                choices = sc1conf[grp == TRUE]$UI,
                                selected = sc1def$grp1),
                    uiOutput("sc1a1sub1.ui"),
                    #---------------------
                    actionButton("sc1a1tog2", "Toggle to further subset cells"),
                    conditionalPanel(
                        condition = "input.sc1a1tog2 % 2 == 1",
                        selectInput("sc1a1sub2_1", "Cell information to subset:",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp2),
                        uiOutput("sc1a1sub2.ui"),
                        #---------------------
                        actionButton("sc1a1tog3", "Toggle to further subset cells"),
                        conditionalPanel(
                            condition = "input.sc1a1tog3 % 2 == 1",
                            selectInput("sc1a1sub3_1", "Cell information to subset:",
                                        choices = sc1conf[grp == TRUE]$UI,
                                        selected = sc1def$grp2),
                            uiOutput("sc1a1sub3.ui")
                        )
                    )
                )
            ), # End of column (3 space)
            column(
                6, actionButton("sc1a1tog4", "Toggle graphics controls"),
                conditionalPanel(
                    condition = "input.sc1a1tog4 % 2 == 1",
                    fluidRow(
                        column(
                            6, sliderInput("sc1a1siz", "Point size:",
                                           min = 0, max = 4, value = 1.25, step = 0.25),
                            checkboxInput("sc1a1max", "unified scale bar", value = TRUE) %>%
                                helper(type = "inline", size = "m", fade = TRUE,
                                       title = "Use global maximal scale",
                                       content = c("Select maximal expression level in the entire dataset.
                                                   This is useful for cross-group comparison")),
                            radioButtons("sc1a1psz", "Plot size:",
                                         choices = c("Small", "Medium", "Large", "Extra Large"),
                                         selected = "Large", inline = TRUE),
                            radioButtons("sc1a1fsz", "Font size:",
                                         choices = c("Extra Small","Small", "Medium", "Large"),
                                         selected = "Medium", inline = TRUE),
                            radioButtons("sc1a1asp1", "Aspect ratio:",
                                         choices = c("Square", "Fixed", "Free"),
                                         selected = "Fixed", inline = TRUE),
                            checkboxInput("sc1a1txt", "Show axis text", value = FALSE)
                        ),
                        column(
                            6,
                            selectInput("sc1a1col1","Select a color palette",
                                        choices=c(rownames(pal.info)[pal.info$category %in% c("seq","div")]),
                                        selected="default"),
                            uiOutput("sc1a1bg.ui"),
                            radioButtons("sc1a1ord1", "Plot order:",
                                         choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                         selected = "Max-1st", inline = TRUE),
                            checkboxInput("sc1a1leg", "Show Legend", value = TRUE),
                            radioButtons("sc1a1legpos", "Legend positions:",
                                         choices = c("top", "right", "bottom"),
                                         selected = "bottom", inline = TRUE),
                        )
                    )
                )
            ) # End of column (6 space)
        ),     # End of fluidRow (4 space)
        fluidRow(
            column(
                12, h4("Gene expression"),
                fluidRow(
                    column(
                        4, selectInput("sc1a1inp1", "Gene name:", choices=NULL) %>%
                            helper(type = "inline", size = "m", fade = TRUE,
                                   title = "Gene expression to colour cells by",
                                   content = c("Select gene to colour cells by gene expression",
                                               paste0("- Gene expression are coloured in a ",
                                                      "White-Red colour scheme which can be ",
                                                      "changed in the plot controls")))
                    ),
                    column(
                        4, selectInput("sc1a1inpsplt", "split by:",
                                       choices = c("no split",sc1conf$ID[!is.na(sc1conf$fID)]),
                                       selected = "no split") %>%
                            helper(type = "inline", size = "m", fade = TRUE,
                                   title = "Cell information to split cells by",
                                   content = c("Select cell information to split cells,
                                                                 equivalent to split.by in Seurat DimPlot",
                                               "- Default is no split"))
                    ),
                    column(
                        4,
                        radioButtons("sc1a1arrange", NULL,
                                     choices = c("auto","1row","1column"),
                                     selected = "auto", inline = TRUE)
                    ),
                ),
                fluidRow(column(12, withSpinner(uiOutput("sc1a1oup1.ui")))),
                fluidRow(
                    column(
                        4,
                        actionButton("sc1a1tog9", "Toggle to show cell numbers / statistics"),
                        conditionalPanel(
                            condition = "input.sc1a1tog9 % 2 == 1",
                        h4("All cells"),
                        dataTableOutput("sc1a1_all.dt"),
                        br(),
                        h4("Selected cells"),
                        dataTableOutput("sc1a1_selec.dt"),
                        style="border-bottom: 2px solid black"
                        )
                    ),    # End of column (1 space)
                    column(6,
                           downloadButton("sc1a1oup1.png", "Download png"),
                           downloadButton("sc1a1oup1.jpeg", "Download jpeg"),
                           div(style="display:inline-block",
                               numericInput("sc1a1oup1.h", "height:", width = "70px",
                                            min = 4, max = 20, value = 13, step = 0.5)),
                           div(style="display:inline-block",
                               numericInput("sc1a1oup1.w", "width:", width = "70px",
                                            min = 4, max = 20, value = 13, step = 0.5))
                    ) # End of column (1 space)
                ),  # End of fluidRow (3 space)
                fluidRow(column(12,
                actionButton("sc1a1tog5", "Toggle Density plot"),
                conditionalPanel(
                    condition = "input.sc1a1tog5 % 2 == 1",
                        selectInput("sc1a1grp", "group by:",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp2),
                    actionButton("sc1a1tog5", "Toggle Density plot graphics control"),
                    conditionalPanel(
                        condition = "input.sc1a1tog5 % 2 == 1",
                        fluidRow(
                            column(
                                6,
                                radioButtons("sc1a1asp2", "Aspect ratio:",
                                             choices = c("Square", "Fixed", "Free"),
                                             selected = "Free", inline = TRUE),
                                selectInput("sc1a1col2", "Colour for density plot",
                                            choices=c("default",rownames(pal.info)[pal.info$category %in% "qual"]),
                                            selected="default"),
                                radioButtons("sc1a1ord2", "Group order:",
                                             choices = c("As-it-is", "Reverse"),
                                             selected = "As-it-is", inline = TRUE),
                                radioButtons("sc1a1dtype", "Density plot type:",
                                             choices = c("density","density_ridges"),
                                             selected = "density_ridges", inline = TRUE)
                            ) # End of column (4 space)
                        ) # End of fluidRow (1 space)
                    ), # End of conditionalPanel (1 space)
                       withSpinner(uiOutput("sc1a1oup2.ui")),
                       downloadButton("sc1a1oup2.png", "Download png"),
                       downloadButton("sc1a1oup2.jpeg", "Download jpeg"),
                       div(style="display:inline-block",
                           numericInput("sc1a1oup1.h", "height:", width = "70px",
                                        min = 4, max = 20, value = 13, step = 0.5)),
                       div(style="display:inline-block",
                           numericInput("sc1a1oup1.w", "width:", width = "70px",
                                        min = 4, max = 20, value = 13, step = 0.5))
                           ) # End of conditionalPanel
                    ), # End of fluidRow (1 space)
                ) # End of conditionalPanel
            )  # End of column
        )# End of fluidRow ( space)
    ),         # End of tab (2 space)
 ### Tab.a2: cellInfo vs geneExpr on dimRed
    tabPanel(
        HTML("CellInfo vs GeneExpr"),
        h4("Cell information vs gene expression on reduced dimensions"),
        "In this tab, users can visualise both cell information and gene ",
        "expression side-by-side on low-dimensional representions.",
        br(),br(),
        fluidRow(
            column(
                3, h4("Dimension Reduction"),
                fluidRow(
                    column(
                        12, selectInput("sc1a2drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                     selected = sc1def$dimred[1]),
                        selectInput("sc1a2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                selected = sc1def$dimred[2]))
                )
            ), # End of column (6 space)
            column(
                3, actionButton("sc1a2tog1", "Toggle to subset cells"),
                conditionalPanel(
                    condition = "input.sc1a2tog1 % 2 == 1",
                    selectInput("sc1a2sub1_1", "Cell information to subset:",
                                            choices = sc1conf[grp == TRUE]$UI,
                                            selected = sc1def$grp1),
                    uiOutput("sc1a2sub1.ui"),
                    actionButton("sc1a2tog2", "Toggle to further subset cells"),
                    conditionalPanel(
                        condition = "input.sc1a2tog2 % 2 == 1",
                        selectInput("sc1a2sub2_1", "Cell information to subset:",
                                                choices = sc1conf[grp == TRUE]$UI,
                                                selected = sc1def$grp2),
                        uiOutput("sc1a2sub2.ui"),
                        actionButton("sc1a2tog3", "Toggle to further subset cells"),
                        conditionalPanel(
                            condition = "input.sc1a2tog3 % 2 == 1",
                            selectInput("sc1a2sub3_1", "Cell information to subset:",
                                        choices = sc1conf[grp == TRUE]$UI,
                                        selected = sc1def$grp2),
                            uiOutput("sc1a2sub3.ui")
                        )
                    )
                )
            ), # End of column (6 space)
            column(
                6, actionButton("sc1a2tog4", "Toggle graphics controls"),
                conditionalPanel(
                    condition = "input.sc1a2tog4 % 2 == 1",
                    fluidRow(
                        column(
                            6, sliderInput("sc1a2siz", "Point size:",
                                                         min = 0, max = 4, value = 1.25, step = 0.25),
                            radioButtons("sc1a2psz", "Plot size:",
                                                     choices = c("Small", "Medium", "Large"),
                                                     selected = "Medium", inline = TRUE),
                            radioButtons("sc1a2fsz", "Font size:",
                                                     choices = c("Small", "Medium", "Large"),
                                                     selected = "Medium", inline = TRUE)
                        ),
                        column(
                            6, radioButtons("sc1a2asp", "Aspect ratio:",
                                                            choices = c("Square", "Fixed", "Free"),
                                                            selected = "Square", inline = TRUE),
                            checkboxInput("sc1a2txt", "Show axis text", value = FALSE)
                        )
                    )
                )
            )    # End of column (6 space)
        ),     # End of fluidRow (3 space)
        fluidRow(
            column(
                6, style="border-right: 2px solid black", h4("Cell information"),
                fluidRow(
                    column(
                        6, selectInput("sc1a2inp1", "Cell information:",
                                                     choices = sc1conf$UI,
                                                     selected = sc1def$meta1) %>%
                            helper(type = "inline", size = "m", fade = TRUE,
                                         title = "Cell information to colour cells by",
                                         content = c("Select cell information to colour cells",
                                                                 "- Categorical covariates have a fixed colour palette",
                                                                 paste0("- Continuous covariates are coloured in a ",
                                                                                "Blue-Yellow-Red colour scheme, which can be ",
                                                                                "changed in the plot controls")))
                    ),# End of column (6 space)
                    column(
                        6, actionButton("sc1a2tog5", "Toggle plot controls"),
                        conditionalPanel(
                            condition = "input.sc1a2tog5 % 2 == 1",
                            selectInput("sc1a2col1", "Select a color palette",
                                        choices=rownames(pal.info),
                                        selected="default"),
                            radioButtons("sc1a2ord1", "Plot order:",
                                                     choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                     selected = "Original", inline = TRUE),
                            radioButtons("sc1a2lab1", "Labels:",
                                                     choices = c("No labels", "black text","black labels", "color labels"),
                                                     selected = "black labels", inline = FALSE),
                            checkboxInput("sc1a2leg1", "Show Legend", value = FALSE),
                            radioButtons("sc1a2legpos1", "Legend positions:",
                                         choices = c("top", "right", "bottom"),
                                         selected = "bottom", inline = TRUE),
                        )
                    )# End of column (6 space)
                ), # End of fluidRow (2 column)
                fluidRow(column(12, withSpinner(uiOutput("sc1a2oup1.ui")))),
                downloadButton("sc1a2oup1.png", "Download png"),
                downloadButton("sc1a2oup1.jpeg", "Download jpeg"), br(),
                div(style="display:inline-block",
                        numericInput("sc1a2oup1.h", "png / jpeg height:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5)),
                div(style="display:inline-block",
                        numericInput("sc1a2oup1.w", "png / jpeg width:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5)), br(),
                actionButton("sc1a2tog9", "Toggle to show cell numbers / statistics"),
                conditionalPanel(
                    condition = "input.sc1a2tog9 % 2 == 1",
                    h4("Cell numbers / statistics"),
                    radioButtons("sc1a2splt", "Split continuous cell info into:",
                                             choices = c("Quartile", "Decile"),
                                             selected = "Decile", inline = TRUE),
                    dataTableOutput("sc1a2.dt")
                )
            ), # End of column (6 space)
            column(
                6, h4("Gene expression"),
                fluidRow(
                    column(
                        6, selectInput("sc1a2inp2", "Gene name:", choices=NULL) %>%
                            helper(type = "inline", size = "m", fade = TRUE,
                                         title = "Gene expression to colour cells by",
                                         content = c("Select gene to colour cells by gene expression",
                                                                 paste0("- Gene expression are coloured in a ",
                                                                                "White-Red colour scheme which can be ",
                                                                                "changed in the plot controls")))
                    ),
                    column(
                        6, actionButton("sc1a2tog6", "Toggle plot controls"),
                        conditionalPanel(
                            condition = "input.sc1a2tog6 % 2 == 1",
                            selectInput("sc1a2col2", "Select a color spectrum",
                                        choices=c(rownames(pal.info)[pal.info$category %in% c("seq","div")]),
                                        selected= "default"),
                            checkboxInput("sc1a2max", "unified scale bar", value = TRUE) %>%
                                helper(type = "inline", size = "m", fade = TRUE,
                                       title = "Use global maximal scale",
                                       content = c("Select maximal expression level in the entire dataset.
                                                   This is useful for cross-group comparison")),
                            radioButtons("sc1a2ord2", "Plot order:",
                                                     choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                     selected = "Max-1st", inline = TRUE),
                            checkboxInput("sc1a2leg2", "Show Legend", value = FALSE),
                            radioButtons("sc1a2legpos2", "Legend positions:",
                                         choices = c("top", "right", "bottom"),
                                         selected = "bottom", inline = TRUE),
                        )
                    )
                ) ,
                fluidRow(column(12, withSpinner(uiOutput("sc1a2oup2.ui")))),
                downloadButton("sc1a2oup2.png", "Download png"),
                downloadButton("sc1a2oup2.jpeg", "Download jpeg"), br(),
                div(style="display:inline-block",
                        numericInput("sc1a2oup2.h", "png / jpeg height:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5)),
                div(style="display:inline-block",
                        numericInput("sc1a2oup2.w", "png / jpeg width:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5))
            )    # End of column (6 space)
        )        # End of fluidRow (4 space)
    ),         # End of tab (2 space)

    ### Tab.a3: cellInfo vs cellInfo on dimRed
    tabPanel(
        HTML("CellInfo vs CellInfo"),
        h4("Cell information vs cell information on dimension reduction"),
        "In this tab, users can visualise two cell informations side-by-side",
        "on low-dimensional representions.",
        br(),br(),
        fluidRow(
            column(
                3, h4("Dimension Reduction"),
                fluidRow(
                    column(
                        12, selectInput("sc1a3drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                     selected = sc1def$dimred[1]),
                        selectInput("sc1a3drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                selected = sc1def$dimred[2]))
                )
            ), # End of column (6 space)
            column(
                3, actionButton("sc1a3tog1", "Toggle to subset cells"),
                conditionalPanel(
                    condition = "input.sc1a3tog1 % 2 == 1",
                    selectInput("sc1a3sub1_1", "Cell information to subset:",
                                            choices = sc1conf[grp == TRUE]$UI,
                                            selected = sc1def$grp1),
                    uiOutput("sc1a3sub1.ui"),
                    #---------------
                    actionButton("sc1a3tog2", "Toggle to further subset cells"),
                    conditionalPanel(
                        condition = "input.sc1a3tog2 % 2 == 1",
                        selectInput("sc1a3sub2_1", "Cell information to subset:",
                                                choices = sc1conf[grp == TRUE]$UI,
                                                selected = sc1def$grp2),
                        uiOutput("sc1a3sub2.ui"),
                        #---------------
                        actionButton("sc1a3tog3", "Toggle to further subset cells"),
                        conditionalPanel(
                            condition = "input.sc1a3tog3 % 2 == 1",
                            selectInput("sc1a3sub3_1", "Cell information to subset:",
                                        choices = sc1conf[grp == TRUE]$UI,
                                        selected = sc1def$grp2),
                            uiOutput("sc1a3sub3.ui")
                        )
                    )
                )
            ), # End of column (6 space)
            column(
                6, actionButton("sc1a3tog4", "Toggle graphics controls"),
                conditionalPanel(
                    condition = "input.sc1a3tog4 % 2 == 1",
                    fluidRow(
                        column(
                            6, sliderInput("sc1a3siz", "Point size:",
                                                         min = 0, max = 4, value = 1.25, step = 0.25),
                            radioButtons("sc1a3psz", "Plot size:",
                                                     choices = c("Small", "Medium", "Large"),
                                                     selected = "Medium", inline = TRUE),
                            radioButtons("sc1a3fsz", "Font size:",
                                                     choices = c("Small", "Medium", "Large"),
                                                     selected = "Medium", inline = TRUE)
                        ),
                        column(
                            6, radioButtons("sc1a3asp", "Aspect ratio:",
                                                            choices = c("Square", "Fixed", "Free"),
                                                            selected = "Square", inline = TRUE),
                            checkboxInput("sc1a3txt", "Show axis text", value = FALSE)
                        )
                    )
                )
            )    # End of column (6 space)
        ),     # End of fluidRow (4 space)
        fluidRow(
            column(
                6, style="border-right: 2px solid black", h4("Cell information 1"),
                fluidRow(
                    column(
                        6, selectInput("sc1a3inp1", "Cell information:",
                                                     choices = sc1conf$UI,
                                                     selected = sc1def$meta1) %>%
                            helper(type = "inline", size = "m", fade = TRUE,
                                         title = "Cell information to colour cells by",
                                         content = c("Select cell information to colour cells",
                                                                 "- Categorical covariates have a fixed colour palette",
                                                                 paste0("- Continuous covariates are coloured in a ",
                                                                                "Blue-Yellow-Red colour scheme, which can be ",
                                                                                "changed in the plot controls")))
                    ),
                    column(
                        6, actionButton("sc1a3tog5", "Toggle plot controls"),
                        conditionalPanel(
                            condition = "input.sc1a3tog5 % 2 == 1",
                            selectInput("sc1a3col1", "Select a color palette",
                                        choices=rownames(pal.info),
                                        selected="default"),
                            radioButtons("sc1a3ord1", "Plot order:",
                                                     choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                     selected = "Original", inline = TRUE),
                            radioButtons("sc1a3lab1", "Labels:",
                                                     choices = c("No labels", "black text","black labels", "color labels"),
                                                     selected = "black labels", inline = FALSE)
                        )
                    )
                ),
                fluidRow(column(12, withSpinner(uiOutput("sc1a3oup1.ui")))),
                downloadButton("sc1a3oup1.png", "Download png"),
                downloadButton("sc1a3oup1.jpeg", "Download jpeg"), br(),
                div(style="display:inline-block",
                        numericInput("sc1a3oup1.h", "png / jpeg height:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5)),
                div(style="display:inline-block",
                        numericInput("sc1a3oup1.w", "png / jpeg width:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5))
            ), # End of column (6 space)
            column(
                6, h4("Cell information 2"),
                fluidRow(
                    column(
                        6, selectInput("sc1a3inp2", "Cell information:",
                                                     choices = sc1conf$UI,
                                                     selected = sc1def$meta2) %>%
                            helper(type = "inline", size = "m", fade = TRUE,
                                         title = "Cell information to colour cells by",
                                         content = c("Select cell information to colour cells",
                                                                 "- Categorical covariates have a fixed colour palette",
                                                                 paste0("- Continuous covariates are coloured in a ",
                                                                                "Blue-Yellow-Red colour scheme, which can be ",
                                                                                "changed in the plot controls")))
                    ),
                    column(
                        6, actionButton("sc1a3tog6", "Toggle plot controls"),
                        conditionalPanel(
                            condition = "input.sc1a3tog6 % 2 == 1",
                            selectInput("sc1a3col2", "Select a color palette",
                                        choices=rownames(pal.info),
                                        selected="default"),
                            radioButtons("sc1a3ord2", "Plot order:",
                                                     choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                     selected = "Original", inline = TRUE),
                            radioButtons("sc1a3lab2", "Labels:",
                                                     choices = c("No labels", "black text","black labels", "color labels"),
                                                     selected = "black labels", inline = FALSE)
                        )
                    )
                ),
                fluidRow(column(12, withSpinner(uiOutput("sc1a3oup2.ui")))),
                downloadButton("sc1a3oup2.png", "Download png"),
                downloadButton("sc1a3oup2.jpeg", "Download jpeg"), br(),
                div(style="display:inline-block",
                        numericInput("sc1a3oup2.h", "png / jpeg height:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5)),
                div(style="display:inline-block",
                        numericInput("sc1a3oup2.w", "png / jpeg width:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5))
            )    # End of column (6 space)
        )        # End of fluidRow (4 space)
    ),         # End of tab (2 space)

    ### Tab.a4: geneExpr vs geneExpr on dimRed
    tabPanel(
        HTML("GeneExpr vs GeneExpr"),
        h4("Gene expression vs gene expression on dimension reduction"),
        "In this tab, users can visualise two gene expressions side-by-side ",
        "on low-dimensional representions.",
        br(),br(),
        fluidRow(
            column(
                3, h4("Dimension Reduction"),
                fluidRow(
                    column(
                        12, selectInput("sc1a4drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                     selected = sc1def$dimred[1]),
                        selectInput("sc1a4drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                selected = sc1def$dimred[2]))
                )
            ), # End of column (6 space)
            column(
                3, actionButton("sc1a4tog1", "Toggle to subset cells"),
                conditionalPanel(
                    condition = "input.sc1a4tog1 % 2 == 1",
                    selectInput("sc1a4sub1_1", "Cell information to subset:",
                                            choices = sc1conf[grp == TRUE]$UI,
                                            selected = sc1def$grp1),
                    uiOutput("sc1a4sub1.ui"),
                    #--------------------
                    actionButton("sc1a4tog2", "Toggle to subset cells"),
                    conditionalPanel(
                        condition = "input.sc1a4tog2 % 2 == 1",
                        selectInput("sc1a4sub2_1", "Cell information to subset:",
                                                choices = sc1conf[grp == TRUE]$UI,
                                                selected = sc1def$grp2),
                        uiOutput("sc1a4sub2.ui"),
                        #--------------------
                        actionButton("sc1a4tog3", "Toggle to subset cells"),
                        conditionalPanel(
                            condition = "input.sc1a4tog3 % 2 == 1",
                            selectInput("sc1a4sub3_1", "Cell information to subset:",
                                        choices = sc1conf[grp == TRUE]$UI,
                                        selected = sc1def$grp2),
                            uiOutput("sc1a4sub3.ui")
                        )
                    )
                )
            ), # End of column (6 space)
            column(
                6, actionButton("sc1a4tog4", "Toggle graphics controls"),
                conditionalPanel(
                    condition = "input.sc1a4tog4 % 2 == 1",
                    fluidRow(
                        column(
                            6, sliderInput("sc1a4siz", "Point size:",
                                                         min = 0, max = 4, value = 1.25, step = 0.25),
                            checkboxInput("sc1a4max", "unified scale bar", value = TRUE) %>%
                                helper(type = "inline", size = "m", fade = TRUE,
                                       title = "Use global maximal scale",
                                       content = c("Select maximal expression level in the entire dataset.
                                                   This is useful for cross-group comparison")),
                            radioButtons("sc1a4psz", "Plot size:",
                                                     choices = c("Small", "Medium", "Large"),
                                                     selected = "Medium", inline = TRUE),
                            radioButtons("sc1a4fsz", "Font size:",
                                                     choices = c("Small", "Medium", "Large"),
                                                     selected = "Medium", inline = TRUE)
                        ),
                        column(
                            6, radioButtons("sc1a4asp", "Aspect ratio:",
                                                            choices = c("Square", "Fixed", "Free"),
                                                            selected = "Square", inline = TRUE),
                            checkboxInput("sc1a4txt", "Show axis text", value = FALSE),
                            checkboxInput("sc1a4leg", "Show Legend", value = TRUE),
                            radioButtons("sc1a4legpos", "Legend positions:",
                                         choices = c("top", "right", "bottom"),
                                         selected = "bottom", inline = TRUE),
                        )
                    )
                )
            )    # End of column (6 space)
        ),     # End of fluidRow (4 space)
        fluidRow(
            column(
                6, style="border-right: 2px solid black", h4("Gene expression 1"),
                fluidRow(
                    column(
                        6, selectInput("sc1a4inp1", "Gene name:", choices=NULL) %>%
                            helper(type = "inline", size = "m", fade = TRUE,
                                         title = "Gene expression to colour cells by",
                                         content = c("Select gene to colour cells by gene expression",
                                                                 paste0("- Gene expression are coloured in a ",
                                                                                "White-Red colour scheme which can be ",
                                                                                "changed in the plot controls")))
                    ),
                    column(
                        6, actionButton("sc1a4tog5", "Toggle plot controls"),
                        conditionalPanel(
                            condition = "input.sc1a4tog5 % 2 == 1",
                            selectInput("sc1a4col1", "Select a color spectrum",
                                        choices=c(rownames(pal.info)[pal.info$category %in% c("seq","div")]),
                                        selected = "default"),
                            radioButtons("sc1a4ord1", "Plot order:",
                                         choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                         selected = "Max-1st", inline = TRUE)
                        )
                    )
                ),
                fluidRow(column(12, withSpinner(uiOutput("sc1a4oup1.ui")))),
                downloadButton("sc1a4oup1.png", "Download png"),
                downloadButton("sc1a4oup1.jpeg", "Download jpeg"), br(),
                div(style="display:inline-block",
                        numericInput("sc1a4oup1.h", "png / jpeg height:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5)),
                div(style="display:inline-block",
                        numericInput("sc1a4oup1.w", "png / jpeg width:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5))
            ), # End of column (6 space))
            column(
                6, h4("Gene expression 2"),
                fluidRow(
                    column(
                        6, selectInput("sc1a4inp2", "Gene name:", choices=NULL) %>%
                            helper(type = "inline", size = "m", fade = TRUE,
                                         title = "Gene expression to colour cells by",
                                         content = c("Select gene to colour cells by gene expression",
                                                                 paste0("- Gene expression are coloured in a ",
                                                                                "White-Red colour scheme which can be ",
                                                                                "changed in the plot controls")))
                    ),
                    column(
                        6, actionButton("sc1a4tog6", "Toggle plot controls"),
                        conditionalPanel(
                            condition = "input.sc1a4tog6 % 2 == 1",
                            selectInput("sc1a4col2", "Select a color spectrum",
                                        choices=c(rownames(pal.info)[pal.info$category %in% c("seq","div")]),
                                        selected = "default"),
                            radioButtons("sc1a4ord2", "Plot order:",
                                         choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                         selected = "Max-1st", inline = TRUE)
                        )
                    )
                ),
                fluidRow(column(12, withSpinner(uiOutput("sc1a4oup2.ui")))),
                downloadButton("sc1a4oup2.png", "Download png"),
                downloadButton("sc1a4oup2.jpeg", "Download jpeg"), br(),
                div(style="display:inline-block",
                        numericInput("sc1a4oup2.h", "png / jpeg height:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5)),
                div(style="display:inline-block",
                        numericInput("sc1a4oup2.w", "png / jpeg width:", width = "138px",
                                                 min = 4, max = 50, value = 8, step = 0.5))
            )    # End of column (6 space)
        )        # End of fluidRow (4 space)
    ),         # End of tab (2 space)

 ### Tab.b1: Gene coexpression plot
 tabPanel(
     HTML("Gene coexpression"),
     h4("Coexpression of two genes on reduced dimensions"),
     "In this tab, users can visualise the coexpression of two genes ",
     "on low-dimensional representions.",
     br(),br(),
     fluidRow(
         column(
             3, h4("Gene Expression on Dimension Reduction"),
             fluidRow(
                 column(
                     12, selectInput("sc1b1drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                     selected = sc1def$dimred[1]),
                     selectInput("sc1b1drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                             selected = sc1def$dimred[2]),
                 selectInput("sc1b1inp1", "Gene 1:", choices=NULL) %>%
                     helper(type = "inline", size = "m", fade = TRUE,
                            title = "Gene expression to colour cells by",
                            content = c("Select gene to colour cells by gene expression",
                                        paste0("- Gene expression are coloured in a ",
                                               "White-Red colour scheme which can be ",
                                               "changed in the plot controls"))),
                 selectInput("sc1b1inp2", "Gene 2:", choices=NULL) %>%
                     helper(type = "inline", size = "m", fade = TRUE,
                            title = "Gene expression to colour cells by",
                            content = c("Select gene to colour cells by gene expression",
                                        paste0("- Gene expression are coloured in a ",
                                               "White-Blue colour scheme which can be ",
                                               "changed in the plot controls")))
                 )
             )
         ), # End of column (6 space)
         column(
             3, actionButton("sc1b1tog0", "Toggle to subset cells"),
             conditionalPanel(
                 condition = "input.sc1b1tog0 % 2 == 1",
                 selectInput("sc1b1sub1_1", "Cell information to subset:",
                                         choices = sc1conf[grp == TRUE]$UI,
                                        selected = sc1def$grp1),
                 uiOutput("sc1b1sub1.ui"),
                 actionButton("sc1b1tog1", "Toggle to further subset cells"),
                 conditionalPanel(
                     condition = "input.sc1b1tog1 % 2 == 1",
                     selectInput("sc1b1sub2_1", "Cell information to subset:",
                                             choices = sc1conf[grp == TRUE]$UI,
                                             selected = sc1def$grp2),
                     uiOutput("sc1b1sub2.ui")
                 )
             )
         ), # End of column (6 space)
         column(
             6, actionButton("sc1b1tog2", "Toggle graphics controls"),
             conditionalPanel(
                 condition = "input.sc1b1tog2 % 2 == 1",
                 fluidRow(
                     column(
                         6, sliderInput("sc1b1siz", "Point size:",
                                                        min = 0, max = 4, value = 1.5, step = 0.25),
                         checkboxInput("sc1b1max", "unified scale bar", value = TRUE) %>%
                             helper(type = "inline", size = "m", fade = TRUE,
                                    title = "Use global maximal scale",
                                    content = c("Select maximal expression level in the entire dataset.
                                                   This is useful for cross-group comparison")),
                         radioButtons("sc1b1psz", "Plot size:",
                                                    choices = c("Small", "Medium", "Large", "Extra Large"),
                                                    selected = "Medium", inline = TRUE),
                         radioButtons("sc1b1fsz", "Font size:",
                                                    choices = c("Extra Small","Small", "Medium", "Large"),
                                                    selected = "Medium", inline = TRUE)
                     ),
                     column(
                         6,
                         radioButtons("sc1b1mulcol","Colour for dimension reduction plot",
                                      choices = c("Bi-Colors", "Tri-Colors"),
                                      selected = "Tri-Colors"),
                         radioButtons("sc1b1col1", "",
                                      choices = c("Red (Gene1); Blue (Gene2); Purple (Both)",
                                                  "Green (Gene1); Blue (Gene2); Red (Both)",
                                                  "Orange (Gene1); Blue (Gene2); Red (Both)",
                                                  "Yellow (Gene1); Blue (Gene2); Red (Both)",
                                                  "Purple (Gene1); Green (Gene2); Red (Both)"),
                                      selected = "Orange (Gene1); Blue (Gene2); Red (Both)"),
                         radioButtons("sc1b1ord1", "Plot order:",
                                      choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                      selected = "Max-1st", inline = TRUE),
                         radioButtons("sc1b1asp", "Aspect ratio:",
                                      choices = c("Square", "Fixed", "Free"),
                                      selected = "Square", inline = TRUE),
                         checkboxInput("sc1b1txt", "Show axis text", value = FALSE)
                     )
                 )
             )
         )    # End of column (6 space)
     ),     # End of fluidRow (4 space)
     fluidRow(
         column(
             8,style="border-right: 2px solid black",
             withSpinner(uiOutput("sc1b1oup1.ui")),
             downloadButton("sc1b1oup1.png", "Download png"),
             downloadButton("sc1b1oup1.jpeg", "Download jpeg"),
             div(style="display:inline-block",
                 numericInput("sc1b1oup1.h", "height:", width = "70px",
                              min = 4, max = 20, value = 10, step = 0.5)),
             div(style="display:inline-block",
                 numericInput("sc1b1oup1.w", "width:", width = "70px",
                              min = 4, max = 20, value = 10, step = 0.5)),
             actionButton("sc1b1tog5", "Toggle Scatter plot"),
             conditionalPanel(
                 condition = "input.sc1b1tog5 % 2 == 1",
                 selectInput("sc1b1grp", "group by:",
                             choices = sc1conf[grp == TRUE]$UI,
                             selected = sc1def$grp2),
                 actionButton("sc1b1tog6", "Toggle Scatter plot graphics controls"),
                 conditionalPanel(
                     condition = "input.sc1b1tog6 % 2 == 1",
                     fluidRow(
                         column(
                             6,
                             selectInput("sc1b1col2", "Colour for scatter plot",
                                         choices=c("default",rownames(pal.info)[pal.info$category %in% "qual"]),
                                         selected="default"),
                             radioButtons("sc1b1ord2", "Group order:",
                                          choices = c("As-it-is", "Reverse"),
                                          selected = "As-it-is", inline = TRUE),
                             radioButtons("sc1b1dtype", "Density plot type:",
                                          choices = c("density","density_ridges"),
                                          selected = "density_ridges", inline = TRUE)
                         ) # End of column (4 space)
                     )
                 ), # End of conditionalPanel
                 withSpinner(uiOutput("sc1b1oup3.ui")),
                 downloadButton("sc1b1oup3.png", "Download png"),
                 downloadButton("sc1b1oup3.jpeg", "Download jpeg"),
                 div(style="display:inline-block",
                     numericInput("sc1b1oup3.h", "height:", width = "70px",
                                  min = 4, max = 20, value = 10, step = 0.5)),
                 div(style="display:inline-block",
                     numericInput("sc1b1oup3.w", "width:", width = "70px",
                                  min = 4, max = 20, value = 10, step = 0.5))
                 ) # End of conditionalPanel
         ), # End of column (9 space)
         column(
             4, uiOutput("sc1b1oup2.ui"),
             downloadButton("sc1b1oup2.png", "Download png"),
             downloadButton("sc1b1oup2.jpeg", "Download jpeg"),
             br(),
             style="border-bottom: 2px solid black",
             actionButton("sc1b1tog3", "All cells numbers / statistics"),
             conditionalPanel(
                 condition = "input.sc1b1tog3 % 2 == 1",
                 h4("All cells"),
                 dataTableOutput("sc1b1_all.dt"),
                 dataTableOutput("sc1b1_all_cor.dt")
             ),
             br(),

             style="border-top: 2px solid black",
             actionButton("sc1b1tog4", "Selected cells numbers / statistics"),
             conditionalPanel(
                 condition = "input.sc1b1tog4 % 2 == 1",
                 h4("Selected cells"),
                 dataTableOutput("sc1b1_selec.dt"),
                 dataTableOutput("sc1b1_selec_cor.dt")
             ) # End of conditionalPanel
         )    # End of column (3 space)
     ),        # End of fluidRow (4 space)
 ),         # End of tab (2 space)

 ### Tab.c1: violinplot / boxplot
 tabPanel(
        HTML("Violinplot"),
     h4("Cell information / gene expression violin plot / box plot"),
     "In this tab, users can visualise the gene expression or continuous cell information ",
     "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).",
     br(),br(),
     fluidRow(
         column(
             3, style="border-right: 2px solid black",
             selectInput("sc1c1inp1", "Cell Info (X-axis):",
                                     choices = sc1conf[grp == TRUE]$UI,
                                     selected = sc1def$grp1) %>%
                 helper(type = "inline", size = "m", fade = TRUE,
                                title = "Cell Info to group cells by",
                                content = c("Select categorical cell information to group cells by",
                                                        "- Single cells are grouped by this categorical covariate",
                                                        "- Plotted as the X-axis of the violin plot / box plot")),
             selectInput("sc1c1inp2", "Cell Info / Gene name (Y-axis):",
                         choices = sc1conf[grp == TRUE]$UI,
                         selected = sc1def$grp1) %>%
                 helper(type = "inline", size = "m", fade = TRUE,
                                title = "Cell Info / Gene to plot",
                                content = c("Select cell info / gene to plot on Y-axis",
                                                        "- Can be continuous cell information (e.g. nUMIs / scores)",
                                                        "- Can also be gene expression")),
             selectInput("sc1c1inp3", "Cell Info to split:",
                         choices = sc1conf[grp == TRUE]$UI,
                         selected = sc1def$grp2),
             radioButtons("sc1c1split", NULL, choices =  c("stack","split"),selected = "stack", inline = TRUE),
             radioButtons("sc1c1typ", "Plot type:",
                                        choices = c("violin", "boxplot","barplot"),
                                        selected = "violin", inline = FALSE),
             checkboxInput("sc1c1sig", "run wilcox test", value = FALSE),
             checkboxInput("sc1c1pts", "Show data points", value = FALSE),
             actionButton("sc1c1tog1", "Toggle to subset cells"),
             conditionalPanel(
                 condition = "input.sc1c1tog1 % 2 == 1",
                 selectInput("sc1c1sub1_1", "Cell information to subset:",
                                         choices = sc1conf[grp == TRUE]$UI,
                                         selected = sc1def$grp1),
                 uiOutput("sc1c1sub1.ui"),
                 #---------------------
                 actionButton("sc1c1tog2", "Toggle to further subset cells"),
                 conditionalPanel(
                     condition = "input.sc1c1tog2 % 2 == 1",
                     selectInput("sc1c1sub2_1", "Cell information to subset:",
                                             choices = sc1conf[grp == TRUE]$UI,
                                             selected = sc1def$grp2),
                     uiOutput("sc1c1sub2.ui"),
                     #----------------------
                     actionButton("sc1c1tog3", "Toggle to further subset cells"),
                     conditionalPanel(
                         condition = "input.sc1c1tog3 % 2 == 1",
                         selectInput("sc1c1sub3_1", "Cell information to subset:",
                                     choices = sc1conf[grp == TRUE]$UI,
                                     selected = sc1def$grp2),
                         uiOutput("sc1c1sub3.ui")
                     )
                 )
             ),
             actionButton("sc1c1tog4", "Toggle graphics controls"),
             conditionalPanel(
                 condition = "input.sc1c1tog4 % 2 == 1",
                 selectInput("sc1c1cols","Select a color palette:",
                             choices=c("default","white",rownames(pal.info)[pal.info$category %in% "qual"]),
                             selected="default"),
                 sliderInput("sc1c1siz", "Data point size",
                             min = 0, max = 8, value = 1, step = 0.25),
                 checkboxInput("sc1c1err", "Show error bar", value = TRUE),
                 selectInput("sc1c1scale","Select violin plot scale:",
                             choices=c("area","count","width"),
                             selected="area"),
                 radioButtons("sc1c1psz", "Plot size:",
                                            choices = c("Small", "Medium", "Large","Extra Large"),
                                            selected = "Medium", inline = TRUE),
                 radioButtons("sc1c1fsz", "Font size:",
                                            choices = c("Extra Small","Small", "Medium", "Large"),
                                            selected = "Medium", inline = TRUE),
                 radioButtons("sc1c1frt", "Rotate x axis label:",
                                            choices = c(0,30,45,90),
                                            selected = 45, inline = TRUE),
                 checkboxInput("sc1c1leg", "Show Legend", value = TRUE),
                 radioButtons("sc1c1legpos", "Legend positions:",
                              choices = c("top", "right", "bottom"),
                              selected = "bottom", inline = TRUE),
                 ),
             br(),
             actionButton("sc1c1tog5", "Toggle p value controls"),
             conditionalPanel(
                 condition = "input.sc1c1tog5 % 2 == 1",
                 sliderInput("sc1c1pcut", "log10 p value cut off",
                             min = 0, max = 300, value = 2, step = 1),
                 sliderInput("sc1c1pvalpos", "Change p value position",
                             min = -1, max = 3, value = 0, step = 0.01),
                 sliderInput("sc1c1ylim", "Change Y upper limit %",
                             min = -1, max = 4, value = 0, step = 0.01),
                 radioButtons("sc1c1plab", "Show p-value label type:",
                              selected = "p",inline = FALSE,
                              choices = c("p","p.signif","p.adj","p.adj.signif")),

             )
         ), # End of column (6 space)
         column(9, withSpinner(uiOutput("sc1c1oup.ui")),
                        downloadButton("sc1c1oup.png", "Download png"),
                        downloadButton("sc1c1oup.jpeg", "Download jpeg"), br(),
                        div(style="display:inline-block",
                                numericInput("sc1c1oup.h", "png / jpeg height:", width = "138px",
                                                         min = 4, max = 50, value = 10, step = 0.5)),
                        div(style="display:inline-block",
                                numericInput("sc1c1oup.w", "png / jpeg width:", width = "138px",
                                                         min = 4, max = 50, value = 10, step = 0.5)),br(),
                actionButton("sc1c1tog6", "Toggle to data summary"),
                conditionalPanel(
                    condition = "input.sc1c1tog6 % 2 == 1",
                    h4("data summary"),
                    dataTableOutput("sc1c1.dt")
                ),
                downloadButton("sc1c1oup1.xlsx", "Download xlsx")
         )    # End of column (6 space)
     )        # End of fluidRow (4 space)
 ),         # End of tab (2 space)

### Tab1.c2: Proportion plot
tabPanel(
    HTML("Proportion plot"),
    h4("Proportion / cell numbers across different cell information"),
    "In this tab, users can visualise the composition of single cells based on one discrete ",
    "cell information across another discrete cell information. ",
    "Usage examples include the library or cellcycle composition across clusters.",
    br(),br(),
    fluidRow(
        column(
            3, style="border-right: 2px solid black",
            selectInput("sc1c2inp1", "Cell Info (X-axis):",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp1) %>%
                helper(type = "inline", size = "m", fade = TRUE,
                             title = "Cell information to plot cells by",
                             content = c("Select categorical cell information to plot cells by",
                                                     "- Plotted as the X-axis of the proportion plot")),
            selectInput("sc1c2norm", "Cell Number / Proportion (Y-axis):",
                        choices = c("Cell Number", "Cell % by sample","Cell % by X-axis","Cell % by Y-axis"),
                        selected = "Cell Number"),
            selectInput("sc1c2inp2", "Cell Info to stack or split:",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp2),
            radioButtons("sc1c2split", NULL, choices =  c("stack","split"),selected = "stack", inline = TRUE),
            radioButtons("sc1c2typ", "Plot type:",
                                     choices = c("barplot","boxplot","piechart","pairedplot","scatterplot"),
                                     selected = "barplot", inline = FALSE),
            radioButtons("sc1c2sig", "Run statistical test", choices =  c("no test","unpaired Wilcoxon","paired Wilcoxon",
                                                                          "t test","chisq"),
                         inline = TRUE),
            checkboxInput("sc1c2pts", "Show data points", value = FALSE),
            actionButton("sc1c2tog1", "Toggle to subset cells"),
            conditionalPanel(
                condition = "input.sc1c2tog1 % 2 == 1",
                selectInput("sc1c2sub1_1", "Cell information to subset:",
                                        choices = sc1conf[grp == TRUE]$UI,
                                        selected = sc1def$grp1),
                uiOutput("sc1c2sub1.ui"),
                actionButton("sc1c2tog2", "Toggle to further subset cells"),
                conditionalPanel(
                    condition = "input.sc1c2tog2 % 2 == 1",
                    selectInput("sc1c2sub2_1", "Cell information to subset:",
                                            choices = sc1conf[grp == TRUE]$UI,
                                            selected = sc1def$grp2),
                    uiOutput("sc1c2sub2.ui"),
                    actionButton("sc1c2tog3", "Toggle to further subset cells"),
                    conditionalPanel(
                        condition = "input.sc1c2tog3 % 2 == 1",
                        selectInput("sc1c2sub3_1", "Cell information to subset:",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp2),
                        uiOutput("sc1c2sub3.ui")
                    )
                )
            ),
            br(),
            actionButton("sc1c2tog4", "Toggle graphics controls"),
            conditionalPanel(
                condition = "input.sc1c2tog4 % 2 == 1",
                selectInput("sc1c2cols","Select a color palette:",
                            choices=c("default",rownames(pal.info)[pal.info$category %in% "qual"]),
                            selected="default"),
                sliderInput("sc1c2yrange", "Change y-axis upper limit",
                            min = -1, max = 2, value = 0, step = 0.01),
                sliderInput("sc1c2siz", "Data point size",
                            min = 0, max = 8, value = 2, step = 0.25),
                sliderInput("sc1c2filt", "Remove low percentage label %:",
                            min = 0, max = 10, value = 1, step = 0.01),
                radioButtons("sc1c2psz", "Plot size:",
                                         choices = c("Small", "Medium", "Large", "Extra Large"),
                                         selected = "Medium", inline = TRUE),
                checkboxInput("sc1c2addline", "Add line", value = TRUE),
                checkboxInput("sc1c2err", "Show error bar", value = TRUE),
                radioButtons("sc1c2fsz", "Font size:",
                                         choices = c("Extra Small","Small", "Medium", "Large"),
                                         selected = "Medium", inline = TRUE),
                radioButtons("sc1c2lsz", "line size:",
                                         choices = c("Extra Small","Small", "Medium", "Large"),
                                         selected = "Medium", inline = TRUE),
                radioButtons("sc1c2lab", "Labels:",
                             choices = c("No %","No labels", "black text","color labels"),
                             selected = "color labels", inline = FALSE),
                checkboxInput("sc1c2flpxy", "Flip X/Y", value = FALSE),
                checkboxInput("sc1c2flpx", "Flip X axis order", value = FALSE),
                radioButtons("sc1c2frt", "Rotate x axis label:",
                                         choices = c(0,30,45,90),
                                         selected = 90, inline = TRUE),
                checkboxInput("sc1c2leg", "Show Legend", value = TRUE),
                radioButtons("sc1c2legpos", "Legend positions:",
                             choices = c("top", "right", "bottom"),
                             selected = "right", inline = TRUE),
                ),
            br(),
            actionButton("sc1c2tog5", "Toggle p value controls"),
            conditionalPanel(
                condition = "input.sc1c2tog5 % 2 == 1",
                sliderInput("sc1c2pcut", "log10 p value cut off",
                            min = 0, max = 300, value = 0, step = 1),
                sliderInput("sc1c2pvalpos", "Increase p value position",
                            min = -1, max = 2, value = 0, step = 0.01),
                radioButtons("sc1c2plab", "Show p-value label type:",
                             selected = "p",inline = FALSE,
                             choices = c("p","p.signif","p.adj","p.adj.signif")),

            )

        ), # End of column (6 space)
        column(9, withSpinner(uiOutput("sc1c2oup.ui")),
                     downloadButton("sc1c2oup.png", "Download png"),
                     downloadButton("sc1c2oup.jpeg", "Download jpeg"), br(),
                     div(style="display:inline-block",
                             numericInput("sc1c2oup.h", "png / jpeg height:", width = "138px",
                                                        min = 4, max = 50, value = 10, step = 0.5)),
                     div(style="display:inline-block",
                             numericInput("sc1c2oup.w", "png / jpeg width:", width = "138px",
                                                        min = 4, max = 50, value = 10, step = 0.5))
        )    # End of column (6 space)
    ),        # End of fluidRow (4 space)
    fluidRow(
        actionButton("sc1c2tog6", "Toggle to show cell numbers "),
        conditionalPanel(
            condition = "input.sc1c2tog6 % 2 == 1",
            h4("Toggle to show cell number and percentage"),
            dataTableOutput("sc1c2.dt")
        )
    )
),         # End of tab (2 space)

    ### Tab1.d1: Bubbleplot/ Heatmap
    tabPanel(
        HTML("Bubbleplot / Heatmap"),
        h4("Gene expression bubbleplot / heatmap"),
        "In this tab, users can visualise the gene expression patterns of ",
        "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(),
        "The normalised expression are averaged, log-transformed and then plotted.",
        br(),br(),
        fluidRow(
            column(
                3, style="border-right: 2px solid black",
                selectInput("sc1d1grp", HTML("Cell information (X-axis)<br />
                                                        Group / colour by:"),
                                        choices = sc1conf[grp == TRUE]$UI,
                                        selected = sc1conf[grp == TRUE]$UI[1]) %>%
                    helper(type = "inline", size = "m", fade = TRUE,
                                 title = "Cell information to group cells by",
                                 content = c("Select categorical cell information to group cells by",
                                                         "- Single cells are grouped by this categorical covariate",
                                                         "- Plotted as the X-axis of the bubbleplot / heatmap")),
                textAreaInput("sc1d1inp", HTML("List of gene names (Y-axis)<br />
                                                                                    (Max 50 genes, separated <br />
                                                                                     by , or ; or newline):"),
                                            height = "200px",
                                            value = paste0(sc1def$genes, collapse = ", ")) %>%
                    helper(type = "inline", size = "m", fade = TRUE,
                                 title = "List of genes to plot on bubbleplot / heatmap",
                                 content = c("Input genes to plot",
                                                         "- Maximum 50 genes (due to ploting space limitations)",
                                                         "- Genes should be separated by comma, semicolon or newline")),

                radioButtons("sc1d1plt", "Plot type:",
                                         choices = c("Bubbleplot", "Heatmap"),
                                         selected = "Bubbleplot", inline = TRUE),
                actionButton(inputId = "sc1d1update", "Generate Plot",icon("sync")),
                checkboxInput("sc1d1scl", "Scale gene expression", value = TRUE),
                checkboxInput("sc1d1row", "Cluster rows (genes)", value = TRUE),
                checkboxInput("sc1d1col", "Cluster columns (samples)", value = FALSE),
                actionButton("sc1d1tog1", "Toggle to subset cells"),
                conditionalPanel(
                    condition = "input.sc1d1tog1 % 2 == 1",
                    selectInput("sc1d1sub1_1", "Cell information to subset:",
                                choices = sc1conf[grp == TRUE]$UI,
                                selected = sc1def$grp1),
                    uiOutput("sc1d1sub1.ui"),
                    #---------------------
                    actionButton("sc1d1tog2", "Toggle to further subset cells"),
                    conditionalPanel(
                        condition = "input.sc1d1tog2 % 2 == 1",
                        selectInput("sc1d1sub2_1", "Cell information to subset:",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp2),
                        uiOutput("sc1d1sub2.ui"),
                        #----------------------
                        actionButton("sc1d1tog3", "Toggle to further subset cells"),
                        conditionalPanel(
                            condition = "input.sc1d1tog3 % 2 == 1",
                            selectInput("sc1d1sub3_1", "Cell information to subset:",
                                        choices = sc1conf[grp == TRUE]$UI,
                                        selected = sc1def$grp2),
                            uiOutput("sc1d1sub3.ui")
                        )
                    )
                ),
                br(),
                actionButton("sc1d1tog4", "Toggle graphics controls"),
                conditionalPanel(
                    condition = "input.sc1d1tog4 % 2 == 1",
                    selectInput("sc1d1cols", "Select a color specturm:",
                                choices=c(rownames(pal.info)[pal.info$category %in% "div"]),
                                selected="Blue-White-Red"),
                    checkboxInput("sc1d1colinv", "Reverse color", value = FALSE),
                    sliderInput("sc1d1siz", "Point size:",
                                            min = 0, max = 32, value = 8, step = 1),
                    radioButtons("sc1d1psz", "Plot size:",
                                             choices = c("Small", "Medium", "Large", "Extra Large"),
                                             selected = "Medium", inline = TRUE),
                    radioButtons("sc1d1fsz", "Font size:",
                                             choices = c("Extra Small","Small", "Medium", "Large"),
                                             selected = "Medium", inline = TRUE),
                    radioButtons("sc1d1frt", "Rotate x axis label:",
                                             choices = c(0,30,45,90),
                                             selected = 90, inline = TRUE),
                    radioButtons("sc1d1asp", "Aspect ratio:",
                                             choices = c("Square", "Fixed", "Free"),
                                             selected = "Square", inline = TRUE),
                checkboxInput("sc1d1leg", "Show Legend", value = TRUE),
                radioButtons("sc1d1legpos", "Legend positions:",
                             choices = c("top", "right", "bottom"),
                             selected = "bottom", inline = TRUE)
                )

            ),
            downloadButton("sc1d1oup.png", "Download png"),
            downloadButton("sc1d1oup.jpeg", "Download jpeg"), br(),
            div(style="display:inline-block",
                    numericInput("sc1d1oup.h", "png / jpeg height:", width = "138px",
                                             min = 4, max = 50, value = 10, step = 0.5)),
            div(style="display:inline-block",
                    numericInput("sc1d1oup.w", "png / jpeg width:", width = "138px",
                                             min = 4, max = 50, value = 10, step = 0.5)
                ), # End of column (6 space)
            column(9, h4(htmlOutput("sc1d1oupTxt")),
                   withSpinner(uiOutput("sc1d1oup.ui")),br(),
                   downloadButton("sc1d1oup1.xlsx", "Download xlsx")
            )    # End of column (6 space)
        )        # End of fluidRow (4 space)
    ),            # End of tab (2 space)

### Tab1.e1: multiple gene expression trajectory
tabPanel(
    HTML("Trajectory"),
    h4("Gene expression trajectory"),
    "In this tab, users can visualise the multiple genes expression",
    "(e.g. multiple genes expression across groups of cells (e.g. libary / clusters).",
    br(),br(),
    fluidRow(
        column(
            3, style="border-right: 2px solid black",
            textAreaInput("sc1e1inp1", HTML("List of gene names <br/>
                                            (Max 50 genes, separated <br/>
                                            by , or ; or newline):"),
                                        height = "200px",
                                        value = paste0(sc1def$genes, collapse = ", ")) %>%
                helper(type = "inline", size = "m", fade = TRUE,
                             title = "List of genes to plot on bubbleplot / heatmap",
                             content = c("Input genes to plot",
                                                     "- Maximum 50 genes (due to ploting space limitations)",
                                                     "- Genes should be separated by comma, semicolon or newline")),
            selectInput("sc1e1inp2", "Cell information (X-axis):",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp2) %>%
                helper(type = "inline", size = "m", fade = TRUE,
                             title = "Cell information to group cells by",
                             content = c("Select categorical cell information to group cells by",
                                                     "- Single cells are grouped by this categorical covariate",
                                                     "- Plotted as the X-axis of the violin plot / box plot")),
            actionButton("sc1e1tog1", "Toggle to subset cells"),
            conditionalPanel(
                condition = "input.sc1e1tog1 % 2 == 1",
                selectInput("sc1e1sub1_1", "Cell information to subset:",
                                        choices = sc1conf[grp == TRUE]$UI,
                                        selected = sc1def$grp1),
                uiOutput("sc1e1sub1.ui"),
                actionButton("sc1e1tog2", "Toggle to further subset cells"),
                conditionalPanel(
                    condition = "input.sc1e1tog2 % 2 == 1",
                    selectInput("sc1e1sub2_1", "Cell information to subset:",
                                            choices = sc1conf[grp == TRUE]$UI,
                                            selected = sc1def$grp2),
                    uiOutput("sc1e1sub2.ui")
                )
            ),
            br(),
            actionButton("sc1e1tog", "Toggle graphics controls"),
            conditionalPanel(
                condition = "input.sc1e1tog % 2 == 1",
                selectInput("sc1e1cols","Select a color palette:",
                            choices=c("default",rownames(pal.info)[pal.info$category %in% "qual"]),
                            selected="default"),
                checkboxInput("sc1e1norm", "Y-axis shows %", value = FALSE),
                checkboxInput("sc1e1log", "Log transform", value = TRUE),
                checkboxInput("sc1e1y_0", "Y axis starts from zero", value = TRUE),
                radioButtons("sc1e1psz", "Plot size:",
                                         choices = c("Small", "Medium", "Large", "Extra Large"),
                                         selected = "Medium", inline = TRUE),
                checkboxInput("sc1e1pts", "Show data points", value = TRUE),
                sliderInput("sc1e1lvls", "conf.int.level:",
                                        min = 0, max = 1, value = 0.3, step = 0.05),
                radioButtons("sc1e1fsz", "Font size:",
                                         choices = c("Extra Small","Small", "Medium", "Large"),
                                         selected = "Medium", inline = TRUE),
                radioButtons("sc1e1lsz", "line size:",
                                         choices = c("Extra Small","Small", "Medium", "Large"),
                                         selected = "Medium", inline = TRUE),
                radioButtons("sc1e1frt", "Rotate x axis label:",
                                         choices = c(0,30,45,90),
                                         selected = 0, inline = TRUE),
                checkboxInput("sc1e1leg", "Show Legend", value = TRUE))
        ), # End of column (6 space)
        column(9, withSpinner(uiOutput("sc1e1oup.ui")),
                     downloadButton("sc1e1oup.png", "Download png"),
                     downloadButton("sc1e1oup.jpeg", "Download jpeg"), br(),
                     div(style="display:inline-block",
                             numericInput("sc1e1oup.h", "png / jpeg height:", width = "138px",
                                                        min = 4, max = 50, value = 8, step = 0.5)),
                     div(style="display:inline-block",
                             numericInput("sc1e1oup.w", "png / jpeg width:", width = "138px",
                                                        min = 4, max = 50, value = 10, step = 0.5))
        )    # End of column (6 space)
    )        # End of fluidRow (4 space)
),         # End of tab (2 space)

### Tab.f1: run Differential analysis
tabPanel(
    HTML("Volcano Plots"),
    h4("Run Differential Analysis and show results in a Volcano Plot"),
    "In this tab, users can find and visualise differentially expressed genes.",
    br(),br(),
    fluidRow(
        column(
            12, actionButton("sc1f1tog1", "Toggle to subset cells"),
            conditionalPanel(
                condition = "input.sc1f1tog1 % 2 == 1",
                selectInput("sc1f1sub1_1", "Cell information to subset:",
                            choices = sc1conf[grp == TRUE]$UI,
                            selected = sc1def$grp1),
                uiOutput("sc1f1sub1.ui"),
                #-----------------------
                actionButton("sc1f1tog2", "Toggle to further subset cells"),
                conditionalPanel(
                    condition = "input.sc1f1tog2 % 2 == 1",
                    selectInput("sc1f1sub2_1", "Cell information to subset:",
                                choices = sc1conf[grp == TRUE]$UI,
                                selected = sc1def$grp2),
                    uiOutput("sc1f1sub2.ui"),
                    #-----------------------
                    actionButton("sc1f1tog3", "Toggle to further subset cells"),
                    conditionalPanel(
                        condition = "input.sc1f1tog3 % 2 == 1",
                        selectInput("sc1f1sub3_1", "Cell information to subset:",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp2),
                        uiOutput("sc1f1sub3.ui")
                    )
                )
            )
        ), # End of column (12 space)
    ),     # End of fluidRow (1 column)
    fluidRow(
        column(
            12, h4("Select group to compare"),
            fluidRow(
                column(3, selectInput("sc1f1ident", "group:",
                                   choices = sc1conf[grp == TRUE]$UI,
                                   selected = NULL)
                    ),# End of column (4 space)
                column(6,
                       uiOutput("sc1f1ident1.ui") %>%
                           helper(type = "inline", size = "m", fade = TRUE,
                                  title = "Select group 1",
                                  content = c("Select group 1, equivalent to 'ident.1' in Seurat FindMarkers")),
                    uiOutput("sc1f1ident2.ui") %>%
                        helper(type = "inline", size = "m", fade = TRUE,
                               title = "Select group 2",
                               content = c("Select group 2, equivalent to 'ident.2' in Seurat FindMarkers.
                                           If NULL, use all other cells for comparison"))
                ),# End of column (4 space)
            ),# End of fluidRow (8 space)
            tags$head(
                HTML(
                    "
          <script>
          var socket_timeout_interval
          var n = 0
          $(document).on('shiny:connected', function(event) {
          socket_timeout_interval = setInterval(function(){
          Shiny.onInputChange('count', n++)
          }, 15000)
          });
          $(document).on('shiny:disconnected', function(event) {
          clearInterval(socket_timeout_interval)
          });
          </script>
          "
                )
            ),
            actionButton(inputId = "sc1f1update", "Run differential analysis",icon("sync")),
            actionButton("sc1f1tog4", "Toggle graphics controls"),
            conditionalPanel(
                condition = "input.sc1f1tog4 % 2 == 1",
                fluidRow(
                    column(
                        6, sliderInput("sc1f1siz", "Point size:",
                                       min = 0, max = 8, value = 3, step = 0.5),
                        sliderInput("sc1f1alpha", "Point transparency:",
                                    min = 0, max = 1, value = 0.9, step = 0.05),
                        radioButtons("sc1f1psz", "Plot size:",
                                     choices = c("Small", "Medium", "Large", "Extra Large"),
                                     selected = "Large", inline = TRUE),
                        radioButtons("sc1f1fsz", "Font size:",
                                     choices = c("Extra Small","Small", "Medium", "Large"),
                                     selected = "Medium", inline = TRUE),
                        radioButtons("sc1f1asp", "Aspect ratio:",
                                     choices = c("Square", "Fixed", "Free"),
                                     selected = "Square", inline = TRUE),
                        checkboxInput("sc1f1leg", "Show Legend", value = FALSE),
                        radioButtons("sc1f1legpos", "Legend positions:",
                                     choices = c("top", "right", "bottom"),
                                     selected = "bottom", inline = TRUE),
                        checkboxGroupInput("sc1f1rmgene", "Remove unwanted genes:", inline = TRUE,
                                           choices = c("Mitochondrial (MT) genes",
                                                       "Ribosomal protein large subunit (RPL)",
                                                       "Ribosomal protein small subunit (RPS)"),
                                           selected = "Mitochondrial (MT) genes")
                    ), # end of column 6
                    column(
                        6,
                        radioButtons("sc1f1cutp", "Select p-value or adjusted p-value:",
                                     choices = c("p_val_adj","p_val"),
                                     selected = "p_val_adj", inline = TRUE),
                        textInput("sc1f1cutpval", "Select p-value cut off:", value = "0.05"),
                        textInput("sc1f1cutfc", "Select log2FC cut off:", value = "0.25"),
                        radioButtons("sc1f1sort", "Order DEGs by:",
                                     choices = c("p_val_adj","p_val", "avg_log2FC"),
                                     selected = "p_val_adj", inline = TRUE),
                        textInput("sc1f1top", "Select top N DEGs:", value = "15"),
                        selectInput("sc1f1cols", "Select a color spectrum",
                                    choices=c(rownames(pal.info)[pal.info$category %in% c("div")]),
                                    selected= "Blue-White-Red"),
                        checkboxInput("sc1f1colinv", "Reverse color", value = FALSE),
                        radioButtons("sc1f1lab1", "Labels for top DE genes:",
                                     choices = c("No labels", "black text","black labels"),
                                     selected = "black text", inline = FALSE),
                        radioButtons("sc1f1lab2", "Labels for manually input genes:",
                                     choices = c("No labels", "red text","red labels"),
                                     selected = "red text", inline = FALSE)
                    ) # end of column 6
                )
            ), # End of column (12 space)
            actionButton("sc1f1tog5", "Toggle to add genes manually"),
            conditionalPanel(
                condition = "input.sc1f1tog5 % 2 == 1",
                textAreaInput("sc1f1inp", HTML("List of gene names (Y-axis)<br />
                                                    (Max 50 genes, separated <br />
                                                   by , or ; or newline):"),
                              height = "50px",
                              value = paste0("", collapse = ", "))
                ), # End of conditionalPanel
            fluidRow(column(12, #textOutput("keepAlive"),
                            uiOutput("sc1f1oup1.ui"),
                            downloadButton("sc1f1oup1.png", "Download png"),
                            downloadButton("sc1f1oup1.jpeg", "Download jpeg"), br(),
                            div(style="display:inline-block",
                                numericInput("sc1f1oup1.h", "png / jpeg height:", width = "138px",
                                             min = 4, max = 50, value = 13, step = 0.5)),
                            div(style="display:inline-block",
                                numericInput("sc1f1oup1.w", "png / jpeg width:", width = "138px",
                                             min = 4, max = 50, value = 13, step = 0.5)), br(),
                            downloadButton("sc1f1oup1.csv", "Download csv"),
                            downloadButton("sc1f1oup1.xlsx", "Download xlsx"), br(),
                            actionButton("sc1f1tog9", "Toggle to show DE results"),
                            conditionalPanel(
                                condition = "input.sc1f1tog9 % 2 == 1",
                                h4("Toggle to show differential analysis results"),
                                withSpinner(dataTableOutput("sc1f1.dt"))
                                )
                            ), # End of column (12 space)
                )        # End of fluidRow (4 space)
    ),         # End of tab (2 space)
    )
),         # End of tab (2 space)

## Tab1.cor1: Correlation network
tabPanel(
    HTML("Correlation Network"),
    h4("Run correlation analysis and show results in a volcano plot or network"),
    "In this tab, users can visualise the gene-gene Correlation",
    "and network relation based correlation.",
    br(),br(),
    fluidRow(
        column(
            12,
            h4("Step One: Select group of cells for correlation and run it"),
            actionButton("sc1n1tog1", "Toggle to subset cells"),
            conditionalPanel(
                condition = "input.sc1n1tog1 % 2 == 1",
                selectInput("sc1n1sub1_1", "Cell information to subset:",
                            choices = sc1conf[grp == TRUE]$UI,
                            selected = sc1def$grp1),
                uiOutput("sc1n1sub1.ui"),
                #-----------------------
                actionButton("sc1n1tog2", "Toggle to further subset cells"),
                conditionalPanel(
                    condition = "input.sc1n1tog2 % 2 == 1",
                    selectInput("sc1n1sub2_1", "Cell information to subset:",
                                choices = sc1conf[grp == TRUE]$UI,
                                selected = sc1def$grp2),
                    uiOutput("sc1n1sub2.ui"),
                    #-----------------------
                    actionButton("sc1n1tog3", "Toggle to further subset cells"),
                    conditionalPanel(
                        condition = "input.sc1n1tog3 % 2 == 1",
                        selectInput("sc1n1sub3_1", "Cell information to subset:",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp2),
                        uiOutput("sc1n1sub3.ui")
                    )
                )
            ),
            sliderInput("sc1n1minexpr", "Select genes expressed in at least % of the cells",min = 0, max = 100, value = 5, step = 1),
            actionButton(inputId = "sc1n1update1", "Run correlation analysis,Click twice",icon("sync")),
        ), # End of column (12 space)
    ),     # End of fluidRow (1 column)
    fluidRow(
        column(
            12,
            h4("Step Two: Correlation volcano plot"),
            fluidRow(column(12, #textOutput("keepAlive"),
                    selectInput("sc1n1inp1", "Gene name:", choices=NULL,
                                selected = sc1def$gene2,
                                multiple=FALSE),
                    actionButton(inputId = "sc1n1update2", "Generate volcano plot",icon("sync")),
                    actionButton("sc1n1tog4", "Toggle graphics controls"),
                    conditionalPanel(
                        condition = "input.sc1n1tog4 % 2 == 1",
                        fluidRow(
                            column(
                                6, sliderInput("sc1n1siz1", "Point size:",
                                               min = 0, max = 8, value = 4, step = 0.25),
                                sliderInput("sc1n1alpha1", "Point transparency:",
                                            min = 0, max = 1, value = 0.8, step = 0.05),
                                radioButtons("sc1n1psz1", "Plot size:",
                                             choices = c("Small", "Medium", "Large", "Extra Large"),
                                             selected = "Large", inline = TRUE),
                                radioButtons("sc1n1fsz1", "Font size:",
                                             choices = c("Extra Small","Small", "Medium", "Large"),
                                             selected = "Medium", inline = TRUE),
                                radioButtons("sc1n1asp1", "Aspect ratio:",
                                             choices = c("Square", "Fixed", "Free"),
                                             selected = "Square", inline = TRUE)
                                #checkboxInput("sc1n1leg", "Show Legend", value = FALSE),
                                #radioButtons("sc1n1legpos", "Legend positions:",
                                #             choices = c("top", "right", "bottom"),
                                #             selected = "bottom", inline = TRUE),
                                #checkboxGroupInput("sc1n1rmgene", "Remove unwanted genes:", inline = TRUE,
                                #                   choices = c("Mitochondrial (MT) genes",
                                #                               "Ribosomal protein large subunit (RPL)",
                                #                               "Ribosomal protein small subunit (RPS)"),
                                #                   selected = "Mitochondrial (MT) genes")
                            ), # end of column 6
                            column(
                                6,
                                radioButtons("sc1n1cutp1", "Select p-value or adjusted p-value:",
                                             choices = c("p_val_adj","p_val"),
                                             selected = "p_val_adj", inline = TRUE),
                                textInput("sc1n1cutpval1", "Select p-value cut off:", value = "0.05"),
                                textInput("sc1n1cutfc1", "Select correlation cut off:", value = "0.05"),
                                textInput("sc1n1top1", "Select top N correlated genes:", value = "10"),
                                #selectInput("sc1n1cols", "Select a color spectrum",
                                #            choices=c(rownames(pal.info)[pal.info$category %in% c("div")]),
                                #            selected= "Blue-White-Red"),
                                selectInput("sc1n1lab1", "Labels for top correlated genes:",
                                            choices=c("black",cList[["Blue-Yellow-Red"]]),
                                            selected= "black"),
                                selectInput("sc1n1lab2", "Labels for manually input genes:",
                                            choices=c(cList[["Blue-Yellow-Red"]],"black"),
                                            selected = "red")
                            ) # end of column 6
                        )
                    ), # End of column (12 space)
                    actionButton("sc1n1tog5", "Toggle to add genes manually"),
                    conditionalPanel(
                        condition = "input.sc1n1tog5 % 2 == 1",
                        textAreaInput("sc1n1inp", HTML("List of gene names (Y-axis)<br />
                                            (Max 50 genes, separated <br />
                                           by , or ; or newline):"),
                                      height = "50px",
                                      value = paste0("", collapse = ", "))
                    ), # End of conditionalPanel
                    uiOutput("sc1n1oup1.ui"),
                    downloadButton("sc1n1oup1.png", "Download png"),
                    downloadButton("sc1n1oup1.jpeg", "Download jpeg"), br(),
                    div(style="display:inline-block",
                        numericInput("sc1n1oup1.h", "png / jpeg height:", width = "138px",
                                     min = 4, max = 50, value = 13, step = 0.5)),
                    div(style="display:inline-block",
                        numericInput("sc1n1oup1.w", "png / jpeg width:", width = "138px",
                                     min = 4, max = 50, value = 13, step = 0.5)), br(),
                    h4("Step Three: Correlation network"),
                    selectInput("sc1n1inp2", "Gene names:", choices=NULL,selected = "EZH1,EZH2",
                                multiple=TRUE),
                    actionButton(inputId = "sc1n1update3", "Generate Network",icon("sync")),
                    actionButton("sc1n1tog6", "Toggle graphics controls"),
                    conditionalPanel(
                        condition = "input.sc1n1tog6 % 2 == 1",
                        fluidRow(
                            column(
                                6, sliderInput("sc1n1linklen", "link Distance:",
                                               min = 50, max = 250, value = 120, step = 10),
                                sliderInput("sc1n1linkwid", "link Width:",
                                            min = 0, max = 4, value = 2, step = 0.25),
                                radioButtons("sc1n1psz2", "Plot size:",
                                             choices = c("Small", "Medium", "Large", "Extra Large"),
                                             selected = "Large", inline = TRUE),
                                radioButtons("sc1n1fsz2", "Font size:",
                                             choices = c("Extra Small","Small", "Medium", "Large"),
                                             selected = "Medium", inline = TRUE),
                                checkboxInput("sc1n1name", "Show gene name", value = TRUE),
                                checkboxInput("sc1n1leg", "Show Legend", value = FALSE)
                                #radioButtons("sc1n1legpos", "Legend positions:",
                                #             choices = c("top", "right", "bottom"),
                                #             selected = "bottom", inline = TRUE),
                                #checkboxGroupInput("sc1n1rmgene", "Remove unwanted genes:", inline = TRUE,
                                #                   choices = c("Mitochondrial (MT) genes",
                                #                               "Ribosomal protein large subunit (RPL)",
                                #                               "Ribosomal protein small subunit (RPS)"),
                                #                   selected = "Mitochondrial (MT) genes")
                            ), # end of column 6
                            column(
                                6,
                                radioButtons("sc1n1cutp2", "Select p-value or adjusted p-value:",
                                             choices = c("p_val_adj","p_val"),
                                             selected = "p_val_adj", inline = TRUE),
                                textInput("sc1n1cutpval2", "Select p-value cut off:", value = "0.05"),
                                textInput("sc1n1cutfc2", "Select correlation cut off:", value = "0.05"),
                                textInput("sc1n1top2", "Select top N correlated genes:", value = "10"),
                                radioButtons("sc1n1posneg", "Select positive and/or negative correlated genes only:",
                                             choices = c("positive only","negative only","both"),
                                             selected = "positive only", inline = TRUE),
                                selectInput("sc1n1col2","Select a color palette",
                                            choices=c("default",rownames(pal.info)[pal.info$category %in% c("div")]),
                                            selected="Blue-White-Red"),
                                checkboxInput("sc1n1col2inv", "Reverse color", value = FALSE),

                            ) # end of column 6
                        )
                    ), # End of column (12 space)
                    forceNetworkOutput("sc1n1oup2.ui",height = "500px"),
                    downloadButton("sc1n1oup2.network", "Download html")
                    #downloadButton("sc1n1oup2.csv", "Download csv"),
                    #downloadButton("sc1n1oup2.xlsx", "Download xlsx"), br(),

            ), # End of column (12 space)
            )        # End of fluidRow (4 space)
        ),         # End of tab (2 space)
    )        # End of fluidRow (4 space)
),            # End of tab (2 space)

### Tab1.t1: TCR plot
tabPanel(
    HTML("TCR plot"),
    h4("T-cell receptor (TCR) and immunoglobulin (Ig) enrichment information"),
    "In this tab, users can visualise the single-cell immune receptor profiling",
    br(),br(),
    fluidRow(
        column(
            3, style="border-right: 2px solid black",
            selectInput("sc1t1inp1", "Cell information to plot (X-axis):",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp2),
            selectInput("sc1t1inp2", "Cell information to group:",
                        choices = sc1conf[grp == TRUE]$UI,
                        selected = sc1def$grp1) %>%
                helper(type = "inline", size = "m", fade = TRUE,
                       title = "Cell information to group / colour cells by",
                       content = c("Select categorical cell information to group / colour cells by",
                                   "- Proportion / cell numbers are shown in different colours")),
            radioButtons("sc1t1typ", "Plot Type:",
                         choices = c("Unique Barplot","Cumulative clone sizes distribution",
                                     "Proportion Barplot","Paired Diversity","unpaired Diversity",
                                     "Paired scatter Clonotype","unpaired scatter Clonotype",
                                     "Novel & Expand & Contract Barplot_1",
                                     "Novel & Expand & Contract Barplot_2"),
                         selected = c("Unique Barplot"), inline = FALSE),
            checkboxInput("sc1t1norm", "Y-axis shows %", value = FALSE),
            radioButtons("sc1t1cloneCall", "Definition of clonotype:",
                         choices = c("TCR/Ig genes", "CDR3 nucleotide","CDR3 amino acid",
                                     "TCR/Ig + CDR3 nucleotide"),
                         selected = "TCR/Ig + CDR3 nucleotide", inline = TRUE),
            actionButton("sc1t1tog1", "Toggle to subset cells"),
            conditionalPanel(
                condition = "input.sc1t1tog1 % 2 == 1",
                selectInput("sc1t1sub1_1", "Cell information to subset:",
                            choices = sc1conf[grp == TRUE]$UI,
                            selected = sc1def$grp1),
                uiOutput("sc1t1sub1.ui"),
                actionButton("sc1t1tog2", "Toggle to further subset cells"),
                conditionalPanel(
                    condition = "input.sc1t1tog2 % 2 == 1",
                    selectInput("sc1t1sub2_1", "Cell information to subset:",
                                choices = sc1conf[grp == TRUE]$UI,
                                selected = sc1def$grp2),
                    uiOutput("sc1t1sub2.ui"),
                    actionButton("sc1t1tog3", "Toggle to further subset cells"),
                    conditionalPanel(
                        condition = "input.sc1t1tog3 % 2 == 1",
                        selectInput("sc1t1sub3_1", "Cell information to subset:",
                                    choices = sc1conf[grp == TRUE]$UI,
                                    selected = sc1def$grp2),
                        uiOutput("sc1t1sub3.ui")
                    )
                )
            ),
            br(),
            actionButton("sc1t1tog4", "Toggle graphics controls"),
            conditionalPanel(
                condition = "input.sc1t1tog4 % 2 == 1",
                sliderInput("sc1t1lvls", "Significant level for Clone change:",
                            min = 0, max = 0.1, value = 0.03, step = 0.001),
                sliderInput("sc1t1siz", "Data point size",
                            min = 0, max = 8, value = 2, step = 0.25),
                selectInput("sc1t1cols","Select a color palette:",
                            choices=c("default",rownames(pal.info)[pal.info$category %in% "qual"]),
                            selected="default"),
                selectInput("sc1t1divindex", "Diversity index",
                            choices = c("Shannon","Simpson","Inv.Simpson","Chao","ACE"),
                            selected = "Shannon"),
                radioButtons("sc1t1psz", "Plot size:",
                             choices = c("Small", "Medium", "Large"),
                             selected = "Medium", inline = TRUE),
                checkboxInput("sc1t1pts", "Show data points", value = TRUE),
                radioButtons("sc1t1fsz", "Font size:",
                             choices = c("Extra Small","Small", "Medium", "Large"),
                             selected = "Medium", inline = TRUE),
                radioButtons("sc1t1lsz", "Line size:",
                             choices = c("Extra Small","Small", "Medium", "Large"),
                             selected = "Medium", inline = TRUE),
                checkboxInput("sc1t1flp", "Flip X/Y", value = FALSE),
                radioButtons("sc1t1frt", "Rotate x axis label:",
                             choices = c(0,30,45,90),
                             selected = 0, inline = TRUE),
                checkboxInput("sc1t1leg", "Show Legend", value = TRUE),
                radioButtons("sc1t1legpos", "Legend positions:",
                             choices = c("top", "right", "bottom"),
                             selected = "bottom", inline = TRUE),
            )

        ), # End of column (6 space)
        column(9, withSpinner(uiOutput("sc1t1oup.ui")),
               downloadButton("sc1t1oup.png", "Download png"),
               downloadButton("sc1t1oup.jpeg", "Download jpeg"), br(),
               div(style="display:inline-block",
                   numericInput("sc1t1oup.h", "png / jpeg height:", width = "138px",
                                min = 4, max = 50, value = 10, step = 0.5)),
               div(style="display:inline-block",
                   numericInput("sc1t1oup.w", "png / jpeg width:", width = "138px",
                                min = 4, max = 50, value = 10, step = 0.5))
        )    # End of column (6 space)
    )        # End of fluidRow (4 space)
),         # End of tab (2 space)
br(),
p("", style = "font-size: 125%;"),
p(em("This webpage was made using "), a("ShinyCell",
                                        href = "https://github.com/nyuhuyang/ShinyCell",target="_blank")),
br(),br(),br(),br(),br()
)))



