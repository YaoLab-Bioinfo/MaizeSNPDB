

library(shinythemes)
library(shinyBS)

shinyUI(
  navbarPage(
    
    title = "MaizeSNPDB", theme = shinytheme("flatly"), 
    windowTitle = "SNP database of 1210 maize lines",
    
    ## About
    tabPanel("About", includeMarkdown("About.md")),
    
    # Genome browser
    tabPanel(
      "Browser",
      
      tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});')),
        tags$style(
          HTML(
            "
            #inputs-table {
            border-collapse: collapse;
            }
            
            #inputs-table td {
            padding: 3px;
            vertical-align: bottom;
            }

            .multicol .shiny-options-group{
                            -webkit-column-count: 2; /* Chrome, Safari, Opera */
            -moz-column-count: 2;    /* Firefox */
            column-count: 2;
            -moz-column-fill: balanced;
            -column-fill: balanced;
            }
            .checkbox{
            margin-top: 0px !important;
            -webkit-margin-after: 1px !important; 
            }

            "
          ) #/ HTML
          ) #/ style
          ), #/ head
      
      
      fluidRow(column(12,
                class = "col-md-5",
                style = "margin: 1px 1px 1px 1px",
                tags$table(id = "inputs-table",
                             style = "width: 100%",
                             tags$tr(
                               tags$td(style = "width: 60%",
                                       textInput("regB", label = h5("Genomic region:",
                                                                    bsButton("q6", label="", icon=icon("question"), style="info", size="small")),
                                                 value = "chr1:29765419-29793053"),
                                       
                                       bsPopover("q6", "A genomic region can be determined by chromosome positions or gene locus. For example, chr7:29506705-29659223 or Zm00001d003489.",
                                                 trigger = "focus")
                               ), #/ column 1
                               tags$td(style = "width: 40%; text-align: right",
                                         div(class = "form-group shiny-input-container",
                                               actionButton("submit1", strong("Go!",
                                                                              bsButton("q7", label="", icon=icon("question"), style="info", size="small")
                                                            ), width = "90%", styleclass = "success"),
                                             bsPopover("q7", "Whenever the genomic region is updated, please click Go!",
                                                       trigger = "focus")
                                         )
                               ) #/ column 2
                             ) #/ tr
                ) #/ table
      )),
      
#      br(),
      
      downloadButton("downloadsnp.txt", "Download genotype data"),
      downloadButton("downloadsnpInfo.txt", "Download SNPs information"),
      downloadButton("downloadGB.pdf", "Download pdf-file"),
      
      withSpinner(plotlyOutput("gbrowser", height = '100%', width = '100%')),
      
      fluidRow(
        column(3,
               sliderInput("GBUP", h4("Upstream:",
                                      bsButton("qg2", label="", icon=icon("question"), style="info", size="small")
                                      ), min = 0, max = 50000, value = 0, ticks = FALSE),
               bsPopover("qg2", "SNPs in the upstream of the specified genomic region will be used.",
                         trigger = "focus"),
               sliderInput("GBDOWN", h4("Downstream:",
                                        bsButton("qg4", label="", icon=icon("question"), style="info", size="small")
                                        ), min = 0, max = 50000, value = 0, ticks = FALSE),
               bsPopover("qg4", "SNPs in the downstream of the specified genomic region will be used.",
                         trigger = "focus")
        ),
        
        column(4,
               p(h4("Select maize lines:",
                    bsButton("qg3", label="", icon=icon("question"), style="info", size="small"))),
               bsPopover("qg3", "Only the chosen maize lines will be used.",
                         trigger = "focus"),
               
               chooserInput("mychooserB", "Available frobs", "Selected frobs", c(),
                            all.acc.cho, size = 10, multiple = TRUE)
        ),
        
        column(5,
               tags$div(align = 'left',
                        class = 'multicol', style = "width: 100%",
                        checkboxGroupInput("GB_mut_group", h4("Mutation types:",
                                                              bsButton("qg1", label="", icon=icon("question"), style="info", size="small")),
                                           choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                       "Intron", "Upstream", "Downstream", "Intergenic",
                                                       "five_prime_UTR","three_prime_UTR",
                                                       "Non_synonymous_start","Non_synonymous_coding",
                                                       "Splice_site_acceptor","Splice_site_donor",
                                                       "Synonymous_stop","Synonymous_coding"
                                                       ),
                                           selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                        "Intron", "Upstream", "Downstream", "Intergenic",
                                                        "five_prime_UTR","three_prime_UTR",
                                                        "Non_synonymous_start","Non_synonymous_coding",
                                                        "Splice_site_acceptor","Splice_site_donor",
                                                        "Synonymous_stop","Synonymous_coding"
                                           )),
                        bsPopover("qg1", "Only SNPs with selected mutation effects will be used.",
                                  trigger = "focus")
               )
        )
      ),
      
      br()
      
    ),
    
    # LDheatmap
    tabPanel(
      "LDheatmap",
      
      sidebarPanel(
        textInput("regL", label = h5("Genomic region:",
                                     bsButton("q5", label="", icon=icon("question"), style="info", size="small")),
                  value = "Zm00001d033673"),
        
        bsPopover("q5", "A genomic region can be determined by chromosome positions or gene locus. For example, chr7:29506705-29659223 or Zm00001d003489.",
                  trigger = "focus"),
        
        actionButton("submit2", strong("Go!",
                                       bsButton("q8", label="", icon=icon("question"), style="info", size="small")
        ), styleclass = "success"),
        bsPopover("q8", "Whenever the genomic region is updated, please click Go!",
                  trigger = "focus"),
        
        br(),
        HTML("<h4><font color='red'>Plot options</font></h4>"),
        
        radioButtons("flip", "Flip the figure", list("FALSE" = 0, "TRUE" = 1)),
        
        conditionalPanel(
          condition = "input.flip==1",
          checkboxInput("LDshowGene", "Show gene model", FALSE),
          conditionalPanel(
            condition = "input.LDshowGene",
            numericInput("ldY", "Y:", value = 72),
            numericInput("ldW", "W:", value = 72)
          )
        ),
        
        conditionalPanel(
          condition = "input.flip==0",
          radioButtons("showText", "Print LD measurements", list("FALSE" =  0, "TRUE" = 1)),
          textInput("ldpos", "Label SNPs:", value = "5, 8")
        ),
        
        radioButtons("ldcol",
          "Color", list("grey.colors(20)" = 1, "heat.colors(20)" = 2)
        ),
        
        numericInput("ldUp", h4("Upstream (kb):",
                        bsButton("ql4", label="", icon=icon("question"), style="info", size="small")
                     ), value = 0),
        bsPopover("ql4", "SNPs in the upstream of the specified genomic region will be used.",
                trigger = "focus"),
        numericInput("ldDown", h4("Downstream (kb):",
                                  bsButton("ql5", label="", icon=icon("question"), style="info", size="small")
        ), value = 0),
        bsPopover("ql5", "SNPs in the downstream of the specified genomic region will be used.",
            trigger = "focus"),

        tags$div(align = 'left',
           class = 'multicol', style = "width: 100%",
           checkboxGroupInput("ld_mut_group", h4("Mutation types:",
                                                 bsButton("ql1", label="", icon=icon("question"), style="info", size="small")),
                            choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                        "Intron", "Upstream", "Downstream", "Intergenic",
                                        "five_prime_UTR","three_prime_UTR",
                                        "Non_synonymous_start","Non_synonymous_coding",
                                        "Splice_site_acceptor","Splice_site_donor",
                                        "Synonymous_stop","Synonymous_coding"
                            ),
                            selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                         "Intron", "Upstream", "Downstream", "Intergenic",
                                         "five_prime_UTR","three_prime_UTR",
                                         "Non_synonymous_start","Non_synonymous_coding",
                                         "Splice_site_acceptor","Splice_site_donor",
                                         "Synonymous_stop","Synonymous_coding"
                            )),
           bsPopover("ql1", "Only SNPs with selected mutation effects will be used.",
                     trigger = "focus")
        ),

        p(h4("Select maize lines:",
            bsButton("ql3", label="", icon=icon("question"), style="info", size="small"))),
        bsPopover("ql3", "Only the chosen maize lines will be used.",
            trigger = "focus"),

        chooserInput("mychooserLD", "Available frobs", "Selected frobs",
              c(), all.acc.cho, size = 10, multiple = TRUE
        ),

        radioButtons("uploadLD", h4("SNP sites to be retained:",
                                    bsButton("ql2", label="", icon=icon("question"), style="info", size="small")), 
                     c("ALL" = "1", "Upload SNP sites file" = "2"), "1"),
        bsPopover("ql2", "A text file with SNP IDs (one ID per row) could be uploaded to screen the SNPs used in the analysis. Or else, all the SNPs in the specifid genomic region will be used.",
              trigger = "focus"),
        conditionalPanel(condition="input.uploadLD == '2'",
                 fileInput("LD.snpsite", NULL, multiple = FALSE)),

        checkboxInput("ldSize", "Adjust plot size", FALSE),
        conditionalPanel(
          condition = "input.ldSize",
          numericInput("ldHeight", "Plot height:", value = 550),
          numericInput("ldWidth", "Plot width:", value = 750)
        )
      ),
      
      mainPanel(
        downloadButton("downloadLD.pdf", "Download pdf-file"),
        downloadButton("downloadLD.svg", "Download svg-file"),
        withSpinner((plotOutput("ldheatmap", height = '100%', width = '100%')))
        
      )
    ),
    
 
    # Nucleotide diversity
    tabPanel(
      "Diversity",
      
      sidebarPanel(
        textInput("regD", label = h5("Genomic region:",
                                     bsButton("q3", label="", icon=icon("question"), style="info", size="small")),
                  value = "Zm00001d033673"),
        
        bsPopover("q3", "A genomic region can be determined by chromosome positions or gene locus. For example, chr7:29506705-29659223 or Zm00001d003489.",
                  trigger = "focus"),
        
        HTML("<h4><font color='red'>Plot options</font></h4>"),
        
        numericInput("snpnumD", h5("Number of SNPs in each window:",
                                   bsButton("qd6", label="", icon=icon("question"), style="info", size="small")
                                   ), value = 10, min = 5, max = 20),
        bsPopover("qd6", "A specified genomic region would be split into non-overlapping window so that each window contains specified number of SNPs. The nucleotide diversity of all maize lines belong to the specified ecotypes in each window would be calculated.",
              trigger = "focus"),
        
        tags$div(align = 'left',
                 class = 'multicol', style = "width: 100%",
                 checkboxGroupInput("div_acc_group", "Ecotypes to calculate diversity:",
                                    choices = c("Improved", "Landrace", "Parviglumis"),
                                    selected = c("Improved", "Parviglumis")) 
        ),
        
        selectInput("nuc_numerator", h5("Numerator ecotype:",
                                        bsButton("qd7", label="", icon=icon("question"), style="info", size="small")
                                        ), choices = 
                      c("Parviglumis", "Improved", "Landrace")),
        bsPopover("qd7", "The nucleotide diversity of maize lines belong to the Numerator ecotype would be divided by the nucleotide diversity of maize lines belong to the Denominator ecotype for comparison.",
              trigger = "focus"),
        selectInput("nuc_denominator", "Denominator ecotype:", choices = 
                      c("Improved", "Landrace", "Parviglumis")),
        
        tags$div(align = 'left',
          class = 'multicol', style = "width: 100%",
          checkboxGroupInput("div_mut_group", h4("Mutation types:",
                                                 bsButton("qd1", label="", icon=icon("question"), style="info", size="small")),
                             choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                         "Intron", "Upstream", "Downstream", "Intergenic",
                                         "five_prime_UTR","three_prime_UTR",
                                         "Non_synonymous_start","Non_synonymous_coding",
                                         "Splice_site_acceptor","Splice_site_donor",
                                         "Synonymous_stop","Synonymous_coding"
                             ),
                             selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                          "Intron", "Upstream", "Downstream", "Intergenic",
                                          "five_prime_UTR","three_prime_UTR",
                                          "Non_synonymous_start","Non_synonymous_coding",
                                          "Splice_site_acceptor","Splice_site_donor",
                                          "Synonymous_stop","Synonymous_coding"
                             )),
          bsPopover("qd1", "Only SNPs with selected mutation effects will be used.",
                    trigger = "focus")
        ),

        
        numericInput("divUp", h4("Upstream (kb):",
                                 bsButton("qd4", label="", icon=icon("question"), style="info", size="small")
        ), value = 20),
        bsPopover("qd4", "SNPs in the upstream of the specified genomic region will be used.",
              trigger = "focus"),
        numericInput("divDown", h4("Downstream (kb):",
                                   bsButton("qd5", label="", icon=icon("question"), style="info", size="small")
        ), value = 20),
        bsPopover("qd5", "SNPs in the downstream of the specified genomic region will be used.",
            trigger = "focus"),

        radioButtons("uploadDIV", h4("SNP sites to be retained:",
                                     bsButton("qd2", label="", icon=icon("question"), style="info", size="small")), 
                     c("ALL" = "1", "Upload SNP sites file" = "2"), "1"),
        bsPopover("qd2", "A text file with SNP IDs (one ID per row) could be uploaded to screen the SNPs used in the analysis. Or else, all the SNPs in the specifid genomic region will be used.",
            trigger = "focus"),
        conditionalPanel(condition="input.uploadDIV == '2'",
                 fileInput("DIV.snpsite", NULL, multiple = FALSE)),

        checkboxInput("divSize", "Adjust plot size", FALSE),
        conditionalPanel(
          condition = "input.divSize",
          numericInput("divHeight", "Plot height:", value = 550),
          numericInput("divWidth", "Plot width:", value = 750)
        ),

        actionButton("submit4", strong("Go!",
                               bsButton("q10", label="", icon=icon("question"), style="info", size="small")
        ), styleclass = "success"),
        conditionalPanel(condition="input.submit4 != '0'", busyIndicator(HTML("<div style='color:red;font-size:30px'>Calculation In progress...</div>"), wait = 0)),
        bsPopover("q10", "Whenever the genomic region or any option is updated, please click Go!!",
          trigger = "focus")
        
      ),
      
      mainPanel(
        downloadButton("downloadDiv.pdf", "Download pdf-file"),
        downloadButton("downloadDiv.svg", "Download svg-file"),
        downloadButton("downloadDiv.txt", "Download TXT-file"),
        plotOutput("diversity", height = '100%', width = '100%')
      )
      
    ),
    
    # Phylogenetic tree
    tabPanel(
      "Phylogenetic",
      
      sidebarPanel(
        textInput("regP", label = h5("Genomic region:",
                                     bsButton("q2", label="", icon=icon("question"), style="info", size="small")),
                  value = "Zm00001d033673"),
        
        bsPopover("q2", "A genomic region can be determined by chromosome positions or gene locus. For example, chr7:29506705-29659223 or Zm00001d003489.",
                  trigger = "focus"),
        
        HTML("<h4><font color='red'>Plot options</font></h4>"),
        numericInput("phyUp", h4("Upstream (kb):",
                                 bsButton("qp4", label="", icon=icon("question"), style="info", size="small")
        ), value = 20),
        bsPopover("qp4", "SNPs in the upstream of the specified genomic region will be used.",
                  trigger = "focus"),
        numericInput("phyDown", h4("Downstream (kb):",
                                   bsButton("qp5", label="", icon=icon("question"), style="info", size="small")
        ), value = 20),
        bsPopover("qp5", "SNPs in the downstream of the specified genomic region will be used.",
                  trigger = "focus"),
        
        tags$div(align = 'left',
                 class = 'multicol', style = "width: 100%",
                 checkboxGroupInput("phy_mut_group", h4("Mutation types:",
                                                        bsButton("qp1", label="", icon=icon("question"), style="info", size="small")),
                                    choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                "Intron", "Upstream", "Downstream", "Intergenic",
                                                "five_prime_UTR","three_prime_UTR",
                                                "Non_synonymous_start","Non_synonymous_coding",
                                                "Splice_site_acceptor","Splice_site_donor",
                                                "Synonymous_stop","Synonymous_coding"
                                    ),
                                    selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                 "Intron", "Upstream", "Downstream", "Intergenic",
                                                 "five_prime_UTR","three_prime_UTR",
                                                 "Non_synonymous_start","Non_synonymous_coding",
                                                 "Splice_site_acceptor","Splice_site_donor",
                                                 "Synonymous_stop","Synonymous_coding"
                                    )),
                 bsPopover("qp1", "Only SNPs with selected mutation effects will be used.",
                           trigger = "focus")
        ),
        
        p(h4("Select maize lines:",
             bsButton("qp3", label="", icon=icon("question"), style="info", size="small"))),
        bsPopover("qp3", "Only the chosen maize lines will be used.",
                  trigger = "focus"),
        
        chooserInput("mychooserPhy", "Available frobs", "Selected frobs",
                     c(), all.acc.cho, size = 10, multiple = TRUE
        ),

        radioButtons("uploadPHY", h4("SNP sites to be retained:",
                                     bsButton("qp2", label="", icon=icon("question"), style="info", size="small")), 
                     c("ALL" = "1", "Upload SNP sites file" = "2"), "1"),
bsPopover("qp2", "A text file with SNP IDs (one ID per row) could be uploaded to screen the SNPs used in the analysis. Or else, all the SNPs in the specifid genomic region will be used.",
          trigger = "focus"),
        conditionalPanel(condition="input.uploadPHY == '2'",
                 fileInput("PHY.snpsite", NULL, multiple = FALSE)),

        checkboxInput("phySize", "Adjust plot size", FALSE),
        conditionalPanel(
          condition = "input.phySize",
          numericInput("phyHeight", "Plot height:", value = 700),
          numericInput("phyWidth", "Plot width:", value = 750)
        ),        

        actionButton("submit5", strong("Go!",
                                       bsButton("q11", label="", icon=icon("question"), style="info", size="small")
        ), styleclass = "success"),
        conditionalPanel(condition="input.submit5 != '0'", busyIndicator(HTML("<div style='color:red;font-size:30px'>Calculation In progress...</div>"), wait = 0)),
        bsPopover("q11", "Whenever the genomic region or any option is updated, please click Go!",
                  trigger = "focus")
      ),
      
      mainPanel(
        downloadButton("downloadPhylo.pdf", "Download pdf-file"),
        downloadButton("downloadPhylo.nwk", "Download Newick-file"),
        plotOutput("phylo", height = '100%', width = '100%')
        
      )
      
    ),
    
    # Accession
    tabPanel(
      "Accession",
      
      sidebarPanel(
        downloadButton("acc.info.txt", "Download information of all accessions"),
        
        p(h4("Select maize lines:",
             bsButton("qa1", label="", icon=icon("question"), style="info", size="small"))),
        bsPopover("qa1", "Only the chosen maize lines will be used.",
                  trigger = "focus"),
        
        chooserInput("mychooserA", "Available frobs", "Selected frobs",
                     c(), all.acc.cho, size = 10, multiple = TRUE
        ),
        
        br(),
        downloadButton("sel.acc.info.txt", "Download information of selected accessions")
      ),
      
      mainPanel(
        h4("Information of selected maize lines"),
        dataTableOutput("mytable1")
        )
      ),
    
    # Bulk download of data
    tabPanel(
      "Download",
      
      mainPanel(
        fluidRow(column(12,
                        class = "col-md-5",
                        style = "margin: 1px 1px 1px 1px",
                        tags$table(id = "inputs-table",
                                   style = "width: 100%",
                                   tags$tr(
                                     tags$td(style = "width: 60%",
                                             textInput("regBB", label = h5("Genomic region:",
                                                                          bsButton("q1", label="", icon=icon("question"), style="info", size="small")),
                                                       value = "chr1:29611303-29639223"),
                                             
                                             bsPopover("q1", "A genomic region can be determined by chromosome positions or gene locus. For example, chr7:29506705-29659223 or Zm00001d003489.",
                                                       trigger = "focus")
                                     ), #/ column 1
                                     tags$td(style = "width: 40%; text-align: right",
                                             div(class = "form-group shiny-input-container",
                                                 actionButton("submit6", strong("Go!",
                                                                                bsButton("q12", label="", icon=icon("question"), style="info", size="small")
                                                 ), width = "90%", styleclass = "success"),
                                                 bsPopover("q12", "Whenever the genomic region or any option is updated, please click Go!",
                                                           trigger = "focus")
                                             )
                                     ) #/ column 2
                                   ) #/ tr
                        ) #/ table
        )),
        
        fluidRow(
          column(7,
                 p(h4("Select maize lines:",
                           bsButton("qdl2", label="", icon=icon("question"), style="info", size="small"))),
                 bsPopover("qdl2", "Only the chosen maize lines will be used.",
                           trigger = "focus"),
                 
                 chooserInput("mychooserD", "Available frobs", "Selected frobs",
                              c(), all.acc.cho, size = 10, multiple = TRUE
                 )),
          
          column(5,
                 tags$div(align = 'left',
                          class = 'multicol', style = "width: 100%",
                          checkboxGroupInput("down_mut_group", h4("Mutation types:",
                                                                  bsButton("qdl1", label="", icon=icon("question"), style="info", size="small")),
                                             choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                         "Intron", "Upstream", "Downstream", "Intergenic",
                                                         "five_prime_UTR","three_prime_UTR",
                                                         "Non_synonymous_start","Non_synonymous_coding",
                                                         "Splice_site_acceptor","Splice_site_donor",
                                                         "Synonymous_stop","Synonymous_coding"
                                             ),
                                             selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                          "Intron", "Upstream", "Downstream", "Intergenic",
                                                          "five_prime_UTR","three_prime_UTR",
                                                          "Non_synonymous_start","Non_synonymous_coding",
                                                          "Splice_site_acceptor","Splice_site_donor",
                                                          "Synonymous_stop","Synonymous_coding"
                                             )),
                          bsPopover("qdl1", "Only SNPs with selected mutation effects will be used.",
                                    trigger = "focus")
          ))
        ),
        
        br(),
        
        downloadButton("bulkdownloadsnpInfo.txt", "Download SNPs information"),
        downloadButton("bulkdownloadsnp.txt", "Download genotype data"),
        downloadButton("bulkdownloadgene.txt", "Download gene annotation"),
        h4("SNPs information in specified genomic region:"),
        withSpinner(dataTableOutput("mytable2"))
        
      )
    ),

  ## Help
  tabPanel("Help", includeMarkdown("README.md"))

  )
)

