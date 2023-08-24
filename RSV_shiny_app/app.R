#RSV WGS QC R shiny app
# 1: TAB 1:
#     Global and local run QC - V1 - technical validation:
#         Input post - sequencing fasta files, variant text files and csv files
#         Output - cleaned fasta files, collated csv file with sample QC metrics
#         Source: 
#           initiate_script.py
#           post_QC_stats.py
# 2: TAB 2:
#     Lineage calling for each sample and control:
#         Input - run rack directory, rack ID
#         Output - result table rendered on screen and appended to the run rack 
#                 topsheet 
# 3: TAB 3:
#     V2 - scientific/clinical validation:
#         Input - sample metadata text file, run rack directory, rack ID
#         Output - multifasta file with fasta files and sequence IDs renamed 
#                 according to RSV nomenclature,
#                 appended xls template file for submission to GISAID
# 4: TAB 4:
#     GenBank submission preparation:
#         Input - sample fasta sequences
#         Output - respective annotated gbk files
#         Source:
#           annotate.py
#           defaults.py
################################################################################
#UKHSA, RVU, Tiina Talts v.01 24.08.2023
################################################################################



# Load Packages
library(shiny)
library(reticulate)
library(ranger)
library(kmer)
library(ape)
library(DT)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(shinyjs)
options(shiny.maxRequestSize = 5*1024^2)
library(rintrojs)
library(shinythemes)
library(magrittr)
library(dplyr)
library(plyr)
library(tibble)
library(parallel)
require(data.table)
library(readxl)
library(ShortRead)
library(tidyr)
library(tidyverse)
library("openxlsx")
library(RDCOMClient)
library(stringr)
library(ggpubr)
library(patchwork)
library(phylotools)
library(msa)
library(bios2mds)

use_python("C:/Users/Tiina.Talts/AppData/Local/r-miniconda/envs/r-reticulate/python.exe")
source_python("initiate_script.py")
ctrlMonitoring_data <- readRDS("data/ctrlMonitoring_data.rds")
connection = saveRDS(ctrlMonitoring_data, "data/temp.rds")
model <- readRDS("data/1.rds")
alignedB.df2 <- read.fasta("data/NC_001781_RSVB_resvidex_ref_G_ectodomain_nt.fas")
ref_fastaB <- readDNAStringSet("data/NC_001781_RSVB_resvidex_ref_G_ectodomain_nt.fas")
alignedA.df2 <- read.fasta("data/NC_038235_RSVA_resvidex_ref_G_ectodomain_nt.fas")
ref_fastaA <- readDNAStringSet("data/NC_038235_RSVA_resvidex_ref_G_ectodomain_nt.fas")


ui <- fluidPage(
  useShinyjs(),
  theme = shinytheme("sandstone"),
  tags$style(HTML( 
  '.navbar .navbar-nav {float: right;}
  ')),
  includeScript("www/script.js"),
  includeCSS("www/style.css"),
  navbarPage(title="RSV WGS post-sequencing & post-pipeline analysis",
    tabPanel("V1",
             fluidPage(
               sidebarLayout(
                 sidebarPanel(width = 2,
                   div(id="all",
                   textInput(inputId = "dataset", 
                             label = "Enter Folder Address:", 
                             width = 300),
                   textInput(inputId = "racknumber", 
                             label = "Enter Rack Number: EFAR_xxxx", 
                             width = 300),
                   hr(),
                   h5(HTML(paste0("<strong>","Enter controls<br/>
                            if ngs yield inadequate<br/>
                           in the following format: WTR_a_EFAR_xxxx<br/>
                           RSVAB_Bxxxxxx_xx_EFAR_xxxx<br/>
                           (Min 0 Max 5)","</strong>"))),
                   textInput(inputId = "ctrl1", label = "1", width = 300),
                   textInput(inputId = "ctrl2", label = "2", width = 300),
                   textInput(inputId = "ctrl3", label = "3", width = 300),
                   textInput(inputId = "ctrl4", label = "4", width = 300),
                   textInput(inputId = "ctrl5", label = "5", width = 300),
                   hr(),
                   submitButton("Submit for Processing"),
                   h3(textOutput("caption")),
                   hr(),
                   hr(),
                   hr(),
                   h5(textOutput("caption3")),
                   hr()
                 )),
                 mainPanel(width = 10,
                   div(id="all",
                   titlePanel("RSV WGS Technical Validation V1:"),
                   tags$style(type="text/css",
                              ".shiny-output-error { visibility: hidden; }",
                              ".shiny-output-error:before { visibility: hidden; }"),
                   fluidRow(
                     column(
                       width = 12,
                       align="center",
                       plotOutput("boxPlott")
                     ),
                     column(
                       width = 12,
                       align="left",
                       h3(textOutput("caption2"))
                     ),
                     column(
                       width = 12,
                       align="center",
                       DT::dataTableOutput("controls")
                     )
                   ),
                   HTML("<hr>"),
                   div(
                     id = "footer", align = "center",
                     "V.01 Developed and created by TT in 2023-06-26")
                 ))
               )
              )
    ),
    tabPanel("RSV LINEAGE",
             fluidPage(
               sidebarLayout(
                  sidebarPanel(width = 2,
                    div(id="all",
                        textInput(inputId = "text9", 
                                  label = "Enter Rack Folder Address:", 
                                  width = 300),
                        textInput(inputId = "user", 
                                  "Please enter Rack as EFAR_xxxx:"),
                        submitButton("RUN")
                    )),
                  mainPanel(width = 10,
                    div(id="all",
                    titlePanel("Identify RSV subtype, genotype, subgenotype and genetic lineage for each query sequence:"),
                    fluidRow(
                      DT::dataTableOutput("table")
                    ),
                  HTML("<hr>"),
                  div(
                    id = "footer",
                    "Developed and created by",
                    a("Marco Cacciabue", 
                      href = "https://sourceforge.net/u/marcocacciabue/profile/"),
                    " and Stephanie Goya",
                    a("Goya et al.", 
                      href = "https://onlinelibrary.wiley.com/doi/full/10.1111/irv.12715"),
                    "Implemented for UKHSA by TT in 2023")
                  )
                 )
               )
              )
             
    ),
 ###############################################################################
    tabPanel("V2",
             fluidPage(
               splitLayout(
                 wellPanel(style = "padding: 15px;background-color:#ffffff",
                 wellPanel(style = "padding: 30px;",
                    div(id="all",
                    titlePanel("PART I: Prepare .fasta and .xls files for GISAID batch upload"),
                    hr(),
                    fileInput(inputId = "v2file1", 
                              label = ("STEP 1: Input file: Text metadata file (FLRSSQ_.txt)"),
                               accept = ".txt"),
                    textInput(inputId = "v2text2", 
                              label = "STEP 2: Enter Rack Folder Address:", 
                              width = 1000),
                    hr(),
                    textInput(inputId = "v2text", 
                              label = "STEP 3: Enter GISAID user name:", 
                              width = 1000),
                    hr(),
                    submitButton("STEP 4: Submit for Processing"),
                    hr(),
                    tags$a(href="https://www.epicov.org/epi3/frontend#266eea", 
                           "Click here to open the GISAID Batch Upload Link in new window", 
                           target="_blank"),
                    hr(),
                    hr(),
                    hr(),
                    h3(textOutput("caption7"))
                 ))),
                 wellPanel(style = "padding: 15px;background-color:#ffffff",
                 wellPanel(style = "padding: 30px;",
                   div(id="all",
                   titlePanel(HTML("PART II: Update Topsheet with GISAID data<br/>and create MOLIS upload file")),
                   hr(),
                   fileInput(inputId = "v2file2", 
                             label = ("STEP 1: Input file: retrieved from GISAID \"Patient Status Metadata\" (.txt or .tsv)"),
                             accept = c(".txt",".tsv")),
                   fileInput(inputId = "v2file3", 
                             label = ("STEP 2: Input file: Topsheet_OUT_LIN_ (.xlsx)"),
                             accept = ".xlsx"),
                   textInput(inputId = "v2fileout", 
                             label = "STEP 3: Enter Rack Folder Address:", 
                             width = 1000),
                   textInput(inputId = "racknumber2", 
                             label = "STEP 4: Enter Rack Number: EFAR_xxxx", 
                             width = 300),
                   hr(),
                   submitButton("STEP 5: Submit for Processing"),
                   h3(textOutput("caption5"))
                 )))
               )
             )
    ),
    tabPanel("GENBANK",
             fluidPage(
               sidebarLayout(
                 sidebarPanel(
                 ),
                 mainPanel(
                   titlePanel("Create GenBank .gbk file")
                 )
               )
             )
    )
    )
)
################################################################################

# Define server logic 
server <- shinyServer(function(input,output,session){
  
 
################################################################################
  #runs the python with submit button:
  displayText1 = eventReactive(input$dataset, {
    x = input$dataset
    y = input$racknumber
    
    start_QCscript(x, y)
    
    
    ssttr = str_split_i(x, "\\\\|[^[:print:]]", -1)
    directory_path = dirname(x)
    
      
    path_csvfile = paste0(directory_path,
                          "\\",
                          "qc_out_", 
                          ssttr, 
                          "\\", 
                          "post_qc_", 
                          y, 
                          ".csv")
    
    df4 = as.data.frame(fread(path_csvfile))
    df.long <- df4 %>% 
      pivot_longer('sequence length':'filtered_mapped_reads', 
                   names_to = 'variable', 
                   values_to = 'value')
    
    negative_value <- filter(df.long, grepl("WTR", MOLIS))
    positive_value <- filter(df.long, grepl("RSVAB", MOLIS))
    
    ctrl = c(input$ctrl1, 
             input$ctrl2, 
             input$ctrl3, 
             input$ctrl4, 
             input$ctrl5)
    if(all(is.na(ctrl))){
      return (NULL)
    }else{
      ctrl0 = as.data.frame(ctrl)
      xx = nrow(ctrl0)
      newColNames4 = c("run_folder", 
                       "factor", 
                       "HQ_FMR")
      ctrl0[newColNames4] <- NA
      ctrl0[1:xx,"run_folder"] = ifelse(is.na(ctrl0$ctrl), 
                                        NA, 
                                        ssttr)
      ctrl0[1:xx,"factor"] = ifelse(grepl("WTR", 
                                          ctrl0$ctrl),
                                    "WTR",
                                    "RSVAB")
      ctrl0[1:xx,"HQ_FMR"] = 0
      ctrl0[!apply(ctrl0 == "", 1, all),]
      ctrl00 = relocate(ctrl0, run_folder, .before=ctrl)
      
    }
    
    
    #enter new data from the run just analysed:
    pos_neg = rbind(negative_value, positive_value)
    pos_neg0 = subset(pos_neg, variable == "filtered_mapped_reads")
    newColNames9 = c("run_folder", "factor")
    pos_neg0[newColNames9] <- NA
    x = nrow(pos_neg0)
    pos_neg0[1:x,"run_folder"] = ssttr
    pos_neg0[1:x,"factor"] = ifelse(grepl("WTR", 
                                          pos_neg0$MOLIS),
                                    "WTR",
                                    "RSVAB")
    pos_neg0$variable <- NULL
    pos_neg0 = relocate(pos_neg0, run_folder, .before=MOLIS)
    pos_neg0 = relocate(pos_neg0, factor, .before=value)
    ctrlMonitoring_data = rbind(ctrlMonitoring_data, 
                                setNames(pos_neg0, 
                                         names(ctrlMonitoring_data)))
    
    if(all(is.na(ctrl))){
      return (NULL)
    }else{
      ctrlMonitoring_data %>% rbind(., ctrl00)
    }

    connection
    
    ctrlMonitoring_dataSub <- subset(ctrlMonitoring_data, 
                                     select = c(factor, HQ_FMR))
    
    ctrl_long <- melt(setDT(ctrlMonitoring_dataSub), 
                      id.vars = "factor", 
                      variable.name = "HQ_FMR")
    
    x_wtr = setDT(ctrl_long)[factor == "WTR", as.numeric(value)]
    x_pos = setDT(ctrl_long)[factor == "RSVAB", as.numeric(value)]
    
    upper_limitneg <- median(x_wtr) + 3 * mad(x_wtr, constant = 1)
    lower_limitpos <- median(x_pos) - 3 * mad(x_pos, constant = 1)
    upper_limitpos <- median(x_pos) + 3 * mad(x_pos, constant = 1)
    outs_pos <- x_pos[x_pos < lower_limitpos | x_pos > upper_limitpos]
    upper_negs <- x_wtr[x_wtr > (median(x_wtr) + 3 * mad(x_wtr, constant = 1))]
    lower_poss <- x_pos[x_pos < (median(x_pos) - 3 * mad(x_pos, constant = 1))]
    outliersRemoved_wtr1 = x_wtr[!x_wtr %in% upper_negs]
    outliersRemoved_pos1 = x_pos[!x_pos %in% lower_poss]
    water_median <- median(as.numeric(outliersRemoved_wtr1))
    pos_median <- median(as.numeric(outliersRemoved_pos1), na.rm = T)
    mad_wtr = mad(outliersRemoved_wtr1, water_median, constant = 1, na.rm = T)
    mad_pos = mad(outliersRemoved_pos1, pos_median, constant = 1, na.rm = T)
    ctrlMonitoring_data1 = subset(ctrlMonitoring_data, run_folder == ssttr)
    newColNames7 = c("Control:", 
                     "Outlier_flag:", 
                     "Median_total:", 
                     "MAD_total:", 
                     "Hampel_filter:", 
                     "PASS/FAIL:")
    ctrlMonitoring_data1[newColNames7] <- NA
    ctrlMonitoring_data1$'run_folder' <- NULL
    xxx = nrow(ctrlMonitoring_data1)
    
    ctrlMonitoring_data1 <- ctrlMonitoring_data1 %>%
      mutate("Rack:" = str_extract(ctrl, "(EFAR_[0-9]{4})"))
    
    ctrlMonitoring_data1 = relocate(ctrlMonitoring_data1, 
                                    "Rack:", 
                                    .before="Outlier_flag:")
    ctrlMonitoring_data1[1:xxx,"Control:"] = fifelse(grepl("WTR", 
                                                           ctrlMonitoring_data1$factor),
                                                     "NEG",
                                                     "POS")
    ctrlMonitoring_data1$factor <- NULL
    ctrlMonitoring_data1[1:xxx,"Outlier_flag:"] = fifelse(grepl("POS", 
                                                                ctrlMonitoring_data1$"Control:") & 
                                                            ctrlMonitoring_data1$HQ_FMR %in% outs_pos | 
                                                            grepl("NEG", ctrlMonitoring_data1$"Control:") & 
                                                            ctrlMonitoring_data1$HQ_FMR %in% upper_negs,
                                                           "yes", "no")
    
    ctrlMonitoring_data1[1:xxx,"Median_total:"] = fifelse(grepl("POS", 
                                                                ctrlMonitoring_data1$"Control:"), 
                                                          pos_median, 
                                                          water_median)
    ctrlMonitoring_data1[1:xxx,"MAD_total:"] = fifelse(grepl("POS", 
                                                             ctrlMonitoring_data1$"Control:"), 
                                                       mad_pos, 
                                                       mad_wtr)
    ctrlMonitoring_data1[1:xxx,"Hampel_filter:"] = fifelse(grepl("POS", 
                                                                 ctrlMonitoring_data1$"Control:"), 
                                                           lower_limitpos, 
                                                           upper_limitneg)
    ctrlMonitoring_data1[1:xxx,"PASS/FAIL:"] = fifelse(grepl("POS", 
                                                             ctrlMonitoring_data1$"Control:") & 
                                                         ctrlMonitoring_data1$HQ_FMR %in% lower_poss & 
                                                         ctrlMonitoring_data1$HQ_FMR < pos_median | 
                                                         grepl("NEG", ctrlMonitoring_data1$"Control:") & 
                                                         ctrlMonitoring_data1$HQ_FMR %in% upper_negs & 
                                                         ctrlMonitoring_data1$HQ_FMR > 1000, 
                                                       "FAIL", "PASS")
    ctrlMonitoring_data2<-ctrlMonitoring_data1[which(ctrlMonitoring_data1$ctrl > ctrlMonitoring_data1$HQ_FMR), ]
    generated_datatbl <- datatable(
        ctrlMonitoring_data2, 
        options = list(
          columnDefs = list(list(className = 'dt-center', targets = 2:8))
          )
        ) %>%
          formatStyle("PASS/FAIL:", color = styleEqual("FAIL", "red")) %>%
          formatStyle("PASS/FAIL:", color = styleEqual("PASS", "green"))

    output$controls <- renderDataTable({
      generated_datatbl
      
    })
    
    g <- ggplot(data = df.long, aes(x = variable, y = value)) +
          geom_boxplot() +
          facet_wrap(facets = ~variable, scales = 'free') +
          geom_dotplot(binaxis="y", stackdir="center", dotsize=0.6) +
          geom_point(data = positive_value, color = "red", size=3) +
          geom_point(data = negative_value, color = "blue", size=3) +
          geom_text(aes(label = MOLIS), 
                    data = positive_value, 
                    color = "red", 
                    vjust = -1) +
          geom_text(aes(label = MOLIS), 
                    data = negative_value, 
                    color = "blue", 
                    vjust = -1)
    
    output$boxPlott <- renderPlot({
      g
    })
    
    output$caption2 <- renderText({
      paste(paste0("FTP_", ssttr, ":"))
    })
    output$caption3 <- renderText({
      now<-format(Sys.time(), "%Y-%m-%d-%H%M%S")
      paste(paste0("Date/Time: ", now))
    })
    
    
  })
  output$caption <- renderText({
    displayText1()
  })
  
  
  
  ##############################################################################
  #Lineage
  ##############################################################################
  
  
  data_reactive<- eventReactive(input$user,{
    directory_path1 = input$text9
    path_tofastafiles = paste0(directory_path1, "\\", "8 - Fastas_QC_PASS")
    
    pattB <- ".+_b.+fas$"
    pattA <- ".+_a.+fas$"
    
    if(length(list.files(path = path_tofastafiles, patt=pattB)) == 0){
      alignedB.df2 <- read.fasta("data/NC_001781_RSVB_resvidex_ref_G_ectodomain_nt.fas")
    }else{
      fastasB <- readFasta(paste0(path_tofastafiles, "\\"), pattB)
      ref_fastaB <- readDNAStringSet("data/NC_001781_RSVB_resvidex_ref_G_ectodomain_nt.fas")
      
      tmpdir_locationB = "temp/combinedB.fas"
      writeFasta(fastasB, tmpdir_locationB)
      writeXStringSet(ref_fastaB, tmpdir_locationB, append = T)
      
      mySequencesB <- readDNAStringSet("temp/combinedB.fas")
      alignmentB <- msaMuscle(mySequencesB, 
                              type = "dna", 
                              gapOpening = 800, 
                              gapExtension = 0.2)
      
      alignCW_as_alignB <- msaConvert(alignmentB, "seqinr::alignment")
      alignCW_as_alignBB <- msaConvert(alignmentB, "bios2mds::align")
      export.fasta(alignCW_as_alignBB, 
                   outfile = "temp/alignmentB.fas", 
                   ncol = 60, 
                   open = "w")
      
      alignedB.df <- data.frame(seq.name = alignCW_as_alignB$nam, 
                                seq.text = alignCW_as_alignB$seq)
      
      num_ntstartB <- nchar(sapply(strsplit(grep("([-]+(?i)CCAATCC(?-i))", 
                                                 alignedB.df$seq.text[grepl("NC_001781", 
                                                                            alignedB.df$seq.name)], 
                                                 value = TRUE), 
                                            "(?i)CCAATCC(?-i)"), "[", 1))
      num_nttailB <- nchar(sapply(strsplit(grep("((?i)AAATTCC(?-i)[-]+)", 
                                                alignedB.df$seq.text[grepl("NC_001781", 
                                                                           alignedB.df$seq.name)], 
                                                value = TRUE), "(?i)AAATTCC(?-i)"), "[", 2))
      num_nttotalB <- nchar(alignedB.df$seq.text[grepl("NC_001781", alignedB.df$seq.name)])
      num_ntendB <- num_ntstartB + (num_nttotalB - num_ntstartB - num_nttailB)
      alignedB.df2 <- mutate(alignedB.df, 
                             seq.text = substring(alignedB.df$seq.text, 
                                                  num_ntstartB+1, 
                                                  num_ntendB))
    }
    
    if(length(list.files(path = path_tofastafiles, patt=pattA)) == 0){
      alignedA.df2 <- read.fasta("data/NC_038235_RSVA_resvidex_ref_G_ectodomain_nt.fas")
    }else{
      fastasA <- readFasta(paste0(path_tofastafiles, "\\"), pattA)
      ref_fastaA <- readDNAStringSet("data/NC_038235_RSVA_resvidex_ref_G_ectodomain_nt.fas")
      
      tmpdir_locationA = "temp/combinedA.fas"
      writeFasta(fastasA, tmpdir_locationA)
      writeXStringSet(ref_fastaA, tmpdir_locationA, append = T)
      
      mySequencesA <- readDNAStringSet("temp/combinedA.fas")
      alignmentA <- msaMuscle(mySequencesA, 
                              type = "dna", 
                              gapOpening = 800, 
                              gapExtension = 0.2)
      
      alignCW_as_alignA <- msaConvert(alignmentA, "seqinr::alignment")
      alignCW_as_alignAA <- msaConvert(alignmentA, "bios2mds::align")
      export.fasta(alignCW_as_alignAA, 
                   outfile = "temp/alignmentA.fas", 
                   ncol = 60, open = "w")
      
      alignedA.df <- data.frame(seq.name = alignCW_as_alignA$nam, 
                                seq.text = alignCW_as_alignA$seq)
      
      num_ntstartA <- nchar(sapply(strsplit(grep("([-]+(?i)GTCCCT(?-i))", 
                                                 alignedA.df$seq.text[grepl("NC_038235", 
                                                                            alignedA.df$seq.name)], 
                                                 value = TRUE), "(?i)GTCCCT(?-i)"), "[", 1))
      num_nttailA <- nchar(sapply(strsplit(grep("((?i)CACGCC(?-i)[-]+)", 
                                                alignedA.df$seq.text[grepl("NC_038235", 
                                                                           alignedA.df$seq.name)], 
                                                value = TRUE), "(?i)CACGCC(?-i)"), "[", 2))
      num_nttotalA <- nchar(alignedA.df$seq.text[grepl("NC_038235", alignedA.df$seq.name)])
      num_ntendA <- num_ntstartA + (num_nttotalA - num_ntstartA - num_nttailA)
      alignedA.df2 <- mutate(alignedA.df, 
                             seq.text = substring(alignedA.df$seq.text, 
                                                  num_ntstartA+1, 
                                                  num_ntendA))
    }
    
    comb_aligned_extr_df = rbind(alignedB.df2, alignedA.df2)
    y <- t(sapply(strsplit(comb_aligned_extr_df[,2],""), tolower))
    names(y) <- comb_aligned_extr_df[,1]
    
    model <- readRDS("data/1.rds")
    
    query_count <- kcount(as.DNAbin(y), k = model$kmer)
    genome_length <- 0
    n_length <- 0
    for(i in 1:length(query_count[,1])){
      k <- as.DNAbin(y)[i]
      k <- as.matrix(k)
      query_count[i,] <- query_count[i,]*model$kmer/(length(k))
      genome_length[i] <- length(k)
      n_length[i] <- round(100*base.freq(k, all = TRUE)[15],2)
    }
    
    calling <- predict(model, query_count)
    #Run the predict method from de Ranger package, retaining the 
    #classification result from each tree in the model (to calculate 
    #a probability value for each classification)
    calling_all <- predict(model, query_count, predict.all = TRUE)
    probability <- rep(0, length(calling_all$predictions[,1]))
    
    for (i in 1:length(calling_all$predictions[,1])) {
      #extract predictions for each query sample in temp vector,
      #count the number of correct predictions and divide by number 
      #of trees to get a probability.
      temp <- calling_all$predictions[i,]
      probability[i] <- sum(temp == which(model$forest$levels == calling$predictions[i]))/model$num.trees
      
    }
    
    N_QC <- (n_length < 0.5)
    Length_QC <- (genome_length > 590) & (genome_length < 700)
    temp1 <- strsplit(as.character(calling$prediction), split="_")
    Subtype = unlist(lapply(temp1, function(l) l[[1]]))
    Clade = paste(Subtype,".",unlist(lapply(temp1, function(l) l[[2]])),sep="")
    temp1 <- strsplit(as.character(Subtype), split="G")
    Subtype = unlist(lapply(temp1, function(l) l[[2]]))
    
    data_out <- data.frame(Label = row.names(query_count), 
                           Subtype = Subtype,
                           Clade = Clade, 
                           Probability = probability,
                           Length = genome_length,
                           Length_QC = Length_QC,
                           N = n_length,
                           N_QC = N_QC)
    
    
    data_out = data_out %>%
      mutate(lin_prob = ifelse( (probability > .8), "confirmed", "probable")) %>%
      unite('Clade', c('Clade', 'lin_prob'), sep=' : ', remove=T)
    
    
    path_topsheetfile = paste0(directory_path1, 
                               "\\", 
                               "10 - Topsheet_MOLIS", 
                               "\\", 
                               "Topsheet_", 
                               input$user, 
                               ".xlsx")
    
    wb3 <- loadWorkbook(file = path_topsheetfile, isUnzipped=F)
    
    sheet8 <- openxlsx::addWorksheet(wb3, "LINEAGE_data")
    textstyle <- openxlsx::createStyle(fontName = "Arial", 
                                       fontSize = 7, 
                                       numFmt = "@", 
                                       halign = "center", 
                                       valign = "center"
    )
    
    writeData(wb3, "LINEAGE_data", data_out)
    
    textcells <- expand.grid(row = 1:200, col = 1:40)
    openxlsx::addStyle(wb3, sheet = "LINEAGE_data", 
                       rows = textcells$row, cols = textcells$col, 
                       style = textstyle)
    
    formula3 <- paste(paste(paste0(paste0("=VLOOKUP(B" , 
                                          seq(2,nrow(data_out)+2 ,1), sep=','), 
                                   ' LINEAGE_data!A:H'),
                            ' "3" ', sep=','),' "FALSE") ', sep=',')
    
    writeFormula(wb3, 
                 sheet = "Topsheet", 
                 x = formula3, 
                 startCol = 23, 
                 startRow = 2, 
                 array = FALSE)
    
    now<-format(Sys.time(), "%Y-%m-%d-%H%M%S")
    outputfile3 <- paste0(directory_path1, 
                          "\\", 
                          "10 - Topsheet_MOLIS", 
                          "\\Topsheet_OUT_LIN_", 
                          input$user, 
                          "_", 
                          now, 
                          ".xlsx")
    
    saveWorkbook(wb3, outputfile3, overwrite = TRUE)
    
    xlApp <- COMCreate("Excel.Application")
    
    #open the same file that you saved earlier with saveWorkbook() method
    xlWbk <- xlApp$Workbooks()$Open(outputfile3)
    
    #save the file and close it
    xlWbk$Save()
    xlWbk$Close()
    
    #quit excel instance    
    xlApp$Quit()
    rm(list= c("xlApp", "xlWbk"))
    unlink("temp/*")
    
    
    list(message = "Done!", data_out = data_out)
  })
  output$text <- renderText({
    data_reactive()$message
  })
  
  table <- reactive({
    data_out <- data_reactive()$data_out
    data_out
  })
  
  output$table <- DT::renderDataTable({
    
    col <- brewer.pal(5, "Blues")
    col2 <- brewer.pal(5, "Reds")
    
    table <- table()
    datatable(table, selection = 'single',
              options = list(
                columnDefs = list(list(className = 'dt-center', targets = 2:8))),
    ) %>%
     formatStyle("Length", "Length_QC",
                  backgroundColor = styleEqual(c(0, 1), c(col[3], col[1]))) %>% 
     formatStyle("N", "N_QC", 
                 backgroundColor = styleEqual(c(0, 1), 
                                              c(col[3], col[1]))) 
    
    
    
  })
  
  
  count_parallel <- function(x, kmer) ({
    library(kmer)
    kcount(x, k = kmer)
  })
  

  
  
  ############################################################################## 
  #V2 PART I:
  ##############################################################################
  
  
  displayText2 = eventReactive(input$v2text, {
    req(input$v2file1)
    metadata = as.data.frame(fread(input=input$v2file1$datapath,encoding = "UTF-8"))
    metadata2 = metadata[!is.na(metadata$'SAMPLE_DT'), ]
    newColNames2 <- c("a", "b")
    metadata3 <- separate(metadata2, MOLIS_number, newColNames2, sep=3, remove=T)
    newColNames3 <- c("c", "d")
    metadata4 <- separate(metadata3, SAMPLE_DT, newColNames3, sep=4, remove=F)
    metadata4[,10] <- toupper(metadata4[,10])
    metadata4[,10] = gsub('#', '', metadata4[,10])
    metadata5 = add_column(metadata4, c = 'hRSV', .after = 'b')
    metadata6 = add_column(metadata5, d = 'England', .after = 'c')
    metadata7 = metadata6 %>% unite(Virus_name, 
                                    c(c.1, FLRSAB, d.1, b, c), 
                                    sep='/', remove=F)
    newColNames4 <- c("e", "f")
    metadata8 <- separate(metadata7, ID, newColNames4, sep='_', remove=T)
    metadata8[,2] = gsub('OM', 'Original', metadata8[,2])
    x = nrow(metadata8)
    metadata9 <- subset(metadata8, select = c(f, Virus_name, FLRSAB, SAMPLE_DT))
    metadata99 = mutate(metadata9, SAMPLE_DT = as.character(SAMPLE_DT), 
                        SAMPLE_DT = as.character(paste(SAMPLE_DT)))
    metadata10 = rename(metadata99, 
                        c('f'='rsv_passage', 
                          'Virus_name'='rsv_virus_name', 
                          'FLRSAB'='rsv_subtype', 
                          'SAMPLE_DT'='rsv_collection_date'))
    newColNames5 <- c("rsv_location", 
                      'submitter',
                      'fn',
                      "rsv_add_location", 
                      'rsv_host', 
                      'rsv_add_host_info',
                      'rsv_sampling_strategy',
                      'rsv_gender',
                      'rsv_patient_age',
                      'rsv_patient_status',
                      'rsv_specimen',
                      'rsv_outbreak',
                      'rsv_last_vaccinated',
                      'rsv_treatment',
                      'rsv_seq_technology',
                      'rsv_assembly_method',
                      'rsv_coverage',
                      'rsv_orig_lab',
                      'rsv_orig_lab_addr',
                      'rsv_provider_sample_id',
                      'rsv_subm_lab',
                      'rsv_subm_lab_addr',
                      'rsv_subm_sample_id',
                      'rsv_authors',
                      'rsv_comment',
                      'comment_type'
    )
    
    directory_path3 = input$v2text2
    path_tofastafiles = paste0(directory_path3, "\\", "8 - Fastas_QC_PASS")
    now<-format(Sys.time(), "%Y-%m-%d-%H%M%S")
    patt <- "([HRE]{1,2}[0-9]{8,9}_OM_[ab]{1}).fas"
    fasta <- readFasta(paste0(path_tofastafiles, "\\"), patt)
    dir_location = "temp/combined.fasta"
    writeFasta(fasta, dir_location)
    
    fastas = readBStringSet(dir_location)
    names(fastas)=str_split_fixed(names(fastas),"\\_",8)[,1]
    names(fastas)=ifelse(names(fastas) %in% metadata8$e, 
                         metadata8$Virus_name, 'FALSE')
    outputfile6 <- paste0(directory_path3, 
                          "\\", "9 - GISAID", 
                          "\\", now, 
                          "_EpiRSV_BulkUpload.fasta")
    writeXStringSet(fastas, outputfile6, format="fasta")
    
    metadata10[newColNames5] <- NA
    
    metadata10[1:x,"submitter"] = input$v2text
    metadata10[1:x,"fn"] = paste0(now, "_EpiRSV_BulkUpload.fasta")
    metadata10[1:x,"rsv_location"] = "Europe / United Kingdom / England"
    metadata10[1:x,"rsv_host"] = "Human"
    metadata10[1:x,"rsv_gender"] = "unknown"
    metadata10[1:x,"rsv_patient_age"] = "unknown"
    metadata10[1:x,"rsv_patient_status"] = "unknown"
    metadata10[1:x,"rsv_seq_technology"] = "Illumina NextSeq"
    metadata10[1:x,"rsv_assembly_method"] = "UKHSA pipeline (BWA-SAMtools-Quasibam)"
    metadata10[1:x,"rsv_orig_lab"] = "Respiratory Virus Unit / Reference Microbiology Services / UK Health Security Agency"
    metadata10[1:x,"rsv_orig_lab_addr"] = "61 Colindale Avenue London NW9 5EQ United Kingdom"
    metadata10[1:x,"rsv_subm_lab"] = "Reference Microbiology Services / UK Health Security Agency"
    metadata10[1:x,"rsv_subm_lab_addr"] = "61 Colindale Avenue London NW9 5EQ United Kingdom"
    metadata10[1:x,"rsv_authors"] = "Zambon M. Talts T. Kele B. Miah S"
    
    metadata11 = relocate(metadata10, submitter, .before=rsv_passage)
    metadata11 = relocate(metadata11, fn, .after=submitter)
    metadata11 = relocate(metadata11, rsv_virus_name, .after=fn)
    metadata11 = relocate(metadata11, rsv_subtype, .after=rsv_virus_name)
    
    wb2 <- loadWorkbook(file = "data\\20210611_EpiRSV_BulkUpload_Template.xlsx", 
                        isUnzipped=F)
    
    
    
    textstyle <- openxlsx::createStyle(numFmt = "@")
    textcells <- expand.grid(row = 1:200, col = 6)
    openxlsx::addStyle(wb2, sheet = "Submissions", 
                       rows = textcells$row, cols = textcells$col, 
                       style = textstyle) 
    
    writeData(wb2,
              "Submissions",
              metadata11,
              colNames = FALSE,
              startRow = 3)
    
    outputfile <- paste0(directory_path3, 
                         "\\", 
                         "9 - GISAID", 
                         "\\", 
                         now, 
                         "_EpiRSV_BulkUpload.xlsx")
    saveWorkbook(wb2, outputfile, overwrite = TRUE)
    path_Excel_File1 <- outputfile
    path_Excel_File_Output <- paste0(directory_path3, 
                                     "\\", 
                                     "9 - GISAID", 
                                     "\\", 
                                     now, 
                                     "_EpiRSV_BulkUpload.xls")
    xlApp <- COMCreate("Excel.Application")
    xlWbk1 <- xlApp$Workbooks()$Open(path_Excel_File1)
    xlWbk1$SaveAs(path_Excel_File_Output, -4143) # saving as .xls
    xlApp$Quit()
    rm(list= c("xlApp", "xlWbk1"))
    unlink(path_Excel_File1)
    unlink(dir_location)
    paste("Completed! Files available in \"GISAID_uploads\" folder.")
  })
  output$caption7 <- renderText({
    displayText2()
  })
  
  
  
  ############################################################################## 
  #V2 PART II:
  ##############################################################################
  
  
  displayText5 = eventReactive(input$racknumber2, {
    req(input$v2file2)
    
    data1<-as.data.frame(fread(input=input$v2file2$datapath, 
                               encoding = "UTF-8"))
    
    req(input$v2file3)
    data2 = as.data.frame(read_excel(path = input$v2file3$datapath, 
                                     sheet = "Topsheet"))
    data22 = data2[!is.na(data2$'Col5'), ]
    newColNames <- c("a", "b", "c", "d", "e")
    data11 <- separate(data1, `Virus name`, newColNames, sep="/", remove=T)
    data33 = data22 %>% full_join(data11, by = join_by('Col5' == 'd'), keep=T)
    df = data33 %>% unite(Virus_name, c(a, b, c, d, e), sep='/', remove=T)
    df$"GISAID virus name:" <- NULL
    df$"GISAID ID:" <- NULL
    wb <- loadWorkbook(file = input$v2file3$datapath, isUnzipped=F)
    
    sheet <- openxlsx::addWorksheet(wb, "GISAID_data")
    textstyle <- openxlsx::createStyle(fontName = "Arial", 
                                       fontSize = 7, 
                                       numFmt = "@", 
                                       halign = "center", 
                                       valign = "center"
    )
  
    writeData(wb, "GISAID_data", df)
    
    textcells <- expand.grid(row = 1:200, col = 1:40)
    openxlsx::addStyle(wb, sheet = "GISAID_data", 
                       rows = textcells$row, cols = textcells$col, 
                       style = textstyle)
    
    formula1 <- paste(paste(paste0(paste0('=VLOOKUP(D' ,
                                          seq(6,nrow(df)+5 ,1), sep=','), 
                                   ' GISAID_data!D:X'),' "21" ', sep=','),
                      ' "FALSE") ', sep=',')
    formula2 <- paste(paste(paste0(paste0('=VLOOKUP(D' ,
                                          seq(6,nrow(df)+5 ,1), sep=','), 
                                   ' GISAID_data!D:Z'),' "23" ', sep=','),
                      ' "FALSE") ', sep=',')
    
    writeFormula(wb, 
                 sheet = "Topsheet", 
                 x = formula1, 
                 startCol = 24, 
                 startRow = 6)
    writeFormula(wb, 
                 sheet = "Topsheet", 
                 x = formula2, 
                 startCol = 25, 
                 startRow = 6)
    
    now<-format(Sys.time(), "%Y-%m-%d-%H%M%S")
    outputfile2 <- paste0(input$v2fileout, 
                          "\\", 
                          "10 - Topsheet_MOLIS", 
                          "\\Topsheet_OUT_GIS_", 
                          input$racknumber2, 
                          "_", 
                          now, 
                          ".xlsx")
    
    saveWorkbook(wb, outputfile2, overwrite = TRUE)
    
    
  
    ############################################################################
    #text MOLIS out
    ############################################################################
    
    
    
    df_molis = df %>% select('Col6', 'Virus_name', 'Accession ID', "Lineage")
    df_molis2 = df_molis[!is.na(df_molis$'Accession ID'), ]
    df_molis3 = df_molis2 %>% unite('RES{FLRSSQ}', 
                                    c('Accession ID', 
                                      'Virus_name', 
                                      "Lineage"), 
                                    sep=' ; ', 
                                    remove=T)
    df_molis4 = rename(df_molis3, 'Col6'='SAMPLEID')
    outputfile3 <- paste0(input$v2fileout, 
                          "\\", 
                          "10 - Topsheet_MOLIS", 
                          "\\MOLIS_OUT_", 
                          input$racknumber2, 
                          "_", 
                          now, 
                          ".txt")
    write_tsv(df_molis4, outputfile3)
    
    paste("Completed! Files available in \"Topsheet\" folder.")
  
    
  
  })
  output$caption5 <- renderText({
    displayText5()
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
