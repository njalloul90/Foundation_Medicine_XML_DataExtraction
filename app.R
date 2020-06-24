# This is an application to extract data from FoundationMedicne xml files
# author: nahed jalloul
# date: 05/2020
library(shiny)
library(XML)
library("methods")
library(DT)
library(stringr)
library(readr)
library(shinythemes)
library(crosstalk)
library(leaflet)
library(patchwork)
library(scales)
library(tidyverse)
library(shinypanels)
library(markdown)
############################################################################################################
# Read foundation xml function
foundationXML_general <- function(inputFile){
    
    fileData <- xmlParse(inputFile)
    
    # Read the text from the file
    xmlText <- paste(readLines(inputFile), "\n", collapse="")
    # Extract General Information
    ReportID <- gsub('^.*<ReportId>\\s*|\\s*</ReportId>.*$', '', xmlText)
    SubmittedDiagnosis <- gsub('^.*<SubmittedDiagnosis>\\s*|\\s*</SubmittedDiagnosis>.*$', '', xmlText)
    Gender <- gsub('^.*<Gender>\\s*|\\s*</Gender>.*$', '', xmlText)
    DOB <- gsub('^.*<DOB>\\s*|\\s*</DOB>.*$', '', xmlText)
    
    if (str_length(DOB)> 11) {
        DOB <- gsub('^.*<DOB xmlns:xsd="http://www.w3.org/2001/XMLSchema">\\s*|\\s*</DOB>.*$', '', xmlText)
    }
    
    SpecSite <- gsub('^.*<SpecSite>\\s*|\\s*</SpecSite>.*$', '', xmlText)
    if (str_length(SpecSite) > 50) {
        SpecSite <- gsub('^.*<SpecFormat>\\s*|\\s*</SpecFormat>.*$', '', xmlText)
    }
    
    CollDate <- gsub('^.*<CollDate>\\s*|\\s*</CollDate>.*$', '', xmlText)
    if (str_length(CollDate) > 11) {
        CollDate <- gsub('^.*<CollDate xmlns:xsd="http://www.w3.org/2001/XMLSchema">\\s*|\\s*</CollDate>.*$', '', xmlText)
    }
    
    ReceivedDate <- gsub('^.*<ReceivedDate>\\s*|\\s*</ReceivedDate>.*$', '', xmlText)
    if (str_length(ReceivedDate) > 11) {
        ReceivedDate <- gsub('^.*<ReceivedDate xmlns:xsd="http://www.w3.org/2001/XMLSchema">\\s*|\\s*</ReceivedDate>.*$', '', xmlText)
    }
    CountryOfOrigin <- gsub('^.*<CountryOfOrigin>\\s*|\\s*</CountryOfOrigin>.*$', '', xmlText)
    if (str_length(CountryOfOrigin) > 50) {
        CountryOfOrigin <- c("NA")
    }
    # Purity estimates
    Pathological_Purity <- as.numeric(gsub('^.*percent-tumor-nuclei="\\s*|\\s*".*$', '', xmlText))
    Computational_Purity <- as.numeric(gsub('^.*purity-assessment="\\s*|\\s*".*$', '', xmlText))
    # Biomarkers
    TMB_Score <- as.numeric(gsub('^.*tumor-mutation-burden score="\\s*|\\s*".*$', '', xmlText))
    MSS <- gsub('^.*microsatellite-instability status="\\s*|\\s*".*$', '', xmlText)
    if (str_length(MSS)>9) {
        MSS = ""
    }
    
    # Short Variants
    ShortVariants <- gsub('^.*<short-variants>\\s*|\\s*</short-variants>.*$', '', xmlText)
    ShortVariants <- str_split(ShortVariants, " </short-variant> ")
    ShortVariants<- ShortVariants[[1]][] # all short variants with attributes
    Genes <- c()
    Allele_Freqs <- c()
    Depths <- c()
    Protein_Change <- c()
    Positions <- c()
    Functional_Effect <- c()
    Mutations <- c()
    Strand <- c()
    for (i in 1: length(ShortVariants)) {
        line <- ShortVariants[i]
        geneName <- gsub('^.*gene=\"\\s*|\\s*\".*$', '', line)
        Genes[i] <- geneName
        alleleFreq <- gsub('^.*allele-fraction="\\s*|\\s*".*$', '', line)
        Allele_Freqs[i] <- as.numeric(alleleFreq)
        depth <- gsub('^.*depth="\\s*|\\s*".*$', '', line)
        Depths[i] <- as.numeric(depth)
        proteinchange <- gsub('^.*protein-effect="\\s*|\\s*".*$', '', line)
        Protein_Change[i] <- proteinchange
        position <- gsub('^.*position="\\s*|\\s*".*$', '', line)
        Positions[i] <- position
        funct_effect <- gsub('^.*functional-effect="\\s*|\\s*".*$', '', line)
        Functional_Effect[i] <- funct_effect
        mutation <- gsub('^.*cds-effect="\\s*|\\s*".*$', '', line)
        Mutations[i] <- mutation
        strand <- gsub('^.*strand="\\s*|\\s*".*$', '', line)
        Strand[i] <- strand
    }
    
    # find short variants that overlap with copy-number-alterations
    # Copy Number Alterations
    CopyNumberAlterations <- gsub('^.*<copy-number-alterations>\\s*|\\s*</copy-number-alterations>.*$', '', xmlText)
    CopyNumberAlterations <- str_split(CopyNumberAlterations, " </copy-number-alteration> ")
    CopyNumberAlterations<- CopyNumberAlterations[[1]][] # all CNAs with attributes
    CNA_geneName <- c()
    CNA_copyNumber <- c()
    CNA_position <- c()
    for (i in 1 : length(CopyNumberAlterations)) {
        line <- CopyNumberAlterations[i]
        cna_geneName <- gsub('^.*gene=\"\\s*|\\s*\".*$', '', line)
        
        cna_copynumber <- gsub('^.*copy-number=\"\\s*|\\s*".*$', '', line)
        
        cna_pos <- gsub('^.*position=\"\\s*|\\s*\".*$', '', line)
        
        if (str_length(cna_copynumber)>5) {
            cna_geneName = ""
            cna_copynumber = ""
            cna_pos = ""
        }
        CNA_geneName[i] <- cna_geneName
        CNA_copyNumber[i] <- cna_copynumber
        CNA_position[i] <- cna_pos
    }
    
    
    var_genes <- data.frame(Genes,Positions)
    cna_gene<- data.frame(CNA_geneName,CNA_copyNumber,CNA_position)
    Ploidy <- rep(2,length(var_genes$Genes)) # initialize ploidy
    
    common <- intersect(var_genes$Genes, cna_gene$CNA_geneName) # find common Vars & CNAs
    position_var_gene <- var_genes[Genes==common,]
    pos_cna_gene <- cna_gene[CNA_geneName==common,]
    Genes_CN = data.frame(Genes,Ploidy)
    
    if (dim(position_var_gene[1]) > 0) {
        # if common variants exist
        Genes_CN = data.frame(Genes,Ploidy)
        Genes_CN[Genes_CN$Genes == common, 2]
        cn <- cna_gene[CNA_geneName == common, 2]
        Genes_CN[Genes_CN$Genes == common, 2] = as.numeric(as.character(cn))
    }
    
    
    
    # create temporary txt file with Gene AlleleFreq Depth and Ploidy for R_allfit.py
    # headers SNV	Allele_freq	Depth	Ploidy PathologicalPurity ComputationalPurity
    PathologicalPurity = rep(Pathological_Purity/100,length(Genes))
    ComputationalPurity = rep(Computational_Purity/100,length(Genes))
    temp_df <- data.frame(Genes_CN$Genes,100*Allele_Freqs,Depths,Genes_CN$Ploidy,PathologicalPurity,
                          ComputationalPurity) # allele freq in %
    colnames(temp_df) <- c("ID","Allele_freq","Depth","Ploidy","PathologicalPurity","ComputationalPurity")
    write.table(temp_df, file = "tabledata.txt", sep = "\t", dec = ".",
                row.names = FALSE, col.names = TRUE)
    
    #pur <- py_run_file("All-FIT2.py")
    
    # access output file
    #mystring <- read_file("tabledata_out.txt")
    
    #allfit_pur <- as.numeric(gsub('^.*tabledata_out\t\\s*|\\s*\t.*$', '', mystring))
    #allfit_pur_CI <- gsub('^.*tabledata_out\t*\t\\s*|\\s*".*$', '', mystring)
    #allfit_pur_CI <- gsub('^.*\t\\s*|\\s*\n.*$', '', allfit_pur_CI)
    #allfit_pur <- rep(100*allfit_pur,length(Genes))
    
    # Define Dataframe: 
    df_GeneralInfo <- data.frame(rep(ReportID,length(ShortVariants)),rep(SubmittedDiagnosis,length(ShortVariants)),
                                 rep(Gender, length(ShortVariants)),rep(DOB,length(ShortVariants)),
                                 rep(SpecSite,length(ShortVariants)),rep(CollDate,length(ShortVariants)),
                                 rep(ReceivedDate,length(ShortVariants)),rep(CountryOfOrigin,length(ShortVariants)),
                                 rep(Pathological_Purity,length(ShortVariants)),rep(Computational_Purity,length(ShortVariants)),
                                 rep(TMB_Score,length(ShortVariants)),rep(MSS,length(ShortVariants)))
    colnames(df_GeneralInfo) <- c("Report_ID","Diagnosis","Gender","DOB","Specimen_Site","Collection_Date",
                                  "Recieved_Date","Country_Of_Origin","Pathological_Purity",
                                  "Computational_Purity","TMB_Score","MS_Status")
    
    df_ShortVariants <- data.frame(Genes_CN$Genes,Protein_Change,Allele_Freqs,Genes_CN$Ploidy,Depths,Positions,Strand,
                                   Mutations,Functional_Effect)
    colnames(df_ShortVariants) <- c("Gene","Protein_Change","Allele_Freq","Ploidy","Depth","Position",
                                    "Strand","Mutation","Functional_Effect")
    
    df <- cbind(df_GeneralInfo,df_ShortVariants)
    df <- df[order(Allele_Freqs),] 
    
    # extract model predictions
    #pathModels <- read.delim("pathologicalPurity_Models_tabledata_out.txt", header = TRUE, sep = "\t", dec = ".")
    #compModdels <- read.delim("computationalPurity_Models_tabledata_out.txt", header = TRUE, sep = "\t", dec = ".")
    #allfitModels <- read.delim("AllFIT_Models_tabledata_out.txt", header = TRUE, sep = "\t", dec = ".")
    
    #df_models <- data.frame(pathModels$Model,compModdels$Model,allfitModels$Model)
    #colnames(df_models) <- c("Pathological Purity Model","Computational Purity Model","All-FIT Purity Model")
    
    #df <- cbind(df,df_models)
    
    contents <- df
    
    return(contents)
}
## -------------------------------------------------------------------------------------------------------
foundationXML_VUS <- function(inputFile){
    fileData <- xmlParse(inputFile)
    
    # Read the text from the file
    xmlText <- paste(readLines(inputFile), "\n", collapse="")
    
    # ID
    ReportID <- gsub('^.*<ReportId>\\s*|\\s*</ReportId>.*$', '', xmlText)
    # VUSs
    VariantProperties <- gsub('^.*<VariantProperties>\\s*|\\s*</VariantProperties>.*$', '', xmlText)
    VariantProperties <- str_split(VariantProperties, " \n ")
    VariantProperties <- VariantProperties[[1]][] # all VUSs with attributes
    VUSs <- matrix(data=NA, nrow=length(VariantProperties), ncol=2)
    for (i in 1 : length(VariantProperties)) {
        line <- VariantProperties[i]
        gName <- gsub('^.*geneName=\"\\s*|\\s*\".*$', '', line)
        variantName <- gsub('^.*variantName=\"\\s*|\\s*\".*$', '', line)
        
        if (str_length(gName)>12) {
            gName = ""
            variantName = ""
        }
        
        VUSs[i,1] <- gName
        VUSs[i,2] <- variantName
    }
    
    # Define Dataframe: 
    df_VUS <- data.frame(VUSs)
    df_IDs <- data.frame(rep(ReportID,length(VUSs)))
    df_VUS <- cbind(df_IDs,df_VUS)
    colnames(df_VUS) <- c("Report_ID","Gene","Protein_Change")
    
    VUS <- df_VUS
    
    return(VUS)
}
## -------------------------------------------------------------------------------------------------------
foundationXML_CNA <- function(inputFile){
    
    fileData <- xmlParse(inputFile)
    
    # Read the text from the file
    xmlText <- paste(readLines(inputFile), "\n", collapse="")
    
    # ID
    ReportID <- gsub('^.*<ReportId>\\s*|\\s*</ReportId>.*$', '', xmlText)
    
    # Copy Number Alterations
    CopyNumberAlterations <- gsub('^.*<copy-number-alterations>\\s*|\\s*</copy-number-alterations>.*$', '', xmlText)
    CopyNumberAlterations <- str_split(CopyNumberAlterations, " </copy-number-alteration> ")
    CopyNumberAlterations<- CopyNumberAlterations[[1]][] # all CNAs with attributes
    
    CNA_geneName <- c()
    CNA_copyNumber <- c()
    CNA_position <- c()
    CNA_type <- c()
    CNA_status <- c()
    CNA_numberOfExons <- c()
    for (i in 1 : length(CopyNumberAlterations)) {
        line <- CopyNumberAlterations[i]
        cna_geneName <- gsub('^.*gene=\"\\s*|\\s*\".*$', '', line)
        
        cna_copynumber <- gsub('^.*copy-number=\"\\s*|\\s*".*$', '', line)
        
        cna_pos <- gsub('^.*position=\"\\s*|\\s*\".*$', '', line)
        
        cna_tp <- gsub('^.*type=\"\\s*|\\s*\".*$', '', line)
        
        cna_stat <- gsub('^.*status=\"\\s*|\\s*\".*$', '', line)
        
        cna_exons <- gsub('^.*number-of-exons=\"\\s*|\\s*\".*$', '', line)
        
        if (str_length(cna_copynumber)>5) {
            cna_geneName = ""
            cna_copynumber = ""
            cna_pos = ""
            cna_tp = ""
            cna_stat = ""
            cna_exons = ""
        }
        CNA_geneName[i] <- cna_geneName
        CNA_copyNumber[i] <- cna_copynumber
        CNA_position[i] <- cna_pos
        CNA_type[i] <- cna_tp
        CNA_status[i] <- cna_stat
        CNA_numberOfExons[i] <- cna_exons
        
    }
    df_IDs <- data.frame(rep(ReportID,length(CNA_geneName)))
    df_CNA <- data.frame(CNA_geneName,CNA_copyNumber,CNA_position,CNA_type,CNA_status,CNA_numberOfExons)
    df_CNA <- cbind(df_IDs,df_CNA)
    colnames(df_CNA) <- c("Report_ID","Gene","Copy_Number","Position","Type","Status","number_Exons")
    
    # Define Dataframe: 
    
    CN_Alterations <- df_CNA
    
    return(CN_Alterations)
}

## -------------------------------------------------------------------------------------------------------
foundationXML_Rearrangements <- function(inputFile){
    
    fileData <- xmlParse(inputFile)
    
    # Read the text from the file
    xmlText <- paste(readLines(inputFile), "\n", collapse="")
    
    # ID
    ReportID <- gsub('^.*<ReportId>\\s*|\\s*</ReportId>.*$', '', xmlText)
    
    # Rearrangements
    Rearrangements <- gsub('^.*<rearrangements>\\s*|\\s*</rearrangements>.*$', '', xmlText)
    Rearrangements <- str_split(Rearrangements, " </rearrangement> ")
    Rearrangements<- Rearrangements[[1]][] # all CNAs with attributes
    
    Rearrangement_targetedGene <- c()
    Rearrangement_description <- c()
    Rearrangement_inFrame <- c()
    Rearrangement_otherGene <- c()
    Rearrangement_Position1 <- c()
    Rearrangement_Position2 <- c()
    Rearrangement_status <- c()
    Rearrangement_type <- c()
    for (i in 1 : length(Rearrangements)) {
        line <- Rearrangements[i]
        ra_targetedGene <- gsub('^.*targeted-gene=\"\\s*|\\s*\".*$', '', line)
        
        ra_description <- gsub('^.*<rearrangement description=\"\\s*|\\s*\".*$', '', line)
        
        ra_inFrame <- gsub('^.*in-frame=\"\\s*|\\s*\".*$', '', line)
        
        ra_otherGene <- gsub('^.*other-gene=\"\\s*|\\s*\".*$', '', line)
        
        ra_pos1 <- gsub('^.*pos1=\"\\s*|\\s*\".*$', '', line)
        
        ra_pos2 <- gsub('^.*pos2=\"\\s*|\\s*\".*$', '', line)
        
        ra_status <- gsub('^.*status=\"\\s*|\\s*\".*$', '', line)
        
        ra_type <- gsub('^.*type=\"\\s*|\\s*\".*$', '', line)
        
        if (str_length(ra_targetedGene)>20) {
            ra_targetedGene = ""
            ra_description = ""
            ra_inFrame = ""
            ra_otherGene = ""
            ra_pos1 = ""
            ra_pos2 = ""
            ra_status = ""
            ra_type = ""
        }
        Rearrangement_targetedGene[i] <- ra_targetedGene
        Rearrangement_description[i] <- ra_description
        Rearrangement_inFrame[i] <- ra_inFrame
        Rearrangement_otherGene[i] <- ra_otherGene
        Rearrangement_Position1[i] <- ra_pos1
        Rearrangement_Position2[i] <- ra_pos2
        Rearrangement_status[i] <- ra_status
        Rearrangement_type[i] <- ra_type
        
    }
    df_IDs <- data.frame(rep(ReportID,length(Rearrangement_targetedGene)))
    df_RA <- data.frame(Rearrangement_targetedGene,Rearrangement_description,
                        Rearrangement_otherGene,Rearrangement_Position1,Rearrangement_Position2,
                        Rearrangement_inFrame,Rearrangement_status,Rearrangement_type)
    df_RA <- cbind(df_IDs,df_RA)
    colnames(df_RA) <- c("Report_ID","Tageted_Gene","Description","Other_Gene","Position1",
                         "Position2","In_Frame","Status","Type")
    
    # Define Dataframe: 
    
    Rearrangements <- df_RA
    
    
    return(Rearrangements)
}

############################################################################################################
# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("cerulean"),

                
                navbarPage("FoundationMedicine Data Extraction",
                           tabPanel("Main",
                # topbar layout
                topbar( fileInput("file1", "Choose File",
                                  multiple = FALSE)),
                # Horizontal line ----
                tags$hr(),
                
                mainPanel(
                    fluidRow(
                        column(12, tabsetPanel(
                            id = 'dataset',
                            # Output: Data file ----
                            tabPanel("Variants", DT::dataTableOutput("contents", height = "500px")),
                            tabPanel("VUS", DT::dataTableOutput("VUS", height = "500px")),
                            tabPanel("CN_Alterations", DT::dataTableOutput("CN_Alterations", height = "500px")),
                            tabPanel("Rearrangements", DT::dataTableOutput("Rearrangements", height = "500px"))
                        ) # tabsetPanel end
                        ), # column end
                        #column(6, plotOutput("scatter1"))
                    ), # fluidRow end
                    
                    # horizontal line
                    tags$hr(),
                    
                    # plots for selcted row
                    fluidRow(
                        column(12,plotOutput("models_purities"))
                    )# fluid row end
                    
                    
                ) # mainPanel end
                           ), # main tab panel end
                tabPanel("Help & Documentaiton",
                         includeMarkdown("README.md")
                )
                ) # navbar end
                
) # fluidPage end



# Define server logic required to draw a histogram
server <- function(input, output) {

    inputData <- reactive({
        req(input$file1)
        inputFile <- input$file1$datapath
        contents <- foundationXML_general(inputFile)
    })
    
    output$contents <- DT::renderDataTable(
        datatable( data = inputData()
                   , extensions = 'Buttons',
                   selection = 'multiple'
                   , options = list( 
                       dom = "Blfrtip"
                       , buttons = 
                           list("copy", list(
                               extend = "collection"
                               , buttons = c("csv", "excel", "pdf")
                               , text = "Download"
                           ) ) # end of buttons customization
                   ) # end of options
        ) # end of datatables
    )
    
    output$VUS <- DT::renderDataTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        
        req(input$file1)
        inputFile <- input$file1$datapath
        VUS <- foundationXML_VUS(inputFile)
        
    },extensions = 'Buttons',
    options = list(dom = 'Blfrtip',
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                   lengthMenu = list(c(10,25,50,-1),
                                     c(10,25,50,"All"))))
    
    output$CN_Alterations <- DT::renderDataTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        req(input$file1)
        inputFile <- input$file1$datapath
        CN_Alterations <- foundationXML_CNA(inputFile)
        
        
    },extensions = 'Buttons',
    options = list(dom = 'Blfrtip',
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                   lengthMenu = list(c(10,25,50,-1),
                                     c(10,25,50,"All"))))
    
    output$Rearrangements <- DT::renderDataTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        req(input$file1)
        inputFile <- input$file1$datapath
        Rearrangements <- foundationXML_Rearrangements(inputFile)
        
        
    },extensions = 'Buttons',
    options = list(dom = 'Blfrtip',
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                   lengthMenu = list(c(10,25,50,-1),
                                     c(10,25,50,"All"))))
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
