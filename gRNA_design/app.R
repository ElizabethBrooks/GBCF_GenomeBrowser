# created: Sheri Sanders
# modified: Elizabeth Brooks
# updated: 19 January 2026

# gRNA design Shiny application using crispRdesignR

## Data Setup ##
## forge and load genome data
## https://www.bioconductor.org/packages//2.7/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
##library(BSgenomeForge)
##BSgenomeForge::forgeBSgenomeDataPkg("./data/BSgenome.Dmagna.LRV0", replace=TRUE)
## once forgeBSgenomeDataPkg is done build the source package in the terminal
##R CMD build <pkgdir>
##R CMD check <tarball>
##R CMD INSTALL <tarball>

## Posit Connect Cloud Setup ##
## generate a manifest.json for posit connect cloud
## run in the console
##setwd("/Users/bamflappy/Repos/GBCF_GenomeBrowser/gRNA_design")
##library(rsconnect)
options(rsconnect.max.bundle.files = 60000)
##writeManifest()

# lists of dependent packages
packageList <- c("BiocManager", "shiny", "seqinr", "kableExtra")
biocList <- c("rtracklayer", "GenomicAlignments", "BSgenome", "crispRdesignR")
# check for any missing packages
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
newBioc <- biocList[!(biocList %in% installed.packages()[,"Package"])]
# install any missing packages
if(length(newPackages)){
  install.packages(newPackages)
}
if(length(newBioc)){
  BiocManager::install(newBioc)
}

#library(shiny)
library(crispRdesignR)
library(seqinr)

installed_genomes <- BSgenome::installed.genomes()
installed_genomes_names <- c()
if (length(installed_genomes) == 0) {
  installed_genomes[length(installed_genomes)+1] <- "no_genomes_installed"
  installed_genomes_names[length(installed_genomes_names)+1] <- "No genomes installed"
  names(installed_genomes) <- installed_genomes_names
} else {
  for (i in 1:length(installed_genomes)) {
    genome_name <- paste(BSgenome::organism(BSgenome::getBSgenome(installed_genomes[i])), " (", metadata(get(installed_genomes[i]))$genome, ")", sep="")
    installed_genomes_names[length(installed_genomes_names)+1] <- genome_name
  }
  names(installed_genomes) <- installed_genomes_names
}

gene_list=list("Dmagna031332-T1","Dmagna034057-T2")
gene_list=list.files("./pre_run/", pattern="*hits.RDS", all.files=TRUE, full.names=FALSE)
gene_list=lapply(X=gene_list, FUN = function(t) gsub(pattern=".hits.RDS", replacement="", x=t, fixed=TRUE))

ui <- fluidPage(
  navbarPage("crispRdesignR",
             tabPanel("sgRNA Finder",
                      titlePanel("sgRNA Finder"),
                      sidebarLayout(
                        sidebarPanel(
                          tags$div(id = "placeholder1"),
                          selectizeInput("genome_select", "Select Genome",
                                      installed_genomes),
                          textInput("gene_select","Enter Gene Name","Dmagna034057-T1"),
                          tags$div(id = "placeholder5"),
                          actionButton("run", "Find sgRNA", icon("paper-plane"))
                        ),
                        mainPanel(
                          tags$div(id = "placeholder3"),
                          DT::dataTableOutput("sgRNA_data"),
                          tags$div(id = "placeholder4"),
                          DT::dataTableOutput("offtarget_data"),
                          titlePanel("About"),
                          column(12, HTML("crispRdesignR designs guide RNA sequences (sgRNA) for Cas9 DNA editing.
                                          To begin, enter a sequence into the sequence box, select a genome to search for
                                          Off-Targets, provide a genome annotation file (.gtf) specific to your genome, and click find sgRNA. <br/><br/> Note about Off-target calling in large genomes: When using a large genome like
                                          Homo sapiens, we reccomend using sequences under 250 base pairs. The time it can take
                                          to search these genomes can be multiple hours if too many sgRNA are generated."))
                          )
                          )
                          )
                          )
                        )

server <- function(input, output) {

  ## Creates a list of reactive values that allows the program to
  ## update only when the action button is pressed
  maindf <- reactiveValues(data = NULL)
  downloadmaindf <- reactiveValues(data = NULL)
  offtargetdf <- reactiveValues(data = NULL)
  downloadofftargetdf <- reactiveValues(data = NULL)

  ## Creates default values for the arguments in the find sgRNA function
  callofftargets <- "yes_off"
  annotateofftargets <- "yes_annotate"
  givenPAM <- "NGG"

  ## Creates a variable that assists with adding UI elements
  n <- reactiveVal(0)

  ## Creates a variable for the gene annotation file
  gtf_datapath <- reactiveVal(0)
  gene_annotation_file <- reactiveVal(0)
  ## Runs the sgRNA_design function when the action button is pressed
  observeEvent(input$run, {
     name=input$'gene_select'

     callofftargets <- "yes_off"
     annotateofftargets <- "yes_annotate"
      # Create a Progress object
          designprogress <- shiny::Progress$new()
          designprogress$set(message = "Preparing gene annotation file", value = 0, detail = "This may take a while")
      # Close the progress when this reactive exits (even if there's an error)
          on.exit(designprogress$close())
          if (callofftargets == "no_off" | annotateofftargets == "no_annotate") {
             annotating <- FALSE
          } else {
             annotating <- TRUE
          }
          if (annotating == FALSE) {
             gene_annotation_file("placeholder")
          }

           filename = paste("./pre_run/",input$'gene_select',".hits.RDS",sep="")
           hits = readRDS(file=filename)

           if ((length(hits) == 0) == FALSE) {
          ## Starts creating the sgRNA table
              filename2 = paste("./pre_run/",input$'gene_select',".int_sgRNA.RDS",sep="")
              int_sgRNA_data = readRDS(filename2)

          ## Adds color to indicate unfavorable GC content
              GCinstance <- unlist(int_sgRNA_data[6])*100
              GCindex <- which(GCinstance >=80 | GCinstance <=30)
              GCinstance_color <- as.character(GCinstance)
              for (G in 1:length(GCinstance)){
                 if (G %in% GCindex){
                    GCinstance_color[G] <- paste('<span style="color:red">', GCinstance_color[G], '<span>', sep = "")
                 } 
              }
          ## Adds color to indicate unfavorable homopolymers
              Homopolymerdetect <- unlist(int_sgRNA_data[7])
              Homopolymerdetect_color <- as.character(Homopolymerdetect)
              for (H in 1:length(Homopolymerdetect)){
                  if (Homopolymerdetect[H] == "TRUE"){
                      Homopolymerdetect_color[H] <- paste('<span style="color:red">', Homopolymerdetect[H], '<span>', sep = "")
                  }
              }
          ## Adds color to indicate Self-Complementarity
              Self_comp_list <- unlist(int_sgRNA_data[8])
              Self_comp_index <- which(Self_comp_list > 0)
              Self_comp_list_color <- as.character(Self_comp_list)
              for (C in 1:length(Self_comp_list)){
                   if (C %in% Self_comp_index){
                      Self_comp_list_color[C] <- paste('<span style="color:red">', Self_comp_list_color[C], '<span>', sep = "")
                   }
              }
              proc_sgRNA_data <- data.frame(int_sgRNA_data[1:5], GCinstance_color, Homopolymerdetect_color, Self_comp_list_color, int_sgRNA_data[9:15])
              colnames(int_sgRNA_data) <- c("sgRNA sequence", "PAM", "Strand", "Start", "End", "GC content",
                                        "Homopolymer", "Self Complementary", "Efficiency Score", "MM0", "MM1", "MM2", "MM3", "MM4", "Note Codes")
              colnames(proc_sgRNA_data) <- c("sgRNA sequence", "PAM", "Strand", "Start", "End", "GC %",
                                         "Homopolymer", "Self Complementary", "Efficiency Score", "MM0", "MM1", "MM2", "MM3", "MM4", "Note Codes")


          ## Adds a title and download button for the sgRNA table to the UI
              n(n()+1)
              if (n() == 1) {
              title = paste("sgRNA Table", sep="")
              insertUI(
                 selector = "#placeholder3",
                 where = "afterEnd",
                 ui = tags$div(id = 'sgRNAdftext',
                            titlePanel(title),
                            column(12, "Note Codes: GC - Unfavorable GC content (=< 80% or => 30%), HP - Homopolymer detected (4 or more consectutive base pairs),
                                   SC - Region of self complementarity detected, LE - Low efficiency score (< 0.5)"),
                            downloadButton("Download_sgRNA", "Download sgRNA")
                 )
              )
              }

          ## Outputs the Table
            maindf$sgRNA_data <- proc_sgRNA_data
            downloadmaindf$sgRNA_data <- int_sgRNA_data
            filename3 = paste("./pre_run/",input$'gene_select',".int_offtarget.RDS",sep="")
            int_offtarget_data = readRDS(filename3)

            ## Adds code to color mismatches red within the off target sequences
            off_offseq <- as.character(unlist(int_offtarget_data[8]))
            off_sgRNAseq <- as.character(unlist(int_offtarget_data[1]))
            for (x in 1:length(off_offseq)) {
              justsgRNA <- off_sgRNAseq[x]
              justoff <- off_offseq[x]
              splitjustsgRNA <- stringr::str_split(justsgRNA, "", simplify = TRUE)
              splitoffsgRNA <- stringr::str_split(justoff, "", simplify = TRUE)
              #mismatches <- which(splitjustsgRNA != splitoffsgRNA)
              splitlistoffsgRNA <- as.list(splitoffsgRNA)
              #if (length(mismatches) != 0){
              #  for (g in length(mismatches):1) {
              #     splitlistoffsgRNA <- append(splitlistoffsgRNA, '</span>', after = mismatches[g])
              #     splitlistoffsgRNA <- append(splitlistoffsgRNA, '<span style="color:red">', after = mismatches[g]-1)
              #  }
              #  off_offseq[x] <- paste(splitlistoffsgRNA, sep="", collapse = "")
              #}
            }

            proc_offtarget_data <- data.frame(int_offtarget_data[1:7], off_offseq, int_offtarget_data[9:12])

          ####added 5/17/23 SS - Creates link to JBrowse (currently hardcoded) for the Off-targets
            library(stringr)
            library(kableExtra)
            url=paste("http://daphnia-db.crc.nd.edu/jbrowse2/?loc=",hits$all_offtarget_info.Chromosome,'&assembly=Dmagna_LRV0_1_genome&sessionTracks=[{%22type%22:%22FeatureTrack%22,%22trackId%22:%22CrisprTargets%22,%22name%22:%22CrisprTargets%22,%22assemblyNames%22:[%22Dmagna_LRV0_1_genome%22],%22adapter%22:{%22type%22:%22FromConfigAdapter%22,%22features%22:[{%22uniqueId%22:%22primer%22,%22refName%22:%22', hits$all_offtarget_info.Chromosome,"%22,%22start%22:", hits$all_offtarget_info.Start,",%22end%22:", hits$all_offtarget_info.End, "}]}}]&tracks=CrisprTargets,Dmagna_LRV0_1_genome.gff", sep="")
            url = str_replace_all(url, "chr","scaffold_")
            url = text_spec("Jbrowse", link = url)

            proc_offtarget_data = cbind(proc_offtarget_data,url)
            colnames(int_offtarget_data) <- c("sgRNA sequence", "Chromosome", "Start", "End", "Mismatches", "Strand", "CFD Scores",
                                            "Off-target sequence", "Gene ID", "Gene Name", "Sequence Type", "Exon Number")

          #####added link to below
            colnames(proc_offtarget_data) <- c("sgRNA sequence", "Chr", "Start", "End", "Mismatches", "Strand", "CFD Scores",
                                            "Off-target sequence", "Gene ID", "Gene Name", "Sequence Type", "Exon Number", "JBlink")

          ## Adds a title and download button for the off-target table to the UI
            if (n() == 1) {
               insertUI(
                selector = "#placeholder4",
                where = "afterEnd",
                ui = tags$div(id = 'sgRNAofftext',
                            titlePanel("Off-target Information"),
                            column(12, "Note: this program may report sequences in the target region as potential off-target sequences"),
                            downloadButton("Download_off", "Download Off-Targets")
                )
               )
            }

            offtargetdf$data <- proc_offtarget_data
            downloadofftargetdf$data <- int_offtarget_data
          } else {
              showModal(modalDialog(
            title = "Error",
            "No sgRNA were generated from sequence"
            ))
          } 
  }
)

  output$Download_sgRNA <- downloadHandler(
    filename = function(){"sgRNA.csv"},
    content = function(file) {
      write.csv(downloadmaindf$sgRNA_data, file, row.names = TRUE)
    }
  )

  output$Download_off <- downloadHandler(
    filename = function(){"Offtarget.csv"},
    content = function(file) {
      write.csv(downloadofftargetdf$data, file, row.names = TRUE)
    }
  )

  ## Reactively outputs an sgRNA table when the function is complete
  output$sgRNA_data <- DT::renderDataTable({maindf$sgRNA_data}, escape = FALSE)
  output$offtarget_data <- DT::renderDataTable({offtargetdf$data}, escape = FALSE)

  ## Add fasta file input to the UI
  observeEvent(input$fasta, {
    if (input$fasta == TRUE) {
      insertUI(
        selector = "#placeholder1",
        where = "afterEnd",
        ui = tags$div(id = 'fastainput',
                      fileInput("fastafile", "Choose fasta file",
                                multiple = FALSE)
        )
      )
    } else {
      removeUI(
        selector = 'div#fastainput',
        multiple = FALSE
      )
    }
  })

  ## Add Additional Options input to the UI

  observeEvent(input$options_toggle, {
    if (input$options_toggle == TRUE) {
      insertUI(
        selector = "#placeholder5",
        where = "afterEnd",
        ui = tags$div(id = 'optionsmenu',
                      column(12, HTML("Warning: Doench score not accurate for custom PAMS")),
                      textInput("customPAM", "Custom PAM (Max 6bp)", value = "NGG"),
                      selectInput("toggle_off_targets", "Call Off-Targets?",
                                  c("Yes" = "yes_off",
                                    "No" = "no_off"),
                                  selected = "yes_off"),
                      selectInput("toggle_off_annotation", "Annotate Off-Targets?",
                                  c("No" = "no_annotate",
                                    "Yes" = "yes_annotate"),
                                  selected = "yes_annotate")
        )
      )
    } else {
      removeUI(
        selector = 'div#optionsmenu',
        multiple = TRUE
      )
    }
  })

}

shinyApp(ui=ui, server=server)
