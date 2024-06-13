# Loading required libraries
library(shiny) # For creating the Shiny web application
library(Biostrings) # For handling biological strings and sequences
library(ggtree) # For visualizing phylogenetic trees
library(ape) # For phylogenetic analysis and tree manipulation
library(msa) # For multiple sequence alignment
library(ggplot2) # For creating plots and visualizations
library(DT) # For interactive data tables

server <- function(input, output) {
  
  # Create a reactive value to store the current output
  current_output <- reactiveVal(NULL) #Because I wanted to show only one output at the time, not all of them at once or  below each other
  
  # Observe event for displaying sequence table when 'table' ("What's inside?") button is clicked
  observeEvent(input$table, {
    current_output(renderDataTable({
      # Make sure a file has been uploaded
      req(input$file)
      tryCatch(
        {
          # Load the sequence data from the uploaded file
          fasta <- readDNAStringSet(input$file$datapath)
          
          # Create a data frame with sequence names and sequences
          df <- data.frame(
            Name = names(fasta),
            Sequence = as.character(fasta) # Convert sequences to character format
          )
          datatable(df, 
                    rownames = FALSE, # Disable rownames
                    # Because I want to have my own row with names so the user can sort with this row!
                    # It surprised me that you cannot sort with column "rownames" :( so I needed to add this line
                    options = list(
                      paging = TRUE, # Enable pagination
                      lengthMenu = list(c(10, 25, 50, 100, -1), c("10","25", "50", "100", "All")), # Number of rows per page, with option for "All"
                      # In case you want to find the one with the lowest or the highest C nucleotide percentage! 
                      pageLength = 10 # Initial number of rows per page
                    ) # Options only accepts values as lists, even if you want to change one option, it's weird and caused some troubles! :D
          )
        },
        error = function(e) {
          # Handle any errors that occur during file reading or processing
          stop(safeError(e))
        }
      )
    }))
  })
  
  # Observe event for performing multiple sequence alignment when 'clustal'( "Clustal Omega dendrogram!")  button is clicked
  observeEvent(input$clustal, {
    current_output(renderPlot({
      # Make sure a file has been uploaded
      req(input$file)
      withProgress(message = 'Generating dendrogram...', value = 0, { # Use withProgress to show a progress indicator
        # I tried with library(shinycssloaders) and withSpinner() but it doesn't work as good because of my approach to outputs as dynamic
        # I only did it here because creating the phylogenetic tree from msa takes the longest! 
        # So you can know it actually work and it didn't crash
        tryCatch(
          {
            # Load the sequence data from the uploaded file
            sequences <- readDNAStringSet(input$file$datapath)
            
            # Increment progress and provide detail message
            incProgress(0.3, detail = "Aligning sequences...")  # Increment progress and provide detail message
            
            # Perform multiple sequence alignment using Clustal Omega
            alignment <- msa(sequences, method = "ClustalOmega")
            
            # Convert alignment to a matrix
            alignment_matrix <- as.matrix(alignment)
            
            # Convert the alignment matrix to a DNAbin object required by ape
            alignment_dnabin <- as.DNAbin(alignment_matrix)
            
            # Compute a distance matrix using raw distance model
            dist_matrix <- dist.dna(alignment_dnabin, model = "raw")
            
            # Increment progress and provide detail message
            incProgress(0.7, detail = "Creating UPGMA tree...")
            
            # Create a UPGMA tree using hclust function 
            hc <- hclust(dist_matrix, method = "average")
            tree <- as.phylo(hc)
            
            # Visualize the tree using ggtree
            ggtree(tree, branch.length = "none") + geom_tiplab() + xlim_tree(14) 
            # the xlim was necessary because the names of, for example, influenza viruses was being cut short!
          },
          error = function(e) {
            # Handle any errors that occur during alignment or tree creation
            stop(safeError(e))
          }
        )
      })
    }))
  })
  
  # Observe event for displaying sequence information when 'info' ("Inform me in general!") button is clicked
  observeEvent(input$info, {
    current_output(renderUI({ # I started with renderPrint but the font was too small so I changed it to renderUI to get bigger text!
      # Ensure a file has been uploaded
      req(input$file)
      tryCatch(
        {
          # Load the sequence data from the uploaded file
          fasta <- readDNAStringSet(input$file$datapath)
          
          # Extract sequence information
          Name <- names(fasta)
          Sequence <- as.character(fasta)
          Length <- nchar(Sequence)
          
          # Calculate number of sequences, the shortest sequence, the longest, and mean sequence lengths
          num_sequences <- length(Name)
          min_length <- min(Length)
          max_length <- max(Length)
          mean_length <- mean(Length)
          
          # Create HTML content with larger text size
          tags$div(
            style = "font-size: 18px;",  # Adjust font size 
            HTML(
              paste(
                "<strong>General sequence information:</strong><br><br>",
                "Number of sequences: ", num_sequences, "<br>",
                "The shortest sequence length: ", min_length, "<br>",
                "The longest sequence length: ", max_length, "<br>",
                "Mean sequence length: ", round(mean_length), "<br>"
              )
            )
          )
        },
        error = function(e) {
          # Handle any errors that occur during information extraction
          stop(safeError(e))
        }
      )
    }))
  })
  
  # Observe event for displaying detailed sequence information s when 'info2' ("More info about the details!") button is clicked
  observeEvent(input$info2, {
    current_output(renderDataTable({
      req(input$file)
      tryCatch(
        {
          # Load the sequence data from the uploaded file
          fasta <- readDNAStringSet(input$file$datapath)
          
          # Extract sequence information
          Name <- names(fasta)
          Sequence <- as.character(fasta)
          
          # Determine if sequences are nucleotide or amino acid based on content
          Type <- ifelse(grepl("[^ACGTacgt]", Sequence), "Amino acid", "Nucleotide")
          
          # Calculate sequence lengths
          Length <- nchar(Sequence)
          
          # Function to calculate nucleotide quantity for a given sequence
          calc_nucleotide_percentages <- function(seq) {
            total <- nchar(seq) # Calculate the total length of the sequence
            A_count <- sum(unlist(strsplit(seq, "")) == "A") # Count occurrences of 'A'
            C_count <- sum(unlist(strsplit(seq, "")) == "C") # Count occurrences of 'C'
            G_count <- sum(unlist(strsplit(seq, "")) == "G") # Count occurrences of 'G'
            T_count <- sum(unlist(strsplit(seq, "")) == "T") # Count occurrences of 'T'
            
            # I guess it is a better way of doing that by using bioconductor and alphabetFrequency method 
            # but I really couldn't figure it out to work with my code :( 
            # And after many tries I was so annoyed I ditched the idea.
            # But it works without bioconductor anyway, so it's still win situation!
            
            # Calculate percentages and round them to two decimal places
            list(
              A_Percentage = round((A_count / total) * 100, digits = 2),
              C_Percentage = round((C_count / total) * 100, digits = 2),
              G_Percentage = round((G_count / total) * 100, digits = 2),
              T_Percentage = round((T_count / total) * 100, digits = 2)
            )
          }
          
          # Initialize vectors to store nucleotide percentages
          A_Percentage <- numeric(length(Sequence))
          C_Percentage <- numeric(length(Sequence))
          G_Percentage <- numeric(length(Sequence))
          T_Percentage <- numeric(length(Sequence))
          
          # Calculate percentages for nucleotide sequences
          for (i in seq_along(Sequence)) {
            if (Type[i] == "Nucleotide") {
              percentages <- calc_nucleotide_percentages(Sequence[i])
              A_Percentage[i] <- percentages$A_Percentage
              C_Percentage[i] <- percentages$C_Percentage
              G_Percentage[i] <- percentages$G_Percentage
              T_Percentage[i] <- percentages$T_Percentage
            } else {
              A_Percentage[i] <- NA
              C_Percentage[i] <- NA
              G_Percentage[i] <- NA
              T_Percentage[i] <- NA
            } # First, I wanted to not show them but if someone has a mixed file (amino-acid and nucleotide) it would be a problem 
          }
          
          # Create the final data frame with sequence information and nucleotide percentages
          df <- data.frame(
            Name = Name, # Sequence names
            Length = Length, # Sequence lengths
            Type = Type, # Sequence types (Amino acid or Nucleotide)
            A_Percentage = A_Percentage, # Percentage of A nucleotides
            C_Percentage = C_Percentage, # Percentage of C nucleotides
            G_Percentage = G_Percentage, # Percentage of G nucleotides
            T_Percentage = T_Percentage  # Percentage of T nucleotides
          )
          
          # Render the table with DT::datatable and customize options
          datatable(df, 
                    rownames = FALSE,  # Disable rownames
                    options = list(
                      paging = TRUE, # Enable pagination
                      lengthMenu = list(c(10, 25, 50, 100, -1), c("10","25", "50", "100", "All")), # Number of rows per page, with option for "All"
                      pageLength = 10 # Initial number of rows per page
                    )
          )
        },
        error = function(e) {
          # Handle any errors that occur during processing
          stop(safeError(e))
        }
      )
    }))
  })
  
  
  # Render the current output dynamically in the UI
  output$dynamic_output <- renderUI({
    current_output()
  })
  
}