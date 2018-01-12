require(grDevices)


plot_blocks <- function(block_filename, min_contig_size = 10000) {
  species_colours <- list("no call" = "gray90",
                          "too close" = "gray70",
                          "pneumoniae" = "red",
                          "quasipneumoniae" = "blue",
                          "variicola" = "green",
                          "oxytoca" = "magenta",
                          "michiganensis" = "cyan")
  
  blocks <- read.delim(block_filename,
                       colClasses=c("character", "numeric", "numeric", "character"))
  contigs <- unique(blocks$contig)
  max_size <- max(blocks$end)

  last_contig_index <- 0
  for (i in seq_along(contigs)){
    contig = contigs[i]
    contig_blocks <- blocks[blocks$contig == contig,]
    contig_size <- max(contig_blocks$end) - min(contig_blocks$start)
    if (contig_size >= min_contig_size) {
      last_contig_index <- i
    }
  }
  contigs <- contigs[1:last_contig_index]
  contig_count <- length(contigs)
  
  par(mar=c(5.1,2.1,4.1,2.1))
  plot(c(0, max_size),
       c(1, -contig_count - 1),
       type = "n", xlab = "", ylab = "", main = block_filename,
       yaxt = "n")
  
  for (i in seq_along(contigs)){
    contig = contigs[i]
    contig_blocks <- blocks[blocks$contig == contig,]
    
    for(j in 1:nrow(contig_blocks)) {
      block <- contig_blocks[j,]
      bottom <- -i + 0.1
      top <- bottom + 0.8
      left <- block$start
      right <- block$end
      
      rect(left, bottom, right, top, col = species_colours[[block$classification]], border = NA)
    }
  }
  
  legend("bottom",
         legend = c("pneumoniae", "quasipneumoniae", "variicola", "oxytoca", "michiganensis", "no call", "too close"),
         fill = c("red", "blue", "green", "magenta", "cyan", "gray90", "gray70"),
         ncol = 4)
}
