#source("afd/R/ALF.R")
#num_trees <- 150
#num_genes <- 5
#seed = NULL
make_sims <- function(num_trees = 50, num_genes = 5, seed = NULL) {
  #run_ALF <- function(simname, nspec, ngenes, mingenelen, subrate, brate, drate, indelrate, dir, 
  #ALF_dir, seed = NULL)
  if(!is.null(seed)){
    set.seed(seed)
  }
  sim_params <- data_frame(nspec = sample(c(150:300), num_trees, replace = TRUE),
                           ngenes = num_genes,
                           mingenelen = 100,
                           subrate = runif(num_trees, 0.01, 2),
                           brate = 0.04,
                           drate = 0.025,
                           indelrate = runif(num_trees, 0.0001, 0.001),
                           rearrange = sample(c(TRUE, FALSE), num_trees, replace = TRUE, prob = c(0.25, 0.75)))
  do_ALF <- function(nspec, ngenes, mingenelen, subrate, brate, drate, indelrate, seed, rearrange) {
    #run_ALF <- function(simname, nspec, ngenes, mingenelen, subrate, brate, drate, indelrate, dir, 
    #ALF_dir, seed = NULL)
    dir <- "temp"
    test <- run_ALF("test", nspec, ngenes, mingenelen, subrate, brate, drate, indelrate, 
                    dir, "/setup_files/ALF_standalone/bin", seed, rearrange)
    
    testdna <- load_ALF(test)
    ALF_delete(test)
    Sys.sleep(3)
    testdna
  }
  sims <- sim_params %>%
    rowwise %>%
    do(nspec = .$nspec, ngenes = .$ngenes, mingenelen = .$mingenelen, subrate = .$subrate, 
       brate = .$brate, drate = .$drate, indelrate = .$indelrate, rearrange = .$rearrange,
       sims = do_ALF(.$nspec, .$ngenes, .$mingenelen, .$subrate, .$brate, .$drate, .$indelrate, seed,
                     .$rearrange))
  sims
}

fasta_it <- function(sim_data) {
  reshape_dna <- function(b, z) {
    labs <- names(z$dna)
    seqs <- sapply(z$dna, function(gg) as.character(gg[[b]]))
    names(seqs) <- labs
    DNAStringSet(seqs)
  }
  #y <- dna$dna_sets[[1]]
  concat_dna <- function(y) {
    labs <- names(y[[1]])
    reord <- lapply(y, function(z) z[labs])
    reord <- do.call(xscat, reord)
    names(reord) <- labs
    reord
  }
  dna <- sim_data %>%
    rowwise %>%
    do(dna_sets = lapply(seq_len(.$ngenes), reshape_dna, z = .$sims)) %>%
    do(dna_sets = concat_dna(.$dna_sets))
  dna
}

#x <- sim_data
align_sims <- function(dna) {
  #z <- x$sims[[1]]
  #b <- 5
  
  align <- function(y) {
    writeXStringSet(y, "temp/tempdna.fasta")
    start <- Sys.time()
    system("/muscle/muscle -in /home/rstudio/afd/temp/tempdna.fasta -out /home/rstudio/afd/temp/alignment.fa",
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    end <- Sys.time()
    aligned <- readDNAMultipleAlignment("/home/rstudio/afd/temp/alignment.fa")
    time_taken <- end - start
    data_frame(alignment = list(aligned), time_taken = list(time_taken))
  }
  alignments <- dna %>%
    rowwise %>%
    do(align(.$dna_sets))
  #alignments <- dna %>%
  #  rowwise %>%
  #  do(alignment = lapply(.$dna_sets, align))
  #alignments <- alignments %>%
  #  rowwise %>%
  #  do(alignment = do.call(xscat, lapply(.$alignment, DNAStringSet)))
  alignments
}

#aligns <- align_data
raxML_it <- function(aligns) {
  #x <- aligns$alignment[[1]]
  get_tree <- function(x) {
    writeXStringSet(DNAStringSet(x), "temp/alignment.fa")
    start <- Sys.time()
    system("/raxml/raxmlHPC-SSE3 -s /home/rstudio/afd/temp/alignment.fa -n tree.phy -m GTRGAMMA -p 1149135 -w /home/rstudio/afd/temp",
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    end <- Sys.time()
    time_taken <- end - start
    tree <- read.tree("/home/rstudio/afd/temp/RAxML_bestTree.tree.phy")
    file.remove(list.files("/home/rstudio/afd/temp", full.names = TRUE))
    data_frame(align_tree = list(tree), time_taken = list(time_taken))
  }
  align_trees <- aligns %>%
    rowwise %>%
    do(get_tree(.$alignment))
  align_trees
}

do_FFP <- function(temppath, k=5, ffppath, outpath)
align_free_it <- function(fasta_dna) {
  dna_sets <- fasta_dna$dna_sets[[1]]
  do_align_free <- function(dna_sets) {
    writeXStringSet(dna_sets, "temp/tempdna.fasta")
    FFP_5 <- do_FFP("temp/tempdna.fasta", ffppath = "/afd-python/ffp-3.19/src", temppath = "temp")
  }
  writeXStringSet(y, "temp/tempdna.fasta")
  align_free <- fasta_dna %>%
    rowwise %>%
    do(do_align_free(dna_sets))
}

