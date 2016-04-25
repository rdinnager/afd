#source("afd/R/ALF.R")
#num_trees <- 150
#num_genes <- 5
#seed = NULL
make_sims <- function(num_trees = 150, num_genes = 5, seed = NULL) {
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
                           rearrange = sample(c(TRUE, FALSE), num_trees, replace = TRUE))
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