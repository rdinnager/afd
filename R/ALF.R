## Generate ALF simulations
library(whisker)
library(Biostrings)
library(ape)
ALF_template <- function() {
  "  webRequest := false;
  uuid := '4e4937bd-70e5-4caf-8521-f4340a4b7e09';
  
  # name of simulation - you may want to change this
  mname := {{{simname}}};
  
  # directories for file storage - you may want to change these
  wdir := '{{{dir}}}';
  dbdir := 'DB/';
  dbAncdir := 'DBancestral/';
  
  # time scale for simulation (PAM is default)
  unitIsPam := true:
    
  # parameters concerning the root genome
  realseed := false;
  protStart := {{ngenes}};
  minGeneLength := {{mingenelen}};
  gammaLengthDist := [2.4019, 133.8063];
  blocksize := 3:
    
  # parameters concerning the species tree
  treeType := 'BDTree';
  birthRate := {{brate}};
  deathRate := {{drate}};
  NSpecies := {{nspec}};
  ultrametric := true;
  mutRate := {{PAMunits}};
  scaleTree := false;
  
  # parameters concerning the substitution models
  substModels := [SubstitutionModel('CPAM')];
  indelModels := [IndelModel({{indelrate}}, ZIPF, [1.821], 50)];
  rateVarModels := [RateVarModel(Gamma, 5, 0.01, 1)];
  modelAssignments := [1]:
  modelSwitchS := [[1]]:
  modelSwitchD := [[1]]:
    
  # parameters concerning gene duplication
  geneDuplRate := 0;
  numberDupl := 0;
  transDupl := 0;
  fissionDupl := 0;
  fusionDupl := 0;
  P_pseudogene := 0;
  P_neofunc := 0;
  P_subfunc := 0;
  
  # parameters concerning gene loss
  geneLossRate := 0;
  
  # parameters concerning LGT
  lgtRate := 0;
  lgtGRate := 0;
  lgtGSize := 0;
  
  # parameters concerning rate heterogeneity among genes
  amongGeneDistr := 'Gamma';
  aGAlpha := 1;"
  
}

#' a function to create an ALF configuration file (*.drw)
#' @param simname Name for the simulation to be simulated
#' @param nspec Number of species to simulate
#' @param ngenes Number of genes to simulate in the genome
#' @param mingenlen Minimum length of a gene under the simulation
#' @param PAMunits Distance from root to tips on simulated phylogeny, in PAM units: time it takes
#' to accumulate 1 changed amino acid out of every 100 amino acids (codon based evolutionary unit)
#' @param brate Birth rate for Birth-Death phylogeny simulation 
#' @param drate Death rate for Birth-Death phylogeny simulation
#' @param indelrate Rate at which insertions and deletion are incorporated into genomes
#' @param dir Directory where the simulated genome should be stored
#' @export
gen_ALF_drw <- function(simname, nspec, ngenes, mingenelen, PAMunits, brate, drate, indelrate, dir) {
  parms <- list(simname = simname, nspec = nspec, ngenes = ngenes, mingenelen = mingenelen,
                PAMunits = PAMunits, brate = brate, drate = drate, indelrate = indelrate,
                dir = dir)
  whisker.render(ALF_template(), parms) 
}

#' a function to create an ALF configuration file (*.drw)
#' @param simname Name for the simulation to be simulated
#' @param nspec Number of species to simulate
#' @param ngenes Number of genes to simulate in the genome
#' @param mingenlen Minimum length of a gene under the simulation
#' @param PAMunits Distance from root to tips on simulated phylogeny, in PAM units: time it takes
#' to accumulate 1 changed amino acid out of every 100 amino acids (codon based evolutionary unit)
#' @param brate Birth rate for Birth-Death phylogeny simulation 
#' @param drate Death rate for Birth-Death phylogeny simulation
#' @param indelrate Rate at which insertions and deletion are incorporated into genomes
#' @param dir Directory where the simulated genome should be stored
#' @param ALF_dir Directory containing the \code{alfsim} executable file
#' @return A character vector containing the path to the directory in which the ALF simulation saves its output.
#' @export
run_ALF <- function(simname, nspec, ngenes, mingenelen, PAMunits, brate, drate, indelrate, dir, 
                    ALF_dir) {
  sdir <- paste0(dir, "/", simname)
  filename <- paste0(dir, "/", simname, ".drw")
  cat(gen_ALF_drw(simname, nspec, ngenes, mingenelen, PAMunits, brate, drate, indelrate, dir), file = filename)
  system(paste0(ALF_dir, "/alfsim ", filename))
  return(sdir)
}

#' Load sequences and trees from an ALF simulation into R
#' @param simdir The directory where the simulation output you wish to load is located. This is returned by
#' \code{\link{run_ALF}}.
#' @param which_specs A vector of integers expressing which species sequences you wish to load
#' @param return_tree Boolean. Returns the underlying tree if TRUE, only the sequences if FALSE
#' @param true_align Boolean. Should the function load the true alignement from the simulation. Currently not used.
#' @return An "ALFsim" object containing the requested simulated sequences and the simulated tree 
#' (if return_tree = TRUE).
#' @import Biostrings ape
#' @export
load_ALF <- function(simdir, which_specs = NULL, return_tree = TRUE, true_align = TRUE) {
  results <- list()
  tl <- list.files(paste0(simdir,"/DB"), pattern = "dna.fa")
  if (is.null(which_specs)) which_specs <- seq_along(tl) 
  labels <- sapply(strsplit(tl, "_", fixed = TRUE), function(x) x[1])[which_specs]
  dna <- lapply(paste0(simdir, "/DB/", tl[which_specs]), readDNAStringSet, format="fasta", use.names=TRUE)
  results[["dna"]] <- dna
  names(results$dna) <- labels
  if (return_tree) {
    tree <- read.tree(paste0(simdir, "/RealTree.nwk"))
    results[["tree"]] <- drop.tip(tree, which(!tree$tip.label %in% labels))
  }
  class(results) <- "ALFsim"
  return(results)
}

#' Function to concatenate all gene sequences into one sequence for each species
#' @param ALF_sim An "ALFsim" object containing the sequence to be concatenated
#' @return An "ALFsim" object with concatenated sequences
#' @export
ALF_cat <- function(ALF_sim) {
  dna <- lapply(ALF_sim$dna, function(x) DNAStringSet(unlist(x)))
  ALF_sim$dna <- dna
  return(ALF_sim)
}

## example
#dir <- "afd/temp"
#test <- run_ALF("test", 50, 5, 400, 100, 0.04, 0.025, 0.0001, dir, "/setup_files/ALF_standalone/bin")

#testdna <- load_ALF(test)
#testcat <- ALF_cat(testdna)
