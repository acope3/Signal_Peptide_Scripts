
library(AnaCoDa)
parameter.file.to.fix <- "../Results/Full_genome_main_ht/chain/run_2/R_objects/parameter.Rda"
load(parameter.file.to.fix)

aminoAcids()
paramBase$grouplist <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y", "Z")

save(list = c("paramBase", "currentMutation", "currentSelection",
              "proposedMutation", "proposedSelection", "model",  
              "mutationPrior", "mutationTrace", "selectionTrace", 
              "synthesisOffsetAcceptRatTrace", "synthesisOffsetTrace", 
              "observedSynthesisNoiseTrace", "withPhi"),
     file=parameter.file.to.fix)