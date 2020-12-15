#' @title peptide property calculation
#' @description calculate various different peptide properties according to the peptide sequence.
#' @return The output from \code{\link{pep_property}}
#' @export
#' @param pep The input peptide sequence
#' @examples
#' pep_property("AAAAAAAAAA")
#' @import Peptides

pep_property<-function(pep){

  Pep_aIndex<-Peptides::aIndex(pep)
  Pep_charge<-Peptides::charge(pep, pH = 13, pKscale = "EMBOSS")
  PeP_pI<-Peptides::pI(pep,pKscale = "EMBOSS")
  Pep_InstabilityIndex<-Peptides::instaIndex(pep)
  Pep_hydrophobicity<-Peptides::hydrophobicity(pep,scale = "Eisenberg")
  Pep_comp<-Peptides::aaComp(pep)
  Pep_comp1 <- data.frame(matrix(nrow=1, ncol=18))
  Pep_property_table<- cbind(pep,
                            Pep_aIndex,
                                          Pep_charge,
                                          PeP_pI,
                                          Pep_InstabilityIndex,
                                          Pep_hydrophobicity,
                                          Pep_comp1)
  normalized_Pep_property_table<-Pep_property_table
  normalized_Pep_property_table["Pep_aIndex"]=(normalized_Pep_property_table["Pep_aIndex"]-0)/(278.5714286-0)
  normalized_Pep_property_table["Pep_charge"]=(normalized_Pep_property_table["Pep_charge"]+6.777702387)/(3.215827564+6.777702387)
  normalized_Pep_property_table["PeP_pI"]=(normalized_Pep_property_table["PeP_pI"]-3.370808255)/(12.51646088-3.370808255)
  normalized_Pep_property_table["Pep_InstabilityIndex"]=(normalized_Pep_property_table["Pep_InstabilityIndex"]+71.1875)/(221.7777778+71.1875)
  normalized_Pep_property_table["Pep_hydrophobicity"]=(normalized_Pep_property_table["Pep_hydrophobicity"]+0.848888889)/(0.713333333+0.848888889)
  normalized_Pep_property_table["X1"]=(normalized_Pep_property_table["X1"]-0)/(44-0)
  normalized_Pep_property_table["X2"]=(normalized_Pep_property_table["X2"]-0)/(47-0)
  normalized_Pep_property_table["X3"]=(normalized_Pep_property_table["X3"]-0)/(22-0)
  normalized_Pep_property_table["X4"]=(normalized_Pep_property_table["X4"]-0)/(10-0)
  normalized_Pep_property_table["X5"]=(normalized_Pep_property_table["X5"]-0)/(38-0)
  normalized_Pep_property_table["X6"]=(normalized_Pep_property_table["X6"]-0)/(23-0)
  normalized_Pep_property_table["X7"]=(normalized_Pep_property_table["X7"]-0)/(11-0)
  normalized_Pep_property_table["X8"]=(normalized_Pep_property_table["X8"]-0)/(6-0)
  normalized_Pep_property_table["X9"]=(normalized_Pep_property_table["X9"]-0)/(8-0)
  normalized_Pep_property_table["X10"]=(normalized_Pep_property_table["X10"]-0)/(95.833-0)
  normalized_Pep_property_table["X11"]=(normalized_Pep_property_table["X11"]-0)/(96.774-0)
  normalized_Pep_property_table["X12"]=(normalized_Pep_property_table["X12"]-0)/(87.5-0)
  normalized_Pep_property_table["X13"]=(normalized_Pep_property_table["X13"]-0)/(71.429-0)
  normalized_Pep_property_table["X14"]=(normalized_Pep_property_table["X14"]-9.091)/(96.552-9.091)
  normalized_Pep_property_table["X15"]=(normalized_Pep_property_table["X15"]-3.448)/(90.909-3.448)
  normalized_Pep_property_table["X16"]=(normalized_Pep_property_table["X16"]-0)/(75-0)
  normalized_Pep_property_table["X17"]=(normalized_Pep_property_table["X17"]-0)/(42.857-0)
  normalized_Pep_property_table["X18"]=(normalized_Pep_property_table["X18"]-0)/(50-0)
  return(normalized_Pep_property_table)

}
