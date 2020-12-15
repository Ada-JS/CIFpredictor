#' @title Predict isopropanol elution concentration
#' @description Predict isopropanol elution concentration for a given peptide sequence
#' @return The output from \code{\link{predict_iso_conc}}
#' @export
#' @param pep_property The input peptide property calculated by function pep_property()
#' @examples
#' seq_property<-pep_property("AAAAAAAAAA")
#' predict_iso_conc(seq_property)

predict_iso_conc<-function(pep_property){

  pep_property[is.na(pep_property)]<-0

  conc=76.61610175+
    pep_property[1,"Pep_aIndex"]*5.271212572+
    pep_property[1,"Pep_charge"]*11.86668818+
    pep_property[1,"PeP_pI"]*0-
    pep_property[1,"Pep_InstabilityIndex"]*0.185983347+
    pep_property[1,"Pep_hydrophobicity"]*8.296614231+
    pep_property[1,"X1"]*0+
    pep_property[1,"X2"]*0+
    pep_property[1,"X3"]*0+
    pep_property[1,"X4"]*0+
    pep_property[1,"X5"]*0-
    pep_property[1,"X6"]*6.819724125-
    pep_property[1,"X7"]*21.04899333+
    pep_property[1,"X8"]*0+
    pep_property[1,"X9"]*0-
    pep_property[1,"X10"]*3.828794549+
    pep_property[1,"X11"]*0+
    pep_property[1,"X12"]*0+
    pep_property[1,"X13"]*0+
    pep_property[1,"X14"]*0.846420227+
    pep_property[1,"X15"]*0+
    pep_property[1,"X16"]*0+
    pep_property[1,"X17"]*0-
    pep_property[1,"X18"]*5.336390334

return(conc)

}
