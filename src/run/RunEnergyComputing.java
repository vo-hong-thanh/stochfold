/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package run;

import constants.ENERGY_CONSTANTS;
import model.RNAPrimarySequence;
import model.secondarystructure.RNASecondaryStructure;

/**
 *
 * @author vot2
 */
public class RunEnergyComputing {
    public static void main(String[] args) {
        String sequence;     
        String dotbracket;
        
//        sequence =   "CGAUUGUAGUAGGGUCUCACCUAUUUAGUG";
////        dotbracket = "..((((.((((((......)))))))))).";
//        dotbracket = "...(.(..(..............).).)..";
        
//        dotbracket = "..((((...((((......))))..)))).";
//        dotbracket = "((.....((((((......))))))...))";
        
//        sequence     = "AACUAAAACAAUUUUUGAAGAACAGUUUCUGUACUUCAUUGGUAUGUAGAGACUUC";
//        dotbracket   = ".......................((((((((((((.....)))))..))))))).."; 
        
//        sequence     = "GUAGUAGGGUCU";
//        dotbracket   = "((...)(...))";
        
//        sequence     = "GGGGGGCCCCCC";
//        dotbracket   = "...((...))..";  

/////////////////////////////////////////////
//        hairpin example
/////////////////////////////////////////////
//        sequence   = "CACAAAAAAAUGUG";
//        dotbracket = "((((......))))"; 
        
//        sequence   = "CACAGGAAAUGUG";
//        dotbracket = "((((.....))))"; 

//        sequence   = "CACCGAAAGGUG";
//        dotbracket = "((((....))))"; 
        
//        sequence   = "CACACCCCCCUGUG";
//        dotbracket = "((((......))))"; 
        
//        sequence   = "CGGGGGAAGUCCG";
//        dotbracket = "((((.....))))"; 

/////////////////////////////////////////////        
        //bulge
/////////////////////////////////////////////
//        sequence   = "GCCCGAAACGGC";
//        dotbracket = "(((.(...))))"; 
        
//        sequence   = "GAACAGAAACUC";
//        dotbracket = "((...(...)))";
/////////////////////////////////////////////    

        //internal
/////////////////////////////////////////////
//        sequence   = "CAGACGAAACGGAGUG";
//        dotbracket = "((..((...))...))"; 
  
//        sequence   = "CAGCGAAACGGAAAGUG";
//        dotbracket = "((.((...)).....))"; 
        
//        sequence   = "CAGACGAAACGGAUG";
//        dotbracket = "((..((...))..))"; 
        
//        sequence   = "AAGCGAAACGGUUA";
//        dotbracket = "((.((...)).))."; 
        
               
        //tRNA - 76nt
        sequence              = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA";
        dotbracket     = "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))...."; 


/////////////////////////////////////////////        
        StringBuffer currDotBracketForm = new StringBuffer(dotbracket.length());
        for (int i = 0; i < sequence.length(); i++) {
            currDotBracketForm.append('.');
        }

        //initialize pairing list
        int[] pairList = new int[sequence.length()];
        //build ring list
        for (int i = 0; i < sequence.length(); i++) {
            pairList[i] = -1;
        }

        //primary sequence
        RNAPrimarySequence primarySequence = new RNAPrimarySequence(sequence);
        
        //build pair list
        StringBuffer copyDotBracketForm = new StringBuffer(dotbracket);
        int ipos, jpos, balance = 0;
        for (ipos = 0; ipos < copyDotBracketForm.length(); ipos++) {
            if (copyDotBracketForm.charAt(ipos) == ')') {
                jpos = ipos;
                copyDotBracketForm.setCharAt(ipos, '.');
                balance++;

                while (true) {
                    ipos--;
                    if (ipos < 0) {
                        break;
                    }

                    if (copyDotBracketForm.charAt(ipos) == '(') {
                        if(!primarySequence.isPairable(ipos, jpos)){
                            throw new RuntimeException(primarySequence.getNucleotide(ipos) + " and " + primarySequence.getNucleotide(jpos) + " do not pair as Watson-Crick or Wooble base-paired rules");
                        }
                        
                        if (jpos - ipos <= ENERGY_CONSTANTS.PAIRING_CONSTRAINT) {
                            throw new RuntimeException("The gap between a base-pair must greater than " + ENERGY_CONSTANTS.PAIRING_CONSTRAINT);
                        }

                        balance--;
                        copyDotBracketForm.setCharAt(ipos, '.');

                        pairList[ipos] = jpos;
                        pairList[jpos] = ipos;

                        ipos = jpos;
                        break;
                    }
                }
            }
        }
        if (balance != 0) {
            throw new RuntimeException("Incorrect dot-bracket form");
        }

        System.out.println("Primary sequence: " + sequence);
        System.out.println("Starting conformation: " + dotbracket);
//        System.out.println("Pairing list");
//        for (int i = 0; i < pairList.length; i++) {
//            if(pairList[i] != -1){
//                System.out.println("Position " + i + " is pairing with " + pairList[i]);
//            }
//        }
               
        //compute energy
        RNASecondaryStructure conformation = new RNASecondaryStructure(new RNAPrimarySequence(sequence), dotbracket);
//        //conformation.enumerateActions();
//        double energy = conformation.evaluateStructureEnergy();
//        
        System.out.println("Energy of structure " + conformation.getCurrentEnergy() + " kcal");
        
//        conformation.secondaryLoopDecomposion();
    }
}
