/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package run;

import model.secondarystructure.RNASecondaryStructure;
import model.RNAPrimarySequence;
import mfe_optimization.ZukerAlgorithm;

/**
 *
 * @author vot2
 */
public class RunZukerAlgorithm {
    public static void main(String[] args) {
//        String sequence     = "AGGGUU";
//        String dotbracket   = "(...).";

//        String sequence   = "AUUGAGCAUAUUCAC";
//        String dotbracket = "..((((....)))).";

//        String sequence   = "CGAUUGUAGUAGGGUCUCACCUAUUUAGUG";
//        String dotbracket = "..((((.((((((......)))))))))).";
//        String dotbracket = "..((((.((((((......)))))))))).";

        String sequence   = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA";
        String dotbracket = "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))...."; 
        
//        String dotbracket = "((....)(...((((...))))......))";
        
        //double simTime = 1;
        RNASecondaryStructure structure = new RNASecondaryStructure(new RNAPrimarySequence(sequence), dotbracket);
        structure.printStructureInfo();
        structure.secondaryLoopDecomposion();
        
        System.out.println("-----------------------");  
        System.out.println("MFE optimization by Zuker algorithm");        
        //zuker algorithm
        ZukerAlgorithm zuker = new ZukerAlgorithm();
        zuker.config(sequence);
        zuker.evaluateOptimalEnergy();
        double optimalEnergy = zuker.getOptimalEnergyBetweenPosition(0, sequence.length() - 1);
        System.out.println("Minimum energy is: " + optimalEnergy);
        
        String optimalStructure = zuker.evaluateOptimalStructure();;
        System.out.println("Minimum energy is: " + optimalStructure);
    }
}
