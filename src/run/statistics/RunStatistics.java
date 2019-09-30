/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package run.statistics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 *
 * @author vot2
 */
public class RunStatistics {
    public static void main(String[] args) throws IOException {
        System.out.println("[Current working directory: " + (new File(".")).getCanonicalPath() + "]");

        BufferedReader screenReader = new BufferedReader(new InputStreamReader(System.in));
        
        String outputStatisticFile;
        
        double time;
        double interval;       
        
        
        ArrayList<String> targetDotBracketList = new ArrayList<>();

//        String targetDotBracket = "..((((....)))).";    
//        targetDotBracketList.add(targetDotBracket);
        

////        //switching molecule - 33 nt
////        String sequence           = "GGCCCCUUUGGGGGCCAGACCCCUAAAGGGGUC";
////        String startingDotBracket = ".................................";        
//        String optimalDotBracket  = "((((((((((((((.....))))))))))))))"; //optimal mfe       
//        String suboptimalDotBracket = "((((((....)))))).((((((....))))))"; //suboptimal mfe   
//        targetDotBracketList.add(optimalDotBracket);
//        targetDotBracketList.add(suboptimalDotBracket);
        
////        //spliced leader from Leptomonas collosoma - 56 nt 
////        String sequence             = "AACUAAAACAAUUUUUGAAGAACAGUUUCUGUACUUCAUUGGUAUGUAGAGACUUC";
////        String startingDotBracket   = "........................................................";  
//        String targetDotBracket1    = ".......................((((((((((((.....)))))..))))))).."; 
//        String targetDotBracket2    = "..((...((((((..(((((.((((...)))).)))))..))).)))..))....."; 
//        targetDotBracketList.add(targetDotBracket1);
//        targetDotBracketList.add(targetDotBracket2);  

//        //tRNA - 76nt
//        String sequence              = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA";
//        String startingDotBracket    = "............................................................................";
        String suboptimalDotBracket  = "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))...."; 
        String optimalDotBracket     = "((((((((.(......((((.((((((..((((...........))))..))))))..))))).))))))))....";
//        String optimalDotBracket   = ".((((((....((((((.((......)).)))))).))))))(((.((.((((.......)))))).)))......";
        targetDotBracketList.add(suboptimalDotBracket);        
        targetDotBracketList.add(optimalDotBracket);      

//        //QBeta - 115nt
////        String sequence             = "GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA";
////        String startingDotBracket   = "...................................................................................................................";
//        String optimalmfe           = "(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).";
//        String suboptimal           = "(((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....)))))))).....";
//        targetDotBracketList.add(optimalmfe);
//        targetDotBracketList.add(suboptimal);  
                
        int numBin;
        int numTarget = targetDotBracketList.size();
        CountingBinData[][] bins;
        
        System.out.print("Enter file name: ");
        outputStatisticFile = screenReader.readLine().trim();
        
        System.out.print("Enter time: ");
        time = Double.parseDouble(screenReader.readLine().trim());
        
        System.out.print("Enter interval: ");
        interval = Double.parseDouble(screenReader.readLine().trim());
        
        numBin = (int)(time / interval);
        bins = new CountingBinData[numTarget][numBin];
        for (int t = 0; t < numTarget; t++) {
            for (int i = 0; i < numBin; i++) {
                bins[t][i] = new CountingBinData();
            }
        }
        
        BufferedReader datafileReader = new BufferedReader(new FileReader(outputStatisticFile));
        String dataLine;
        StringTokenizer tokens;
        
        //skip fire line
        dataLine = datafileReader.readLine();
        
        //read data
        double totalRuntime = 0;
        
        //read data
        while( (dataLine = datafileReader.readLine()) != null){
            if(dataLine.trim().equals(""))
                continue;
            tokens = new StringTokenizer(dataLine);
            
            //each line: folding time - predicted structure - found - runtime
            double foldingTime = Double.parseDouble(tokens.nextToken());
            String predictedStructure = tokens.nextToken();
            boolean found = Boolean.parseBoolean(tokens.nextToken());
            double runTime = Double.parseDouble(tokens.nextToken());
                        
            if(found){
                for (int t = 0; t < numTarget; t++) {
                    if (targetDotBracketList.get(t).equals(predictedStructure) ) {
                        int binId = (int) (foldingTime / interval);
                        bins[t][binId].addItem(runTime);
                    }
                }
            }
            
            totalRuntime += runTime;
        }

        System.out.println("Total runtime: " + totalRuntime);
        System.out.println("Total bins: " + numBin); 
        for (int t = 0; t < numTarget; t++) {
            System.out.println("Target structure: " + targetDotBracketList.get(t));
            for(int i = 0; i < numBin; i++){
                System.out.println("bin " + i + " (time < " + ((i+1)*interval) +") : " + bins[t][i].getCount() + " items"); 
            }
        }
    }
}
