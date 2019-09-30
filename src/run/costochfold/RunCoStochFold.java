/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package run.costochfold;

import run.statistics.CountingBinData;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;
import simulator.cotranscription.CoStochFold;
import utility.DataWriter;

/**
 *
 * @author vot2
 */
public class RunCoStochFold {
    public static void main(String[] args) throws Exception{
        System.out.println("[Current working directory: " + (new File(".")).getCanonicalPath() + "]");

        BufferedReader screenReader = new BufferedReader(new InputStreamReader(System.in));
        
        ArrayList<String> targetDotBracketList = new ArrayList<>();

        StringBuilder outputDataFile = new StringBuilder("CoStochFold_");
        
//        //simple hairpin - 15 nt
//        String sequence           = "ACUGAUCGUAGUCAC"; //"AUUGAGCAUAUUCAC";  /*efficient sequence*/ // "ACUGAUCGUAGUCAC"; /*inefficient sequence*/
//        String startingDotBracket = "...............";
//        String targetDotBracket   = "..((((....)))).";    
//        targetDotBracketList.add(targetDotBracket);
//        outputDataFile.append("Hairpin");

//        //small 20 nt 
//        String sequence             = "UUGCUAAGCAACCAUUGGUU";
//        String startingDotBracket   = "....................";  
//        String targetDotBracket1    = "..(((((.......)))))."; 
//        String targetDotBracket2    = ".........((((...))))"; 
//        String targetDotBracket3    = "((((...))))........."; 
//        targetDotBracketList.add(targetDotBracket1);
//        targetDotBracketList.add(targetDotBracket2);  
//        targetDotBracketList.add(targetDotBracket3);
//        outputDataFile.append("20nt");

//        //switching molecule - 33 nt
//        String sequence           = "GGCCCCUUUGGGGGCCAGACCCCUAAAGGGGUC";
//        String startingDotBracket = ".................................";        
//        String optimalDotBracket  = "((((((((((((((.....))))))))))))))"; //optimal mfe       
//        String suboptimalDotBracket = "((((((....)))))).((((((....))))))"; //suboptimal mfe   
//        targetDotBracketList.add(optimalDotBracket);
//        targetDotBracketList.add(suboptimalDotBracket);
//        outputDataFile.append("Switching");

//        //spliced leader from Leptomonas collosoma - 56 nt 
//        String sequence             = "AACUAAAACAAUUUUUGAAGAACAGUUUCUGUACUUCAUUGGUAUGUAGAGACUUC";
//        String startingDotBracket   = "........................................................";  
//        String targetDotBracket1    = ".......................((((((((((((.....)))))..))))))).."; 
//        String targetDotBracket2    = "..((...((((((..(((((.((((...)))).)))))..))).)))..))....."; 
//        targetDotBracketList.add(targetDotBracket1);
//        targetDotBracketList.add(targetDotBracket2);  
//        outputDataFile.append("spliced-leader");
       
        //tRNA - 76nt
        String sequence              = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA";
        String startingDotBracket    = "............................................................................";
        String suboptimalDotBracket  = "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))...."; 
        String optimalDotBracket     = "((((((((.(......((((.((((((..((((...........))))..))))))..))))).))))))))....";
//        String optimalDotBracket   = ".((((((....((((((.((......)).)))))).))))))(((.((.((((.......)))))).)))......";
        targetDotBracketList.add(suboptimalDotBracket);        
        targetDotBracketList.add(optimalDotBracket);              
        outputDataFile.append("tRNA");         

//        //QBeta - 115nt
//        String sequence             = "GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA";
//        String startingDotBracket   = "...................................................................................................................";
//        String optimalmfe           = "(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..))).";
//        String suboptimal           = "(((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....)))))))).....";
//        targetDotBracketList.add(optimalmfe);
//        targetDotBracketList.add(suboptimal);  
//        outputDataFile.append("QBeta");
        
        int startingFoldingLength = 1;
        
        System.out.print("Enter transcriptional time: ");        
        double transcriptionalInterval = Double.parseDouble(screenReader.readLine().trim());// 5; //1; //0.2; //0.02; //0.005; // s
        
        ////////////////////////////////////
        //run simulation
        double simTime = 10000;
        
        int numRuns = 10000;
        
        outputDataFile.append("-" + transcriptionalInterval + "_" + numRuns + ".txt");
        DataWriter dataWriter = new DataWriter(outputDataFile.toString());
        dataWriter.writeLine("Time\tPrediction\tFound\tRuntime");
        
        CoStochFold coFold;
        for(int i = 0; i < numRuns; i++){       
            System.out.print("run @" + i  +": ");
            coFold = new CoStochFold();
            coFold.config(sequence, startingFoldingLength, targetDotBracketList, simTime, transcriptionalInterval, dataWriter);
            coFold.runSim();        
        }
        
        dataWriter.flush();
        dataWriter.close();    
        
        ///////////////////////////////////////
        //run statistics
        double interval = 10;

        int numBin;
        int numTarget = targetDotBracketList.size();
        CountingBinData[][] bins;

        numBin = (int) (simTime / interval);
        bins = new CountingBinData[numTarget][numBin];
        for (int t = 0; t < numTarget; t++) {
            for (int i = 0; i < numBin; i++) {
                bins[t][i] = new CountingBinData();
            }
        }

        BufferedReader datafileReader = new BufferedReader(new FileReader(outputDataFile.toString()));
        String dataLine;
        StringTokenizer tokens;

        //skip fire line
        dataLine = datafileReader.readLine();
        //read data
        double totalRuntime = 0;
        while( (dataLine = datafileReader.readLine()) != null){
            if(dataLine.trim().equals(""))
                continue;
            tokens = new StringTokenizer(dataLine);

            //each line: folding time - predicted structure - expected structure - found - runtime
            double foldingTime = Double.parseDouble(tokens.nextToken());
            String predictedStructure = tokens.nextToken();
            boolean found = Boolean.parseBoolean(tokens.nextToken());
            double runTime = Double.parseDouble(tokens.nextToken());

            if(found){
                for (int t = 0; t < targetDotBracketList.size(); t++) {
                    if (targetDotBracketList.get(t).equals(predictedStructure) ) {
                        int binId = (int) (foldingTime / interval);
                        bins[t][binId].addItem(runTime);
                    }
                }
            }
            totalRuntime += runTime;
        }

        String outputStatisticsFile = "(Statistics)" + outputDataFile;
        DataWriter statsWriter = new DataWriter(outputStatisticsFile);

        statsWriter.writeLine("Total run: " + numRuns);
        statsWriter.writeLine("Total runtime: " + totalRuntime);
        statsWriter.writeLine("Total bins: " + numBin);

        for (int t = 0; t < targetDotBracketList.size(); t++) {
            statsWriter.writeLine("Target structure: " + targetDotBracketList.get(t));
            for (int i = 0; i < numBin; i++) {
                statsWriter.writeLine("bin " + i + " (time < " + ((i + 1) * interval) + ") : " + bins[t][i].getCount());
            }
        }
        
        statsWriter.flush();
        statsWriter.close();
    }
}
