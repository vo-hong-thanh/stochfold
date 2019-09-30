/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package simulator.cotranscription;

import constants.SETTINGS;
import java.util.ArrayList;
import java.util.Random;
import model.Action;
import model.ActionList;
import model.RNAPrimarySequence;
import model.secondarystructure.RNASecondaryStructure;
import utility.Cache;
import utility.ComputingMachine;
import utility.DataWriter;

/**
 *
 * @author vot2
 */
public class CoStochFold {
    private Random rand = new Random();

    private Cache<String, ActionList> cacheStructureActions = new Cache<>(SETTINGS.MAX_STOREED_STRUCTURE_WITH_ACTIONS);

    private RNASecondaryStructure structure;
    private ArrayList<String> targetDotBracketList;
    
    private double simTime;    
     
    private double currentTime = 0;
    
    private double transcriptionalInterval;
    private double transcriptionalTime;
    
    private DataWriter dataWriter;
    
    public void config(String _primarySequence, int startingFoldingLength, ArrayList<String> _tergetDotBracketList, double _simTime, double _transcriptionalInterval, DataWriter _dataWriter) {        
        this.structure = new RNASecondaryStructure(new RNAPrimarySequence(_primarySequence), startingFoldingLength);
        this.targetDotBracketList = _tergetDotBracketList;
        this.simTime = _simTime;
        this.transcriptionalInterval = _transcriptionalInterval;
        this.transcriptionalTime = _transcriptionalInterval;
        this.dataWriter = _dataWriter;
    }

    public void runSim() throws Exception{
//        System.out.println("-----------------------------");
        System.out.println("Naive CoStochFold - stochastic co-transcriptional folding of RNA");        
        
        ActionList actionList = null;
        
        long startTimeCounting = System.currentTimeMillis();
        
//        int step = 1;
        do {
//            System.out.println("-------------- STEP = " + step + " ------------------");
//            step++;
//            structure.printStructureInfo();
            
            if(targetDotBracketList.contains(structure.getId())){                
                dataWriter.write(currentTime + "\t" + structure.getId() + "\t" + (targetDotBracketList.contains(structure.getId())));
                break;
            }
            
            //get structure from cacheStructureActions
            if (cacheStructureActions.containsKey(structure.getId())) {
//                System.out.println("Lookup actions from cache");
                actionList = cacheStructureActions.get(structure.getId());
            } else {
                actionList = structure.enumerateActions(false);                
                cacheStructureActions.put(structure.getId(), actionList);
            }

            double partialRateSum = 0;
            double totalRateSum = actionList.getTotalRate();
            
            //generate time
            double delta = ComputingMachine.computeTime(rand, totalRateSum / ComputingMachine.computeRate(structure.getCurrentEnergy()));
            
//            System.out.print("current time: "+ currentTime + " (random) delta t: " + delta);
            
            //update time
            currentTime = currentTime + delta;
            
//            System.out.println(" => updated time: "+ currentTime);
            
            if( (!structure.isFullSequence()) && (currentTime > transcriptionalTime) ){               
                currentTime = transcriptionalTime;
                structure.increaseFoldingLength();                
                
//                System.out.println("@Time " + currentTime + ": transcription event " + " => new structure " + structure.getId() + " (" + structure.getCurrentEnergy() +" kcal)");            

                
                transcriptionalTime += transcriptionalInterval;
                
                continue;
            }
            
            if(currentTime > simTime){                
                dataWriter.write(simTime + "\t" + structure.getId() + "\t" + (targetDotBracketList.contains(structure.getId())));
                break;
            }
            
//           //select action
            Action action = null;            
            double searchValue = rand.nextDouble() * totalRateSum;
            
//            System.out.println("total rate: " + totalRateSum + " =>  (random) search value " + searchValue);
            for (int i = 0; i < actionList.getNumAction(); i++) {
                action = actionList.getAction(i);

                partialRateSum += action.getRate();
//                System.out.println(" partial rate sum " + partialRateSum);
                if (partialRateSum >= searchValue) {
//                    System.out.println(" select action " + action.toReducedString());
                    break;
                }
            }

            //update action
            structure.applyAction(action);

//            System.out.println("@Time " +currentTime + ": folding event " + action.toReducedString() + " => new structure " + structure.getId() + " (" + structure.getCurrentEnergy() +" kcal)");
            
        } while (currentTime <= simTime);
        
        long endTimeCounting = System.currentTimeMillis();
        long cpuTime = endTimeCounting - startTimeCounting;
        
        dataWriter.writeLine("\t" + cpuTime);
    }
}
