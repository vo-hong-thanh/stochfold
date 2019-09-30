/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package simulator.openstrand;

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
public class StochFold {
    private Random rand = new Random();

    private Cache<String, ActionList> cacheStructureActions = new Cache<>(SETTINGS.MAX_STOREED_STRUCTURE_WITH_ACTIONS);

    private RNASecondaryStructure structure;
    private double simTime;
    private ArrayList<String> targetDotBracketList;
    
    private DataWriter dataWriter;

    public void config(String _primarySequence, String _startingStructure, ArrayList<String> _tergetDotBracketList, double _simTime, DataWriter _dataWriter) {        
        this.structure = new RNASecondaryStructure(new RNAPrimarySequence(_primarySequence), _startingStructure);
        this.targetDotBracketList = _tergetDotBracketList;
        this.simTime = _simTime;
        this.dataWriter = _dataWriter;
    }

    public void runSim() throws Exception{
//        System.out.println("-----------------------------");
        System.out.println("StochFold - stochastic folding of RNA");
        System.out.println("--- Initial condition ---");
        structure.printStructureInfo();   
        
        double currentTime = 0;
        ActionList actionList = null;
        
        long timeComputeNextStructures = 0;
        long timeSelectNextStructure = 0;
        long totalRuntime = 0;
        
        System.out.println("--- Simulation run ---");
        long startRuntimeCounting = System.currentTimeMillis();

//        long step = 0;
        do {
//            step++;
//            System.out.println("-------------- STEP = " + step + " ------------------");
//            structure.printStructureInfo();
            
            if(targetDotBracketList.contains(structure.getId())){                
                dataWriter.write(currentTime + "\t" + structure.getId() + "\t" + (targetDotBracketList.contains(structure.getId())));
                break;
            }
                        
            long startComputeNextStructures = System.currentTimeMillis();
            //get structure from cacheStructureActions
            if (cacheStructureActions.containsKey(structure.getId())) {
//                System.out.println("Lookup actions from cache");
                actionList = cacheStructureActions.get(structure.getId());
            } else {
                actionList = structure.enumerateActions(false);
                cacheStructureActions.put(structure.getId(), actionList);
            }
            long endComputeNextStructures = System.currentTimeMillis();
            timeComputeNextStructures += endComputeNextStructures - startComputeNextStructures;
            
            double totalRateSum = actionList.getTotalRate();
            
            //generate time
            double delta = ComputingMachine.computeTime( rand, totalRateSum / ComputingMachine.computeRate(structure.getCurrentEnergy()) );

//            System.out.print("current time: "+ currentTime + "(random) delta t: " + delta);
            
            //update time
            currentTime = currentTime + delta;
            
//            System.out.println(" => updated time: "+ currentTime);
            
            if(currentTime > simTime){                
                dataWriter.write(simTime + "\t" + structure.getId() + "\t" + (targetDotBracketList.contains(structure.getId())));
                break;
            }
            
            long startSelectNextStructure = System.currentTimeMillis(); 
            //select action
            Action action = null;

            double partialRateSum = 0;            
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
            long endSelectNextStructure = System.currentTimeMillis(); 
            timeSelectNextStructure += endSelectNextStructure - startSelectNextStructure;
            
            //update action
            structure.applyAction(action);
            
//            System.out.println("@Time " +currentTime + ": folding event " + action.toString() + " => " + structure.getId() + " (" + structure.getCurrentEnergy() +" kcal)");
            
        } while (currentTime <= simTime);        
        long endRuntimeCounting = System.currentTimeMillis();
        totalRuntime = endRuntimeCounting - startRuntimeCounting;
        
        dataWriter.writeLine("\t" + totalRuntime + "\t" + timeComputeNextStructures + "\t" + timeSelectNextStructure);
        System.out.println("-----------------------------");
    }
}
