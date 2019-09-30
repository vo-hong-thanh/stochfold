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
public class DelayedStochFold {
    private Random rand = new Random();

    private Cache<String, ActionList> cacheStructureActions = new Cache<>(SETTINGS.MAX_STOREED_STRUCTURE_WITH_ACTIONS);
    private Cache<String, Double> cacheStructureEnergys = new Cache<>(SETTINGS.MAX_STORED_STRUCTURE_WITH_ENEGRY);
    
    private RNASecondaryStructure structure;
    private ArrayList<String> tergetDotBracketList;
    private double simTime;
    
    private DataWriter dataWriter;

    public void config(String _primarySequence, String _startingStructure, ArrayList<String> _tergetDotBracketList , double _simTime, DataWriter _dataWriter) {        
        this.structure = new RNASecondaryStructure(new RNAPrimarySequence(_primarySequence), _startingStructure);
        structure.setOptimalEnergy();
        structure.initializeEnergyCache(cacheStructureEnergys);        
        this.tergetDotBracketList = _tergetDotBracketList;
        this.simTime = _simTime;       
        this.dataWriter = _dataWriter;
    }

    public void runSim() throws Exception{
//        System.out.println("-----------------------------");
        System.out.println("DelayedStochFold - efficeint stochastic folding of RNA");
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
        long trial = 0;
        long acceptedStep = 0;
    
        do {
//            step++;
//            System.out.println("-------------- STEP = " + step + " ------------------");
//            structure.printStructureInfo();
            
            if(tergetDotBracketList.contains(structure.getId())){                
                dataWriter.write(currentTime + "\t" + structure.getId() + "\t" + (tergetDotBracketList.contains(structure.getId())));
                break;
            }
            
            long startComputeNextStructures = System.currentTimeMillis();
            //get structure from cacheStructureActions
            if (cacheStructureActions.containsKey(structure.getId())) {
//                System.out.println("Lookup actions from cache");
                actionList = cacheStructureActions.get(structure.getId());
            } else {
                actionList = structure.enumerateActions(true);
                cacheStructureActions.put(structure.getId(), actionList);
            }
            long endComputeNextStructures = System.currentTimeMillis();
            timeComputeNextStructures += endComputeNextStructures - startComputeNextStructures;
           
            double searchValue;
            double partialRateSum = 0;            
            
            int numActions = actionList.getNumAction();
            int actionIndex = -1;
            Action action = null;            
            boolean accept = false;
            
            double currentStructureRate = ComputingMachine.computeRate(structure.getCurrentEnergy());
            do{
                trial++;
                double totalOverestimatedRateSum = actionList.getTotalRate();

                //generate time
                double delta = ComputingMachine.computeTime(rand, totalOverestimatedRateSum / currentStructureRate );

    //            System.out.print("current time: "+ currentTime + "(random) delta t: " + delta);

                //update time
                currentTime = currentTime + delta;

    //            System.out.println(" => updated time: "+ currentTime);

                if(currentTime > simTime){                
                    dataWriter.write(simTime + "\t" + structure.getId() + "\t" + (tergetDotBracketList.contains(structure.getId())));
                    break;
                }

                long startSelectNextStructure = System.currentTimeMillis(); 

                //select action                
                searchValue = rand.nextDouble() * totalOverestimatedRateSum;
                double overestimatedRate = 0;                
//                System.out.println("total overestimatedRate: " + totalOverestimatedRateSum + " =>  (random) search value " + searchValue);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////       Linear Search          //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                partialRateSum = 0; 
//                for (actionIndex = 0; actionIndex < numActions; actionIndex++) {
//                    action = actionList.getAction(actionIndex);
//                    overestimatedRate = action.getRate();
//                    partialRateSum += overestimatedRate;
//    //                System.out.println(" partial overestimatedRate sum " + partialRateSum);
//                    if (partialRateSum >= searchValue) {
//    //                    System.out.println(" select action " + action.toReducedString());
//
//                        //rejection-based test
//                        if(action.checkIsEnergyComputed()){
//    //                        System.out.println(" accepted because has energy is already computed");
//                            accept = true;
//                        }
//                        else{
//                            double energy = structure.evaluateEnergyByAction(action);
//
//                            action.setEnergy(energy);
//                            double trueRate = action.getRate();
//
//    //                        System.out.println(" energy = " + energy + " => rate = " + trueRate);
//
//                            double randValue = rand.nextDouble();
//
//    //                        System.out.println(" random value = " + randValue + ", current estmated rate " + overestimatedRate);
//
//                            if(trueRate >= randValue*overestimatedRate){
//    //                            System.out.println(" => accepted by rejection test");
//                                accept = true;
//                            }
//                            actionList.updateTotalRate(trueRate - overestimatedRate);
//                        }
//
//                        //update action
//                        if(accept){
//                            acceptedStep++;  
//                            structure.applyAction(action);
//            //                System.out.println("@Time " + currentTime + ": folding event " + action.toString() + " => " + structure.getId() + " (" + structure.getCurrentEnergy() +" kcal)");
//                        }            
//            //            else{
//            //                System.out.println(" a rejection");
//            //            }
//                        break;
//                    }
//                }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////       Cache-based search          /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////                
                if(partialRateSum < searchValue){
                    //sum-up
//                    System.out.println("partialRateSum = " + partialRateSum  + " < searchValue => sum-up method from current actionIndex = " + actionIndex);
                    for(actionIndex++; actionIndex < numActions; actionIndex++){
                        action = actionList.getAction(actionIndex);
                        overestimatedRate = action.getRate();
                        partialRateSum += overestimatedRate;                    

                        if(partialRateSum >= searchValue){
//                            System.out.println(" => found actionIndex  = " + actionIndex);
                            break;
                        }
                    }
                }else{
                    //chop-down
//                    System.out.println("partialRateSum = " + partialRateSum  + " >= searchValue => chop-down method from current actionIndex = " + actionIndex);
                    for(; actionIndex >= 0; actionIndex--){
                        action = actionList.getAction(actionIndex);
                        overestimatedRate = action.getRate();                           
                        if(partialRateSum - overestimatedRate < searchValue){
//                            System.out.println(" => found actionIndex  = " + actionIndex);
                            break;
                        }
                        partialRateSum -= overestimatedRate;
                    }
                }

                //rejection-based test
                if(action.checkIsEnergyComputed()){
//                    System.out.println(" accepted because has energy is already computed");
                    accept = true;
                }
                else{
                    double energy = structure.evaluateEnergyByAction(action);
                    action.setEnergy(energy);
                    double trueRate = action.getRate();

//                    System.out.println(" energy = " + energy + " => rate = " + trueRate);

                    double randValue = rand.nextDouble();

//                    System.out.println(" random value = " + randValue + ", current estmated rate " + overestimatedRate);

                    if(trueRate >= randValue*overestimatedRate){
//                        System.out.println(" => accepted by rejection test");
                        accept = true;
                    }
                    double diff = trueRate - overestimatedRate;
                    
                    actionList.updateTotalRate(trueRate - overestimatedRate);
                    partialRateSum += diff;
                }

                //update action
                if(accept){
                    acceptedStep++;
                    structure.applyAction(action);
//                    System.out.println("@Time " + currentTime + ": folding event " + action.toString() + " => " + structure.getId() + " (" + structure.getCurrentEnergy() +" kcal)");
                }            
//                else{
//                    System.out.println(" a rejection");
//                }                       
               
                long endSelectNextStructure = System.currentTimeMillis(); 
                timeSelectNextStructure += endSelectNextStructure - startSelectNextStructure;
            }while(!accept);            
        } while (currentTime <= simTime);        
        long endRuntimeCounting = System.currentTimeMillis();
        totalRuntime = endRuntimeCounting - startRuntimeCounting;
        
        dataWriter.writeLine("\t" + totalRuntime + "\t" + timeComputeNextStructures + "\t" + timeSelectNextStructure + "\t" + trial + "\t" + acceptedStep);
        System.out.println("-----------------------------");
    }
}
