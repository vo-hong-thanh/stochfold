/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package model.secondarystructure;

import constants.ACTION_CODE;
import constants.ENERGY_CONSTANTS;
import constants.RINGNODE_CODE;
import java.util.ArrayList;
import constants.LOOP_TYPE;
import model.Action;
import model.ActionList;
import model.RNAPrimarySequence;
import utility.Cache;
import utility.ComputingMachine;

/**
 *
 * @author vot2
 */
public class RNASecondaryStructure {
    private RNAPrimarySequence primarySequence;

    private int currentFoldingLength;
    private StringBuffer currentDotBracketForm;
    private double currentStructureEnergy;

    private int[] pairList;

    private RingNode[] ringList;
    private RingNode root;

    private double optimalEnergy;
    
    private Cache<String, Double> cacheStructureEnergys = null;

    public RNASecondaryStructure(RNAPrimarySequence _primarySequence, String _currDotBracketForm) {
        //initialize information
        initializeStructure(_primarySequence);

        //dot-bracket of secondary structure
        currentFoldingLength = primarySequence.getSequenceLength();
        currentDotBracketForm = new StringBuffer(currentFoldingLength);
        for (int i = 0; i < currentFoldingLength; i++) {
            currentDotBracketForm.append('.');
        }

        //build ring list
        buildRingList();

        //build internal representation
        buildTree(_currDotBracketForm);

        //compute energy of structure
        currentStructureEnergy = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);

//        //print info        
//        printStructureInfo();         
    }

    public RNASecondaryStructure(RNAPrimarySequence _primarySequence, int _startingFoldingLength) {
        //initialize information
        initializeStructure(_primarySequence);

        //dot-bracket of secondary structure
        currentFoldingLength = _startingFoldingLength;
        currentDotBracketForm = new StringBuffer(currentFoldingLength);
        for (int i = 0; i < currentFoldingLength; i++) {
            currentDotBracketForm.append('.');
        }

        //build ring list
        buildRingList();

        //compute energy of structure
        currentStructureEnergy = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);

//        //printInfo
//        printStructureInfo();          
    }    

    private void initializeStructure(RNAPrimarySequence _primarySequence) {
        //primary squence
        primarySequence = _primarySequence;

        //pairing list
        pairList = new int[primarySequence.getSequenceLength() + 1];
        for (int i = 0; i <= primarySequence.getSequenceLength(); i++) {
            pairList[i] = -1;
        }

        //ring list
        ringList = new RingNode[primarySequence.getSequenceLength() + 1];
        for (int i = 0; i <= primarySequence.getSequenceLength(); i++) {
            ringList[i] = new RingNode();
        }

        //sentinel nodes
        //root
        root = new RingNode();
        root.setPosition(-1);
        root.setType(RINGNODE_CODE.ROOT);
        root.setDown(ringList[primarySequence.getSequenceLength()]);

        //last element (dummy node)
        ringList[primarySequence.getSequenceLength()].setPosition(primarySequence.getSequenceLength());
        ringList[primarySequence.getSequenceLength()].setNext(ringList[0]);
        ringList[primarySequence.getSequenceLength()].setUp(root);
        ringList[primarySequence.getSequenceLength()].setType(RINGNODE_CODE.UNKNOWN);
    }

    private void buildRingList() {
        //build ring elements
        for (int i = 0; i < currentFoldingLength; i++) {
            ringList[i].setPosition(i);
            ringList[i].setType(RINGNODE_CODE.UNPAIRED);

            //link to previous node
            if (i == 0) {
                ringList[i].setPrev(ringList[primarySequence.getSequenceLength()]);
            } else {
                ringList[i].setPrev(ringList[i - 1]);
            }

            //link to next node
            if (i == currentFoldingLength - 1) {
                ringList[i].setNext(ringList[primarySequence.getSequenceLength()]);
                //set link of the last node   
                ringList[primarySequence.getSequenceLength()].setPrev(ringList[i]);
            } else {
                ringList[i].setNext(ringList[i + 1]);
            }
        }
    }

    private void buildTree(String dotBracketForm) {
        StringBuffer copyDotBracketForm = new StringBuffer(dotBracketForm);
        int ipos, jpos, balance = 0;

        for (ipos = 0; ipos < dotBracketForm.length(); ipos++) {
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
                        if (!primarySequence.isPairable(ipos, jpos)) {
                            throw new RuntimeException(primarySequence.getNucleotide(ipos) + " and " + primarySequence.getNucleotide(jpos) + " do not pair as Watson-Crick or Wooble base-paired rules");
                        }

                        if (jpos - ipos <= ENERGY_CONSTANTS.PAIRING_CONSTRAINT) {
                            throw new RuntimeException("The gap between a base-pair must be greater than " + ENERGY_CONSTANTS.PAIRING_CONSTRAINT);
                        }

                        balance--;
                        copyDotBracketForm.setCharAt(ipos, '.');
                        do_close_bp(ringList[ipos], ringList[jpos]);
                        ipos = jpos;
                        break;
                    }
                }
            }
        }

        if (balance != 0) {
            throw new RuntimeException("Incorrect dot-bracket form");
        }
    }

    public boolean isFullSequence() {
        return currentFoldingLength == primarySequence.getSequenceLength();
    }

    public String getId() {
        return currentDotBracketForm.toString();
    }

    public double getCurrentEnergy() {
        return currentStructureEnergy;
    }

    @Override
    public int hashCode() {
        return currentDotBracketForm.hashCode();
    }

    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }
        if (!(o instanceof RNASecondaryStructure)) {
            return false;
        }

        RNASecondaryStructure c = (RNASecondaryStructure) o;
        return c.currentDotBracketForm.equals(this.currentDotBracketForm);
    }

    public void initializeEnergyCache(Cache<String, Double> _cacheStructureEnergys) {
        cacheStructureEnergys = _cacheStructureEnergys;
    }
    
    private boolean isEnergyCached(String structureId ) {
        if(cacheStructureEnergys != null && cacheStructureEnergys.containsKey(structureId))
            return true;
        return false;
    }

    private double lookupEnergyCache(String structureId) {
        return cacheStructureEnergys.get(structureId);
    }
    
    private double getEnergyofStructure() {
        double energy = 0.0;
        energy = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);
        //put in cache
        cacheStructureEnergys.put(getId(), energy);
          
        return energy;
    }    
    
    public void setOptimalEnergy(){        
//        //get optimal energy
        optimalEnergy = ComputingMachine.computeOptimalEnergy(primarySequence);
        
        System.out.println("Optimal energy: " + optimalEnergy + " => rate " + ComputingMachine.computeRate(optimalEnergy));
    }
    
    public void printStructureInfo() {
//        System.out.println("-------------------");
        System.out.println("Structure info:");
        System.out.println(" + Current folding length: " + currentFoldingLength);
        System.out.println(" + Nucleotide sequence: " + primarySequence.getFoldingSequence(currentFoldingLength));
        System.out.println(" + Current conformation: " + currentDotBracketForm.toString());        
        System.out.println(" + Energy of structure " + currentStructureEnergy + " (kcal)");
//        System.out.println("-------------------");

//        //internal representation
//        System.out.println("Internal representation:");
//        System.out.println("Root " + root.toString());
//                
//        System.out.println("Ringlist");
//        for(int i = 0; i < ringList.length; i++){
//            if(ringList[i].getPosition() == -1) {
//                System.out.println("Position " + i + " : Skip (not yet fold to this position)");  
//                continue;
//            }
//            System.out.println(ringList[i]);
//            
//            if(i == primarySequence.getSequenceLength()){
//                System.out.println(" - Last node is MULTILOOP_A dummy");
//                continue;
//            }            
//
//            if(pairList[i] != -1){
//                System.out.println(" - Pairing with node at position " + pairList[i]);   
//            }else{
//                System.out.println(" - Not pairing yet");   
//            }
//        }
    }

    public void increaseFoldingLength() {
//        System.out.println("[Function - Increase folding length]");

//        System.out.println("Before increasing length");
//        printStructureInfo();
        int foldingPosition = currentFoldingLength;
        currentFoldingLength++;

        currentDotBracketForm.append('.');
        ringList[foldingPosition].setPosition(foldingPosition);
        ringList[foldingPosition].setType(RINGNODE_CODE.UNPAIRED);

        RingNode up = ringList[foldingPosition - 1].getUp();
        if (up == null) {
            ringList[foldingPosition - 1].setNext(ringList[foldingPosition]);
            ringList[foldingPosition].setPrev(ringList[foldingPosition - 1]);

        } else {
            up.setNext(ringList[foldingPosition]);
            ringList[foldingPosition].setPrev(up);
        }

        //set link of the last node   
        ringList[foldingPosition].setNext(ringList[primarySequence.getSequenceLength()]);
        ringList[primarySequence.getSequenceLength()].setPrev(ringList[foldingPosition]);

//        System.out.println("-------------------------");
//        System.out.println("After increasing length");
//        printStructureInfo();
    }

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    //////////////////  Enuneration possible conformation ///////
    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    public ActionList enumerateActions(boolean willDelayComputing) {
//        System.out.println("[Function - Enumerate conformations by applying elementary actions]");
//        System.out.println(" Enumerate action for new nucleotide sequence: " + primarySequence.getFoldingSequence(currentFoldingLength));        

        ActionList actionList = new ActionList();
        int i;

        //insert base-pair for root
//        System.out.println(" + insert base-pair at utmost level");
        form_bp(root, actionList, willDelayComputing);
        for (i = 0; i < currentFoldingLength; i++) {
            if (pairList[i] > i) {
                /* insert base-pair for neighbours */
//                System.out.println(" + insert base-pair");
                form_bp(ringList[i], actionList, willDelayComputing);

                /* delete base-pair neighbours */
//                System.out.println(" + delete base-pair");
                break_bp(ringList[i], actionList, willDelayComputing);

                /* shift pair neighbour */
//                System.out.println(" + shift base-pair");
                shift_bp(ringList[i], actionList, willDelayComputing);
            }
        }

//        System.out.println("all actions by enumerate-action:");
//        if(actionList.getAllActions().length == 0){
//            System.out.println(" null ");
//        }else{                
//            for(Action MULTILOOP_A : actionList.getAllActions()){
//                System.out.println(" - " + MULTILOOP_A);
//            }
//        }
        return actionList;
    }

    public ActionList enumerateMoreActions(boolean willDelayComputing) {
//        System.out.println("[Function - Enumerate more conformations with extending sequence]");
//        System.out.println(" Enumerate more action for new nucleotide sequence: " + primarySequence.getFoldingSequence(currentFoldingLength));        
        ActionList actionList = new ActionList();

//        System.out.println(" + insert more base-pair at utmost level");        
        form_more_bp(root, actionList, willDelayComputing);
//        int i;        
//        for (i = 0; i < currentFoldingLength; i++) {
//            if (pairList[i] > i) {
////                System.out.println(" + shift more base-pair");
//                shift_more_bp(ringList[i]);
//            }
//        }
//        System.out.println(" + shift more base-pair");
        shift_more_bp(root, actionList, willDelayComputing);

//        System.out.println("all actions by enumerate-more-actions:");
//        if(actionList.getAllActions().length == 0){
//            System.out.println(" null ");
//        }else{                
//            for(Action MULTILOOP_A : actionList.getAllActions()){
//                System.out.println(" - " + MULTILOOP_A);
//            }
//        }
        return actionList;
    }

    public void applyAction(Action action) {
//        System.out.println("apply action " + action);

//        ArrayList<Integer> actionInfo = action.getInfo();
//        int i = actionInfo.get(0);
//        int j = actionInfo.get(1);
        int i = action.getInfoI();
        int j = action.getInfoJ();

        RingNode help_rli;
        RingNode help_rlj;

        switch (action.getActionType()) {
            case ACTION_CODE.ADDING_BASE_PAIR:
                /* insert */
                do_close_bp(ringList[i], ringList[j]);
                break;
            case ACTION_CODE.REMOVING_BASE_PAIR:
                /* delete */
                do_open_bp(ringList[i]);
                break;
            case ACTION_CODE.SHIFTING_HEAD_BASE_PAIR:
                /* shift head*/
                RingNode old_rli = ringList[j].getUp();
                do_open_bp(old_rli);

                if (i < j) {
                    help_rli = ringList[i];
                    help_rlj = ringList[j];
                } else {
                    help_rli = ringList[j];
                    help_rlj = ringList[i];
                }

                do_close_bp(help_rli, help_rlj);
                break;

            case ACTION_CODE.SHIFTING_TAIL_BASE_PAIR:
                /* shift tail*/
                do_open_bp(ringList[i]);

                if (i < j) {
                    help_rli = ringList[i];
                    help_rlj = ringList[j];
                } else {
                    help_rli = ringList[j];
                    help_rlj = ringList[i];
                }

                do_close_bp(help_rli, help_rlj);

                break;
        }
        currentStructureEnergy = action.getEnergy();
    }

    public double evaluateEnergyByAction(Action action) {
//        System.out.println("evaluate energy og structure with action " + action);

//        ArrayList<Integer> actionInfo = action.getInfo();
//        int i = actionInfo.get(0);
//        int j = actionInfo.get(1);
        double energy = 0.0;

        int i = action.getInfoI();
        int j = action.getInfoJ();

        RingNode help_rli;
        RingNode help_rlj;

        switch (action.getActionType()) {
            case ACTION_CODE.ADDING_BASE_PAIR:
                /* insert */
                do_close_bp(ringList[i], ringList[j]);
                
                energy = getEnergyofStructure();
//                energy = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);                
                
                do_open_bp(ringList[i]);
                break;
            case ACTION_CODE.REMOVING_BASE_PAIR:
                /* delete */
                do_open_bp(ringList[i]);

                energy = getEnergyofStructure();                
//                energy = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);

                do_close_bp(ringList[i], ringList[j]);
                break;
            case ACTION_CODE.SHIFTING_HEAD_BASE_PAIR:
                /* shift head*/
                RingNode old_rli = ringList[j].getUp();
                do_open_bp(old_rli);

                if (i < j) {
                    help_rli = ringList[i];
                    help_rlj = ringList[j];
                } else {
                    help_rli = ringList[j];
                    help_rlj = ringList[i];
                }

                do_close_bp(help_rli, help_rlj);

                energy = getEnergyofStructure();                
//                energy = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);

                do_open_bp(help_rli);

                do_close_bp(old_rli, ringList[j]);

                break;

            case ACTION_CODE.SHIFTING_TAIL_BASE_PAIR:
                /* shift tail*/
                RingNode old_rlj = ringList[i].getDown();
                do_open_bp(ringList[i]);

                if (i < j) {
                    help_rli = ringList[i];
                    help_rlj = ringList[j];
                } else {
                    help_rli = ringList[j];
                    help_rlj = ringList[i];
                }

                do_close_bp(help_rli, help_rlj);

                energy = getEnergyofStructure();                
//                energy = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);

                do_open_bp(help_rli);

                do_close_bp(ringList[i], old_rlj);
                break;
        }
        return energy;
    }

    private void form_bp(RingNode root, ActionList actionList, boolean willDelayComputing) {
//        System.out.println("  + insert base-pair starting from " + (root.getPosition() == -1? "root" : "node at position " + root.getPosition()));
        Action action;
        double EoT;

        RingNode stop, rli, rlj;

        /* loop ringlist over all possible i positions */
        stop = root.getDown();
        for (rli = stop.getNext(); rli != stop; rli = rli.getNext()) {
            /* potential i-position is already paired */
            if (rli.getType() == RINGNODE_CODE.PAIRED) {
                continue;
            }

            /* loop ringlist over all possible j positions */
            for (rlj = rli.getNext(); rlj != stop; rlj = rlj.getNext()) {
                /* potential j-position is already paired */
                if (rlj.getType() == RINGNODE_CODE.PAIRED) {
                    continue;
                }

                /* base pair must enclose at least some base defined by paring constraint */
                if (rlj.getPosition() - rli.getPosition() <= ENERGY_CONSTANTS.PAIRING_CONSTRAINT) {
                    continue;
                }

                /* if i-j can form MULTILOOP_A base pair*/
                if (primarySequence.isPairable(rli.getPosition(), rlj.getPosition())) {
                    boolean isEnergyComputed;     
//                    System.out.println(" - examine to create a base-pair (" + rli.getPosition() + " , " + rlj.getPosition() + ")");
                        /* close the base bair*/
                    do_close_bp(rli, rlj);

                    if (!willDelayComputing) {
                        /* evaluate energy of the structure */
                        EoT = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);
                        isEnergyComputed = true;
                    }    
                    else {
                        String structureId = getId();
                        if(isEnergyCached(structureId)){
                            EoT = lookupEnergyCache(structureId);
                            isEnergyComputed = true;
                        }
                        else{
                            /* evaluate energy of the structure */
                            EoT = optimalEnergy;
                            isEnergyComputed = false;
                        }
                    }
                    /* open the base pair again... */
                    do_open_bp(rli);                    

                    /* put the move and the enegy of the structure into the action list */
//                    ArrayList<Integer> actionInfo = new ArrayList();
//                    actionInfo.add(rli.getPosition());
//                    actionInfo.add(rlj.getPosition());
//                    action = new Action(ACTION_CODE.ADDING_BASE_PAIR, EoT, actionInfo);
                    action = new Action(ACTION_CODE.ADDING_BASE_PAIR, EoT, rli.getPosition(), rlj.getPosition(), isEnergyComputed);
                    actionList.addAction(action);

//                    System.out.println("  => create action " + action.toString());
                }
            }
        }
    }

    private void form_more_bp(RingNode root, ActionList actionList, boolean willDelayComputing) {
//        System.out.println("  + insert more base-pair starting from " + (root.getPosition() == -1? "root" : "node at position " + root.getPosition()));
        Action action;
        double EoT;

        RingNode stop, rli, rlj;

        rlj = ringList[currentFoldingLength - 1];

        /* loop ringlist over all possible i positions */
        stop = root.getDown();
        for (rli = stop.getNext(); rli != stop; rli = rli.getNext()) {
            /* potential i-position is already paired */
            if (rli.getType() == RINGNODE_CODE.PAIRED) {
                continue;
            }

            /* base pair must enclose at least some base defined by paring constraint */
            if (rlj.getPosition() - rli.getPosition() <= ENERGY_CONSTANTS.PAIRING_CONSTRAINT) {
                continue;
            }

            /* if i-j can form MULTILOOP_A base pair*/
            if (primarySequence.isPairable(rli.getPosition(), rlj.getPosition())) {
//                System.out.println(" - examine to create a base-pair (" + rli.getPosition() + " , " + rlj.getPosition() + ")");
                boolean isEnergyComputed;                
                /* close the base bair*/
                do_close_bp(rli, rlj);
                if (!willDelayComputing) {
                    /* evaluate energy of the structure */
                    EoT = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);
                    isEnergyComputed = true;
                }

                else {
                    String structureId = getId();
                    if(isEnergyCached(structureId)){
                        EoT = lookupEnergyCache(structureId);
                        isEnergyComputed = true;
                    }
                    else{
                        /* evaluate energy of the structure */
                        EoT = optimalEnergy;
                        isEnergyComputed = false;
                    }
                }

                /* open the base pair again... */
                do_open_bp(rli);
                
                /* put the move and the enegy of the structure into the action list */
//                ArrayList<Integer> actionInfo = new ArrayList();
//                actionInfo.add(rli.getPosition());
//                actionInfo.add(rlj.getPosition());
//                action = new Action(ACTION_CODE.ADDING_BASE_PAIR, EoT, actionInfo);

                action = new Action(ACTION_CODE.ADDING_BASE_PAIR, EoT, rli.getPosition(), rlj.getPosition(), isEnergyComputed);
                actionList.addAction(action);

//                System.out.println("   => create action " + action.toString());
            }
        }
    }

    private void break_bp(RingNode rli, ActionList actionList, boolean willDelayComputing) {
//        System.out.println("  + delete base-pair for node position " + rli.getPosition());
        Action action;

        double EoT;

        RingNode rlj, r;
        rlj = rli.getDown();
        
        boolean isEnergyComputed;        
        
        /* open the base bair*/
        do_open_bp(rli);
//        System.out.println(" - examine to delete base-pair (" + rli.getPosition() + " , " + rlj.getPosition() + ")");
        if (!willDelayComputing) {
            /* evaluate energy of the structure */
            
            EoT = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);
            isEnergyComputed = true;
        }
        else {
            String structureId = getId();
            if(isEnergyCached(structureId)){
                EoT = lookupEnergyCache(structureId);
                isEnergyComputed = true;
            }
            else{
                /* evaluate energy of the structure */
                EoT = optimalEnergy;
                isEnergyComputed = false;
            }
        }

        /* close the base pair again... */
        do_close_bp(rli, rlj);        
        
        
        /* put the move and the enegy of the structure into the action list */
//        ArrayList<Integer> actionInfo = new ArrayList();
//        actionInfo.add(rli.getPosition());
//        actionInfo.add(rlj.getPosition());
//        action = new Action(ACTION_CODE.REMOVING_BASE_PAIR, EoT, actionInfo);
        action = new Action(ACTION_CODE.REMOVING_BASE_PAIR, EoT, rli.getPosition(), rlj.getPosition(), isEnergyComputed);
        
        actionList.addAction(action);

//        System.out.println("   => create action " + action.toString());
    }

    private void shift_bp(RingNode rli, ActionList actionList, boolean willDelayComputing) {
//        System.out.println("  + shift base-pair for node position " + rli.getPosition());
        Action action;

        double EoT;
        int x;
        RingNode rlj, stop, help_rli, help_rlj;

        stop = rli.getDown();

        /* examine interior loop of bp(ij); (.......)
           i of j move                      ->   <- */
        for (rlj = stop.getNext(); rlj != stop; rlj = rlj.getNext()) {
            /* prevent shifting to paired position */
            if (rlj.getType() == RINGNODE_CODE.PAIRED || rlj.getType() == RINGNODE_CODE.ROOT) {
                continue;
            }
            /* j-position of base pair shifts to k position (ij)->(ik) i<k<j */
            if ((rlj.getPosition() - rli.getPosition() > ENERGY_CONSTANTS.PAIRING_CONSTRAINT)
                    && (primarySequence.isPairable(rli.getPosition(), rlj.getPosition()))) {
                boolean isEnergyComputed;                
                
                /* open original basepair */
                do_open_bp(rli);

                /* close shifted version of original basepair */
                do_close_bp(rli, rlj);

//                System.out.println(" - examine to (internal) shift base-pair (" + rli.getPosition() + " , " + stop.getPosition() + ") -> (" + rli.getPosition() + " , " + rlj.getPosition() + ")");
                                        
                if (!willDelayComputing) {
                    /* evaluate energy of the structure */

                    EoT = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);
                    isEnergyComputed = true;
                }
                else {
                    String structureId = getId();
                    if(isEnergyCached(structureId)){
                        EoT = lookupEnergyCache(structureId);
                        isEnergyComputed = true;
                    }
                    else{
                        /* evaluate energy of the structure */
                        EoT = optimalEnergy;
                        isEnergyComputed = false;
                    }
                }
                    
                /* open shifted basepair */
                do_open_bp(rli);

                /* restore original basepair */
                do_close_bp(rli, stop);
                
                /* put the move and the enegy of the structure into the action list */
//                ArrayList<Integer> actionInfo = new ArrayList();
//                actionInfo.add(rli.getPosition());
//                actionInfo.add(rlj.getPosition());
//                action = new Action(ACTION_CODE.SHIFTING_TAIL_BASE_PAIR, EoT, actionInfo);
                action = new Action(ACTION_CODE.SHIFTING_TAIL_BASE_PAIR, EoT, rli.getPosition(), rlj.getPosition(), isEnergyComputed);
                actionList.addAction(action);

//                System.out.println("   => create action " + action.toString());
            }
            /* i-position of base pair shifts to position k (ij)->(kj) i<k<j */
            if ((stop.getPosition() - rlj.getPosition() > ENERGY_CONSTANTS.PAIRING_CONSTRAINT)
                    && (primarySequence.isPairable(stop.getPosition(), rlj.getPosition()))) {
                boolean isEnergyComputed;
                
                /* open original basepair */
                do_open_bp(rli);

                /* close shifted version of original basepair */
                do_close_bp(rlj, stop);

//                System.out.println(" - examine to (internal) shift base-pair (" + rli.getPosition() + " , " + stop.getPosition() + ") -> (" + rlj.getPosition() + " , " + stop.getPosition() + ")");

                if (!willDelayComputing) {
                    /* evaluate energy of the structure */

                    EoT = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);
                    isEnergyComputed = true;
                }
                else {
                    String structureId = getId();
                    if(isEnergyCached(structureId)){
                        EoT = lookupEnergyCache(structureId);
                        isEnergyComputed = true;
                    }
                    else{
                        /* evaluate energy of the structure */
                        EoT = optimalEnergy;
                        isEnergyComputed = false;
                    }
                }

                /* open shifted basepair */
                do_open_bp(rlj);

                /* restore original basepair */
                do_close_bp(rli, stop);
                
//              /* put the move and the enegy of the structure into the neighbour list */
//                ArrayList<Integer> actionInfo = new ArrayList();
//                actionInfo.add(rlj.getPosition());
//                actionInfo.add(stop.getPosition());
//                action = new Action(ACTION_CODE.SHIFTING_HEAD_BASE_PAIR, EoT, actionInfo);
                action = new Action(ACTION_CODE.SHIFTING_HEAD_BASE_PAIR, EoT, rlj.getPosition(), stop.getPosition(), isEnergyComputed);
                actionList.addAction(action);

//                System.out.println("   => create action " + action.toString());
            }
        }
        /* examine exterior loop of bp(ij);   (.......)
           i or j moves                     <-         -> */
        for (rlj = rli.getNext(); rlj != rli; rlj = rlj.getNext()) {
            if ((rlj.getType() == RINGNODE_CODE.PAIRED) || (rlj.getType() == RINGNODE_CODE.UNKNOWN)) {
                continue;
            }
            x = rlj.getPosition() - rli.getPosition();
            if (x < 0) {
                x = -x;
            }
            /* j-position of base pair shifts to position k */
            if ((x > ENERGY_CONSTANTS.PAIRING_CONSTRAINT) && (primarySequence.isPairable(rli.getPosition(), rlj.getPosition()))) {
                if (rli.getPosition() < rlj.getPosition()) {
                    help_rli = rli;
                    help_rlj = rlj;
                } else {
                    help_rli = rlj;
                    help_rlj = rli;
                }
                
                boolean isEnergyComputed;
                
                /* open original basepair */
                do_open_bp(rli);

                /* close shifted version of original basepair */
                do_close_bp(help_rli, help_rlj);

//                System.out.println(" - examine to (external) shift base-pair (" + rli.getPosition() + " , " + stop.getPosition() + ") -> (" + rli.getPosition() + " , " + rlj.getPosition() + ")");

                if (!willDelayComputing) {
                    /* evaluate energy of the structure */

                    EoT = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);
                    isEnergyComputed = true;
                }
                else {
                    String structureId = getId();
                    if(isEnergyCached(structureId)){
                        EoT = lookupEnergyCache(structureId);
                        isEnergyComputed = true;
                    }
                    else{
                        /* evaluate energy of the structure */
                        EoT = optimalEnergy;
                        isEnergyComputed = false;
                    }
                }

                /* open shifted base pair */
                do_open_bp(help_rli);

                /* restore original basepair */
                do_close_bp(rli, stop);
                
                /* put the move and the enegy of the structure into the neighbour list */
//                ArrayList<Integer> actionInfo = new ArrayList();
//                actionInfo.add(rli.getPosition());
//                actionInfo.add(rlj.getPosition());
//                action = new Action(ACTION_CODE.SHIFTING_TAIL_BASE_PAIR, EoT, actionInfo);
                action = new Action(ACTION_CODE.SHIFTING_TAIL_BASE_PAIR, EoT, rli.getPosition(), rlj.getPosition(), isEnergyComputed);
                actionList.addAction(action);

//                System.out.println("   => create action " + action.toString());
            }

            x = rlj.getPosition() - stop.getPosition();
            if (x < 0) {
                x = -x;
            }
            /* i-position of base pair shifts to position k */
            if ((x > ENERGY_CONSTANTS.PAIRING_CONSTRAINT) && (primarySequence.isPairable(stop.getPosition(), rlj.getPosition()))) {
                if (stop.getPosition() < rlj.getPosition()) {
                    help_rli = stop;
                    help_rlj = rlj;
                } else {
                    help_rli = rlj;
                    help_rlj = stop;
                }
                
                boolean isEnergyComputed;                
                /* open original basepair */
                do_open_bp(rli);

                /* close shifted version of original basepair */
                do_close_bp(help_rli, help_rlj);

//                System.out.println(" - examine to (external) shift base-pair (" + rli.getPosition() + " , " + stop.getPosition() + ") -> (" + rlj.getPosition() + " , " + stop.getPosition() + ")");
                if (!willDelayComputing) {
                    /* evaluate energy of the structure */
                    EoT = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);
                    isEnergyComputed = true;
                }
                else {
                    String structureId = getId();
                    if(isEnergyCached(structureId)){
                        EoT = lookupEnergyCache(structureId);
                        isEnergyComputed = true;
                    }
                    else{
                        /* evaluate energy of the structure */
                        EoT = optimalEnergy;
                        isEnergyComputed = false;
                    }
                }                    
                    
                /* open shifted basepair */
                do_open_bp(help_rli);
                /* restore original basepair */
                do_close_bp(rli, stop);
                
                /* put the move and the enegy of the structure into the neighbour list */
//                ArrayList<Integer> actionInfo = new ArrayList();
//                actionInfo.add(rlj.getPosition());
//                actionInfo.add(stop.getPosition());
//                action = new Action(ACTION_CODE.SHIFTING_HEAD_BASE_PAIR, EoT, actionInfo);
                action = new Action(ACTION_CODE.SHIFTING_HEAD_BASE_PAIR, EoT, rlj.getPosition(), stop.getPosition(), isEnergyComputed);
                actionList.addAction(action);

//                System.out.println("   => create action " + action.toString());
            }
        }
    }

    private void shift_more_bp(RingNode root, ActionList actionList, boolean willDelayComputing) {
//        System.out.println("  + shift more base-pair starting from " + (root.getPosition() == -1? "root" : "node at position " + root.getPosition()));
        Action action;

        double EoT;
        RingNode rlj, stop, rli;

        RingNode rlk = ringList[currentFoldingLength - 1];

        stop = root.getDown();

        /* only examine exterior loop of bp(ij);   (.......)
           i or j moves                     <-         -> */
        for (rli = stop.getNext(); rli != stop; rli = rli.getNext()) {
            if (rli.getType() == RINGNODE_CODE.PAIRED) {
                rlj = rli.getDown();

                /* j-position of base pair shifts to position k */
                if ((rlk.getPosition() - rli.getPosition() > ENERGY_CONSTANTS.PAIRING_CONSTRAINT) && (primarySequence.isPairable(rli.getPosition(), rlk.getPosition()))) {
                    boolean isEnergyComputed;
                    
                    /* open original basepair */
                    do_open_bp(rli);

                    /* close shifted version of original basepair */
                    do_close_bp(rli, rlk);

//                    System.out.println(" - examine to (external) shift base-pair (" + rli.getPosition() + " , " + rlj.getPosition() + ") -> (" + rli.getPosition() + " , " + rlk.getPosition() + ")");

                    if (!willDelayComputing) {
                        /* evaluate energy of the structure */
                        EoT = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);
                        isEnergyComputed = true;
                    }
                    else {
                        String structureId = getId();
                        if(isEnergyCached(structureId)){
                            EoT = lookupEnergyCache(structureId);
                            isEnergyComputed = true;
                        }
                        else{
                            /* evaluate energy of the structure */
                            EoT = optimalEnergy;
                            isEnergyComputed = false;
                        }
                    }          
                    
                    /* open shifted base pair */
                    do_open_bp(rli);

                    /* restore original basepair */
                    do_close_bp(rli, rlj);
                    
                    /* put the move and the enegy of the structure into the neighbour list */
//                    ArrayList<Integer> actionInfo = new ArrayList();
//                    actionInfo.add(rli.getPosition());
//                    actionInfo.add(rlk.getPosition());
//                    action = new Action(ACTION_CODE.SHIFTING_TAIL_BASE_PAIR, EoT, actionInfo);
                    action = new Action(ACTION_CODE.SHIFTING_TAIL_BASE_PAIR, EoT, rli.getPosition(), rlk.getPosition(), isEnergyComputed);
                    actionList.addAction(action);

//                    System.out.println("   => create action " + action.toString());
                }

                /* i-position of base pair shifts to position k */
                if ((rlk.getPosition() - rlj.getPosition() > ENERGY_CONSTANTS.PAIRING_CONSTRAINT) && (primarySequence.isPairable(rlj.getPosition(), rlk.getPosition()))) {
                    boolean isEnergyComputed;
                    
                    /* open original basepair */
                    do_open_bp(rli);

                    /* close shifted version of original basepair */
                    do_close_bp(rlj, rlk);

//                    System.out.println(" - examine to (external) shift base-pair (" + rli.getPosition() + " , " + rlj.getPosition() + ") -> (" + rlj.getPosition() + " , " + rlk.getPosition() + ")");

                    if (!willDelayComputing) {
                        /* evaluate energy of the structure */
                        EoT = ComputingMachine.computeStructureEnergy(primarySequence, pairList, currentFoldingLength);
                        isEnergyComputed = true;
                    }
                    else {
                        String structureId = getId();
                        if(isEnergyCached(structureId)){
                            EoT = lookupEnergyCache(structureId);
                            isEnergyComputed = true;
                        }
                        else{
                            /* evaluate energy of the structure */
                            EoT = optimalEnergy;
                            isEnergyComputed = false;
                        }
                    }

                    /* open shifted basepair */
                    do_open_bp(rlj);

                    /* restore original basepair */
                    do_close_bp(rli, rlj);
                    
                    /* put the move and the enegy of the structure into the neighbour list */
//                    ArrayList<Integer> actionInfo = new ArrayList();
//                    actionInfo.add(rlk.getPosition());
//                    actionInfo.add(rlj.getPosition());
//                    action = new Action(ACTION_CODE.SHIFTING_HEAD_BASE_PAIR, EoT, actionInfo);
                    action = new Action(ACTION_CODE.SHIFTING_HEAD_BASE_PAIR, EoT, rlk.getPosition(), rlj.getPosition(), isEnergyComputed);
                    actionList.addAction(action);

//                    System.out.println("   => create action " + action.toString());
                }

            }
        }
    }

    private void do_close_bp(RingNode i, RingNode j) {
        RingNode jNext;

        currentDotBracketForm.setCharAt(i.getPosition(), '(');
        currentDotBracketForm.setCharAt(j.getPosition(), ')');

        //pair list
        pairList[i.getPosition()] = j.getPosition();
        pairList[j.getPosition()] = i.getPosition();

        //change tree representation
        jNext = j.getNext();

        i.setType(RINGNODE_CODE.PAIRED);
        i.setDown(j);
        i.getNext().setPrev(j);

        j.setType(RINGNODE_CODE.PAIRED);
        j.setUp(i);
        j.getNext().setPrev(i);
        j.setNext(i.getNext());

        i.setNext(jNext);
    }

    private void do_open_bp(RingNode i) {
        RingNode node;

        //change to dots
        currentDotBracketForm.setCharAt(i.getPosition(), '.');
        currentDotBracketForm.setCharAt(i.getDown().getPosition(), '.');

        //pair list
        pairList[i.getPosition()] = -1;
        pairList[i.getDown().getPosition()] = -1;

        //change tree representation
        node = i.getNext();

        i.setType(RINGNODE_CODE.UNPAIRED);
        i.getDown().setType(RINGNODE_CODE.UNPAIRED);
        i.setNext(i.getDown().getNext());
        i.getNext().setPrev(i);
        node.setPrev(i.getDown());
        i.getDown().setNext(node);
        i.setDown(null);
        node.getPrev().setUp(null);
    }

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    //////////////////  decomposition of structure into loops //////////////
    /////////////////////////////////////////////////////////////    
    /////////////////////////////////////////////////////////////
    public void secondaryLoopDecomposion() {
        System.out.println("--- decomposion of the current structure ---");
        ArrayList<SecondaryLoop> components = new ArrayList<>();
        for (int i = 0; i < currentFoldingLength; i++) {
            if (pairList[i] > i) {
                composingLoop(i, components);
                i = pairList[i];
            }
        }

        for (SecondaryLoop l : components) {
            System.out.println(l);
        }
    }

    private void composingLoop(int i, ArrayList<SecondaryLoop> components) {
        SecondaryLoop l = new SecondaryLoop();
        boolean isSimpleLoop = true;

        int j, p, q;

        j = pairList[i];

        p = i;
        q = j;
        while (p < q) {
            while (pairList[++p] == -1);
            while (pairList[--q] == -1);

            /* MULTILOOP_A hairpin or multi-loop */
            if ((pairList[q] != p) || (p > q)) {
                isSimpleLoop = false;
                break;
            }

            int n1 = p - i - 1;
            int n2 = j - q - 1;

            int nl, ns;

            if (n1 > n2) {
                nl = n1;
                ns = n2;
            } else {
                nl = n2;
                ns = n1;
            }

            if (nl == 0) {
                /* stack */
//                System.out.println(" stack (" + i + " < " + p + " , " + q + " > "+ j + ")");
                l.addSackIndex(i);
            } else if (ns == 0) {
                /* bulge */
//                System.out.println(" => create bulge (" + i + " < " + p + " , " + q + " > "+ j + ")");
                l.setLoopType(LOOP_TYPE.BULGE);
                l.addLoopIndex(i);
                l.addLoopIndex(p);
                components.add(l);
                l = new SecondaryLoop();

            } else {
                /* interior loop */
//                System.out.println(" => create interior (" + i + " < " + p + " , " + q + " > "+ j + ")");
                l.setLoopType(LOOP_TYPE.INTERNAL);
                l.addLoopIndex(i);
                l.addLoopIndex(p);
                components.add(l);
                l = new SecondaryLoop();
            }

            //loopStructure(i, j, p, q);            
            i = p;
            j = q;
        }

        if (!isSimpleLoop) {
            /* p,q don't pair must have found hairpin or multiloop */
            if (p > q) {
                /* hairpin */
//                System.out.println(" => create hairpin enclosed by (" + i + ", "+ j + ")");
                l.setLoopType(LOOP_TYPE.HAIRPIN);
                l.addLoopIndex(i);
                l.addLoopIndex(p);
                components.add(l);
                return;
            }

            /* multiloop */
//            System.out.println(" => create multi-loop enclosed by (" + i + ", "+ j + ")");
            l.setLoopType(LOOP_TYPE.MULTILOOP);
            l.addLoopIndex(i);
            while (p < j) {
                /* add up the contributions of the substructures of the multiloop */
                l.addLoopIndex(p);
                composingLoop(p, components);
                p = pairList[p];
                while (pairList[++p] == -1);
            }
            components.add(l);
        }
    }    
}
