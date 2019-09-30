/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package model;

import constants.ACTION_CODE;
import utility.ComputingMachine;

/**
 *
 * @author vot2
 */
public class Action {
    private final int actionType;

    private double energy;
    private double rate;

    private final int i;
    private final int j;
    
    private boolean isEnergycomputed;

//    private ArrayList<Integer> info = null;
    public Action(int _actionType, double _energy, int _i, int _j, boolean _isEnergyComputed) {
        actionType = _actionType;
        energy = _energy;

        rate = ComputingMachine.computeRate(energy);

        i = _i;
        j = _j;
        
        isEnergycomputed = _isEnergyComputed;
    }

//    public Action(int _actionType, double _energy, ArrayList<Integer> _sources) {
//        actionType = _actionType;
//        energy = _energy;
//        info = _source;
//    }
    public int getActionType() {
        return actionType;
    }

    public void setEnergy(double _energy) {
        energy = _energy;
        rate = ComputingMachine.computeRate(energy);
        isEnergycomputed = true;
    }

    public boolean checkIsEnergyComputed(){
        return isEnergycomputed;
    }
    
    public double getEnergy() {
        return energy;
    }

    public double getRate() {
        return rate;
//        return ComputingMachine.computeRate(energy);
    }

    public int getInfoI() {
        return i;
    }

    public int getInfoJ() {
        return j;
    }

//    public ArrayList<Integer> getInfo() {
//        return info;
//    }
    public String toReducedString() {
        StringBuilder actionString = new StringBuilder();
        actionString.append("(");

        switch (actionType) {
            case ACTION_CODE.ADDING_BASE_PAIR:
                actionString.append("Form base-pair: " + i + " - " + j);
                break;
            case ACTION_CODE.REMOVING_BASE_PAIR:
                actionString.append("Delete base-pair: " + i + " x " + j);
                break;
            case ACTION_CODE.SHIFTING_HEAD_BASE_PAIR:
                actionString.append("Shift head base-pair: " + i + " - " + j);
                break;
            case ACTION_CODE.SHIFTING_TAIL_BASE_PAIR:
                actionString.append("Shift tail-pair: " + i + " - " + j);
                break;
            case ACTION_CODE.MUTATING_NUCLEOTIDE:
                actionString.append("Not yet implement the mutate action");
                break;
            default:
                actionString.append("Unknown action");
                break;
        }
        actionString.append(")");
        return actionString.toString();
    }

    public String toString() {
        StringBuilder actionString = new StringBuilder();
        actionString.append("(");

        switch (actionType) {
            case ACTION_CODE.ADDING_BASE_PAIR:
                actionString.append("Form base-pair: " + i + " - " + j);
                break;
            case ACTION_CODE.REMOVING_BASE_PAIR:
                actionString.append("Delete base-pair: " + i + " x " + j);
                break;
            case ACTION_CODE.SHIFTING_HEAD_BASE_PAIR:
                actionString.append("Shift head base-pair: " + i + " - " + j);
                break;
            case ACTION_CODE.SHIFTING_TAIL_BASE_PAIR:
                actionString.append("Shift tail-pair: " + i + " - " + j);
                break;
            case ACTION_CODE.MUTATING_NUCLEOTIDE:
                actionString.append("Not yet implement the mutate action");
                break;
            default:
                actionString.append("Unknown action");
                break;
        }

        actionString.append(") resulting in new structure with energy " + energy + " (kcal), transition rate " + getRate());
        return actionString.toString();
    }
}
