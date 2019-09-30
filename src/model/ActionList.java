/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package model;

import java.util.ArrayList;

/**
 *
 * @author vot2
 */
public class ActionList {
    private ArrayList<Action> actions;
    private double totalRate;
    public ActionList(){
        actions = new ArrayList();
        totalRate = 0;
    }
    
    public void addAction(Action action){
        actions.add(action);
        totalRate += action.getRate();
    }
    
    public void addAllActions(ActionList otherActionList){
        for(Action action : otherActionList.getAllActions()){
            actions.add(action);
            totalRate += action.getRate();
        }        
    }
    
    public Action getAction(int index){
        return actions.get(index);
    }
    
    public Action[] getAllActions(){
        int size = actions.size();
        Action[] list = new Action[size];
        for(int index = 0; index < size; index++){
            list[index] = actions.get(index);
        }
        return list;        
    }
    
    public void updateTotalRate(double amount){
        totalRate += amount;
    }
    
    public double getTotalRate(){
        return totalRate;
    }
    
    public int getNumAction(){
        return actions.size();
    }
}
