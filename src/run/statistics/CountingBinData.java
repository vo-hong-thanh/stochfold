/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package run.statistics;

import java.util.ArrayList;

/**
 *
 * @author vot2
 */
public class CountingBinData{
    private int count;    
    private ArrayList<Double> dataItem;
    
    public CountingBinData(){
        count = 0;
        dataItem = new ArrayList<Double>();
    }  
    
    public void addItem(double runTime){
        count++;
        dataItem.add(runTime);
    }
    
    public int getCount(){
        return count;
    }    
}

