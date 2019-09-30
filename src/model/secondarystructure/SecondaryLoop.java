/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package model.secondarystructure;

import java.util.ArrayList;

/**
 *
 * @author vot2
 */
public class SecondaryLoop {
    private String loopType;
    private ArrayList<Integer> stemIndexList;
    private ArrayList<Integer> loopIndexList;

    public SecondaryLoop() {
        stemIndexList = new ArrayList<>();
        loopIndexList = new ArrayList<>();
    }

    public void setLoopType(String _loopType) {
        loopType = _loopType;
    }

    public String getLooptype() {
        return loopType;
    }

    public void addSackIndex(int i) {
        stemIndexList.add(i);
    }

    public void addLoopIndex(int i) {
        loopIndexList.add(i);
    }

    @Override
    public String toString() {
        StringBuilder loopInfo = new StringBuilder(loopType);
        loopInfo.append(" (StemInfo: ");
        for (int i : stemIndexList) {
            loopInfo.append(i + " ");
        }
        loopInfo.append(" )");

        loopInfo.append(" (LoopInfo : ");
        for (int i : loopIndexList) {
            loopInfo.append(i + " ");
        }
        loopInfo.append(" )");
        return loopInfo.toString();
    }
}
