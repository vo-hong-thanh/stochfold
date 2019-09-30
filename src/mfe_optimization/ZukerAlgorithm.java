/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mfe_optimization;

import model.RNAPrimarySequence;
import constants.ENERGY_CONSTANTS;
import utility.ComputingMachine;

/**
 *
 * @author vot2
 */
public class ZukerAlgorithm {
    RNAPrimarySequence primarySequence;
    StringBuilder dotBracket;

    //minimum folding energy of all non-empty folding sequence i -> j
    double[][] W;
    int[][] traceW;

    //minimum folding energy of all non-empty folding sequence i -> j containing the base-pair (i, j)
    double[][] V;
    int[][] traceV;
    int[][] traceVL;
    int[][] traceVR;

    //minimum folding energy of all non-empty folding sequence i -> j that is part of MULTILOOP_A multiloop
    double[][] pMu;
    int[][] tracepMu;
    
    public void config(String _primarySequence) {
        primarySequence = new RNAPrimarySequence(_primarySequence);

        dotBracket = new StringBuilder();
        for (int i = 0; i < _primarySequence.length(); i++) {
            dotBracket.append(".");
        }

//        System.out.println("Sequence");
//        for (int x = 0; x < primarySequence.getSequenceLength(); x++) {
//            System.out.println("position " + x + " = " + primarySequence.getNucleotide(x));
//        }
    }
    
    public void config(RNAPrimarySequence _primarySequence) {
        primarySequence = _primarySequence;

        dotBracket = new StringBuilder();
        for (int i = 0; i < _primarySequence.getSequenceLength(); i++) {
            dotBracket.append(".");
        }

//        System.out.println("Sequence");
//        for (int x = 0; x < primarySequence.getSequenceLength(); x++) {
//            System.out.println("position " + x + " = " + primarySequence.getNucleotide(x));
//        }
    }

    public void printWorkingMatrix() {
        int foldingLength = primarySequence.getSequenceLength();

        System.out.println("W matrix");
        for (int i = 0; i < foldingLength; i++) {
            for (int j = 0; j < foldingLength; j++) {
                System.out.print("W[" + i + "][" + j + "] = " + W[i][j] + " \t");
            }
            System.out.println();
        }

        System.out.println("V matrix");
        for (int i = 0; i < foldingLength; i++) {
            for (int j = 0; j < foldingLength; j++) {
                System.out.print("V[" + i + "][" + j + "] = " + V[i][j] + " \t");
            }
            System.out.println();
        }
    }

    public void printTraceMatrix() {
        int foldingLength = primarySequence.getSequenceLength();

        System.out.println("traceW");
        for (int i = 0; i < foldingLength; i++) {
            for (int j = 0; j < foldingLength; j++) {
                System.out.print("traceW[" + i + "][" + j + "] = " + traceW[i][j] + " \t");
            }
            System.out.println();
        }

        System.out.println("traceV");
        for (int i = 0; i < foldingLength; i++) {
            for (int j = 0; j < foldingLength; j++) {
                System.out.print("traceV[" + i + "][" + j + "] = " + traceV[i][j] + " \t");
            }
            System.out.println();
        }
        
        System.out.println("tracepMu");
        for (int i = 0; i < foldingLength; i++) {
            for (int j = 0; j < foldingLength; j++) {
                System.out.print("tracepMu[" + i + "][" + j + "] = " + tracepMu[i][j] + " \t");
            }
            System.out.println();
        }
    }
    
    /////////////////////////////////////////////////////////////////
    ///////  Zucker algorithm for computing optimal energy  /////////
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////   
    public void evaluateOptimalEnergy() {        
        int foldingLength = primarySequence.getSequenceLength();
        int pairingConstraint = ENERGY_CONSTANTS.PAIRING_CONSTRAINT;

        //initialization of W and V
        W = new double[foldingLength][foldingLength];
        V = new double[foldingLength][foldingLength];
        pMu = new double[foldingLength][foldingLength];
        
        traceW = new int[foldingLength][foldingLength];
        traceV = new int[foldingLength][foldingLength];
        traceVL = new int[foldingLength][foldingLength];
        traceVR = new int[foldingLength][foldingLength];
        tracepMu = new int[foldingLength][foldingLength];

        for (int i = 0; i < foldingLength; i++) {
            for (int j = 0; j < foldingLength; j++) {
                W[i][j] = 0;
                V[i][j] = Double.MAX_VALUE;
                pMu[i][j] = Double.MAX_VALUE;

                traceW[i][j] = 0;
                traceV[i][j] = 0;
                traceVL[i][j] = 0;
                traceVR[i][j] = 0;
                tracepMu[i][j] = 0;
            }
        }

        //W[i][j] = INF with j - pairingConstraint < i < j 
        for (int i = 0; i < foldingLength; i++) {
            for (int k = 0; k <= pairingConstraint; k++) {
                int j = i + k;
                if(j < foldingLength){
                    W[i][j] = Double.MAX_VALUE;
                    V[i][j] = Double.MAX_VALUE;                     
                }
            }
        }

//        printWorkingMatrix();

        //dynamic update W and V
        double min = Double.MAX_VALUE;
        int i, j;

        for (int d = pairingConstraint + 1; d < foldingLength; d++) {
//            System.out.println("---- d = " + d);
            for (j = d; j < foldingLength; j++) {
                i = j - d;
//                System.out.println("Consider: i = " + i + ", j = " + j);

                min = Double.MAX_VALUE;
                if (min > W[i + 1][j]) {
                    min = W[i + 1][j];
                    
//                    System.out.println(" unpaired => move i up. Set new min energy " + min);
                    traceW[i][j] = -1;
                }

                if (min > W[i][j - 1]) {
                    min = W[i][j - 1];
                    
//                    System.out.println(" unpaired => move j down. Set new min energy " + min);
                    traceW[i][j] = -2;
                }

                double val;
                if (primarySequence.isPairable(i, j)) {
//                    System.out.println(" pairable, calculate energy of new possible structures ");
                    val = computeV(i, j);
                    
//                    System.out.print(" => the structure with smallest energy = " + val);
//                    System.out.print(" and tracebackV[" + i + "][" + j + "] = " + traceV[i][j]);
                    if (min > val) {
                        min = val;
//                        System.out.println(". Set new min energy " + min);
                        
                        traceW[i][j] = -3;
                    }
//                    else{
//                        System.out.println();
//                    }
                }

                for (int k = i + pairingConstraint + 1; k < j - pairingConstraint; k++) {
//                    System.out.print(" separate i = " + i + ", j = " + j + " by k = " + k);
                    val = W[i][k] + W[k + 1][j];
//                    System.out.print(" => energy = " + val);
                    if (min > val) {
                        min = val;
//                        System.out.println(". Set new min energy " + min + ")");
                        traceW[i][j] = k;
                    }
//                    else{
//                        System.out.println();
//                    }
                }
                W[i][j] = min;

//                System.out.println("W[" + i + "][" + j + "] = " + W[i][j]);                
//                System.out.println("V[" + i + "][" + j + "] = " + V[i][j]);
//                
//                System.out.println("tracebackW[" + i + "][" + j + "] = " + traceW[i][j]);
            }
//            System.out.println("-------------------------------------------------");
        }

//        System.out.println("Minimum energy is: " + W[0][foldingLength - 1]);
    }

    public double getOptimalEnergy()
    {
        return W[0][primarySequence.getSequenceLength() - 1];
    }
    
    public double getOptimalEnergyBetweenPosition(int startingPosition, int endingPosition)
    {
        return W[startingPosition][endingPosition];
    }    
    
    private double computeV(int i, int j) {
        double min = Double.MAX_VALUE;
        double val;

        val = hairpin(i, j);
        if (min > val) {
            min = val;            
            traceV[i][j] = -1;
//            System.out.println(" reset the best structure so far to hairpin with energy " + min);
        }

        val = stack(i, j);
        if (min > val) {
            min = val;            
            traceV[i][j] = -2;
            
//            System.out.println(" reset the best structure so far to stack structure with energy " + min);
        }

        val = loop(i, j);
        if (min > val) {
            min = val;
            traceV[i][j] = -3;
            
//            System.out.println(" reset the best structure so far to bulge/interior loop with energy " + min);
        }

        val = multiloop(i, j);
        if (min > val) {
            min = val;
            traceV[i][j] = -4;
            
//            System.out.println(" reset the best structure so far to multiloop with energy " + min);
        }

        V[i][j] = min;               
        return min;
    }

    private double hairpin(int i, int j) {
//        System.out.println("  - call hairpin for (" + i + " , " + j + ")");
        double eHairpin = ComputingMachine.hairpinEnergy(primarySequence, i, j);
//        System.out.println("  => hairpin energy: " + eHairpin);
        return eHairpin;
    }
    
    private double stack(int i, int j) {
//        System.out.println("  - call stack for (" + i + " , " + j + ")");
        int p = i + 1;
        int q = j - 1;

        double eStack = Double.MAX_VALUE;
        if (primarySequence.isPairable(p, q)) {
            eStack = ComputingMachine.stackEnergy(primarySequence, i, j, p, q);
//            System.out.print("  => energy of stack = " + eStack + " and V[" + p + "][" + q +"] = " + V[p][q]);
            eStack += V[p][q];
//            System.out.println(" thus total energy of stack structure " + eStack);
        }
        return eStack;
    }

    private double loop(int i, int j) {
//        System.out.println("  - call bulge/interior loop for (" + i + " , " + j + ")");
        int pairingConstraint = ENERGY_CONSTANTS.PAIRING_CONSTRAINT;
        double min = Double.MAX_VALUE;

        for (int p = i + 1; p < j - pairingConstraint; p++) {
            for (int q = p + pairingConstraint + 1; q < j; q++) {
                if ( (p == i + 1) && (q == j - 1) ){
                    //a stack - already computed
                    continue;
                }

                double eLoop = Double.MAX_VALUE;
                if (primarySequence.isPairable(p, q)) {
//                    System.out.println("    + consider bulge/interior loop for (" + i + " , " + p + " , " + q + " , " + j + ")");
                    eLoop = ComputingMachine.loopEnergy(primarySequence, i, j, p, q);
//                    System.out.print("    => energy of bulge/interior loop = " + eLoop + " and V[" + p + "][" + q +"] = " + V[p][q]);
                    eLoop += V[p][q];                                        
//                    System.out.println(" thus total energy of loop structure " + eLoop);
                }

                if (min > eLoop) {
                    min = eLoop;                    
                    traceVL[i][j] = p;
                    traceVR[i][j] = q;
                }
            }
        }

        return min;
    }

    private double multiloop(int i, int j) {
//        System.out.println("  - call multiloop for (" + i + " , " + j + ")");
        int pairingConstraint = ENERGY_CONSTANTS.PAIRING_CONSTRAINT;
        double min = Double.MAX_VALUE;

        for (int k = i + pairingConstraint + 1; k < j - pairingConstraint - 1; k++) {
            if (primarySequence.isPairable(i + 1, k) && (primarySequence.isPairable(k + 1, j - 1) && ((j - k - 2) > pairingConstraint))  ){
//                System.out.println("    + consider multiloop with parts (" + (i + 1) + " , " + k + ") and (" + (k + 1) + ", " + (j - 1) + ")");

                double eMultiloop = Double.MAX_VALUE;
                double p1 = partMultiloop(i + 1, k);
                double p2 = partMultiloop(k + 1, j - 1);
                if (p1 < Double.MAX_VALUE && p2 < Double.MAX_VALUE) {
                    eMultiloop = p1 + p2 + ENERGY_CONSTANTS.MULTILOOP_A;
                }
//                System.out.println("    => multiloop energy: " + eMultiloop);
                if (min > eMultiloop) {
                    min = eMultiloop;
                    tracepMu[i][j] = k;
                }
            }
        }

        return min;
    }

    public double partMultiloop(int i, int j) {
//        System.out.println("    - call part of multiloop enclosing by (" + i + " , " + j + ")");
        double min = Double.MAX_VALUE;
        
        double pluI = pMu[i + 1][j] + ENERGY_CONSTANTS.MULTILOOP_B;
//        System.out.println("      move i up with energy " + pluI);
        if (min > pluI) {            
            min = pluI;
            tracepMu[i][j] = -1;
        }

        double minusr = pMu[i][j - 1] + ENERGY_CONSTANTS.MULTILOOP_B;
//        System.out.println("      move j down with energy " + minusr);
        if (min > minusr) {
            min = minusr;   
            tracepMu[i][j] = -2;
        }

        double interior = V[i][j] + ENERGY_CONSTANTS.MULTILOOP_C;
//        System.out.println("      V[" + i + "][" + j +"] has energy = " + V[i][j] + " thus total interior structure energy = " + interior);
        if (min > interior) {
            min = interior;
            tracepMu[i][j] = -3;
        }
        
        for (int k = i + 1; k < j; k++) {            
            double s = pMu[i][k] + pMu[k + 1][j];
//            System.out.println("      divide part of multiloop at k = " + k + " with energy " + s);
            if (min > s) {
                min = s;   
                tracepMu[i][j] = k;
            }
        }
        
//        System.out.println("     => min part of multiloop energy " + min);
        pMu[i][j] = min;
        return min;
    }
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    ///////  Trace back algorithm for computing optimal structure  //
    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////
    public String evaluateOptimalStructure() {
        System.out.println("--- calculate optimal structure ---");        
//        System.out.println("trace matrices:");
//        printTraceMatrix();
        
        int foldingLength = primarySequence.getSequenceLength();

        backtrackW(0, foldingLength - 1);
        //System.out.println("Optimal structure is: " + dotBracket);
        return dotBracket.toString();
    }

    private void backtrackW(int i, int j) {
        if (traceW[i][j] == 0) {
            return;
        } else if (traceW[i][j] == -1) {//i up
            backtrackW(i + 1, j);
        } else if (traceW[i][j] == -2) {//j down
            backtrackW(i, j - 1);
        } else if (traceW[i][j] == -3) {//pair i,j
            backtrackV(i, j);
        } else { //separate structure
            int k = traceW[i][j];
            
            backtrackW(i, k);
            backtrackW(k + 1, j);
        }

    }

    private void backtrackV(int i, int j) {
        if (traceV[i][j] == 0) {
            return;
        } else if (traceV[i][j] == -1) {//hairpin
            dotBracket.setCharAt(i, '(');
            dotBracket.setCharAt(j, ')');
        } else if (traceV[i][j] == -2) {//stack
            dotBracket.setCharAt(i, '(');
            dotBracket.setCharAt(j, ')');
            
            backtrackV(i + 1, j - 1);            
        } else if (traceV[i][j] == -3) {//bulge/interior loop
            dotBracket.setCharAt(i, '(');
            dotBracket.setCharAt(j, ')');
            
            int p = traceVL[i][j];
            int q = traceVR[i][j];
            
            backtrackV(p, q);            
        } else {//multiloop
            dotBracket.setCharAt(i, '(');
            dotBracket.setCharAt(j, ')');
            
            int k = tracepMu[i][j];
            
            backtrackpMu(i + 1, k);        
            backtrackpMu(k + 1, j - 1);
        }
    }
    
    private void backtrackpMu(int i, int j) {
        if (tracepMu[i][j] == 0) {
            return;
        } else if (tracepMu[i][j] == -1) {//move i up
            backtrackpMu(i + 1, j);
        } else if (tracepMu[i][j] == -2) {//move j down
            backtrackpMu(i, j - 1);            
        } else if (tracepMu[i][j] == -3) {//interior
            backtrackpMu(i, j);            
        } else {//multiloop
            int k = tracepMu[i][j];
            
            backtrackpMu(i, k);        
            backtrackpMu(k + 1, j);
        }
    }
    /////////////////////////////////////////////////////////////////
}
