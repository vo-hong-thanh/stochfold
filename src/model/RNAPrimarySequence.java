/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package model;

import constants.NUCLEOTIDE_CODE;

/**
 *
 * @author vot2
 */
public class RNAPrimarySequence {

    private String nucleotides;
    private int sequenceLength;

    public RNAPrimarySequence(String _nucleotides) {
        nucleotides = _nucleotides;
        sequenceLength = nucleotides.length();
    }

    public int getSequenceLength() {
        return sequenceLength;
    }

    public String getFoldingSequence(int foldinglength) {
        return nucleotides.substring(0, foldinglength);
    }

    @Override
    public String toString() {
        return nucleotides;
    }

    public char getNucleotide(int i) {
        return nucleotides.charAt(i);
    }

    public String getSubsequence(int i, int j) {
        return nucleotides.substring(i, j + 1);
    }

    public boolean isPairable(int i, int j) {
        if (i >= sequenceLength || j >= sequenceLength) {
            return false;
        }

        boolean pairable = false;
        //A-U
        if ((nucleotides.charAt(i) == NUCLEOTIDE_CODE.A && nucleotides.charAt(j) == NUCLEOTIDE_CODE.U)
                || (nucleotides.charAt(i) == NUCLEOTIDE_CODE.U && nucleotides.charAt(j) == NUCLEOTIDE_CODE.A)) {
            pairable = true;
        }

        //G-C
        if ((nucleotides.charAt(i) == NUCLEOTIDE_CODE.G && nucleotides.charAt(j) == NUCLEOTIDE_CODE.C)
                || (nucleotides.charAt(i) == NUCLEOTIDE_CODE.C && nucleotides.charAt(j) == NUCLEOTIDE_CODE.G)) {
            pairable = true;
        }

        //G-U: wooble
        if ((nucleotides.charAt(i) == NUCLEOTIDE_CODE.G && nucleotides.charAt(j) == NUCLEOTIDE_CODE.U)
                || (nucleotides.charAt(i) == NUCLEOTIDE_CODE.U && nucleotides.charAt(j) == NUCLEOTIDE_CODE.G)) {
            pairable = true;
        }

        return pairable;
    }
}
