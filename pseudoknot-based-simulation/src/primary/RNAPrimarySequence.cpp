// cpp file for the primary sequence of the RNA strand 

#include "RNAPrimarySequence.hpp"

RNAPrimarySequence::RNAPrimarySequence(std::string pSeq) : nucleotides(pSeq) ,sequenceLength(pSeq.size()) , foldingLength(pSeq.size()) { }
RNAPrimarySequence::RNAPrimarySequence(std::string pSeq,int flen) : nucleotides(pSeq) ,sequenceLength(pSeq.size()) , foldingLength(flen) { }


bool RNAPrimarySequence::isPairable(int i,int j) {
    if(i>foldingLength || j > foldingLength) return false;
    if(std::abs(j-i) <4) return false;
    if (nucleotides[i] =='A'){
        if( nucleotides[j]=='U') return true;
        else return false;
    } else if (nucleotides[i] =='G'){
        if( nucleotides[j]=='U' || nucleotides[j]=='C') return true;
        else return false;
    } else if (nucleotides[i] =='U'){
        if( nucleotides[j]=='A' || nucleotides[j]=='G') return true;
        else return false;
    } else if (nucleotides[i] =='C'){
        if( nucleotides[j]=='G') return true;
        else return false;
    }
    return false;
}

const std::string & RNAPrimarySequence::getString() {
    return nucleotides;
}

int RNAPrimarySequence::getLength() {
    return sequenceLength;
}

int RNAPrimarySequence::getFoldingLength() {
    return foldingLength;
}
bool RNAPrimarySequence::increaseFoldingLength() {
    if(foldingLength >= sequenceLength) return false;
    else{
        foldingLength++;
        return true;
    }
}