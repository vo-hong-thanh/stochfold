#pragma once 
#include <string>
// hearedfile for the primary sequence of the RNA strand




class RNAPrimarySequence {

    public:
        //constructor, takes string with 5' end as first character
        RNAPrimarySequence(std::string pSeq,int flen);
        RNAPrimarySequence(std::string pSeq);
        //returns true if i is pairable with j according to Watson-Crick or Wobble
        bool isPairable(int i, int j);

        const std::string & getString();
        bool increaseFoldingLength();

        int getLength();
        int getFoldingLength();
    private:
        std::string nucleotides;
        int sequenceLength;
        int foldingLength;

};