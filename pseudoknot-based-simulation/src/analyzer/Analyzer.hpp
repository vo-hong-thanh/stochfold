
#pragma once
#include <string>
#include "secondary/RNASecondaryStructure.hpp"

class Analyzer {
    public:
    void countAnalysis( double starttme, double interval,double maxtime , int thr,std::vector<std::string> dotBrackets,std::string ifile,std::string ofile);
    void countAnalysis(    double starttme, double interval,double maxtime , int thr,int mostFreq, std::string ifile,std::string ofile);
    void timeAnalysis( double starttme, double interval,double maxtime , int thr,std::vector<std::string> dotBrackets,std::string ifile,std::string ofile);
    std::pair<int,double>  neighboranalysis(RNASecondaryStructure sec);
    void pathAnalysis( int thr, std::string ifile);

};