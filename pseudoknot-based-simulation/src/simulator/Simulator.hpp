// header file for the simulator class modelling the RNA folding

#pragma once
#include "../secondary/RNASecondaryStructure.hpp"

class Simulator {
    public:
        Simulator();

        //running the actual simulation 
        void runSim(double duration, RNASecondaryStructure secondary,std::string file ,bool cotranscription = false, double transSpeed=5); // 40 somewhat arbitrary chosen

        //running the simulation with defined number of steps
        void simSteps(int steps, RNASecondaryStructure secondary,std::string file ,bool cotranscription, double transSpeed=5);
        //getResult() or runSim modified to return the result 

    private:
        //results attribute ?

};