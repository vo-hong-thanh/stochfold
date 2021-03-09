

#pragma once
#include "../constants/Constants.hpp"
#include <string>

namespace EnergyCalculator {

    double stemPairEnergy( const char i,const char j, const char p, const char q);

    double multiloopEnergy( int unpaired, int pairs);

    double endPenalty(const char i, const char j);

     double terminalMismatchEnergy(const char i, const char j, const char x, const char y);

     

    double hairpinSpecialEnergy(std::string str);

    double hairpinTetraEnergy(std::string str);

    double internalLoop11Energy(const char i, const char j, const char p, const char q, const char x, const char y);

    double danglingEnergy(const char i, const char j, const char x, bool prime3);

    double internalmmEnergy(const char i, const char j, const char p, const char q);


};