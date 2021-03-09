// header file for actions I.e elmentary steps the RNA structure can take (adding, deleting and shifting a base pairing)

#pragma once

enum ActionType {
    Addition,
    Deletion
};

class Action {

    public:
        Action(int i, int j,ActionType _type,double energy);

        double getEnergy();

        ActionType getType();
        friend bool operator<(const Action& lhs,const Action& rhs);
        friend bool operator==(const Action& lhs,const Action& rhs);
        int i;
        int j;
    private: 

        ActionType type;
        double energy;
        
};