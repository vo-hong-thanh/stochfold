
#include "Action.hpp"


Action::Action(int _i, int _j,ActionType _type, double _energy) : i(_i), j(_j), type(_type), energy(_energy) { }

ActionType Action::getType(){
    return type;
}

double Action::getEnergy(){
    return energy;
}
bool operator<(const Action& lhs,const Action& rhs){
    if( lhs.i < lhs.j){
        if(rhs.i < rhs.j){
            if(lhs.i < rhs.i) return true;
            else if (lhs.i == rhs.i && lhs.j < rhs.j) return true;
            else return false;
            
        }else{
            if(lhs.i < rhs.j) return true;
            else if (lhs.i == rhs.j && lhs.j < rhs.i) return true;
            else return false;
        }
    }else{
        if(rhs.i < rhs.j){
            if(lhs.j < rhs.i) return true;
            else if (lhs.j == rhs.i && lhs.i < rhs.j) return true;
            else return false;
            
        }else{
            if(lhs.j < rhs.j) return true;
            else if (lhs.j == rhs.j && lhs.i < rhs.i) return true;
            else return false;
        }
    }
    
        
}

bool operator==(const Action& lhs,const Action& rhs){
    return (lhs.i==rhs.i && lhs.j == rhs.j) || (lhs.i==rhs.j && lhs.j == rhs.i);
}