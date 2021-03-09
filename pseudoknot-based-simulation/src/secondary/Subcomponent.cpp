
#include "Subcomponent.hpp"
#include <iostream>






Subcomponent::Subcomponent(SCType _type, std::vector<int> defBases) : type(_type), defB(defBases) { }

const std::vector<int> & Subcomponent::getDefiningBases() {
    return defB;
}

void Subcomponent::addChild(std::shared_ptr<Subcomponent> child) {
    children.push_back(child);
    // std::cout<<"adding subcomponent of type: "<<child->getType()<<std::endl; //debugging
}

SCType Subcomponent::getType(){
    return type;
}

const  std::vector<std::shared_ptr<Subcomponent>> & Subcomponent::getChildren(){
    return children;
}
const  std::vector<std::shared_ptr<Subcomponent>>  Subcomponent::getChildren(int i , int j){
    std::vector<std::shared_ptr<Subcomponent>> ret;
    for( auto c : children){
        if(c->getLeft()>=i && c->getRight() <=j) ret.push_back(c);
    }
    return ret;
}
    int Subcomponent::getLeft() const {
        if(type == Multiloop ){
            return defB.empty() ?children[0]->getDefiningBases()[0] : std::min(defB[0],children[0]->getDefiningBases()[0]);
        }else return defB[0];

    }
    int Subcomponent::getRight(){
        if(type == Multiloop){
            return std::max(defB.back(),children.back()->getDefiningBases().back());
        }else return defB.back();
    }

    bool operator==( std::shared_ptr<Subcomponent> lhs, std::shared_ptr<Subcomponent> rhs){
        if(lhs->getType() != rhs->getType()) return false;
        auto bases1 = lhs->getDefiningBases();
        auto bases2 = rhs->getDefiningBases();
        if (bases1.size() != bases2.size()) return false;

        for ( int i=0; i<bases1.size() ; i++){
            if(bases1[i]!= bases2[i]) return false;
        }
        return true;
    }