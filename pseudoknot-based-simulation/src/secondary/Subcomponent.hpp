// class for all subcomponents: Interior, Bulge Exrerior, Multiloop, Hairpin , Pseudoknot, Stem, Unpaired

#pragma once 

#include <vector>
#include <memory>

enum SCType{
    Root,
    Exterior,
    Stem,
    Interior,
    Bulge,
    Hairpin, 
    Multiloop,
    Pseudoknot,  //will be removed 
    PseudoknotH,
    PseudoknotK,
    PseudoknotL, // not working 
    PseudoknotM  // not working 
};



class Subcomponent {
    public:
        Subcomponent(SCType _type, std::vector<int> defBases);

       const std::vector<int> & getDefiningBases();

        void addChild(std::shared_ptr<Subcomponent> child);

        SCType getType();

        int getLeft() const;
        int getRight();
        //get all children
        const  std::vector<std::shared_ptr<Subcomponent>> & getChildren();
        //get all children in range [i,j]
        const  std::vector<std::shared_ptr<Subcomponent>> getChildren(int i , int j);

        friend bool operator==( std::shared_ptr<Subcomponent> lhs, std::shared_ptr<Subcomponent> rhs);

    private:

        SCType type;
        //a vector of those basis indices to perfectly define the substructure with the help of bp[] 
        std::vector<int> defB;
        std::vector<std::shared_ptr<Subcomponent>> children;
    
};


namespace std {

  template <>
  struct hash<Subcomponent>
  {
    std::size_t operator()(const Subcomponent& k) const
    {
      return (k.getLeft());
    }
  };

}

class SCLess {
    public:

    bool operator()(const std::shared_ptr<Subcomponent> f,const std::shared_ptr<Subcomponent> s) const {
        if(f->getType() == Exterior && s->getType() != Exterior) return true;
        else if (s->getType() == Exterior && f->getType() != Exterior) return false;
        else return (f->getLeft() < s->getLeft());
    }
};