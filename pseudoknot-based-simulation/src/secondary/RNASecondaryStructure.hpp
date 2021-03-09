#pragma once 
#include <string>
#include <vector>
#include "Action.hpp"
#include "../primary/RNAPrimarySequence.hpp"
#include "Subcomponent.hpp"
#include <functional>
#include <unordered_map>
#include <list>
#include <tuple>

//header file for secondary structure of the RNA strand

class RNASecondaryStructure {

    public:
        //constructor takes dot-bracket form (and pseudoknt numbers?) as input, initializes the internal data structure
        RNASecondaryStructure(RNAPrimarySequence prima, std::string dotBracket, std::vector<int> psnum);

        // creates a vector of all possible actions that can be done for the structure 
        //done by travering tree of subcomponents
        std::vector<Action> enumerateActions();
        // enumerating the actions in the correct order
        std::vector<Action> enumerateActionsOrder();
        //only used in cotranscriptional folding
        // enumerates all actions with the newly added base
        std::vector<Action> additonalActions();

        // calculates the current energy of the sructure and updates the energy value
        double evaluateEnergy();
        // calculates the energy of substructure subroot included
        double evaluateSubstructureEnergy(std::shared_ptr<Subcomponent> subRoot );
        // returns the cached energy should always be same as evaluateEnergy 
        double getEnergy();
        
        //decomposes the structure to loops; hairpin, interior, exterior, multi, bulge and pseudoknot
        void docomposeStructure(std::shared_ptr<Subcomponent> root =nullptr,int start =0 ,int end =-1, std::function<bool(int)> region = [](int in){ return true;});
        // Takes action as input and changes the secondary structure and its energy depending on the  action
        void executeAction(Action & ac);
        //gives last evaluated dot bracket /unnecessary
        const std::string & getDB();
        //evaluates new dotbracket
        const std::string & getDotBracket();

        bool increaseFoldingLength();

        int getFoldingLength();

        std::string getSCTree();

        std::string getSequence();

        void inititEffAction();


    private:
        RNAPrimarySequence primary;
        std::string dotBracket;
        //decomposed form of secondary structure
        std::shared_ptr<Subcomponent> root;
        std::shared_ptr<Subcomponent> oldRoot;
        //base pair function -1 if unpaired and index of pair otherwise
        std::vector<int> bp;

        double energy;
       ///////////////////////////////////////////
        std::vector<std::vector<std::vector<Action>> > currentActions;
        std::list<std::list<double>> currentSCFlux;
 
    

        //dim2 list
        // one element in dim2 list corresponds to a set of actions between two sc (or within a sc)
        // minimum length of 1 elemnt ,the first element is always a within sc 
        // linked list where each element has a vector of actions and sum of their rates
        //std::list<std::tuple<std::vector<Action>,double>> dim2;

        //dim1 list
        //one element in dim1 list corresponds to all actions associated to a sc
        //an element consists of a dim2 list and their sum of rates
        std::list<std::tuple<std::list<std::tuple<std::vector<Action>,double>>,double>> dim1;

        //std::vector<std::tuple<std::shared_ptr<Subcomponent>,std::list<std::tuple<std::vector<Action>> >::iterator >> 
        
        //sc -> { [ (sc,dim2ptr)... ], [sc...], dim1ptr}
        std::unordered_map<Subcomponent,std::tuple<
            std::vector<std::tuple<std::shared_ptr<Subcomponent>,std::list<std::tuple<std::vector<Action>,double> >::iterator >>
            ,std::vector<std::shared_ptr<Subcomponent>>
            ,std::list<std::tuple<std::list<std::tuple<std::vector<Action>,double>>,double>>::iterator 
            >
        ,std::hash<Subcomponent> > map3;

        //////////////////////////////////////////
        double componentEnergy(std::shared_ptr<Subcomponent> comp, int ii = -1, int jj= -1);

        double pkEnergy(std::shared_ptr<Subcomponent> comp);

        double pkEnergyDP09(std::shared_ptr<Subcomponent> comp);

        double terminalEnergy(const int i ,const int j);
        //enumerates all actions inside comp and all
        const std::vector<Action>  subcomponentActions(std::shared_ptr<Subcomponent> comp,const std::vector<std::shared_ptr<Subcomponent>> siblings);
        const std::vector<Action> SCInnerActions(std::shared_ptr<Subcomponent> comp);
        const std::vector<Action> SCCrossingActions(std::shared_ptr<Subcomponent> SCfrom,std::shared_ptr<Subcomponent>SCto  );
        //energy of an action calculated by redecomposition
        double rebuidEnergy(int i, int j);
        //enumerates actions betveen first and second, which will result in a pseudoknot
        std::vector<Action> makePseudo(std::shared_ptr<Subcomponent> first ,std::shared_ptr<Subcomponent> second);

        bool canFormPseudo(std::shared_ptr<Subcomponent> sc1,std::shared_ptr<Subcomponent> sc2 );
        ///////////////////////////////////////////
        void updateActions(Action lastAction);

        //saves inner actions to  dim1 and updates
        void saveSCActions(std::shared_ptr<Subcomponent> sc);

        //saves crossing actions to dim1 and updates map3 correspondingly
        //assumes that sc and asc are already saved (using saveSCActions)
        void saveCrossingActions(std::shared_ptr<Subcomponent> sc, std::shared_ptr<Subcomponent> asc);
        /////////////////////////////////////////////
};