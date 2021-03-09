//cpp file for the secondary structure of RNA

#include "RNASecondaryStructure.hpp"
#include "CRTree.hpp"
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>
#include "../utility/EnergyCalc.hpp"
#include <deque>
#include <unordered_set>

RNASecondaryStructure::RNASecondaryStructure(RNAPrimarySequence prima, std::string _dotBracket,std::vector<int> psnum) :  primary(prima),dotBracket(_dotBracket) {
    assert( prima.getLength()== (int)_dotBracket.size());
    //convert dotbracket to bp vector
    
    std::vector<int> stack;
    
    //  ...(((...[[[...[[[)))...]]]...]]]
    //            0     1        0     1
    std::vector<std::vector<int>> pstack(psnum.size()/2);
    int i;
    int pi =0;
    for( int e=0; e < (int)_dotBracket.size(); e++){
       switch (_dotBracket[e])
       {
        case '(':
            bp.push_back(-1);
            stack.push_back(e);
            break;

        case ')':
            bp.push_back(-1);
            i=stack.back();
            stack.pop_back();
            bp[i]=e;
            bp[e]=i;
            
            break;
        case '[':
                    while(_dotBracket[e]=='['){

                        bp.push_back(-1);
                        pstack[psnum[pi]].push_back(e);
                        e++;
                    }
                    pi++;
                    e--;
                    break;
        case ']':
                    while(_dotBracket[e]==']'){

                        bp.push_back(-1);
                        i=pstack[psnum[pi]].back();
                        pstack[psnum[pi]].pop_back();
                        bp[i]=e;
                        bp[e]=i;
        
                        
                        e++;
                    }
                    pi++;
                    e--;
                   
                    break;
        default:
            bp.push_back(-1);
            break;
        }

    }
    assert(stack.empty());
 for( int e=0; e < (int) _dotBracket.size(); e++){
        if(bp[e] > -1){
            bool NOT_WC_OR_WOBBLE = (primary.isPairable(e,bp[e]));
            if (!NOT_WC_OR_WOBBLE) std::cout<<e<<" is not pairable with "<<bp[e] <<std::endl;
            assert(NOT_WC_OR_WOBBLE);
        }

    
    }

    
    root = std::make_shared<Subcomponent>( Subcomponent(Root,std::vector<int>{-1,primary.getLength()}));
    docomposeStructure();
    evaluateEnergy();
}
const std::string & RNASecondaryStructure::getDB(){
    return dotBracket;
}

void RNASecondaryStructure::docomposeStructure(std::shared_ptr<Subcomponent> froot, int start ,int end, std::function<bool(int)> region) {

    if(!froot){
        root = std::make_shared<Subcomponent>( Subcomponent(Root,std::vector<int>{-1,primary.getLength()}));
        froot=root;}
    if(end ==-1) end = primary.getLength();

    // building the clsed region tree
    CRTree T(Node(-1,bp.size()));
    std::vector<std::tuple<int,int>> stack;

    // std::cout<<"startin decomposition"<< bp.size() << std::endl;
    for( int base = start; base< end; base++){
        if(region(base)){// just test
            if(base< bp[base] ) stack.push_back(std::make_tuple(base,bp[base]));
            else if (-1 < bp[base] && bp[base] <base){
                int E = base;
                while(!stack.empty() && std::get<0>(stack.back()) > bp[base]){
                    E= std::max(E,std::get<1>(stack.back()));
                    stack.pop_back();
                }
                if(!stack.empty()) std::get<1>(stack.back()) = std::max(E,std::get<1>(stack.back()));

            }
            if(!stack.empty() && base == std::get<1>(stack.back())){
                const int i = std::get<0>(stack.back());
                const int j = std::get<1>(stack.back());
                stack.pop_back();
                Node n=Node(i,j);
                // std::cout<<"adding node "<<i<<"  "<< j << std::endl;
                T.addToTree(n); 
            } 
        }
    }
    
    //using the closed regien tree to make subcomponent tree 

    std::vector< std::tuple<std::shared_ptr<Node>,std::shared_ptr<Subcomponent>> > stack2;
    for(int i=0;i<(int) T.getRoot()->children.size(); i++){
        stack2.push_back( std::make_tuple(T.getRoot()->children[i],froot));
    }
    
    

    while(!stack2.empty()){
        auto currentstart = stack2.back();
        stack2.pop_back();
        auto currCRTS = std::get<0>(currentstart);
        auto currSCTS = std::get<1>(currentstart);

        // determine and add subcomponents

        

            //interior, bulge or stem
            if(currCRTS->children.size() ==1 && currCRTS->left == bp[currCRTS->right]){ 
                auto child = currCRTS->children[0];
                while(currCRTS->children.size() ==1  && currCRTS->children[0]->left == bp[currCRTS->children[0]->right]){ // child->left ==bp[child->right] if fails "multi loop" with pseudochild
                    child = currCRTS->children[0];
                    //interior
                    if(currSCTS->getType()==Stem && currCRTS->left+1 < child->left && currCRTS->right-1 > child->right  ){
                        auto comp=Subcomponent(Interior,std::vector<int> {currCRTS->left +1,child->left -1, child->right +1, currCRTS->right -1});
                        auto SCChild = std::make_shared<Subcomponent>(comp);
                        currSCTS->addChild(SCChild);
                        currSCTS = SCChild;
                        currCRTS = child;
                    }
                    //bulge
                    else if(currSCTS->getType()==Stem && (currCRTS->left+1 < child->left) != (currCRTS->right-1 > child->right)){
                        int k,l;
                        if (currCRTS->left +1 == child->left) {
                            l=currCRTS->right -1;
                            k=child->right +1;
                        }else{
                            k=currCRTS->left +1;
                            l=child->left -1;
                        }
                        auto comp=Subcomponent(Bulge,std::vector<int> { k, l});
                        auto SCChild = std::make_shared<Subcomponent>(comp);
                        currSCTS->addChild(SCChild);
                        currSCTS = SCChild;
                        currCRTS = child;
                    }
                    //stem
                    else {
                        auto start = currCRTS;
                        auto end = currCRTS; //child;
                        while( currCRTS->left+1 == child->left && currCRTS->right-1 == child->right && child->left == bp[child->right]){
                            end = child;
                            currCRTS = child;
                            if(currCRTS->children.empty()) break;
                            else child = child->children[0];
                            
                        }
                        auto comp=Subcomponent(Stem,std::vector<int> {start->left ,end->left, end->right, start->right});
                        auto SCChild = std::make_shared<Subcomponent>(comp);
                        currSCTS->addChild(SCChild);
                        currSCTS = SCChild;
                        currCRTS = end;

                    }
                }
                

            }
            //hairpin loop, or Multi loop (or lone pair before)
            if(currCRTS->left == bp[currCRTS->right]) {

                //lone pair before multiloop or hairpin
                if(currSCTS->getType()!=Stem){

                        auto comp=Subcomponent(Stem,std::vector<int> {currCRTS->left ,currCRTS->left, currCRTS->right, currCRTS->right});
                        auto SCChild = std::make_shared<Subcomponent>(comp);
                        currSCTS->addChild(SCChild);
                        currSCTS = SCChild;
                }

                //hairpin loop
                if(currCRTS->children.empty()){
                    auto comp=Subcomponent(Hairpin,std::vector<int> {currCRTS->left+1 ,currCRTS->right -1});
                    auto SCChild = std::make_shared<Subcomponent>(comp);
                    currSCTS->addChild(SCChild);
                }
                //multi loop
                else{
                    // make multi loop
                    std::vector<int> defmulti;
                    int i=currCRTS->left +1;
                    int childn = 0;;
                    int siz =currCRTS->children.size();
                    
                    while(currCRTS->right >i ){
                        if(childn < siz && i==currCRTS->children[siz-1-childn]->left){
                            i= currCRTS->children[siz-1-childn]->right +1;
                            childn++;
                        }else {
                            defmulti.push_back(i);
                            i++;
                        }
                    }
                    
                        auto comp=Subcomponent(Multiloop,defmulti);
                        auto SCChild = std::make_shared<Subcomponent>(comp);
                        currSCTS->addChild(SCChild);

                    //add children of multi loop to stack2 for continuing 
                    for( int i=0; i<(int) currCRTS->children.size();i++){
                        stack2.push_back( std::make_tuple(currCRTS->children[i] ,SCChild));
                    }
                }
                
            }
            // pseudoknot 
            else {
                //calculate type of pseudoknot 
                std::vector<int> pvec{currCRTS->left};
                int siz = currCRTS->children.size();
                int childn =0;
                
                int curr=currCRTS->left;
                int currend=bp[curr];
                int f1= bp[curr];
                int f2= -2;
                for(int i=currCRTS->left; i<=currCRTS->right; i++){
                   assert(i<primary.getLength() && i>=0);
                    if(childn < siz && i==currCRTS->children[siz-1-childn]->left){
                        i= currCRTS->children[siz-1-childn]->right;
                        childn++;
                        
                    }else if(bp[i]>currend){
                        pvec.push_back(curr);
                        curr=i;
                        pvec.push_back(i);
                        currend = bp[i];
                        f2=bp[i];
                    }else if(bp[i]>i){
                        curr=i;
                    }  else if(bp[i] != -1){
                        
                            pvec.push_back(curr);
                            pvec.push_back(i);
                            i = f1;
                            curr = f1;
                            f1=f2;
                            f2= -2;
                        
                    }
                }
                pvec.push_back(currCRTS->right);

                /* int i = currCRTS->left;
                int siz = currCRTS->children.size();
                int childn =0;
                int paired =1;
                //end of current arc
                int currend=bp[i];
                while( i<currCRTS->right ){
                    if(childn < siz && i==currCRTS->children[siz-1-childn]->left){
                        i= currCRTS->children[siz-1-childn]->right +1;
                        childn++;
                    }else {   
                     if(bp[i]>currend){
                            paired++;
                            int pairedidxb = currend;
                            currend=bp[i];
                            i=pairedidxb;
                        }
                        i++;
                    }
                }
                // only H and K type allowed (for now)
                assert(paired <4);
                SCType typ = paired==2 ? PseudoknotH : PseudoknotK; */
                assert(pvec.size()==8 || pvec.size() == 12);
                SCType typ = pvec.size()==8 ? PseudoknotH : PseudoknotK; 
                
                auto comp=Subcomponent(typ,pvec);
                auto SCChild = std::make_shared<Subcomponent>(comp);
                currSCTS->addChild(SCChild);
                //add children of pseudoknot to stack2 for continuing 
                for( int i=0; i<(int) currCRTS->children.size();i++){
                    stack2.push_back( std::make_tuple(currCRTS->children[i] ,SCChild));
                }                
                
            }
        

        
  
    }
    //add  exterior
    if(root==froot)
    {


        std::shared_ptr<Subcomponent> SCChild1 = nullptr;
        std::shared_ptr<Subcomponent> SCChild2 = nullptr;
        // "end" exterior loop
        if(!root->getChildren().empty()){

            int end =root->getChildren()[0]->getDefiningBases().back() +1;
            if(end < primary.getLength()) {
                auto comp=Subcomponent(Exterior,std::vector<int> {end,primary.getLength()-1});
                SCChild1 = std::make_shared<Subcomponent>(comp);
            }
            int start = root->getChildren().back()->getDefiningBases()[0] -1;   
            if(start > -1 ){
                auto comp=Subcomponent(Exterior,std::vector<int> {0,start});
                SCChild2 = std::make_shared<Subcomponent>(comp);         
            }
        } else {
            auto comp=Subcomponent(Exterior,std::vector<int> {0,primary.getLength()-1});
                SCChild1 = std::make_shared<Subcomponent>(comp); 
        }
    

        // "middle" exterior
        std::vector<std::shared_ptr<Subcomponent>> middle;
        int siz =  root->getChildren().size()-1;
        for(int i=0;i< siz; i++){
            int start =root->getChildren()[i+1]->getDefiningBases().back() +1;
            int end  = root->getChildren()[i]->getDefiningBases()[0] -1;
            if(start <= end){
                auto comp=Subcomponent(Exterior,std::vector<int> {start,end});
                auto SCChild = std::make_shared<Subcomponent>(comp);      
                //root->addChild(SCChild);
                middle.push_back(SCChild);
            }
            
        }

        //adding the "end" exterior
        if(SCChild1) root->addChild(SCChild1);
        
        
        for(auto c : middle){
            root->addChild(c);
        } 
        
        if(SCChild2) root->addChild(SCChild2);
    }

    // std::cout<<"decdone"<<std::endl;
}
const std::string & RNASecondaryStructure::getDotBracket(){
    std::vector<int> rb {primary.getLength()};
    for(int i=0; i<= std::min(primary.getFoldingLength(),primary.getLength()-1); i++){
        if(bp[i]==-1){
            dotBracket[i]='.';
        }else if(i<bp[i] ){
            if(bp[i]<rb.back()) {
                dotBracket[i]='(';
                dotBracket[bp[i]]=')';
                rb.push_back(bp[i]);
            }
            else {
                dotBracket[i]='[';
                dotBracket[bp[i]]=']';
            }
        }else if (dotBracket[i]==')'){
            rb.pop_back();
        }
    }
    for(int i =primary.getFoldingLength()+1;i<primary.getLength();i++ ){
        dotBracket[i]='.';
    }
    return dotBracket;
}

double RNASecondaryStructure::evaluateEnergy(){
    energy =evaluateSubstructureEnergy(root);
    return energy;
}

double RNASecondaryStructure::evaluateSubstructureEnergy(std::shared_ptr<Subcomponent> subRoot ) {
    double energyTot = 0;
    std::vector<std::shared_ptr<Subcomponent>> stack;
    stack.push_back(subRoot);
    while(!stack.empty()){

        auto curr = stack.back();
        stack.pop_back();
        double e =componentEnergy(curr);
        energyTot += e;
        // std::cout<<"type "<< curr->getType()<< " energy "<< e<< std::endl;

        for(auto c :  curr->getChildren()){
            stack.push_back(c);
        }

        
    }

    return energyTot;
}



double RNASecondaryStructure::componentEnergy(std::shared_ptr<Subcomponent> comp, int ii, int jj) {
    // here a as lambda instead of regular vector so that action energy can be calculated NOTE not yet in pkenergy 
    auto tmp2= bp;
    bool moi = ii!=- 1;
    if(moi){
        if(bp[ii]==jj){
            bp[jj]=-1;
            bp[ii]=-1;
        }else{
            bp[ii]=jj;
            bp[jj]=ii;
        }
    }
    //exactly same as bp
    std::function<int(int)> bpf = [&bp=bp,ii,jj](int i)->int{ 
        if(true ||ii==-1) return bp[i];
        else if (i==ii ){
            if ( bp[i] > -1) return -1;
            else return jj;

        }else if (i==jj ){
            if ( bp[i] > -1) return -1;
            else return ii;

        }else return bp[i];
    };

    SCType ty= comp->getType();
    double energy = 0;
    auto vec = comp->getDefiningBases();


    switch (ty)
    {
    case Stem:
        { 
            int i = 0;
            while( i+vec[0]< vec[1]){
                char ii = primary.getString()[vec[0]+i];
                char pp = primary.getString()[vec[0]+i+1];
                char jj = primary.getString()[vec[3]-i];
                char qq = primary.getString()[vec[3]-i-1];
                energy += EnergyCalculator::stemPairEnergy(ii,jj,pp,qq);
                i++;
            }

            // AU or GU ends  
            energy += EnergyCalculator::endPenalty(primary.getString()[vec[0]],primary.getString()[vec[3]]);
            energy += EnergyCalculator::endPenalty(primary.getString()[vec[1]],primary.getString()[vec[2]]);
        }
        break;

    case Interior:
        {
            int left = vec[1] -vec[0] +1;
            int right = vec[3] -vec[2] +1;
            int diff= std::abs(right-left);
            int n= right+left;
            assert(n>1);
            if (n == 2) {
                energy+= EnergyCalculator::internalLoop11Energy(primary.getString()[vec[0]-1],primary.getString()[vec[3]+1], primary.getString()[vec[1]+1],primary.getString()[vec[2]-1], primary.getString()[vec[0]],primary.getString()[vec[3]]);
                if(EnergyCalculator::endPenalty(primary.getString()[vec[0]-1],primary.getString()[vec[3]+1])>0.0) energy-=  EnergyConstants::AU_end;
                if(EnergyCalculator::endPenalty(primary.getString()[vec[1]+1],primary.getString()[vec[2]-1])>0.0) energy-=  EnergyConstants::AU_end;
                break;
                return energy;
            }


            energy+= (n <= EnergyConstants::maxLoopSize) ? (EnergyConstants::internal_loop[n]) : (EnergyConstants::internal_loop[EnergyConstants::maxLoopSize] + EnergyConstants::lxc37 * std::log((n) / (double) EnergyConstants::maxLoopSize));
            energy+= diff*EnergyConstants::internalAsymmetry;
            //special internal end penalty
            //   AU/GU closure 
            if(EnergyCalculator::endPenalty(primary.getString()[vec[0]-1],primary.getString()[vec[3]+1])>0.0) energy+= EnergyConstants::internalAUGUClosure - EnergyConstants::AU_end;
            if(EnergyCalculator::endPenalty(primary.getString()[vec[1]+1],primary.getString()[vec[2]-1])>0.0) energy+= EnergyConstants::internalAUGUClosure - EnergyConstants::AU_end;

            energy += EnergyCalculator::internalmmEnergy(primary.getString()[vec[0]],primary.getString()[vec[3]],primary.getString()[vec[1]],primary.getString()[vec[2]]);
        
        }

        break;

    case Bulge:
        {
            int n = vec[1] -vec[0] +1;
            if( n==1) {
                // check if special C only in 2004
                //if((primary.getString()[vec[0]] == 'C') &&(primary.getString()[vec[0]-1] == 'C' || primary.getString()[vec[0]+1] == 'C' ) )   energy += EnergyConstants::bulge[0];
                // add stacking energy 
                int i,j,p,q;
                if (vec[0]-1 < bpf(vec[0]-1)){
                    i = vec[0]-1;
                    j = bpf(vec[0]-1);
                    p = vec[0]+1;
                    q = bpf(vec[0]+1);
                     
                }else{
                    q = vec[0]-1;
                    p = bpf(vec[0]-1);
                    j = vec[0]+1;
                    i = bpf(vec[0]+1);

                }
                char ii = primary.getString()[i];
                char pp = primary.getString()[p];
                char jj = primary.getString()[j];
                char qq = primary.getString()[q];
                energy += EnergyCalculator::stemPairEnergy(ii,jj,pp,qq);
                // No end penalty if "stacking"
                energy -= EnergyCalculator::endPenalty(primary.getString()[vec[0]-1],primary.getString()[bpf(vec[0]-1)]);
                energy -= EnergyCalculator::endPenalty(primary.getString()[vec[1]+1],primary.getString()[bpf(vec[1]+1)]);

            } 
            
        
           energy += (n <= EnergyConstants::maxLoopSize) ? (EnergyConstants::bulge[n]) : (EnergyConstants::bulge[EnergyConstants::maxLoopSize] + EnergyConstants::lxc37 * std::log(n / (double) EnergyConstants::maxLoopSize));

        }
        break;

    case Multiloop:
   
        {
            int unpaired = vec.size();


            int paired=1; //closing pair 

            //closing pair helix
            int i= comp->getLeft() -1;//std::min(vec[0] -1 , comp->getChildren()[0]->getDefiningBases()[0] -1 );
            
            energy += terminalEnergy(bpf(i),i);

            for( auto c : comp->getChildren()){
        
                SCType t = c->getType();
                if(t == PseudoknotH || t == PseudoknotK){
                    // adding pseudochild paired bases 
                    paired += 2; // t -6;
                    // adding childs initialisation energy
                    energy += EnergyConstants::pseudoknot[t-8][1] -EnergyConstants::pseudoknot[t-8][0];

                    //DP09
                    // energy += EnergyConstants::pseudoknotDP09[1]-EnergyConstants::pseudoknotDP09[0];

                    //pseudo helices terminal energy
                    //ONLY DANGLE
                    // energy += terminalEnergy(c->getDefiningBases()[0],bpf(c->getDefiningBases()[0]));
                    // energy += terminalEnergy(bpf(c->getDefiningBases().back()),c->getDefiningBases().back());

                }
                else {
                    energy += terminalEnergy(c->getDefiningBases()[0],c->getDefiningBases().back());
                    paired++;
                }
                
            }


            energy += EnergyCalculator::multiloopEnergy(unpaired,paired);
        }
        break;

    case Exterior:
        {
        }
        
        break;

    case Pseudoknot:
        {
            std::cout<< "type generic Pseudoknot, Energy will not be calculated" <<std::endl;
            energy = pkEnergy(comp);
        }
        break;

    case PseudoknotH:
        energy = pkEnergy(comp);
        // energy = pkEnergyDP09(comp);
        break;

    case PseudoknotK:
        energy = pkEnergy(comp);
        // energy = pkEnergyDP09(comp);
        break;

    case Hairpin:
        {
            int n = vec[1] -vec[0] +1;
            assert(n>2);

            auto str = primary.getString().substr(vec[0]-1,n+2);

            // if special hairpin return special hairpin energy ONLY IN 2004
            // double spec = EnergyCalculator::hairpinSpecialEnergy(str);
            //if(spec >0.0) return spec;

            

            //loop size energy
            energy += (n <= EnergyConstants::maxLoopSize) ? (EnergyConstants::hairpin[n]) : (EnergyConstants::hairpin[EnergyConstants::maxLoopSize] + EnergyConstants::lxc37 * std::log(n / (double) EnergyConstants::maxLoopSize));
            
            auto loopstr = str.substr(1,n);
            // all c penalty
            bool all_C_Loop = true;
            for (auto c : loopstr) {
                if (c != 'C') {
                    all_C_Loop = false;
                    break;
                }
            }

            if (all_C_Loop) {
                if (n == 3) energy += EnergyConstants::hairpinC3;
                else energy +=  n * EnergyConstants::hairpinAllCA + EnergyConstants::hairpinAllCB;      
            }

            if(n>3){
                // mismatch
                energy += EnergyCalculator::terminalMismatchEnergy(str[0],str.back(),loopstr[0],loopstr.back());
                // special GU closure
                if (str[0]== 'G' && str.back() == 'U') energy += EnergyConstants::hairpinGUclosure;
                // UU or GA first mismatch
                if((loopstr[0] == 'U' && loopstr.back() == 'U') || (loopstr[0] == 'G' && loopstr.back() == 'A')) energy += EnergyConstants::hairpinUUGAmm;

                // GG fist mismatch  ONLY IN 2004
                //if(loopstr[0] == 'G' && loopstr.back() == 'G') energy += EnergyConstants::hairpinGGmm;

            }

            if(n==4) energy += EnergyCalculator::hairpinTetraEnergy(str);
           
        }
        break;

        case Root:
        {
            for(auto c :  comp->getChildren()){
                SCType typ= c->getType();
                if(typ==Exterior) { }
                else if(typ==Stem){
                    
                   energy += terminalEnergy(c->getDefiningBases()[0],c->getDefiningBases().back());
                }else if(typ==PseudoknotH || typ == PseudoknotK){
                    
                    // energy += terminalEnergy(c->getDefiningBases()[0],bpf(c->getDefiningBases()[0]));
                    // energy += terminalEnergy(bpf(c->getDefiningBases().back()),c->getDefiningBases().back());
                    //ONLY DANGLE 
                }else{
                    assert(typ==Stem ||typ==PseudoknotH || typ== PseudoknotK);
                }
            }
        }
        break;


    default:
        std::cout<<"WARNING Undefined component energy = 0"<<std::endl;
        break;
    }

    if(ii!=-1){
        if(bp[ii]==jj){
            bp[jj]=-1;
            bp[ii]=-1;
        }else{
            bp[ii]=jj;
            bp[jj]=ii;
        }
    }
    assert(tmp2 == bp);
    return energy;
}


double RNASecondaryStructure::pkEnergy(std::shared_ptr<Subcomponent> comp){
//remenber part in the multiloop
    SCType ty= comp->getType();
    double energy = 0;
    auto vec = comp->getDefiningBases();
    assert(ty ==Pseudoknot || ty==PseudoknotH || ty==PseudoknotK);

    int i= vec[0];
    //number of next child 
    int childn=0;
    // size of children
    int siz = comp->getChildren().size();
    //number of unpaired bases in pseudoknot: all unpaired which do not belong to a child
    int unpaired=0;
    //number of base pairs in pseudoknot (2 if H-type and 3 if K-type)
    int paired = ty==PseudoknotH ? 2 : 3;

    
    //calculationg number of unpaired 
    // while(vec.back() >i ){
    //         if(childn < siz && i==comp->getChildren()[childn]->getDefiningBases()[0]){
    //             i= comp->getChildren()[childn]->getDefiningBases().back() +1;
    //             childn++;
    //         }else if(bp[i] < i){ 
    //             if(bp[i] ==-1){
    //                 unpaired++;
                
    //             }
    //             i++;
    //         }else {
    //             i++;
    //         }
    //     }
   // int typ= paired-2;

   int typ = ty - 8;

   // int pkchildren= std::count(comp->getChildren().begin(),comp->getChildren().end(), [](std::shared_ptr<Subcomponent> c){c->getType() ==Pseudoknot;}) ;
    int eend = ty==PseudoknotH ? 5 : 9;
    for(int idx=1; idx<=eend; idx++){
        unpaired += vec[idx+1] -vec[idx] -1;
    }

    for( auto c : comp->getChildren()){
        unpaired-= c->getRight() - c->getLeft() +1;
        SCType t = c->getType();
        if(t == PseudoknotH || t == PseudoknotK){
            // adding pseudochild paired bases 
            paired += t -6;
            // adding childs initialisation energy
            energy += EnergyConstants::pseudoknot[t-8][2] -EnergyConstants::pseudoknot[t-8][0];
            //MISMATCH ENERGY HERE !!!
            //look implementation from multiloop
        }
        else {
            //Mismatch energy of children
            energy += terminalEnergy(c->getDefiningBases()[0],c->getDefiningBases().back());
            paired++;
        }
        
    }

    
    //stacking energy
    std::vector<std::shared_ptr<Subcomponent>> pseudostem;
        
    //pseudoknots  consists of pure stems 
    //stems added to pseudoknot energy
    auto comp2=Subcomponent(Stem,std::vector<int> {vec[0],vec[1],vec[4],vec[5]});
    auto SCChild = std::make_shared<Subcomponent>(comp2);      
    pseudostem.push_back(SCChild);

    if(ty==PseudoknotH){
        comp2=Subcomponent(Stem,std::vector<int> {vec[2],vec[3],vec[6],vec[7]});
        SCChild = std::make_shared<Subcomponent>(comp2);      
        pseudostem.push_back(SCChild);
    }else{
        comp2=Subcomponent(Stem,std::vector<int> {vec[2],vec[3],vec[8],vec[9]});
        SCChild = std::make_shared<Subcomponent>(comp2);      
        pseudostem.push_back(SCChild);
        comp2=Subcomponent(Stem,std::vector<int> {vec[6],vec[7],vec[10],vec[11]});
        SCChild = std::make_shared<Subcomponent>(comp2);      
        pseudostem.push_back(SCChild);
    }
        
    
    for(auto e : pseudostem){
        energy += componentEnergy(e);
    }

    // pseudoknot initialisation energy
    energy += EnergyConstants::pseudoknot[typ][0];
    energy += unpaired * EnergyConstants::pseudoknot[typ][4] + paired * EnergyConstants::pseudoknot[typ][3];

    return energy;

}


double RNASecondaryStructure::pkEnergyDP09(std::shared_ptr<Subcomponent> comp){
//remenber part in the multiloop
    SCType ty= comp->getType();
    double energy = 0;
    auto vec = comp->getDefiningBases();
    assert(ty ==Pseudoknot || ty==PseudoknotH || ty==PseudoknotK);

    int i= vec[0];
    //number of next child 
    int childn=0;
    // size of children
    int siz = comp->getChildren().size();
    //number of unpaired bases in pseudoknot: all unpaired which do not belong to a child
    int unpaired=0;
    //number of base pairs in pseudoknot ( nested substructures)
    int paired = 0;

    
    //calculationg number of unpaired 
    // while(vec.back() >i ){
    //         if(childn < siz && i==comp->getChildren()[childn]->getDefiningBases()[0]){
    //             i= comp->getChildren()[childn]->getDefiningBases().back() +1;
    //             childn++;
    //         }else if(bp[i] < i){ 
    //             if(bp[i] ==-1){
    //                 unpaired++;
                
    //             }
    //             i++;
    //         }else {
    //             i++;
    //         }
    //     }
   // int typ= paired-2;

   int typ = ty - 8;

   // int pkchildren= std::count(comp->getChildren().begin(),comp->getChildren().end(), [](std::shared_ptr<Subcomponent> c){c->getType() ==Pseudoknot;}) ;
    int eend = ty==PseudoknotH ? 5 : 9;
    for(int idx=1; idx<=eend; idx++){
        unpaired += vec[idx+1] -vec[idx] -1;
    }

    for( auto c : comp->getChildren()){
        unpaired-= c->getRight() - c->getLeft() +1;
        SCType t = c->getType();
        if(t == PseudoknotH || t == PseudoknotK){
            // adding pseudochild paired bases 
            paired += t -6;
            // adding childs initialisation energy
            energy += EnergyConstants::pseudoknotDP09[2]-EnergyConstants::pseudoknotDP09[0];
            //MISMATCH ENERGY HERE !!!
            //look implementation from multiloop
        }
        else {
            //Mismatch energy of children
            energy += terminalEnergy(c->getDefiningBases()[0],c->getDefiningBases().back());
            paired++;
        }
        
    }

    
    //stacking energy
    std::vector<std::shared_ptr<Subcomponent>> pseudostem;
        
    //pseudoknots  consists of pure stems 
    //stems added to pseudoknot energy
    auto comp2=Subcomponent(Stem,std::vector<int> {vec[0],vec[1],vec[4],vec[5]});
    auto SCChild = std::make_shared<Subcomponent>(comp2);      
    pseudostem.push_back(SCChild);

    //minimum of 2 bands
    energy += EnergyConstants::pseudoknotDP09[3]*2;

    if(ty==PseudoknotH){
        comp2=Subcomponent(Stem,std::vector<int> {vec[2],vec[3],vec[6],vec[7]});
        SCChild = std::make_shared<Subcomponent>(comp2);      
        pseudostem.push_back(SCChild);
    }else{
        //type K has 3 bands
        energy += EnergyConstants::pseudoknotDP09[3];

        comp2=Subcomponent(Stem,std::vector<int> {vec[2],vec[3],vec[8],vec[9]});
        SCChild = std::make_shared<Subcomponent>(comp2);      
        pseudostem.push_back(SCChild);
        comp2=Subcomponent(Stem,std::vector<int> {vec[6],vec[7],vec[10],vec[11]});
        SCChild = std::make_shared<Subcomponent>(comp2);      
        pseudostem.push_back(SCChild);
    }
        
    
    for(auto e : pseudostem){
        energy += componentEnergy(e)* EnergyConstants::pseudoknotDP09[9];
    }

    // pseudoknot initialisation energy
    energy += EnergyConstants::pseudoknotDP09[0];
    energy += unpaired * EnergyConstants::pseudoknotDP09[4] + paired * EnergyConstants::pseudoknotDP09[5];

    return energy;

}

//dangling energies disabled !!!!!!
double RNASecondaryStructure::terminalEnergy( const int i ,const int j){
// i.j i is on 5' and j on 3'
/*    terminal
        | 
        v
    5' jX 3'
       |
    3' iY 5'
NOTE closing pair of multiloop  will be inserted so that j<i 
 */   
//  assert(i== bp[j]); // Can not be asserted because will fail when enumerating actions
 assert(j>=0 && j<primary.getLength());
 assert(i>=0 && i<primary.getLength());
    double energy = 0;
    
    if(j == primary.getLength()-1 && i == 0){
        return 0;
    }
    else if (i==0){
        // dangling 3'
        // energy = EnergyCalculator::danglingEnergy(primary.getString()[i],primary.getString()[j],primary.getString()[j+1],true);
    }
    else if(j == primary.getLength()-1){
        // dangling 5'
        // energy = EnergyCalculator::danglingEnergy(primary.getString()[i],primary.getString()[j],primary.getString()[i-1],false);
    }
    else{
        // mismatch
        if(bp[i-1] == -1 && bp[j+1]== -1){
            energy = EnergyCalculator::terminalMismatchEnergy(primary.getString()[i],primary.getString()[j],primary.getString()[i-1],primary.getString()[j+1]);
        } 
        // dangling 3'
        else if (bp[j+1] == -1){
            // energy = EnergyCalculator::danglingEnergy(primary.getString()[i],primary.getString()[j],primary.getString()[j+1],true);

        }
        // dangling 5'
        else if(bp[i-1] ==-1){
            // energy = EnergyCalculator::danglingEnergy(primary.getString()[i],primary.getString()[j],primary.getString()[i-1],false);
        }

        
    }
    return energy;
}

void RNASecondaryStructure::executeAction(Action & ac){
    assert(ac.j >=0 && ac.j < primary.getLength());
    assert(ac.i >=0 && ac.i < primary.getLength());
    // assert(ac.i < ac.j);
    assert(primary.isPairable(ac.i,ac.j));
    assert((ac.getType() == Addition && bp[ac.i]==-1 && bp[ac.j]==-1) || (ac.getType()==Deletion && bp[ac.i]==ac.j));
    //saving old tree
    oldRoot = root;
    if(bp[ac.i]==ac.j){
        bp[ac.j]=-1;
        bp[ac.i]=-1;
    }else{
        bp[ac.i]=ac.j;
        bp[ac.j]=ac.i;
    }
    
    docomposeStructure();
    double ee =energy + ac.getEnergy();
    double ee2 = evaluateSubstructureEnergy(root);
    assert((float)ee ==(float) ee2 ); 
    energy=ee;
    
}
// OLD enumerates in  wrong order
std::vector<Action> RNASecondaryStructure::enumerateActions(){

    //stack for "root" subcomponents i.e. Root, Multiloop or PseudoknotH/K
    std::vector<std::shared_ptr<Subcomponent>> rStack;
    //stack for siblings of current subcomponent
    // std::vector<std::shared_ptr<Subcomponent>> sStack;
    // parent of the current subcomponent
    //std::shared_ptr<Subcomponent> parent;
    //current subcomponent
    std::shared_ptr<Subcomponent> current;
    //parent for subcomponents in sStack
    std::shared_ptr<Subcomponent> rParent;
    
    std::vector<Action> allActions;

    rStack.push_back(root);

    while(!rStack.empty()){
        
        rParent =rStack.back();
        rStack.pop_back();

        // for pseudoknots siblings will be divided into regions 
        std::vector<std::vector<std::shared_ptr<Subcomponent>> >sRegions;

        if(rParent->getType() == PseudoknotK ){
            for(int idx=1;idx<=9;idx +=2){
                auto tmp = rParent->getChildren(rParent->getDefiningBases()[idx],rParent->getDefiningBases()[idx+1]);
                sRegions.push_back(tmp);
            }
        }
        else if( rParent->getType() == PseudoknotH){
            for(int idx=1;idx<=5;idx +=2){
                auto tmp = rParent->getChildren(rParent->getDefiningBases()[idx],rParent->getDefiningBases()[idx+1]);
                sRegions.push_back(tmp);
            }
        }
        else{

            std::vector<std::shared_ptr<Subcomponent>> sStack;
            for( auto c : rParent->getChildren()){
                sStack.push_back(c);
            }
            sRegions.push_back(sStack);
        }
        
        

        for( auto sStack : sRegions)
        while(!sStack.empty()){
            //Exterior Stem or Pseudoknot enumerate action for this SC
            current=sStack.back();
            sStack.pop_back();
            if(current->getType()==Exterior){
                // inside and pseudoneighbours
                // siblings (including descendants) exteriors 
                auto tmp=subcomponentActions(current,sStack);
                allActions.insert(allActions.end(), tmp.begin(), tmp.end());

            }else if(current->getType()==Stem){
                //all delete options
                auto tmp=subcomponentActions(current,std::vector<std::shared_ptr<Subcomponent>> { });
                allActions.insert(allActions.end(), tmp.begin(), tmp.end());
            }
            // NOTE Exterior loops are always evaluated before pseudo or stem 
            // psudoknots cant form base pairs with other siblings
            // so no need to look for base pairs to siblings
            else if(current->getType()==PseudoknotH){
               // inside and descendants
                auto tmp=subcomponentActions(current,std::vector<std::shared_ptr<Subcomponent>> { });
                allActions.insert(allActions.end(), tmp.begin(), tmp.end());
               // add to rStack
               rStack.push_back(current);

            }else if(current->getType()==PseudoknotK){
                // inside and descendants
                auto tmp=subcomponentActions(current,std::vector<std::shared_ptr<Subcomponent>> { });
                allActions.insert(allActions.end(), tmp.begin(), tmp.end());
                // add to rStack
                rStack.push_back(current);

            }

            //first subcomponent level untill hairpin loop or multiloop 
            //can not be pseudoknot 
            if(current->getType() == Stem){

                current = current->getChildren()[0];
                if(current->getType()==Multiloop){
                        auto tmp=subcomponentActions(current,sStack);
                        allActions.insert(allActions.end(), tmp.begin(), tmp.end());
                        rStack.push_back(current);
                        
                    }

                    else{
                        auto tmp=subcomponentActions(current,sStack);
                        allActions.insert(allActions.end(), tmp.begin(), tmp.end());

                        while(current->getChildren().size()==1){
                        
                            current = current->getChildren()[0];

                            if(current->getType()==Multiloop){
                                auto tmp=subcomponentActions(current,std::vector<std::shared_ptr<Subcomponent>> { });
                                allActions.insert(allActions.end(), tmp.begin(), tmp.end());
                                rStack.push_back(current);
                                break;
                            }

                            else{
                                auto tmp=subcomponentActions(current,std::vector<std::shared_ptr<Subcomponent>> { });
                                allActions.insert(allActions.end(), tmp.begin(), tmp.end());
                            } 


                                        
                                    
                        }
                    } 

                
            }
        }

    }
    
    return allActions;
}

// new version of enumerating action with FIFO queue not STACK
std::vector<Action> RNASecondaryStructure::enumerateActionsOrder(){

    //stack for "root" subcomponents i.e. Root, Multiloop or PseudoknotH/K
    std::vector<std::shared_ptr<Subcomponent>> rStack;
    //queue for "root" subcomponents i.e. Root, Multiloop or PseudoknotH/K
    std::deque<std::shared_ptr<Subcomponent>> rQueue;
    
    
    //current subcomponent
    std::shared_ptr<Subcomponent> current;
    //parent for subcomponents in sStack
    std::shared_ptr<Subcomponent> rParent;
    

    std::vector<Action> allActions;

    rQueue.push_back(root);

    while(!rQueue.empty()){
        
        rParent =rQueue.front();
        rQueue.pop_front();

        // for pseudoknots siblings will be divided into regions 
        std::vector<std::deque<std::shared_ptr<Subcomponent>> >sRegions;

        if(rParent->getType() == PseudoknotK ){
            for(int idx=1;idx<=9;idx +=2){
                auto tmp = rParent->getChildren(rParent->getDefiningBases()[idx],rParent->getDefiningBases()[idx+1]);
                std::deque<std::shared_ptr<Subcomponent> > q(tmp.begin(),tmp.end());
                sRegions.push_back(q);
            }
        }
        else if( rParent->getType() == PseudoknotH){
            for(int idx=1;idx<=5;idx +=2){
                auto tmp = rParent->getChildren(rParent->getDefiningBases()[idx],rParent->getDefiningBases()[idx+1]);
                std::deque<std::shared_ptr<Subcomponent> > q(tmp.begin(),tmp.end());
                sRegions.push_back(q);
            }
        }
        else{

            std::deque<std::shared_ptr<Subcomponent>> sQueue;
            for( auto c : rParent->getChildren()){
                sQueue.push_back(c);
            }
            if(rParent->getType()==Root)  std::reverse(sQueue.begin(),sQueue.end());
            sRegions.push_back(sQueue);
        }
        
        

        for( auto sQueue : sRegions)
        while(!sQueue.empty()){
            //Exterior Stem or Pseudoknot enumerate action for this SC
            current=sQueue.front();
            sQueue.pop_front();
            if(current->getType()==Exterior){
                // inside and pseudoneighbours
                // siblings (including descendants) exteriors 
                auto tmp=subcomponentActions(current,{sQueue.begin(), sQueue.end()});
                allActions.insert(allActions.end(), tmp.begin(), tmp.end());

            }else if(current->getType()==Stem){
                //all delete options
                auto tmp=subcomponentActions(current,std::vector<std::shared_ptr<Subcomponent>> { });
                allActions.insert(allActions.end(), tmp.begin(), tmp.end());
            }
            // NOTE Exterior loops are always evaluated before pseudo or stem 
            // psudoknots cant form base pairs with other siblings
            // so no need to look for base pairs to siblings
            else if(current->getType()==PseudoknotH){
               // inside and descendants
                auto tmp=subcomponentActions(current,std::vector<std::shared_ptr<Subcomponent>> { });
                allActions.insert(allActions.end(), tmp.begin(), tmp.end());
               // add to rQueue
               rQueue.push_back(current);

            }else if(current->getType()==PseudoknotK){
                // inside and descendants
                auto tmp=subcomponentActions(current,std::vector<std::shared_ptr<Subcomponent>> { });
                allActions.insert(allActions.end(), tmp.begin(), tmp.end());
                // add to rQueue
                rQueue.push_back(current);

            }

            //first subcomponent level untill hairpin loop or multiloop 
            //can not be pseudoknot 
            if(current->getType() == Stem){

                current = current->getChildren()[0];
                if(current->getType()==Multiloop){
                        auto tmp=subcomponentActions(current,{sQueue.begin(), sQueue.end()});
                        allActions.insert(allActions.end(), tmp.begin(), tmp.end());
                        rQueue.push_back(current);
                        
                    }

                    else{
                        auto tmp=subcomponentActions(current,{sQueue.begin(), sQueue.end()});
                        allActions.insert(allActions.end(), tmp.begin(), tmp.end());

                        while(current->getChildren().size()==1){
                        
                            current = current->getChildren()[0];

                            if(current->getType()==Multiloop){
                                auto tmp=subcomponentActions(current,std::vector<std::shared_ptr<Subcomponent>> { });
                                allActions.insert(allActions.end(), tmp.begin(), tmp.end());
                                rQueue.push_back(current);
                                break;
                            }

                            else{
                                auto tmp=subcomponentActions(current,std::vector<std::shared_ptr<Subcomponent>> { });
                                allActions.insert(allActions.end(), tmp.begin(), tmp.end());
                            } 


                                        
                                    
                        }
                    } 

                
            }
        }

    }
    
    return allActions;
}

void RNASecondaryStructure::saveSCActions(std::shared_ptr<Subcomponent> sc){


    //adding the inner action to dim1
    std::list<std::tuple<std::vector<Action>,double>> dim2;
    auto innerA = SCInnerActions(sc);
    double dim2RateSum=0;
    for( auto a : innerA){
        dim2RateSum +=std::exp(-a.getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) );
    }

    dim2.push_back(std::make_tuple(innerA,dim2RateSum));
    dim1.push_back(std::make_tuple(dim2,dim2RateSum));

    // adding sc to map3
    std::vector<std::tuple<std::shared_ptr<Subcomponent>,std::list<std::tuple<std::vector<Action>,double> >::iterator >> first;
    std::vector<std::shared_ptr<Subcomponent>> second{sc};
    std::list<std::tuple<std::list<std::tuple<std::vector<Action>,double>>,double>>::iterator third = dim1.end() ;
    std::advance(third,-1);
    map3.insert({ *sc , std::make_tuple(first,second,third) });
}


void RNASecondaryStructure::saveCrossingActions(std::shared_ptr<Subcomponent> sc, std::shared_ptr<Subcomponent> asc){

    

    //add crossing actions and rate sum (dim2 elements to dim1)
    auto crossA= SCCrossingActions(sc,asc);
    double dim2RateSum=0;
    for( auto a : crossA){
        dim2RateSum +=std::exp(-a.getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) );
    }
    std::get<0>(* std::get<2>(map3[*sc])).push_back(std::make_tuple(crossA,dim2RateSum));
    
    //updating map3
    //asc
    std::list<std::tuple<std::vector<Action>,double> >::iterator dim2ptr= std::get<0>(* std::get<2>(map3[*sc])).end();
    std::advance(dim2ptr,-1);
    std::get<0>(map3[*asc]).push_back(std::make_tuple(sc,dim2ptr));
    //sc
    std::get<1>(map3[*sc]).push_back(asc);

    // updating dim1 rate sum
    std::get<1>(* std::get<2>(map3[*sc])) += dim2RateSum;
}

void RNASecondaryStructure::inititEffAction(){



    //stack for "root" subcomponents i.e. Root, Multiloop or PseudoknotH/K
    std::vector<std::shared_ptr<Subcomponent>> rStack;
    //queue for "root" subcomponents i.e. Root, Multiloop or PseudoknotH/K
    std::deque<std::shared_ptr<Subcomponent>> rQueue;
    
    
    //current subcomponent
    std::shared_ptr<Subcomponent> current;
    //parent for subcomponents in sStack
    std::shared_ptr<Subcomponent> rParent;
    

    

    rQueue.push_back(root);

    while(!rQueue.empty()){
        
        rParent =rQueue.front();
        rQueue.pop_front();

        // for pseudoknots siblings will be divided into regions 
        std::vector<std::deque<std::shared_ptr<Subcomponent>> >sRegions;

        if(rParent->getType() == PseudoknotK ){
            for(int idx=1;idx<=9;idx +=2){
                auto tmp = rParent->getChildren(rParent->getDefiningBases()[idx],rParent->getDefiningBases()[idx+1]);
                std::deque<std::shared_ptr<Subcomponent> > q(tmp.begin(),tmp.end());
                sRegions.push_back(q);
            }
        }
        else if( rParent->getType() == PseudoknotH){
            for(int idx=1;idx<=5;idx +=2){
                auto tmp = rParent->getChildren(rParent->getDefiningBases()[idx],rParent->getDefiningBases()[idx+1]);
                std::deque<std::shared_ptr<Subcomponent> > q(tmp.begin(),tmp.end());
                sRegions.push_back(q);
            }
        }
        else{

            std::deque<std::shared_ptr<Subcomponent>> sQueue;
            for( auto c : rParent->getChildren()){
                sQueue.push_back(c);
            }
            if(rParent->getType()==Root)  std::reverse(sQueue.begin(),sQueue.end());
            sRegions.push_back(sQueue);
        }
        
        

        for( auto sQueue : sRegions)
        while(!sQueue.empty()){
            //Exterior Stem or Pseudoknot enumerate action for this SC
            current=sQueue.front();
            sQueue.pop_front();
            if(current->getType()==Exterior){
                // initial inner actions save
                saveSCActions(current);
                // save crossing actions
                for(auto asc : sQueue){
                    if(asc->getType() != Stem){
                        if(map3.find(*asc)==map3.end()){
                            saveSCActions(asc);
                        }
                        saveCrossingActions(current,asc);
                    }else{
                        auto tmpasc = asc->getChildren()[0];
                        if(map3.find(*tmpasc)==map3.end()){
                            saveSCActions(tmpasc);
                        }
                        saveCrossingActions(current,tmpasc);
                    }
                }



            }else if(current->getType()==Stem){
                saveSCActions(current);
            }
            // NOTE Exterior loops are always evaluated before pseudo or stem 
            // psudoknots cant form base pairs with other siblings
            // so no need to look for base pairs to siblings
            else if(current->getType()==PseudoknotH){
               // inside and descendants
                saveSCActions(current);
                for( auto asc : current->getChildren()){
                    if(asc->getType() != Stem){  
                        if(map3.find(*asc)==map3.end()){
                            saveSCActions(asc);
                        }
                        saveCrossingActions(current,asc);
                    }else{
                        auto tmpasc = asc->getChildren()[0];
                        if(map3.find(*tmpasc)==map3.end()){
                            saveSCActions(tmpasc);
                        }
                        saveCrossingActions(current,tmpasc);
                    }
                }

               // add to rQueue
               rQueue.push_back(current);

            }else if(current->getType()==PseudoknotK){
                // inside and descendants
                saveSCActions(current);
                for( auto asc : current->getChildren()){
                    if(asc->getType() != Stem){
                        if(map3.find(*asc)==map3.end()){
                            saveSCActions(asc);
                        }
                        saveCrossingActions(current,asc);
                    }else{
                        auto tmpasc = asc->getChildren()[0];
                        if(map3.find(*tmpasc)==map3.end()){
                            saveSCActions(tmpasc);
                        }
                        saveCrossingActions(current,tmpasc);
                    }
                }
                // add to rQueue
                rQueue.push_back(current);

            }

            //first subcomponent level untill hairpin loop or multiloop 
            //can not be pseudoknot 
            if(current->getType() == Stem){

                current = current->getChildren()[0];
                if(current->getType()==Multiloop){
                    // saving inside
                    saveSCActions(current);
                    // saving crossterm to children
                    for( auto asc : current->getChildren()){
                        if(asc->getType() != Stem){
                            if(map3.find(*asc)==map3.end()){
                                saveSCActions(asc);
                            }
                            saveCrossingActions(current,asc);
                        }else{
                            auto tmpasc = asc->getChildren()[0];
                            if(map3.find(*tmpasc)==map3.end()){
                                saveSCActions(tmpasc);
                            }
                            saveCrossingActions(current,tmpasc);
                        }
                    }
                    //saving crossterms to siblings
                    for(auto asc : sQueue){
                        if(asc->getType() != Stem){
                            if(map3.find(*asc)==map3.end()){
                                saveSCActions(asc);
                            }
                            saveCrossingActions(current,asc);
                        }else{
                            auto tmpasc = asc->getChildren()[0];
                            if(map3.find(*tmpasc)==map3.end()){
                                saveSCActions(tmpasc);
                            }
                            saveCrossingActions(current,tmpasc);
                        }
                    }

                    rQueue.push_back(current);
                        
                }

                else{

                    saveSCActions(current);
                    //grand child
                    auto tmpasc = current->getChildren()[0]->getChildren()[0];
                    if(map3.find(*tmpasc)==map3.end()){
                        saveSCActions(tmpasc);
                    }
                    saveCrossingActions(current,tmpasc);



                    while(current->getChildren().size()==1){
                    
                        current = current->getChildren()[0];

                        if(current->getType()==Multiloop){
                            saveSCActions(current);
                            // saving crossterm to children
                            for( auto asc : current->getChildren()){
                                if(asc->getType() != Stem){
                                    if(map3.find(*asc)==map3.end()){
                                        saveSCActions(asc);
                                    }
                                    saveCrossingActions(current,asc);
                                }else{
                                    auto tmpasc = asc->getChildren()[0];
                                    if(map3.find(*tmpasc)==map3.end()){
                                        saveSCActions(tmpasc);
                                    }
                                    saveCrossingActions(current,tmpasc);
                                }
                            }
                            //saving crossterms to siblings
                            for(auto asc : sQueue){
                                if(asc->getType() != Stem){
                                    if(map3.find(*asc)==map3.end()){
                                        saveSCActions(asc);
                                    }
                                    saveCrossingActions(current,asc);
                                }else{
                                    auto tmpasc = asc->getChildren()[0];
                                    if(map3.find(*tmpasc)==map3.end()){
                                        saveSCActions(tmpasc);
                                    }
                                    saveCrossingActions(current,tmpasc);
                                }
                            }
                            rQueue.push_back(current);
                            break;
                        }

                        else{
                            //interior bulge or hairpin
                            saveSCActions(current);
                            //grand child
                            if(current->getType()!= Hairpin){
                                auto tmpasc = current->getChildren()[0]->getChildren()[0];
                                if(map3.find(*tmpasc)==map3.end()){
                                    saveSCActions(tmpasc);
                                }
                                saveCrossingActions(current,tmpasc);
                            }
                        } 


                                    
                                
                    }
                } 

                
            }
        }

    }
    
   
}

const std::vector<Action>  RNASecondaryStructure::subcomponentActions(std::shared_ptr<Subcomponent> comp,const std::vector<std::shared_ptr<Subcomponent>> siblings){
//sibling and descendant pseudoknot action enum changed from while loop so that only pure H-type and K-type pseudoknots form
    
    std::vector<Action> actions;
    

    switch (comp->getType())
        {
        case Exterior:
            
            {
                double ene = 0;
                ene = -componentEnergy(root);
                //inside
                for( int i= comp->getLeft(); i<=comp->getRight(); i++){
                    for( int j=i+4; j<=comp->getRight(); j++){
                        if(primary.isPairable(i,j)){
                            ene = -componentEnergy(root);
                            auto ste =std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j}));
                            ene += componentEnergy(ste,i,j);
                            ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Hairpin,std::vector<int> {i+1,j-1})),i,j);
                            //ene += terminalEnergy(i,j);
                            auto tmpRoot = std::make_shared<Subcomponent>(Subcomponent(Root,std::vector<int> {-1,primary.getLength()}));
                            for( auto c : root->getChildren()){
                                tmpRoot->addChild(c);
                            }
                            tmpRoot->addChild(ste);
                            double moi =componentEnergy(tmpRoot,i,j);
                            ene += moi;
                            actions.push_back(Action(i,j,Addition,ene));
                            double ene2 =rebuidEnergy(i,j);
                            assert((float)ene==(float)ene2);
                        }
                    }
                }
                for ( auto sib : siblings){
                    
                    auto curr = sib;
                    if(curr->getType() == Exterior){
                        for (int  i= comp->getLeft(); i<=comp->getRight(); i++) {
                            for(int j= curr->getLeft(); j<= curr->getRight(); j++){
                                if(primary.isPairable(i,j)){
                                ene =rebuidEnergy(i,j);
                                actions.push_back(Action(i,j,Addition,ene));
                                }
                            }
                        }
                    }
                    else if(curr->getType() == PseudoknotH ){
                            //get inside pseudo
                        //exterior before pseudoknot
                        if(comp->getLeft()< curr->getLeft()){
                            for(int i=comp->getLeft();i<= comp->getRight(); i++){
                                auto childrenInside = curr->getChildren(curr->getDefiningBases()[1],curr->getDefiningBases()[2]);
                                for(int j= curr->getDefiningBases()[2]-1;j>curr->getDefiningBases()[1];j--){
                                    if( !childrenInside.empty() && j== childrenInside.back()->getRight()){
                                        j=childrenInside.back()->getLeft();
                                        childrenInside.pop_back();
                                    }else{
                                        if(primary.isPairable(i,j)){
                                            ene =rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }
                                }
                                assert(childrenInside.empty());
                            }
                            
                            if(comp->getRight() + 1 == curr->getLeft()){
                                int i = comp->getRight();
                                int j = curr->getDefiningBases()[5]+1;
                                if(bp[j]==-1 && primary.isPairable(i,j)){
                                    ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }
                            }
                        }
                        else{
                            for(int j=comp->getLeft();j<= comp->getRight(); j++){
                                auto childrenInside = curr->getChildren(curr->getDefiningBases()[5],curr->getDefiningBases()[6]);
                                for(int i= curr->getDefiningBases()[6]-1;i>curr->getDefiningBases()[5];i--){
                                    if( !childrenInside.empty() && i== childrenInside.back()->getRight()){
                                        i=childrenInside.back()->getLeft();
                                        childrenInside.pop_back();
                                    }else{
                                        if(primary.isPairable(i,j)){
                                            ene =rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }
                                }
                                assert(childrenInside.empty());
                            }

                            if(comp->getLeft() -1 == curr->getRight()){
                                int i= curr->getDefiningBases()[2] -1;
                                int j= comp->getLeft();
                                if(bp[i]==-1 && primary.isPairable(i,j)){
                                    ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }

                            }
                        }
                    }
                    else if(curr->getType() == PseudoknotK ){
                        //get inside pseudo
                        if(comp->getRight() + 1 == curr->getLeft()){
                                int i = comp->getRight();
                                int j = curr->getDefiningBases()[5]+1;
                                if(bp[j]==-1 && primary.isPairable(i,j)){
                                    ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }
                            }
                        if(comp->getLeft() -1 == curr->getRight()){
                                int i= curr->getDefiningBases()[6] -1;
                                int j= comp->getLeft();
                                if(bp[i]==-1 && primary.isPairable(i,j)){
                                    ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }

                            }
                    }
                    else{
                        /* while(curr->getType() != Multiloop && !curr->getChildren().empty() ){
                            curr=curr->getChildren()[0];
                            if(curr->getType() != Stem){
                                auto tmp =makePseudo(comp,curr,siblings,sib,nullptr);
                                actions.insert(actions.end(),tmp.begin(),tmp.end());
                            }
                            
                        } */

                        
                        // Skipping first stem
                        curr=curr->getChildren()[0];
                        // makeing pseudo with next loop
                        auto tmp =makePseudo(comp,curr);
                        actions.insert(actions.end(),tmp.begin(),tmp.end());
                    }
                    
                }

                
            }
                

        break;
        case Stem:
        {
            //purely inside
            for( int i=comp->getLeft()+1; i< comp->getDefiningBases()[1]; i++){
                double ene=0;
                ene -= componentEnergy(comp);
                ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {comp->getLeft(), i-1, bp[i] +1, comp->getRight()})),i,bp[i]);
                ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> { i+1, comp->getDefiningBases()[1], comp->getDefiningBases()[2], bp[i]-1})),i,bp[i]);
                ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Interior,std::vector<int> {i,i,bp[i],bp[i]})),i,bp[i]);
                actions.push_back(Action(i,bp[i],Deletion,ene));

            }
            //border cases (~optimizable)
            // first pair (and last if lone)
            {
                 
                int i=comp->getLeft();
                int j= bp[i];
                double ene = rebuidEnergy(i,j);
                actions.push_back(Action(i,bp[i],Deletion,ene));

            }
            //last pair
            if(comp->getLeft() != comp->getDefiningBases()[1])
            {
                int i=comp->getDefiningBases()[1];
                int j= bp[i];
                double ene = rebuidEnergy(i,j);
                actions.push_back(Action(i,bp[i],Deletion,ene));
            }



        }
            
        break;
        case Interior:
        {

            //purely inside
            //upper side
            for( int i= comp->getLeft(); i<=comp->getDefiningBases()[1]; i++){
                for( int j=i+4; j<=comp->getDefiningBases()[1]; j++){
                    if(primary.isPairable(i,j)){
                        double ene = -componentEnergy(comp);
                        auto s =std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j}));
                        ene += componentEnergy(s,i,j);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Hairpin,std::vector<int> {i+1,j-1})),i,j);
                        std::vector<int> mvec;
                        for( int idx=comp->getLeft(); idx<i; idx++) mvec.push_back(idx);
                        for( int idx=j+1; idx<=comp->getDefiningBases()[1]; idx++) mvec.push_back(idx);
                        for( int idx=comp->getDefiningBases()[2]; idx<=comp->getRight(); idx++) mvec.push_back(idx);
                        
                        auto multi = std::make_shared<Subcomponent>(Subcomponent(Multiloop,mvec)); 
                        multi->addChild(s);
                        multi->addChild(comp->getChildren()[0]);
                        ene += componentEnergy(multi,i,j);
                        actions.push_back(Action(i,j,Addition,ene));

                    }
                }
            }
            //bottom side
            for( int i= comp->getDefiningBases()[2]; i<=comp->getRight(); i++){
                for( int j=i+4; j<=comp->getRight(); j++){
                    if(primary.isPairable(i,j)){
                        double ene = -componentEnergy(comp);
                        auto s =std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j}));
                        ene += componentEnergy(s,i,j);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Hairpin,std::vector<int> {i+1,j-1})),i,j);
                        std::vector<int> mvec;
                        for( int idx=comp->getLeft(); idx<=comp->getDefiningBases()[1]; idx++) mvec.push_back(idx);
                        for( int idx=comp->getDefiningBases()[2]; idx<i; idx++) mvec.push_back(idx);
                        for( int idx=j+1; idx<=comp->getRight(); idx++) mvec.push_back(idx);
                        
                        
                        auto multi = std::make_shared<Subcomponent>(Subcomponent(Multiloop,mvec)); 
                        multi->addChild(comp->getChildren()[0]);
                        multi->addChild(s);
                        ene += componentEnergy(multi,i,j);
                        actions.push_back(Action(i,j,Addition,ene));

                    }
                }
            }
            //both sides 
            for( int i=comp->getLeft()+1; i< comp->getDefiningBases()[1]; i++){
                for( int j=comp->getDefiningBases()[2]+1; j <comp->getRight(); j++){
                    if(primary.isPairable(i,j)){
                        double ene=0;
                        ene -= componentEnergy(comp);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Interior,std::vector<int> {comp->getLeft(), i-1, j +1, comp->getRight()})),i,j);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Interior,std::vector<int> { i+1, comp->getDefiningBases()[1], comp->getDefiningBases()[2], j-1})),i,j);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j})),i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                    }
                }
                

            }

            // border cases (optimizable) 
            {
            int j = comp->getDefiningBases()[2];
            for( int i=comp->getLeft(); i<= comp->getDefiningBases()[1]; i++){
                if(primary.isPairable(i,j)){
                double ene =  rebuidEnergy(i,j);
                actions.push_back(Action(i,j,Addition,ene));
                }
            }
            j=comp->getRight();
            for( int i=comp->getLeft(); i<= comp->getDefiningBases()[1]; i++){
                if(primary.isPairable(i,j)){
                double ene =  rebuidEnergy(i,j);
                actions.push_back(Action(i,j,Addition,ene));
                }
            }
            }
            {
            int i = comp->getLeft();
            for( int j=comp->getDefiningBases()[2] +1; j< comp->getRight(); j++){
                if(primary.isPairable(i,j)){
                double ene =  rebuidEnergy(i,j);
                actions.push_back(Action(i,j,Addition,ene));
                }
            }
            i=comp->getDefiningBases()[1];
             for( int j=comp->getDefiningBases()[2] +1; j< comp->getRight(); j++){
                if(primary.isPairable(i,j)){
                double ene = rebuidEnergy(i,j);
                actions.push_back(Action(i,j,Addition,ene));
                }
            }
            }
            
        
            

            //PSEUDOKNOTS
            //descendants
            auto curr = comp->getChildren()[0];

            /* while(curr->getType() != Multiloop && !curr->getChildren().empty() ){
                curr=curr->getChildren()[0];
                if(curr->getType() != Stem){
                        auto tmp =makePseudo(comp,curr,siblings,nullptr,nullptr);
                    actions.insert(actions.end(),tmp.begin(),tmp.end());
                }
                            
            } */

            // Skipping first stem
            curr=curr->getChildren()[0];
            // makeing pseudo with next loop
            auto tmp =makePseudo(comp,curr);
            actions.insert(actions.end(),tmp.begin(),tmp.end());

            //siblings
            for( auto sib : siblings){
                auto curr = sib;
                if(curr->getType() == PseudoknotK ){
                    //do nothing
                    }
                    else if(curr->getType() == PseudoknotH ){
                        //do nothing
                    }
                    else{
                        /* while(curr->getType() != Multiloop && !curr->getChildren().empty() ){
                            curr=curr->getChildren()[0];
                            if(curr->getType() != Stem){
                                 auto tmp =makePseudo(comp,curr,siblings,sib,nullptr);
                                actions.insert(actions.end(),tmp.begin(),tmp.end());
                            }
                            
                        } */

                        // Skipping first stem
                        curr=curr->getChildren()[0];
                        // makeing pseudo with next loop
                        auto tmp =makePseudo(comp,curr);
                        actions.insert(actions.end(),tmp.begin(),tmp.end());
                    }
            }
        }
        break;
        case Bulge:
        {
            for( int i= comp->getLeft(); i<=comp->getRight(); i++){
                for( int j=i+4; j<=comp->getRight(); j++){
                    if(primary.isPairable(i,j)){
                        double ene = 0;
                        auto s =std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j}));
                        ene += componentEnergy(s,i,j);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Hairpin,std::vector<int> {i+1,j-1})),i,j);
                        std::vector<int> mvec;
                        
                        for( int idx=comp->getLeft(); idx<i; idx++) mvec.push_back(idx);
                        for( int idx=j+1; idx<=comp->getRight(); idx++) mvec.push_back(idx);
                        
                        
                        auto multi = std::make_shared<Subcomponent>(Subcomponent(Multiloop,std::vector<int> {i+1,j-1}));
                        multi->addChild(comp->getChildren()[0]);
                        multi->addChild(s);
                        ene += componentEnergy(multi,i,j);
                        actions.push_back(Action(i,j,Addition,ene));

                    }
                }
            }

            //PSEUDOKNOTS
            //descendants
            auto curr = comp->getChildren()[0];
            /* while(curr->getType() != Multiloop && !curr->getChildren().empty() ){
                curr=curr->getChildren()[0];
                if(curr->getType() != Stem){
                        auto tmp =makePseudo(comp,curr,siblings,nullptr,nullptr);
                    actions.insert(actions.end(),tmp.begin(),tmp.end());
                }
                            
            } */

            // Skipping first stem
            curr=curr->getChildren()[0];
            // makeing pseudo with next loop
            auto tmp =makePseudo(comp,curr);
            actions.insert(actions.end(),tmp.begin(),tmp.end());

            //siblings
            for( auto sib : siblings){
                auto curr = sib;
                if(curr->getType() == PseudoknotK ){
                        //do nothing
                    }
                    else if(curr->getType() == PseudoknotH ){
                        //do nothing
                    }
                    else{
                        /* while(curr->getType() != Multiloop && !curr->getChildren().empty() ){
                            curr=curr->getChildren()[0];
                            if(curr->getType() != Stem){
                                 auto tmp =makePseudo(comp,curr,siblings,sib,nullptr);
                                actions.insert(actions.end(),tmp.begin(),tmp.end());
                            }
                            
                        } */

                        // Skipping first stem
                        curr=curr->getChildren()[0];
                        // makeing pseudo with next loop
                        auto tmp =makePseudo(comp,curr);
                        actions.insert(actions.end(),tmp.begin(),tmp.end());
                    }
            }
        }
        break;
        case Hairpin:
        {
            //pure inside 
        
         for( int i= comp->getLeft()+1; i<comp->getRight(); i++){
            for( int j=i+4; j<comp->getRight(); j++){
                if(primary.isPairable(i,j)){
                    double ene=0;
                    ene -= componentEnergy(comp);
                    ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Interior,std::vector<int> {comp->getLeft(), i-1, j +1, comp->getRight()})),i,j);
                    ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Hairpin,std::vector<int> { i+1,  j-1})),i,j);
                    ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j})),i,j);
                    
                    actions.push_back(Action(i,j,Addition,ene));
                }

            }

            
         }   

         //border cases (optimizable)
        {
            int i= comp->getLeft();
            for( int j=i+4; j<=comp->getRight(); j++){
                if(primary.isPairable(i,j)){
                    double ene=  rebuidEnergy(i,j);
                    actions.push_back(Action(i,j,Addition,ene));

                }
            }
            int j = comp->getRight();
            for( int i= comp->getLeft() +1; i<=j-4; i++){
                if(primary.isPairable(i,j)){
                    double ene= rebuidEnergy(i,j);
                    actions.push_back(Action(i,j,Addition,ene));

                }
            }
        }

        //PSEUDOKNOTS
            //siblings
            for( auto sib : siblings){
                auto curr = sib;
                if(curr->getType() == PseudoknotK ){
                        //do nothing
                    }
                    else if(curr->getType() == PseudoknotH ){
                        //do nothing
                    }
                    else{
                        /* while(curr->getType() != Multiloop && !curr->getChildren().empty() ){
                            curr=curr->getChildren()[0];
                            if(curr->getType() != Stem){
                                 auto tmp =makePseudo(comp,curr,siblings,sib,nullptr);
                                actions.insert(actions.end(),tmp.begin(),tmp.end());
                            }
                            
                        } */

                        // Skipping first stem
                        curr=curr->getChildren()[0];
                        // makeing pseudo with next loop
                        auto tmp =makePseudo(comp,curr);
                        actions.insert(actions.end(),tmp.begin(),tmp.end());
                    }
            }

        }
        break;
        case Multiloop:
        {
            // (~~optimizable)
            for( int i : comp->getDefiningBases()){
                for( int j : comp->getDefiningBases()){
                    if( i<j && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));

                    }
                }
            }


            //PSEUDOKNOTS
            //descendants
            for(auto curr : comp->getChildren()){
                /* while(curr->getType() != Multiloop && !curr->getChildren().empty() ){
                            curr=curr->getChildren()[0];
                            if(curr->getType() != Stem){
                                 auto tmp =makePseudo(comp,curr,siblings,nullptr,nullptr);
                                actions.insert(actions.end(),tmp.begin(),tmp.end());
                            }
                            
                        } */
                        // Skipping first stem
                        if(curr->getType() == Stem){
                            curr=curr->getChildren()[0];
                            // makeing pseudo with next loop
                            auto tmp =makePseudo(comp,curr);
                            actions.insert(actions.end(),tmp.begin(),tmp.end());

                        } else if( curr->getType() == PseudoknotH){
                            //from left side of multiloop to 1-2 making K-type
                            int idx=0;
                            while(idx <comp->getDefiningBases().size() && comp->getDefiningBases()[idx] < curr->getLeft()){
                                int i =comp->getDefiningBases()[idx];
                                auto childrenInside = curr->getChildren(curr->getDefiningBases()[1],curr->getDefiningBases()[2]);
                                for(int j= curr->getDefiningBases()[2]-1;j>curr->getDefiningBases()[1];j--){
                                    if( !childrenInside.empty() && j== childrenInside.back()->getRight()){
                                        j=childrenInside.back()->getLeft();
                                        childrenInside.pop_back();
                                    }else{
                                        if(primary.isPairable(i,j)){
                                            double ene =rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }
                                }
                                
                                assert(childrenInside.empty());
                                idx++;
                            }
                            //adding outer base pair to stem 
                            if(idx > 0 && comp->getDefiningBases()[idx-1] ==curr->getLeft()-1){ 
                                int i = comp->getDefiningBases()[idx-1];
                                int j = curr->getDefiningBases()[5]+1;
                                if(bp[j]==-1 && primary.isPairable(i,j)){
                                    double ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }
                            }

                            //adding outer base pair to stem 
                            if(idx <comp->getDefiningBases().size() && comp->getDefiningBases()[idx] == curr->getRight()+1){ 
                                int i= curr->getDefiningBases()[2] -1;
                                int j= comp->getDefiningBases()[idx];
                                if(bp[i]==-1 && primary.isPairable(i,j)){
                                    double ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }
                            }


                            //from Right side of multiloop to 5-6 making K-type
                            while(idx <comp->getDefiningBases().size()){
                                int j =comp->getDefiningBases()[idx];
                                auto childrenInside = curr->getChildren(curr->getDefiningBases()[5],curr->getDefiningBases()[6]);
                                for(int i= curr->getDefiningBases()[6]-1;i>curr->getDefiningBases()[5];i--){
                                    if( !childrenInside.empty() && i== childrenInside.back()->getRight()){
                                        i=childrenInside.back()->getLeft();
                                        childrenInside.pop_back();
                                    }else{
                                        if(primary.isPairable(i,j)){
                                            double ene =rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }
                                }
                                assert(childrenInside.empty());

                                idx++;
                            }
                            




                        }else if( curr->getType() == PseudoknotK){
                            
                            int idx=0;
                            while(idx <comp->getDefiningBases().size() && comp->getDefiningBases()[idx] < curr->getLeft()){
                                idx++;
                            }
                            //adding outer base pair to stem 
                            if(idx > 0 && comp->getDefiningBases()[idx-1] ==curr->getLeft()-1){ 
                                int i = comp->getDefiningBases()[idx-1];
                                int j = curr->getDefiningBases()[5]+1;
                                if(bp[j]==-1 && primary.isPairable(i,j)){
                                    double ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }
                            }

                            //adding outer base pair to stem 
                            if(idx <comp->getDefiningBases().size() && comp->getDefiningBases()[idx] == curr->getRight()+1){ 
                                int i= curr->getDefiningBases()[6] -1;
                                int j= comp->getDefiningBases()[idx];
                                if(bp[i]==-1 && primary.isPairable(i,j)){
                                    double ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }
                            }
                        }
            }
            
            //siblings
            for( auto sib : siblings){
                auto curr = sib;
                if(curr->getType() == PseudoknotK ){
                        //do nothing
                    }
                    else if(curr->getType() == PseudoknotH ){
                        //do nothing
                    }
                    else{
                        /* while(curr->getType() != Multiloop && !curr->getChildren().empty() ){
                            curr=curr->getChildren()[0];
                            if(curr->getType() != Stem){
                                 auto tmp =makePseudo(comp,curr,siblings,sib,nullptr);
                                actions.insert(actions.end(),tmp.begin(),tmp.end());
                            }
                            
                        } */
                        // Skipping first stem
                        curr=curr->getChildren()[0];
                        // makeing pseudo with next loop
                        auto tmp =makePseudo(comp,curr);
                        actions.insert(actions.end(),tmp.begin(),tmp.end());
                    }
            }
        }
        break;
        case PseudoknotH:
        {
            // all iside pseudo
            //inside 1-2, 3-4, 5-6
            for(int idx=1;idx<=5;idx +=2){
                auto childrenInside = comp->getChildren(comp->getDefiningBases()[idx],comp->getDefiningBases()[idx+1]);
                for(int j= comp->getDefiningBases()[idx+1]-1;j>comp->getDefiningBases()[idx];j--){
                    if( !childrenInside.empty() && j== childrenInside.back()->getRight()){
                        j=childrenInside.back()->getLeft();
                        childrenInside.pop_back();
                    }else{
                        auto childrenInside2=childrenInside;
                        for( int i = j;i>comp->getDefiningBases()[idx];i--){
                            

                        
                            if(!childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                i=childrenInside2.back()->getLeft();
                                childrenInside2.pop_back();
                            }else if(primary.isPairable(i,j)){
                                double ene=  rebuidEnergy(i,j);
                                actions.push_back(Action(i,j,Addition,ene));
                            }
                        }
                        assert(childrenInside2.empty());
                        
                    }
                }
                assert(childrenInside.empty());
            }
            // 1-2,3-4 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[1]+1;
                int j = comp->getDefiningBases()[4]-1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }
            //3-4,5-6 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[3]+1;
                int j = comp->getDefiningBases()[6]-1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }

            //deletions (ONLY outermost/innermost of bands)
            for(int idx=0; idx<4; idx++){
                double ene = rebuidEnergy(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]]);
                actions.push_back(Action(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]],Deletion,ene));
            }

            // make pseudo ( of children) 
            for(int idx=1;idx<=5;idx +=2){
                auto childrenInside = comp->getChildren(comp->getDefiningBases()[idx],comp->getDefiningBases()[idx+1]);
                if(!childrenInside.empty())
                for(int j= comp->getDefiningBases()[idx+1]-1;j>comp->getDefiningBases()[idx];j--){
                    
                    if( !childrenInside.empty() && j== childrenInside.back()->getRight()){
                        j=childrenInside.back()->getLeft();
                        childrenInside.pop_back();
                    }else{
                        for( auto c : comp->getChildren(comp->getDefiningBases()[idx],comp->getDefiningBases()[idx+1])){
                            if(c->getType() == Stem){
                                auto curr=c->getChildren()[0];
                                // makeing pseudo with next loop
                                if(curr->getType() == Multiloop){
                                    for(int i : curr->getDefiningBases()){
                                        if(primary.isPairable(i,j)){
                                            double ene=  rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }

                                    }

                                }else if(curr->getType() == Interior){
                                    for(int i= curr->getLeft(); i<=curr->getDefiningBases()[1]; i++){
                                        if(primary.isPairable(i,j)){
                                            double ene=  rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }
                                    for(int i=curr->getDefiningBases()[2]; i<= curr->getRight(); i++){
                                        if(primary.isPairable(i,j)){
                                            double ene=  rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }

                                }else{
                                    for(int i= curr->getLeft(); i<=curr->getRight(); i++){
                                        if(primary.isPairable(i,j)){
                                            double ene=  rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }

                                }

                            }else if(c->getType() == PseudoknotH){
                                if(j< c->getLeft()){
                                    auto childrenInside2 = c->getChildren(c->getDefiningBases()[1],c->getDefiningBases()[2]);
                                    for(int i= c->getDefiningBases()[2]-1;i>c->getDefiningBases()[1];i--){
                                        if( !childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                            i=childrenInside2.back()->getLeft();
                                            childrenInside2.pop_back();
                                        }else{
                                            if(primary.isPairable(i,j)){
                                                double ene =rebuidEnergy(j,i);
                                                actions.push_back(Action(j,i,Addition,ene));
                                            }
                                        }
                                    }
                                
                                    assert(childrenInside2.empty());
                                    if(j==c->getLeft()-1){
                                        int i = c->getDefiningBases()[5]+1;
                                        if(bp[i]==-1 && primary.isPairable(i,j)){
                                            double ene =rebuidEnergy(j,i);
                                            actions.push_back(Action(j,i,Addition,ene));
                                        }
                                    }
                                }else{
                                    auto childrenInside2 = c->getChildren(c->getDefiningBases()[5],c->getDefiningBases()[6]);
                                    for(int i= c->getDefiningBases()[6]-1;i>c->getDefiningBases()[5];i--){
                                        if( !childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                            i=childrenInside2.back()->getLeft();
                                            childrenInside2.pop_back();
                                        }else{
                                            if(primary.isPairable(i,j)){
                                                double ene =rebuidEnergy(i,j);
                                                actions.push_back(Action(i,j,Addition,ene));
                                            }
                                        }
                                    }
                                    assert(childrenInside2.empty());
                                    if(j==c->getRight()+1){
                                        int i = c->getDefiningBases()[2]-1;
                                        if(bp[i]==-1 && primary.isPairable(i,j)){
                                            double ene =rebuidEnergy(j,i);
                                            actions.push_back(Action(j,i,Addition,ene));
                                        }
                                    }
                                }

                            }else if(c->getType() == PseudoknotK){
                                if(j==c->getLeft()-1){
                                    int i = c->getDefiningBases()[5]+1;
                                    if(bp[i]==-1 && primary.isPairable(i,j)){
                                        double ene =rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }else if(j==c->getRight()+1){
                                    int i = c->getDefiningBases()[6]-1;
                                    if(bp[i]==-1 && primary.isPairable(i,j)){
                                        double ene =rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }
                                
                                
                                
                            }
                        }
                        
                        
                    }
                }
                assert(childrenInside.empty());
            }
            
        }
        break;
        case PseudoknotK:
        {
            // all iside pseudo
            //inside 1-2, 3-4, 5-6, 7-8 ,9-10
            for(int idx=1;idx<=9;idx +=2){
                auto childrenInside = comp->getChildren(comp->getDefiningBases()[idx],comp->getDefiningBases()[idx+1]);
                for(int j= comp->getDefiningBases()[idx+1]-1;j>comp->getDefiningBases()[idx];j--){
                    if(!childrenInside.empty() && j== childrenInside.back()->getRight()){
                        j=childrenInside.back()->getLeft();
                        childrenInside.pop_back();
                    }else{
                        auto childrenInside2=childrenInside;
                        for( int i = j;i>comp->getDefiningBases()[idx];i--){
                            

                        
                            if(!childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                i=childrenInside2.back()->getLeft();
                                childrenInside2.pop_back();
                            }else if(primary.isPairable(i,j)){
                                double ene=  rebuidEnergy(i,j);
                                actions.push_back(Action(i,j,Addition,ene));
                            }
                        
                        }
                        assert(childrenInside2.size() == 0);
                    
                    }
                }
                assert(childrenInside.size() == 0);
            }

            //deletions (ONLY outermost/innermost of bands)
            for(int idx=0; idx<4; idx++){
                double ene = rebuidEnergy(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]]);
                actions.push_back(Action(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]],Deletion,ene));
            }for(int idx=6; idx<8; idx++){
                double ene = rebuidEnergy(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]]);
                actions.push_back(Action(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]],Deletion,ene));
            }



            // 1-2 ,3-4 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[1]+1;
                int j = comp->getDefiningBases()[4]-1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }
            //3-4,7-8 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[3]+1;
                int j = comp->getDefiningBases()[8]-1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }
            //7-8,9-10 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[7]+1;
                int j = comp->getDefiningBases()[10]-1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }
            //1-2,9-10 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[2]-1;
                int j = comp->getDefiningBases()[9]+1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }

            // make pseudo ( of children) 
            for(int idx=1;idx<=9;idx +=2){
                auto childrenInside = comp->getChildren(comp->getDefiningBases()[idx],comp->getDefiningBases()[idx+1]);
                if(!childrenInside.empty())
                for(int j= comp->getDefiningBases()[idx+1]-1;j>comp->getDefiningBases()[idx];j--){
                    
                    if( !childrenInside.empty() && j== childrenInside.back()->getRight()){
                        j=childrenInside.back()->getLeft();
                        childrenInside.pop_back();
                    }else{
                        for( auto c : comp->getChildren(comp->getDefiningBases()[idx],comp->getDefiningBases()[idx+1])){
                            if(c->getType() == Stem){
                                auto curr=c->getChildren()[0];
                                // makeing pseudo with next loop
                                if(curr->getType() == Multiloop){
                                    for(int i : curr->getDefiningBases()){
                                        if(primary.isPairable(i,j)){
                                            double ene=  rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }

                                    }

                                }else if(curr->getType() == Interior){
                                    for(int i= curr->getLeft(); i<=curr->getDefiningBases()[1]; i++){
                                        if(primary.isPairable(i,j)){
                                            double ene=  rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }
                                    for(int i=curr->getDefiningBases()[2]; i<= curr->getRight(); i++){
                                        if(primary.isPairable(i,j)){
                                            double ene=  rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }

                                }else{
                                    for(int i= curr->getLeft(); i<=curr->getRight(); i++){
                                        if(primary.isPairable(i,j)){
                                            double ene=  rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }

                                }

                            }else if(c->getType() == PseudoknotH){
                                if(j< c->getLeft()){
                                    auto childrenInside2 = c->getChildren(c->getDefiningBases()[1],c->getDefiningBases()[2]);
                                    for(int i= c->getDefiningBases()[2]-1;i>c->getDefiningBases()[1];i--){
                                        if( !childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                            i=childrenInside2.back()->getLeft();
                                            childrenInside2.pop_back();
                                        }else{
                                            if(primary.isPairable(i,j)){
                                                double ene =rebuidEnergy(j,i);
                                                actions.push_back(Action(j,i,Addition,ene));
                                            }
                                        }
                                    }
                                
                                    assert(childrenInside2.empty());
                                    if(j==c->getLeft()-1){
                                        int i = c->getDefiningBases()[5]+1;
                                        if(bp[i]==-1 && primary.isPairable(i,j)){
                                            double ene =rebuidEnergy(j,i);
                                            actions.push_back(Action(j,i,Addition,ene));
                                        }
                                    }
                                }else{
                                    auto childrenInside2 = c->getChildren(c->getDefiningBases()[5],c->getDefiningBases()[6]);
                                    for(int i= c->getDefiningBases()[6]-1;i>c->getDefiningBases()[5];i--){
                                        if( !childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                            i=childrenInside2.back()->getLeft();
                                            childrenInside2.pop_back();
                                        }else{
                                            if(primary.isPairable(i,j)){
                                                double ene =rebuidEnergy(i,j);
                                                actions.push_back(Action(i,j,Addition,ene));
                                            }
                                        }
                                    }
                                    assert(childrenInside2.empty());
                                    if(j==c->getRight()+1){
                                        int i = c->getDefiningBases()[2]-1;
                                        if(bp[i]==-1 && primary.isPairable(i,j)){
                                            double ene =rebuidEnergy(j,i);
                                            actions.push_back(Action(j,i,Addition,ene));
                                        }
                                    }
                                }

                            }else if(c->getType() == PseudoknotK){
                                if(j==c->getLeft()-1){
                                    int i = c->getDefiningBases()[5]+1;
                                    if(bp[i]==-1 && primary.isPairable(i,j)){
                                        double ene =rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }else if(j==c->getRight()+1){
                                    int i = c->getDefiningBases()[6]-1;
                                    if(bp[i]==-1 && primary.isPairable(i,j)){
                                        double ene =rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }
                                
                                
                                
                            }
                        }
                        
                        
                    }
                }
                assert(childrenInside.empty());
            }

        }
        break;
        
        

    }
    return actions;

}

const std::vector<Action> RNASecondaryStructure::SCInnerActions(std::shared_ptr<Subcomponent> comp){
    
    std::vector<Action> actions;
    

    switch (comp->getType())
        {
        case Exterior:
            
            {
                double ene = 0;
                ene = -componentEnergy(root);
                //inside
                for( int i= comp->getLeft(); i<=comp->getRight(); i++){
                    for( int j=i+4; j<=comp->getRight(); j++){
                        if(primary.isPairable(i,j)){
                            ene = -componentEnergy(root);
                            auto ste =std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j}));
                            ene += componentEnergy(ste,i,j);
                            ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Hairpin,std::vector<int> {i+1,j-1})),i,j);
                            //ene += terminalEnergy(i,j);
                            auto tmpRoot = std::make_shared<Subcomponent>(Subcomponent(Root,std::vector<int> {-1,primary.getLength()}));
                            for( auto c : root->getChildren()){
                                tmpRoot->addChild(c);
                            }
                            tmpRoot->addChild(ste);
                            double moi =componentEnergy(tmpRoot,i,j);
                            ene += moi;
                            actions.push_back(Action(i,j,Addition,ene));
                            double ene2 =rebuidEnergy(i,j);
                            assert((float)ene==(float)ene2);
                        }
                    }
                }
                

                
            }
                

        break;
        case Stem:
        {
            //purely inside
            for( int i=comp->getLeft()+1; i< comp->getDefiningBases()[1]; i++){
                double ene=0;
                ene -= componentEnergy(comp);
                ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {comp->getLeft(), i-1, bp[i] +1, comp->getRight()})),i,bp[i]);
                ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> { i+1, comp->getDefiningBases()[1], comp->getDefiningBases()[2], bp[i]-1})),i,bp[i]);
                ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Interior,std::vector<int> {i,i,bp[i],bp[i]})),i,bp[i]);
                actions.push_back(Action(i,bp[i],Deletion,ene));

            }
            //border cases (~optimizable)
            // first pair (and last if lone)
            {
                 
                int i=comp->getLeft();
                int j= bp[i];
                double ene = rebuidEnergy(i,j);
                actions.push_back(Action(i,bp[i],Deletion,ene));

            }
            //last pair
            if(comp->getLeft() != comp->getDefiningBases()[1])
            {
                int i=comp->getDefiningBases()[1];
                int j= bp[i];
                double ene = rebuidEnergy(i,j);
                actions.push_back(Action(i,bp[i],Deletion,ene));
            }



        }
            
        break;
        case Interior:
        {

            //purely inside
            //upper side
            for( int i= comp->getLeft(); i<=comp->getDefiningBases()[1]; i++){
                for( int j=i+4; j<=comp->getDefiningBases()[1]; j++){
                    if(primary.isPairable(i,j)){
                        double ene = -componentEnergy(comp);
                        auto s =std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j}));
                        ene += componentEnergy(s,i,j);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Hairpin,std::vector<int> {i+1,j-1})),i,j);
                        std::vector<int> mvec;
                        for( int idx=comp->getLeft(); idx<i; idx++) mvec.push_back(idx);
                        for( int idx=j+1; idx<=comp->getDefiningBases()[1]; idx++) mvec.push_back(idx);
                        for( int idx=comp->getDefiningBases()[2]; idx<=comp->getRight(); idx++) mvec.push_back(idx);
                        
                        auto multi = std::make_shared<Subcomponent>(Subcomponent(Multiloop,mvec)); 
                        multi->addChild(s);
                        multi->addChild(comp->getChildren()[0]);
                        ene += componentEnergy(multi,i,j);
                        actions.push_back(Action(i,j,Addition,ene));

                    }
                }
            }
            //bottom side
            for( int i= comp->getDefiningBases()[2]; i<=comp->getRight(); i++){
                for( int j=i+4; j<=comp->getRight(); j++){
                    if(primary.isPairable(i,j)){
                        double ene = -componentEnergy(comp);
                        auto s =std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j}));
                        ene += componentEnergy(s,i,j);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Hairpin,std::vector<int> {i+1,j-1})),i,j);
                        std::vector<int> mvec;
                        for( int idx=comp->getLeft(); idx<=comp->getDefiningBases()[1]; idx++) mvec.push_back(idx);
                        for( int idx=comp->getDefiningBases()[2]; idx<i; idx++) mvec.push_back(idx);
                        for( int idx=j+1; idx<=comp->getRight(); idx++) mvec.push_back(idx);
                        
                        
                        auto multi = std::make_shared<Subcomponent>(Subcomponent(Multiloop,mvec)); 
                        multi->addChild(comp->getChildren()[0]);
                        multi->addChild(s);
                        ene += componentEnergy(multi,i,j);
                        actions.push_back(Action(i,j,Addition,ene));

                    }
                }
            }
            //both sides 
            for( int i=comp->getLeft()+1; i< comp->getDefiningBases()[1]; i++){
                for( int j=comp->getDefiningBases()[2]+1; j <comp->getRight(); j++){
                    if(primary.isPairable(i,j)){
                        double ene=0;
                        ene -= componentEnergy(comp);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Interior,std::vector<int> {comp->getLeft(), i-1, j +1, comp->getRight()})),i,j);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Interior,std::vector<int> { i+1, comp->getDefiningBases()[1], comp->getDefiningBases()[2], j-1})),i,j);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j})),i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                    }
                }
                

            }

            // border cases (optimizable) 
            {
            int j = comp->getDefiningBases()[2];
            for( int i=comp->getLeft(); i<= comp->getDefiningBases()[1]; i++){
                if(primary.isPairable(i,j)){
                double ene =  rebuidEnergy(i,j);
                actions.push_back(Action(i,j,Addition,ene));
                }
            }
            j=comp->getRight();
            for( int i=comp->getLeft(); i<= comp->getDefiningBases()[1]; i++){
                if(primary.isPairable(i,j)){
                double ene =  rebuidEnergy(i,j);
                actions.push_back(Action(i,j,Addition,ene));
                }
            }
            }
            {
            int i = comp->getLeft();
            for( int j=comp->getDefiningBases()[2] +1; j< comp->getRight(); j++){
                if(primary.isPairable(i,j)){
                double ene =  rebuidEnergy(i,j);
                actions.push_back(Action(i,j,Addition,ene));
                }
            }
            i=comp->getDefiningBases()[1];
             for( int j=comp->getDefiningBases()[2] +1; j< comp->getRight(); j++){
                if(primary.isPairable(i,j)){
                double ene = rebuidEnergy(i,j);
                actions.push_back(Action(i,j,Addition,ene));
                }
            }
            }
            
        
        
        }
        break;
        case Bulge:
        {
            for( int i= comp->getLeft(); i<=comp->getRight(); i++){
                for( int j=i+4; j<=comp->getRight(); j++){
                    if(primary.isPairable(i,j)){
                        double ene = 0;
                        auto s =std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j}));
                        ene += componentEnergy(s,i,j);
                        ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Hairpin,std::vector<int> {i+1,j-1})),i,j);
                        std::vector<int> mvec;
                        
                        for( int idx=comp->getLeft(); idx<i; idx++) mvec.push_back(idx);
                        for( int idx=j+1; idx<=comp->getRight(); idx++) mvec.push_back(idx);
                        
                        
                        auto multi = std::make_shared<Subcomponent>(Subcomponent(Multiloop,std::vector<int> {i+1,j-1}));
                        multi->addChild(comp->getChildren()[0]);
                        multi->addChild(s);
                        ene += componentEnergy(multi,i,j);
                        actions.push_back(Action(i,j,Addition,ene));

                    }
                }
            }

            
        }
        break;
        case Hairpin:
        {
            //pure inside 
        
         for( int i= comp->getLeft()+1; i<comp->getRight(); i++){
            for( int j=i+4; j<comp->getRight(); j++){
                if(primary.isPairable(i,j)){
                    double ene=0;
                    ene -= componentEnergy(comp);
                    ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Interior,std::vector<int> {comp->getLeft(), i-1, j +1, comp->getRight()})),i,j);
                    ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Hairpin,std::vector<int> { i+1,  j-1})),i,j);
                    ene += componentEnergy(std::make_shared<Subcomponent>(Subcomponent(Stem,std::vector<int> {i,i,j,j})),i,j);
                    
                    actions.push_back(Action(i,j,Addition,ene));
                }

            }

            
         }   

         //border cases (optimizable)
        {
            int i= comp->getLeft();
            for( int j=i+4; j<=comp->getRight(); j++){
                if(primary.isPairable(i,j)){
                    double ene=  rebuidEnergy(i,j);
                    actions.push_back(Action(i,j,Addition,ene));

                }
            }
            int j = comp->getRight();
            for( int i= comp->getLeft() +1; i<=j-4; i++){
                if(primary.isPairable(i,j)){
                    double ene= rebuidEnergy(i,j);
                    actions.push_back(Action(i,j,Addition,ene));

                }
            }
        }

        

        }
        break;
        case Multiloop:
        {
            // (~~optimizable)
            for( int i : comp->getDefiningBases()){
                for( int j : comp->getDefiningBases()){
                    if( i<j && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));

                    }
                }
            }


            
        }
        break;
        case PseudoknotH:
        {
            // all iside pseudo
            //inside 1-2, 3-4, 5-6
            for(int idx=1;idx<=5;idx +=2){
                auto childrenInside = comp->getChildren(comp->getDefiningBases()[idx],comp->getDefiningBases()[idx+1]);
                for(int j= comp->getDefiningBases()[idx+1]-1;j>comp->getDefiningBases()[idx];j--){
                    if( !childrenInside.empty() && j== childrenInside.back()->getRight()){
                        j=childrenInside.back()->getLeft();
                        childrenInside.pop_back();
                    }else{
                        auto childrenInside2=childrenInside;
                        for( int i = j;i>comp->getDefiningBases()[idx];i--){
                            

                        
                            if(!childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                i=childrenInside2.back()->getLeft();
                                childrenInside2.pop_back();
                            }else if(primary.isPairable(i,j)){
                                double ene=  rebuidEnergy(i,j);
                                actions.push_back(Action(i,j,Addition,ene));
                            }
                        }
                        assert(childrenInside2.empty());
                        
                    }
                }
                assert(childrenInside.empty());
            }
            // 1-2,3-4 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[1]+1;
                int j = comp->getDefiningBases()[4]-1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }
            //3-4,5-6 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[3]+1;
                int j = comp->getDefiningBases()[6]-1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }

            //deletions (ONLY outermost/innermost of bands)
            for(int idx=0; idx<4; idx++){
                double ene = rebuidEnergy(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]]);
                actions.push_back(Action(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]],Deletion,ene));
            }

            
            
        }
        break;
        case PseudoknotK:
        {
            // all iside pseudo
            //inside 1-2, 3-4, 5-6, 7-8 ,9-10
            for(int idx=1;idx<=9;idx +=2){
                auto childrenInside = comp->getChildren(comp->getDefiningBases()[idx],comp->getDefiningBases()[idx+1]);
                for(int j= comp->getDefiningBases()[idx+1]-1;j>comp->getDefiningBases()[idx];j--){
                    if(!childrenInside.empty() && j== childrenInside.back()->getRight()){
                        j=childrenInside.back()->getLeft();
                        childrenInside.pop_back();
                    }else{
                        auto childrenInside2=childrenInside;
                        for( int i = j;i>comp->getDefiningBases()[idx];i--){
                            

                        
                            if(!childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                i=childrenInside2.back()->getLeft();
                                childrenInside2.pop_back();
                            }else if(primary.isPairable(i,j)){
                                double ene=  rebuidEnergy(i,j);
                                actions.push_back(Action(i,j,Addition,ene));
                            }
                        
                        }
                        assert(childrenInside2.size() == 0);
                    
                    }
                }
                assert(childrenInside.size() == 0);
            }

            //deletions (ONLY outermost/innermost of bands)
            for(int idx=0; idx<4; idx++){
                double ene = rebuidEnergy(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]]);
                actions.push_back(Action(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]],Deletion,ene));
            }for(int idx=6; idx<8; idx++){
                double ene = rebuidEnergy(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]]);
                actions.push_back(Action(comp->getDefiningBases()[idx],bp[comp->getDefiningBases()[idx]],Deletion,ene));
            }



            // 1-2 ,3-4 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[1]+1;
                int j = comp->getDefiningBases()[4]-1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }
            //3-4,7-8 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[3]+1;
                int j = comp->getDefiningBases()[8]-1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }
            //7-8,9-10 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[7]+1;
                int j = comp->getDefiningBases()[10]-1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }
            //1-2,9-10 addition to pseudoknotted stem 
            {
                int i = comp->getDefiningBases()[2]-1;
                int j = comp->getDefiningBases()[9]+1;
                if(bp[i]== -1 &&bp[j]== -1 && primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));
                }
            }

        }
        break;
        
        

    }
    return actions;
}

const std::vector<Action> RNASecondaryStructure::SCCrossingActions(std::shared_ptr<Subcomponent> SCfrom,std::shared_ptr<Subcomponent> SCto  ){

    std::vector<Action> actions; 

    switch (SCfrom->getType())
    {
    //exterior-child of Root
    // -E
    // -Ph
    // -Pk
    case Exterior:
        {
        double ene =0;
        if(SCto->getType() == Exterior){
                        for (int  i= SCfrom->getLeft(); i<=SCfrom->getRight(); i++) {
                            for(int j= SCto->getLeft(); j<= SCto->getRight(); j++){
                                if(primary.isPairable(i,j)){
                                ene =rebuidEnergy(i,j);
                                actions.push_back(Action(i,j,Addition,ene));
                                }
                            }
                        }
                    }
                    else if(SCto->getType() == PseudoknotH ){
                            //get inside pseudo
                        //exterior before pseudoknot
                        if(SCfrom->getLeft()< SCto->getLeft()){
                            for(int i=SCfrom->getLeft();i<= SCfrom->getRight(); i++){
                                auto childrenInside = SCto->getChildren(SCto->getDefiningBases()[1],SCto->getDefiningBases()[2]);
                                for(int j= SCto->getDefiningBases()[2]-1;j>SCto->getDefiningBases()[1];j--){
                                    if( !childrenInside.empty() && j== childrenInside.back()->getRight()){
                                        j=childrenInside.back()->getLeft();
                                        childrenInside.pop_back();
                                    }else{
                                        if(primary.isPairable(i,j)){
                                            ene =rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }
                                }
                                assert(childrenInside.empty());
                            }
                            
                            if(SCfrom->getRight() + 1 == SCto->getLeft()){
                                int i = SCfrom->getRight();
                                int j = SCto->getDefiningBases()[5]+1;
                                if(bp[j]==-1 && primary.isPairable(i,j)){
                                    ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }
                            }
                        }
                        else{
                            for(int j=SCfrom->getLeft();j<= SCfrom->getRight(); j++){
                                auto childrenInside = SCto->getChildren(SCto->getDefiningBases()[5],SCto->getDefiningBases()[6]);
                                for(int i= SCto->getDefiningBases()[6]-1;i>SCto->getDefiningBases()[5];i--){
                                    if( !childrenInside.empty() && i== childrenInside.back()->getRight()){
                                        i=childrenInside.back()->getLeft();
                                        childrenInside.pop_back();
                                    }else{
                                        if(primary.isPairable(i,j)){
                                            ene =rebuidEnergy(i,j);
                                            actions.push_back(Action(i,j,Addition,ene));
                                        }
                                    }
                                }
                                assert(childrenInside.empty());
                            }

                            if(SCfrom->getLeft() -1 == SCto->getRight()){
                                int i= SCto->getDefiningBases()[2] -1;
                                int j= SCfrom->getLeft();
                                if(bp[i]==-1 && primary.isPairable(i,j)){
                                    ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }

                            }
                        }
                    }
                    else if(SCto->getType() == PseudoknotK ){
                        //get inside pseudo
                        if(SCfrom->getRight() + 1 == SCto->getLeft()){
                                int i = SCfrom->getRight();
                                int j = SCto->getDefiningBases()[5]+1;
                                if(bp[j]==-1 && primary.isPairable(i,j)){
                                    ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }
                            }
                        if(SCfrom->getLeft() -1 == SCto->getRight()){
                                int i= SCto->getDefiningBases()[6] -1;
                                int j= SCfrom->getLeft();
                                if(bp[i]==-1 && primary.isPairable(i,j)){
                                    ene =rebuidEnergy(i,j);
                                    actions.push_back(Action(i,j,Addition,ene));
                                }
             

                            }
                    }
        }
        break;
    case Multiloop:
        //multiloop-child
        // -Ph
        // -Pk
        {
        if( SCto->getType() == PseudoknotH){
            //from left side of multiloop to 1-2 making K-type
            int idx=0;
            while(idx <SCfrom->getDefiningBases().size() && SCfrom->getDefiningBases()[idx] < SCto->getLeft()){
                int i =SCfrom->getDefiningBases()[idx];
                auto childrenInside = SCto->getChildren(SCto->getDefiningBases()[1],SCto->getDefiningBases()[2]);
                for(int j= SCto->getDefiningBases()[2]-1;j>SCto->getDefiningBases()[1];j--){
                    if( !childrenInside.empty() && j== childrenInside.back()->getRight()){
                        j=childrenInside.back()->getLeft();
                        childrenInside.pop_back();
                    }else{
                        if(primary.isPairable(i,j)){
                            double ene =rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));
                        }
                    }
                }
                
                assert(childrenInside.empty());
                idx++;
            }
            //adding outer base pair to stem 
            if(idx > 0 && SCfrom->getDefiningBases()[idx-1] ==SCto->getLeft()-1){ 
                int i = SCfrom->getDefiningBases()[idx-1];
                int j = SCto->getDefiningBases()[5]+1;
                if(bp[j]==-1 && primary.isPairable(i,j)){
                    double ene =rebuidEnergy(i,j);
                    actions.push_back(Action(i,j,Addition,ene));
                }
            }

            //adding outer base pair to stem 
            if(idx <SCfrom->getDefiningBases().size() && SCfrom->getDefiningBases()[idx] == SCto->getRight()+1){ 
                int i= SCto->getDefiningBases()[2] -1;
                int j= SCfrom->getDefiningBases()[idx];
                if(bp[i]==-1 && primary.isPairable(i,j)){
                    double ene =rebuidEnergy(i,j);
                    actions.push_back(Action(i,j,Addition,ene));
                }
            }


            //from Right side of multiloop to 5-6 making K-type
            while(idx <SCfrom->getDefiningBases().size()){
                int j =SCfrom->getDefiningBases()[idx];
                auto childrenInside = SCto->getChildren(SCto->getDefiningBases()[5],SCto->getDefiningBases()[6]);
                for(int i= SCto->getDefiningBases()[6]-1;i>SCto->getDefiningBases()[5];i--){
                    if( !childrenInside.empty() && i== childrenInside.back()->getRight()){
                        i=childrenInside.back()->getLeft();
                        childrenInside.pop_back();
                    }else{
                        if(primary.isPairable(i,j)){
                            double ene =rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));
                        }
                    }
                }
                assert(childrenInside.empty());

                idx++;
            }
                            




        }else if( SCto->getType() == PseudoknotK){
            
            int idx=0;
            while(idx <SCfrom->getDefiningBases().size() && SCfrom->getDefiningBases()[idx] < SCto->getLeft()){
                idx++;
            }
            //adding outer base pair to stem 
            if(idx > 0 && SCfrom->getDefiningBases()[idx-1] ==SCto->getLeft()-1){ 
                int i = SCfrom->getDefiningBases()[idx-1];
                int j = SCto->getDefiningBases()[5]+1;
                if(bp[j]==-1 && primary.isPairable(i,j)){
                    double ene =rebuidEnergy(i,j);
                    actions.push_back(Action(i,j,Addition,ene));
                }
            }

            //adding outer base pair to stem 
            if(idx <SCfrom->getDefiningBases().size() && SCfrom->getDefiningBases()[idx] == SCto->getRight()+1){ 
                int i= SCto->getDefiningBases()[6] -1;
                int j= SCfrom->getDefiningBases()[idx];
                if(bp[i]==-1 && primary.isPairable(i,j)){
                    double ene =rebuidEnergy(i,j);
                    actions.push_back(Action(i,j,Addition,ene));
                }
            }
        }
        }
        break;

    //pseudoknotH-(grand)child
    // -M
    // -I
    // -H
    // -B
    // -Pk
    // -Ph
    case PseudoknotH:
        {
        int reg1= -1;
        int reg2;
        for(int i =1; i<6; i+=2){

            if(SCfrom->getDefiningBases()[i]<SCto->getLeft()  && SCfrom->getDefiningBases()[i+1]>SCto->getRight()){
                reg1 = SCfrom->getDefiningBases()[i];
                reg2 = SCfrom->getDefiningBases()[i+1];
                break;
            }
        }

        auto childrenInside = SCfrom->getChildren(reg1,reg2);
                if(!childrenInside.empty())
                for(int j= reg2-1;j>reg1;j--){
                    
                    if( !childrenInside.empty() && j== childrenInside.back()->getRight()){
                        j=childrenInside.back()->getLeft();
                        childrenInside.pop_back();
                    }else{
                        
                            
                                
                            if(SCto->getType() == Multiloop){
                                for(int i : SCto->getDefiningBases()){
                                    if(primary.isPairable(i,j)){
                                        double ene=  rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }

                                }

                            }else if(SCto->getType() == Interior){
                                for(int i= SCto->getLeft(); i<=SCto->getDefiningBases()[1]; i++){
                                    if(primary.isPairable(i,j)){
                                        double ene=  rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }
                                for(int i=SCto->getDefiningBases()[2]; i<= SCto->getRight(); i++){
                                    if(primary.isPairable(i,j)){
                                        double ene=  rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }

                            }else if(SCto->getType() ==Hairpin || SCto->getType() == Bulge){
                                for(int i= SCto->getLeft(); i<=SCto->getRight(); i++){
                                    if(primary.isPairable(i,j)){
                                        double ene=  rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }

                            }

                            else if(SCto->getType() ==PseudoknotH){
                                if(j< SCto->getLeft()){
                                    auto childrenInside2 = SCto->getChildren(SCto->getDefiningBases()[1],SCto->getDefiningBases()[2]);
                                    for(int i= SCto->getDefiningBases()[2]-1;i>SCto->getDefiningBases()[1];i--){
                                        if( !childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                            i=childrenInside2.back()->getLeft();
                                            childrenInside2.pop_back();
                                        }else{
                                            if(primary.isPairable(i,j)){
                                                double ene =rebuidEnergy(j,i);
                                                actions.push_back(Action(j,i,Addition,ene));
                                            }
                                        }
                                    }
                                
                                    assert(childrenInside2.empty());
                                    if(j==SCto->getLeft()-1){
                                        int i = SCto->getDefiningBases()[5]+1;
                                        if(bp[i]==-1 && primary.isPairable(i,j)){
                                            double ene =rebuidEnergy(j,i);
                                            actions.push_back(Action(j,i,Addition,ene));
                                        }
                                    }
                                }else{
                                    auto childrenInside2 = SCto->getChildren(SCto->getDefiningBases()[5],SCto->getDefiningBases()[6]);
                                    for(int i= SCto->getDefiningBases()[6]-1;i>SCto->getDefiningBases()[5];i--){
                                        if( !childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                            i=childrenInside2.back()->getLeft();
                                            childrenInside2.pop_back();
                                        }else{
                                            if(primary.isPairable(i,j)){
                                                double ene =rebuidEnergy(i,j);
                                                actions.push_back(Action(i,j,Addition,ene));
                                            }
                                        }
                                    }
                                    assert(childrenInside2.empty());
                                    if(j==SCto->getRight()+1){
                                        int i = SCto->getDefiningBases()[2]-1;
                                        if(bp[i]==-1 && primary.isPairable(i,j)){
                                            double ene =rebuidEnergy(j,i);
                                            actions.push_back(Action(j,i,Addition,ene));
                                        }
                                    }
                                }

                            }else if(SCto->getType() == PseudoknotK){
                                if(j==SCto->getLeft()-1){
                                    int i = SCto->getDefiningBases()[5]+1;
                                    if(bp[i]==-1 && primary.isPairable(i,j)){
                                        double ene =rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }else if(j==SCto->getRight()+1){
                                    int i = SCto->getDefiningBases()[6]-1;
                                    if(bp[i]==-1 && primary.isPairable(i,j)){
                                        double ene =rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }
                                
                                
                                
                            }
                        
                        
                        
                    }
                }
                assert(childrenInside.empty());
        }
        break;   
    
      //pseudoknotK-(grand)child
    // -M
    // -I
    // -H
    // -B
    // -Pk
    // -Ph
    case PseudoknotK:
        {
        int reg1= -1;
        int reg2;
        for(int i =1; i<9; i+=2){

            if(SCfrom->getDefiningBases()[i]<SCto->getLeft()  && SCfrom->getDefiningBases()[i+1]>SCto->getRight()){
                reg1 = SCfrom->getDefiningBases()[i];
                reg2 = SCfrom->getDefiningBases()[i+1];
                break;
            }
        }

        auto childrenInside = SCfrom->getChildren(reg1,reg2);
                if(!childrenInside.empty())
                for(int j= reg2-1;j>reg1;j--){
                    
                    if( !childrenInside.empty() && j== childrenInside.back()->getRight()){
                        j=childrenInside.back()->getLeft();
                        childrenInside.pop_back();
                    }else{
                        
                            
                                
                            if(SCto->getType() == Multiloop){
                                for(int i : SCto->getDefiningBases()){
                                    if(primary.isPairable(i,j)){
                                        double ene=  rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }

                                }

                            }else if(SCto->getType() == Interior){
                                for(int i= SCto->getLeft(); i<=SCto->getDefiningBases()[1]; i++){
                                    if(primary.isPairable(i,j)){
                                        double ene=  rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }
                                for(int i=SCto->getDefiningBases()[2]; i<= SCto->getRight(); i++){
                                    if(primary.isPairable(i,j)){
                                        double ene=  rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }

                            }else if(SCto->getType() ==Hairpin || SCto->getType() == Bulge){
                                for(int i= SCto->getLeft(); i<=SCto->getRight(); i++){
                                    if(primary.isPairable(i,j)){
                                        double ene=  rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }

                            }

                            else if(SCto->getType() ==PseudoknotH){
                                if(j< SCto->getLeft()){
                                    auto childrenInside2 = SCto->getChildren(SCto->getDefiningBases()[1],SCto->getDefiningBases()[2]);
                                    for(int i= SCto->getDefiningBases()[2]-1;i>SCto->getDefiningBases()[1];i--){
                                        if( !childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                            i=childrenInside2.back()->getLeft();
                                            childrenInside2.pop_back();
                                        }else{
                                            if(primary.isPairable(i,j)){
                                                double ene =rebuidEnergy(j,i);
                                                actions.push_back(Action(j,i,Addition,ene));
                                            }
                                        }
                                    }
                                
                                    assert(childrenInside2.empty());
                                    if(j==SCto->getLeft()-1){
                                        int i = SCto->getDefiningBases()[5]+1;
                                        if(bp[i]==-1 && primary.isPairable(i,j)){
                                            double ene =rebuidEnergy(j,i);
                                            actions.push_back(Action(j,i,Addition,ene));
                                        }
                                    }
                                }else{
                                    auto childrenInside2 = SCto->getChildren(SCto->getDefiningBases()[5],SCto->getDefiningBases()[6]);
                                    for(int i= SCto->getDefiningBases()[6]-1;i>SCto->getDefiningBases()[5];i--){
                                        if( !childrenInside2.empty() && i== childrenInside2.back()->getRight()){
                                            i=childrenInside2.back()->getLeft();
                                            childrenInside2.pop_back();
                                        }else{
                                            if(primary.isPairable(i,j)){
                                                double ene =rebuidEnergy(i,j);
                                                actions.push_back(Action(i,j,Addition,ene));
                                            }
                                        }
                                    }
                                    assert(childrenInside2.empty());
                                    if(j==SCto->getRight()+1){
                                        int i = SCto->getDefiningBases()[2]-1;
                                        if(bp[i]==-1 && primary.isPairable(i,j)){
                                            double ene =rebuidEnergy(j,i);
                                            actions.push_back(Action(j,i,Addition,ene));
                                        }
                                    }
                                }

                            }else if(SCto->getType() == PseudoknotK){
                                if(j==SCto->getLeft()-1){
                                    int i = SCto->getDefiningBases()[5]+1;
                                    if(bp[i]==-1 && primary.isPairable(i,j)){
                                        double ene =rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }else if(j==SCto->getRight()+1){
                                    int i = SCto->getDefiningBases()[6]-1;
                                    if(bp[i]==-1 && primary.isPairable(i,j)){
                                        double ene =rebuidEnergy(i,j);
                                        actions.push_back(Action(i,j,Addition,ene));
                                    }
                                }
                                
                                
                                
                            }
                        
                        
                        
                    }
                }
                assert(childrenInside.empty());
        }
        break; 
    default:
        actions = makePseudo(SCfrom,SCto);
        break;
    }
    

    

    
return actions;
  
}

double RNASecondaryStructure::rebuidEnergy(int i, int j){
    
                
    double ene=0;
    auto tmp2= bp;
    
    auto tmp = root;
    ene -= evaluateSubstructureEnergy(root);

    if(bp[i]==j){
        bp[j]=-1;
        bp[i]=-1;
    }else{
        bp[i]=j;
        bp[j]=i;
    }
    
    docomposeStructure();
    ene += evaluateSubstructureEnergy(root);
    if(bp[i]==j){
        bp[j]=-1;
        bp[i]=-1;
    }else{
        bp[i]=j;
        bp[j]=i;
    }
    
    
    root=tmp;
    docomposeStructure();
    assert(bp==tmp2);
    return ene;
}

std::vector<Action> RNASecondaryStructure::makePseudo(std::shared_ptr<Subcomponent> first ,std::shared_ptr<Subcomponent> second ){
    assert(first->getType() != Stem && second->getType() != Stem );
    std::vector<Action> actions;
    //pseudoknots already considered in component action
    assert(second->getType() != PseudoknotK && second->getType() != PseudoknotH);
    assert(first->getType() != PseudoknotK && first->getType() != PseudoknotH);

    switch (first->getType())
    {
    case Exterior:
        {
            //double ene =0;
            

             if (second->getType() == Interior) {
                double ene= 0; //- componentEnergy(root);

                for(int j= first->getLeft(); j<=first->getRight(); j++){
                    for(int i= second->getLeft(); i<= second->getDefiningBases()[1]; i++){
                        if(  primary.isPairable(i,j)){
                                double ene=  rebuidEnergy(i,j);
                                actions.push_back(Action(i,j,Addition,ene));
                        }
                    }           
                    for(int i= second->getDefiningBases()[2]; i<=second->getRight() ; i++){   
                        if(  primary.isPairable(i,j)){
                                double ene=  rebuidEnergy(i,j);
                                actions.push_back(Action(i,j,Addition,ene));
                        }
                    }
                }
                
 
            }else if (second->getType() == Multiloop){
                // optimizable
                for(int j= first->getLeft(); j<=first->getRight(); j++){
                    for(int i : second->getDefiningBases()){
                        if(  primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                    }
                }
            }
            //hairpin or bulge
            else{
                // optimizable
                for(int j= first->getLeft(); j<=first->getRight(); j++){
                    for(int i = second->getLeft();i<=second->getRight();i++){
                        if( primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                    }
                }

            }
        }
        break;
    case Interior:
    //all optimizable
    {
        if(second->getType()==Interior){
            for(int i= first->getLeft(); i<= first->getDefiningBases()[1]; i++){
                for(int j= second->getLeft(); j<= second->getDefiningBases()[1]; j++){
                    if(  primary.isPairable(i,j)){
                            double ene= rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }
                for(int j= second->getDefiningBases()[2]; j<=second->getRight() ; j++){
                    if(  primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }

            }
            for(int i= first->getDefiningBases()[2]; i<=first->getRight() ; i++){
                for(int j= second->getLeft(); j<= second->getDefiningBases()[1]; j++){
                    if( primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }
                for(int j= second->getDefiningBases()[2]; j<=second->getRight() ; j++){
                    if(  primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }
            }

        }else if(second->getType() == Multiloop){

            for(int i= first->getLeft(); i<= first->getDefiningBases()[1]; i++){
                for(int j : second->getDefiningBases()){
                    if(  primary.isPairable(i,j)){
                            double ene= rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }
                

            }
            for(int i= first->getDefiningBases()[2]; i<=first->getRight() ; i++){
                for(int j : second->getDefiningBases()){
                    if(  primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }
                
            }
        }
        //bulge or hairpin (should NOT be exterior ( exterior would work))
        else{
            for(int i= first->getLeft(); i<= first->getDefiningBases()[1]; i++){
                for(int j = second->getLeft(); j<= second->getRight();j++){
                    if( primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }
                

            }
            for(int i= first->getDefiningBases()[2]; i<=first->getRight() ; i++){
                for(int j = second->getLeft(); j<= second->getRight();j++){
                    if(  primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }
                
            }

        }
    }
    break;

    case Multiloop:
    //all optimizable
    {
        if( second->getType() == PseudoknotK){

        }else if (second->getType() ==  PseudoknotH){

        }else if (second->getType() ==  Interior){

            for(int j : first->getDefiningBases()){
                for(int i= second->getLeft(); i<= second->getDefiningBases()[1]; i++){
                    if(  primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }           
                for(int i= second->getDefiningBases()[2]; i<=second->getRight() ; i++){   
                    if(  primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                
                
                }
            }
            
        }else if (second->getType() ==  Multiloop){

            for(int i : first->getDefiningBases()){
                for(int j : second->getDefiningBases()){
                    if( primary.isPairable(i,j)){
                            double ene= rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }
                
            }
        }else{
            for(int i : first->getDefiningBases()){
                for(int j = second->getLeft(); j<= second->getRight(); j++){
                    if( primary.isPairable(i,j)){
                            double ene= rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }
                
            }

        }
    }
    break;
    
    case PseudoknotH:
    {
        if( second->getType() == PseudoknotK){

        }else if (second->getType() ==  PseudoknotH){

        }else if (second->getType() ==  Interior){

        }else if (second->getType() ==  Multiloop){

        }else{
            
        }
    }
    break;
    
    case PseudoknotK:
    {
        if( second->getType() == PseudoknotK){

        }else if (second->getType() ==  PseudoknotH){

        }else if (second->getType() ==  Interior){


        }else if (second->getType() ==  Multiloop){

        }else{
            
        }
    }
    break;
    //hairpin or bulge
    default:
    {
        if (second->getType() ==  Interior){
            for(int j= first->getLeft(); j<=first->getRight(); j++){
                for(int i= second->getLeft(); i<= second->getDefiningBases()[1]; i++){
                    if( primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                }           
                for(int i= second->getDefiningBases()[2]; i<=second->getRight() ; i++){   
                    if(  primary.isPairable(i,j)){
                            double ene=  rebuidEnergy(i,j);
                            actions.push_back(Action(i,j,Addition,ene));

                        }
                
                
                }
            }
        }else if (second->getType() == Multiloop){
            // optimizable
            for(int j= first->getLeft(); j<=first->getRight(); j++){
                for(int i : second->getDefiningBases()){
                    if(  primary.isPairable(i,j)){
                        double ene=  rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));

                    }
                }
            }
        }
        //hairpin or bulge
        else{
            // optimizable
            for(int j= first->getLeft(); j<=first->getRight(); j++){
                for(int i = second->getLeft();i<=second->getRight();i++){
                    if(  primary.isPairable(i,j)){
                        double ene= rebuidEnergy(i,j);
                        actions.push_back(Action(i,j,Addition,ene));

                    }
                }
            }

        }
    }
        break;
    }
    return actions;
}

bool RNASecondaryStructure::increaseFoldingLength(){
    return primary.increaseFoldingLength();
}

int RNASecondaryStructure::getFoldingLength(){
    return primary.getFoldingLength();
}

std::vector<Action> RNASecondaryStructure::additonalActions() {

    auto Ebase = std::make_shared<Subcomponent>(Subcomponent(Exterior,std::vector<int> { primary.getFoldingLength()}));

    return subcomponentActions(Ebase,root->getChildren());
    
}

std::string getString(std::shared_ptr<Subcomponent> comp){
    std::string ret ="";
    auto curr=comp;
    while(curr->getChildren().size() ==1){
        if(curr->getType()==Stem){
            ret += 'S';
        }else if(curr->getType()==Interior){
            ret += 'I';
        }else if(curr->getType()==Bulge){
            ret += 'B';
        }else break;

        curr= curr->getChildren()[0];
    }
    if(curr->getChildren().size() ==0){
         if(curr->getType()==Hairpin){
            ret += 'H';
        }else if(curr->getType()==PseudoknotH  ){
            ret += "Ph";
        }else if(curr->getType()==PseudoknotK  ){
            ret += "Pk";
        }
    }else{
        if(curr->getType()==Multiloop){
            ret += "M(";
            for(auto c :curr->getChildren()){
                ret += getString(c);
            }
            ret +=')';
        }else if(curr->getType()==PseudoknotH  ){
            ret += "Ph(";
            auto vec = curr->getDefiningBases();
            for(int idx=1;idx<=5;idx +=2){
                ret +='(';
                if(curr->getChildren(vec[idx],vec[idx+1]).size()==0){
                    // ret += 'L';
                }else{
                    for(auto c :curr->getChildren(vec[idx],vec[idx+1])){
                        ret += getString(c);
                    }
                }
                ret+=')';
            }
            ret+= ')';
        }else if(curr->getType()==PseudoknotK  ){
            ret += "Pk(";
            auto vec = curr->getDefiningBases();
            for(int idx=1;idx<=9;idx +=2){
                ret +='(';
                if(curr->getChildren(vec[idx],vec[idx+1]).size()==0){
                    // ret += 'L';
                }else{
                    for(auto c :curr->getChildren(vec[idx],vec[idx+1])){
                        ret += getString(c);
                    }
                }
                ret +=')';
            }
            ret+= ')';
        }
    }
    return ret;
}

std::string RNASecondaryStructure::getSCTree(){
    std::string ret="";
    int curRight=-1;  
    std::vector<std::shared_ptr<Subcomponent>> rStack;
    std::vector<int> extRight;
    for(auto sc : root->getChildren()){
        if(sc->getType() == Exterior){
            int r = sc->getRight();
            extRight.push_back(r);
        }else {
            rStack.push_back(sc);
        }
    }
    // std::reverse(rStack.begin(),rStack.end());
    std::sort(extRight.begin(),extRight.end(), std::greater<int>() );
    while(/* curRight<primary.getLength()-1 */ !rStack.empty()){
        auto curr = rStack.back();
        if(extRight.empty() || curr->getRight()<extRight.back()){
            ret += getString(curr);
            rStack.pop_back();
        }else{
            ret +='E';
            extRight.pop_back();
        }
    }
    if(! extRight.empty()){
        ret +="E";
    }
    return ret;
}

std::string RNASecondaryStructure::getSequence(){
    return primary.getString();
}

//must be tested
std::vector<std::shared_ptr<Subcomponent>> SCWithin(std::shared_ptr<Subcomponent> root ,int base1, int base2)
{

    //stack for "root" subcomponents i.e. Root, Multiloop or PseudoknotH/K
    std::vector<std::shared_ptr<Subcomponent>> rStack;
    //queue for "root" subcomponents i.e. Root, Multiloop or PseudoknotH/K
    std::deque<std::shared_ptr<Subcomponent>> rQueue;
    
    
    //current subcomponent
    std::shared_ptr<Subcomponent> current;
    //parent for subcomponents in sStack
    std::shared_ptr<Subcomponent> rParent;
    

   std::vector<std::shared_ptr<Subcomponent>> within;

    rQueue.push_back(root);

    while(!rQueue.empty()){
        
        rParent =rQueue.front();
        rQueue.pop_front();

        // for pseudoknots siblings will be divided into regions 
        std::vector<std::deque<std::shared_ptr<Subcomponent>> >sRegions;

        if(rParent->getType() == PseudoknotK ){
            for(int idx=1;idx<=9;idx +=2){
                auto tmp = rParent->getChildren(rParent->getDefiningBases()[idx],rParent->getDefiningBases()[idx+1]);
                std::deque<std::shared_ptr<Subcomponent> > q(tmp.begin(),tmp.end());
                sRegions.push_back(q);
            }
        }
        else if( rParent->getType() == PseudoknotH){
            for(int idx=1;idx<=5;idx +=2){
                auto tmp = rParent->getChildren(rParent->getDefiningBases()[idx],rParent->getDefiningBases()[idx+1]);
                std::deque<std::shared_ptr<Subcomponent> > q(tmp.begin(),tmp.end());
                sRegions.push_back(q);
            }
        }
        else{

            std::deque<std::shared_ptr<Subcomponent>> sQueue;
            for( auto c : rParent->getChildren()){
                sQueue.push_back(c);
            }
            if(rParent->getType()==Root)  std::reverse(sQueue.begin(),sQueue.end());
            sRegions.push_back(sQueue);
        }
        
        

        for( auto sQueue : sRegions)
        while(!sQueue.empty()){
            //Exterior Stem or Pseudoknot enumerate action for this SC
            current=sQueue.front();
            sQueue.pop_front();


            if(current->getType()==PseudoknotH || current->getType()==PseudoknotK){ 
                rQueue.push_back(current);
            }
            //add subcomponent to list if it is within the region
            if(!(current->getLeft() >base2 || current->getRight() < base1) ) within.push_back(current);

            //first subcomponent level untill hairpin loop or multiloop 
            //cannot be pseudoknot 
            if(current->getType() == Stem){

                current = current->getChildren()[0];
                if(current->getType()==Multiloop){
                    //add subcomponent to list if it is within the region
                    if(!(current->getLeft() >base2 || current->getRight() < base1) ) within.push_back(current);
                    rQueue.push_back(current);
                        
                }
                else{
                    //add subcomponent to list if it is within the region
                    if(!(current->getLeft() >base2 || current->getRight() < base1) ) within.push_back(current);

                    while(current->getChildren().size()==1){
                    
                        current = current->getChildren()[0];

                        if(current->getType()==Multiloop){
                            rQueue.push_back(current);
                            break;
                        }

                        //add subcomponent to list if it is within the region
                        if(!(current->getLeft() >base2 || current->getRight() < base1) ) within.push_back(current);


                                    
                                
                    }
                } 

                
            }
        }

    }
    return within;
}
//must be tested
bool RNASecondaryStructure::canFormPseudo(std::shared_ptr<Subcomponent> sc1,std::shared_ptr<Subcomponent> sc2 ){
    assert(sc1->getType()!=Stem);
    assert(sc2->getType()!=Stem);
    //other is grandparent of other

    if(sc1->getLeft()<sc2->getLeft() && sc1->getRight()> sc2->getRight() ){
        std::shared_ptr<Subcomponent> parent;
        for(auto c : sc1->getChildren()){
            if(c->getLeft()<sc2->getLeft() && c->getRight()> sc2->getRight()){
                if(c==sc2 ) return true;  //not working
                else if(c->getType()==Stem && c->getChildren()[0]==sc2) return true;
                else return false;
            }
        }
    }else if(sc2->getLeft()<sc1->getLeft() && sc2->getRight()> sc1->getRight() ){
        std::shared_ptr<Subcomponent> parent;
        for(auto c : sc2->getChildren()){
            if(c->getLeft()<sc1->getLeft() && c->getRight()> sc1->getRight()){
                if(c==sc1) return true;
                else if(c->getType()==Stem && c->getChildren()[0]==sc1) return true;
                else return false;
            }
        }
    }

    //find common ancestor 
    auto curr = root;
    bool bothInside =true;
    while(bothInside){
        bothInside=false;
        for(auto c : curr->getChildren()){
            if(c->getLeft()<= std::min(sc1->getLeft(),sc2->getLeft()) && c->getRight()>= std::max(sc1->getRight(),sc2->getRight())){
                curr = c;
                bothInside=true;
                break;
            }
        }
    }
    // chek if common grandparent (or exterior or pseudo)
    if(curr->getType() ==Root){
        //if exterior
        if(sc1->getType()==Exterior){
            if(sc2->getType()==Exterior) return true;
            for(std::shared_ptr<Subcomponent> c : root->getChildren()){
                if(c->getLeft()<=sc2->getLeft() && c->getRight()>= sc2->getRight()){
                    if(c->getType()==PseudoknotH && c==sc2) return true;
                    else if(c->getType()==PseudoknotK && c==sc2){
                        if(sc1->getLeft()==sc2->getRight() +1 ||sc1->getRight() +1 == sc2->getLeft())return true;
                        else return false;
                    }else if(c->getType()==Stem){
                        if(c->getChildren()[0]==sc2) return true;
                        else return false;
                    }
                }
            }
            return false;
        }else if(sc2->getType()==Exterior){
            
            for(std::shared_ptr<Subcomponent> c : root->getChildren()){
                if(c->getLeft()<=sc1->getLeft() && c->getRight()>= sc1->getRight()){
                    if(c->getType()==PseudoknotH && c==sc1) return true;
                    else if(c->getType()==PseudoknotK && c==sc1){
                        if(sc2->getLeft()==sc1->getRight() +1 ||sc2->getRight() +1 == sc1->getLeft())return true;
                        else return false;
                    }else if(c->getType()==Stem){
                        if(c->getChildren()[0]==sc1) return true;
                        else return false;
                    }
                }
            }
            return false;
        }
        else{
            bool possible=false;
            for(std::shared_ptr<Subcomponent> c : curr->getChildren()){
                if(c->getLeft()<=sc1->getLeft() && c->getRight()>= sc1->getRight()){
                    if(c->getType()==Stem){
                        if(c->getChildren()[0]==sc1) possible= true;
                        else return false;
                    }else return false;
                }
               
                if(c->getLeft()<=sc2->getLeft() && c->getRight()>= sc2->getRight()){
                    if(c->getType()==Stem){
                        if(c->getChildren()[0]==sc2 && possible) return true;
                        else return false;
                    }
                }
            }
        }

    }else if( curr->getType()==Multiloop){
        bool possible=false;
            for(std::shared_ptr<Subcomponent> c : curr->getChildren()){
                if(c->getLeft()<=sc1->getLeft() && c->getRight()>= sc1->getRight()){
                    if(c->getType()==Stem){
                        if(c->getChildren()[0]==sc1) possible= true;
                        else return false;
                    }else return false;
                }
               
                if(c->getLeft()<=sc2->getLeft() && c->getRight()>= sc2->getRight()){
                    if(c->getType()==Stem){
                        if(c->getChildren()[0]==sc2 && possible) return true;
                        else return false;
                    }
                }
            }
    }else if(curr->getType()==PseudoknotH){
        int reg1= -1;
        int reg2;
        for(int i =1; i<6; i+=2){

        
            if(curr->getDefiningBases()[i]<sc1->getLeft() && curr->getDefiningBases()[i]<sc2->getLeft() && curr->getDefiningBases()[i+1]>sc1->getRight()&& curr->getDefiningBases()[i+1]>sc2->getRight()){
                reg1 = curr->getDefiningBases()[i];
                reg2 = curr->getDefiningBases()[i+1];
                break;
            }
        }
        if(reg1 == -1) return false;
        else{
            bool possible=false;
            for(std::shared_ptr<Subcomponent> c : curr->getChildren(reg1,reg2)){
                if(c->getLeft()<=sc1->getLeft() && c->getRight()>= sc1->getRight()){
                    if(c->getType()==Stem){
                        if(c->getChildren()[0]==sc1) possible= true;
                        else return false;
                    }else return false;
                }
               
                if(c->getLeft()<=sc2->getLeft() && c->getRight()>= sc2->getRight()){
                    if(c->getType()==Stem){
                        if(c->getChildren()[0]==sc2 && possible) return true;
                        else return false;
                    }
                }
            }
        }


    }else if(curr->getType()==PseudoknotK){
        int reg1= -1;
        int reg2;
        for(int i =1; i<8; i+=2){

            if(curr->getDefiningBases()[i]<sc1->getLeft() && curr->getDefiningBases()[i]<sc2->getLeft() && curr->getDefiningBases()[i+1]>sc1->getRight()&& curr->getDefiningBases()[i+1]>sc2->getRight()){
                reg1 = curr->getDefiningBases()[i];
                reg2 = curr->getDefiningBases()[i+1];
                break;
            }
        }
        if(reg1 == -1) return false;
        else{
            bool possible=false;
            for(std::shared_ptr<Subcomponent> c : curr->getChildren(reg1,reg2)){
                if(c->getLeft()<=sc1->getLeft() && c->getRight()>= sc1->getRight()){
                    if(c->getType()==Stem){
                        if(c->getChildren()[0]==sc1) possible= true;
                        else return false;
                    }else return false;
                }
               
                if(c->getLeft()<=sc2->getLeft() && c->getRight()>= sc2->getRight()){
                    if(c->getType()==Stem){
                        if(c->getChildren()[0]==sc2 && possible) return true;
                        else return false;
                    }
                }
            }
        }

    }
    return false;
}




void RNASecondaryStructure::updateActions( Action lastAction){ 

   //Finding old and new sc for the update 

    auto oldSCsa = SCWithin(oldRoot,std::min(lastAction.i,lastAction.j),std::max(lastAction.i,lastAction.j));
    auto newSCsa = SCWithin(oldRoot,std::min(lastAction.i,lastAction.j),std::max(lastAction.i,lastAction.j));
    SCLess le;
    std::sort(oldSCsa.begin(),oldSCsa.end(),le);
    std::sort(newSCsa.begin(),newSCsa.end(),le);
    // subcomponents that got deleted with the last action
    std::vector<std::shared_ptr<Subcomponent>> oldSC;
    // subcomponents that got generated with last action
    std::vector<std::shared_ptr<Subcomponent>> newSC;

    std::set_difference(oldSCsa.begin(),oldSCsa.end(),newSCsa.begin(),newSCsa.end(),std::inserter(oldSC, oldSC.begin()),le);
    std::set_difference(newSCsa.begin(),newSCsa.end(),oldSCsa.begin(),oldSCsa.end(),std::inserter(newSC, newSC.begin()),le);
    
    // Make new SC associated actions
    std::vector<std::shared_ptr<Subcomponent>> newCrossterms;
    for( auto sc :oldSC){
        
        if(sc->getType()!=Stem){
            if(std::get<1>(map3[*sc]).size()>1){
                newCrossterms.insert(newCrossterms.end(),std::get<1>(map3[*sc]).begin()+1,std::get<1>(map3[*sc]).end()); 
            }
        }
    }
    std::sort(newCrossterms.begin(),newCrossterms.end(),le);
    newCrossterms.erase( std::unique( newCrossterms.begin(), newCrossterms.end() ), newCrossterms.end() );

    // set of all possible new crossterms for all newSC (filtered for old sc and stems), has to be individually filtered for each newSC
    std::vector<std::shared_ptr<Subcomponent>> newCrosstermsFiltered;
    std::set_difference(newCrossterms.begin(),newCrossterms.end(),oldSC.begin(),oldSC.end(),std::inserter(newCrosstermsFiltered ,newCrosstermsFiltered.begin()),le);
    
    //vector of crossterms filtered for each newSC (in same order as sc in newSC)
    std::vector<std::vector<std::shared_ptr<Subcomponent>> > newActionSC;
    
    for( const std::shared_ptr<Subcomponent> sc: newSC){
        std::vector<std::shared_ptr<Subcomponent>> tmp;
        // tmp.push_back(sc);
        std::copy_if(newSC.begin(),newSC.end(),std::back_inserter(tmp),[sc,le](const std::shared_ptr<Subcomponent> a){ return a->getType()!= Stem && le(sc,a);});
        std::copy_if(newCrosstermsFiltered.begin(),newCrosstermsFiltered.end(),std::back_inserter(tmp),[sc,le](const std::shared_ptr<Subcomponent> a){ return  le(sc,a);});
        if(tmp.size()>0)
        std::remove_if(tmp.begin()/* +1 */,tmp.end(),[this,sc](const std::shared_ptr<Subcomponent> a){return !canFormPseudo(sc,a);});
        newActionSC.push_back(tmp);
    }


    //doing the actual update

    //delete all old actions and add those SC ,that the deleted actions are associated to, to upd 
    std::unordered_set<std::shared_ptr<Subcomponent>> upd;

   

    //delete old sc
    for( auto sc : oldSC){
        for( auto tup : std::get<0>( map3[*sc]) ){
            //add SC ((dim1) that the crossterm action between SC and sc are associated to)  to update set for addition of new actions with the new Scs
            // upd.insert( *std::get<0>(tup));
            upd.insert( std::get<0>(tup));
            // delete actions associated to SC that are sc and subtract the total flux of dim2 elem from dim1 elem
            auto dim1elem = *std::get<2>(map3[*std::get<0>(tup)]);
            std::get<1>(dim1elem) -= std::get<1>(*std::get<1>(tup));
            std::get<0>(dim1elem).erase(std::get<1>(tup));

            // delet sc From map3[SC][1]
            for(auto sc1it = std::get<1>(map3[*std::get<0>(tup)]).begin(); sc1it!=std::get<1>(map3[*std::get<0>(tup)]).end(); sc1it++ ){
                if(*sc1it==sc){
                    std::get<1>(map3[*std::get<0>(tup)]).erase(sc1it);
                    break;
                }
            }
        }
        for(int idx =1;idx< std::get<1>( map3[*sc]).size(); idx++){
            for( auto tupit= std::get<0>(map3[*std::get<1>( map3[*sc])[idx]]).begin(); tupit!= std::get<0>(map3[*std::get<1>( map3[*sc])[idx]]).end(); tupit++ ){
                if(std::get<0>(*tupit) ==sc){
                    std::get<0>(map3[*std::get<1>( map3[*sc])[idx]]).erase(tupit);
                    break;
                }
            }
        }
        dim1.erase(std::get<2>(map3[*sc]));
        map3.erase(*sc);

    }
    //add new to map3 and dim1
    // add actions within sc

    for( auto idx =0; idx< newSC.size(); idx++){
        auto sc = newSC[idx];
        //adding the inner action to dim1
        saveSCActions(sc);
    }

    //add crossterm actions of newSC

    for( auto idx =0; idx< newSC.size(); idx++){
        auto sc = newSC[idx];
        /* double dim1RateSum =0;
        for(auto asc : newActionSC[idx]){
            //add crossing actions and rate sum (dim2 elements to dim1)
            auto crossA= SCCrossingActions(sc,asc);
            double dim2RateSum=0;
            for( auto a : crossA){
                dim2RateSum +=std::exp(-a.getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) );
            }
            std::get<0>(* std::get<2>(map3[*sc])).push_back(std::make_tuple(crossA,dim2RateSum));
            dim1RateSum += dim2RateSum;
            //updating map3
            //asc
            std::list<std::tuple<std::vector<Action>,double> >::iterator dim2ptr= std::get<0>(* std::get<2>(map3[*sc])).end();
            std::advance(dim2ptr,-1);
            std::get<0>(map3[*asc]).push_back(std::make_tuple(sc,dim2ptr));
            //sc
            std::get<1>(map3[*sc]).push_back(asc);
        }
        // updating dim1 rate sum
        std::get<1>(* std::get<2>(map3[*sc])) += dim1RateSum; */

        for(auto asc : newActionSC[idx]){
            saveCrossingActions(sc,asc);
        }
        

    }

    //Update sc in upd with crossterm actions to newSC
    for(auto sc : upd){

        std::vector<std::shared_ptr<Subcomponent>> tmp;
        std::copy_if(newSC.begin(),newSC.end(),std::back_inserter(tmp),[sc,le](const std::shared_ptr<Subcomponent> a){ return a->getType()!= Stem && le(sc,a);});
        std::remove_if(tmp.begin()/* +1 */,tmp.end(),[this,sc](const std::shared_ptr<Subcomponent> a){return !canFormPseudo(sc,a);});
        
        /* double deltaDim1RateSum =0;
        for(auto asc : tmp){
            //add crossing actions and rate sum (dim2 elements to dim1)
            auto crossA= SCCrossingActions(sc,asc);
            double dim2RateSum=0;
            for( auto a : crossA){
                dim2RateSum +=std::exp(-a.getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) );
            }
            std::get<0>(* std::get<2>(map3[*sc])).push_back(std::make_tuple(crossA,dim2RateSum));
            deltaDim1RateSum += dim2RateSum;

            //updating map3
            //asc
            std::list<std::tuple<std::vector<Action>,double> >::iterator dim2ptr= std::get<0>(* std::get<2>(map3[*sc])).end();
            std::advance(dim2ptr,-1);
            std::get<0>(map3[*asc]).push_back(std::make_tuple(sc,dim2ptr));
            //sc
            std::get<1>(map3[*sc]).push_back(asc);
        }
        // updating dim1 rate sum
        std::get<1>(* std::get<2>(map3[*sc])) += deltaDim1RateSum; */

        for(auto asc : tmp){
            saveCrossingActions(sc,asc);
        }
        
        
    }

    // be happy 
}