
#include "CRTree.hpp"

CRTree::CRTree(Node r):  root(std::make_shared<Node>(r))  { }

std::shared_ptr<Node> CRTree::getLLB() {
    std::shared_ptr<Node> llb = root;
    for( auto it : root->children ){
        if(it->left > llb->left) llb = it;
    }
    return llb;
}

void CRTree::addToTree(Node & ijr) {
    std::shared_ptr<Node> ij= std::make_shared<Node>(ijr);
    std::shared_ptr<Node> ab = getLLB();
    
    while (ab->left > ij->left) {
        
        ij->children.push_back(ab);      
        root->children.pop_back();       
        ab=getLLB();

    }
    

     
    root->children.push_back(ij);

}

std::shared_ptr<Node> CRTree::getRoot() {
    return root;
}

