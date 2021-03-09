//class for the closed region tree to divide the secondary struture in loops

#pragma once
#include "Node.hpp"

class CRTree {

    public:
        CRTree(Node r);

        void addToTree(Node & n);

        std::shared_ptr<Node> getLLB();
        std::shared_ptr<Node> getRoot();
    private:

        std::shared_ptr<Node> root;
};