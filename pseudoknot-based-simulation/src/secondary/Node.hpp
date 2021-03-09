
#pragma once 
#include <memory>
#include <vector>


class Node {
    public:
        Node( int i, int j);
        
        const int left;
        const int right;
        std::vector<std::shared_ptr<Node>> children;


};