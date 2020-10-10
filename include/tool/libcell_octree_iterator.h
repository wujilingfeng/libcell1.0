#ifndef LIBCELL_OCTREE_ITERATOR_H_
#define LIBCELL_OCTREE_ITERATOR_H_
#include<tool/libcell_tools_octree.h>

typedef struct Octree_Trav{
    OC_Node* tree;
    OC_Node* it;
    void* prop;
}Octree_Trav;

void Octree_Trav_init(Octree_Trav* oct);
Octree_Trav octree_begin(OC_Node* ocn);
Octree_Trav octree_end(OC_Node*ocn);

Octree_Trav operator++(Octree_Trav& trav);
Octree_Trav operator++(Octree_Trav& trav,int);
bool operator!=(const Octree_Trav&trav1,const Octree_Trav&trav2);
#endif