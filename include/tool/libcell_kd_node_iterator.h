#ifndef LIBCELL_KD_NODE_ITERATOR_H_
#define LIBCELL_KD_NODE_ITERATOR_H_

#include<tool/libcell_tools_kd_node.h>


typedef struct KD_Trav{
    KD_Node*tree;
    KD_Node* it;
    //void*(*first)(struct RB_Trav*);
    //void*(*second)(struct RB_Trav*);
    void*prop;

}KD_Trav;
void KD_Trav_init(KD_Trav*);
KD_Trav kd_tree_begin(KD_Node* kn);
KD_Trav kd_tree_end(KD_Node* kn);

KD_Trav kd_tree_rbegin(KD_Node*kn);
KD_Trav operator++(KD_Trav& trav);
KD_Trav operator++(KD_Trav& trav,int);
bool operator!=(const KD_Trav&,const KD_Trav&);
#endif