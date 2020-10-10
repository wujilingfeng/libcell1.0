#ifndef LIBCELL_TOOLS_OCTREE_H_
#define LIBCELL_TOOLS_OCTREE_H_
#include<Mesh_IO/Mesh_IO.h>
#include "math.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct OC_Node{
    struct OC_Node* parent;
    struct OC_Node *children[8];
    double loc_min[3],loc_max[3]; 
    Node* value; 
    void* prop;
}OC_Node;

void OC_Node_init(OC_Node* ocn);

void oc_node_free(OC_Node*ocn);

//对叶子结点再分割
//输入的ocn是一个叶子结点
//返回含有顶点体素的个数
int oc_node_divide_one_leaf(OC_Node* ocn);

//八叉树中的顶点个数
int oc_node_vertices_size(OC_Node* ocn);


//含有顶点的体素
int oc_node_voxel_size(OC_Node*ocn);

//八叉树的叶子结点个数
int oc_node_leaf_size(OC_Node*ocn);

//输入八叉树的某个结点，对它以下的所有含有顶点的体素划分
//返回划分后含有顶点的体素个数
int oc_node_divide_all_leaves(OC_Node* ocn);


//int oc_tree_size();
#ifdef __cplusplus
}
#endif
#endif