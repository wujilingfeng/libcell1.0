#ifndef LIBCELL_KD_TREE_H_
#define LIBCELL_KD_TREE_H_

#include<Mesh_IO/Mesh_IO.h>
#include "math.h"
#ifdef __cplusplus
extern "C"{
#endif
#ifndef SAFE_SQRT
#define SAFE_SQRT(x) x>=0?sqrt(x):sqrt(-x)
#endif
#ifndef SAFE_FREE
#define SAFE_FREE(x) if(x!=NULL){free(x);}
#endif
#define KD_BROTHER(x) x->parent==NULL?NULL:(x->parent->left==x?x->parent->right:x->parent->left)
#define KD_IS_INTERSECT_WITH_SPHERE(m,r,k) r<0?1:(fabs(k->value->point[k->dim]-m->point[k->dim])>=r?0:1)
//KD树用来储存反对称张量的指标 

typedef struct KD_Node{
    struct KD_Node *left;
    struct KD_Node *right;
    struct KD_Node* parent;
    int dim; 
    template_v* value;
    void * prop;
}KD_Node;

void KD_Node_init(KD_Node* kn);




KD_Node* create_kd_tree(template_v** values,int size,int dim);


double two_vertex_distance(Vertex* v1,Vertex *v2);
Node* kd_tree_find_nearest_sphere(Vertex* mu,double r,KD_Node* kn);

//寻找排除vs点集且离mu最近的点
KD_Node* kd_tree_find_nearest_and_except_vs(Vertex* mu,KD_Node* kn,Vertex** vs,int size);
//度量空间负数表示无穷
//寻找距离mu最近点

KD_Node* kd_tree_find_nearest(Vertex* mu,KD_Node* kn);
int kd_node_size(KD_Node* kn);
void free_kdnode(KD_Node* kn);
KD_Node* create_kd_tree_from_mesh(Mesh* m);
KD_Node* create_kd_tree_from_RBTree(std::map<int,template_v*> &mp);
/*KD_Node*  create_kd_tree(Mesh * m)
{
    
    Vertex** values=(Vertex**)malloc(sizeof(Vertex*)*m->num_v(m));
    int i=0;
    for(auto vit=m->vertices.begin();vit!=m->vertices.end();vit++)
    {

        values[i]=vit->second; 
        i++;
    } 
    
    KD_Node* kn=create_kd_node(values,i,0);
    free(values);
    
    return kn; 
}*/

#ifdef __cplusplus
}
#endif
#endif
