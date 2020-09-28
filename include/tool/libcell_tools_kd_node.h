#ifndef LIBCELL_KD_TREE_H_
#define LIBCELL_KD_TREE_H_

#include<Mesh_IO/Mesh_IO.h>
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
void KD_Node_init(KD_Node* kn)
{
    kn->left=NULL;
    kn->right=NULL;
    kn->parent=NULL;
    kn->dim=0;
    //Vertex_init(kn->value);
    kn->value=NULL;
    kn->prop=NULL;
}


KD_Node* create_kd_tree(template_v** values,int size,int dim);

/*

double two_vertex_distance(Vertex* v1,Vertex *v2)
{
    if(v1==NULL||v2==NULL)
    {
        printf("cuowu\n");
    }
    double re=0;
    for(int i=0;i<v1->point_size;i++)
    {   
        re+=(v1->point[i]-v2->point[i])*(v1->point[i]-v2->point[i]);
    } 
    re=SAFE_SQRT(re);
    return re; 
}
*/
/*
Node* kd_tree_find_nearest_sphere(Vertex* mu,double r,KD_Node* kn)
{
    if(kn==NULL)
    {
        return NULL;
    }
    Node* n=NULL; 
    KD_Node* knn=kn;
    int d=0,point_size=mu->point_size;
    while(knn!=NULL)
    {
       // printf("id:%d\n",knn->value->id);
        if(knn!=kn)
        {
            n=node_overlying(n,knn);
        }
        if(knn->value==NULL)
        {
            printf("节点的value为NULL\n");
            return NULL;
        }
        d=knn->dim%point_size;
        if(mu->point[d]>=knn->value->point[d])
        {
            knn=knn->right;
        }  
        else
        {
            knn=knn->left;
        }
    }
    int flag=0;
    Node* re=NULL;
   
    if(n==NULL)
    {
        if(two_vertex_distance(mu,kn->value)<r)
        {
            re=node_overlying(re,kn->value);
        }    
        return re;
    }
    Node* nit=n;
    KD_Node* kn1=NULL;
    kn1=(KD_Node*)(nit->value);
    if(two_vertex_distance(mu,kn1->value)<r)
    {
        flag=1;
        re=node_overlying(re,kn1->value);
    }

   // nit=(Node*)(nit->Next);
    Node* n1=NULL,*n2=NULL;
    for(;nit!=NULL;nit=(Node*)(nit->Next))
    {
        kn1=(KD_Node*)(nit->value);
        kn1=KD_BROTHER(kn1);
        n1=NULL;
        if(kn1==NULL)
        {
            continue;
        }
        if(flag==0)
        {
            n1=kd_tree_find_nearest_sphere(mu,r,kn1);
                

        }
        else if(flag==1)
        {
            if(KD_IS_INTERSECT_WITH_SPHERE(mu,r,kn1->parent)==1)
            {
                n1=kd_tree_find_nearest_sphere(mu,r,kn1);
            }

        }
        if(n1!=NULL)
        {
            flag=1;
            n2=node_reverse(re);
            if(re==NULL)
            {
                re=n1;
            }
            else
            {
                n2->Next=n1;
                n1->Prev=n2;
            } 
        }
    }
    free_node(n);
    return re;

}*/
/*
static int is_existence_in_vs(Vertex* v,Vertex** vs,int size)
{
    int n=size-1;
    while(n>=0)
    {
        if(v->id==vs[n]->id)
        {
            return 1;
        }
        n--;
    }
    return 0;
}*/
/*
//寻找排除vs点集且离mu最近的点
KD_Node* kd_tree_find_nearest_and_except_vs(Vertex* mu,KD_Node* kn,Vertex** vs,int size)
{
    if(kn==NULL)
    {
        return NULL;
    }
    Node* n=NULL; 
    KD_Node* knn=kn;
    int d=0,point_size=mu->point_size;
    while(knn!=NULL)
    {
       // printf("id:%d\n",knn->value->id);
        if(knn!=kn)
        {
            n=node_overlying(n,knn);
        }
        if(knn->value==NULL)
        {
            printf("节点的value为NULL\n");
            return NULL;
        }
        d=knn->dim%point_size;
        if(mu->point[d]>=knn->value->point[d])
        {
            knn=knn->right;
        }  
        else
        {
            knn=knn->left;
        }
    }
    if(n==NULL)
    {
        return kn;
    }
    double distance1=0,distance2=0;
    KD_Node*re=NULL,*kn1=NULL,*kn2=NULL;
    Node* nit=n;
    re=(KD_Node*)(nit->value);
    if(is_existence_in_vs(re->value,vs,size)==1)
    {
        distance1=-1;
    }
    else
    {
    
       distance1=two_vertex_distance(re->value,mu);

    }
    //nit=(Node*)(nit->Next);
    for(;nit!=NULL;nit=(Node*)(nit->Next))
    {
        kn1=(KD_Node*)(nit->value);
        kn1=KD_BROTHER(kn1);
        if(kn1==NULL)
        {
            continue;
        }
        if(KD_IS_INTERSECT_WITH_SPHERE(mu,distance1,kn1->parent)==1)
        {
            kn2=kd_tree_find_nearest_and_except_vs(mu,kn1,vs,size);
            if(is_existence_in_vs(kn2->value,vs,size)==1)
            {
                distance2=-1;
            }
            else
            {
                distance2=two_vertex_distance(mu,kn2->value);
            }
            if(RELATIONSHIP_PARTIAL_ORDER(distance1,distance2)==1)
            {
                distance1=distance2;
                re=kn2;
            }
        }


    }
    free_node(n);
   // printf("re:%lf\n",distance1);  
    return re;

}
*/
//度量空间负数表示无穷
//寻找距离mu最近点
/*
KD_Node* kd_tree_find_nearest(Vertex* mu,KD_Node* kn)
{
    if(kn==NULL)
    {
        return NULL;
    }
    Node* n=NULL; 
    KD_Node* knn=kn;
    int d=0,point_size=mu->point_size;
    while(knn!=NULL)
    {
       // printf("id:%d\n",knn->value->id);
        if(knn!=kn)
        {
            n=node_overlying(n,knn);
        }
        if(knn->value==NULL)
        {
            printf("节点的value为NULL\n");
            return NULL;
        }
        d=knn->dim%point_size;
        if(mu->point[d]>=knn->value->point[d])
        {
            knn=knn->right;
        }  
        else
        {
            knn=knn->left;
        }
    }
    if(n==NULL)
    {
        return kn;
    }
    double distance1=0,distance2=0;
    KD_Node*re=NULL,*kn1=NULL,*kn2=NULL;
    Node* nit=n;
    re=(KD_Node*)(nit->value);
    distance1=two_vertex_distance(re->value,mu);
  //  nit=(Node*)(nit->Next);
    for(;nit!=NULL;nit=(Node*)(nit->Next))
    {
        kn1=(KD_Node*)(nit->value);
        kn1=KD_BROTHER(kn1);
        if(kn1==NULL)
        {
            continue;
        }
        if(KD_IS_INTERSECT_WITH_SPHERE(mu,distance1,kn1->parent)==1)
        {
            kn2=kd_tree_find_nearest(mu,kn1);
            distance2=two_vertex_distance(mu,kn2->value);
            if(RELATIONSHIP_PARTIAL_ORDER(distance1,distance2)==1)
            {
                distance1=distance2;
                re=kn2;
            }
        }


    }
    free_node(n);
   // printf("re:%lf\n",distance1);  
    return re;

}*/
/*
void free_kdnode(KD_Node* kn)
{
    if(kn==NULL)
    {
        return;
    }
    free_kdnode(kn->left);
    free_kdnode(kn->right);
    free(kn);
}*/
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
