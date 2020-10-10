#include<tool/libcell_tools_kd_node.h>

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

static void sort_vertex(template_v** values,int size,int dim)
{
    if(size<=1)
    {
        return;
    }
    template_v* key=values[0];
    int low=0,high=size-1;
    int flag=1;
    template_v* temp=NULL;
    int d=dim%key->point_size;
    while(low<high)
    {
        if(flag==1)
        {
            if(values[high]->point[d]<key->point[d])
            {
                temp=values[high];
                values[high]=values[low];
                values[low]=temp;
                low++;
                flag=0;
            }
            else
            {
                high--;
            }
        }
        else
        {
            if(values[low]->point[d]>key->point[d])
            {
                temp=values[low];
                values[low]=values[high];
                values[high]=temp;
                high--;
                flag=1;
            }
            else
            {
                low++;
            }

        }
    
    }
    sort_vertex(values,low,d);
    sort_vertex(values+low+1,size-low-1,d); 

}

static KD_Node* create_kd_node(template_v** values,int size,int dim)
{
    if(size==0)
    {
        return NULL;
    }
    sort_vertex(values,size,dim);
 
    int wei=size/2;
    int d=dim%values[wei]->point_size;

    for(;(wei-1)>=0;wei--)
    {
        if(values[wei]->point[d]!=values[wei-1]->point[d])
        {
            break;
        }
    }
     if(wei==0)
    {
        wei++;
    } 
    KD_Node* re=(KD_Node*)malloc(sizeof(KD_Node));
    KD_Node_init(re);
    re->dim=d;

    
    re->value=values[wei]; 
      
    if(size==1)
    {
        //printf("here4 wei:%d\n",wei);
        re->value=values[0];
       // printf("here4 wei:%d\n",wei);
        return re;
    }
    
    re->left=create_kd_node(values,wei,dim+1);
    re->right=create_kd_node(values+wei,size-wei,dim+1);

    re->left->parent=re;
    re->right->parent=re;

     
    return re;
}
KD_Node* create_kd_tree(template_v** values,int size,int dim)
{
   KD_Node* re=create_kd_node(values,size,dim);
   return re;
}
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

}
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
}

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

}

int kd_node_size(KD_Node* kn)
{
    if(kn==NULL)
    {
        return 0;
    }
    if(kn->left==NULL&& kn->right==NULL)
    {
        return 1;
    }



    int re=0;
    re=re+kd_node_size(kn->left);
    re=re+kd_node_size(kn->right);

    return re; 

}
void free_kdnode(KD_Node* kn)
{
    if(kn==NULL)
    {
        return;
    }
    free_kdnode(kn->left);
    free_kdnode(kn->right);
    free(kn);
}



KD_Node*  create_kd_tree_from_mesh(Mesh * m)
{ 
    if(m->num_v(m)<=0)
    {
        return NULL;
    }
    template_v** values=(template_v**)malloc(sizeof(template_v*)*m->num_v(m));
    int i=0;
    for(auto vit=m->vertices.begin();vit!=m->vertices.end();vit++)
    {

        values[i]=vit->second; 
        i++;
    } 
    
    KD_Node* kn=create_kd_node(values,i,0);
    free(values);
    
    return kn; 
}
KD_Node* create_kd_tree_from_RBTree(std::map<int,template_v*> &mp)
{
    if(mp.size()<=0)
    {
        return NULL;
    }
    template_v** values=(template_v**)malloc(sizeof(template_v*)*mp.size());
    int i=0;
    for(auto vit=mp.begin();vit!=mp.end();vit++)
    {
        values[i]=vit->second;
        i++;
    }
    KD_Node* kn=create_kd_node(values,i,0);
    free(values);
    return kn;


}