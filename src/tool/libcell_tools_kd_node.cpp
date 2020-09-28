#include<tool/libcell_tools_kd_node.h>


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