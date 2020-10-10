#include <tool/libcell_tools_octree.h>
void OC_Node_init(OC_Node*ocn)
{
    ocn->parent=NULL;
    memset(ocn->children,0,sizeof(OC_Node*)*8);
    ocn->loc_min[0]=0;ocn->loc_min[1]=0;ocn->loc_min[2]=0;
    ocn->loc_max[0]=0;ocn->loc_max[1]=0;ocn->loc_max[2]=0;
    ocn->value=NULL;
    ocn->prop=NULL; 
}
void oc_node_free(OC_Node*ocn)
{
    if(ocn==NULL)
    {
        return ;
    } 
    for(int i=0;i<8;i++)
    {
        oc_node_free(ocn->children[i]);
    }
    
    free_node(ocn->value);
    free(ocn);
    return;
}

int oc_node_divide_one_leaf(OC_Node* ocn)
{
    if(ocn==NULL)
    {
        return 0;
    }
    if(ocn->children[0]!=NULL)
    {
        printf("oc_node_divide_one_leaf is not leaf\n");
        return 0;
    }
    double steps[3]={0,0,0};
    for(int i=0;i<3;i++)
    {
        steps[i]=(ocn->loc_max[i]-ocn->loc_min[i])/2.0;
    }
    if(steps[0]<=0||steps[1]<=0||steps[2]<=0)
    {
        printf("steps is zero\n");
        return 0;
    }
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            for(int k=0;k<2;k++)
            {
                ocn->children[i*4+j*2+k]=(OC_Node*)malloc(sizeof(OC_Node));
                OC_Node_init(ocn->children[i*4+j*2+k]);
                ocn->children[i*4+j*2+k]->parent=ocn;
                ocn->children[i*4+2*j+k]->loc_min[0]=ocn->loc_min[0]+steps[0]*i;
                ocn->children[i*4+2*j+k]->loc_min[1]=ocn->loc_min[1]+steps[1]*j; 
                ocn->children[i*4+2*j+k]->loc_min[2]=ocn->loc_min[2]+steps[2]*k;
                ocn->children[i*4+2*j+k]->loc_max[0]=ocn->children[i*4+2*j+k]->loc_min[0]+steps[0]; 
                ocn->children[i*4+2*j+k]->loc_max[1]=ocn->children[i*4+2*j+k]->loc_min[1]+steps[1]; 
                ocn->children[i*4+2*j+k]->loc_max[2]=ocn->children[i*4+2*j+k]->loc_min[2]+steps[2]; 
            }

        }
    }
    for(Node* nit=(Node*)(ocn->value);nit!=NULL;nit=(Node*)(nit->Next))
    {
        template_v* vv=(template_v*)(nit->value);
        int i=0,j=0,k=0;
        i=(int)((vv->point[0]-ocn->loc_min[0])/steps[0]);  
        j=(int)((vv->point[1]-ocn->loc_min[1])/steps[1]);
        k=(int)((vv->point[2]-ocn->loc_min[2])/steps[2]);
        if(i>=2)
        {
            i=1;
        }    
        if(j>=2)
        {
            j=1;
        }
        if(k>=2)
        {
            printf("k:%d\n",k);
            k=1;
        }
        Node* n=(Node*)(ocn->children[i*4+j*2+k]->value);
        n=node_overlying(n,vv);
        ocn->children[i*4+2*j+k]->value=n; 
    } 
    free_node((Node*)(ocn->value));
    ocn->value=NULL;

    int re=0;
    for(int i=0;i<8;i++)
    {
        if(ocn->children[i]->value!=NULL)
        {
            re++;
        }
    } 

    return re;
}


int oc_node_vertices_size(OC_Node* ocn)
{
    int re=0;
    if(ocn==NULL)
    {
        return re;
    }
    if(ocn->children[0]==NULL)
    {
        return node_size((Node*)(ocn->value)); 
    } 
    for(int i=0;i<8;i++)
    {
        re+=oc_node_vertices_size(ocn->children[i]);

    }
    return re;
}
int oc_node_voxel_size(OC_Node* ocn)
{
    int re=0;
    if(ocn==NULL)
    {
        return re;
    } 
    if(ocn->children[0]==NULL)
    {
        if(ocn->value!=NULL)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    for(int i=0;i<8;i++)
    {
        re+=oc_node_voxel_size(ocn->children[i]);
    }
    return re;
}

int oc_node_leaf_size(OC_Node* ocn)
{
    int re=0;
    if(ocn==NULL)
    {
        return re;
    } 

    if(ocn->children[0]==NULL)
    {
        return 1;
    }
    for(int i=0;i<8;i++)
    {
        re+=oc_node_leaf_size(ocn->children[i]);
    } 
    return re;

}
int oc_node_divide_all_leaves(OC_Node* ocn)
{
    int re=0;
    if(ocn==NULL)
    {
        return re;
    }
    if(ocn->children[0]==NULL)
    {
        int sum=node_size((Node*)(ocn->value));
        if(sum>1)
        {
            return  oc_node_divide_one_leaf(ocn);

        }
        else if(sum==1)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    for(int i=0;i<8;i++)
    {
        re+=oc_node_divide_all_leaves(ocn->children[i]);
    }
    return re;
}

