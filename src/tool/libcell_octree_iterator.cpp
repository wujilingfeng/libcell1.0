#include<tool/libcell_octree_iterator.h>
void Octree_Trav_init(Octree_Trav* oct)
{
    oct->tree=NULL;
    oct->it=NULL;
    oct->prop=NULL;
}
Octree_Trav octree_begin(OC_Node* ocn)
{
    Octree_Trav oct;
    Octree_Trav_init(&oct);
    oct.tree=ocn;
    OC_Node* it=ocn;
    while(it->children[0]!=NULL)
    {
        it=it->children[0];
    } 
    oct.it=it;

    return oct;  
}
Octree_Trav octree_end(OC_Node*ocn)
{
    Octree_Trav oct;
    Octree_Trav_init(&oct);

    oct.tree=ocn;

    return  oct;
}
Octree_Trav operator++(Octree_Trav& trav)
{
    if(trav.it==NULL)
    {
        return trav;
    }
    OC_Node* it=trav.it;
    while(it->parent!=NULL&&it->parent->children[7]==it)
    {
        it=it->parent;
    } 
    if(it->parent==NULL)
    {
        trav.it=NULL;
    }
    else
    {
        for(int i=0;i<8;i++)
        {
            if(it->parent->children[i]==it)
            {
                it=it->parent->children[i+1];
                break;
            }
        }
        while(it->children[0]!=NULL)
        {
            it=it->children[0];
        } 
        trav.it=it;

    }
    return trav;

}

Octree_Trav operator++(Octree_Trav& trav,int)
{
    if(trav.it==NULL)
    {
        return trav;
    }
    Octree_Trav trav1=trav;
    OC_Node* it=trav.it;
    while(it->parent!=NULL&&it->parent->children[7]==it)
    {
        it=it->parent;
    } 
    if(it->parent==NULL)
    {
        trav.it=NULL;
    }
    else
    {
        for(int i=0;i<8;i++)
        {
            if(it->parent->children[i]==it)
            {
                it=it->parent->children[i+1];
                break;
            }
        }
        while(it->children[0]!=NULL)
        {
            it=it->children[0];
        } 
        trav.it=it;

    }
    return trav1;

}
bool operator!=(const Octree_Trav&trav1,const Octree_Trav&trav2)
{
    if(trav1.tree!=trav2.tree)
    {
        return true;
    }
    if(trav1.it!=trav2.it)
    {
        return true;
    }
    return false;
}

