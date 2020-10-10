#include<tool/libcell_kd_node_iterator.h>

void KD_Trav_init(KD_Trav*trav)
{
    trav->tree=NULL;
    trav->it=NULL;
    trav->prop=NULL;
}
KD_Trav kd_tree_begin(KD_Node* kn)
{
    KD_Node* it=kn;
    while(it->left!=NULL)
    {
        it=it->left;
    }
    KD_Trav trav;
    KD_Trav_init(&trav);
    trav.it=it;
    trav.tree=kn;

    return trav;
}

KD_Trav kd_tree_rbegin(KD_Node*kn)
{
    KD_Node* it=kn;
    while(it->right!=NULL)
    {
        it=it->right;
    }
    KD_Trav trav;
    KD_Trav_init(&trav);
    trav.it=it;
    trav.tree=kn;

    return trav;

}
KD_Trav kd_tree_end(KD_Node* kn)
{
    KD_Trav trav;
    KD_Trav_init(&trav);
    trav.tree=kn;
    return trav;
}


KD_Trav operator++(KD_Trav& trav)
{
    if(trav.it==NULL)
    {
        return trav;
    }
    KD_Node* it=trav.it;
    while(it->parent!=NULL&&it->parent->right==it)
    {
        it=it->parent;
    }
    //trav已经是终点
    if(it->parent==NULL)
    {
        trav.it=NULL; 
    }
    else
    {  
        it=it->parent->right;
        while(it->left!=NULL)
        {
            it=it->left;
        }
        trav.it=it;
    }
    return trav;

}
KD_Trav operator++(KD_Trav& trav,int)
{
    if(trav.it==NULL)
    {
        return trav;
    }
    KD_Trav trav1=trav;
    KD_Node* it=trav.it;
    while(it->parent!=NULL&&it->parent->right==it)
    {
        it=it->parent;
    }
    //trav已经是终点
    if(it->parent==NULL)
    {
        trav.it=NULL;
    }
    else
    {
        it=it->parent->right;
        while(it->left!=NULL)
        {
            it=it->left;
        }
        trav.it=it;
    }  

   return trav1;

}
bool operator!=(const KD_Trav&trav1,const KD_Trav&trav2)
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