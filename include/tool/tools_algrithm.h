//一些小算法
#include<Mesh/_Func_.h>
#include<Tensors/Tensors_Operation.h>
#include<tools/mpf_tools.h>
#ifdef MANIFOLD_REQUIRE
//face要保证是单形
Node* Mesh_adjust_nature_halffaces(Node*node)
{
    Node* re_n=NULL,*re_n1=NULL;
    if(node==NULL)
    {
    
        return NULL;
    }
    int *temp_int=(int*)malloc(sizeof(int));
template_hf* temp_hf=(template_hf*)(node->value);

    re_n=node_overlying(re_n,(void*)temp_hf);

    re_n1=re_n;
    
    //((template_f*)(temp_hf->face))->prop=(void*)temp_int;
    while(re_n1!=NULL)
    {
push:

        temp_hf=(template_hf*)(re_n1->value);

        temp_hf->prop=(void*)temp_int;
        for(int i=0;i<temp_hf->vertices_size;i++)
        {
        


            template_v*temp_v=(template_v*)(temp_hf->vertices[i]);
            for(Node n_it=*node;*n_it!=NULL;n_it++)
            {
                if(((template_hf*)(n_it.value))->prop==temp_int)
                {
                    continue;
                
                }
                template_hf* temp_hf1=(template_hf*)(n_it.value);
                for(int j=0;j<temp_hf1->vertices_size;j++)
                {
                    if(temp_hf1->vertices[j]==temp_v)
                    {
                        re_n=node_overlying(re_n,(void*)temp_hf1);

                        re_n1=re_n;
                      //  re_n1->traits=(void*)temp_int;
                        goto push;
                    
                
                
                    }
            
                }
             
            }   
    
        }
 
        re_n1=(Node*)(re_n1->Next); 
    }
    free(temp_int);
   

        return re_n;

}
template_v* Mesh_division_face(template_m*m,template_f*f,double *p)
{
    if(f==NULL)
    {
    
        return NULL;
    }
    int p_size,fv_size;
    fv_size=f->vertices_size;
    p_size=((template_v*)(f->vertices[0]))->point_size;

    template_v* v1=m->create_vertexv(m,p,p_size),*temp_v[2][fv_size];
    temp_v[0][0]=NULL;
    temp_v[1][0]=NULL;
    Node* node1[2];
    node1[0]=NULL;
    node1[1]=NULL;

    for(int l=0;l<2;l++)
    {
        if(f->halffaces[l]->cell==NULL)
        {
            continue;
        }
        for(int i=0;i<fv_size;i++)
        {
            temp_v[l][i]=(template_v*)(f->halffaces[l]->vertices[i]);
            printf("%d ",temp_v[l][i]->id);
    
        }
        node1[l]=node_copy(((template_c*)(f->halffaces[l]->cell))->halffaces);
        node1[l]=node_delete_value(node1[l],(void*)(f->halffaces[l]));
        
    }

    m->delete_face(m,*f,true);
    for(int l=0;l<2;l++)
    {
        if(temp_v[l][0]==NULL)
        {
        
            continue;
        }
        template_v* temp_v1[fv_size][fv_size];
        for(int i=0;i<fv_size;i++)
        {
            temp_v1[i][0]=v1;
            int k=1;
            for(int j=0;j<fv_size;j++)
            {
                if(j==i)
                {
                    continue;
                }
                temp_v1[i][k]=temp_v[l][j];
                k++;
            }
            if(i%2==1&&fv_size>=2)
            {
                template_v* temp_v3=temp_v1[i][0];
                temp_v1[i][0]=temp_v1[i][1];
                temp_v1[i][1]=temp_v3;
        
            
            }
        
        }
        for(int i=0;i<fv_size;i++)
        {
        
        
            template_f*temp_f=m->create_facev(m,temp_v1[i],fv_size);
            template_hf* temp_hf=m->create_halfface(m,temp_f,temp_v[i],fv_size);

            node1[l]=node_overlying(node1[l],(void*)temp_hf); 
        }
        Node *m_node;
         m_node=Mesh_adjust_nature_halffaces(node1[l]);
        //
      //  m_node=node1[l];
        int n_size=node_size(m_node);
        
        template_hf* temp_hf1[n_size];
        int i=0;
        for(Node n_it=*(m_node);*n_it!=NULL;n_it++)
        {
            temp_hf1[i]=(template_hf*)(n_it.value);
        
        i++;
        }
        printf("n_szie:%d\n",n_size);
        m->create_cellf(m,temp_hf1,n_size);
        m->simplex=0;


    }
    return v1;
    
}
//division_size是分割的几多变形(几多面体..)
//由于face是单形，又因为欧拉定理(亏格),所以由face组成的流形的边界的“边”数是face个数的函数
template_v * Mesh_division_cell(template_c* c,template_hf* hf,double *p,int p_size,int division_size)
{
    template_v* v0=NULL;
    if(c==NULL)
    {
        return v0;

    }

    



    return v0;

}

int Mesh_is_local_min_max(template_m*m,template_v* v)
{
    int min=1,max=1;
    Node *vvit=m->vv_begin(m,*v);
    for(Node vit=*vvit;*vit!=NULL;vit++)
    {
        if(*((double*)(v->prop))>*((double*)(((template_v*)(*vit))->prop)))
        {
            min=0;

        }
        else if(*((double*)(v->prop))<*((double*)(((template_v*)(*vit))->prop)))
        {
            max=0;
        }
    }
    free_node(vvit);
    if(min==1)
    {

        return 0;
    }
    else if(max==1)
    {
        return 1;
    }
    else
    {
        return 2;
    }
    return 2;

}


#endif
