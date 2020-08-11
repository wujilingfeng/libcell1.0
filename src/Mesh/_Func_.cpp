#include<Mesh/_Func_.h>

template_v Mesh_get_vertex(struct Mesh* own,int id)
{
	template_v v;
    Vertex_init_(&v);
    auto it=own->vertices.find(id);
    if(it!=own->vertices.end())
    {
        v=*(it->second);
//return *(it->second);
    }
    return v;
//return NULL;
}
template_c Mesh_get_cell(struct Mesh* own,int id)
{
	template_c c;
    Cell_init_(&c);
    auto it=own->cells.find(id);
    if(it!=own->cells.end())
    {
        c= *(it->second);
    }
    return c;
//return NULL;

}
template_f Mesh_get_face(struct Mesh*own,int id)
{
    template_f f;
    Face_init_(&f);
    auto it=own->faces.find(id);
    if(it!=own->faces.end())
    {
        f=*(it->second);
    }
    return f;
}
template_f* Mesh_get_facev(struct Mesh*own,template_v** temp_v,int size)
{
    Node* node=temp_v[0]->faces;
    bool b=false;
    while(node!=NULL)
    {   b=true;
        for(int i=1;i<size;i++)
        {
            if(node_find(temp_v[i]->faces,node->value)==NULL)
            {
                b=false;
                break;
            }

        }
        if(b)
        {

            break;
        }
        node=(Node*)(node->Next);
    }

    if(b)
    {
        return (template_f*)(node->value);
    }

    return NULL;
}
template_c* Mesh_get_cellv(struct Mesh*own,template_v** temp_v,int size)
{
    Node* node=temp_v[0]->cells;
    bool b=false;
    while(node!=NULL)
    {
        b=true;
        for(int i=1;i<size;i++)
        {
            if(node_find(temp_v[i]->cells,node->value)==NULL)
            {
                b=false;
                break;
            }
        
        }
        if(b)
        {
            break;
        
        }
        node=(Node*)node->Next;
    
    }
    if(b)
    {
        return (template_c*)(node->value);
    
    }
    return NULL;

}
template_c* Mesh_get_cellf(struct Mesh* own,template_f** temp_f,int size)
{
    bool b=false;
    int i=0;
    for(i=0;i<2;i++)
    {   if(temp_f[0]->halffaces[i]->cell==NULL)
        {
            //b=false;
            continue;
        
        }
        b=true;
        for(int j=1;j<size;j++)
        {
            if(temp_f[j]->halffaces[0]->cell==temp_f[0]->halffaces[i]->cell||temp_f[j]->halffaces[1]->cell==temp_f[0]->halffaces[i]->cell)
            {
                continue;
            
            
            }
            else
            {
                b=false;
                break;
            }
        
        
        }
        if(b)
        {
        break;
        }
    
    
    }
    if(b)
    {
        return (template_c*)temp_f[0]->halffaces[i]->cell;
    }
    return NULL;


}
int Mesh_num_c(struct Mesh* own)
{

    return own->cells.size();
}
int Mesh_num_f(struct Mesh* own)
{

    return own->faces.size();
}
int Mesh_num_v(struct Mesh* own)
{
	return own->vertices.size();
}
int Mesh_num_hf(struct Mesh* own)
{

    return own->halffaces.size();
}
template_v* Mesh_create_vertex(struct Mesh* own)
{
    template_v* v=(template_v*)malloc(sizeof(template_v));

    Vertex_init_(v);
    v->id=own->vertex_id++;
    own->vertices[v->id]=v;
    if(own->init_v_prop!=NULL)
    {
        own->init_v_prop(v);
    }
    return v;
}
template_v* Mesh_create_vertexv(struct Mesh*own,double* p,int size)
{

    template_v*v= Mesh_create_vertex(own);
    v->point_size=size;
   v->point=(double*)malloc(sizeof(double)*size);
    for(int i=0;i<size;i++)
    {
        v->point[i]=p[i];
    
    }
    return v;
}
template_f * Mesh_create_face(struct Mesh* own)
{
    template_f* f=(template_f*)malloc(sizeof(template_f));
    Face_init_(f);
    f->halffaces[0]->cell=NULL;f->halffaces[1]->cell=NULL;
    f->id=own->face_id++;
    own->faces[f->id]=f;
    if(own->init_f_prop!=NULL)
    {
        own->init_f_prop(f);
    }
    return f;

}
template_c* Mesh_create_cell(struct Mesh* own)
{
    template_c *c=(template_c*)malloc(sizeof(template_c));

    Cell_init_(c);
    c->id=own->cell_id++;
    own->cells[c->id]=c;
    if(own->init_c_prop!=NULL)
    {
        own->init_c_prop(c);
    }

    return c;

}
//如果返回NULL,意味创建失败
template_hf * Mesh_create_halfface(struct Mesh* own,template_f* f,template_v** temp_v,int size)
{
    if(f==NULL)
    {
        return 0;
    }
    template_hf *hf=NULL;
    if(f->halffaces[0]->vertices==NULL)
    {
        hf=f->halffaces[0];
    }
    else if(f->halffaces[1]->vertices==NULL)
    {

        hf=f->halffaces[1];
    }
    else
    {
        if(f->halffaces[0]->cell==NULL)
        {
            hf=f->halffaces[0];
        }
        else if(f->halffaces[1]->cell==NULL)
        {
            hf=f->halffaces[1];
        }
        else
        {

            printf("流形的网格拓扑错误\n"); 
            for(int i=0;i<size;i++)
            {
                printf("%d ",temp_v[i]->id);
            }
            printf("\n");
            return NULL;
        }
    }
    hf->id=own->halfface_id++;
    own->halffaces[hf->id]=hf;
    if(hf->vertices==NULL)
    {
        hf->vertices=(template_v**)malloc(sizeof(template_v*)*size);
    }
    for(int i=0;i<size;i++)
    {
        hf->vertices[i]=temp_v[i];
    }
    hf->vertices_size=size;
    if(own->init_hf_prop!=NULL)
    {
        own->init_hf_prop(hf);
    }

    return hf;
}
//以下程序表明点到face的遍历没有不能保证顺(逆)时针(当cell是曲面)
template_f* Mesh_create_facev(struct Mesh*own,template_v** temp_v,int size)
{
	template_f*f=Mesh_get_facev(own,temp_v,size);
	if(f==NULL)
	{
	    f=Mesh_create_face(own);
	}
	else
	{
	    return f;
	}
	f->vertices_size=size;
	f->vertices=(template_v**)malloc(sizeof(template_v*)*size);
    for(int i=0;i<size;i++)
    {

        f->vertices[i]=temp_v[i];
        temp_v[i]->faces=node_overlying(temp_v[i]->faces,(void*)f);

    }
    return f;
}
//要判断是否已经创建相应的cell
template_c* Mesh_create_cellf(struct Mesh* own,template_hf**temp_hf,int size)
{
    template_c* c=Mesh_create_cell(own);
    
    Node* node=NULL,*node_v=NULL;
    for(int i=0;i<size;i++)
    {
        temp_hf[i]->cell=c;
        node=node_overlying(node,(void*)temp_hf[i]);
        for(int j=0;j<temp_hf[i]->vertices_size;j++)
        {
            if(node_find(node_v,temp_hf[i]->vertices[j])==NULL)
            {
                node_v=node_overlying(node_v,temp_hf[i]->vertices[j]); 
            
            } 
        } 
    }
    c->halffaces=node;

    Node* temp_nodev=node_v;
    node_v=node_reverse(node_v);
    int le=node_size(node_v);
    c->vertices=(template_v**)malloc(sizeof(template_v*)*le);
    c->vertices_size=le;
    int i=0;
    while(node_v!=NULL)
    {       c->vertices[i]=(template_v*)(node_v->value);

        
        ((template_v*)(node_v->value))->cells=node_overlying(((template_v*)(node_v->value))->cells,(void*)c);
        node_v=(Node*)(node_v->Prev);
        i++;
    }
    free_node(temp_nodev);
    return c;
}
//应该写 Mesh_get_cellv，并用在下面判断(但速度会变慢)
//以下程序也不能保证点到cell的遍历是顺(逆)时针(当cell表示面时)
//创建单形使用的函数
//如果返回NULL意味创建失败
template_c* Mesh_create_cellv(struct Mesh* own,template_v** temp_v,int size)
{
    template_hf** hfs=(template_hf**)malloc(sizeof(template_hf*)*size);

    template_v** temp_v1=(template_v**)malloc(sizeof(template_v*)*(size-1));
    for(int i=0;i<size;i++)
    {
        int temp_i=0;
        for(int j=0;j<size;j++)
        {
            if(j!=i)
            {

                temp_v1[temp_i]=temp_v[j];
                temp_i++;
            }
        }
        if((size-i-1)%2==0)
        {}
        else if(size>2)
        {
            template_v* vv=temp_v1[0];
            temp_v1[0]=temp_v1[1];
            temp_v1[1]=vv;
        }
        template_f* f=Mesh_create_facev(own,temp_v1,size-1);
        hfs[i]=Mesh_create_halfface(own,f,temp_v1,size-1);
        if(hfs[i]==NULL)
        {
           free(temp_v1);free(hfs);
           return NULL; 
        }
    }


    template_c* c=Mesh_create_cell(own);
    c->vertices=(template_v**)malloc(sizeof(template_v*)*size);
    c->vertices_size=size;


    for(int i=0;i<size;i++)
    {
        c->vertices[i]=temp_v[i];

        temp_v[i]->cells=node_overlying(temp_v[i]->cells,(void*)c);
//#ifdef SIMPLEX_REQUIRE   
     
    

        hfs[i]->cell=c;
        c->halffaces=node_overlying(c->halffaces,(void*)(hfs[i]));
//#endif
    }
    free(hfs);
    free(temp_v1);
    return c;
}
bool Mesh_delete_vertex(struct Mesh*own,const template_v &v_,bool burning)
{
    template_v*v;
    auto it=own->vertices.find(v_.id);
    if(it!=own->vertices.end())
    {
        own->vertices.erase(it);
        v=it->second;
    }
    else
    {
        return false;
    }

    template_c*c;
	Node* node=node_copy(v->cells),*node1,*node2;
	node1=node;
	while(node1!=NULL)
	{   
        node2=(Node*)(node1->Next);
	    c=(template_c*)(node1->value);
	    Mesh_delete_cell(own,*c,true);
	    node1=node2;
	
    }
	free_node(node);


#ifndef MANIFOLD_REQUIRE
	node=node_copy(v->faces);
	template_f*f;
	node1=node;
	while(node1!=NULL)
    {   node2=(Node*)(node1->Next);
		f=(template_f*)(node1->value);
        Mesh_delete_face(own,*f,true);
	    node1=node2;
    }
    free_node(node);

#endif
	if(own->free_v_prop!=NULL)
    {
        own->free_v_prop(v);
    }

    if(burning)
    {
 	    free_Vertex(v);       
	}

    return true;
}
bool Mesh_delete_halfface(struct Mesh*own,const template_hf& f_,bool burning)
{
    auto it=own->halffaces.find(f_.id);
    if(it!=own->halffaces.end())
    {
        if(it->second->cell!=NULL)
        {

            Mesh_delete_cell(own,*((template_c*)(it->second->cell)),burning);

            it->second->cell=NULL;

        }
#ifndef MANIFOLD_REQUIRE

        if(it->second->vertices!=NULL)
        {
            free(it->second->vertices);
            it->second->vertices=NULL;
        }
        it->second->vertices_size=0;
        own->halffaces.erase(it);
        if(own->free_hf_prop!=NULL)
        {
            own->free_hf_prop(it->second);
        }
#endif
    }
    else
    {
        return false;
    }
    return true;
}
bool Mesh_delete_face(struct Mesh*own,const template_f &f_,bool burning)
{
    template_f*f;
    auto it=own->faces.find(f_.id);
    if(it!=own->faces.end())
    {
        f=it->second;
        template_hf hf1=*(f->halffaces[0]),hf2=*(f->halffaces[1]);
        Mesh_delete_halfface(own,hf1,burning);

        Mesh_delete_halfface(own,hf2,burning);
#ifndef MANIFOLD_REQUIRE

        own->faces.erase(it);
#endif
    }
    else
    {
        return false;
    }
   /* if(f->halffaces[0]->cell!=NULL)
    {
	    Mesh_delete_cell(own,*((template_c*)f->halffaces[0]->cell),burning);
//Cell_remove_f((template_c*)f->cells[0],f);
    }
    if(f->halffaces[1]->cell!=NULL)
    {
//Cell_remove_f((template_c*)f->cells[1],f);
        Mesh_delete_cell(own,*((template_c*)f->halffaces[1]->cell),burning);
    }*/
#ifndef MANIFOLD_REQUIRE
    for(int i=0;i<f->vertices_size;i++)
    {
        Vertex_remove_f((template_v*)f->vertices[i],f);

    }
    if(own->free_f_prop!=NULL)
    {
        own->free_f_prop(f);
    }
    if(burning)
    {
        free_Face(f);
    }
#endif
    return true;
}
bool Mesh_delete_cell(struct Mesh* own,const template_c &c_,bool burning)
{
    template_c* c;
    auto it=own->cells.find(c_.id);

    if(it!=own->cells.end())
    {

        c=it->second;
        own->cells.erase(it);
    }
    else
    {
        return false;
    }
    //handle with face
    template_f*f;
    Node* node=c->halffaces,*node1;
    while(node!=NULL)
    {
        node1=(Node*)(node->Next);
      //  f=(template_f*)(((template_hf*)node->value)->face);
        //face要甩掉cell
        //Face_remove_c(f,c);
        Halfface_remove_c((template_hf*)(node->value),c);

        auto hfit=own->halffaces.find(((template_hf*)(node->value))->id);
        if(hfit!=own->halffaces.end())
        {
            own->halffaces.erase(hfit);
        }
        node=node1;
    }
    template_v*v;
    for(int i=0;i<c->vertices_size;i++)
    {
        v=(template_v*)c->vertices[i];
        Vertex_remove_c(v,c);
    }
#ifdef MANIFOLD_REQUIRE
    Node* node_f=c->halffaces;
    while(node_f!=NULL)
    {
        f=(template_f*)((template_hf*)(node_f->value))->face;
        if(f->halffaces[0]->cell==NULL&&f->halffaces[1]->cell==NULL)
        {
            for(int i=0;i<f->vertices_size;i++)
            {
                Vertex_remove_f((template_v*)f->vertices[i],f);

            }
  
            auto hfit=own->halffaces.find(f->halffaces[0]->id);
            if(hfit!=own->halffaces.end())
            {
                own->halffaces.erase(hfit);   
            }
            hfit=own->halffaces.find(f->halffaces[1]->id);
            if(hfit!=own->halffaces.end())
            {
                own->halffaces.erase(hfit);               
            }

  //          Mesh_delete_halfface(own,*(f->halffaces[0]),burning);
  //          Mesh_delete_halfface(own,*(f->halffaces[1]),burning);
            auto fit=own->faces.find(f->id);
            if(fit!=own->faces.end())
            {
                own->faces.erase(fit);
                if(own->free_f_prop!=NULL)
                {
                    own->free_f_prop(f);
                }
                if(burning)
                {
                    free_Face(f);
                }
            }

        }
        node_f=(Node*)(node_f->Next);

    }
    //如果要严格删除孤立点，则取消注释一下代码
/*    for(int i=0;i<c->vertices_size;i++)
    {
	    v=(template_v*)(c->vertices[i]);
        if(v->cells==NULL)
        {
	        auto vit=own->vertices.find(v->id);
	        if(vit!=own->vertices.end())
	        {
	            own->vertices.erase(vit);
	        }
	        free_Vertex(v);
        }

    }
  */  
#endif
    if(own->free_c_prop!=NULL)
    {
        own->free_c_prop(c);
    }
    if(burning)
    {
        free_Cell(c);

    }
    return true;
}
bool Mesh_face_is_boundary(struct Mesh* own,template_f *f)
{

    return Mesh_nface_is_boundary(own,*f);
}
bool Mesh_nface_is_boundary(struct Mesh* own,template_f &f)
{
    if(f.halffaces[0]->cell==NULL||f.halffaces[1]->cell==NULL)
    {   return true;}
    return false;
}
bool Mesh_vertex_is_boundary(struct Mesh* own,template_v &v)
{
    Node* node=v.faces;
    while(node!=NULL)
    {
        if(Mesh_nface_is_boundary(own,*((template_f*)(node->value))))
        {   return true;}
        node=(Node*)(node->Next);
    }
    return false;
}
//给它一个初始面，他会返回boundary的node链.这种寻找边界的方式.
//返回半边
//用红黑树RB_Tree做标记，代替prop
Node* Mesh_node_of_boundary_face(struct Mesh* own,template_f *f_)
{
    RB_Tree* tree=(RB_Tree*)malloc(sizeof(RB_Tree));
    RB_Tree_init_int(tree);
    RB_int rbt,*rbt1=NULL;
    if(f_==NULL)
    {
        return NULL;
    }
    Node* node=NULL,*node3=NULL;
    template_f *f;

    template_hf* hf;
    auto it=own->faces.find(f_->id);
    if(it==own->faces.end())
    {   printf("没找到\n");
        return NULL;
    }
    else
    {
	//printf("it\r\n");
        f=it->second;
    }
    if(!Mesh_nface_is_boundary(own,*f_))
    {return NULL;}
    else
    {
        //printf("shi\n");
//printf("boundary\r\n");
        node=node_overlying(node,(void*)f);
        
        if(f->halffaces[1]->cell==NULL)
        {
            hf=f->halffaces[1];
            node3=node_overlying(node3,(void*)hf);
        }
        if(f->halffaces[0]->cell==NULL)
        {
            hf=f->halffaces[0];
            node3=node_overlying(node3,(void*)hf);
        }
//node1->value=(void*)f;
        
//node1->value=
    }
   // 深度优先的遍历
//   printf("shendu\n");
   //int* temp_faces=(int*)malloc(sizeof(int));
   //*temp_faces=1;
    template_v *v;
    Node* node2;
    while(node!=NULL)
    {   
        push:
        f=(template_f*)node->value;
//        printf("f : %d\n",f->id);
        if(f->halffaces[0]->cell!=NULL)
        {
            hf=f->halffaces[0];
        }
        else
        {
            hf=f->halffaces[1];
        }
        rbt.key=f->id;
        tree->insert(tree,&rbt);
        //f->prop=(void*)temp_faces;
        for(int i=0;i<hf->vertices_size;i++)
        {
            v=(template_v*)hf->vertices[i];
            node2=v->faces;

            while(node2!=NULL)
            {
                if(Mesh_nface_is_boundary(own,*((template_f*)(node2->value))))
                {
                    /*if(((template_f*)(node2->value))->prop==NULL)
                    {
                        node=node_overlying(node,node2->value);
                       void* temp_value=NULL;
                        if(((template_f*)node2->value)->halffaces[0]->cell==NULL)
                        {
                           temp_value=(void*)(((template_f*)node2->value)->halffaces[0]);
                            node3=node_overlying(node3,temp_value);
                        }
                        else if(((template_f*)node2->value)->halffaces[1]->cell==NULL)
                        {

                            temp_value=(void*)(((template_f*)node2->value)->halffaces[1]);
                             node3=node_overlying(node3,temp_value);
                        }
                       
                        goto push;
                    
                    }*/
                    rbt.key=((template_f*)(node2->value))->id;
                    rbt1=(RB_int*)tree->find(tree,&rbt);
                    if(rbt1==NULL)
                    {
                        node=node_overlying(node,node2->value);
                       void* temp_value=NULL;
                        if(((template_f*)node2->value)->halffaces[0]->cell==NULL)
                        {
                           temp_value=(void*)(((template_f*)node2->value)->halffaces[0]);
                            node3=node_overlying(node3,temp_value);
                        }
                        else if(((template_f*)node2->value)->halffaces[1]->cell==NULL)
                        {

                            temp_value=(void*)(((template_f*)node2->value)->halffaces[1]);
                             node3=node_overlying(node3,temp_value);
                        }
                       
                        goto push;
                    
                    }
                }
                node2=(Node*)node2->Next;

            }
        }
    
        Node *node1=node;
        node=(Node*)node->Next;
        
        free(node1);
        node1=NULL;
    }
    //free(temp_faces);
    RB_Tree_free(tree);
    return node3;
}
static Node* Mesh_one_dim_boundary(struct Mesh* own)
{
    Node* node=NULL;
    for(auto it=own->faces.begin();it!=own->faces.end();it++)
    {
        if(it->second->halffaces[0]->cell==NULL)
        {

            node=node_overlying(node,it->second->halffaces[0]);
        }
        else if(it->second->halffaces[1]->cell==NULL)
        {
            node=node_overlying(node,it->second->halffaces[1]);
        }

    }
    return node;
}
//
void Mesh_external_cell_init_(struct Mesh* own)
{
    
    free_node(own->external_cell.halffaces);
    own->external_cell.halffaces=NULL;
//own->external_cell=(template_c*)malloc(sizeof(template_c));
#ifdef MANIFOLD_REQUIRE
    if(own->dimension==1)
    {
          own->external_cell.halffaces= Mesh_one_dim_boundary(own);
//            printf("external_cell size:%d\n",node_size(own->external_cell.halffaces));

          return;
    }
#endif
    for(auto it=own->faces.begin();it!=own->faces.end();it++)
    {   if(Mesh_nface_is_boundary(own,*(it->second)))
        {

            own->external_cell.halffaces=Mesh_node_of_boundary_face(own,it->second);

            break;
        }
    }
    template_hf* hf1,*hf2;
    Node* node=own->external_cell.halffaces;
  //  printf("external_cell size:%d\n",node_size(node));
    if(node==NULL)
    {

        return;
    }
    //reset_f_prop(own);
    while(node!=NULL)
    {

        hf1=(template_hf*)(node->value);
        if(hf1->vertices!=NULL)
        {
		    free(hf1->vertices);
		    hf1->vertices=NULL;
	    }
	
         hf2=Mesh_s_opposite_halfface(hf1);
         hf1->vertices=(template_v**)malloc(sizeof(template_v*)*hf2->vertices_size);

         hf1->vertices_size=hf2->vertices_size;
         for(int i=0;i<hf2->vertices_size;i++)
         {
              hf1->vertices[i]=hf2->vertices[i];
//              printf("external cell h_f id:%d ",((template_v*)(hf1->vertices[i]))->id);
         }
//	printf("\n");
        if(hf1->vertices_size>=2)
        {
              template_v * temp=hf1->vertices[0];
              hf1->vertices[0]=hf1->vertices[1];
              hf1->vertices[1]=temp;
         }
         node=(Node*)(node->Next);
    }
}
template_hf Mesh_opposite_halfface(template_hf &hf)
{
    template_f*f=(template_f*)hf.face;
    if(f->halffaces[0]->cell==hf.cell)
    {
        return *(f->halffaces[1]);
    }
    else
    {
        return *(f->halffaces[0]);

    }
}
template_hf* Mesh_s_opposite_halfface(template_hf*hf)
{
    template_f*f=(template_f*)(hf->face);
     if(f->halffaces[0]==hf)
     {
        return f->halffaces[1];
     }
     else if(f->halffaces[1]==hf)
     {
        return f->halffaces[0];
     }
     return NULL;
}   
Node* Mesh_non_manifold_vertices(struct Mesh* own)
{
    Node* re=NULL;
    for(auto iter=own->vertices.begin();iter!=own->vertices.end();iter++)
    {
        if(iter->second->cells==NULL)
        {
            re=node_overlying(re,iter->second);
        
        } 
    }
    return re;
}
void reset_c_prop(struct Mesh* mesh)
{
    for(auto cit=mesh->cells.begin();cit!=mesh->cells.end();cit++)
    {
	
    cit->second->prop=NULL;

    }
}
void reset_v_prop(struct Mesh* mesh)
{
    for(auto vit=mesh->vertices.begin();vit!=mesh->vertices.end();vit++)
    {
        vit->second->prop=NULL;    
    
    }

}
void reset_hf_prop(struct Mesh* mesh)
{
    for(auto hfit=mesh->halffaces.begin();hfit!=mesh->halffaces.end();hfit++)
    {
        hfit->second->prop=NULL;
        
    }

}
void reset_f_prop(struct Mesh* mesh)
{
    for(auto fit=mesh->faces.begin();fit!=mesh->faces.end();fit++)
    {
    
        fit->second->prop=NULL;
    
    }

}
iterator_v Mesh_fv_begin(struct Mesh* own,const template_f& f_)
{

    iterator_v it;
    iterator_v_init(&it);
    it.value=f_.vertices;
    return it;

}

iterator_v Mesh_fv_end(struct Mesh* own,const template_f& f_)
{
    iterator_v it;
    iterator_v_init(&it);
    it.value=f_.vertices;
    it.i=f_.vertices_size;
    return it;

}
iterator_v Mesh_hfv_begin(struct Mesh* own,const template_hf &hf)
{
    iterator_v it;
    iterator_v_init(&it);
    it.value=hf.vertices;
    return it;

}
iterator_v Mesh_hfv_end(struct Mesh* own,const template_hf &hf)
{
    iterator_v it;
    iterator_v_init(&it);
    it.value=hf.vertices;
    it.i=hf.vertices_size;

    return it;
}
iterator_v Mesh_cv_begin(struct Mesh* own,const template_c& c)
{
    iterator_v it;
    iterator_v_init(&it);
    it.value=c.vertices;
    return it;
}
iterator_v Mesh_cv_end(struct Mesh* own,const template_c& c)
{
    iterator_v it;
    iterator_v_init(&it);
    it.value=c.vertices;
    it.i=c.vertices_size;
    return it;
}


iterator_f Mesh_vf_begin(struct Mesh* own,const template_v& v)
{
    iterator_f iff;
    iterator_f_init(&iff);
    if(v.faces==NULL)
    {
        return iff;
    }
    //可以通过node->traits记录起点增加速度
    if(v.faces->traits==NULL)
    {
        v.faces->traits=(void*)(node_reverse(v.faces));
    }
    iff.node=*((Node*)(v.faces->traits));
    return iff;

}
iterator_f Mesh_vf_end(struct Mesh* own,const template_v& v)
{
    iterator_f iff;
    iterator_f_init(&iff);
    return iff;

}
iterator_hf Mesh_chf_begin(struct Mesh* own,const template_c& c)
{

    iterator_hf iff;
    iterator_hf_init(&iff);
    if(c.halffaces==NULL)
    {
        return iff;
    }
    if(c.halffaces->traits==NULL)
    {
        c.halffaces->traits=(void*)(node_reverse(c.halffaces));
    }
    iff.node=*((Node*)(c.halffaces->traits));

    return iff;
}
iterator_hf Mesh_chf_end(struct Mesh* own,const template_c&c)
{
    iterator_hf iff;
    iterator_hf_init(&iff);
    return iff;

}
iterator_c Mesh_vc_begin(struct Mesh* own,const template_v&v)
{
    iterator_c it;
    iterator_c_init(&it);
    if(v.cells==NULL)
    {
        return it;
    }
    it.node=*(v.cells);
    return it;

}
iterator_c Mesh_vc_end(struct Mesh* own,const template_v&v)
{

    iterator_c it;
    iterator_c_init(&it);
    return it;

}
Node* Mesh_vv_begin(struct Mesh* own,const template_v&v)
{
    if(own->simplex!=1)
    {
        printf("is not simplex\n");
        return NULL;
    }
    Node* node=NULL;
    for(iterator_f fit=Mesh_vf_begin(own,v);fit!=Mesh_vf_end(own,v);fit++)
    {
        for(auto vit=Mesh_fv_begin(own,*fit);vit!=Mesh_fv_end(own,*fit);vit++)
        {
            if((*vit).id!=v.id)
            {

                if(node_find(node,(void*)quote(vit))==NULL)
                {
                    node=node_overlying(node,(void*)quote(vit));
                }
            }
        }
    }
    return node;


}
Node* Mesh_vv_end(struct Mesh* own,const template_v&)
{

    return NULL;
}
void Mesh_printself(struct Mesh* own)
{
    printf("****************************\n\n");
    for(auto it=own->vertices.begin();it!=own->vertices.end();it++)
    {
        printf("id :%d \n",it->second->id);
        for(int j=0;j<it->second->point_size;j++)
        {
            printf("%lf  ",it->second->point[j]); 
        }
        printf("\n");
    }
    printf("faces\n");
    for(auto it=own->faces.begin();it!=own->faces.end();it++)
    {
       printf("id:%d\n",it->second->id);
       if(it->second->prop!=NULL)
       {
            printf("prop ");
       }
       for(int j=0;j<it->second->vertices_size;j++)
       {
           printf("id :%d  ",((template_v*)(it->second->vertices[j]))->id);
       
       }
       printf("\n");
    }
    printf("cell\n");

    for(auto it=own->cells.begin();it!=own->cells.end();it++)
    {
        printf("cell id:%d\n",it->second->id);
        for(int j=0;j<it->second->vertices_size;j++)
        {
            printf("id %d  ",((template_v*)(it->second->vertices[j]))->id);
        }
        printf("\n");
    }
 printf("****************************\n\n");

}
Node* Mesh_intersection_two_faces(struct Mesh*own,template_f*f0,template_f*f1)
{
    Node * node=NULL;
    int flag=1;
    for(int i=0;i<f0->vertices_size;i++)
    {
        flag=1;
        for(int j=0;j<f1->vertices_size;j++)
        {
            if((f0->vertices[i])==f1->vertices[j])
            {
                flag=0;
                break;
            } 
        }
        if(flag==0)
        {
            Node* node=node_overlying(node,f0->vertices[i]);
        }
    }
    return node;
}
void Mesh_free(struct Mesh* own)
{
    if(own->external_cell.halffaces!=0)
    {
//        printf("external\n");
        free_node(own->external_cell.halffaces);
        own->external_cell.halffaces=NULL;
    }

    for(auto iter=own->cells.begin();iter!=own->cells.end();iter++)
    {
//        printf("cells\n");
        free_Cell(iter->second);
    }
    own->cells.clear();
    for(auto iter=own->faces.begin();iter!=own->faces.end();iter++)
    {
//        printf("faces\n");
        free_Face(iter->second);
    }
    own->faces.clear();
    own->halffaces.clear();
    for(auto iter=own->vertices.begin();iter!=own->vertices.end();iter++)
    {
//        printf("vertices\n");
        free_Vertex(iter->second);
    }
    own->vertices.clear();
    /*for(auto iter=own->halffaces.begin();iter!=own->vertices.end();iter++)
    {
        
    }*/

}
static void default_free_v_prop(template_v*v)
{
    if(v->user_prop!=NULL)
    {
        free(v->user_prop);
    }
}
static void default_free_hf_prop(template_hf*hf)
{
    if(hf->user_prop!=NULL)
    {
        free(hf->user_prop);
    }
}
static void default_free_f_prop(template_f*f)
{
    if(f->user_prop!=NULL)
    {
        free(f->user_prop);
    }
}
static void default_free_c_prop(template_c*c)
{
    if(c->user_prop!=NULL)
    {
        free(c->user_prop);
    }
}

void Mesh_init(struct Mesh* own)
{
    Cell_init_(&(own->external_cell));
    own->cell_id=0;own->vertex_id=0;own->face_id=0;own->halfface_id=0;
    own->get_vertex=Mesh_get_vertex;
    own->get_cell=Mesh_get_cell;
    own->get_face=Mesh_get_face;
    own->get_facev=Mesh_get_facev;
    own->get_cellv=Mesh_get_cellv;
    own->get_cellf=Mesh_get_cellf;
    own->num_v=Mesh_num_v;
    own->num_c=Mesh_num_c;
    own->num_f=Mesh_num_f;
    own->num_hf=Mesh_num_hf;
    own->create_vertex=Mesh_create_vertex;
    own->create_vertexv=Mesh_create_vertexv;
    own->create_cell=Mesh_create_cell;
    own->create_face=Mesh_create_face;
    own->create_cellv=Mesh_create_cellv;
    own->create_cellf=Mesh_create_cellf;
    own->create_facev=Mesh_create_facev;
    own->create_halfface=Mesh_create_halfface;
    own->delete_cell=Mesh_delete_cell;
    own->delete_vertex=Mesh_delete_vertex;
    own->delete_face=Mesh_delete_face;
    own->vertex_is_boundary=Mesh_vertex_is_boundary;
    own->face_is_boundary=Mesh_nface_is_boundary;
    own->node_of_boundary_face=Mesh_node_of_boundary_face;
    own->external_cell_init_=Mesh_external_cell_init_;
    own->opposite_halfface=Mesh_opposite_halfface;
    own->s_opposite_halfface=Mesh_s_opposite_halfface;
    own->cv_begin=Mesh_cv_begin;
    own->cv_end=Mesh_cv_end;
    own->fv_begin=Mesh_fv_begin;
    own->fv_end=Mesh_fv_end;
    own->hfv_begin=Mesh_hfv_begin;
    own->hfv_end=Mesh_hfv_end;
    own->vc_begin=Mesh_vc_begin;
    own->vc_end=Mesh_vc_end;
    own->vf_begin=Mesh_vf_begin;
    own->vf_end=Mesh_vf_end;
    own->chf_begin=Mesh_chf_begin;
    own->chf_end=Mesh_chf_end;
    own->vv_begin=Mesh_vv_begin;
    own->vv_end=Mesh_vv_end;
    own->intersection_two_faces=Mesh_intersection_two_faces;
    own->non_manifold_vertices=Mesh_non_manifold_vertices;
    own->printself=Mesh_printself;
    own->init_v_prop=NULL;own->init_c_prop=NULL;own->init_f_prop=NULL;own->init_hf_prop=NULL;
    own->free_v_prop=default_free_v_prop;
    own->free_c_prop=default_free_c_prop;
    own->free_f_prop=default_free_f_prop;
    own->free_hf_prop=default_free_hf_prop;
  //  own->free_v_prop=NULL;
   // own->free_c_prop=NULL;
   // own->free_f_prop=NULL;
   // own->free_hf_prop=NULL;
    own->prop=NULL;
    MeshT_init(&(own->traits));
}

