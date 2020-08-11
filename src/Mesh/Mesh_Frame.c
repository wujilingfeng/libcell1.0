#include<Mesh/Mesh_Frame.h>
void Vertex_init_(Vertex* pv)
{

    pv->cells=NULL;
    pv->faces=NULL;
    pv->id=-1;
    pv->point=NULL;
    pv->point_size=0;
    pv->prop=NULL;
    pv->user_prop=NULL;
	VertexT_init(&(pv->traits));
}
void Cell_init_(Cell *pv)
{
	pv->vertices=NULL;
	pv->halffaces=NULL;
	pv->id=-1;
    pv->vertices_size=0;
    pv->prop=NULL;
    pv->user_prop=NULL;
	CellT_init(&(pv->traits));
}

void Face_init_(Face *f)
{
	f->prop=NULL;
	f->user_prop=NULL;
    f->vertices=NULL;
	f->id=-1;
    f->halffaces[0]=(template_hf*)malloc(sizeof(template_hf));
    f->halffaces[1]=(template_hf*)malloc(sizeof(template_hf));
    HalfFace_init_(f->halffaces[0]);
    HalfFace_init_(f->halffaces[1]);
    f->halffaces[0]->face=f;
    f->halffaces[1]->face=f;
	 f->vertices_size=0;
	FaceT_init(&(f->traits));
 
}
void HalfFace_init_(HalfFace *hf)
{
    hf->vertices=NULL;
    hf->vertices_size=0;

    hf->cell=NULL;
    hf->face=NULL;
    hf->prop=NULL;
    hf->user_prop=NULL;
	hf->id=-1;
	HalfT_init(&(hf->traits));
}
void free_HalfFace(HalfFace *hf)
{
	if(hf->vertices!=NULL)
	{
        	free(hf->vertices);
            hf->vertices=NULL;
	}
	free(hf);
	hf=NULL;

}
void free_Vertex(Vertex*v)
{   
    if(v->point!=NULL)
    {
	
    		free(v->point);
	}
	if(v->faces!=NULL)
	{
		free_node(v->faces);
	}
//free_node(v->faces);
	if(v->cells!=NULL)
	{
		free_node(v->cells);
	}
    	free(v);
    	v=NULL;
}
void free_Face(Face*f)
{
	//printf("here1\n");
	if(f==NULL)
	{
        printf("cuowu\n");
        return;
    }
    if(f->vertices!=NULL)
    {
        free(f->vertices);
    }
    free_HalfFace(f->halffaces[0]);
    free_HalfFace(f->halffaces[1]);
    free(f);
    f=NULL;
}
void free_Cell(Cell*c)
{
    if(c->vertices!=NULL)
    {
        free(c->vertices);
    }
    if(c->halffaces!=NULL)
    {
        free_node(c->halffaces);
    }
    free(c);
    c=NULL;
}
void Halfface_remove_c(HalfFace*hf,Cell*c)
{
	if(hf==NULL||c==NULL)
	{return;}
	if(hf->cell==c)
	{
		hf->cell=NULL;
	}
}
void Face_remove_c(Face*f ,Cell*c)
{
    if(c==NULL)
    {
        return ;
    }

    if((f->halffaces[0])->cell==c)
    {
        ((template_hf*)f->halffaces[0])->cell=NULL;
    }
    else if(((template_hf*)f->halffaces[1])->cell==c)
    {
        ((template_hf*)f->halffaces[1])->cell=NULL;
    }
    else
    {
        return;
    }
}
void Vertex_remove_f(Vertex*v,Face*f)
{
    Node* node=node_find(v->faces,(void*)f);
    if(node==NULL)
    {
        return;
    }
    Node* node1=(Node*)(node->Prev),*node2=(Node*)(node->Next);
    if(node1!=NULL)
    {
        node1->Next=(void*)node2;
    }
    else
    {
        v->faces=node2;
    }
    if(node2!=NULL)
    {
        node2->Prev=(void*)node1;
    }
    free(node);
    node=NULL;
}

void Cell_remove_hf(Cell*c,HalfFace*f)
{
    Node*node=node_find(c->halffaces,(void*)f);
    if(node==NULL)
    {
        return;
    }
    Node* node1=(Node*)(node->Prev),*node2=(Node*)(node->Next);
    if(node1!=NULL)
    {
        node1->Next=(void*)node2;
    }
    else
    {
        c->halffaces=node2;
    }
    if(node2!=NULL)
    {
        node2->Prev=(void*)node1;
    }   
    free(node);
    node=NULL;
}
void Vertex_remove_c(Vertex*v,Cell*c)
{
    Node* node=node_find(v->cells,(void*)c);
    if(node==NULL)
    {
        return;
    }
    Node* node1=(Node*)(node->Prev),*node2=(Node*)(node->Next);
    if(node1!=NULL)
    {
        node1->Next=(void*)node->Next;
    }
    else
    {
        v->cells=node2;
    }
    if(node2!=NULL)
    {
        node2->Prev=(void*)node1;
    }
    free(node);
    node=NULL;

}


