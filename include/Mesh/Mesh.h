#ifndef LIB_CELL_MESH_H
#define LIB_CELL_MESH_H
#include<map>
#include<tools_rbtree.h>
#include "lib_cell_Iterator.h"
//度量空间的偏序关系宏，负数表示正无穷
//返回x>=y的关系
#define RELATIONSHIP_PARTIAL_ORDER(x,y) x<0?1:(y<0?0:(x>=y?1:-1)) 
typedef struct Mesh{

    template_c external_cell;
	std::map<int,template_c*> cells;
	std::map<int, template_f*> faces;
	std::map<int,template_v*>vertices;
    std::map<int,template_hf*> halffaces;

	int cell_id,vertex_id,face_id,halfface_id;

    template_v  (*get_vertex)(struct Mesh*,int);
    template_c (*get_cell)(struct Mesh*,int);
    template_f (*get_face)(struct Mesh*,int);
	template_f*(*get_facev)(struct Mesh*,template_v**,int);
	template_c*(*get_cellv)(struct Mesh*,template_v**,int);
	template_c* (*get_cellf)(struct Mesh*,template_f**,int);
	int (*num_v)(struct Mesh*);
	int (*num_c)(struct Mesh*);
	int (*num_f)(struct Mesh*);
    int (*num_hf)(struct Mesh*);
	template_v *(*create_vertex)(struct Mesh*);
    template_v* (*create_vertexv)(struct Mesh*,double*,int);
	template_c *(*create_cell)(struct Mesh*);
    template_f*(*create_face)(struct Mesh*);
    
    template_c *(*create_cellv)(struct Mesh*,template_v**,int);
    template_c *(*create_cellf)(struct Mesh*,template_hf**,int);
    template_f *(*create_facev)(struct Mesh*,template_v**,int);
    template_hf*(*create_halfface)(struct Mesh*,template_f*,template_v**,int);
//	void (*rearrange_id)(struct Mesh*)=mesh_rearrange_id;
    bool (*delete_cell)(struct Mesh*,const template_c&,bool);       
	bool (*delete_vertex)(struct Mesh*,const template_v&,bool); 
    bool (*delete_face)(struct Mesh*,const template_f&,bool);
    bool (*vertex_is_boundary)(struct Mesh*,template_v&);
    bool (*face_is_boundary)(struct Mesh*,template_f&);
    Node* (*node_of_boundary_face)(struct Mesh*,template_f*);
    void (*external_cell_init_)(struct Mesh*);
    template_hf (*opposite_halfface)(template_hf&);
	template_hf*(*s_opposite_halfface)(template_hf*);
	iterator_v (*cv_begin)(struct Mesh*,const template_c&);
    iterator_v (*fv_begin)(struct Mesh*,const template_f&);
    iterator_v (*fv_end)(struct Mesh*,const template_f&);
	iterator_v (*hfv_begin)(struct Mesh*,const template_hf&);
	iterator_v (*hfv_end)(struct Mesh*,const template_hf&);
	iterator_v (*cv_end)(struct Mesh*,const template_c&);
    iterator_c (*vc_begin)(struct Mesh*,const template_v&);
    iterator_c (*vc_end)(struct Mesh*,const template_v& );
	iterator_f (*vf_begin)(struct Mesh*,const template_v&);
	iterator_f (*vf_end)(struct Mesh*,const template_v&);
    iterator_hf (*chf_begin)(struct Mesh*,const template_c&);
    iterator_hf (*chf_end)(struct Mesh*,const template_c&);
//下面也要换成iteractor_v
    Node* (*vv_begin)(struct Mesh*,const template_v&);
	Node* (*vv_end)(struct Mesh*,const template_v&);
	
	void (*reset_c_prop)(struct Mesh*);
	void (*reset_v_prop)(struct Mesh*);
    void (*reset_f_prop)(struct Mesh*);
	void (*reset_hf_prop)(struct Mesh*);
	void (*init_v_prop)(template_v*);
	void (*init_c_prop)(template_c*);
	void (*init_f_prop)(template_f*);
	void (*init_hf_prop)(template_hf*);
	void (*free_v_prop)(template_v*);
	void (*free_c_prop)(template_c*);
	void (*free_f_prop)(template_f*);
	void (*free_hf_prop)(template_hf*);
	Node* (*intersection_two_faces)(struct Mesh*,template_f*,template_f*);
	Node* (*non_manifold_vertices)(struct Mesh*);
	void (*printself)(struct Mesh*);
	int simplex;
	void* prop;
	MeshT traits;
#ifdef MANIFOLD_REQUIRE
    int dimension;
#endif
        	
}Mesh;
#endif
