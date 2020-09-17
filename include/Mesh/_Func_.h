
#ifndef LIB_CELL_FUNC_H
#define LIB_CELL_FUNC_H
#include "Mesh.h"

#ifdef __cplusplus
extern "C"{
#endif
template_v  Mesh_get_vertex(struct Mesh*,int);
template_c  Mesh_get_cell(struct Mesh*,int);
template_f Mesh_get_face(struct Mesh*,int);
template_f* Mesh_get_facev(struct Mesh*,template_v**,int);
template_c* Mesh_get_cellv(struct Mesh*,template_v**,int);
template_c* Mesh_get_cellf(struct Mesh*,template_f**,int);


int Mesh_num_v(struct Mesh*);
int Mesh_num_c(struct Mesh*);
int Mesh_num_f(struct Mesh*);
int Mesh_num_hf(struct Mesh*);

template_v * Mesh_create_vertex(struct Mesh*);
template_v* Mesh_create_vertexv(struct Mesh*,double*,int);
template_c * Mesh_create_cell(struct Mesh*);
template_f * Mesh_create_face(struct Mesh*);
template_c* Mesh_create_cellv(struct Mesh*,template_v**,int);
template_c* Mesh_create_cellf(struct Mesh*,template_hf**,int size);
template_f * Mesh_create_facev(struct Mesh*,template_v**,int);
template_hf* Mesh_create_halfface(struct Mesh*,template_f*,template_v**,int);
bool Mesh_delete_vertex(struct Mesh*,const template_v&,bool burning);
bool Mesh_delete_face(struct Mesh*,const template_f&,bool burning);
bool Mesh_delete_halfface(struct Mesh*,const template_hf&,bool burning);
bool Mesh_delete_cell(struct Mesh*,const template_c&,bool burning);
bool Mesh_vertex_is_boundary(struct Mesh*,template_v&);
bool Mesh_nface_is_boundary(struct Mesh* own,template_f &f);
bool Mesh_face_is_boundary(struct Mesh*,template_f*);
Node* Mesh_node_of_boundary_face(struct Mesh*,template_f*);
void Mesh_external_cell_init_(struct Mesh*);
template_hf Mesh_opposite_halfface(template_hf&);
template_hf* Mesh_s_opposite_halfface(template_hf*);
//void Mesh_cell_division(struct Mesh*,template_c*);
void reset_c_prop(struct Mesh*);
void reset_v_prop(struct Mesh*);
void reset_hf_prop(struct Mesh*);
void reset_f_prop(struct Mesh*);
iterator_v Mesh_fv_begin(struct Mesh*,const template_f& );
iterator_v Mesh_fv_end(struct Mesh*,const template_f&);
iterator_v Mesh_hfv_begin(struct Mesh*,const template_hf&);
iterator_v Mesh_hfv_end(struct Mesh*,const template_hf&);
iterator_v Mesh_cv_begin(struct Mesh*,const template_c&);
iterator_v Mesh_cv_end(struct Mesh*,const template_c&);

//iterator_f
iterator_f Mesh_vf_begin(struct Mesh*,const template_v&);
iterator_f Mesh_vf_end(struct Mesh*,const template_v&);
iterator_hf Mesh_chf_begin(struct Mesh*,const template_c&);
iterator_hf Mesh_chf_end(struct Mesh*,const template_c&);
//iterator_c

iterator_c Mesh_vc_begin(struct Mesh*,const template_v&);
iterator_c Mesh_vc_end(struct Mesh*,const template_v&);
Node* Mesh_vv_begin(struct Mesh*,const template_v&);
Node* Mesh_vv_end(struct Mesh*,const template_v&);
void Mesh_printself(struct Mesh*);
Node* Mesh_intersection_two_faces(struct Mesh*,template_f*,template_f*);
bool Mesh_two_cell_is_connected(struct Mesh* own,template_c* c1,template_c* c2);
Node* Mesh_non_manifold_vertices(struct Mesh*);
bool Mesh_is_manifold_vertices(struct Mesh*own,template_v*v);
void Mesh_free(struct Mesh*);
void Mesh_init(struct Mesh*);


#ifdef __cplusplus
}
#endif
//bool mesh_delete_vertex(struct Mesh*,template_v *);
//bool mesh_delete_cell(struct Mesh*,template_c *);
//void mesh_rearrange_id(struct Mesh*);
#endif
