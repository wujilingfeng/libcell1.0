#ifndef LIB_CELL_ITERATOR_H
#define LIB_CELL_ITERATOR_H
#define iterator_v lib_cell_iterator_vertices
#define iterator_f lib_cell_iterator_faces
#define iterator_c lib_cell_iterator_cells
#define quote lib_cell_quote
#include "Mesh_Frame.h"
#include<iterator/cstruct_iterator.h>
#include<map>
typedef struct iterator_v{
template_v** value;
int i;
}iterator_v;
void iterator_v_init(iterator_v*);
iterator_v operator++(iterator_v&); 
iterator_v operator++(iterator_v&,int);
template_v* quote(iterator_v&);
template_v operator*(iterator_v& );
bool operator!=(const iterator_v&,const iterator_v&);
typedef struct iterator_f{
Node node;
}iterator_f;
void iterator_f_init(iterator_f*);
iterator_f operator++(iterator_f&);
iterator_f operator++(iterator_f&,int);
template_f* quote(iterator_f&);
template_f operator*(iterator_f&);
bool operator!=(const iterator_f&,const iterator_f&);
typedef struct iterator_c{
Node node;
}iterator_c;
void iterator_c_init(iterator_c*);
iterator_c operator++(iterator_c&);
iterator_c operator++(iterator_c&,int);
template_c operator*(iterator_c&);
template_c* quote(iterator_c&);
bool operator!=(const iterator_c&,const iterator_c&);
typedef struct iterator_hf
{Node node;}iterator_hf;
void iterator_hf_init(iterator_hf*);
iterator_hf operator++(iterator_hf&);
iterator_hf operator++(iterator_hf&,int);
template_hf* quote(iterator_hf&);
template_hf operator*(iterator_hf&);
bool operator!=(const iterator_hf&,const iterator_hf&);
//auto Mesh_vertices_begin(struct Mesh*);
//auto Mesh_veritces_end(struct Mesh*);

#endif
