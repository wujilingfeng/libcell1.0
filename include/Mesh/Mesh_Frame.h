#ifndef MESH_FRAME_H
#define MESH_FRAME_H
#include "Cell_Traits.h"
//给了一堆有顺序的点p0 p1 ..pn 那么定义的反对称张量是(p1-p0)/\(p2-p0)/\...(pn-p0)
#ifdef __cplusplus
extern "C"{
#endif
typedef struct Vertex{
    int id;
    double *point;
    int point_size;
   
    Node * cells;
    Node *faces;
    void *prop,*user_prop;
	VertexT traits;
}Vertex;

typedef struct HalfFace{
    int id;
    struct Vertex** vertices;
    int vertices_size;

    struct Cell* cell;
    struct Face* face;
    void *prop,*user_prop;
    HalfT traits;
}HalfFace;

typedef struct Face{
    int id;
    struct Vertex**vertices;
    int vertices_size;
    struct HalfFace* halffaces[2];
//int halffaces_size;

    void *prop,*user_prop;
    FaceT traits;
}Face;

typedef struct Cell{
    int id;

    struct Vertex**vertices;
    int vertices_size;
    Node* halffaces;
    void *prop,*user_prop;
    CellT traits;
}Cell;
void Vertex_init_(Vertex*);
void Face_init_(Face*);
void Cell_init_(Cell*);
void HalfFace_init_(HalfFace*);

void free_Vertex(Vertex*);
void free_Face(Face*);
void free_HalfFace(HalfFace*);
void free_Cell(Cell*);
void Halfface_remove_c(HalfFace*,Cell*);
void Face_remove_c(Face*,Cell*);
void Vertex_remove_f(Vertex*,Face*);
void Vertex_remove_c(Vertex*v,Cell*c);
void Cell_remove_f(Cell*,Face*);
#ifdef __cplusplus
}
#endif
#endif
