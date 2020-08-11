#ifndef LIB_CELL_TRAITS_H
#define LIB_CELL_TRAITS_H

//表示我要处理的是单形cell(要弃用这个宏)
#define SIMPLEX_REQUIRE
//表示我要按照流形的标准来处理
#define MANIFOLD_REQUIRE
#define template_v Vertex
#define template_f Face
#define template_c Cell
#define template_hf HalfFace

#define template_m Mesh
#include<stdlib.h>
#include <tools_node.h>
#include<stdio.h>
#ifdef SIMPLEX_REQUIRE

#endif
#ifdef __cplusplus
extern "C"{
#endif
typedef struct VertexT{}VertexT;
typedef struct CellT{}CellT;
typedef struct MeshT{}MeshT;
typedef struct FaceT{}FaceT;
typedef struct HalfT{}HalfT;
void VertexT_init(VertexT* vt);
void CellT_init(CellT* ct);
void MeshT_init(MeshT*mt);
void FaceT_init(FaceT*ft);
void HalfT_init(HalfT*ht);
#ifdef __cplusplus
}
#endif
#endif
