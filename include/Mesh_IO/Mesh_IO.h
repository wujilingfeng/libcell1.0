
#ifndef LIBCELL_MESH_IO_H
#define LIBCELL_MESH_IO_H
#include<stdio.h>
#include<Mesh/_Func_.h>
#include<string.h>
#define _ReadOff_   LibCell_ReadOff_
#define _ReadCell_  LibCell_ReadCell_
#define _WriteCell_ LibCell_WriteCell_
#define _ReadArray_ LibCell_ReadArray_
//当off文件描述单形是可以用_ReadOff_    需要在函数内部加判断给出提示
#ifdef __cplusplus
extern "C"{
#endif
//只支持三角形网格off文件的读取，tools_formats文件可以将off文件转为cell文件
void _ReadOff_(template_m* ,char const * ,int);
//float **用来储存纹理
Node* _ReadM_(template_m*,char const *);


void _ReadCell_(template_m*,char const *);
void _WriteCell_(template_m*,char const *);
//
void _ReadArray_(template_m*,double**,int**,int **,int*,int,int);
#ifdef __cplusplus
}
#endif
#endif

