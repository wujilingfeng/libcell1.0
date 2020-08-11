#include<Mesh/_Func_.h>
#include<string.h>
#ifndef LIBCELL_TOOLS_VIEW_H_
#define LIBCELL_TOOLS_VIEW_H_
#ifdef __cplusplus
extern "C"{
#endif
void get_data_from_1dim_cell(Mesh*m,float **v,unsigned int **e);
//二维流形可视化
void get_data_from_2dim_cell(Mesh* m,float**v,unsigned int**f);

void get_data_from_cell(Mesh* m,float**v,unsigned int**f);
#ifdef __cplusplus
}
#endif
#endif

