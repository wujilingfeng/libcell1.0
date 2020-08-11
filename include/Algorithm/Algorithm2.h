//基于高精度的张量库算法,优点是概念上的统一（符合数学,张量天然稠密性）,算法简洁。高精度库保证精度不会出错
#ifndef LIB_CELL_ALGORITHM1_H
#define LIB_CELL_ALGORITHM1_H
#include<math.h>
#include<Mesh/_Func_.h>
//#include<Tensors/Antisymmetric_Tensor.h>
#include<Tensors/Tensors_Operation.h>
#include<tools/mpf_tools.h>
int libcell_factorial(int n)
{
    int result=1;
    for(int i=2;i<=n;i++)
    {
        result=result*i;
    }

    return result;
}
//把一个二维数组(矩阵)转化为反对称张量(它的模长也等价求矩阵广义行列式)
Tensor* Anti_tensor_mpf_from_v(Tensors_Algebra_System*tas,double**M,int rows,int cols)
{
    Tensor* re=tas->T_create();
    for(int j=0;j<cols;j++)
    {
        re->insert(tas->as,re,&j,1,tas->copy_from_double(M[0][j]));
    }
    for(int i=1;i<rows;i++)
    {   
        Tensor *t=tas->T_create();
        for(int j=0;j<cols;j++)
        {

            t->insert(tas->as,t,&j,1,tas->copy_from_double(M[i][j]));
        }
        Tensor*t1=Tensor_Wedge_(tas,re,t);
        tas->T_free(tas,t);tas->T_free(tas,re);
        re=t1;
    }
    return re;
}
//求一些点集对应的反对称张量
Tensor* Anti_tensor_mpf_from_point(Tensors_Algebra_System*tas,double**M,int rows,int cols)
{
    Tensor* re=tas->T_create();
    for(int j=0;j<cols;j++)
    {
        if((M[1][j]-M[0][j])==0)
        {
            continue;
        }
        //printf("d:%lf\n",M[1][j]-M[0][j]);
        re->insert(tas->as,re,&j,1,tas->copy_from_double(M[1][j]-M[0][j]));
    }
    for(int i=2;i<rows;i++)
    {   
        Tensor *t=tas->T_create();
        for(int j=0;j<cols;j++)
        {

            if((M[i][j]-M[0][j])==0)
            {
                continue;
            }
            //printf("d:%lf\n",M[i][j]-M[0][j]);
            t->insert(tas->as,t,&j,1,tas->copy_from_double(M[i][j]-M[0][j]));
        }
        Tensor*t1=Tensor_Wedge_(tas,re,t);
        tas->T_free(tas,t);
      //  tensor_mpf_print_self(re);
        tas->T_free(tas,re);
        re=t1;
    }
    return re;
}
//输入rows个点组成的单形，输出体积
__mpf_struct* area_simplex(Tensors_Algebra_System*tas,double**M,int rows,int cols)
{
    Tensor*t=Anti_tensor_mpf_from_point(tas,M,rows,cols);
    __mpf_struct* re=(__mpf_struct*)(tas->T_norm(tas,t));
    tas->T_free(tas,t);
 //   mpf_t re1;
  //  mpf_inits(re1,NULL);
  //  mpf_sqrt(re1,re); 
    mpf_sqrt(re,re);
    //int lr=mpf_cmp(re1,re);
    mpf_div_ui(re,re,libcell_factorial(rows-1));

    return re;
}
static bool libcell_is_same_dir(Tensors_Algebra_System*tas,Tensor*t1,Tensor*t2)
{
    __mpf_struct * re=(__mpf_struct*)(tas->T_inner_product(tas,t1,t2));
    //gmp_printf("re:%.Ff\n",re);
    int re1=mpf_cmp_ui(re,0);
    mpf_clear(re);
    free(re);
    if(re1>=0)
    {
        return true;
    }
  //  printf("shir\n");
    return false;

}
//增量凸包算法

//以给定的反对称张量作为正方向
//应该改为bool返回值的函数，用来接受是否成功
bool increasing_convex_hull(Tensors_Algebra_System*tas,Tensor*t,template_m*mesh,template_v*v)
{
    RB_Tree* tree=(RB_Tree*)malloc(sizeof(RB_Tree));
    RB_Tree_init_int(tree);
    RB_int rbt,*rbt1=NULL;
//    int rows=v->point_size+1,cols=v->point_size;
    /*
    if(mesh->simplex!=)
    {
        return;
    }
    */
   // printf("increasing once\n");
    int rows=t->order(t)+1,cols=v->point_size;
    //printf("rows:%d cols:%d\n",rows,cols);
    double**M=(double**)malloc(sizeof(double*)*rows);
    for(int i=0;i<rows;i++)
    {
        M[i]=(double*)malloc(sizeof(double)*cols);
    }
    for(int j=0;j<cols;j++)
    {
        M[rows-1][j]=v->point[j];
    }
    Node* node=NULL,*node2=NULL;
    //int *temp_int=(int*)malloc(sizeof(int)),*temp_int1=(int*)malloc(sizeof(int));
    //*temp_int=1;*temp_int1=0;
    for(auto cit=mesh->cells.begin();cit!=mesh->cells.end();cit++)
    {
        for(int i=0;i<rows-1;i++)
        {
            for(int j=0;j<cols;j++)
            {

                M[i][j]=((template_v*)(cit->second->vertices[i]))->point[j];
            }

        }
        Tensor* t1=Anti_tensor_mpf_from_point(tas,M,rows,cols);
     //   printf("cid:%d\n",cit->first);
      //  tensor_mpf_print_self(t1);
        if(!libcell_is_same_dir(tas,t,t1))
        {
        //    printf("shi\n");
            node=node_overlying(node,(void*)cit->second);
            //cit->second->prop=(void*)temp_int;
            rbt.key=cit->second->id;
            tree->insert(tree,&rbt);
            node2=node_overlying(node2,(void*)cit->second);
            tas->T_free(tas,t1);
            break;

        }
        else
        {

            //cit->second->prop=(void*)temp_int1;
            rbt.key=cit->second->id;
            tree->insert(tree,&rbt);
            tas->T_free(tas,t1);
        }

    } 
    template_c*c0=NULL,*c1=NULL;
    template_hf*hf1=NULL;
    //广度优先
    while(node!=NULL)
    {
        Node* temp_node=node_copy(node);
        free_node(node);
        node=NULL;
        for(Node *n_it=temp_node;n_it!=NULL;n_it=(Node*)(n_it->Next))
        {
            c0=(template_c*)(n_it->value);
            for(auto chf_it=mesh->chf_begin(mesh,*c0);chf_it!=mesh->chf_end(mesh,*c0);chf_it++)
            {
                hf1=mesh->s_opposite_halfface(quote(chf_it));
                c1=(template_c*)(hf1->cell);
                if(c1==NULL)
                {continue;} 
                rbt.key=c1->id;
                rbt1=(RB_int*)(tree->find(tree,&rbt));
                if(rbt1!=NULL)
                {
                    continue;
                } 
                else
                {
                    for(int i=0;i<rows-1;i++)
                    {
                        for(int j=0;j<cols;j++)
                        {
                            M[i][j]=((template_v*)(c1->vertices[i]))->point[j];
                        }
                    }
                    Tensor*t1=Anti_tensor_mpf_from_point(tas,M,rows,cols);
          //          printf("cid:%d\n",c1->id);
          //          tensor_mpf_print_self(t1);
                    if(!libcell_is_same_dir(tas,t,t1))
                    {
          //              printf("shi\n");
                        node=node_overlying(node,c1);
                        node2=node_overlying(node2,c1);
                      //  c1->prop=(void*)temp_int;
                    }
                    else
                    {
                     //   c1->prop=(void*)temp_int1;
                    }
                    tas->T_free(tas,t1);
                    rbt.key=c1->id;
                    tree->insert(tree,&rbt);

                }
            }

        } 
        free_node(temp_node);
    }
    RB_Tree_free(tree); 
    //reset_c_prop(mesh);
    Node* node1=node2;
    while(node2!=NULL)
    {
    //    printf("delete id:%d\n",((template_c*)(node2->value))->id);
        mesh->delete_cell(mesh,*((template_c*)(node2->value)),true);
        node2=(Node*)(node2->Next); 
    }
    //mesh->printself(mesh);
    free_node(node1);
    node1=NULL;
    mesh->external_cell_init_(mesh);
    template_v**temp_v=(template_v**)malloc(sizeof(template_v*)*(rows-1));
    temp_v[rows-2]=v;
    
    for(Node* hfit=(mesh->external_cell.halffaces);hfit!=NULL;hfit=(Node*)(hfit->Next))
    {
       // printf("external size:%d\n",((template_hf*)(hfit->value))->vertices_size);
        for(int i=0;i<((template_hf*)(hfit->value))->vertices_size;i++)
        {

            temp_v[i]=((template_hf*)(hfit->value))->vertices[i];
        }
        //解决halfface退化为点的情况
        if(rows==3)
        {
            hf1=mesh->s_opposite_halfface((template_hf*)(hfit->value));
                template_v* temp_vv[2];

            if(hf1->cell->vertices[0]==temp_v[0])
            {
                temp_vv[0]=temp_v[1];temp_vv[1]=temp_v[0];
                c0=(mesh->create_cellv(mesh,temp_vv,rows-1));
                
            } 
            else
            {
                c0=mesh->create_cellv(mesh,temp_v,rows-1);
            }
        }
        else
        {

            c0=mesh->create_cellv(mesh,temp_v,rows-1);
        }
        if(c0==NULL)
        {
            for(int i=0;i<rows;i++)
            {
                free(M[i]);
            }
            free(M);
   // free(temp_int);
   // free(temp_int1);
            free(temp_v);
            return false;
        }
    }
   // mesh->printself(mesh);
    for(int i=0;i<rows;i++)
    {
        free(M[i]);
    }
    free(M);
   // free(temp_int);
   // free(temp_int1);
    free(temp_v);
    return true;
}

//mesh是只包含顶点集合的mesh
//t是凸包用到的反对称张量作为方向
bool mesh_createconvex(Tensors_Algebra_System*tas,Tensor* t,template_m*m)
{
    if(m->num_v(m)<=0)
    {
        return false;
    }
    m->simplex=1;
#ifdef MANIFOLD_REQUIRE
    m->dimension=t->order(t)-1;
#endif
    auto iter=m->vertices.begin();
    m->simplex=1;m->dimension=t->order(t)-1;
    int rows=m->dimension+2,cols=iter->second->point_size;
    double **VV=(double**)malloc(sizeof(double*)*(rows));
    for(int i=0;i<rows;i++)
    {
        VV[i]=(double*)malloc(sizeof(double)*(cols));
    }
    template_v** temp_v=(template_v**)malloc(sizeof(template_v*)*rows); 
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<cols;j++)
        {
             VV[i][j]=iter->second->point[j];
        }
        temp_v[i]=iter->second;
        iter++;
    }
    int *index=(int*)malloc(sizeof(int)*(rows));
    for(int i=0;i<rows;i++)
    {
        index[i]=i;
    }
    Tensor* t1=Anti_tensor_mpf_from_point(tas,VV,rows,cols);
    //tensor_mpf_print_self(t1);
    if(!libcell_is_same_dir(tas,t,t1))
    {
        index[0]=1;
        index[1]=0;
    }
    tas->T_free(tas,t1);
   
    template_v **temp_v1=(template_v**)malloc(sizeof(template_v*)*cols);
    template_v*v0=NULL;
    template_c*c0=NULL;
    for(int i=0;i<rows;i++)
    {
        int k=0;
        for(int j=0;j<rows;j++)
        {
            if(j!=i)
            {
                temp_v1[k]=temp_v[index[j]];
                k++;
            }

        }
        if((rows-1-i)%2!=0)
        {
            v0=temp_v1[0];
            temp_v1[0]=temp_v1[1];
            temp_v1[1]=v0;
        }
        c0=m->create_cellv(m,temp_v1,cols);
        if(c0==NULL)
        {
            free(index);
            for(int i=0;i<rows;i++)
            {
                free(VV[i]);
            }
            free(VV);
            free(temp_v);free(temp_v1);
            printf("mesh createconvex\n");
            return false;
        }
    }
    free(index);
    for(int i=0;i<rows;i++)
    {
        free(VV[i]);
    }
    free(VV);
    //printf("beigin\n");
    for(;iter!=m->vertices.end();iter++)
    {
        if(!increasing_convex_hull(tas,t,m,iter->second))
        {
            free(temp_v);free(temp_v1);
            return false;
        } 
    }
    free(temp_v);free(temp_v1);
    return true;
}
//给一个张量的代数系统，和反对称张量做方向，计算点集的凸包
bool  from_v_createconvex(Tensors_Algebra_System* tas,Tensor* t,template_m* mesh,double** VV,int rows,int cols)
{
    mesh->simplex=1;
    template_c*c0=NULL;
    int dim=t->order(t);
#ifdef MANIFOLD_REQUIRE
    mesh->dimension=dim-1;
#endif
    //***初始化一个单形
    for(int i=0;i<dim+1;i++)
    {
        mesh->create_vertexv(mesh,VV[i],cols);
    }
    int *index=(int*)malloc(sizeof(int)*(dim+1));
    for(int i=0;i<dim+1;i++)
    {
        index[i]=i;
    } 
    Tensor* t1=Anti_tensor_mpf_from_point(tas,VV,dim+1,cols);
    //tensor_mpf_print_self(t1);
    if(!libcell_is_same_dir(tas,t,t1))
    {
        index[0]=1;
        index[1]=0;
    }
    tas->T_free(tas,t1);
    //*****
    template_v **temp_v1=(template_v**)malloc(sizeof(template_v*)*(dim+1));
    template_v*v0=NULL;
    for(int i=0;i<dim+1;i++)
    {
        int k=0;
        for(int j=0;j<dim+1;j++)
        {
            if(j!=i)
            {
                temp_v1[k]=mesh->vertices.find(index[j])->second;
                k++;
            }
        }
        if((dim-i)%2!=0)
        {
            v0=temp_v1[0];
            temp_v1[0]=temp_v1[1];
            temp_v1[1]=v0;
        }
        c0=mesh->create_cellv(mesh,temp_v1,dim);
        if(c0==NULL)
        {
            return false;
        }
    }
    free(index);
    free(temp_v1);
    mesh->printself(mesh);
    for(int i=dim+1;i<rows;i++)
    {
        template_v* v=mesh->create_vertexv(mesh,VV[i],cols);
        if(!increasing_convex_hull(tas,t,mesh,v))
        {
            return false;
        }
    }
    return true;

}
//求cell的体积
__mpf_struct* compute_simplex_cell_volume(Tensors_Algebra_System*tas,template_m *own,const template_c*c)
{
    if(own->simplex!=1)
    {
        printf("is not simplex\n");
        return NULL;
    }
    int rows=c->vertices_size;

    int cols=((template_v*)(c->vertices[0]))->point_size; 
    double**M=(double**)malloc(sizeof(double*)*rows);
    for(int i=0;i<rows;i++)
    {
       M[i]=(double*)malloc(sizeof(double)*cols);        
    }
    int i=0;
    for(auto cv_iter=own->cv_begin(own,*c);cv_iter!=own->cv_end(own,*c);cv_iter++)
    {
        for(int j=0;j<(*cv_iter).point_size;j++) 
        {
            M[i][j]=(*cv_iter).point[j];
        }
        i++; 
    }
    __mpf_struct* re=area_simplex(tas,M,rows,cols);
    for(i=0;i<rows;i++)
    {
        free(M[i]);
    }
    free(M);
    return re;
}

//借助凸包算法的剖分算法，注意和凸包算法的不一样
//t是这些点所在子空间对应的反对称张量

bool convex_subdivision(Tensors_Algebra_System*tas,Tensor*t,Mesh* mesh,double **VV,int rows,int cols)
{ 
    Tensor* t1=NULL,*t2=NULL;
    int k=0;
    int size=0;
    if(t->order(t)==cols)
    {
        size=cols+1; 
    }
    else
    {
        size=cols;
    }
    double *temp_v=(double*)malloc(sizeof(double)*size);
    memset(temp_v,0,sizeof(double)*size);
    //保证tas->as->elements->size和size一致
    for(unsigned int i=0;i<tas->as->elements->size;i++)
    {
        t1=tas->T_create();
        k=i;
        t1->insert(tas->as,t1,&k,1,tas->copy_from_double(1));
        t2=Tensor_Wedge_(tas,t,t1);
        if(t2->value->size!=0)
        {
            temp_v[i]=10000;
            tas->T_free(tas,t1);
            t1=NULL;
            break;
        }

        tas->T_free(tas,t2);
        tas->T_free(tas,t1);
    }
    //printf("t2\n");
    //tensor_mpf_print_self(t2);
    //Mesh mesh;
    //Mesh_init(&mesh);
    template_v*v0=mesh->create_vertexv(mesh,temp_v,size);
    temp_v[size-1]=0;
    //memset(temp_v,0,sizeof(double)*size);
    for(int i=0;i<rows;i++)
    {
        if(size!=cols)
        {
            for(int j=0;j<cols;j++)
            {
                temp_v[j]=VV[i][j];
            }
            mesh->create_vertexv(mesh,temp_v,size);
        }
        else
        {
            mesh->create_vertexv(mesh,VV[i],cols);
        }
    } 
    bool re=mesh_createconvex(tas,t2,mesh);
    tas->T_free(tas,t2);
    free(temp_v);
    if(re){mesh->delete_vertex(mesh,*v0,true);}
    return re;
}
//cols背景空间
//给出mesh空间的反对称张量t
bool delauny_subdivision(Tensors_Algebra_System* tas,Tensor*t,Mesh*mesh,double**VV,int rows,int cols)
{

    Tensor* t1=NULL,*t2=NULL;
    int k=0;
    int size=0;
    if(t->order(t)==cols)
    {
        size=cols+1; 
    }
    else
    {
        size=cols;
    }
    double *temp_v=(double*)malloc(sizeof(double)*size);
    memset(temp_v,0,sizeof(double)*size);
    //保证tas->as->elements->size和size一致
    //这里应该是寻找法向量t1
    //然后将点延向量t1提升
    for(unsigned int i=0;i<tas->as->elements->size;i++)
    {
        t1=tas->T_create();
        k=i;
        t1->insert(tas->as,t1,&k,1,tas->copy_from_double(1));
        t2=Tensor_Wedge_(tas,t,t1);
        if(t2->value->size!=0)
        {
            temp_v[i]=10000;
            tas->T_free(tas,t1);
            t1=NULL;
            break;
        }

        tas->T_free(tas,t2);
        tas->T_free(tas,t1);
    }
    double temp_d=0;
    template_v* v0=mesh->create_vertexv(mesh,temp_v,size);
    //memset(temp_v,0,sizeof(double)*size);
    temp_v[k]=0;
    //temp_v[size-1]=0;
    for(int i=0;i<rows;i++)
    {
        if(size!=cols)
        {
            temp_d=0;
            for(int j=0;j<cols;j++)
            {
                temp_d+=VV[i][j]*VV[i][j];
                temp_v[j]=VV[i][j];
            }
            temp_v[k]=temp_d/2.0;
            mesh->create_vertexv(mesh,temp_v,size);
        }
        else
        {
            temp_d=0;
            for(int j=0;j<cols;j++)
            {
                temp_v[j]=VV[i][j];
                temp_d+=VV[i][j]*VV[i][j];
            }
            temp_v[k]+=temp_d/2.0;
            mesh->create_vertexv(mesh,temp_v,cols);
        }
    }
    //tensor_mpf_print_self(t2);
    //printf("k:%d\n",k);
    bool re= mesh_createconvex(tas,t2,mesh);
    if(!re)
    {
    }
    else
    {

        mesh->delete_vertex(mesh,*v0,true);
    }
    tas->T_free(tas,t2);

    free(temp_v);
     
    for(auto iter=mesh->vertices.begin();iter!=mesh->vertices.end();iter++)
    {
        iter->second->point[k]=0;
    } 
    return re;
     
}
//t是mesh所在子空间对应的反对称张量
//借助了凸包算法
//求rows个点所围凸包的面积
//似乎这样的算法不稳定
__mpf_struct* compute_convex_area(Tensors_Algebra_System*tas,Tensor* t,double **VV,int rows,int cols)
{
    Mesh mesh; 
    Mesh_init(&mesh);
    __mpf_struct* re=(__mpf_struct*)malloc(sizeof(__mpf_struct)),*re1=NULL;
    mpf_inits(re,NULL);
    mpf_set_ui(re,0);
    //mesh.simplex=1;
    //还有一种方案是
    /*
    convex_subdivision(tas,t,&mesh,VV,rows,cols);
    */

    if(!convex_subdivision(tas,t,&mesh,VV,rows,cols))
    {
        return re;
    }
    for(auto c_it=mesh.cells.begin();c_it!=mesh.cells.end();c_it++)
    {
        re1=compute_simplex_cell_volume(tas,&mesh,c_it->second);
        //printf("libo\n");
        mpf_add(re,re,re1);
        mpf_clear(re1);
        free(re1);
        //printf("end\n");
    }
    Mesh_free(&mesh);
    return re;
}
//输入同上
__mpf_struct* compute_delauny_convex_area(Tensors_Algebra_System*tas,Tensor* t,double **VV,int rows,int cols)
{
    Mesh mesh; 
    Mesh_init(&mesh);
    __mpf_struct* re=(__mpf_struct*)malloc(sizeof(__mpf_struct)),*re1=NULL;
    mpf_inits(re,NULL);
    mpf_set_ui(re,0);
    //mesh.simplex=1;
    //还有一种方案是
    /*
    convex_subdivision(tas,t,&mesh,VV,rows,cols);
    */
    if(!delauny_subdivision(tas,t,&mesh,VV,rows,cols))
    {
        return re;
    }
    for(auto c_it=mesh.cells.begin();c_it!=mesh.cells.end();c_it++)
    {
        re1=compute_simplex_cell_volume(tas,&mesh,c_it->second);
        mpf_add(re,re,re1);
        mpf_clear(re1);
        free(re1);
    }
    Mesh_free(&mesh);
    return re; 
}
/*double area_simplex_double(double**M,int rows,int cols)
{
    //rows respect point size
    if(rows>(cols+1))
    {
        printf("points are too much\r\n");
        return 0;
    }
    Eigen::MatrixXd A(rows-1,cols);
    for(int i=0;i<(rows-1);i++)
    {
        for(int j=0;j<cols;j++)
        {
            A.coeffRef(i,j)=M[i][j]-M[rows-1][j];

        }
    }
    double test=(A*A.transpose()).determinant();
    if(test<0)
    {
        test=-test;
       // printf("xiaoyu 0  :%lf\n",test);
    }
    double re=sqrt(test)/factorial(rows-1);

    if(isnan(re)==1)
    {
        printf("nan:%lf\n",re);
    }
    return re;
}
*/




#endif
