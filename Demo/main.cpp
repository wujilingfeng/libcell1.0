#include<stdio.h>
#include<Mesh/Mesh_Frame.h>
#include<Mesh/lib_cell_Iterator.h>
#include<Mesh_IO/Mesh_IO.h>
//#include<Algorithm/Algorithm1.h>
#include<Algorithm/Algorithm2.h>
#include<tool/libcell_tools_view.h>
//#include<Algorithm/Harmonic.h>
#include <Algorithm/Conformal.h>
void test_delauny() 
{
    srand((unsigned)time(0));
    Mesh mesh;
    Mesh_init(&mesh);
    double **v=(double**)malloc(sizeof(double*)*2000);
    for(int i=0;i<2000;i++)
    {
        v[i]=(double*)malloc(sizeof(double)*4);
    }
 /*   
    for(int i=0;i<700;i++)
    {
        double r=1,delta=(rand()%2000)/1000.0-1,theta=(rand()%1000)/1000.0;
        //theta=0.5;
        v[i][0]=r*sin(theta*M_PI)*cos(delta*M_PI);
        v[i][1]=r*sin(theta*M_PI)*sin(delta*M_PI);
        v[i][2]=r*cos(theta*M_PI); 
    }
*/
  /*  for(int i=0;i<2000;i++)
    {
        double r=(rand()%1000000)/1000000.0,delta=(rand()%2000000)/1000000.0-1,theta=(rand()%1000000)/1000000.0;
        theta=0.5;
        v[i][0]=r*sin(theta*M_PI)*cos(delta*M_PI);
        v[i][1]=r*sin(theta*M_PI)*sin(delta*M_PI);
        v[i][2]=r*cos(theta*M_PI);
          
    }*/
    for(int i=0;i<2000;i++)
    {
        //theta=0.5;
        v[i][0]=(rand()%20000000)/10000000.0-1;
        v[i][1]=(rand()%20000000)/10000000.0-1;
        v[i][2]=(rand()%20000000)/10000000.0-1;
        v[i][3]=(rand()%20000000)/10000000.0-1;
        v[i][4]=(rand()%20000000)/10000000.0-1;
    }
    for(int i=0;i<2000;i++)
    {
        double sum=0;
        for(int j=0;j<5;j++)
        {
            sum+=v[i][j]*v[i][j];
        }
        sum=sqrt(sum);
        if(sum>=1)
        for(int j=0;j<5;j++)
        {

            v[i][j]/=sum;
        } 
    }
    Tensors_Algebra_System*tas=(Tensors_Algebra_System*)malloc(sizeof(Tensors_Algebra_System));
    Tensors_Algebra_System_mpf_init(tas,6);//一般tas的维数要大一位
    Tensor*t=tas->T_create();
    int ids[5]={0,1,2,3,4};
    t->insert(tas->as,t,ids,5,tas->copy_from_double(1));
    tensor_mpf_print_self(t);
    //convex_subdivision(tas,t,&mesh,v,1110,2);
    if(!delauny_subdivision(tas,t,&mesh,v,10,5))
    {
        printf("liboodfsdferro\n");
    }
    //from_v_createconvex(tas,t,&mesh,v,1200,3);

    //_WriteCell_(&mesh,"surface.cell");
    //mesh.printself(&mesh);
    //要记得删除非流形点
    Node *nmv=mesh.non_manifold_vertices(&mesh);
    for(auto nit=nmv;nit!=NULL;nit=(Node*)(nit->Next))
    {
        mesh.delete_vertex(&mesh,*((template_v*)(nit->value)),true);
    }
    free_node(nmv);
    _WriteCell_(&mesh,"delauny_subdivision5.cell");
    Tensors_Algebra_System_free(tas);
    Mesh_free(&mesh);
/*    from_v_createdelauny_simplex(&mesh,v,2000,3);
    for(auto it=mesh.vertices.begin();it!=mesh.vertices.end();it++)
    {
        it->second->point_size--; 
    }

    Node*re=mesh.isolate_vertices(&mesh);
    Node* n_it=re;
    while(re!=NULL)
    {
        mesh.delete_vertex(&mesh,*((template_v*)(re->value)),true);
        re=(Node*)(re->Next);
    }
    free_node(n_it);
    _WriteCell_(&mesh,"delauny_sphere2.cell");
    
    Mesh_free(&mesh);
*/
    //Mesh_viewer_world mw;
    //Mesh_viewer_world_init(&mw);
    /*char ch[]="faces";
    Node* n=mw.create_something(&mw,ch);
    Mesh_viewer_something*ms=(Mesh_viewer_something*)(n->value);
    Mesh_viewer_faces *mf=(Mesh_viewer_faces*)(ms->evolution);
    mf->color_rows=mesh.num_v(&mesh);
    mf->random_color(mf);
    //mf->normal_rows=mesh.num
    get_data_from_3dim_cell(&mesh,&(mf->Data),&(mf->Data_index),&(mf->Data_rows),&(mf->Data_index_rows));
    mf->normal_rows=mf->Data_rows;
*/
    /*char sp_name[]="faces";
    Node* n=mw.create_something(&mw,sp_name);
    auto ms=(Mesh_viewer_something*)(n->value);
    auto mf=(Mesh_viewer_faces*)(ms->evolution);
    mf->color_rows=mesh.num_c(&mesh);
    mf->random_color(mf);
//    get_data_from_2dim_cell(&mesh,&(mf->Data),&(mf->Data_index));
    mf->normal_rows=mesh.num_c(&mesh);
    mf->Data_rows=mesh.num_v(&mesh);
    mf->Data_index_rows=mesh.num_c(&mesh);
    free_node(n);*/
    //test_camera_and_intera(&mw);
    //Mesh_viewer_opengl_Interpreter moi;
    //Mesh_viewer_opengl_Interpreter_init(&moi);
    //moi.world=&mw;
    //moi.routine_show(&moi);
    
}
void test_area()
{
    Mesh mesh;
    Mesh_init(&mesh);
    double **v=(double**)malloc(sizeof(double*)*8);
    for(int i=0;i<8;i++)
    {
        v[i]=(double*)malloc(sizeof(double)*3);
    }
    v[0][0]=1;v[0][1]=-1;v[0][2]=-1;
    v[1][0]=1;v[1][1]=-1;v[1][2]=1;
    v[2][0]=1;v[2][1]=1;v[2][2]=1;
    v[3][0]=1;v[3][1]=1;v[3][2]=-1;
    v[4][0]=-1;v[4][1]=-1;v[4][2]=-1;
    v[5][0]=-1;v[5][1]=-1;v[5][2]=1;
    v[6][0]=-1;v[6][1]=1;v[6][2]=1;
    v[7][0]=-1;v[7][1]=1;v[7][2]=-1;

    Tensors_Algebra_System*tas=(Tensors_Algebra_System*)malloc(sizeof(Tensors_Algebra_System));
    Tensors_Algebra_System_mpf_init(tas,4);
    Tensor*t=tas->T_create();
    int ids[3]={0,1,2};
    t->insert(tas->as,t,ids,2,tas->copy_from_double(1));
    tensor_mpf_print_self(t);
    convex_subdivision(tas,t,&mesh,v,4,3);

    auto re=compute_convex_area(tas,t,v,4,3);
    gmp_printf("re:%.Ff\n",re);

    
}
void test_convex()
{
    Mesh mesh;
    Mesh_init(&mesh);
    _ReadCell_(&mesh,"dual_topo_mesh.cell");
    Node*n=NULL;
    for(auto cit=mesh.cells.begin();cit!=mesh.cells.end();cit++)
    {
        n=node_overlying(n,cit->second);
    }
    for(Node*nit=n;nit!=NULL;nit=(Node*)(nit->Next))
    {
     //   printf("once\n");
        mesh.delete_cell(&mesh,*((template_c*)(nit->value)),true);
    }
    Tensors_Algebra_System*tas=(Tensors_Algebra_System*)malloc(sizeof(Tensors_Algebra_System));
    Tensors_Algebra_System_mpf_init(tas,4);
    Tensor*t=tas->T_create();
    int ids[4]={0,1,2,3};
    t->insert(tas->as,t,ids,4,tas->copy_from_double(1));
    mesh_createconvex(tas,t,&mesh);   
    printf("f:%d c:%d\n",mesh.num_f(&mesh),mesh.num_c(&mesh));
}
static bool zhuan_is_same_dir(Tensors_Algebra_System*tas,Tensor*t,template_v** v,int rows)
{
    int cols=v[0]->point_size;
    double **M=(double**)malloc(sizeof(double*)*rows);
    for(int i=0;i<rows;i++)
    {
        M[i]=(double*)malloc(sizeof(double)*cols);
        for(int j=0;j<cols;j++)
        {
            M[i][j]=v[i]->point[j]; 
        }
    } 
    printf("libo:%d %d \n",rows,cols); 

    Tensor* t1=Anti_tensor_mpf_from_point(tas,M,rows,cols);
  printf("liboend\n"); 
    for(int i=0;i<rows;i++)
    {
        free(M[i]);
        //M[i]=(double*)malloc(sizeof(double)*cols);
    }
    free(M);
    bool re=libcell_is_same_dir(tas,t,t1);
    tas->T_free(tas,t1);    
    return re;
}
bool adjust_3dim_cell_vertices_order(Mesh*mesh)
{
    if(mesh->simplex==1)
    {
        return true;
    }
    int dim=mesh->dimension;
    Tensors_Algebra_System*tas=(Tensors_Algebra_System*)malloc(sizeof(Tensors_Algebra_System));
    Tensors_Algebra_System_mpf_init(tas,dim);
    double ** M=(double**)malloc(sizeof(double*)*dim);
    for(int i=0;i<dim;i++)
    {
        M[i]=(double*)malloc(sizeof(double)*dim);
    }
    template_v**temp_v=(template_v**)malloc(sizeof(template_v*)*dim);
    template_v* v=NULL;
    bool flag=false;
    template_hf*hf=NULL;
    for(auto iter=mesh->halffaces.begin();iter!=mesh->halffaces.end();iter++)
    {
        hf=iter->second;
        for(int i=0;i<dim;i++)
        {
            for(int j=0;j<dim;j++)
            {
                M[i][j]=hf->vertices[i]->point[j];

                printf("%lf ",M[i][j]);
            }
            printf("\n");
        }
        Tensor*t=Anti_tensor_mpf_from_point(tas,M,dim,dim);
        //printf("here\n");
        //tensor_mpf_print_self(t);
       // break;
        printf("vsize:%d\n",iter->second->vertices_size);
        for(int i=dim;i<hf->vertices_size;i++)
        {
            printf("i:%d\n",i);
            temp_v[0]=hf->vertices[i-1];
            temp_v[1]=hf->vertices[0];

            temp_v[2]=hf->vertices[i];
            printf("hre\n");   
            flag=zhuan_is_same_dir(tas,t,temp_v,dim);
            

            printf("he\n");
            if(!flag)
            {
                printf("h\n");
                //break;
            }
            else
            {
                printf("hre\n");      
                for(int j=1;j<i;j++)
                {
                    temp_v[0]=hf->vertices[j-1];
                    temp_v[1]=hf->vertices[j];
                    temp_v[2]=hf->vertices[i];

                    printf("hre\n");
                    flag=zhuan_is_same_dir(tas,t,temp_v,dim);

                    printf("ib\n"); 
                    if(!flag)
                    {
                        v=hf->vertices[i];
                        for(int k=i;k>j;k--)
                        {
                            hf->vertices[k]=hf->vertices[k-1];
                        }
                        hf->vertices[j]=v;
                        break;
                    }
                }
            }
        }
        //printf("once\n");
        //tensor_mpf_print_self(t);
        //printf("onceend\n");
        tas->T_free(tas,t);
        
    }
    for(int i=0;i<dim;i++)
    {
        free(M[i]);
    }
    free(M); 
    return true;
}
void test_cell()
{
   Mesh mesh;
   Mesh_init(&mesh);
    _ReadCell_(&mesh,"out_mesh.cell");
    Tensors_Algebra_System* tas=(Tensors_Algebra_System*)malloc(sizeof(Tensors_Algebra_System));
    Tensors_Algebra_System_mpf_init(tas,3);
    
    adjust_3dim_cell_vertices_order(&mesh);
    for(auto iter=mesh.faces.begin();iter!=mesh.faces.end();iter++)
    {
       template_hf* hf1=iter->second->halffaces[0];
       template_hf* hf2=iter->second->halffaces[1];
       if(hf1->cell==NULL)
       {
           for(int i=0;i<hf2->vertices_size;i++)
           {
                int j=hf2->vertices_size-i-1;
                if(j>i)
                {
                    template_v* v=hf2->vertices[i];
                    hf2->vertices[i]=hf2->vertices[j];
                    hf2->vertices[j]=v;
                }

           } 
       }
       if(hf2->cell==NULL)
       {
            for(int i=0;i<hf1->vertices_size;i++)
            {
                int j=hf1->vertices_size-i-1;
                if(j>i)
                {
                    template_v* v=hf1->vertices[i];
                    hf1->vertices[i]=hf1->vertices[j];
                    hf1->vertices[j]=v;
                }

           } 
       }
        
    }
    _WriteCell_(&mesh,"out_mesh_adjust.cell");
    //Tensor* t=tas->T_create();
    //int id;
    //id=0;
    //t->insert(tas->as,t,&id,1,tas->copy_from_double(0));
    //id=1;
    //t->insert(tas->as,t,&id,1,tas->copy_from_double(0));
    //id=2;
    //t->insert(tas->as,t,&id,1,tas->copy_from_double(1));
    //Tensor*t1=Hodge_Anti_tensor_(tas,t);
    //Mesh mesh1;
    //Mesh_init(&mesh1);
    //from_v_createconvex(tas,t1,&mesh1,);
         
      
}
void test_harmonic()
{
    Mesh mesh;
    Mesh_init(&mesh);

}

void test_read_M()
{
    //Mesh mesh;
    //Mesh_init(&mesh);
    //Node* n_c=_ReadM_(&mesh,"bull1.m");

    //printf("%d %d \n",mesh.num_c(&mesh),mesh.num_v(&mesh));
   /* for(auto it=mesh.vertices.begin();it!=mesh.vertices.end();it++)
    {
        double *temp=(double*)(n_c->value);

        for(int i=0;i<3;i++)
        {
            
        }
        n_c=(Node*)(n_c->Prev);
    }*/ 
}
int main(int argc,char**argv)
{

    test_read_M();
 /*   mpf_set_default_prec(200);

    Mesh mesh,mesh1;
    Mesh_init(&mesh);Mesh_init(&mesh1);
  
    _ReadCell_(&mesh,"cube.cell");
    template_v v0=mesh.get_vertex(&mesh,0);
    template_f f0=mesh.get_face(&mesh,0);
    //mesh.delete_vertex(&mesh,v0,true);
    mesh.delete_face(&mesh,f0,true);
    double M[3][3],M1[3][2];
    M[0][0]=21.0;M[0][1]=-2.41;M[0][2]=-212.97;
    M[1][0]=1.001;M[1][1]=-0.51;M[1][2]=-2.4;
    M[2][0]=2.9;M[2][1]=92.25;M[2][2]=12.23;
    M1[0][0]=1;M1[0][1]=0;
    M1[1][0]=1;M1[1][1]=0;
    M1[2][0]=0;M1[2][1]=1;
    Tensors_Algebra_System*tas=(Tensors_Algebra_System*)malloc(sizeof(Tensors_Algebra_System));
    Tensors_Algebra_System_mpf_init(tas,3);
    Tensors_Algebra_System* tas1=(Tensors_Algebra_System*)malloc(sizeof(Tensors_Algebra_System));
    Tensors_Algebra_System_double_init(tas1,3);
    Tensor* t=tas1->T_create_matrix(tas1,(double*)M,3,3);
    //tas1->print_T(t);
    Tensor* t1=tas1->T_create_matrix(tas1,(double*)M1,3,2);
    tas1->print_T(t1);
    Tensor* t2=W_Tensor_Contraction(tas1,t,t1,1,0);
    tas1->print_T(t2);

    //tensor_double_print_self(t); 
    
    //test_area();
    //test_convex();
    test_delauny();
    */
    //test_cell();
    //printf("%f\n",area_simplex_double(M,3,3));
    //Tensor*t=Anti_tensor_mpf_from_v(tas,M,3,3);
    //tensor_mpf_print_self(t);
    //gmp_printf("%.Ff\n",(__mpf_struct*)tas->T_norm(tas,t));
     
   // mesh.printself(&mesh);
 //   printf("v:%d,f:%d ,c:%d\n",mesh1.num_v(&mesh1),mesh1.num_f(&mesh1),mesh1.num_c(&mesh1));

   /* for(auto it=mesh.cells.begin();it!=mesh.cells.end();it++)
    {
        //printf("d");
    
    }*/
   // _ReadCell_(&mesh,"hand.cell");
    printf("zuiend\n");
    return 0;
}
