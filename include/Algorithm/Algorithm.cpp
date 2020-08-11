#include<Algorithm/Algorithm.h>
//
double* temp_normal_vector(Eigen::MatrixXd& A)
{
    int dim=A.rows();
    double * re=(double*)malloc(sizeof(double)*(dim+1));
    Eigen::MatrixXd temp;
    temp.resize(dim,dim);
    int l=0;
    for(int i=0;i<dim+1;i++)
    {
        for(int k=0;k<dim;k++)
        {
            l=0;
        
            for(int j=0;j<dim+1;j++)
            {
                if(j!=i)
                {
                    temp.coeffRef(k,l)=A.coeff(k,j);
                    l++;
                }
            
            }
        }
    
        re[i]=temp.determinant();
        if((dim-i)%2==1)
        {
            re[i]=-re[i];
        
        } 
    }
    return re;
}

//#include<Cell_Traits.h>
double  factorial(int n)
{
	int result=1;
    for(int i=2;i<=n;i++)
    {
        result=result*i;
    }
    return (double)result;
}
//rows代表点数，cols代表背景空间，比如3,7表示7维背景空间下三角形(单形))
double area_simplex(double**M,int rows,int cols)
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
//rows代表点数，cols代表背景空间
double ori_area_simplex(double ** M,int rows,int cols)
{
    if((rows-1)!=cols)
    {
        printf("rows is not corresponse to cols\r\n");
        return 0;
    }
    Eigen::MatrixXd A(rows,rows);
    for(int i=0;i<rows;i++)
    {
        for(int j=0;j<rows;j++)
        {
            if(j!=0)
            {
                A.coeffRef(i,j)=M[i][j-1];
            }
            else
            {
                A.coeffRef(i,j)=1;
            }
        }
    }
    return A.determinant()/factorial(rows-1);

}
Eigen::MatrixXd orthogonal_complementation(Eigen::MatrixXd A)
{
    //Eigen::EigenSolver<Eigen::MatrixXd> egs(A);
    //Eigen::LLT<Eigen::MatrixXd> llt;
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A);
    Eigen::MatrixXd k=lu_decomp.kernel();
   //lu_decomp.image(A) 
    //std::cout<<A.lu().matrixLU();

    return k;
}

void set_cell_normal(template_m*mesh)
{
    Eigen::MatrixXd A;
    Eigen::VectorXd temp_A;
    int n=(mesh->vertices.begin())->second->point_size;
    A.resize(n,n);
    A.setZero(); 
    for(auto cit=mesh->cells.begin();cit!=mesh->cells.end();cit++)
    {
        int length=0;    
        if(mesh->simplex==1)
        {
            length=cit->second->vertices_size;
        }
        else
        {
            length=n;
        
        }
        for(int i=1;i<length;i++)
        {
            template_v* temp_v=(template_v*)cit->second->vertices[0],*temp_v1=(template_v*)cit->second->vertices[i];
            for(int j=0;j<n;j++)
            {
                A.coeffRef(i-1,j)=temp_v1->point[j]-temp_v->point[j];   
            }
                        
        } 
        double *nor=(double*)malloc(sizeof(double)*n),sum=0;
        //Eigen::MatrixXd k=A.adjoint();
        //Eigen::MatrixXd temp_A(n-1,n-1);
        for(int i=0;i<n;i++)
        {
            temp_A=A.col(i);
            A.col(i)=A.col(0);
            nor[i]=-A.block(0,1,n-1,n-1).determinant();
            A.col(i)=temp_A;
            sum+=nor[i]*A.coeff(0,i);
            
           // printf("%lf  ",nor[i]);
        }
       // printf("sum:%lf\n",sum);
       // printf("\n");
        nor[0]=-nor[0];
        cit->second->prop=(void*)nor;
    }

}
#ifdef MANIFOLD_REQUIRE
void increasing_convex_hull(template_m* mesh,template_v* v)
{
 //   printf("begin convex hull\n");
    if(mesh->simplex!=1)
    {
        return;
    }
	int rows=v->point_size+1,cols=v->point_size;
	double **M=(double**)malloc(sizeof(double*)*rows);
	for(int i=0;i<rows;i++)
	{
	    M[i]=(double*)malloc(sizeof(double)*cols);	
	}
	for(int j=0;j<cols;j++)
	{
	    M[rows-1][j]=v->point[j];
	
	}
    Node* node=0,*node2=0;
    
//find first
    int *temp_int=(int*)malloc(sizeof(int)),*temp_int1=(int*)malloc(sizeof(int));
    *temp_int=1;*temp_int1=0;
    for(auto cit=mesh->cells.begin();cit!=mesh->cells.end();cit++)
    {
        for(int i=0;i<rows-1;i++)
	    {
	        for(int j=0;j<cols;j++)
	        {	
                M[i][j]=((template_v*)(cit->second->vertices[i]))->point[j];		
	        }
	    }
       	if(ori_area_simplex(M,rows,cols)<0)
        {
            node=node_overlying(node,(void*)cit->second);
            cit->second->prop=(void*)temp_int;
            node2=node_overlying(node2,(void*)cit->second);

            break;
        }
        else
        {
            cit->second->prop=(void*)temp_int1;
        } 
    }
    template_c* c0=NULL,*c1=NULL;
    template_hf*hf1=NULL;
    //广度优先
    while(node!=NULL)
    {
        Node * temp_node=node_copy(node);
        free_node(node);
        node=NULL;
        for(Node node_it=*temp_node;*node_it!=NULL;node_it++)
        {
            c0=(template_c*)(node_it.value);
            for(auto chf_it=mesh->chf_begin(mesh,*c0);chf_it!=mesh->chf_end(mesh,*c0);chf_it++)
            {
                hf1=mesh->s_opposite_halfface(quote(chf_it));
                c1=(template_c*)(hf1->cell);
		        if(c1==NULL)
		        {continue;}
                
                if(c1->prop!=NULL)
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
       	            if(ori_area_simplex(M,rows,cols)<0)
                    {
     //                       printf("push cell id:%d\n",(*vc_it).id); 
                        node=node_overlying(node,c1);
                        node2=node_overlying(node2,c1); 
                        c1->prop=(void*)temp_int; 
                    }
                    else
                    {
                        c1->prop=(void*)temp_int1; 
                    }

                }
            }
        }
        free_node(temp_node);
    }
    Node *node1=node2;
    while(node2!=NULL)
    {
  //      printf("delete cell id:%d\n",((template_c*)(node2->value))->id);
//        int myid=((template_c*)(node2->value))->id;
	    mesh->delete_cell(mesh,*((template_c*)(node2->value)),true);
	
/*        if(myid==5)
        {
            template_v* myv[2];
            auto it=mesh->vertices.find(1);
            myv[0]=it->second;
            auto it1=mesh->vertices.find(3);
            myv[1]=it1->second;
            auto it2=mesh->vertices.find(4);
            template_f* myf=mesh->get_facev(mesh,myv,2);
            if(mesh->face_is_boundary(mesh,myf))
            {
                if(myf->prop!=0)
                {
                    printf("prop ");
                }
                printf("f:id %dshi1\n",myf->id);
            }
            myv[1]=it2->second;
            myf=mesh->get_facev(mesh,myv,2);
            if(mesh->face_is_boundary(mesh,myf))
            {
                if(myf->prop!=0)
                {
                    printf("prop ");
                }

                printf("f:id %dshi2\n",myf->id);
            }
            myv[0]=it1->second;
            myf=mesh->get_facev(mesh,myv,2);

            if(mesh->face_is_boundary(mesh,myf))
            {
                if(myf->prop!=0)
                {
                    printf("prop ");
                }

                printf("f:id %dshi3\n",myf->id);
            }
        }*/
        
	    node2=(Node*)(node2->Next);
//delete_cell
    }
    free_node(node1);
    node1=0;
    reset_c_prop(mesh);
    /*for(auto cit=mesh->cells.begin();cit!=mesh->cells.end();cit++)
    {
        if((cit->second->prop)==(void*)temp_int)
        {
        
            mesh->delete_cell(mesh,*(cit->second),true); 
        }
        else
        {
        } 
        (cit->second->prop)=NULL;
    
    }*/
    mesh->external_cell_init_(mesh);
    template_v* temp_v[cols];
    temp_v[cols-1]=v;
    for(Node* hfit=(mesh->external_cell.halffaces);hfit!=NULL;hfit=(Node*)(hfit->Next))
    {   if(((template_hf*)(hfit->value))->cell!=NULL)
        {
            printf("here\n");
        }
        for(int i=0;i<((template_hf*)(hfit->value))->vertices_size;i++)
        {
            temp_v[i]=(template_v*)((template_hf*)(hfit->value))->vertices[i];
        }
        mesh->create_cellv(mesh,temp_v,cols); 
    }
    for(int i=0;i<rows;i++)
    {
        free(M[i]);
    }
    free(M);
    free(temp_int);
    free(temp_int1);
   // printf("begin convex hull\n");
}
//z只包含顶点集的mesh
void mesh_createconvex(template_m*m)
{
    if(m->num_v(m)<=0)
    {
        return;
    }
    auto iter=m->vertices.begin();
    m->simplex=1;m->dimension=iter->second->point_size-1;
    int rows=m->dimension+2,cols=m->dimension+1;
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
    //int i=0;
    int *index=(int*)malloc(sizeof(int)*(rows));
    for(int i=0;i<rows;i++)
    {
        index[i]=i;
    }
    if(ori_area_simplex(VV,rows,cols)<0)
    {
        index[0]=1;
        index[1]=0;   
    }
    template_v **temp_v1=(template_v**)malloc(sizeof(template_v*)*cols);
    template_v*v0=NULL;
    for(int i=0;i<rows;i++)
    {
        int k=0;
//      template_v* temp_v[cols];
        for(int j=0;j<rows;j++)
        {
            if(j!=i)
            {
                temp_v1[k]=temp_v[index[j]];
                k++;
            }

        }
        if((cols-i)%2!=0)
        {
            v0=temp_v1[0];
            temp_v1[0]=temp_v1[1];
            temp_v1[1]=v0;
        }
        m->create_cellv(m,temp_v1,cols);
    }
    //free(temp_v);
    free(index);
    for(int i=0;i<rows;i++)
    {
        free(VV[i]);
    }
    free(VV);
    for(;iter!=m->vertices.end();iter++)
    {
        increasing_convex_hull(m,iter->second); 
    }
    free(temp_v);free(temp_v1);
}

void  from_v_createconvex(template_m* mesh,double** VV,int rows,int cols)
{
    mesh->simplex=1;
#ifdef MANIFOLD_REQUIRE
    mesh->dimension=cols-1;
#endif
    for(int i=0;i<cols+1;i++)
    {
        mesh->create_vertexv(mesh,VV[i],cols);

    }
    int *index=(int*)malloc(sizeof(int)*(cols+1));
    for(int i=0;i<cols+1;i++)
    {
        index[i]=i;
    }
    if(ori_area_simplex(VV,cols+1,cols)<0)
    {
        index[0]=1;
        index[1]=0;   
    }
  template_v **temp_v=(template_v**)malloc(sizeof(template_v*)*cols);

    for(int i=0;i<cols+1;i++)
    {
        int k=0;
//        template_v* temp_v[cols];
        for(int j=0;j<cols+1;j++)
        {
            if(j!=i)
            {
                temp_v[k]=mesh->vertices.find(index[j])->second;
                k++;
            }
        }
        if((cols-i)%2!=0)
        {
            template_v* temp_v1=temp_v[0];
            temp_v[0]=temp_v[1];
            temp_v[1]=temp_v1;
        }
        mesh->create_cellv(mesh,temp_v,cols);
    }
    free(temp_v);
    free(index);
	//mesh->printself(mesh);
    for(int i=cols+1;i<rows;i++)
    {
        template_v* v=mesh->create_vertexv(mesh,VV[i],cols);
        increasing_convex_hull(mesh,v);
	    //mesh->printself(mesh);
    }
}
//
void from_v_createdelauny_simplex(template_m*mesh,double**VV,int rows,int cols)
{
#ifdef MANIFOLD_REQUIRE
    mesh->dimension=cols;
#endif
    mesh->simplex=1;
    double **vv=(double**)malloc(sizeof(double*)*(rows+1));
    for(int i=0;i<=rows;i++)
    {
        vv[i]=(double*)malloc(sizeof(double)*(cols+1));
    }
    memset(vv[0],0,sizeof(double)*cols);
    vv[0][cols]=1000.0;
    double sum=0;
    for(int i=0;i<rows;i++)
    {
        sum=0;
        for(int j=0;j<cols;j++)
        {
            vv[i+1][j]=VV[i][j];
            sum+=(VV[i][j]*VV[i][j])/2.0;
        }
        vv[i+1][cols]=sum;
    }
    from_v_createconvex(mesh,vv,rows+1,cols+1);
    auto v_it=mesh->vertices.begin();
 /*   for(auto it=mesh->vc_begin(mesh,*(v_it->second));it!=mesh->vc_end(mesh,*(v_it->second));it++)
    {
        mesh->delete_cell(mesh,*(it),true); 
    }*/
    mesh->delete_vertex(mesh,*(v_it->second),true);
   
    for(int i=0;i<=rows;i++)
    {
        free(vv[i]);
    }
    free(vv);
}
double compute_simplex_cell_volume(template_m *own,const template_c*c)
{
    if(own->simplex!=1)
    {
        printf("is not simplex\n");
        return -1;
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
    double re=area_simplex(M,rows,cols);
    for(i=0;i<rows;i++)
    {
        free(M[i]);
    }
    free(M);
    return re;
}
//求rows个点所围凸包的面积，由于没有反对称张量，所以要限制背景空间
double compute_convex_area(double **v,int rows,int cols,int dim)
{
    if((rows-1)<dim)
    {
        printf("rows is not correct \n");
        return -1; 
    }
//因为没有反对称张量库，所以只能限制背景空间
    if(dim!=(cols-1))
    {
        printf("cols is not  correspond to dim\n");
        return -1;
    }
    double **vv=(double**)malloc(sizeof(double*)*(rows+1));
    for(int i=0;i<(rows+1);i++)
    {
        vv[i]=(double*)malloc(sizeof(double)*(cols));
    }
    memset(vv[0],0,sizeof(double)*(cols-1));   
    vv[0][cols-1]=10000.0;
    for(int i=1;i<(rows+1);i++)
    {
        for(int j=0;j<cols;j++)
        {
            vv[i][j]=v[i-1][j];
        }
    }
    Mesh mesh; 

    A_Mesh_init(&mesh);
    //mesh.simplex=1;
    from_v_createconvex(&mesh,vv,rows+1,cols);
    auto iter=mesh.vertices.begin();
    mesh.delete_vertex(&mesh,*(iter->second),true);
   // mesh.printself(&mesh);
 //   mesh.printself(&mesh);
    double re=0;
    for(auto c_it=mesh.cells.begin();c_it!=mesh.cells.end();c_it++)
    {
        re+=mesh.cell_volume(&mesh,c_it->second);
    }

    Mesh_free(&mesh);
    for(int i=0;i<(rows+1);i++)
    {
        free(vv[i]);
    }
    free(vv);
    return re;
}
bool libcell_face_is_convex(template_m*m,template_f*f)
{
    if(m->face_is_boundary(m,f))
    {
        return true;
    }
    template_c*c0=(template_c*)(f->halffaces[0]->cell),*c1=(template_c*)(f->halffaces[1]->cell);
    template_v*v0=NULL;
    int flag=1;
    for(auto cv_it=m->cv_begin(m,*c0);cv_it!=m->cv_end(m,*c0);cv_it++)
    {
        flag=1;
        for(auto cv_it1=m->cv_begin(m,*c1);cv_it1!=m->cv_end(m,*c1);cv_it1++)
        {
            if(quote(cv_it)==quote(cv_it1))
            {
                flag=0;
                break;
            }
        
        }
        if(flag==1)
        {
            v0=(template_v*)quote(cv_it);
            break;
        }
    
    }
    int rows=c0->vertices_size+1,cols=c0->vertices_size;
    double **VV=(double**)malloc(sizeof(double*)*rows);
    for(int i=0;i<rows;i++)
    {
        VV[i]=(double*)malloc(sizeof(double)*cols); 
    }
    for(int j=0;j<cols;j++)
    {
    
        VV[rows-1][j]=v0->point[j];
    }
   // VV[rows-1][cols-1]=*((double*)(v0->prop));
    for(int i=0;i<rows-1;i++)
    {
        for(int j=0;j<cols;j++)
        {
            VV[i][j]=((template_v*)(c1->vertices[i]))->point[j];
        }
      //  VV[i][cols-1]=*((double*)(((template_v*)(c1->vertices[i]))->prop));
     }
    if(ori_area_simplex(VV,rows,cols)<0)
    {
    
        return false;
    }
    return true;

}

//概念不对，被弃用,在二维有用
bool libcell_is_flip(template_m*mesh,template_f*f)
{
    int rows=f->vertices_size;
    template_v* v=((template_v*)(f->vertices[0]));
    //因为没有反对称张量所以限制背景空间
    if(rows!=v->point_size)
    {
        printf("libcell_is_filp:background dim is not right\n");
        return false;
    }
    template_c* c0=0,*c1=0;template_hf*hf0=0,*hf1=0;
    
    hf0=f->halffaces[0];
    
    hf1=f->halffaces[1];
    
    c0=(template_c*)(hf0->cell);c1=(template_c*)(hf1->cell);
    if(c0==0||c1==0)
    {
        return false;
    }
    template_v* v1=0;
    int flag=1;
    for(auto iter=mesh->cv_begin(mesh,*c1);iter!=mesh->cv_end(mesh,*c1);iter++)
    {
        flag=1;
        for(auto iter1=mesh->hfv_begin(mesh,*hf1);iter1!=mesh->hfv_end(mesh,*hf1);iter1++)
        {
            if((*iter).id==(*iter1).id)
            {   
                flag=0;
                break;
            }
        
        }
        if(flag==1)
        {
            v1=quote(iter);
            break;
        }
    }
    if(v1==0)
    {
        printf("here\n");
        return false;
    }
    double **M=(double**)malloc(sizeof(double*)*(rows+1));
    for(int i=0;i<=rows;i++)
    {
        M[i]=(double*)malloc(sizeof(double)*rows);
    }
    for(int j=0;j<rows;j++)
    {
        M[rows][j]=v1->point[j];
    }
    
    for(auto iter=mesh->chf_begin(mesh,*c0);iter!=mesh->chf_end(mesh,*c0);iter++)
    {
        for(int i=0;i<(*iter).vertices_size;i++)
        {
            v=(template_v*)((*iter).vertices[i]);
            for(int j=0;j<(v->point_size);j++)
            {
                M[i][j]=v->point[j];
            }
        }
        double temp_d=ori_area_simplex(M,rows+1,rows);
        if((*iter).id==hf0->id)
        {
            if(temp_d>0)
            {
                return false;
            }
        }
        else
        {
            if(temp_d<0)
            {
                return true;
            }
        }
    }
    return true;    

}
//暂时对n维背景空间的n维流形判断
void A_Mesh_init(struct Mesh*own)
{
    Mesh_init(own);
    own->cell_volume=compute_simplex_cell_volume;
}

#endif
