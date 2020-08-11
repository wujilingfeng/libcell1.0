#ifndef LIB_CELL_HARMONIC_H
#define LIB_CELL_HARMONIC_H
#include<math.h>
#include<Mesh/_Func_.h>
#include<IterativeLinearSolvers>
#include<Sparse>
#include<Dense>
//#include<MatOp/SparseGenMatProd.h>
//#include<Tensors/Antisymmetric_Tensor.h>
#include<Tensors/Tensors_Operation.h>
#include<tools/double_tools.h>

#ifndef M_PI
#define M_PI 3.1415926
#endif
typedef struct Ham_v_prop
{
	int index;
	double ux,uy;
}Ham_v_prop;
typedef struct Ham_f_prop
{
	double weight;
}Ham_f_prop;
/*static void Ham_v_prop_init(template_v*v)
{
	Ham_v_prop* hvp=(Ham_v_prop*)malloc(sizeof(Ham_v_prop));
	hvp->index=0;
	v->user_prop=hvp;
}*/
/*static void  Ham_v_prop_free(template_v*v)
{
	Ham_v_prop* hvp=(Ham_v_prop*)(v->user_prop);
	free(hvp);
}*/
static double cotan(double *p0,double* p1,double* p2)
{
	double q1[3],q2[3];
	for(int i=0;i<3;i++)
	{
			q1[i]=p1[i]-p0[i];
			q2[i]=p2[i]-p0[i];
	}

	double temp=0,re=0;
	temp=q1[1]*q2[2]-q1[2]*q2[1];
	re+=temp*temp;
	temp=-q1[0]*q2[2]+q1[2]*q2[0];
	re+=temp*temp;
	temp=q1[0]*q2[1]-q1[1]*q2[0];
	re+=temp*temp;
	re=sqrt(re);
	temp=0;
	for(int i=0;i<3;i++)
	{
		temp+=q1[i]*q2[i];
	}	

	return temp/re; 
}
double get_cotan_from_halfface(Mesh* mesh,template_hf* hf)
{
	double re=0;
	template_c* c=hf->cell;
	if(c==NULL)
	{
		return re;
	}
	template_v* vs[3];
	for(int i=0;i<c->vertices_size;i++)
	{
		bool flag=false;
		for(int j=0;j<hf->vertices_size;j++)
		{
			if(c->vertices[i]==hf->vertices[j])
			{
					flag=true;
			}
		}
		if(!flag)
		{
			vs[0]=c->vertices[i];
			break;
		}
	}
	vs[1]=hf->vertices[0];vs[2]=hf->vertices[1];

	re=cotan(vs[0]->point,vs[1]->point,vs[2]->point); 
	return re;
}



class Harmonic
{
public:
	Harmonic(Mesh *m):mesh(m)
	{
		tas=(Tensors_Algebra_System*)malloc(sizeof(Tensors_Algebra_System));
		//Tensors_Algebra_System_double_init(tas,10);
		//int sum=0;
		for(auto vit=mesh->vertices.begin();vit!=mesh->vertices.end();vit++)
		{
			if(vit->second->user_prop!=NULL)
			{
				free(vit->second->user_prop);
				vit->second->user_prop=NULL;
			}
			vit->second->user_prop=malloc(sizeof(Ham_v_prop));
		}
		for(auto fit=mesh->faces.begin();fit!=mesh->faces.end();fit++)
		{
			if(fit->second->user_prop!=NULL)
			{
				free(fit->second->user_prop);
				fit->second->user_prop=NULL;
			}
			fit->second->user_prop=malloc(sizeof(Ham_f_prop));
		}
		//mesh->init_v_prop=Ham_v_prop_init;
		//mesh->free_v_prop=Ham_v_prop_free;
		LL=NULL;
		//Tensors_Algebra_System_double_init(tas,)
	}
	~Harmonic()
	{
		if(LL!=NULL)
		{
			tas->T_free(tas,LL);
		}
		free(tas);
	}
	void adjust_standard_loc()
	{
		auto vit=mesh->vertices.begin();
		double mins[3]={vit->second->point[0],vit->second->point[1],vit->second->point[2]},maxs[3]={mins[0],mins[1],mins[2]};
		for(;vit!=mesh->vertices.end();vit++)
		{
			for(int i=0;i<3;i++)
			{
				if(vit->second->point[i]>maxs[i])
				{
					maxs[i]=vit->second->point[i];

				}
				if(vit->second->point[i]<mins[i])
				{
					mins[i]=vit->second->point[i];
				}
			}
		}

		double scales[3];		
		for(int i=0;i<3;i++)
		{
		//	printf("%lf %lf\n",mins[i],maxs[i]);
			center[i]=(maxs[i]+mins[i])/2.0;
			scales[i]=0.99/(maxs[i]-mins[i]);
		}
		scale=scales[2];
		if(scales[0]<scales[1]&&scales[0]<scales[2])
		{
			scale=scales[0];
		}
		else if(scales[1]<scales[0]&&scales[1]<scales[2])
		{
			scale=scales[1];
		}
		for(vit=mesh->vertices.begin();vit!=mesh->vertices.end();vit++)
		{
			for(int i=0;i<3;i++)
			{

				vit->second->point[i]=(vit->second->point[i]-center[i])*scale+0.5;
			
			}
		}	

	}
	void init()
	{
		int i=0,j=0;
		Ham_v_prop * temp=NULL;
		for(auto vit=mesh->vertices.begin();vit!=mesh->vertices.end();vit++)
		{
			temp=(Ham_v_prop*)(vit->second->user_prop);
			temp->ux=0;temp->uy=0;
			if(mesh->vertex_is_boundary(mesh,*(vit->second)))
			{
				temp->index=j;
				j++;
			}
			else
			{
				temp->index=i;
				i++;
			}				
		}
		//printf("%d %d \n",i,j);	
		//Tensors_Algebra_System_double_init(tas,i);
		H.resize(i,i);
		H.setZero();
	}
	void init_weight()
	{
		Ham_f_prop* temp1=NULL;
		for(auto fit=mesh->faces.begin();fit!=mesh->faces.end();fit++)
		{
			temp1=(Ham_f_prop*)(fit->second->user_prop);
			temp1->weight=0;
			for(int i=0;i<2;i++)
			{
				if(fit->second->halffaces[i]->cell!=NULL)
				{
					temp1->weight+=get_cotan_from_halfface(mesh,fit->second->halffaces[i]);	
				}	
			}
		}
	}
	void construct_LLE()
	{
		template_v* vs[2];
		Eigen::MatrixXd v(H.rows(),mesh->num_v(mesh)-H.rows());
		v.setZero();
		LLE=v;
		Ham_v_prop* temp=NULL,*temp0=NULL;Ham_f_prop*temp1=NULL;	
		for(auto fit=mesh->faces.begin();fit!=mesh->faces.end();fit++)
		{
			if(!mesh->face_is_boundary(mesh,*(fit->second)))
			{
				vs[0]=fit->second->vertices[0];
				vs[1]=fit->second->vertices[1];
				temp1=(Ham_f_prop*)(fit->second->user_prop);
				for(int i=0;i<2;i++)
				{
					if(!mesh->vertex_is_boundary(mesh,*vs[i]))
					{
						temp=(Ham_v_prop*)(vs[i]->user_prop);
						temp0=(Ham_v_prop*)(vs[(i+1)%2]->user_prop);
						H.coeffRef(temp->index,temp->index)+=temp1->weight;
						if(!mesh->vertex_is_boundary(mesh,*vs[(i+1)%2]))
						{
							H.coeffRef(temp->index,temp0->index)=-temp1->weight;
						}
						else
						{
							v.coeffRef(temp->index,temp0->index)=temp1->weight;
						}
					}
				}	
			}
		}	
		//printf("here\n");
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>llt;
		llt.compute(H);
		//printf("here\n");
		for(int j=0;j<v.cols();j++)
		{
			LLE.col(j)=llt.solve(v.col(j));
		}

	}
	void construct_LL()
	{
		//Tensor* t=tas->T_create();
		template_v* vs[2];
		int rows=tas->as->elements->size,cols=mesh->num_v(mesh)-rows;
		double *A,*B;
		A=(double*)malloc(sizeof(double)*rows*rows);
		memset(A,0,sizeof(double)*rows*rows);
		B=(double*)malloc(sizeof(double)*rows*cols);
		memset(B,0,sizeof(double)*rows*cols);
		Ham_v_prop* temp=NULL,*temp0=NULL;Ham_f_prop*temp1=NULL;	
		for(auto fit=mesh->faces.begin();fit!=mesh->faces.end();fit++)
		{
			if(!mesh->face_is_boundary(mesh,*(fit->second)))
			{
				vs[0]=fit->second->vertices[0];
				vs[1]=fit->second->vertices[1];
				temp1=(Ham_f_prop*)(fit->second->user_prop);
				for(int i=0;i<2;i++)
				{
					if(!mesh->vertex_is_boundary(mesh,*vs[i]))
					{
						temp=(Ham_v_prop*)(vs[i]->user_prop);
						temp0=(Ham_v_prop*)(vs[(i+1)%2]->user_prop);
						A[temp->index*rows+temp->index]+=temp1->weight;
						if(!mesh->vertex_is_boundary(mesh,*vs[(i+1)%2]))
						{
							A[temp->index*rows+temp0->index]=-temp1->weight;
						}
						else
						{
							B[temp->index*cols+temp0->index]=temp1->weight;
						}
						//int k=vs[i]->
					}
				}
			
			}
		}	
		Tensor* t1,*t2,*t3;
		t1=tas->T_create_matrix(tas,A,rows,rows);
		t3=Tensor_inverse(tas,t1);
		t2=tas->T_create_matrix(tas,B,rows,cols);
		LL=W_Tensor_Contraction(tas,t3,t2,1,0);

		tas->T_free(tas,t1);tas->T_free(tas,t2);tas->T_free(tas,t3);

		free(A);free(B);
	}
	void map_squre()
	{
		mesh->external_cell_init_(mesh);

		int size=node_size(mesh->external_cell.halffaces);
		int step=size/4;
		double dir[4][2];
		dir[0][0]=1;dir[0][1]=0;	
		dir[1][0]=0;dir[1][1]=1;	
		dir[2][0]=-1;dir[2][1]=0;	
		dir[3][0]=0;dir[3][1]=-1;
		int i=1,k=0;

	//	printf("%d\n",step);
		template_v* v0=NULL,*v1=NULL;
		Ham_v_prop* temp=NULL,*temp1=NULL;
		int steps[4];
		steps[0]=step;steps[1]=2*step;steps[2]=3*step;steps[3]=size;
		for(Node* nit=mesh->external_cell.halffaces;nit!=NULL;nit=(Node*)(nit->Next))
		{
			v0=((template_hf*)(nit->value))->vertices[1];
			v1=((template_hf*)(nit->value))->vertices[0];
			//v=(template_v*)(nit->value);
			temp=(Ham_v_prop*)(v0->user_prop);
			temp1=(Ham_v_prop*)(v1->user_prop);
			
			if(k==3)
			{
				double temp_step=size-step*3;
				temp1->ux=temp->ux+dir[k][0]/(double)temp_step;
				temp1->uy=temp->uy+dir[k][1]/(double)temp_step;
			}
			else
			{
				temp1->ux=temp->ux+dir[k][0]/(double)step;

				temp1->uy=temp->uy+dir[k][1]/(double)step;
			}
		
			//printf("v1:%d v2:%d,ux:%lf uy:%lf\n",v0->id,v1->id,temp1->ux,temp1->uy);
			if(i==steps[0]||i==steps[1]||i==steps[2]||i==steps[3])
			{
				k++;
			}
			i++;
		}	
		temp1->ux=0;temp1->uy=0;
	}
	
	void map()
	{
		Eigen::VectorXd ux(LLE.cols()),uy(LLE.cols()),x(LLE.cols()),y(LLE.cols());
		Ham_v_prop* temp=NULL;
		for(auto iter=mesh->vertices.begin();iter!=mesh->vertices.end();iter++)
		{
			if(mesh->vertex_is_boundary(mesh,*(iter->second)))
			{
				temp=(Ham_v_prop*)(iter->second->user_prop);
				ux.coeffRef(temp->index)=temp->ux;
				uy.coeffRef(temp->index)=temp->uy;

			}
		}	
		x=LLE*ux;y=LLE*uy;

		for(auto iter=mesh->vertices.begin();iter!=mesh->vertices.end();iter++)
		{
			if(!mesh->vertex_is_boundary(mesh,*(iter->second)))
			{
				temp=(Ham_v_prop*)(iter->second->user_prop);
				temp->ux=x.coeff(temp->index);
				//printf("ux:%lf\n",temp->ux);
				temp->uy=y.coeff(temp->index);

			}
		}	

	}
	void get_data_uv(Mesh* m,float**v,unsigned int**f)
	{
#ifdef MANIFOLD_REQUIRE
    	if(m->dimension!=2)
    	{
        	printf("can't draw\n");
        	return ;
    	}
#endif
    	RB_Tree* tree=(RB_Tree*)malloc(sizeof(RB_Tree));
    	RB_Tree_init_int(tree);
    	RB_int rbt,*rbt1=NULL;
    //要重排vertices id
    	if((*v)!=0)
    	{
        	free((*v));
    	}
    	*v=(float*)malloc(sizeof(float)*m->num_v(m)*3);
    	memset(*v,0,sizeof(float)*m->num_v(m)*3);
    	int f_size=0;
    	for(auto iter=m->cells.begin();iter!=m->cells.end();iter++)
    	{	
        	f_size+=((iter->second)->vertices_size+1);  
    	}

		if((*f)!=0)
		{free((*f));
		}
    //int *temp_v_id=(int*)malloc(sizeof(int)*m->num_v(m));
    	*f=(unsigned int*)malloc(sizeof(unsigned int)*f_size);
    	memset(*f,0,sizeof(unsigned int)*f_size);
    	int i=0;
    	Ham_v_prop* temp=NULL;	
    	int * temp_ids=(int*)malloc(sizeof(int)*m->num_v(m));
    	for(auto iter=m->vertices.begin();iter!=m->vertices.end();iter++)
    	{       
    		temp=(Ham_v_prop*)(iter->second->user_prop);
    		(*v)[i*3+0]=temp->ux;
    		(*v)[i*3+1]=temp->uy;
        	//int * temp_id=(int*)malloc(sizeof(int));
        	temp_ids[i]=i;
        	rbt.key=iter->second->id;
        	rbt.value=(void*)(&temp_ids[i]);
        	tree->insert(tree,&rbt);
        //iter->second->prop=(void*)(&temp_v_id[i]);
//       (*v)[i*3+2]=(float)temp_z/2.0;
        	i++; 
    	}
    	i=0;
    	for(auto iter=m->cells.begin();iter!=m->cells.end();iter++)
    	{
        	(*f)[i]=iter->second->vertices_size;
        	i++;
        	for(auto iter1=m->cv_begin(m,*(iter->second));iter1!=m->cv_end(m,*(iter->second));iter1++)
        	{
        		rbt.key=(*iter1).id;
        		rbt1=(RB_int*)(tree->find(tree,&rbt));
            	(*f)[i]=*((int*)((rbt1->value)));
            	i++;
        	}

    	}
    	free(temp_ids);
    	RB_Tree_free(tree); 
   // free(temp_v_id);
    //reset_v_prop(m);
	}
	void get_edges_data(Mesh* m,float**v,unsigned int ** f)
	{
#ifdef MANIFOLD_REQUIRE
    	if(m->dimension!=2)
    	{
        	printf("can't draw\n");
        	return ;
    	}
#endif
    	RB_Tree* tree=(RB_Tree*)malloc(sizeof(RB_Tree));
    	RB_Tree_init_int(tree);
    	RB_int rbt,*rbt1=NULL;
    	int *temp_ids=(int*)malloc(sizeof(int)*m->num_v(m));

    	if((*v)!=0)
    	{
        	free((*v));
    	}
    	*v=(float*)malloc(sizeof(float)*m->num_v(m)*3);
    	memset(*v,0,sizeof(float)*m->num_v(m)*3);
    	int f_size=m->num_f(m);

		if((*f)!=0)
		{free((*f));
		}
		*f=(unsigned int*)malloc(sizeof(unsigned int)*f_size*2);
    	memset(*f,0,sizeof(unsigned int)*f_size*2);
    	int i=0;
    	Ham_v_prop* temp=NULL;	
    	for(auto iter=m->vertices.begin();iter!=m->vertices.end();iter++)
    	{
    		temp=(Ham_v_prop*)(iter->second->user_prop);

    		
    		(*v)[i*3+0]=temp->ux;
    		(*v)[i*3+1]=temp->uy;
    	
			temp_ids[i]=i;
        	rbt.key=iter->second->id;
        	rbt.value=(void*)(&temp_ids[i]);
        	tree->insert(tree,&rbt);
        	i++;
    	}
    	i=0;
    	for(auto iter=m->faces.begin();iter!=m->faces.end();iter++)
    	{
    		for(int j=0;j<2;j++)
    		{
        		rbt.key=iter->second->vertices[j]->id;
        		rbt1=(RB_int*)(tree->find(tree,&rbt));
            	(*f)[i]=*((int*)((rbt1->value)));
            	i++;
    		}

    	}
    	free(temp_ids);
    	RB_Tree_free(tree);	
    	//me->Data_index_rows=mesh->num_f(mesh);
    
    	//me->Data_rows=me->Data_index_rows*2; 
	}

public:
	Mesh *mesh;
	Tensors_Algebra_System *tas;
	Tensor *LL;
	Eigen::SparseMatrix<double> H;
	Eigen::MatrixXd V,LLE,F;
	double scale;
	double center[3];
// Tensors_Algebra_System*tas=(Tensors_Algebra_System*)malloc(sizeof(Tensors_Algebra_System));
 //   Tensors_Algebra_System_mpf_init(tas,6);//一般tas的维数要大一位
};
class Geometry_image
{
public:
	Geometry_image(Mesh *m):mesh(m)
	{
		image_tensor=NULL;
		test_data=NULL;
		image_tensor1=NULL;
		image_tensor2=NULL;
	}
	~Geometry_image()
	{
		if(image_tensor!=NULL)
		{
			free(image_tensor);
		}
		if(test_data!=NULL)
		{
			free(test_data);
		}
		if(image_tensor1!=NULL)
		{
			free(image_tensor1);
		}
		if(image_tensor2!=NULL)
		{
			free(image_tensor2);
		}
	}
	void WriteBMP(unsigned char*img,const char* filename,int w,int h)
	{
#ifdef WIN32
    FILE *fp;
    fopen_s(&fp,filename,"r");
#else
    FILE *fp=fopen(filename,"wb");
#endif
    	int l=(w*3+3)/4*4;
     	int bmi[]= {l*h+54,0,54,40,w,h,1|3*8<<16,0,l*h,0,0,100,0};
     	//FILE *fp = fopen(filename,"wb");
     	fprintf(fp,"BM");
     	fwrite(&bmi,52,1,fp);
     	fwrite(img,1,l*h,fp);
     	fclose(fp);
	}
	void double_image(int w,int h)
	{
		image_tensor1=(unsigned char*)malloc(sizeof(unsigned char)*w*h*3);
		image_tensor2=(unsigned char*)malloc(sizeof(unsigned char)*w*h*3);	
		for(int i=0;i<w;i++)
		{
			for(int j=0;j<h;j++)
			{
				for(int k=0;k<3;k++)
				{
					int temp=round(test_data[(j*w+i)*3+k]*(255*255-1));
					int w1=temp/255;
					int w2=temp-w1*255;
					image_tensor1[(j*w+i)*3+k]=w1;
					image_tensor2[(j*w+i)*3+k]=w2;
//					int ww1=image_tensor1[(j*w+i)*3+k];
//					int ww2=image_tensor2[(j*w+i)*3+k];
//					//printf("%lf %lf\n",test_data[(j*w+i)*3+k],(ww1*255+ww2)/(double)(255*255-1));				
				}
			}
		}
	}
	/*void WriteBMP1(unsigned char*img,const char* filename,int w,int h)
	{

    	int l=(w*3+3)/4*4;
     	int bmi[]= {l*h+54,0,54,40,w,h,1|3*8<<16,0,l*h,0,0,100,0};
     	FILE *fp = fopen(filename,"wb");
     	fprintf(fp,"BM");
     	fwrite(&bmi,52,1,fp);
     	fwrite(img,1,l*h,fp);
     	fclose(fp);
	}*/
	void double_image2mesh(unsigned char *data1,unsigned char *data2,int w,int h,Mesh* m)
	{

		Mesh_init(m);
		m->simplex=1;
#ifdef MANIFOLD_REQUIRE
		m->dimension=2;
#endif
		double p[3];
		printf("hereli %d %d \n",w,h);
		template_v**indexs=(template_v**)malloc(sizeof(template_v*)*w*h);
		printf("hereli\n");
		for(int i=0;i<w;i++)
		{
			for(int j=0;j<h;j++)
			{
				for(int k=0;k<3;k++)
				{
					p[k]=(data1[(j*w+i)*3+k]*255+data2[(j*w+i)*3+k])/(double)(255.0*255.0-1.0);
					//p[k]=test_data[(j*w+i)*3+k];
				}
				template_v*v=m->create_vertexv(m,p,3);
				indexs[j*w+i]=v;	
			}
		}
		template_v* vs[3];
		for(int i=0;i<w-1;i++)
		{
			for(int j=0;j<h-1;j++)
			{
				vs[0]=indexs[j*w+i];
				vs[2]=indexs[(j+1)*w+i];
				vs[1]=indexs[(j+1)*w+i+1];

				m->create_cellv(m,vs,3);
				vs[0]=indexs[j*w+i];
				vs[2]=indexs[(j+1)*w+i+1];
				vs[1]=indexs[j*w+i+1];
				m->create_cellv(m,vs,3);	
			}
		}
		free(indexs);	
	}
	void  image2mesh(unsigned char * data,int w,int h,Mesh* m)
	{
		Mesh_init(m);
		m->simplex=1;
#ifdef MANIFOLD_REQUIRE
		m->dimension=2;
#endif
		double p[3];
		printf("hereli %d %d \n",w,h);
		template_v**indexs=(template_v**)malloc(sizeof(template_v*)*w*h);
		printf("hereli\n");
		for(int i=0;i<w;i++)
		{
			for(int j=0;j<h;j++)
			{
				for(int k=0;k<3;k++)
				{
						//int w1=image_tensor1[(j*w+i)*3+k];
					//int w2=image_tensor2[(j*w+i)*3+k];
					//p[k]=(w1*255+w2)/(double)(255*255-1);
					p[k]=data[(j*w+i)*3+k]/255.0;
					//p[k]=test_data[(j*w+i)*3+k];
				}
				template_v*v=m->create_vertexv(m,p,3);
				indexs[j*w+i]=v;	
			}
		}
		template_v* vs[3];
		for(int i=0;i<w-1;i++)
		{
			for(int j=0;j<h-1;j++)
			{
				vs[0]=indexs[j*w+i];
				vs[2]=indexs[(j+1)*w+i];
				vs[1]=indexs[(j+1)*w+i+1];

				m->create_cellv(m,vs,3);
				vs[0]=indexs[j*w+i];
				vs[2]=indexs[(j+1)*w+i+1];
				vs[1]=indexs[j*w+i+1];
				m->create_cellv(m,vs,3);	
			}
		}
		free(indexs);	

	}
	void regularize_uv()
	{
		double max_x=0,min_x=1,max_y=0,min_y=1;
		Ham_v_prop* temp=NULL;		
		for(auto it=mesh->vertices.begin();it!=mesh->vertices.end();it++)
		{

				temp=(Ham_v_prop*)(it->second->user_prop);
				if(temp->uy>max_y)
				{
					max_y=temp->uy;
				}
				if(temp->uy<min_y)
				{
					min_y=temp->uy;
				}
				if(temp->ux>max_x)
				{
					max_x=temp->ux;
				}
				if(temp->ux<min_x)
				{
					min_x=temp->ux;
				}
		}
		printf("%lf %lf %lf %lf\n",max_x,min_x,max_y,min_y);	

	}
	//边界调整
	void map_squre(int w,int h)
	{
		mesh->external_cell_init_(mesh);
		template_v* v0=NULL,*v1=NULL,*vv1=NULL,*vv2=NULL;	
		Ham_v_prop* temp=NULL,*temp1=NULL;
		int x1,x2,y1,y2;	
		for(Node* nit=mesh->external_cell.halffaces;nit!=NULL;nit=(Node*)(nit->Next))
		{
			v0=((template_hf*)(nit->value))->vertices[1];
			v1=((template_hf*)(nit->value))->vertices[0];
			//v=(template_v*)(nit->value);
			temp=(Ham_v_prop*)(v0->user_prop);
			temp1=(Ham_v_prop*)(v1->user_prop);
			x1=round(temp->ux*(w-1));x2=round(temp1->ux*(w-1));
			y1=round(temp->uy*(h-1));y2=round(temp1->uy*(h-1));
			double loc=0;
			double lambda=0;
			if(x1!=x2)
			{
				vv1=v0;vv2=v1;
				int in1=x1,in2=x2;
				if(x2>x1)
				{

				}
				else
				{
					vv1=v1;vv2=v0;
					in1=x2,in2=x1;
				}

				for(int i=in1;i<=in2&&i<w;i++)
				{
					temp=(Ham_v_prop*)(vv1->user_prop);temp1=(Ham_v_prop*)(vv2->user_prop);

					lambda=(i/(double)(w-1)-temp->ux)/(temp1->ux-temp->ux);
					if(lambda>=0&&lambda<=1)
					{

						for(int k=0;k<3;k++)
						{
							loc=(lambda*vv2->point[k]+(1-lambda)*vv1->point[k]);
							int temp_i=round(loc*254.0);
						//if(fabs(loc-test_data[(y1*w+i)*3+k])>0.1)
						//{
							/*if(y1==4||y1==5)
							{
								printf("i:%d j:%d %lf %lf %d %d %lf %lf\n",i,y1,loc,test_data[(y1*w+i)*3+k],in1,in2,vv1->point[k],vv2->point[k]);
							}*/
						//}
							test_data[(y1*w+i)*3+k]=loc;
							image_tensor[(y1*w+i)*3+k]=temp_i;
						}
					}	
					//image_tensor[(y1*w+i)*3];
					//printf("\n");
				}
			}
			if(y1!=y2)
			{
				vv1=v0;vv2=v1;
				int in1=y1,in2=y2;
				if(y2>y1)
				{

				}
				else
				{
					vv1=v1;vv2=v0;
					in1=y2,in2=y1;
				}
				for(int i=in1;i<=in2&&i<h;i++)
				{
					temp=(Ham_v_prop*)(vv1->user_prop);temp1=(Ham_v_prop*)(vv2->user_prop);

					lambda=(i/(double)(h-1)-temp->uy)/(temp1->uy-temp->uy);
					if(lambda>=0&&lambda<=1)
					{


						for(int k=0;k<3;k++)
						{
							loc=(lambda*vv2->point[k]+(1-lambda)*vv1->point[k]);
							int temp_i=round(loc*254.0);
						//if(fabs(loc-test_data[(i*w+x1)*3+k])>0.1)
						//{
							/*if(i==4)
							{
								printf("%lf %lf\n",i/(double)(h-1),temp->uy);
								printf("i:%d j:%d %lf %lf %d %d %lf %lf %lf\n",x1,i,loc,test_data[(i*w+x1)*3+k],in1,in2,vv1->point[k],vv2->point[k],lambda);
							}*/
						//}
							test_data[(i*w+x1)*3+k]=loc;
							image_tensor[(i*w+x1)*3+k]=temp_i;
						}	
					}
					//printf("\n");
					//image_tensor[(y1*w+i)*3];
				}

			}
			//printf("%d %d %d %d\n",x1,x2,y1,y2);	
			
		}	
		//printf("%d\n",(int)(0.999));

	}
	void regular_sample(int w,int h)
	{
		regularize_uv();
		image_tensor=(unsigned char*)malloc(sizeof(unsigned char)*w*h*3);
		test_data=(float*)malloc(sizeof(float)*w*h*3);
		memset(test_data,0,sizeof(float)*w*h*3);
		memset(image_tensor,0,sizeof(unsigned char)*w*h*3);
//		printf("libo\n");
		double max_y=0,min_y=1;
		double max_x=0,min_x=1;
		Eigen::MatrixXd A(3,3);
		A.coeffRef(2,0)=1;A.coeffRef(2,1)=1;A.coeffRef(2,2)=1;	
		Eigen::VectorXd b(3);
		b.coeffRef(2)=1;
		Ham_v_prop* temp=NULL;
		int max_y_,min_y_,max_x_,min_x_;
		for(auto cit=mesh->cells.begin();cit!=mesh->cells.end();cit++)
		{
			//printf("begin\n");
			max_y=0;min_y=1;
			max_x=0;min_x=1;
			int max_y_index=0,min_y_index=0;
			int max_x_index=0,min_x_index=0;
			for(int j=0;j<cit->second->vertices_size;j++)
			{
				temp=(Ham_v_prop*)(cit->second->vertices[j]->user_prop);
				if(temp->uy>max_y)
				{
					max_y_index=j;
					max_y=temp->uy;
				}
				if(temp->uy<min_y)
				{
					min_y_index=j;
					min_y=temp->uy;
				}
				if(temp->ux>max_x)
				{
					max_x_index=j;
					max_x=temp->ux;
				}
				if(temp->ux<min_x)
				{
					min_x_index=j;
					min_x=temp->ux;
				}
				A.coeffRef(0,j)=temp->ux;
				A.coeffRef(1,j)=temp->uy;
				//printf("%lf %lf\n",temp->ux,temp->uy);
			}
			min_y_=((int)(min_y*h));
			max_y_=((int)(max_y*h));	
			min_x_=((int)(min_x*w));
			max_x_=((int)(max_x*w));	

			//printf("%lf %lf %lf %lf\n",min_y,max_y,min_x,max_x);	
			if(min_y_==max_y_)
			{
				for(int i=0;i<3;i++)
				{
					A.coeffRef(1,i)=0;
					if(i!=min_x_index&&i!=max_x_index)
					{
						A.coeffRef(1,i)=1;
					}
				}
				//A.coeffRef(1,0)=1;A.coeffRef(1,1)=0;A.coeffRef(1,2)=0;
			}	
			else if(min_x_==max_x_)
			{
				//b.coeff(1)=0;
				for(int i=0;i<3;i++)
				{
					A.coeffRef(0,i)=0;
					if(i!=min_y_index&&i!=max_y_index)
					{
						A.coeffRef(0,i)=1;
					}
				}
			}
			{
				Eigen::MatrixXd Ain=A.inverse(),x,y;
				//printf("%d %d %d %d\n",min_y_,max_y_,min_x_,max_x_);		
				for(int i=min_x_;i<=max_x_&&i<w&&i>=0;i++)
				{

					for(int j=min_y_;j<=max_y_&&j<h&&j>=0;j++)
					{
					//	printf("libo\n");
						b.coeffRef(0)=i/(double)(w-1);
						b.coeffRef(1)=j/(double)(h-1);
						//bool flag= false;
						if(min_y_==max_y_)
						{
						
							b.coeffRef(1)=0;
							//x=Ain*b;
							//y=A*x;
						//	flag=true;
						}
						else if(min_x_==max_x_)
						{
							b.coeffRef(0)=0;
						//	flag=true;
						}
						else
						{

						}
					//	printf("b %lf %lf\n",b.coeff(0),b.coeff(1));	
						x=Ain*b;
						y=A*x;
					//	printf("i:%d j:%d jiance x:%lf %lf %lf b:%lf %lf,A:%lf %lf %lf %lf\n",i,j,x.coeff(0),x.coeff(1),x.coeff(2),b.coeff(0),b.coeff(1),A.coeff(1,0),A.coeff(1,1),A.coeff(1,2),A.coeff(1,0));
						if(x.coeff(0)>=-0.00001&&x.coeff(1)>=-0.00001&&x.coeff(2)>=-0.00001)
						{
					//		printf("here\n");
							double loc[3]={0,0,0};
							for(int k=0;k<3;k++)
							{
								loc[0]+=x.coeff(k)*cit->second->vertices[k]->point[0];
								loc[1]+=x.coeff(k)*cit->second->vertices[k]->point[1];	
								loc[2]+=x.coeff(k)*cit->second->vertices[k]->point[2];
							}
							//printf("here\n");
							for(int k=0;k<3;k++)
							{
								int temp_i=round(loc[k]*254);

			//					printf("i:%d j%d libo%d\n",i,j,(j*w+i)*3+k);
								if((j*w+i)*3+k>=30000)
								{
									printf("hauile\n");
									return;
								}
								test_data[(j*w+i)*3+k]=loc[k];
								//printf("%lf %d ",loc[k],temp_i);
								image_tensor[(j*w+i)*3+k]=temp_i;
								//printf("%d ",image_tensor[(i*w+j)*3+k]);	
							}
				//			printf("\n");
						}							
					}
				
				}
				//printf("end\n");
			}
		}
	//解决没有扫描的点
		/*int sum1=0;
		for(int i=0;i<h;i++)
		{
			for(int j=0;j<w;j++)
			{

				if(image_tensor[(i*w+j)*3+0]==0&&image_tensor[(i*w+j)*3+1]==0&&image_tensor[(w*i+j)*3+2]==0)
				{
					printf("i:%d j:%d\n",i,j);
					sum1++;
					double sum=0;
					if(i+1<h)
					{
						for(int k=0;k<3;k++)
						{
							image_tensor[(i*w+j)*3+k]+=image_tensor[((i+1)*w+j)*3+k];
						}
						sum++;
					}
					if(i-1>=0)
					{
						for(int k=0;k<3;k++)
						{
							image_tensor[(i*w+j)*3+k]+=image_tensor[((i-1)*w+j)*3+k];
						}
						sum++;
					}
					if(j+1<w)
					{
						for(int k=0;k<3;k++)
						{
							image_tensor[(i*w+j)*3+k]+=image_tensor[((i)*w+j+1)*3+k];
						}
						sum++;
					}
					if(j-1>=0)
					{
						for(int k=0;k<3;k++)
						{
							image_tensor[(i*w+j)*3+k]+=image_tensor[((i)*w+j-1)*3+k];
						}
						sum++;
					}
					for(int k=0;k<3;k++)
					{
						//printf("sum:%lf\n",sum);
						image_tensor[(i*w+j)*3+k]=image_tensor[(i*w+j)*3+k]/(double)sum;
					}
				}
			}
			
		}	
		
		printf("sum1:%d\n",sum1);	*/
	//解决边界不对称

			
		map_squre(w,h);
		double_image(w,h);	

	}

	float *test_data;
	
//protected:
	unsigned char *image_tensor;
	unsigned char *image_tensor1,*image_tensor2;
	Mesh* mesh;
};



#endif
