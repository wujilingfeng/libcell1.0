//#undef LibCell_IO
//一次性读取文件到内存会提高速度
#include<Mesh_IO/Mesh_IO.h>
//#include<sys/sysinfo.h>
//num[3]
void _ReadArray_(template_m*m,double**V,int**F,int**C,int*num,int back_dim,int manifold_dim)
{
#ifdef MANIFOLD_REQUIRE
    m->dimension=manifold_dim;
#endif
	int n=back_dim;
    for(int i=0;i<num[0];i++)
    {
        template_v* v=m->create_vertex(m);
	    v->point_size=n;
        v->point=(double*)malloc(sizeof(double)*n);
        for(int j=0;j<n;j++)
	    {   
            v->point[j]=V[i][j];   
	    }
    }
    int cols;
	if(F==NULL||num[1]==0)
	{
		m->simplex=1;
        for(int i=0;i<num[2];i++)
        {    
            cols=C[i][0];
            template_v** temp_v=(template_v**)malloc(sizeof(template_v*)*cols);
            //template_v* temp_v[cols];
            for(int j=0;j<cols;j++)
	        {
                int id=C[i][j+1];
		        auto it=m->vertices.find(id);
		        template_v* v;
		        if(it!=m->vertices.end())
		        {
		            v=it->second;
		        }
		        else
                {
                    printf("can t find this vertex: %d\r\n",id);
		            return;
		        }
		        temp_v[j]=v;
	        }
	        m->create_cellv(m,temp_v,cols);
            free(temp_v);

            
        }

	}
	else
	{
        m->simplex=0;
        for(int i=0;i<num[1];i++)
        {
            cols=F[i][0];
            template_v**temp_v=(template_v**)malloc(sizeof(template_v*)*cols);
            //template_v* temp_v[cols];
                
            for(int j=0;j<cols;j++)
            {
                int id=F[i][j+1];
                auto it=m->vertices.find(id);
                template_v* v;
                if(it!=m->vertices.end())
                {
                    v=it->second;
                    
                }
                else
                {
                    printf("can t find this vertex\r\n");
                    return ;
                    
                }
                temp_v[j]=v;
                    
                
                
            }
            template_f* temp_f=m->create_facev(m,temp_v,cols);
            m->create_halfface(m,temp_f,temp_v,cols);  
            free(temp_v);
        }
        for(int i=0;i<num[2];i++)
        {
                cols=C[i][0];
                template_hf**temp_hf=(template_hf**)malloc(sizeof(template_hf*)*cols); 
                //template_hf* temp_hf[cols];
                for(int j=0;j<cols;j++)
                {
                    int id=C[i][j+1];
                    auto it=m->halffaces.find(id);
                    template_hf * hf;
                    if(it!=m->halffaces.end())
                    {
                        hf=it->second;
                    }
                    else
                    {
                        printf("can t find this halfface\n");
                        return ;
                    
                    }
                
                    temp_hf[j]=hf;
                
                
                }
                m->create_cellf(m,temp_hf,cols);
                free(temp_hf);
        }

    
    }
}

void _ReadCell_(template_m*m,char const*filename)
{
#ifdef WIN32
    FILE *infile;
    fopen_s(&infile,filename,"r");
#else
    FILE *infile=fopen(filename,"r");
#endif
    if(!infile)
	{
	    printf("cant't open this file\r\n");
	    return;
	}
    fseek(infile,0,SEEK_SET);
    char off[20];
    fscanf(infile,"%s",off);
    if(strcmp("CELL",off)==0)
    {
        //printf("cell read\n");
        fscanf(infile,"%s",off);
        int n;
        fscanf(infile,"%d",&n);
        fscanf(infile,"%s",off);
        int manifold_dim;
        fscanf(infile,"%d",&manifold_dim);
#ifdef MANIFOLD_REQUIRE
        m->dimension=manifold_dim;
#endif
        fscanf(infile,"%s",off);
        int simplex;
        int num[3];
        fscanf(infile,"%d",&simplex);m->simplex=simplex;
        fscanf(infile,"%d %d %d",&num[0],&num[1],&num[2]);
        for(int i=0;i<num[0];i++)
        {
                template_v* v=m->create_vertex(m);
	            v->point_size=n;
                v->point=(double*)malloc(sizeof(double)*n);
                for(int j=0;j<n;j++)
	            {   
                    fscanf(infile,"%lf",&(v->point[j]));
	            }


        
        }
        if(simplex==1)
        {
            int cols=0,id=0;
            for(int i=0;i<num[2];i++)
            {
                fscanf(infile,"%d",&cols);
                template_v** temp_v=(template_v**)malloc(sizeof(template_v*)*cols); 
                //template_v* temp_v[cols];
                for(int j=0;j<cols;j++)
	            {
                    fscanf(infile,"%d",&id);
		            auto it=m->vertices.find(id);
		            template_v* v;
		            if(it!=m->vertices.end())
		            {
		                v=it->second;
		            }
		            else
                    {
                        printf("can t find this vertex: %d\r\n",id);
		                return;
		            }
		            temp_v[j]=v;
	            }
	            m->create_cellv(m,temp_v,cols);
                free(temp_v);
            
            }
             
        }
        else
        {   
		    m->simplex=0;

//		printf("simplex:%d\n",simplex);
            int cols,id;
            for(int i=0;i<num[1];i++)
            {
                fscanf(infile,"%d",&cols);
                //template_v* temp_v[cols];
                template_v**temp_v=(template_v**)malloc(sizeof(template_v*)*cols);
                for(int j=0;j<cols;j++)
                {
                    fscanf(infile,"%d",&id);
                    auto it=m->vertices.find(id);
                    template_v* v;
                    if(it!=m->vertices.end())
                    {
                        v=it->second;
                    
                    }
                    else
                    {
                        printf("can t find this vertex\r\n");
                        return ;
                    
                    }
                    temp_v[j]=v;
                    
                
                
                }
                template_f* temp_f=m->create_facev(m,temp_v,cols);
                m->create_halfface(m,temp_f,temp_v,cols); 
                free(temp_v);
            
            
            
            }
            for(int i=0;i<num[2];i++)
            {
                fscanf(infile,"%d",&cols);
                //temp_hf放到循环外面
                //template_hf* temp_hf[cols];
                template_hf**temp_hf=(template_hf**)malloc(sizeof(template_hf*)*cols);
                for(int j=0;j<cols;j++)
                {
                    fscanf(infile,"%d",&id);
                    auto it=m->halffaces.find(id);
                    template_hf * hf;
                    if(it!=m->halffaces.end())
                    {
                        hf=it->second;
                    }
                    else
                    {
                        printf("can t find this halfface\n");
                        return ;
                    
                    }
                
                    temp_hf[j]=hf;
                
                
                }
                m->create_cellf(m,temp_hf,cols);
                free(temp_hf);
            }
        }    
    }

}
Node* _ReadM_(template_m*m,char const *filename)
{
    m->simplex=1;
#ifdef MANIFOLD_REQUIRE
    m->dimension=2;
#endif
    #ifdef WIN32
    FILE *infile;
    fopen_s(&infile,filename,"r");
#else
    FILE *infile=fopen(filename,"r");
#endif
    if(!infile)
    {
        printf("cant't open this file\r\n");
        return NULL;
    }
    fseek(infile,0,SEEK_SET);
   // char *MM=(char*)malloc(sizeof(char)*20);
   // memset(MM,0,sizeof(char)*20);
    char MM[20];
    char  temp_m[110];
    Node * n_c=NULL;
    RB_Tree* tree=(RB_Tree*)malloc(sizeof(RB_Tree));
    RB_Tree_init_int(tree);
    RB_int rbt,*rbt1=NULL;
    while(fscanf(infile,"%s",MM)>0)
    {
        int temp_id;
        fscanf(infile,"%d",&temp_id);
        if(MM[0]=='#')
        {
            fgets(temp_m,110,infile);
        }
        if(strcmp("Vertex",MM)==0)
        {
            template_v* v=m->create_vertex(m);
            rbt.key=temp_id;
            rbt.value=v;
            tree->insert(tree,&rbt);
            v->point_size=3;
            v->point=(double*)malloc(sizeof(double)*3);
            for(int i=0;i<3;i++)
            {
                fscanf(infile,"%lf",&(v->point[i]));
                
            }
            fgets(temp_m,110,infile);  
            if(temp_m[2]=='r'&&temp_m[3]=='g'&&temp_m[4]=='b')
            {
                double *temp_d=(double*)malloc(sizeof(double)*3);
                n_c=node_overlying(n_c,temp_d); 
                sscanf(&(temp_m[7]),"%lf %lf %lf ",&(temp_d[0]),&(temp_d[1]),&(temp_d[2]));

              //  printf("shi:%lf %lf %lf\n",temp_d[0],temp_d[1],temp_d[2]); 
            }

            //printf("%s\n",temp_m);
        }
        else if(strcmp("Face",MM)==0)
        {
            template_v * vs[3];
            int id=0;
            for(int i=0;i<3;i++)
            {   
                fscanf(infile,"%d",&id);

                rbt.key=id;
                rbt1=(RB_int*)(tree->find(tree,&rbt));
                if(rbt1==NULL)
                {
                    printf("cant find vertex\n");
                    return n_c;
                }
                else
                {
                   vs[i]=(template_v*)(rbt1->value); 
                }
                /*auto it=m->vertices.find(id-1);
                if(it!=m->vertices.end())
                {
                    vs[i]=it->second;
                } 
                else
                {
                    printf("cant find vertex\n");
                    break;
                }*/
                
            }
            
        //    printf("sum%d\n",sum);
            m->create_cellv(m,vs,3);
            //printf("cid:%d\n",c->id);
        }
        else
        {
            break;
        }
        //printf("sum\n");
        //printf("%s \n",MM);
    }
    RB_Tree_free(tree);
    return node_reverse(n_c); 
}
//writecell 
void _WriteCell_(template_m*m,char const* filename)
{
#ifdef WIN32
    FILE *outfile;
    fopen_s(&outfile,filename,"w");
#else
    FILE *outfile=fopen(filename,"w");
#endif
   if(!outfile)
   {
        printf("cant open this file\n");
        return ;
   }
   fprintf(outfile,"CELL\n");
   int back_ground=m->vertices.begin()->second->point_size;
   fprintf(outfile,"background_dim= %d\n",back_ground);
#ifdef MANIFOLD_REQUIRE
   fprintf(outfile,"manifold_dim= %d\n",m->dimension);
#else
   fprintf(outfile,"manifold_dim= 0\n");
#endif
   fprintf(outfile,"simplex= %d\n",m->simplex);
   int num[3];
   num[0]=m->num_v(m);
   num[1]=m->num_hf(m);
   num[2]=m->num_c(m);
    printf("herere libo ***********************************\n"); 
   if(m->simplex==1)
   {

        fprintf(outfile,"%d %d %d\n",num[0],0,num[2]);
        printf("numc_%d\n",m->num_c(m));
        int *temp_v_id=(int*)malloc(sizeof(int)*num[0]);
        int temp_id=0;
        for(auto v_it=m->vertices.begin();v_it!=m->vertices.end();v_it++)
        {
            temp_v_id[temp_id]=temp_id;
            v_it->second->prop=(void*)(&temp_v_id[temp_id]);
            for(int j=0;j<back_ground;j++)
            {
                fprintf(outfile,"%.11lf ",v_it->second->point[j]);
            }
            fprintf(outfile,"\n");
            temp_id++; 
        }
        int i=0;
        for(auto c_it=m->cells.begin();c_it!=m->cells.end();c_it++)
        { 
            fprintf(outfile,"%d ",c_it->second->vertices_size);
            for(auto cv_it=m->cv_begin(m,*(c_it->second));cv_it!=m->cv_end(m,*(c_it->second));cv_it++)
            {

                //printf("i:%d %d  ",i,*((int*)((*cv_it).prop)));
                fprintf(outfile,"%d ",*((int*)((*cv_it).prop))); 
            }
            fprintf(outfile,"\n");
            i++;
        }
        free(temp_v_id);
        reset_v_prop(m);

   }
   else
   {
       int temp_sum=0;
       for(auto c_it=m->cells.begin();c_it!=m->cells.end();c_it++)
       {
           temp_sum+=node_size(c_it->second->halffaces);
            
       }
       num[1]=temp_sum;
        fprintf(outfile,"%d %d %d\n",num[0],num[1],num[2]);
        int *temp_v_id=(int*)malloc(sizeof(int)*num[0]);
        int temp_id=0;
        for(auto v_it=m->vertices.begin();v_it!=m->vertices.end();v_it++)
        {
            temp_v_id[temp_id]=temp_id;
            v_it->second->prop=(void*)(&temp_v_id[temp_id]);
            for(int j=0;j<back_ground;j++)
            {
                fprintf(outfile,"%lf ",v_it->second->point[j]);
            }
            fprintf(outfile,"\n");
            temp_id++; 
        }

        int *temp_hf_id=(int*)malloc(sizeof(int)*num[1]);
        temp_id=0;
        for(auto c_it=m->cells.begin();c_it!=m->cells.end();c_it++)
        {
            for(auto chf_it=m->chf_begin(m,*(c_it->second));chf_it!=m->chf_end(m,*(c_it->second));chf_it++)
            {
                temp_hf_id[temp_id]=temp_id;
                quote(chf_it)->prop=(void*)(&temp_hf_id[temp_id]);
                fprintf(outfile,"%d ",(*chf_it).vertices_size);
                for(auto hfv_it=m->hfv_begin(m,*chf_it);hfv_it!=m->hfv_end(m,*chf_it);hfv_it++)
                {
                    fprintf(outfile,"%d ",*((int*)(quote(hfv_it)->prop))); 
                }

                fprintf(outfile,"\n");
                temp_id++;
            }
        }


        for(auto c_it=m->cells.begin();c_it!=m->cells.end();c_it++)
        {
            int temp_size=node_size(c_it->second->halffaces);
            fprintf(outfile,"%d ",temp_size);
            for(auto chf_it=m->chf_begin(m,*c_it->second);chf_it!=m->chf_end(m,*c_it->second);chf_it++)
            {
                fprintf(outfile,"%d ",*((int*)(quote(chf_it)->prop)));
            }
            fprintf(outfile,"\n");
        }

        free(temp_v_id);
        free(temp_hf_id);
        reset_v_prop(m);
        reset_hf_prop(m);

   }
   

}
//只有off文件描述单形时才有效（我认为off文件还能描述n维流形）
void _ReadOff_(template_m *m,char const * filename, int n)
{
#ifdef WIN32
	FILE* infile;
	fopen_s(&infile,filename,"r");
#else
	FILE* infile=fopen(filename,"r");
#endif
    if(!infile)
	{
	    printf("cant't open this file\r\n");
	    return;
	}
    m->dimension=2;
    m->simplex=1; 
    //fseek(infile,0,SEEK_END);
    //int len=ftell(infile);
    //char*source=(char*)malloc(len+1);
    //fread(source,1,len,infile);
    //free(source);
    fseek(infile,0,SEEK_SET);
	char off[6];
	fscanf(infile,"%s",off);
    //printf("off\r\n");
    if(strcmp("OFF",off)==0)
    {
            
	    int num[3];
        if(fscanf(infile,"%d %d %d",&num[0],&num[1],&num[2])==3)
        {
	        int rows=0,cols=0;
	        while(rows<num[0]&&!feof(infile))
	        {
	            template_v* v=m->create_vertex(m);
	            v->point_size=n;
                v->point=(double*)malloc(sizeof(double)*n);
                for(int i=0;i<n;i++)
	            {   
                    fscanf(infile,"%lf",&(v->point[i]));
	            }
	   
	            rows++;
	        }
	   	    //printf("num c%d",m->num_c(m));
		    rows=0;

	        int id;

	        while(!feof(infile)&&rows<num[1])
	        {
                fscanf(infile,"%d",&cols);
 
                template_v**temp_v=(template_v**)malloc(sizeof(template_v*)*cols);
                for(int i=0;i<cols;i++)
	            {
                    fscanf(infile,"%d",&id);
		            auto it=m->vertices.find(id);
		            template_v* v;
		            if(it!=m->vertices.end())
		            {
		                v=it->second;
		            }
		            else
                    {
                        printf("can t find this vertex :%d\n",id);
		                return;
		            }
		            temp_v[i]=v;
	            }

	            m->create_cellv(m,temp_v,cols);

                free(temp_v);
                rows++; 
	        }


	   //printf("%d\r\n",m->num_f(m));
	    }
        else
        {
	   
	        printf("this file is not complete\r\n"); 

	        return;
	   }
    }
    else
	{
	    printf("this is not off\r\n");
	    return;
	
	}
}
