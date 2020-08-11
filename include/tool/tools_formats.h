#include<stdio.h>
#include<stdlib.h>
#include<string.h>
static double myfscanf_num(FILE * infile,char* b)
{
    double re=0,temp_d=0;
    unsigned int fac=1;
    int dot=0;
    char a;
    fscanf(infile,"%c",&a);
    *b=a;
    
    while(a==' ')
    {
        fscanf(infile,"%c",&a);
        *b=a;
    }
    while((a>='0'&&a<='9')||a=='.')
    {
        if(a=='.')
        {
            if(dot==0)
            {
                dot=1;
            }
            else
            {
                printf("cuowu\n");
                return 0;
            }
        }
        else
        {
            temp_d=a-'0';
            if(dot==0)
            {
                re*=10;
                re+=temp_d; 

            }
            else
            {
                fac*=10;
                temp_d/=fac;
                re+=temp_d;

            }
        }

        fscanf(infile,"%c",&a);

        *b=a;
    }    
    return re;
}
static void myfgets(FILE* infile)
{
    char a;
    do{
            
        fscanf(infile,"%c",&a);
    }while(a!='\n'&&!feof(infile));
}
void objtooff(char const * filename)
{
    int filenamesize=strlen(filename);
    char *outfilename=(char*)malloc(sizeof(char)*(filenamesize+2));
    strcpy(outfilename,filename);
    strcpy(&(outfilename[filenamesize-3]),"off");
    printf("%s\n",outfilename);
#ifdef WIN32
    FILE *infile,*outfile;
    fopen_s(&infile,filename,"r");
    fopen_s(&outfile,outfilename,"w");
#else
    FILE * infile=fopen(filename,"r");
    FILE *outfile=fopen(outfilename,"w");
#endif
    if(!infile)
    {
        printf("can't open this file %s\n",filename);
        return;
    }
    if(!outfile)
    {
        printf("can't open this file %s\n",filename);
        return;
    }
    char a;
    int v_rows=0,f_rows=0;
    while(!feof(infile))
    {
         fscanf(infile,"%c",&a);
        if(a=='\r')
        {
             fscanf(infile,"%c",&a);
        }
        if(a=='#')
        {
            myfgets(infile); 
        }
        else if(a=='v')
        {
            fscanf(infile,"%c",&a);
            if(a==' ')
            {
               v_rows++;
            }
            
            myfgets(infile);
            
        }
        else if(a=='f')
        {
            f_rows++;
            myfgets(infile);
        }
        else if(a=='\n')
        {
        }
        else
        {
            
           myfgets(infile);
        } 
    }
    fseek(outfile,0,SEEK_SET);
    fprintf(outfile,"OFF\n");
    fprintf(outfile,"%u %u %d\n",v_rows,f_rows,0);

    fseek(infile,0,SEEK_SET);
    unsigned int f_ids[10];
    while(!feof(infile))
    {
        double v_p;
        unsigned int f_id;
        fscanf(infile,"%c",&a);
        if(a=='\r')
        {
             fscanf(infile,"%c",&a);
        }
        if(a=='#')
        {
            myfgets(infile); 
        }
        else if(a=='v')
        {
            fscanf(infile,"%c",&a);
            if(a==' ')
            {
                for(int i=0;i<3;i++)
                {
                
                    fscanf(infile,"%lf ",&v_p);
                    //printf("%lf ",v_p);
                    fprintf(outfile,"%lf ",v_p);
                }
              //  v_rows++;
                fprintf(outfile,"\n");
            }
            else
            {
            
                myfgets(infile);
            }
            
        }
        else if(a=='f')
        {

            int cols=0;
            do{
                f_id=myfscanf_num(infile,&a);
                if(a=='\n')
                {
                    break;
                }
                f_ids[cols]=f_id;
                //fprintf(outfile,"%u ",f_id);
                cols++;
                //printf("%u\n",f_id); 
                while(a=='/')
                {
                    f_id=myfscanf_num(infile,&a); 
                }

            }while(a!='\n');
            //f_rows++;
            fprintf(outfile,"%u ",cols);

            for(int i=0;i<cols;i++)
            {
                fprintf(outfile,"%u ",f_ids[i]-1);
            
            }
            fprintf(outfile,"\n");
        }
        else if(a=='\n')
        {
        }
        else
        {
            
           myfgets(infile);
        }
    
    }
    
    printf("v_rows: %u,f_rows: %u here\n",v_rows,f_rows);
   free(outfilename); 
}

void offtocell(char const* filename)
{
    int filenamesize=strlen(filename);
    char *outfilename=(char*)malloc(filenamesize+2);
    for(int i=0;i<filenamesize;i++)
    {
        outfilename[i]=filename[i];
    
    }
    outfilename[filenamesize-3]='c';
    outfilename[filenamesize-2]='e';
    outfilename[filenamesize-1]='l';
    outfilename[filenamesize]='l';
    outfilename[filenamesize+1]='\0';
    printf("%s\n",outfilename);
#ifdef WIN32
    FILE* infile,*outfile;
    fopen_s(&infile,filename,"r");
    fopen_s(&outfile,outfilename,'w');
#else
    FILE* infile=fopen(filename,"r");
    FILE* outfile =fopen(outfilename,"w");
#endif
    if(!infile)
    {
        printf("can't open this file %s\n",filename);
        return ;
    
    }
    if(!outfile)
    {
    
        printf("can't open this file %s\n",outfilename);
        return;
    
    }
	free(outfilename);
    char off[6];
    fscanf(infile,"%s",off);
    int num_el[3],background_dim=0;int simplex=1;
    if(strcmp("OFF",off)==0)
    {
        //fwrite("CELL\n",1,5,ourfile);
        fscanf(infile,"%d %d %d",&num_el[0],&num_el[2],&num_el[1]);
        char a;int flag=0;
        for(fscanf(infile,"%c",&a);a!='\n'||background_dim==0;fscanf(infile,"%c",&a))
        {
            if(a==' ')
            {
                if(flag==1)
                {
                    background_dim++;
    
                }
                flag=0;

            }
            else
            {
                flag=1;

            }
            
        printf("a :%c\n",a);
        }
        if(flag==1){
        background_dim++;
        }
       // printf("background_dim:%d\n",background_dim);
        //fscanf(infile,"%c",&a);
        double b=0;
        for(int i=1;i<num_el[0];i++)
        {
            for(int j=0;j<background_dim;j++)
            {
        
                fscanf(infile,"%lf ",&b);
                //printf("%lf\n",b);
        
            }
        
        
        }
        int max_min[2];
        //这里的num_face是半面的个数
        int num_face=0; 
        int face_index[(background_dim+2)*num_el[2]];
        int face_temp=0;

        for(int i=0;i<num_el[2];i++)
        {   
            fscanf(infile,"%d",&max_min[0]);
            num_face+=max_min[0];
            face_index[face_temp]=max_min[0];
            face_temp++;

            if(max_min[0]>3)
            {
                simplex=0;
               // break;
            
            }
            for(int j=0;j<max_min[0];j++)
            {
            
            //    fscanf(infile,"%d",&max_min[1]);
                fscanf(infile,"%d",&face_index[face_temp]);
                face_temp++;

            }
            
        
        }
        //以上分析获取数据，以下写数据
        //printf("simplex: %d\n",simplex);
        fprintf(outfile,"CELL\n");
        fprintf(outfile,"background_dim= %d\n",background_dim);
        fprintf(outfile,"manifold_dim= %d\n",2);
        fprintf(outfile,"simplex= %d\n",simplex);
        if(simplex==1)
        {
            fprintf(outfile,"%d %d %d\n",num_el[0],0,num_el[2]); 
        
        }
        else
        {
        
        
            fprintf(outfile,"%d %d %d\n",num_el[0],num_face,num_el[2]);
        }
       
        fseek(infile,0,SEEK_SET);
        fscanf(infile,"%s %d %d %d",off,&max_min[1],&max_min[1],&max_min[1]);
         
        for(int i=0;i<num_el[0];i++)
        {
            for(int j=0;j<background_dim;j++)
            {
                fscanf(infile,"%lf",&b);
                fprintf(outfile,"%lf ",b);
            
            }
            fprintf(outfile,"\n"); 
        }
        /*for(int i=0;i<num_el[2];i++)
        {   
            fscanf(infile,"%d",&max_min[0]);
            face_index[face_temp]=max_min[0];
            face_temp++;
            for(int j=0;j<max_min[0];j++)
            {
            
                fscanf(infile,"%d",&face_index[face_temp]);
                face_temp++;
            
            }
            
        
        }*/

        if(simplex==0)
        {
            int i=0;
            while(i<face_temp)
            {
                int tmp_num=face_index[i];
                for(int j=0;j<tmp_num;j++)
                {
                    fprintf(outfile,"%d %d %d\n",2,face_index[i+(j)%tmp_num+1],face_index[i+(j+1)%tmp_num+1]);  
                    
               
                }
                i=i+tmp_num+1;
            }
            i=0;int tmp_id=0;
            while(i<face_temp)
            {
                int tmp_num=face_index[i];
                fprintf(outfile,"%d ",tmp_num);
                for(int j=0;j<tmp_num;j++)
                {
                    fprintf(outfile,"%d ",tmp_id);
                    tmp_id++;
                     
                }
                fprintf(outfile,"\n");
                i=i+tmp_num+1;
            }
           
        }
        else
        {
         
            int i=0;
            while(i<face_temp)
            {
                int tmp_num=face_index[i];
                fprintf(outfile,"%d ",tmp_num);
                for(int j=0;j<tmp_num;j++)
                {
                    fprintf(outfile,"%d ",face_index[i+j+1]);  
                    
               
                }
                fprintf(outfile,"\n");
                i=i+tmp_num+1;
            }
 

        }

       /* 
               for(int i=0;i<num_el[1];i++)
        {
            for(int j=0;j<)
            {
            
            
            }
        }
        */
        /*fwrite("CELL\n",1,5,outfile);
        char temp_s[100];
        sprintf(temp_s,"background_dim=%d\n",background_dim);
        fwrite(temp_s,1,5,outfile);*/
    }
}
void meshtocell(char const * filename)
{
    int filenamesize=strlen(filename);
    char *outfilename=(char*)malloc(filenamesize+1);
    for(int i=0;i<filenamesize;i++)
    {
        outfilename[i]=filename[i];
    
    }
    outfilename[filenamesize-4]='c';
    outfilename[filenamesize-3]='e';
    outfilename[filenamesize-2]='l';
    outfilename[filenamesize-1]='l';
    outfilename[filenamesize]='\0';

#ifdef WIN32
    FILE* infile,*outfile;
    fopen_s(&infile,filename,"r");
    fopen_s(&outfile,filename,"w");
#else
    FILE *outfile=fopen(outfilename,"w");
    FILE *infile=fopen(filename,"r");
#endif
    if(!infile||!outfile)
    {
        printf("cant open this file\r\n");
        return ;

    
    }
    fseek(infile,0,SEEK_SET);
    char mesh[20];
    fscanf(infile,"%s",mesh);
    double *vertices;
    int *tetrahedra;

//    printf("daozheli\n");
    int flag=0,dimension=3; 
    int num[3];
    while(!feof(infile))
    {
        if(strcmp("Dimension",mesh)==0)
        {
            fscanf(infile,"%d",&dimension);
            
        }
        else if(strcmp("#",mesh)==0)
        {    
            while(mesh[0]!='\n')
            {
                fscanf(infile,"%c",mesh);
                
            }
            
        }
        else if(strcmp("Vertices",mesh)==0)
        {
            fscanf(infile,"%d",&num[0]);
            printf("num0: %d\n",num[0]);
            double temp_d;
            vertices=(double*)malloc(sizeof(double)*dimension*num[0]);

            for(int i=0;i<num[0];i++)
            {
                for(int j=0;j<dimension;j++)
                {
                    fscanf(infile,"%lf ",&vertices[i*dimension+j]);

                
                }
                fscanf(infile,"%lf ",&temp_d);
            
            
            }
            flag++;
            
        }
        else if(strcmp("Tetrahedra",mesh)==0)
        {
                
            fscanf(infile,"%d",&num[1]);
            printf("num1: %d\n",num[1]);
            tetrahedra=(int*)malloc(sizeof(int)*4*num[1]);
            int temp_i,temp_i1;
            for(int i=0;i<num[1];i++)
            {
                for(int j=0;j<4;j++)
                {
                    fscanf(infile, "%d ",&temp_i1);
                     
		            tetrahedra[i*4+j]=temp_i1-1;
                }
                fscanf(infile,"%d ",&temp_i);
            }
            flag++;
            
        }
        if(flag==2)
        {    
            break;
        }
            
        
        fscanf(infile,"%s",mesh);
    }
    fprintf(outfile,"CELL\n");
    fprintf(outfile,"background_dim= %d\n",dimension);
    fprintf(outfile,"manifold_dim= %d\n",3);
    fprintf(outfile,"simplex= %d\n",1);
    fprintf(outfile,"%d %d %d\n",num[0],0,num[1]);
    for(int i=0;i<num[0];i++)
    {
        for(int j=0;j<dimension;j++)
        {
        
            fprintf(outfile,"%lf ",vertices[i*dimension+j]);
        }
        fprintf(outfile,"\n");

    
    }
    for(int i=0;i<num[1];i++)
    {
        fprintf(outfile,"%d ",4);
        for(int j=0;j<4;j++)
        {
        
            fprintf(outfile,"%d ",tetrahedra[i*4+j]);
        }
        fprintf(outfile,"\n");

    
    }
     
        /*for(int i=0;i<4;i++)
        {
            fscanf(infile,"%s",mesh);
        }

        fscanf(infile,"%d",&num[0]);
        printf("vertices,%d",num[0]);
        V.resize(num[0],3);

        V.setZero();
        for(int i=0;i<num[0];i++)
        {
            for(int j=0;j<3;j++)
            {
                fscanf(infile,"%lf",&vertices[i][j]);
                V.coeffRef(i,j)=vertices[i][j];
            
            }
            fscanf(infile,"%s",mesh); 
        }
        
        fscanf(infile,"%s",mesh);
        fscanf(infile,"%d",&num[1]);
        printf("num1%d\n",num[1]);
        int temp_int;
        F.resize(num[1],3);
        F.setZero();
        int faces[num[1]][3];
        for(int i=0;i<num[1];i++)
        {
            for(int j=0;j<3;j++)
            {
                 fscanf(infile,"%d",&temp_int);
                 F.coeffRef(i,j)=temp_int;
                 faces[i][j]=temp_int;
            }
            fscanf(infile,"%d",&temp_int);
            
        
        }
        fprintf(outfile,"%d %d %d\n",num[0],num[1],0);
        for(int i=0;i<num[0];i++)
        {
            for(int j=0;j<3;j++)
            {
            
                fprintf(outfile,"%lf ",vertices[i][j]);

            }
            fprintf(outfile,"\n");
        }
        for(int i=0;i<num[1];i++)
        {
        
            for(int j=0;j<3;j++)
            {
                fprintf(outfile,"%d ",faces[i][j]);
            }
            fprintf(outfile,"\n");
        }

        printf("vertices:%d num[1]:%d\n",num[0],num[1]);

        */
    
    free(vertices);
    free(tetrahedra);
}

