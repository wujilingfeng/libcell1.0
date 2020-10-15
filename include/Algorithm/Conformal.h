#ifndef LIB_CELL_CONFORMAL_H
#define LIB_CELL_CONFORMAL_H
#include "Harmonic.h"

class Conformal_module
{
public:
    Conformal_module()
    {

    }
    ~Conformal_module()
    {

    }
    void init_vs()
    {
        Mesh* mesh=ham->mesh;
        mesh->external_cell_init_(mesh);

        int size=node_size(mesh->external_cell.halffaces);
        int step=size/4;
        
        int i=1,k=0;

    //  printf("%d\n",step);
        int steps[4];
        steps[0]=0;steps[1]=step;steps[2]=2*step;steps[3]=3*step;
        for(Node* nit=mesh->external_cell.halffaces;nit!=NULL;nit=(Node*)(nit->Next))
        {
            if(i==steps[0]||i==steps[1]||i==steps[2]||i==steps[3])
            {
                v_keys[k]=((template_hf*)(nit->value))->vertices[1];
                k++;
            }
            i++;
        }   
    }
    void neumman_boundary()
    {
        Eigen::MatrixXd temp_H,temp_M;
        Mesh* mesh=ham->mesh;
        mesh->external_cell_init_(mesh);
        int rows=node_size(mesh->external_cell.halffaces);  
        temp_H.resize(rows,mesh->num_v(mesh));
        temp_H.setZero();
        temp_M.resize(mesh->num_v(mesh),rows);
        temp_M.setZero();

        Eigen::VectorXd vx(rows);
        vx.setZero();
        int i=0,j=0;
        Ham_v_prop * temp=NULL;
        template_v* v=NULL;
        template_v* vs1[2],*vs2[2];
        vs1[0]=v_keys[0];vs1[1]=v_keys[1];
        vs2[0]=v_keys[2];vs2[1]=v_keys[3];
        int value_=0;
        bool flag=false;
        for(Node* nit=mesh->external_cell.halffaces;nit!=NULL;nit=(Node*)(nit->Next))
        {
            template_v* v1=(((template_hf*)(nit->value))->vertices[1]);
            temp=(Ham_v_prop*)(v1->user_prop);
            i=temp->index;
            if(v1==vs1[0]||v1==vs2[0])
            {
                flag=true;
            }
            if(v1==vs1[1]||v1==vs2[1])
            {
                flag=false;
                value_++;
            }
            if(flag)
            {
                temp_H.coeffRef(i,i)=value_;
                vx.coeffRef(i)=0;
            }   
            else
            {
                Node* vvit=mesh->vv_begin(mesh,*(v1));
                Node* temp_vvit=vvit;
                for(;vvit!=mesh->vv_end(mesh,*(v1));vvit=(Node*)(vvit->Next))
                {
                    v=(template_v*)(vvit->value);
                    temp=(Ham_v_prop*)(v->user_prop);
                    j=temp->index;
                    if(!mesh->vertex_is_boundary(mesh,*v))
                    {
                        j+=rows;
                    }
                    template_v* vs[2];
                    vs[0]=v1;vs[1]=v;
                    template_f* f=mesh->get_facev(mesh,vs,2);
                    Ham_f_prop* temp1=(Ham_f_prop*)(f->user_prop);
                    temp_H.coeffRef(i,j)=-temp1->weight;
                    temp_H.coeffRef(i,i)+=temp1->weight;    
                }       
                free_node(temp_vvit);           
            }           
        }
        
        temp_M.block(0,0,rows,rows)=Eigen::MatrixXd::Identity(rows, rows);
        temp_M.block(rows,0,mesh->num_v(mesh)-rows,rows)=ham->LLE;
        Eigen::MatrixXd temp_MM=temp_H*temp_M;

        Eigen::SparseMatrix<double> MM;
        MM.resize(rows,rows);
        MM.setZero();
        for(i=0;i<rows;i++)
        {
            for(j=0;j<rows;j++)
            {
                MM.coeffRef(i,j)=temp_MM.coeffRef(i,j);
            }
        }
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>llt;
        llt.compute(MM);
        Eigen::VectorXd x=llt.solve(vx);


    }
template_v* v_keys[4];
Harmonic * ham;
};
#endif