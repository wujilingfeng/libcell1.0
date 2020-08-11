#include<Mesh/lib_cell_Iterator.h>
void iterator_v_init(iterator_v* it)
{
    it->i=0;
    it->value=NULL;
}
iterator_v operator++(iterator_v& vv_)
{

    (vv_.i++);
    return vv_;
}
iterator_v operator++(iterator_v&vv_,int)
{
    iterator_v vv=vv_;
    vv_.i++;
    return vv;


}
template_v* quote(iterator_v &vv_)
{
    return vv_.value[vv_.i];
}
template_v operator*(iterator_v &vv_)
{
    return *quote(vv_);

}
bool operator!=(const iterator_v&v_1,const iterator_v& v_2)
{
    if(v_1.value!=v_2.value)
    {
        return true;
    }
    if(v_1.i!=v_2.i)
    {
        return true;
    }
    return false;

}
void iterator_f_init(iterator_f*it)
{
    Node_init(&(it->node));
}
iterator_f operator++(iterator_f& f_)
{
    f_.node--;
    return f_;
}
iterator_f operator++(iterator_f&f_,int)
{
    iterator_f f1=f_;
    f_.node--;
    return f1;

}
template_f* quote(iterator_f&f_)
{

    return ((template_f*)(*(f_.node)));

}
template_f operator*(iterator_f&f_)
{

    return *quote(f_);
}
bool operator!=(const iterator_f& f_1,const iterator_f& f_2)
{
    if(f_1.node.value!=f_2.node.value)
    {
        return true;
    }
    return false;
}
void iterator_c_init(iterator_c*it)
{
    Node_init(&(it->node));
}
iterator_c operator++(iterator_c& c_)
{
    c_.node++;
    return c_;
}
iterator_c operator++(iterator_c &c_,int)
{
    iterator_c c_1=c_;
    c_.node++;
    return c_1;

}

template_c operator*(iterator_c&c_)
{

    return *quote(c_);
}
template_c* quote(iterator_c&c_)
{
    return (template_c*)(*(c_.node));

}
bool operator!=(const iterator_c& c_1,const iterator_c& c_2)
{
    if(c_1.node.value!=c_2.node.value)
    {
        return true;
    }
    return false;
}
void iterator_hf_init(iterator_hf* hf_)
{
    Node_init(&(hf_->node));

}
iterator_hf operator++(iterator_hf& it)
{
    it.node--;
    return it;

}
iterator_hf operator++(iterator_hf& it,int)
{
    iterator_hf hf_=it;
    it.node--;
    return hf_;

}
template_hf* quote(iterator_hf&it)
{
    return (template_hf*)(*(it.node));

}
template_hf operator*(iterator_hf& it)
{
    return *quote(it);

}
bool operator!=(const iterator_hf& it1,const iterator_hf&it2)
{
    if(it1.node.value!=it2.node.value)
    {
        return true;
    }
    return false;
}
