#include "GUI.h"

get_color(float gval, color_table* tbl,float min, float max){
    int i;
    float val;
    val=gval;
    if (val>max) val=max-0.00001;
    if (val<min) val=min+0.00001;

    double alpha;

    if ((max-min)>0.00001){
        alpha=(val-min)/(max-min)*(tbl->x_num);
        i=tbl->xt[(int)(alpha)];
        // alpha=alpha-i;
    }else{
        alpha=0.0;
        i=2;
    }
    alpha=(alpha-tbl->x[i])*1.0/(tbl->x[i+1]-tbl->x[i]);

    VEC3 res;
    res.x=tbl->r[i]*(1.0-alpha)+alpha*tbl->r[i+1];
    res.y=tbl->g[i]*(1.0-alpha)+alpha*tbl->g[i+1];
    res.z=tbl->b[i]*(1.0-alpha)+alpha*tbl->b[i+1];

    glColor3f(res.x,res.y,res.z);
    //glColor3f(alpha,alpha,-alpha);
}
