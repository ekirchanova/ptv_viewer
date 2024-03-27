#include "IO.h"

/*
void save_tecplot(char* filename,int i_max, int j_max,int k_max,double dr)
{
    VEC3 pos;
    VEC3 v;
    FILE*    f=fopen(filename, "w");
    fprintf(f, "TITLE =\"  \" \n VARIABLES=\"X\" \n \"Y\" \n \"Z\" \n \"U\"\n \"V\"\n \"W\"\n zone T=\"  \" \n");
    fprintf(f," I=%d J=%d K=%d F=POINT \n", i_max, j_max, k_max);
    int i,j,k;
    for ( k=-k_max/2;k<k_max/2;k++){
        printf ("saving ,k=%d \n",k);
        for ( j=-j_max/2; j<j_max/2; j++){
            for ( i=-i_max/2; i<i_max/2; i++){
                pos.x = i*dr;
                pos.y = j*dr;
                pos.z = k*dr;
                v = get_v(pos);
                fprintf(f, "%e  %e  %e %e %e %e \n",
                        pos.x,pos.y,pos.z,
                        v.x,v.y,v.z);
            }

        }
    }
    fclose(f);
}
*/
