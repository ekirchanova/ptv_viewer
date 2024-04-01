
#include <stdlib.h>
#include  <math.h>
#include <stdio.h>
#include <time.h>

#include <GUI/GUI.h>


char dir[1024];
int main(int argc, char** argv){
   //  calc_radials();

    char fname[1024];
    sprintf(dir,"%s",__FILE__);
    dir[strlen(dir)-9]=0; //hack to get source directory



    GUI_init(argc,argv);
}
