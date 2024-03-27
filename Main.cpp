
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

    sprintf(fname,"%s/test_in/ImperxCamera0_18_10_123/test.txt",dir);

    printf("Program started, sorce dir is %s \n fname id %s \n", dir,fname);
    FILE* f = fopen(fname,"r");
    char rrr[1024];
    fscanf(f,"%s",rrr);
    printf("the file content is:\n %s \n",rrr);
    fclose(f);

    GUI_init(argc,argv);
}
