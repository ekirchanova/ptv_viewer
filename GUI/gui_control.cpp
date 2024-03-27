#include "GUI.h"
#include "IO/IO.h"
#include <stdlib.h>
#include <stdio.h>

float rx=0.0;
float ry=0.0;
int mx0,my0;
int rotate=1;
float rx0=0;
float ry0=0;

float t_x=0.0;
float t_y=0.0;
float t_z=-2.4;

float sc =1.0;
int mode =0;
int screen_num =0;
void kb(unsigned char key, int x, int y){ //keyboard event
    if (key==' '){
     printf("key pressed \n");
    }

    if (key=='P'){
        mode=0;
     printf("perspective projection is set \n");
     glMatrixMode (GL_PROJECTION);
     glLoadIdentity ();
     gluPerspective(60, W_HEIGHT*1.0/W_WIDTH, 0.0001, 150);
     glMatrixMode (GL_MODELVIEW);

    }

    if (key=='O'){
        mode=0;
     printf("orthogonal projection is set \n");
     glMatrixMode (GL_PROJECTION);
     glLoadIdentity ();
     glOrtho(w_x0, w_x1, W_HEIGHT*w_y0/W_WIDTH,W_HEIGHT*w_y1/W_WIDTH, -10.0, 10.0);
     glMatrixMode (GL_MODELVIEW);
    }

    if(key=='q'){
        t_z-=0.01;
        printf("t_z = %f \n",t_z);
    }

    if(key=='e'){
        t_z+=0.01;
        printf("t_z = %f \n",t_z);
    }
    if(key=='w'){
        t_y+=0.01;
        printf("t_y = %f \n",t_y);
    }
    if(key=='s'){
        t_y-=0.01;
        printf("t_y = %f \n",t_y);
    }
    if(key=='d'){
        t_x+=0.01;
        printf("t_x = %f \n",t_x);
    }
    if(key=='a'){
        t_x-=0.01;
        printf("t_x = %f \n",t_x);
    }


    if (key=='S'){
        char fname[4096];
         sprintf(fname,"%s/test_out/%d.bmp",dir,screen_num);
     printf("saving screenshot to %s \n",fname);
     screenshot_bmp(fname, W_WIDTH, W_HEIGHT);
     screen_num++;
    }

    if (key>='1' && key <='7'){


      active_camera =key-'1';
     printf("%d cam view \n",active_camera);
      mode =1;

      glMatrixMode (GL_PROJECTION);
      glLoadIdentity ();
      glLoadMatrixd(cs.projM.m);
      glMatrixMode (GL_MODELVIEW);
      glLoadMatrixd(cs.modelvMatrices[active_camera].m);
    }





    if (key=='C'){
      GLfloat m_matrix[16];
      glGetFloatv (GL_MODELVIEW_MATRIX, m_matrix);

      GLfloat p_matrix[16];
      glGetFloatv (GL_PROJECTION_MATRIX, p_matrix);

      printf("CURR_PROJ: \n%f %f %f %f \n ",p_matrix[0],p_matrix[4],p_matrix[8],p_matrix[12]);
      printf("%f %f %f %f \n ",p_matrix[1],p_matrix[5],p_matrix[9],p_matrix[13]);
      printf("%f %f %f %f \n ",p_matrix[2],p_matrix[6],p_matrix[10],p_matrix[14]);
      printf("%f %f %f %f \n ",p_matrix[3],p_matrix[7],p_matrix[11],p_matrix[15]);

      GLdouble *projm = cs.projM.m;
      GLdouble *modelvm = cs.modelvMatrices[active_camera].m;
      printf("OPENCV_PROJ: \n%f %f %f %f \n ",projm[0],projm[4],projm[8],projm[12]);
      printf("%f %f %f %f \n ",projm[1],projm[5],projm[9],projm[13]);
      printf("%f %f %f %f \n ",projm[2],projm[6],projm[10],projm[14]);
      printf("%f %f %f %f \n ",projm[3],projm[7],projm[11],projm[15]);


      printf("\n CURR_MODELVIEW: \n%f %f %f %f \n ",m_matrix[0],m_matrix[4],m_matrix[8],m_matrix[12]);
      printf("%f %f %f %f \n ",m_matrix[1],m_matrix[5],m_matrix[9],m_matrix[13]);
      printf("%f %f %f %f \n ",m_matrix[2],m_matrix[6],m_matrix[10],m_matrix[14]);
      printf("%f %f %f %f \n ",m_matrix[3],m_matrix[7],m_matrix[11],m_matrix[15]);

      printf("\n OPENCV_MODELVIEW: \n%f %f %f %f \n ",modelvm[0],modelvm[4],modelvm[8],modelvm[12]);
      printf("%f %f %f %f \n ",modelvm[1],modelvm[5],modelvm[9], modelvm[13]);
      printf("%f %f %f %f \n ",modelvm[2],modelvm[6],modelvm[10],modelvm[14]);
      printf("%f %f %f %f \n ",modelvm[3],modelvm[7],modelvm[11],modelvm[15]);

    }

    glutPostRedisplay();
}

void skb(int key, int x, int y){ //special keys event
    if (key == GLUT_KEY_UP){
        sc*=1.1;
       // printf("special key pressed \n");
    }
    if (key == GLUT_KEY_DOWN){
        sc/=1.1;
       // printf("special key pressed \n");
    }

    if (key == GLUT_KEY_RIGHT){
        active_camera++;
        if (active_camera>=7) active_camera =6;
       // printf("special key pressed \n");
    }

    if (key == GLUT_KEY_LEFT){
        active_camera--;
        if (active_camera<0) active_camera =0;
       // printf("special key pressed \n");
    }
    glutPostRedisplay();
}

void m_m(int x,int y){ //mouse move
    if (rotate==1){
        rx=rx0+0.5*(x-mx0);
        ry=ry0+0.5*(y-my0);
    }
    glutPostRedisplay();
}

void m_d(int button, int state,int x, int y){  //mouse down
    if (state==GLUT_UP)
    {
        rotate=0;
        rx0=rx;
        ry0=ry;
    }
    if (state==GLUT_DOWN)
    {
        rotate=1;
        mx0=x;
        my0=y;

    }    
    glutPostRedisplay();
}
