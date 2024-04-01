#include "GUI.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

void gL4x4MatrixInverse(OpenGLMatrix *M, OpenGLMatrix *Inv){
    double  det;
    int i;

    double *m = M->m;
    double *inv =Inv->m;

    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

 //   if (det == 0)
   //     return false;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        inv[i] = inv[i] * det;

  //  return true;
}


void gL4x4MatrixMult(OpenGLMatrix *Ma, OpenGLMatrix *Mb, OpenGLMatrix *Mres){
    for (int i=0;i<16;i++){
        int row = i%4;
        int col = i/4;
        Mres->m[i]=0.0;
        for (int j=0; j<4;j++){
            Mres->m[i] += Ma->m[row+4*j]*Mb->m[col*4 + j];
         //   printf("m %d %d %f %f \n",row+4*j,col*4 + j,Ma.m[row+4*j],Mb.m[col*4 + j]);
        }

    }
 //   for (int i=0;i<16;i++)
   //      printf("m %d %f \n",i,Mres.m[i]);
}

void printMatr(OpenGLMatrix *Ma){
   // for (int i=0;i<16;i++)
     //    printf("m %d %f \n",i,Ma.m[i]);
 GLdouble * p_matrix = Ma->m;
    printf("%f %f %f %f \n ",p_matrix[0],p_matrix[4],p_matrix[8],p_matrix[12]);
    printf("%f %f %f %f \n ",p_matrix[1],p_matrix[5],p_matrix[9],p_matrix[13]);
    printf("%f %f %f %f \n ",p_matrix[2],p_matrix[6],p_matrix[10],p_matrix[14]);
    printf("%f %f %f %f \n ",p_matrix[3],p_matrix[7],p_matrix[11],p_matrix[15]);
}

void drawCameraFrustum(camera_sequence *c, int active_image){


    OpenGLMatrix M,InvM, zNear, zFar,zNear_w,zFar_w;
    gL4x4MatrixMult(&(c->projM),&(c->modelvMatrices[active_image]),&M);
   // printf("\n pXm \n");
   // printMatr(&M);
    gL4x4MatrixInverse(&M,&InvM);



    zNear.m[0] = -1.0; zNear.m[4] =  1.0; zNear.m[8]  =  1.0;  zNear.m[12] = -1.0;  //x
    zNear.m[1] = -1.0; zNear.m[5] = -1.0; zNear.m[9]  =  1.0;  zNear.m[13] =  1.0;  //y
    zNear.m[2] = -1.0; zNear.m[6] = -1.0; zNear.m[10] = -1.0;  zNear.m[14] = -1.0;  //z
    zNear.m[3] =  1.0; zNear.m[7] =  1.0; zNear.m[11] =  1.0;  zNear.m[15] =  1.0;  //w

    zFar.m[0] = -1.0; zFar.m[4] =  1.0; zFar.m[8]  =  1.0;  zFar.m[12] = -1.0;  //x
    zFar.m[1] = -1.0; zFar.m[5] = -1.0; zFar.m[9]  =  1.0;  zFar.m[13] =  1.0;  //y
    zFar.m[2] =  1.0; zFar.m[6] =  1.0; zFar.m[10] =  1.0;  zFar.m[14] =  1.0;  //z
    zFar.m[7] =  1.0; zFar.m[3] =  1.0; zFar.m[11] =  1.0;  zFar.m[15] =  1.0;  //w

    gL4x4MatrixMult(&InvM,&zNear,&zNear_w);
    gL4x4MatrixMult(&InvM,&zFar,&zFar_w);

    printf("projection: \n");
    printMatr(&(c->projM));
    printf("\n modelview: \n");
    printMatr(&(c->modelvMatrices[active_image]));
    printf("\n pXm \n");
    printMatr(&M);
    printf("\n invM \n");
    printMatr(&InvM);
    printf("\n zNearw \n");
    printMatr(&zNear_w);

    printf("\n zFarw \n");
    printMatr(&zFar_w);


    //dividing by weight
    for (int i=0;i<4;i++){
        for (int j=0;j<3;j++){
          zNear_w.m[i*4+j] = zNear_w.m[i*4+j]/zNear_w.m[i*4+3];
          zFar_w.m[i*4+j] = zFar_w.m[i*4+j]/zFar_w.m[i*4+3];
        }
        zNear_w.m[i*4+3] = 1.0;
        zFar_w.m[i*4+3] = 1.0;
    }
    glLineWidth(2.0);
    glBegin(GL_LINE_LOOP);
    glColor3f(0.0,0.5,0.0);
    glVertex3f(zNear_w.m[0], zNear_w.m[1], zNear_w.m[2]);
    glVertex3f(zNear_w.m[4], zNear_w.m[5], zNear_w.m[6]);
    glVertex3f(zNear_w.m[8], zNear_w.m[9], zNear_w.m[10]);
    glVertex3f(zNear_w.m[12],zNear_w.m[13],zNear_w.m[14]);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glColor3f(0.0,0.5,0.0);
    glVertex3f(zFar_w.m[0], zFar_w.m[1], zFar_w.m[2]);
    glVertex3f(zFar_w.m[4], zFar_w.m[5], zFar_w.m[6]);
    glVertex3f(zFar_w.m[8], zFar_w.m[9], zFar_w.m[10]);
    glVertex3f(zFar_w.m[12],zFar_w.m[13],zFar_w.m[14]);
    glEnd();

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, c->texture_indices[active_image]);
    glColor3f(1.0,1.0,1.0);


    glBegin(GL_QUADS);
     glTexCoord2f(0.0f, 0.0f); glVertex3f(zFar_w.m[0], zFar_w.m[1], zFar_w.m[2]);    // Низ лево
     glTexCoord2f(1.0f, 0.0f); glVertex3f(zFar_w.m[4], zFar_w.m[5], zFar_w.m[6]);    // Низ право
     glTexCoord2f(1.0f, 1.0f); glVertex3f(zFar_w.m[8], zFar_w.m[9], zFar_w.m[10]);    // Верх право
     glTexCoord2f(0.0f, 1.0f); glVertex3f(zFar_w.m[12],zFar_w.m[13],zFar_w.m[14]);    // Верх лево
    glEnd();


    glDisable(GL_TEXTURE_2D);

    /*draw_points(c->imagePoints[active_image], 0.0, 1.0, 0.0, 0.0);
    printf("found %d points\n",c->imagePoints[active_image].size());
    draw_points(c-> undistortImagePoints[active_image], 0.0, 0.0, 1.0, 0.0);
    printf("found %d points\n",c->imagePoints[active_image].size());*/

    glBegin(GL_LINES);
    glColor3f(0.0,0.5,0.0);
    glVertex3f(zNear_w.m[0], zNear_w.m[1], zNear_w.m[2]);
    glVertex3f(zFar_w.m[0], zFar_w.m[1], zFar_w.m[2]);

    glVertex3f(zNear_w.m[4], zNear_w.m[5], zNear_w.m[6]);
    glVertex3f(zFar_w.m[4], zFar_w.m[5], zFar_w.m[6]);

    glVertex3f(zNear_w.m[8], zNear_w.m[9], zNear_w.m[10]);
    glVertex3f(zFar_w.m[8], zFar_w.m[9], zFar_w.m[10]);

    glVertex3f(zNear_w.m[12],zNear_w.m[13],zNear_w.m[14]);
    glVertex3f(zFar_w.m[12],zFar_w.m[13],zFar_w.m[14]);
    glEnd();


}
