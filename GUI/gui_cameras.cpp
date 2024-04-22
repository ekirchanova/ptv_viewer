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

//findning nearest point between two lines, return sqaured distance between two lines
double findNearestPoint(double l1_p1[3], double l1_p2[3], double l2_p1[3], double l2_p2[3], double res[3]){
    double n1[3],n2[3]; //two normals of the first line

    double dl_1[3],dl_2[3],d_l2_l1[3];
    SUB(dl_1,l1_p2,l1_p1);
    SUB(dl_2,l2_p2,l2_p1);
    n1[0] = 1.0; n1[1] = 1.0; n1[2] = 1.0; //just a test vector;
    CROSS(n2,dl_1,n1);
    double zl = 1.0/sqrt(DOT(n2,n2));
    SMULT(n2,zl,n2);
    CROSS(n1,n2,dl_1);
    zl = 1.0/sqrt(DOT(n1,n1));
    SMULT(n1,zl,n1);
    //now we got two norlmals

    //the second line will be parametrized with l1_p1 + dl_2*t and we'll search for t nearest to line1
    SUB(d_l2_l1, l2_p1, l1_p1);
    double dot_l_n1,dot_l_n2, dot_dl2_n1,dot_dl2_n2;
    dot_l_n1 = DOT(d_l2_l1,n1);
    dot_l_n2 = DOT(d_l2_l1,n2);
    dot_dl2_n1 = DOT(dl_2,n1);
    dot_dl2_n2 = DOT(dl_2,n2);

    //now let's calculate t for minimal distance:
    double t = -(dot_l_n1*dot_dl2_n1 + dot_l_n2*dot_dl2_n2) / (dot_dl2_n1*dot_dl2_n1 + dot_dl2_n2*dot_dl2_n2);
    SMULT(res,t,dl_2);
    SUM(res, res, d_l2_l1); //relative vector to the l1_p1
    double dx,dy;
    dx= DOT(res,n1);
    dy= DOT(res,n2);
    SUM(res, res, l1_p1); //adding l1_p1 to get the absolute coordinate
    double l2_min =  dx*dx + dy*dy;
    return l2_min;
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
   /* glBegin(GL_LINE_LOOP);
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
    glEnd();*/

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

   /* glBegin(GL_LINES);
>>>>>>> 3cb46f5ac1ffecdffba37af6f8701c7761dc660d
    glColor3f(0.0,0.5,0.0);
    glVertex3f(zNear_w.m[0], zNear_w.m[1], zNear_w.m[2]);
    glVertex3f(zFar_w.m[0], zFar_w.m[1], zFar_w.m[2]);

    glVertex3f(zNear_w.m[4], zNear_w.m[5], zNear_w.m[6]);
    glVertex3f(zFar_w.m[4], zFar_w.m[5], zFar_w.m[6]);

    glVertex3f(zNear_w.m[8], zNear_w.m[9], zNear_w.m[10]);
    glVertex3f(zFar_w.m[8], zFar_w.m[9], zFar_w.m[10]);

    glVertex3f(zNear_w.m[12],zNear_w.m[13],zNear_w.m[14]);
    glVertex3f(zFar_w.m[12],zFar_w.m[13],zFar_w.m[14]);
    glEnd();*/

    glLineWidth(1.0);
    if (active_image==0)
        glColor3f(0,1,0);
    else
        glColor3f(1,1,1);

    glPointSize(4);
    glDisable(GL_DEPTH_TEST);
    glBegin(GL_LINES);
  /*  for (int i=0; i<c->imagePoints[active_image].size();i+=4){
=======
    for (int i=counter; i<c->imagePoints[active_image].size();i+=19){
>>>>>>> 3cb46f5ac1ffecdffba37af6f8701c7761dc660d
        double a,b,dx0[3],dy0[3],dx1[3],dy1[3],x0[3],x1[3];
        a=c->imagePoints[active_image][i].x*1.0/c->images[active_image].size().width;
        b=1.0-(c->imagePoints[active_image][i].y*1.0/c->images[active_image].size().height);
        dx1[0] = zFar_w.m[4] - zFar_w.m[0];
        dx1[1] = zFar_w.m[5] - zFar_w.m[1];
        dx1[2] = zFar_w.m[6] - zFar_w.m[2];

        dy1[0] = (zFar_w.m[8] - zFar_w.m[4]);
        dy1[1] = (zFar_w.m[9] - zFar_w.m[5]);
        dy1[2] = (zFar_w.m[10]- zFar_w.m[6]);

        dx0[0] = zNear_w.m[4] - zNear_w.m[0];
        dx0[1] = zNear_w.m[5] - zNear_w.m[1];
        dx0[2] = zNear_w.m[6] - zNear_w.m[2];

        dy0[0] = (zNear_w.m[8] - zNear_w.m[4]);
        dy0[1] = (zNear_w.m[9] - zNear_w.m[5]);
        dy0[2] = (zNear_w.m[10]- zNear_w.m[6]);

        x1[0] =  zFar_w.m[0];
        x1[1] =  zFar_w.m[1];
        x1[2] =  zFar_w.m[2];

        x0[0] =  zNear_w.m[0];
        x0[1] =  zNear_w.m[1];
        x0[2] =  zNear_w.m[2];

        double r0[3],dr[3],r1[3];
        for(int j=0;j<3;++j){
            r0[j] = x0[j] + dx0[j]*a + dy0[j]*b;
            r1[j] = x1[j] + dx1[j]*a + dy1[j]*b;
            dr[j] = x1[j] - x0[j] + (dx1[j]-dx0[j])*a + (dy1[j]-dy0[j])*b;
        }
        glVertex3f(r0[0],r0[1],r0[2]);
        double dz = fabs(r0[2]/dr[2]);
       // glVertex3f(r1[0],r1[1],r1[2]);
        glVertex3f(r0[0]+dr[0]*dz,r0[1]+dr[1]*dz,r0[2]+dr[2]*dz);
        // glVertex3f(x1[0]+ dx1[0]*a + dy1[0]*b, x1[1]+ dx1[1]*a + dy1[1]*b, x1[2]+ dx1[2]*a + dy1[2]*b);
    }

    glEnd();*/

    /*glBegin(GL_POINTS);
    glColor3f(0.0,0.0  ,0.0);
    for (int i=0; i<c->undistortPoints[active_image].size();i+=4){
        double a,b,dx0[3],dy0[3],dx1[3],dy1[3],x0[3],x1[3];
        a=c->undistortPoints[active_image][i].x*1.0/c->images[active_image].size().width;
        b=1.0-(c->undistortPoints[active_image][i].y*1.0/c->images[active_image].size().height);
        dx1[0] = zFar_w.m[4] - zFar_w.m[0];
        dx1[1] = zFar_w.m[5] - zFar_w.m[1];
        dx1[2] = zFar_w.m[6] - zFar_w.m[2];

        dy1[0] = (zFar_w.m[8] - zFar_w.m[4]);
        dy1[1] = (zFar_w.m[9] - zFar_w.m[5]);
        dy1[2] = (zFar_w.m[10]- zFar_w.m[6]);

        dx0[0] = zNear_w.m[4] - zNear_w.m[0];
        dx0[1] = zNear_w.m[5] - zNear_w.m[1];
        dx0[2] = zNear_w.m[6] - zNear_w.m[2];

        dy0[0] = (zNear_w.m[8] - zNear_w.m[4]);
        dy0[1] = (zNear_w.m[9] - zNear_w.m[5]);
        dy0[2] = (zNear_w.m[10]- zNear_w.m[6]);

        x1[0] =  zFar_w.m[0];
        x1[1] =  zFar_w.m[1];
        x1[2] =  zFar_w.m[2];

        x0[0] =  zNear_w.m[0];
        x0[1] =  zNear_w.m[1];
        x0[2] =  zNear_w.m[2];

        double r0[3],dr[3];
        for(int j=0;j<3;++j){
            r0[j] = x0[j] + dx0[j]*a + dy0[j]*b;
            dr[j] = x1[j] - x0[j] + (dx1[j]-dx0[j])*a + (dy1[j]-dy0[j])*b;
        }
        glVertex3f(r0[0],r0[1],r0[2]);
        double dz = fabs(r0[2]/dr[2]);
        glVertex3f(r0[0]+dr[0]*dz,r0[1]+dr[1]*dz,r0[2]+dr[2]*dz);
        // glVertex3f(x1[0]+ dx1[0]*a + dy1[0]*b, x1[1]+ dx1[1]*a + dy1[1]*b, x1[2]+ dx1[2]*a + dy1[2]*b);
    }

    glEnd();*/
    glEnable(GL_DEPTH_TEST);

}


void drawIntersections_min(camera_sequence *c){

    std::vector<VEC3> points;
    std::vector<VEC3> p0_0;
    std::vector<VEC3> p0_1;

    std::vector<VEC3> p1_0;
    std::vector<VEC3> p1_1;

    OpenGLMatrix M,InvM, zNear, zFar,zNear_w[2],zFar_w[2];

    for (int nn = 0; nn<2; ++nn){
        gL4x4MatrixMult(&(c->projM),&(c->modelvMatrices[nn]),&M);
        gL4x4MatrixInverse(&M,&InvM);

        zNear.m[0] = -1.0; zNear.m[4] =  1.0; zNear.m[8]  =  1.0;  zNear.m[12] = -1.0;  //x
        zNear.m[1] = -1.0; zNear.m[5] = -1.0; zNear.m[9]  =  1.0;  zNear.m[13] =  1.0;  //y
        zNear.m[2] = -1.0; zNear.m[6] = -1.0; zNear.m[10] = -1.0;  zNear.m[14] = -1.0;  //z
        zNear.m[3] =  1.0; zNear.m[7] =  1.0; zNear.m[11] =  1.0;  zNear.m[15] =  1.0;  //w

        zFar.m[0] = -1.0; zFar.m[4] =  1.0; zFar.m[8]  =  1.0;  zFar.m[12] = -1.0;  //x
        zFar.m[1] = -1.0; zFar.m[5] = -1.0; zFar.m[9]  =  1.0;  zFar.m[13] =  1.0;  //y
        zFar.m[2] =  1.0; zFar.m[6] =  1.0; zFar.m[10] =  1.0;  zFar.m[14] =  1.0;  //z
        zFar.m[7] =  1.0; zFar.m[3] =  1.0; zFar.m[11] =  1.0;  zFar.m[15] =  1.0;  //w

        gL4x4MatrixMult(&InvM,&zNear,&(zNear_w[nn]));
        gL4x4MatrixMult(&InvM,&zFar,&(zFar_w[nn]));

        //dividing by weight
        for (int i=0;i<4;i++){
            for (int j=0;j<3;j++){
                zNear_w[nn].m[i*4+j] = zNear_w[nn].m[i*4+j]/zNear_w[nn].m[i*4+3];
                zFar_w[nn].m[i*4+j] = zFar_w[nn].m[i*4+j]/zFar_w[nn].m[i*4+3];
            }
            zNear_w[nn].m[i*4+3] = 1.0;
            zFar_w[nn].m[i*4+3] = 1.0;
        }
    }
    glDisable(GL_DEPTH_TEST);


    VEC3 p_min_0_0, p_min_1_0,p_min_0_1,p_min_1_1;

    for (int nn=counter; nn<c->imagePoints[0].size();nn+=19){
        double a,b,dx0[3],dy0[3],dx1[3],dy1[3],x0[3],x1[3];
        a=c->imagePoints[0][nn].x*1.0/c->images[0].size().width;
        b=1.0-(c->imagePoints[0][nn].y*1.0/c->images[0].size().height);
        dx1[0] = zFar_w[0].m[4] - zFar_w[0].m[0];
        dx1[1] = zFar_w[0].m[5] - zFar_w[0].m[1];
        dx1[2] = zFar_w[0].m[6] - zFar_w[0].m[2];
        dy1[0] = zFar_w[0].m[8] - zFar_w[0].m[4];
        dy1[1] = zFar_w[0].m[9] - zFar_w[0].m[5];
        dy1[2] = zFar_w[0].m[10]- zFar_w[0].m[6];

        dx0[0] = zNear_w[0].m[4] - zNear_w[0].m[0];
        dx0[1] = zNear_w[0].m[5] - zNear_w[0].m[1];
        dx0[2] = zNear_w[0].m[6] - zNear_w[0].m[2];
        dy0[0] = zNear_w[0].m[8] - zNear_w[0].m[4];
        dy0[1] = zNear_w[0].m[9] - zNear_w[0].m[5];
        dy0[2] = zNear_w[0].m[10]- zNear_w[0].m[6];

        x1[0] =  zFar_w[0].m[0];
        x1[1] =  zFar_w[0].m[1];
        x1[2] =  zFar_w[0].m[2];

        x0[0] =  zNear_w[0].m[0];
        x0[1] =  zNear_w[0].m[1];
        x0[2] =  zNear_w[0].m[2];

        double l1_p0[3],l1_p1[3];
        for(int j=0;j<3;++j){
            l1_p0[j] = x0[j] + dx0[j]*a + dy0[j]*b;
            l1_p1[j] = x1[j] + dx1[j]*a + dy1[j]*b;
        }
        double l2_min = 1e10;
        double nearest[3];

        for (int i=0; i<c->imagePoints[1].size();i+=1){
            a=c->imagePoints[1][i].x*1.0/c->images[1].size().width;
            b=1.0-(c->imagePoints[1][i].y*1.0/c->images[1].size().height);
            dx1[0] = zFar_w[1].m[4] - zFar_w[1].m[0];
            dx1[1] = zFar_w[1].m[5] - zFar_w[1].m[1];
            dx1[2] = zFar_w[1].m[6] - zFar_w[1].m[2];
            dy1[0] = zFar_w[1].m[8] - zFar_w[1].m[4];
            dy1[1] = zFar_w[1].m[9] - zFar_w[1].m[5];
            dy1[2] = zFar_w[1].m[10]- zFar_w[1].m[6];

            dx0[0] = zNear_w[1].m[4] - zNear_w[1].m[0];
            dx0[1] = zNear_w[1].m[5] - zNear_w[1].m[1];
            dx0[2] = zNear_w[1].m[6] - zNear_w[1].m[2];
            dy0[0] = zNear_w[1].m[8] - zNear_w[1].m[4];
            dy0[1] = zNear_w[1].m[9] - zNear_w[1].m[5];
            dy0[2] = zNear_w[1].m[10]- zNear_w[1].m[6];

            x1[0] =  zFar_w[1].m[0];
            x1[1] =  zFar_w[1].m[1];
            x1[2] =  zFar_w[1].m[2];

            x0[0] =  zNear_w[1].m[0];
            x0[1] =  zNear_w[1].m[1];
            x0[2] =  zNear_w[1].m[2];

            double l2_p0[3],l2_p1[3];
            for(int j=0;j<3;++j){
                l2_p0[j] = x0[j] + dx0[j]*a + dy0[j]*b;
                l2_p1[j] = x1[j] + dx1[j]*a + dy1[j]*b;
            }
            double ne[3];
            double l2 = findNearestPoint(l1_p0,l1_p1,l2_p0,l2_p1,ne);
            if (l2<l2_min){
                l2_min = l2;
                nearest[0] = ne[0];
                nearest[1] = ne[1];
                nearest[2] = ne[2];
                p_min_0_0.x=l1_p0[0];
                p_min_0_0.y=l1_p0[1];
                p_min_0_0.z=l1_p0[2];

                p_min_0_1.x=l1_p1[0];
                p_min_0_1.y=l1_p1[1];
                p_min_0_1.z=l1_p1[2];

                p_min_1_0.x=l2_p0[0];
                p_min_1_0.y=l2_p0[1];
                p_min_1_0.z=l2_p0[2];

                p_min_1_1.x=l2_p1[0];
                p_min_1_1.y=l2_p1[1];
                p_min_1_1.z=l2_p1[2];
            }

        }
        VEC3 pp;
        pp.x = nearest[0];
        pp.y = nearest[1];
        pp.z = nearest[2];
        points.push_back(pp);
        p0_0.push_back(p_min_0_0);
        p0_1.push_back(p_min_0_1);

        p1_0.push_back(p_min_1_0);
        p1_1.push_back(p_min_1_1);
    }

     glPointSize(7);
    glColor3f(1,0,1);
    glBegin(GL_POINTS);
    for (int i=0; i<points.size(); ++i)
    glVertex3f(points[i].x,points[i].y,points[i].z);
    glEnd();

    glBegin(GL_LINES);
    for (int i=0; i<points.size(); ++i){
        glColor3f(0,0.5,0);
    glVertex3f(p0_0[i].x, p0_0[i].y, p0_0[i].z);
    glVertex3f(p0_1[i].x, p0_1[i].y, p0_1[i].z);
    glColor3f(0.5,0.5,0.5);
    glVertex3f(p1_0[i].x, p1_0[i].y, p1_0[i].z);
    glVertex3f(p1_1[i].x, p1_1[i].y, p1_1[i].z);
    }
    glEnd();

    glEnable(GL_DEPTH_TEST);


}



void drawIntersections(camera_sequence *c){

    std::vector<VEC3> points;
    std::vector<VEC3> p0_0;
    std::vector<VEC3> p0_1;

    std::vector<VEC3> p1_0;
    std::vector<VEC3> p1_1;

    OpenGLMatrix M,InvM, zNear, zFar,zNear_w[2],zFar_w[2];


    for (int nn = 0; nn<2; ++nn){
        gL4x4MatrixMult(&(c->projM),&(c->modelvMatrices[nn]),&M);
        gL4x4MatrixInverse(&M,&InvM);

        zNear.m[0] = -1.0; zNear.m[4] =  1.0; zNear.m[8]  =  1.0;  zNear.m[12] = -1.0;  //x
        zNear.m[1] = -1.0; zNear.m[5] = -1.0; zNear.m[9]  =  1.0;  zNear.m[13] =  1.0;  //y
        zNear.m[2] = -1.0; zNear.m[6] = -1.0; zNear.m[10] = -1.0;  zNear.m[14] = -1.0;  //z
        zNear.m[3] =  1.0; zNear.m[7] =  1.0; zNear.m[11] =  1.0;  zNear.m[15] =  1.0;  //w

        zFar.m[0] = -1.0; zFar.m[4] =  1.0; zFar.m[8]  =  1.0;  zFar.m[12] = -1.0;  //x
        zFar.m[1] = -1.0; zFar.m[5] = -1.0; zFar.m[9]  =  1.0;  zFar.m[13] =  1.0;  //y
        zFar.m[2] =  1.0; zFar.m[6] =  1.0; zFar.m[10] =  1.0;  zFar.m[14] =  1.0;  //z
        zFar.m[7] =  1.0; zFar.m[3] =  1.0; zFar.m[11] =  1.0;  zFar.m[15] =  1.0;  //w

        gL4x4MatrixMult(&InvM,&zNear,&(zNear_w[nn]));
        gL4x4MatrixMult(&InvM,&zFar,&(zFar_w[nn]));

        //dividing by weight
        for (int i=0;i<4;i++){
            for (int j=0;j<3;j++){
                zNear_w[nn].m[i*4+j] = zNear_w[nn].m[i*4+j]/zNear_w[nn].m[i*4+3];
                zFar_w[nn].m[i*4+j] = zFar_w[nn].m[i*4+j]/zFar_w[nn].m[i*4+3];
            }
            zNear_w[nn].m[i*4+3] = 1.0;
            zFar_w[nn].m[i*4+3] = 1.0;
        }
    }
    glDisable(GL_DEPTH_TEST);


    VEC3 p_min_0_0, p_min_1_0,p_min_0_1,p_min_1_1;

    for (int nn=counter; nn<c->imagePoints[0].size();nn+=19){
        double a,b,dx0[3],dy0[3],dx1[3],dy1[3],x0[3],x1[3];
        a=c->imagePoints[0][nn].x*1.0/c->images[0].size().width;
        b=1.0-(c->imagePoints[0][nn].y*1.0/c->images[0].size().height);
        dx1[0] = zFar_w[0].m[4] - zFar_w[0].m[0];
        dx1[1] = zFar_w[0].m[5] - zFar_w[0].m[1];
        dx1[2] = zFar_w[0].m[6] - zFar_w[0].m[2];
        dy1[0] = zFar_w[0].m[8] - zFar_w[0].m[4];
        dy1[1] = zFar_w[0].m[9] - zFar_w[0].m[5];
        dy1[2] = zFar_w[0].m[10]- zFar_w[0].m[6];

        dx0[0] = zNear_w[0].m[4] - zNear_w[0].m[0];
        dx0[1] = zNear_w[0].m[5] - zNear_w[0].m[1];
        dx0[2] = zNear_w[0].m[6] - zNear_w[0].m[2];
        dy0[0] = zNear_w[0].m[8] - zNear_w[0].m[4];
        dy0[1] = zNear_w[0].m[9] - zNear_w[0].m[5];
        dy0[2] = zNear_w[0].m[10]- zNear_w[0].m[6];

        x1[0] =  zFar_w[0].m[0];
        x1[1] =  zFar_w[0].m[1];
        x1[2] =  zFar_w[0].m[2];

        x0[0] =  zNear_w[0].m[0];
        x0[1] =  zNear_w[0].m[1];
        x0[2] =  zNear_w[0].m[2];

        double l1_p0[3],l1_p1[3];
        for(int j=0;j<3;++j){
            l1_p0[j] = x0[j] + dx0[j]*a + dy0[j]*b;
            l1_p1[j] = x1[j] + dx1[j]*a + dy1[j]*b;
        }
        double l2_min = 1e10;
        double nearest[3];

        for (int i=0; i<c->imagePoints[1].size();i+=1){
            a=c->imagePoints[1][i].x*1.0/c->images[1].size().width;
            b=1.0-(c->imagePoints[1][i].y*1.0/c->images[1].size().height);
            dx1[0] = zFar_w[1].m[4] - zFar_w[1].m[0];
            dx1[1] = zFar_w[1].m[5] - zFar_w[1].m[1];
            dx1[2] = zFar_w[1].m[6] - zFar_w[1].m[2];
            dy1[0] = zFar_w[1].m[8] - zFar_w[1].m[4];
            dy1[1] = zFar_w[1].m[9] - zFar_w[1].m[5];
            dy1[2] = zFar_w[1].m[10]- zFar_w[1].m[6];

            dx0[0] = zNear_w[1].m[4] - zNear_w[1].m[0];
            dx0[1] = zNear_w[1].m[5] - zNear_w[1].m[1];
            dx0[2] = zNear_w[1].m[6] - zNear_w[1].m[2];
            dy0[0] = zNear_w[1].m[8] - zNear_w[1].m[4];
            dy0[1] = zNear_w[1].m[9] - zNear_w[1].m[5];
            dy0[2] = zNear_w[1].m[10]- zNear_w[1].m[6];

            x1[0] =  zFar_w[1].m[0];
            x1[1] =  zFar_w[1].m[1];
            x1[2] =  zFar_w[1].m[2];

            x0[0] =  zNear_w[1].m[0];
            x0[1] =  zNear_w[1].m[1];
            x0[2] =  zNear_w[1].m[2];

            double l2_p0[3],l2_p1[3];
            for(int j=0;j<3;++j){
                l2_p0[j] = x0[j] + dx0[j]*a + dy0[j]*b;
                l2_p1[j] = x1[j] + dx1[j]*a + dy1[j]*b;
            }
            double ne[3];
            double l2 = findNearestPoint(l1_p0,l1_p1,l2_p0,l2_p1,ne);
            if (l2<tres){
                l2_min = l2;
                nearest[0] = ne[0];
                nearest[1] = ne[1];
                nearest[2] = ne[2];
                p_min_0_0.x=l1_p0[0];
                p_min_0_0.y=l1_p0[1];
                p_min_0_0.z=l1_p0[2];

                p_min_0_1.x=l1_p1[0];
                p_min_0_1.y=l1_p1[1];
                p_min_0_1.z=l1_p1[2];

                p_min_1_0.x=l2_p0[0];
                p_min_1_0.y=l2_p0[1];
                p_min_1_0.z=l2_p0[2];

                p_min_1_1.x=l2_p1[0];
                p_min_1_1.y=l2_p1[1];
                p_min_1_1.z=l2_p1[2];

                VEC3 pp;
                pp.x = nearest[0];
                pp.y = nearest[1];
                pp.z = nearest[2];
                points.push_back(pp);
                p0_0.push_back(p_min_0_0);
                p0_1.push_back(p_min_0_1);

                p1_0.push_back(p_min_1_0);
                p1_1.push_back(p_min_1_1);
            }

        }

    }

     glPointSize(7);
    glColor3f(1,0,1);
    glBegin(GL_POINTS);
    for (int i=0; i<points.size(); ++i)
    glVertex3f(points[i].x,points[i].y,points[i].z);
    glEnd();

    glBegin(GL_LINES);
    for (int i=0; i<points.size(); ++i){
        glColor3f(0,0.5,0);
    glVertex3f(p0_0[i].x, p0_0[i].y, p0_0[i].z);
    glVertex3f(p0_1[i].x, p0_1[i].y, p0_1[i].z);
    glColor3f(0.5,0.5,0.5);
    glVertex3f(p1_0[i].x, p1_0[i].y, p1_0[i].z);
    glVertex3f(p1_1[i].x, p1_1[i].y, p1_1[i].z);
    }
    glEnd();

    glEnable(GL_DEPTH_TEST);


}
