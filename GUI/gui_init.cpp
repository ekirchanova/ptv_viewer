#include "GUI.h"
#include <stdlib.h>
#include <stdio.h>
#include "image/image.h"
unsigned int W_HEIGHT=900;
unsigned int W_WIDTH=900;
double w_x0=-1.3;
double w_x1=1.3;

double w_y0=-1.3;
double w_y1=1.3;
std::vector<VEC3> tracers;
std::vector<VEC3> vort;

void GUI_init(int argc, char **argv){  //window initialization
    //test directory
    printf("GUI dir = %s \n", dir);
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(W_WIDTH,W_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("PTV_viewer");
    glutDisplayFunc(display);
    glutReshapeFunc(resize);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    glutSpecialFunc(skb);

    //loading textures
    char fname[8][1024];
   // sprintf(fname,"%s/test_in/ImperxCamera0_18_10_123/image0.bmp",dir);
    // sprintf(fname,"image0.bmp",dir);
    //printf("textures0 before = %d \n",textures[0]);
    //textures[0] = LoadBMPTexture(fname);
    //printf("textures0 after = %d \n",textures[0]);

    sprintf(fname[0],"%s/test_out/0.bmp",dir);
    //sprintf(fname[0],"%s/test_in/3_1.bmp",dir);
    sprintf(fname[1],"%s/test_out/1.bmp",dir);
    sprintf(fname[2],"%s/test_out/2.bmp",dir);
    sprintf(fname[3],"%s/test_out/3.bmp",dir);
    sprintf(fname[4],"%s/test_out/4.bmp",dir);
    sprintf(fname[5],"%s/test_out/5.bmp",dir);
    sprintf(fname[6],"%s/test_out/6.bmp",dir);
    sprintf(fname[7],"%s/test_out/7.bmp",dir);

    cs.dz = 1.0;
    cs.fov =0.45;

    cs.f[0] = 0.0; cs.f[1] = 0.0; cs.f[2] = -1.0;
    cs.u[0] = 0.0; cs.u[1] = 1.0; cs.u[2] = 0.0;
    cs.r[0] = 1.0; cs.r[1] = 0.0; cs.r[2] = 0.0;

    for (int i=0;i<1;i++){

    prepareCalibImageCV(fname[i],&(cs));

     }
    loadTextureToVideoMemory(&(cs));

    try_calibration(&(cs));
    glEnable(GL_TEXTURE);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(w_x0, w_x1, W_HEIGHT*w_y0/W_WIDTH,W_HEIGHT*w_y1/W_WIDTH, -10.0, 10.0);
    glMatrixMode (GL_MODELVIEW);


    double r_max=0.5;
    for (int i=0; i<500; i++)
    {
        VEC3 tracer;
        tracer.x = (2.0*(rand()*1.0/RAND_MAX) - 1.0)*r_max;
        tracer.y = (2.0*(rand()*1.0/RAND_MAX) - 1.0)*r_max;
        tracer.z = (2.0*(rand()*1.0/RAND_MAX) - 1.0)*r_max;
        tracers.push_back(tracer);
    }

    int n_vort =20;
    double r_vort =0.5;
    for (int i=0; i<n_vort; i++)
    {
        VEC3 tracer;
        tracer.x = 0.5;
        tracer.y = r_vort * sin(2*i*M_PI/(n_vort-1));
        tracer.z = r_vort * cos(2*i*M_PI/(n_vort-1));
        vort.push_back(tracer);
    }

    glutMainLoop();
}

void resize(int w, int h){ //on resize event
    glViewport(0, 0, w, h);
    glClearColor (1.0, 1.0, 1.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    W_HEIGHT  = h; W_WIDTH = w;
    glOrtho(w_x0, w_x1, W_HEIGHT*w_y0/W_WIDTH,W_HEIGHT*w_y1/W_WIDTH, -10.0, 10.0);
    glMatrixMode (GL_MODELVIEW);
}

