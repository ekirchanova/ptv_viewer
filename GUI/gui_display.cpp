#include "GUI.h"

int redr = 0; //for automatic redrawing
int active_camera =0;
void display(void){ //main display loop
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glLoadIdentity();

    if (mode == 0){
    glTranslatef(t_x,t_y,t_z);
    glRotatef(ry,1.0,0,0);
    glRotatef(rx,0,1.0,0);
    glScalef(sc,sc,sc);
    }
    else //looking from cameras' view
        glLoadMatrixd(cs.modelvMatrices[active_camera].m);

    glColor3f(1,0,0);

    glBegin(GL_LINE_LOOP);
    glVertex2f(w_x0,w_y0);
    glVertex2f(w_x0,w_y1);
    glVertex2f(w_x1,w_y1);
    glVertex2f(w_x1,w_y0);
    glEnd();

     draw_origin(0.1);
     draw_target();


    if(mode==0)
   drawCameraFrustum(&cs,active_camera);


    glutSwapBuffers();
    if (redr==1) glutPostRedisplay();
}
