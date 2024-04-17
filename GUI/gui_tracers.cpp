#include "GUI/GUI.h"

void get_current_screen_basis(double up[3], double right[3], double forward[3]){
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetIntegerv(GL_VIEWPORT, viewport);

    //1. convert line to screen space
    double wScr,hScr;
    wScr = viewport[2];
    hScr = viewport[3];
    double c_x, c_y; //center of the screen coordinates
    c_x = viewport[2] / 2;
    c_y = viewport[3] / 2;


    double r0[3];
    double r1[3];
    gluUnProject(c_x, c_y, 0.0, modelview, projection, viewport, &(r0[0]), &(r0[1]), &(r0[2]));
    gluUnProject(c_x, c_y+1, 0.0, modelview, projection, viewport, &(r1[0]), &(r1[1]), &(r1[2]));
    SUB(up,r1,r0);
    gluUnProject(c_x+1, c_y, 0.0, modelview, projection, viewport, &(r1[0]), &(r1[1]), &(r1[2]));
    SUB(right,r1,r0);
    double zlen = 1.0/sqrt(DOT(up,up));
    SMULT(up,zlen,up);
    zlen = 1.0/sqrt(DOT(right,right));
    SMULT(right,zlen,right);
    CROSS(forward,up,right);
}

//tracer coordinates, their radius, and number of triangles in both radial and azimuthal directions
void draw_tracers(std::vector<VEC3> & tracers,double rad,int n_rad,int n_phi){
      glDisable(GL_LIGHTING);
   double up[3],right[3],forward[3];
    get_current_screen_basis(up,right,forward);

    for (int i=0; i<tracers.size(); ++i){
        glBegin(GL_TRIANGLE_FAN);
         glColor3f(1.0,1.0,1.0);
         glVertex3f(tracers[i].x,tracers[i].y,tracers[i].z);
         for (int ph=0; ph<n_phi; ++ph){
             double phi = ph*2.0*M_PI/(n_phi-1);
             glColor3f(0.0,0.0,0.0);
           glVertex3f(tracers[i].x + rad*(up[0]*sin(phi) + right[0]*cos(phi)),
                      tracers[i].y + rad*(up[1]*sin(phi) + right[1]*cos(phi)),
                      tracers[i].z + rad*(up[2]*sin(phi) + right[2]*cos(phi)));
         }
       glEnd();
    }


  /*  for (int i=0; i<tracers.size(); ++i){
        glBegin(GL_LINES);
         glColor3f(1.0,0.0,0.0);
         glVertex3f(tracers[i].x,tracers[i].y,tracers[i].z);

         glVertex3f(tracers[i].x +    rad*right[0],
                      tracers[i].y +  rad*right[1],
                      tracers[i].z +  rad*right[2]);

         glColor3f(0.0,1.0,0.0);
         glVertex3f(tracers[i].x,tracers[i].y,tracers[i].z);

         glVertex3f(tracers[i].x +    rad*up[0],
                      tracers[i].y +  rad*up[1],
                      tracers[i].z +  rad*up[2]);

       glEnd();
    }*/

 /*  glBegin(GL_POINTS);
    for (int i=0; i<tracers.size(); ++i){
         glColor3f(0.0,1.0,0.0);
         glVertex3f(tracers[i].x,tracers[i].y,tracers[i].z);
    }
       glEnd();*/
}

void get_velocity_from_vortex_line(std::vector<VEC3> & vortex, double r[3], double vel[3]){
    vel[0]=0.0;
    vel[1]=0.0;
    vel[2]=0.0;

    double dr[3];
    double dl[3];
    double v[3];

    for (int i=0; i<vortex.size()-1; ++i){
        dl[0] = vortex[i+1].x - vortex[i].x;
        dl[1] = vortex[i+1].y - vortex[i].y;
        dl[2] = vortex[i+1].z - vortex[i].z;

        dr[0] = 0.5*(vortex[i+1].x + vortex[i].x) - r[0];
        dr[1] = 0.5*(vortex[i+1].y + vortex[i].y) - r[1];
        dr[2] = 0.5*(vortex[i+1].z + vortex[i].z) - r[2];
        double zd3r = sqrt(DOT(dr,dr));
        zd3r = 1.0/(zd3r*zd3r*zd3r+0.001);
        CROSS(v,dl,dr);
        SMULT(v,zd3r,v);
        SUM(vel,v,vel);
    }
}

void evolve_tracers(std::vector<VEC3> & vortex, std::vector<VEC3> & tracers, double dt)
{
    static VEC3 vels[500];

    for (int i=0; i<tracers.size(); ++i){
        double vel[3];
        double r[3];
        r[0] = tracers[i].x; r[1] = tracers[i].y; r[2] = tracers[i].z;
        get_velocity_from_vortex_line(vortex,r,vel);
        vels[i].x = vel[0]; vels[i].y = vel[1]; vels[i].z = vel[2];
    }

    for (int i=0; i<tracers.size(); ++i){
        tracers[i].x += vels[i].x*dt;
        tracers[i].y += vels[i].y*dt;
        tracers[i].z += vels[i].z*dt;
    }
}
