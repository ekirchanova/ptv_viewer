#include "GUI.h"
#include <math.h>
void draw_sphere(double rad)
{
    //glEnable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(1,1,1,0.65);
    glutSolidSphere(rad,
                    40,40);

    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);

    glColor3f(1,0,0);
    glBegin(GL_LINE_LOOP);
    int i=0;
    for (i=0;i<100;i++)
        glVertex3f(rad*sin(i*2*3.14159265/100),0,rad*cos(i*2*3.14159265/100));
    glEnd();
}
void draw_cube(double rad)
{
    glColor3f(0.5,0.,0);
    glBegin(GL_LINE_LOOP);
    glVertex3f(-rad,rad,rad);
    glVertex3f(-rad,-rad,rad);
    glVertex3f( rad,-rad,rad);
    glVertex3f( rad,rad,rad);
    glEnd();

    glColor3f(0.0,0.5,0);
    glBegin(GL_LINE_LOOP);
    glVertex3f(-rad,rad,-rad);
    glVertex3f(-rad,-rad,-rad);
    glVertex3f( rad,-rad,-rad);
    glVertex3f( rad,rad,-rad);
    glEnd();

    glBegin(GL_LINES);
    glColor3f(0.5,0.,0);
    glVertex3f(-rad,rad,rad);
    glColor3f(0.0,0.5,0);
    glVertex3f(-rad,rad,-rad);

    glColor3f(0.5,0.,0);
    glVertex3f(-rad,-rad,rad);
    glColor3f(0.0,0.5,0);
    glVertex3f(-rad,-rad,-rad);

    glColor3f(0.5,0.,0);
    glVertex3f( rad,-rad,rad);
    glColor3f(0.0,0.5,0);
    glVertex3f( rad,-rad,-rad);

    glColor3f(0.5,0.,0);
    glVertex3f( rad,rad,rad);
    glColor3f(0.0,0.5,0);
    glVertex3f( rad,rad,-rad);
    glEnd();
}

void draw_cone(VEC3 p1, VEC3 p2,double r1,double r2,int close_p1, int close_p2, int J_NUM,VEC3 col)
{
    double l2,l, alpha;
    VEC3 n, dx, cross_l, cross_u, parr[2][10], nrm[10];
    n.x = 1.0/sqrt(3.0); n.y = n.x; n.z = n.x;

    //get cone surface with some vector algebra
    dx.x = p2.x - p1.x;
    dx.y = p2.y - p1.y;
    dx.z = p2.z - p1.z;

    l=sqrt(dx.x*dx.x + dx.y*dx.y + dx.z*dx.z);
    dx.x /= l;
    dx.y /= l;
    dx.z /= l;
    //first normal:
    //both vectors have unit lengths
    cross_l.z =  n.x*dx.y - n.y*dx.x;
    cross_l.y = -n.x*dx.z + n.z*dx.x;
    cross_l.x =  n.y*dx.z - n.z*dx.y;

    l2=sqrt(cross_l.x*cross_l.x + cross_l.y*cross_l.y + cross_l.z*cross_l.z);

    cross_l.x /= l2;
    cross_l.y /= l2;
    cross_l.z /= l2;
    //binormal:
    //dx and cross_l both unit and orthogonal
    cross_u.z = -(cross_l.x*dx.y - cross_l.y*dx.x);
    cross_u.y = -(-cross_l.x*dx.z + cross_l.z*dx.x);
    cross_u.x = -(cross_l.y*dx.z - cross_l.z*dx.y);


    for (int j=0; j <= J_NUM; j++)
    {
        double s_n,c_s;
        alpha=2*M_PI*j*1.0/J_NUM;
        s_n=sin(alpha);
        c_s=cos(alpha);

        n.x= cross_u.x*c_s - s_n*cross_l.x;
        n.y= cross_u.y*c_s - s_n*cross_l.y;
        n.z= cross_u.z*c_s - s_n*cross_l.z;

        parr[0][j].x = p1.x + n.x*r1;
        parr[0][j].y = p1.y + n.y*r1;
        parr[0][j].z = p1.z + n.z*r1;
        parr[1][j].x = p2.x + n.x*r2;
        parr[1][j].y = p2.y + n.y*r2;
        parr[1][j].z = p2.z + n.z*r2;
        //local cone surface normal:
        double nl = sqrt((r1-r2)*(r1-r2) + l*l);
        nrm[j].x = (n.x*l + dx.x*(r1-r2)) / nl;
        nrm[j].y = (n.y*l + dx.y*(r1-r2)) / nl;
        nrm[j].z = (n.z*l + dx.z*(r1-r2)) / nl;
    }

    //start drawing
    GLfloat col0[4],col1[4];
    col0[0]=col.x; col0[1]=col.y; col0[2]=col.z; col0[3] =1.0;
    col1[0]=1; col1[1]=1; col1[2]=1; col1[3] =1.0;

    glEnable(GL_COLOR_MATERIAL);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, col0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, col1);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 2.5f);
    glColor4fv(col0);
    glBegin(GL_TRIANGLE_STRIP);
    for (int j=0;j<J_NUM;j++)
    {
        glNormal3f(nrm[j].x, nrm[j].y, nrm[j].z);
        glVertex3f(parr[0][j].x, parr[0][j].y, parr[0][j].z);
        glVertex3f(parr[1][j].x, parr[1][j].y, parr[1][j].z);
    }
    glNormal3f(nrm[0].x, nrm[0].y, nrm[0].z);
    glVertex3f(parr[0][0].x, parr[0][0].y, parr[0][0].z);
    glVertex3f(parr[1][0].x, parr[1][0].y, parr[1][0].z);
    glEnd();

    if (close_p1){
        glBegin(GL_TRIANGLE_FAN);
        glNormal3f(-dx.x, -dx.y, -dx.z);
        glVertex3f(p1.x,p1.y,p1.z);
        for (int j=0;j<J_NUM;j++)
            glVertex3f(parr[0][j].x, parr[0][j].y, parr[0][j].z);
        glVertex3f(parr[0][0].x, parr[0][0].y, parr[0][0].z);
        glEnd();
    }

    if (close_p2){
        glBegin(GL_TRIANGLE_FAN);
        glNormal3f(dx.x, dx.y, dx.z);
        glVertex3f(p2.x,p2.y,p2.z);
        for (int j=0;j<J_NUM;j++)
            glVertex3f(parr[1][j].x, parr[1][j].y, parr[1][j].z);
        glVertex3f(parr[1][0].x, parr[1][0].y, parr[1][0].z);
        glEnd();
    }


}

void draw_arrow(VEC3 p1, VEC3 p2, VEC3 col){
    VEC3 dx, pm;
    dx.x = p2.x - p1.x;
    dx.y = p2.y - p1.y;
    dx.z = p2.z - p1.z;

    double l=sqrt(dx.x*dx.x + dx.y*dx.y + dx.z*dx.z);
    /*dx.x /= l;
    dx.y /= l;
    dx.z /= l;*/
    double a=0.5;
    pm.x = p1.x + dx.x*a;
    pm.y = p1.y + dx.y*a;
    pm.z = p1.z + dx.z*a;
    double r0  = l*0.1;
    draw_cone(p1,pm,r0,r0,1,0,5,col);
    draw_cone(pm,p2,2.0*r0,1e-8,1,0,5,col);
}

void draw_target(){
    int ni,nj,n_phi;
    ni = 19; nj = 19; n_phi = 20;
    double dr, r0, rc;
    r0 = 0.01;
    dr = 5.1*r0;
    rc = 1.5*r0;
    glColor3f(1.0,0.0,0.0);
    for (int i=0; i < ni; ++i)
        for (int j=0; j < nj; ++j){
            glBegin(GL_TRIANGLE_FAN);
                    glVertex3f(dr*(i - (ni-1)*0.5) , dr*(j - (nj-1)*0.5) , 0.0);
            for (int n=0; n <= n_phi ; ++n){
                if (!(i == (ni-1)/2 && j == (nj-1)/2))
                    glVertex3f(dr*(i - (ni-1)*0.5) + r0*cos(n*2.0*M_PI/n_phi), dr*(j - (nj-1)*0.5) + r0*sin(n*2.0*M_PI/n_phi), 0.0);
                else
                    glVertex3f(dr*(i - (ni-1)*0.5) + rc*cos(n*2.0*M_PI/n_phi), dr*(j - (nj-1)*0.5) + rc*sin(n*2.0*M_PI/n_phi), 0.0);
            }
            glEnd();
        }
}

void draw_origin(double s){
    glBegin(GL_LINES);
    glColor3f(1., 0., 0.);
    glVertex3f(0., 0., 0.);
    glVertex3f(s, 0., 0.);
    glColor3f(0., 1., 0.);
    glVertex3f(0., 0., 0.);
    glVertex3f(0., s, 0.);
    glColor3f(0., 0., 1.);
    glVertex3f(0., 0., 0.);
    glVertex3f(0., 0., s);
    glEnd();
}


void draw_points(const std::vector<cv::Point2f>& points, double z, GLdouble red, GLdouble green, GLdouble blue)
{
    size_t size = points.size();

    printf("found %d points\n",size);
    glPointSize(5);
    glBegin(GL_POINTS);
    glColor3d(red,green,blue);
    for(size_t i = 0; i < size; ++i)
    {
        glVertex3d(points[i].x, points[i].y, z);
        printf("points %i x : %f y :%d\n", i, points[i].x, points[i].y );

    }
    glEnd();
}
