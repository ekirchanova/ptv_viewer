#ifndef GUI_H
#define GUI_H

#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>

#include "data/basic.h"

#define DATA_OFFSET_OFFSET 0x000A
#define WIDTH_OFFSET 0x0012
#define HEIGHT_OFFSET 0x0016
#define BITS_PER_PIXEL_OFFSET 0x001C
#define HEADER_SIZE 14
#define INFO_HEADER_SIZE 40
#define NO_COMPRESION 0
#define MAX_NUMBER_OF_COLORS 0
#define ALL_COLORS_REQUIRED 0

typedef unsigned int int32;
typedef short int16;
typedef unsigned char byte;

typedef struct
{
    int num;
    float r[255],g[255],b[255];
    int x[255];
    int x_num;
    int xt[1024];
    float dx;
}color_table;


extern double tres;
extern int counter;
extern int active_camera;
extern char dir[1024];

extern unsigned int W_HEIGHT; //window height
extern unsigned int W_WIDTH; //window height
extern double w_x0; //screen corners space coordinates
extern double w_x1; //screen corners space coordinates
extern double w_y0; //screen corners space coordinates
extern double w_y1;//screen corners space coordinates


extern float rx; //camera angle1
extern float ry; //camera angle2
extern int mx0,my0;//previous mouse coordinates

extern float rx0; //prev camera angle
extern float ry0; //
extern float t_x; //view translation vector
extern float t_y; //
extern float t_z;
extern float sc;

extern GLuint textures[1024];
extern camera_sequence cs;

extern std::vector<VEC3> tracers;
extern std::vector<VEC3> vort;

extern int mode;
void loadTextureToVideoMemory(camera_sequence *c);
int loadBmpImage(const char * fileName, camera_sequence *c);
GLuint LoadBMPTexture( const char * filename ); //textures

void GUI_init(int argc, char **argv); //window initialization
void resize(int w, int h); //on resize event
void kb(unsigned char key, int x, int y); //keyboard event
void skb(int key, int x, int y); //special keys event
void m_d(int button, int state,int x, int y);  //mouse down event
void m_m(int x,int y); //mouse move

void display(void); //main dispaly redraw event


void gL4x4MatrixInverse(OpenGLMatrix *M, OpenGLMatrix *Inv);
void drawCameraFrustum(camera_sequence *c, int active_image);
void drawIntersections(camera_sequence *c);

void draw_target();
void draw_origin(double s);
void draw_tracers(std::vector<VEC3> & tracers,double rad,int n_rad,int n_phi);
void evolve_tracers(std::vector<VEC3> & vortex, std::vector<VEC3> & tracers, double dt);
#endif
