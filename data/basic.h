#ifndef BASIC_H
#define BASIC_H

/*
#define X 0
#define Y 1
#define Z 2*/

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define CROSS(dest,v1,v2) \
    dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
    dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
    dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) \
    dest[0]=v1[0]-v2[0]; \
    dest[1]=v1[1]-v2[1]; \
    dest[2]=v1[2]-v2[2];

#define SUM(dest,v1,v2) \
    dest[0]=v1[0]+v2[0]; \
    dest[1]=v1[1]+v2[1]; \
    dest[2]=v1[2]+v2[2];

#define SMULT(dest,s,v) \
    dest[0]=(s)*v[0]; \
    dest[1]=(s)*v[1]; \
    dest[2]=(s)*v[2];


#include <vector>
#include <opencv2/core/mat.hpp>

typedef struct{
    GLdouble m[16];
}OpenGLMatrix;

typedef struct{
    //int w,h; //width ang height
    double focal_length, fov; //fov is  x_g = z_l*fov*x_l, y_g=z_l*fov*y_l
    double c[3]; //word coordinates of center of focal plane
    double f[3],u[3],r[3];//forward up and right vectors basis. all should have unit length
    double dz; //center of the image coordinate offset in local coordinates
    //int n; //number of snapshots
   /* unsigned char * data; //rgb data
    unsigned int texture_index;
    double points[100024][3];
    int pointSize;
    double centers[10000][3];
    int centers_size;*/
    std::vector<unsigned int> texture_indices;
    std::vector<cv::Mat> rvects;
    std::vector<cv::Mat> tvects;
    cv::Mat cameraMatrix;
    cv::Mat distCoeffs;
    std::vector<std::vector<cv::Point2f> > imagePoints;
    std::vector<std::vector<cv::Point2f> > undistortImagePoints;
    std::vector<std::vector<cv::Point3f> > objectPoints;
    std::vector<cv::Mat> images;
    OpenGLMatrix projM;
    std::vector<OpenGLMatrix> modelvMatrices;
}camera_sequence;
#endif
