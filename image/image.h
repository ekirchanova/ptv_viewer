#ifndef IMAGE_H
#define IMAGE_H

#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>

#include "data/basic.h"

int prepareCalibImageCV(const char * fileName, camera_sequence *c);
void try_calibration(camera_sequence *c);

void convertOpenCVModelViewMatrixToOpenGL(cv::Mat rvec, cv::Mat tvec, GLdouble* GLModelV );



void convertOpenCVProjectionMatrixToOpenGL(cv::Mat cameraMatrix, double zN,double zF,GLdouble* GLProjM);


#endif
