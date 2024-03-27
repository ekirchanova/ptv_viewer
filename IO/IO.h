#ifndef IO_H
#define IO_H

#include <stdio.h>
#include "GUI/GUI.h"
//bmp

//ppm
void screenshot_ppm(const char *filename, unsigned int width,
                           unsigned int height, GLubyte *pixels) ;
void screenshot_bmp(char *fileName, unsigned int width, unsigned int height);

#endif
