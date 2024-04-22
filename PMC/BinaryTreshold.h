#pragma once

#include <vector>
#include "windows.h"
//#include "opencv/cv.h"

unsigned int charsToInt(char* b, unsigned length);

class BinaryTreshold
{
	double m_TrshldLevel;
	int m_partSize;
public:
	BinaryTreshold() : m_TrshldLevel(40), m_partSize(10) {}
	virtual void SetParams(double trshldLevel, int partSize) {m_TrshldLevel = trshldLevel;m_partSize = partSize;}
	virtual bool GetMap(IplImage *img, std::vector<CvPoint2D32f> &particles);
};