#pragma once

#include <vector>
#include "windows.h"
//#include "cv.h"

class ParticleMaskCorrelation
{
	double m_CorrelationCoefficientThreshold;
	int m_WindowSide;
	double m_GaussRadius;
public:
	ParticleMaskCorrelation() :	m_CorrelationCoefficientThreshold(0.7), m_WindowSide(5), m_GaussRadius(1) {}
	virtual void SetParams(double correlationCoefficientThreshold, int windowSide, double gaussRadius);
	virtual bool GetMap(IplImage *img, std::vector<CvPoint2D32f> &particles);
};