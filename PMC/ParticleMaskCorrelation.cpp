#include "stdafx.h"
#include "ParticleMaskCorrelation.h"
#include <math.h>
#include <assert.h>
#include "DomainDetection.h"

struct ParValue
{
	double Val;
	tstring Dim;
};

void ParticleMaskCorrelation::SetParams(double correlationCoefficientThreshold, int windowSide, double gaussRadius)
{
	m_CorrelationCoefficientThreshold = correlationCoefficientThreshold;
	m_WindowSide = windowSide;
	m_GaussRadius = gaussRadius;
}

class CoefficientResolver
{
	std::vector<double> gausscoeffs;
	const std::vector<DWORD>& data;
	int width, height;
	int side;
public:
	CoefficientResolver(const std::vector<DWORD> &d, int w, int h, double particleradius, int s = 0):
		width(w), height(h), data(d), side(s)
	{
		if (!side)
		{
			side = (int)((particleradius * 1.5 + 0.5)*2.0);
		}
		int len = side*side;
		//gausscoeffs.reserve(len);
		gausscoeffs.resize(len);
		double x0, y0, square_sigma = particleradius*particleradius, dx, dy, g_av = 0.0;
		x0 = y0 = double(side - 1)/2.0;
		for (int i = 0; i < len; ++i)
		{
			dx = ((i % side) - x0);
			dy = ((i / side) - y0);
			gausscoeffs[i] = exp(-(dx*dx + dy*dy)/2*square_sigma);
			g_av += gausscoeffs[i];
		}
		g_av /= double(len);
		for (int i = 0; i < len; ++i)
		{
			gausscoeffs[i] -= g_av;
		}
	}

	DWORD getdata(int x, int y)
	{
		if (x >= 0 && x < width && y >=0 && y < height)
		{
			return data[y * width + x];
		}
		return 0;
	}

	double operator() (int x, int y)
	{
		static double delta = 1e-10;
		int len = side*side;
		double i_av = 0.0;
		int ai, aj;
		int aj1 = y - side/2, aj2 = y - side/2 + side,
			ai1 = x - side/2, ai2 = x - side/2 + side;
		for ( aj = aj1; aj < aj2; ++aj )
		{
			for ( ai = ai1; ai < ai2; ++ai )
			{
				i_av += getdata(ai, aj);
			}
		}
		i_av /= double(len);
		double up_sum = 0.0, bottom_sum1 = 0.0, bottom_sum2 = 0.0, di;
		for ( aj = aj1; aj < aj2; ++aj )
		{
			for ( ai = ai1; ai < ai2; ++ai )
			{
				di = getdata(ai, aj) - i_av;
				double &c = gausscoeffs[(aj - aj1) * side + (ai - ai1)];
				up_sum += di * c;
				bottom_sum1 += di*di;
				bottom_sum2 += c*c;
			}
		}
		double bottom = bottom_sum1*bottom_sum2;
		if (bottom == 0.0)
		{
			return delta;
		}
		else
		{
			double res = up_sum/sqrt(bottom);
			return (res > 0) ? res : delta;
		}
	}
};

//////////////////////////////////////////////////////////////////////////
//	I1	I2	I3
//	I4	I5	I6
//	I7	I8	I9
//////////////////////////////////////////////////////////////////////////
class PeakResolver
{
	int width, height;
	std::vector<double> & coeffs;
public:
	PeakResolver(int w, int h, std::vector<double> & _c):
		width(w), height(h), coeffs(_c)
	{}
	bool GetPeak(int indx, float &x_pos, float &y_pos)
	{
		int x = indx % width,
			y = indx / width;  
		if ((x == 0) || (x == width - 1) ||
			(y == 0) || (y == height - 1))
		{
			return false;
		}
		int up_indx = indx - width, 
			down_indx = indx + width, 
			left_indx = indx - 1, 
			right_indx = indx + 1;
		float	up_val = log(coeffs[up_indx]), 
				down_val = log(coeffs[down_indx]), 
				left_val = log(coeffs[left_indx]), 
				right_val = log(coeffs[right_indx]),
				peak_val = log(coeffs[indx]);
		x_pos = float(x) + ((left_val - right_val)/(2*left_val - 4*peak_val + 2*right_val));
		y_pos = float(y) + ((up_val - down_val)/(2*up_val - 4*peak_val + 2*down_val));
		return true;
	}
};

bool ParticleMaskCorrelation::GetMap(IplImage *img, std::vector<CvPoint2D32f> &particles)
{
	int i, j;

	int height = img->height;
	int width = img->width;
	int index;

	std::vector<DWORD> data(width * height);

	for ( j = 0; j < height; ++j )
	{
		for ( i = 0; i < width; ++i )
		{
			index = j * width + i;
			data[index] = (DWORD)img->imageDataOrigin[index];
		}
	}

	std::vector<double> coefficient_field;

	coefficient_field.resize(width * height);

	CoefficientResolver cr(data, width, height, m_GaussRadius, m_WindowSide);

	for ( j = 0; j < height; ++j )
	{
		for ( i = 0; i < width; ++i )
		{
			coefficient_field[j * width + i] = cr(i,j);
		}
	}

	Field f(width, height);

	for ( j = 0; j < height; ++j )
	{
		for ( i = 0; i < width; ++i )
		{
			if (coefficient_field[j * width + i] >= m_CorrelationCoefficientThreshold)
			{
				f.data[j * width + i] = _psExist;
			} 
		}
	}
/*
	std::vector<BYTE> gray_field(coefficient_field.size());
	for ( j = 0; j < height; ++j )
	{
		for ( i = 0; i < width; ++i )
		{
			gray_field[j * width + i] = (coefficient_field[j * width + i] + 1.0)*255.0/2.0;
		}
	}
	DWORD dwHigh;
	HANDLE hFile = CreateFile("d:\\temp\\coefficient.raw", GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	WriteFile(hFile, (void*)(gray_field._Myfirst), gray_field.size(), &dwHigh, NULL);
	CloseHandle(hFile);
*/
	std::vector<int> particleIndxs;
	CvPoint2D32f prtcl;
	PeakResolver pr(width, height, coefficient_field);
	int datasize = width * height;

	for ( i = 0; i < datasize; ++i )
	{
		particleIndxs.resize(0);
		ObtainParticle( f, i, particleIndxs );
		int length = (int)particleIndxs.size();
		if (length)
		{
			double max_val = coefficient_field[particleIndxs[0]];
			int imax = particleIndxs[0];
			for ( j = 1; j < length; ++j )
			{
				if (max_val < coefficient_field[particleIndxs[j]])
				{
					max_val = coefficient_field[particleIndxs[j]];
					imax = particleIndxs[j];
				}
			}
			if (pr.GetPeak(imax, prtcl.x, prtcl.y))
			{
				prtcl.y = img->height - prtcl.y;
				particles.push_back(prtcl);
			}
		}
	}

	return true;
}