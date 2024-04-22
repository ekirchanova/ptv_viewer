#include "stdafx.h"
#include "BinaryTreshold.h"
#include "DomainDetection.h"
#include <math.h>

bool BinaryTreshold::GetMap(IplImage *img, std::vector<CvPoint2D32f> &particles)
{
	int width	= img->width;
	int height	= img->height;
	Field f(width, height);

	int i, j;

	int maxIntensityValue = 1 << img->depth - 1;
	unsigned int TrshldLevel = m_TrshldLevel;//(unsigned int)( maxIntensityValue * m_TrshldLevel + 0.5);

	/*for ( i = 0; i < f.height; ++i )
	{
	for ( j = 0; j < f.width; ++j )
	{
	DWORD ghgh = (DWORD)img->imageDataOrigin[ i * f.width + j ];
	if ( ghgh > (DWORD)TrshldLevel )
	{
	f.data[i*f.width + j] |= _psExist;
	}
	}
	}*/
	unsigned int PixelSize = img->depth/CHAR_BIT;
	for ( i = 0; i < f.height; ++i )
	{
		for ( j = 0; j < f.width; ++j )
		{
			DWORD ghgh = charsToInt(&img->imageDataOrigin[ i * f.width*PixelSize + j*PixelSize ],PixelSize);
			if ( ghgh > (DWORD)TrshldLevel )
			{
				f.data[i*f.width + j] |= _psExist;
			}
		}
	}

	std::vector<int> particleIndxs;
	CvPoint2D32f prtcl;
	int datasize = f.width * f.height;
	int length;
	double dbllength;
	for ( i = 0; i < datasize; ++i )
	{
		particleIndxs.resize(0);
		ObtainParticle( f, i, particleIndxs );
		length = (int)particleIndxs.size();
		if ( length > m_partSize )
		{
			dbllength = length;
			prtcl.x = 0.0;
			prtcl.y = 0.0;
			for ( j = 0; j < length; ++j )
			{
				prtcl.x += ((float)(particleIndxs[j] % f.width )) / dbllength;
				prtcl.y += ((float)(particleIndxs[j] / f.width )) / dbllength;
			}
			prtcl.y = img->height - 1 - prtcl.y;
			particles.push_back(prtcl);
		}
	}
	return true;
}

unsigned int charsToInt(char* b, unsigned length)
{
	int val = 0;
	int j = 0;
	for (int i = length-1; i >= 0; --i)
	{
		val += (b[i] & 0xFF) << (8*j);
		++j;
	}
	return val;
}