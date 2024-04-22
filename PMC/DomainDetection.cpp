#include "stdafx.h"
#include "DomainDetection.h"

void ObtainParticle( Field & f, int indx, std::vector<int>& particleIndxs)
{
	if ( 
		!(f.data[indx] & _psExist) ||
		(f.data[indx]& _psAnalyzed)
		)
	{
		return;
	}

	int datasize = (int)f.data.size();
	std::vector<int> toscan, toscannext;
	toscan.push_back(indx);
	f.data[indx] |= _psAnalyzed;
	int toscansize = (int)toscan.size();
	int neighb_indx;
	while (toscansize)
	{
		for ( int i = 0; i < toscansize; ++i )
		{
			int & currindx = toscan[i];
			neighb_indx = currindx - f.width;
			if ( (neighb_indx >= 0) && 
				(f.data[neighb_indx] & _psExist) && 
				!(f.data[neighb_indx] & _psAnalyzed)	)
			{
				f.data[neighb_indx] |= _psAnalyzed;
				toscannext.push_back(neighb_indx);
			}
			neighb_indx = currindx + f.width;
			if ( (neighb_indx < datasize) &&
				(f.data[neighb_indx] & _psExist) && 
				!(f.data[neighb_indx] & _psAnalyzed)	)
			{
				f.data[neighb_indx] |= _psAnalyzed;
				toscannext.push_back(neighb_indx);
			}
			neighb_indx = currindx - 1;
			if ( (currindx%f.width != 0) && 
				(f.data[neighb_indx] & _psExist) && 
				!(f.data[neighb_indx] & _psAnalyzed)	)
			{
				f.data[neighb_indx] |= _psAnalyzed;
				toscannext.push_back(neighb_indx);
			}
			neighb_indx = currindx + 1;
			if ( (currindx%f.width != (f.width - 1)) && 
				(f.data[neighb_indx] & _psExist) && 
				!(f.data[neighb_indx] & _psAnalyzed)	)
			{
				f.data[neighb_indx] |= _psAnalyzed;
				toscannext.push_back(neighb_indx);
			}
		}
		particleIndxs.insert(particleIndxs.end(), toscan.begin(), toscan.end());
		toscan.resize(0);
		toscan = toscannext;
		toscannext.resize(0);
		toscansize = (int)toscan.size();
	}
}
