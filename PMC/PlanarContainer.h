#pragma once

#include <set>
#include <vector>
#include <limits>

class Particle
{
public:
	double x;
	double y;
};

bool operator < (const Particle& p1, const Particle& p2);

//typedef uint unsigned int;

template<class Target>
class PlanarContainer
{
	int m_BWidth, m_BHeight;
	double m_width, m_height;
	double m_side;
	static double reminder;

	std::vector<Target> m_targets;
	std::vector<std::vector<int> > m_buckets;

	int GetBucketIndx( double/*&*/ x, double/*&*/ y );
	void GetInvolvedBuckets( double& x, double& y, double& radius, std::vector<int>& bucketIndxs );
	bool IsBucketOverlap( int indx, double x, double y, double radius );

public:
	class NoTarget : public std::exception
	{
	public:
		NoTarget():std::exception("No target."){}
	};

	PlanarContainer();
	PlanarContainer(double w, double h, double s);
	~PlanarContainer(void);
	void Init(double w, double h, double s);
	double GetWidth() {return m_width;};
	double GetHeight() {return m_height;};
	void Add( Target& point );
//	double GetDensity( double x, double y );
	double GetOverlappedBucketsArea (double& x, double& y, double& radius )
	{
		std::vector<int> bucketIndxs;
		GetInvolvedBuckets(x, y, radius, bucketIndxs);
		return bucketIndxs.size()*m_side*m_side;
	}
	void GetInvolvedBucketsTargets (double x, double y, double radius, std::vector<int>& bounded )
	{
		bounded.resize(0);
		std::vector<int> bucketsindxs;
		GetInvolvedBuckets(x, y, radius, bucketsindxs);
		size_t length = bucketsindxs.size();
		for (size_t i = 0; i < length; ++i)
		{
			int &indx = bucketsindxs[i];
			if ( indx >= (int)m_buckets.size() || indx < 0	)
			{
				continue;
			}
			std::vector<int>& curr = m_buckets[indx];
			bounded.insert(bounded.end(), curr.begin(),curr.end());
		}
	}
	const std::vector<int>& GetBucket( double x, double y );
	const std::vector<int>& GetBucket( int indx )
	{
		return m_buckets[indx];
	}
	int GetBucketCount( )
	{
		return m_buckets.size();
	}
	size_t size()
	{
		return m_targets.size();
	}

	__declspec(property(get=get_targets)) const Target targets[];
	const Target& get_targets(int i) 
	{
		if ( i < 0 || i >= (int)m_targets.size() )
		{
			throw std::exception();
		}
		return m_targets[i];
	}

	__declspec(property(get=get_targets_non_const)) Target targets_non_const[];
	Target& get_targets_non_const(int i) 
	{
		if ( i < 0 || i >= (int)m_targets.size() )
		{
			throw std::exception();
		}
		return m_targets[i];
	}

	int GetNearest( double x, double y, double rad = DBL_MAX );
	void GetCircleBounded( double x, double y, double radius, std::vector<int>& bounded );
};

#include <math.h>

inline double distance(double x1, double y1, double x2, double y2)
{
	return sqrt( ( x1 - x2 )*( x1- x2 ) + ( y1 - y2 )*( y1 - y2 ));
}

inline bool pointinside(double rcLeft, double rcRight, double rcBottom, double rcTop, double x, double y)
{
	return ( ( x < rcRight ) && ( x > rcLeft) && ( y < rcTop ) && ( y > rcBottom) );
}

template<class Target>
double PlanarContainer<Target>::reminder = .999999999999;

template<class Target>
PlanarContainer<Target>::PlanarContainer(double w, double h, double s)
{
	Init( w, h, s);
}

template<class Target>
PlanarContainer<Target>::PlanarContainer():
	m_width(), m_height(), m_BHeight(), m_BWidth(), m_side()
{}

template<class Target>
PlanarContainer<Target>::~PlanarContainer(void){};

template<class Target>
void PlanarContainer<Target>::Init(double w, double h, double s)
{
	m_buckets.clear();
	m_width = w;
	m_height = h;
	m_side = s;
	m_BWidth = (int)(w/m_side + reminder);
	m_BHeight = (int)(h/m_side + reminder);
	m_buckets.resize(m_BWidth * m_BHeight);
}
/*
template<class Target>
inline double PlanarContainer<Target>::GetDensity( double x, double y )
{

}
*/

template<class Target>
inline int PlanarContainer<Target>::GetBucketIndx( double/*&*/ x, double/*&*/ y )
{
	int ix = (int)( x/m_side ),
		iy = (int)( y/m_side );
	return ( iy *m_BWidth + ix );
}

template<class Target>
void PlanarContainer<Target>::Add( Target& point )
{
	m_targets.push_back(point);
	m_buckets[GetBucketIndx(point.x, point.y)].push_back((int)m_targets.size()-1);
}

template<class Target>
const std::vector<int>& PlanarContainer<Target>::GetBucket( double x, double y )
{
	return m_buckets[GetBucketIndx(x, y)];
}

template<class Target>
int PlanarContainer<Target>::GetNearest( double x, double y, double rad )
{
	int indx = GetBucketIndx( x, y);
	double rcframeleft, rcframeright, rcframetop, rcframebottom;
	int xi = indx % m_BWidth, yi = indx / m_BWidth;
	rcframeleft = xi * m_side;
	rcframeright = xi * m_side + m_side;
	rcframebottom = yi * m_side;
	rcframetop = yi * m_side + m_side;

	std::set<double> radiusArray;
	radiusArray.insert(abs(rcframetop - y));
	radiusArray.insert(abs(y - rcframebottom));
	radiusArray.insert(abs(rcframeright - x));
	radiusArray.insert(abs(x - rcframeleft));

	std::set<double>::iterator ptrrad;
	std::vector<int> bounded;

	for( ptrrad = radiusArray.begin(); ptrrad != radiusArray.end() && bounded.empty(); ++ptrrad)
	{
		if (*ptrrad > rad)
		{
			throw NoTarget();
		}
		GetCircleBounded( x, y, *ptrrad, bounded);
	}

	double currRadius = *(radiusArray.begin()) + m_side;
	bool bEnd = false;
	while (bounded.empty() && !bEnd)
	{
		if ((x - currRadius) < 0)
		{
			currRadius = x;
			bEnd = true;
		}
		if ((x + currRadius) > m_width)
		{
			currRadius = m_width - x;
			bEnd = true;
		}
		if ((y - currRadius) < 0)
		{
			currRadius = y;
			bEnd = true;
		}
		if ((y + currRadius) > m_height)
		{
			currRadius = m_height - y;
			bEnd = true;
		}
		if (currRadius > rad)
		{
			currRadius = rad;
			bEnd = true;
		}
		GetCircleBounded( x, y, currRadius, bounded);
		currRadius += m_side;
	}

	if (bEnd && bounded.empty())
	{
		throw NoTarget();
	}

	std::vector<int>::iterator ptr = bounded.begin();
	int res = *ptr;
	Target t0 = m_targets[res];
	double minRadius = distance( t0.x, t0.y, x, y );

	for (; ptr != bounded.end(); ++ptr)
	{
		const Target& t = m_targets[*ptr];
		double currrad = distance( t.x, t.y, x, y );
		if ( currrad < minRadius )
		{
			res = *ptr;
			minRadius = currrad;
		}
	}
	return res;
}

template<class Target>
void PlanarContainer<Target>::GetCircleBounded( double x, double y, double radius, std::vector<int> &bounded )
{
	bounded.resize(0);
	//std::vector<int> previous;
	std::vector<int> bucketsindxs;
	GetInvolvedBuckets(x, y, radius, bucketsindxs);
	size_t length = bucketsindxs.size();
	//std::cout << "ln = "<< length << "\n";
	for (size_t i = 0; i < length; ++i)
	{
		int &indx = bucketsindxs[i];
		if ( indx >= (int)m_buckets.size() || indx < 0	)
		{
			continue;
		}
		std::vector<int>& curr = m_buckets[indx];
		//previous.insert(previous.end(), curr.begin(),curr.end());
		bounded.insert(bounded.end(), curr.begin(),curr.end());
	}

	//std::copy(previous.begin(), previous.end(), std::back_inserter(bounded));

	//std::cout << "pr = " <<(int)previous.size() << "\n";
	//Sleep(50);

	//std::vector<int>::iterator ptr;

	//for (ptr = previous.begin(); ptr != previous.end(); ++ptr)
	//{
	//	const Target& t = m_targets[*ptr];
	//	if ( distance( t.x, t.y, x, y ) < radius )
	//	{
	//		bounded.push_back(*ptr);
	//	}
	//}
}

template<class Target>
void PlanarContainer<Target>::GetInvolvedBuckets( double& x, double& y, double& radius, std::vector<int>& bucketIndxs )
{
	double frameleft, frameright, frametop, framebottom;
	std::vector<int> to_search;
	frameleft = x - radius - m_side;
	frameright = x + radius;
	framebottom = y - radius;
	frametop = y + radius;

	int leftbottom, rightbottom, lefttop, righttop;
	leftbottom = GetBucketIndx(frameleft,framebottom);
	rightbottom = GetBucketIndx(frameright,framebottom);
	lefttop = GetBucketIndx(frameleft,frametop);
	righttop = GetBucketIndx(frameright,frametop);

	for ( int i = leftbottom; i <= righttop; ++i )
	{
		if ( i >= (int)m_buckets.size() || i < 0	)
		{
			continue;
		}
		to_search.push_back(i);
		if( ( (i - rightbottom) % m_BWidth ) == 0 )
		{
			i = leftbottom + m_BWidth*((i - rightbottom) / m_BWidth + 1);
		}
	}

	bucketIndxs.resize(0);
	int length = (int)to_search.size();
	for (int i = 0; i < length; ++i)
	{
		if ( IsBucketOverlap( to_search[i], x, y, radius))
		{
			bucketIndxs.push_back( to_search[i] );
		}
	}
};

template<class Target>
inline bool PlanarContainer<Target>::IsBucketOverlap( int indx, double x, double y, double radius )
{
	double rcframeleft, rcframeright, rcframetop, rcframebottom;
	int xi = indx % m_BWidth, yi = indx / m_BWidth;
	rcframeleft = xi * m_side;
	rcframeright = xi * m_side + m_side;
	rcframebottom = yi * m_side;
	rcframetop = yi * m_side + m_side;

	if ( pointinside( rcframeleft, rcframeright, rcframebottom, rcframetop, x, y))
	{
		return true;
	}

	if (
		( distance( x, y, rcframeleft, rcframebottom) < radius ) || 
		( distance( x, y, rcframeleft, rcframetop) < radius ) || 
		( distance( x, y, rcframeright, rcframebottom) < radius ) || 
		( distance( x, y, rcframeright, rcframetop) < radius )
		)
	{
		return true;
	}

	double sphfrmleft, sphfrmright, sphfrmtop, sphfrmbottom;
	sphfrmleft = x - radius;
	sphfrmright = x + radius;
	sphfrmbottom = y - radius;
	sphfrmtop = y + radius;

	if (
		pointinside( rcframeleft, rcframeright, rcframebottom, rcframetop, sphfrmleft, y) || 
		pointinside( rcframeleft, rcframeright, rcframebottom, rcframetop, sphfrmright, y) || 
		pointinside( rcframeleft, rcframeright, rcframebottom, rcframetop, x, sphfrmbottom) || 
		pointinside( rcframeleft, rcframeright, rcframebottom, rcframetop, x, sphfrmtop)
		)
	{
		return true;
	}

	return false;
}
