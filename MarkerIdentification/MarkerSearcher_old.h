// MarkerSearcher.h: interface for the MarkerSearcher class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MARKERSEARCHER_H__22E323BC_8FEC_4551_B633_1F13962EE173__INCLUDED_)
#define AFX_MARKERSEARCHER_H__22E323BC_8FEC_4551_B633_1F13962EE173__INCLUDED_

#include <vector>
#include <functional>
#include <algorithm>
#include "..\Data\ImageMap\IImageMap.h"
#include "..\Data\Coordinates\Coordinates.h"
#include "CalibrationTarget.h"
#include "Marker.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

typedef std::vector< std::vector< Marker > > MarkerField;
typedef std::vector< std::vector< Marker > >::iterator MarkerFieldIt;
typedef std::vector< Marker > MarkerColumn;
typedef std::vector< Marker >::iterator MarkerColumnIt;

/**
���������������� ������:
	- ������������� �������� ������������� �������
	- ������������� ��������������� �������� ������� ��������� (������� �������, ������ �������)
	- ����� ��������� ������� �� ����������� � ������������� ���������
	- ����� ������������� �������������� ������������ ��������
	- ������������� ��������� �������� �� ����������� ��������� ����������� ������������� �������
	- ���������� ����������� ���� ������������� �������� �� ��������
*/
class MarkerSearcher
{
public:
	MarkerSearcher() : m_pCT(0)
	{
		m_Z = 0.0f;
		m_xOrient = xRight;
		m_yOrient = yUp;
		m_Q = 5.0f;
		m_minDotArea = 100;
		m_maxDotArea = 10000;
		m_minDotCount = 25;
		m_OriginAreaThreshold = 1.1f;
		m_boundaryPart = 0.0f;
		
		//////////////////////////////////////////////////////////////////////////
		m_id = s_id++;
	};
	MarkerSearcher::MarkerSearcher(
		CalibrationTarget &CT,
		float Z,
		xImgCoordOrient xOrient,
		yImgCoordOrient yOrient,
		float Q,
		int minDotArea,
		int maxDotArea,
		int minDotCount,
		float OriginAreaThreshold,
		float boundaryPart);
	virtual ~MarkerSearcher();
	
	/**����� �������� �� �����������*/
	/**� ���������� ����������� ����� ��������!!!*/
	int Search(IImageMap *pI);
	
	MarkerField& get_MarkerSet() {return m_markerSet;};

	const Position& get_Origin() const {return m_Origin;};

	std::vector< CoordPair >& get_CoordPairSet() {return m_cpSet;}

	const ModelRegion& get_ModelRegion() const {return m_r;}	

	const xImgCoordOrient& get_xOrient() const {return m_xOrient;}

	const yImgCoordOrient& get_yOrient() const {return m_yOrient;}

	/**������ ����������� ���� �������� ��������*/
	/**����� ����������� ��������� ����������!!!*/
	/**�������� ������ ��� ������������� ��������!!!*/
	static int Intersection(MarkerSearcher &leftMS, MarkerSearcher &rightMS);

	/**����������� ����� ������� � ��������� ����������� ��� ���� �������� ��������*/
	static ModelRegion CommonModelRegion(const ModelRegion &, const ModelRegion &);

private:	
	/**�������� �������� ����������� � �����������*/
	int InitTransform(IImageMap *pI);
		
	/**������������ ����������� � ������� ��������*/
	int MarkerScan(IImageMap *pI);

	/**��������� ������������� �������*/
	int ClipRectangle();

	/**���������� ������������ ����� ������������ �������� �� ����������� � � ������*/
	int BindImageMarkersToModel();
	
	MarkerField m_markerSet;
	std::vector< CoordPair > m_cpSet;

	CalibrationTarget *m_pCT;
	Position m_Origin;
	Marker m_OriginMarker;
	float m_Z;

	/**������ ������������� (��� ����� �.�. ���� Z-����������), ����������� ������� �� �������� � ��������� ������� ���������*/
	ModelRegion m_r;	

	//////////////////////////////////////////////////////////////////////////
	// ���������� ������ ���������� �� �������� ������

	// ���������� ������������� �������� ��� ������� �������
	float m_Q;

	// ����������� ���������� �������� � ������� �������
	int m_minDotArea;

	// ������������ ���������� �������� � ������� �������
	int m_maxDotArea;

	// ����� ����������� (����������) ����������� �������� ����������� (����� �� ������ ����� � ������ �����������)
	float m_boundaryPart;

	// ����������� ����� ��������, ��� ������� ������ ����� ��������
	int m_minDotCount;

	// ����������� �� �������� ������� ����������� �������
	float m_OriginAreaThreshold;

	//////////////////////////////////////////////////////////////////////////
	// ����������� ���������� ���� ������� ���������

	// ��� x
	xImgCoordOrient m_xOrient;

	// ��� y
	yImgCoordOrient m_yOrient;

	//////////////////////////////////////////////////////////////////////////
	
	// ������ � ������ ����� ��� �������
	CString _debug_f_name_;

	// id ����� �������� ������ ������ �� ����� � ���������
	static int s_id;
	int m_id;
	//

};

/**
������� ��� ���������� ������ ���� �������
*/
class positionAverageF : public std::unary_function< Position, void >
{
public:
	positionAverageF() : iSum(0), jSum(0), numPositions(0) {}

	void operator () (const Position& p)
	{
		++numPositions;
		iSum += p.m_i;
		jSum += p.m_j;
	}

	ImageCoord result() const
	{
		return ImageCoord((float)iSum / numPositions, (float)jSum / numPositions);
	}
private:
	size_t numPositions;
	int iSum;
	int jSum;
};

/**
������� ��� ��������� ��������� ��������
*/
class dCompareMarkerF : public std::binary_function< Marker, Marker, bool >
{
public:
	dCompareMarkerF() {}
	bool operator () (const Marker &m1, const Marker &m2) const
	{
		return m1.get_d() < m2.get_d();
	}

};

/**
������� ��� ������ ������� ������������� �������� � ���� �� ��������
*/
class findMaxDMarkerF : public std::unary_function< MarkerColumn, void >
{
public:
	findMaxDMarkerF() : m_maxD(0.0), m_i(0) {}

	findMaxDMarkerF(const findMaxDMarkerF& _F) : m_Origin(_F.m_Origin), m_maxD(_F.m_maxD), m_i(_F.m_i), m_OriginMarker(_F.m_OriginMarker) {};

	void operator = (const findMaxDMarkerF& _F)
	{
		m_Origin = _F.m_Origin;
		m_maxD = _F.m_maxD;
		m_i = _F.m_i;
		m_OriginMarker = _F.m_OriginMarker;
	};

	void operator () (const MarkerColumn& mc)
	{		
		std::vector< Marker >::const_iterator mc_it = std::max_element(mc.begin(), mc.end(), dCompareMarkerF());
		const Marker &m = *mc_it;

		if ( m.get_d() > m_maxD ) {

			m_maxD = m.get_d();
			m_Origin.m_i = m_i;
			m_Origin.m_j = std::distance(mc.begin(), mc_it);
			m_OriginMarker = m;

		}

		m_i++;
	}

	Position get_Origin() const	{return m_Origin;}

	Marker get_OriginMarker() const {return m_OriginMarker;}

private:
	Position m_Origin;
	float m_maxD;
	int m_i;
	Marker m_OriginMarker;

};

/**
������� ��� ��������� � ���������� ������� � ������� ���������
*/
class yLessMarkerF : public std::unary_function< Marker, bool >
{
public:
	yLessMarkerF(float yThreshold) : m_yThreshold(yThreshold) {}
	bool operator () (const Marker &m) const
	{
		return m.get_y() < m_yThreshold;
	}

private:
	float m_yThreshold;

};

/**
������� ��� ��������� � ���������� ������� � ������ ���������
*/
class yGreaterMarkerF : public std::unary_function< Marker, bool >
{
public:
	yGreaterMarkerF(float yThreshold) : m_yThreshold(yThreshold) {}
	bool operator () (const Marker &m) const
	{
		return m.get_y() > m_yThreshold;
	}

private:
	float m_yThreshold;

};

/**
������� ��������������� �������� �� ������������� ���������
*/
class yEqualMarkerF : public std::binary_function< Marker, Marker, bool >
{
public:
	yEqualMarkerF() {}
	bool operator () (const Marker &m1, const Marker &m2) const
	{
		return ( m1.get_y() > ( m2.get_y() - m2.get_d() ) ) && ( m1.get_y() < ( m2.get_y() + m2.get_d() ) );
	}
};

/**
������� ��������������� �������� �� ��������� _Pr
*/
template<typename _Pr>
class sizeEqualVertMarkerLineF : public std::binary_function< const MarkerColumn&, const MarkerColumn&, bool >
{
public:
	sizeEqualVertMarkerLineF(const _Pr& _F) : m_F(_F) {}
	bool operator () (const MarkerColumn& mc1, const MarkerColumn& mc2)
	{
		return std::count_if(mc1.begin(), mc1.end(), m_F) < std::count_if(mc2.begin(), mc2.end(), m_F);
	}

private:
	const _Pr &m_F;
	
};

/**
������� ��� ��������� ��������������� ���� ������������ ����� �������� �� ���������� ���������
*/
class compareSizeVertMarkerLineF : public std::binary_function< MarkerColumn, MarkerColumn, bool >
{
public:
	compareSizeVertMarkerLineF() {}

	bool operator () (const MarkerColumn& mc1, const MarkerColumn& mc2) const
	{
		// ATLTRACE("max (%d, %d)\n", mc1.size(), mc2.size());
		return mc1.size() < mc2.size();
	}

};

/**
������� ��� ��������� ���������� ���� ������������ ����� �������� �� �������
*/
class equalSizeVertLineF : public std::binary_function< MarkerColumn, MarkerColumn, bool >
{
public:
	equalSizeVertLineF() {}

	bool operator () (const MarkerColumn& mc1, const MarkerColumn& mc2) const
	{
		// ATLTRACE("equal (%d, %d)\n", mc1.size(), mc2.size());
		return mc1.size() == mc2.size();
	}

};

/**
������� ��� �������� �������� � ������������ �������, ������� ������� �� �������������� ������
*/
class rangeDeleteMarkersF : public std::unary_function< Marker, bool >
{
public:
	rangeDeleteMarkersF(float k_down, float b_down, float k_up, float b_up) : m_k_down(k_down), m_b_down(b_down), m_k_up(k_up), m_b_up(b_up) {}

	bool operator () (const Marker&m) const
	{
		return ( m.get_y() < ( m_k_down * m.get_x() + m_b_down - m.get_d() ) ) || ( m.get_y() > ( m_k_up * m.get_x() + m_b_up + m.get_d() ) );
	}

private:
	float m_k_down;
	float m_b_down;
	float m_k_up;
	float m_b_up;

};

/**
������� ��� �������� �������� � ������������ ��������, ������� ������������ ������ ��������� _Pr
*/
template<typename _Pr>
class condDeleteMarkersF : public std::unary_function< MarkerColumn, void >
{
public:
	condDeleteMarkersF(const _Pr &_F) : m_F(_F) {}

	void operator () (MarkerColumn& mc)
	{		
		mc.erase(std::remove_if(mc.begin(), mc.end(), m_F), mc.end());
	}

private:
	const _Pr &m_F;

};

/**
������� ��� ���������� �������� ������������� ��������� ����� �������
*/
class compareMarkerPositionF : public std::binary_function< Marker, Marker, bool >
{
public:
	compareMarkerPositionF(int h) : m_h(h) {}

	bool operator () (const Marker &m1, const Marker &m2) const
	{
		return ( m1.get_x() * m_h + m1.get_y() ) < ( m2.get_x() * m_h + m2.get_y() );
	}

private:
	int m_h;

};

/**
������� ��� ���������� �������� �� ������������� ���������
*/
class compareMarkerYF : public std::binary_function< Marker, Marker, bool >
{
public:
	compareMarkerYF() {}

	bool operator () (const Marker &m1, const Marker &m2) const
	{
		return m1.get_y() < m2.get_y();
	}

};

/**
������� ��� ������� ������� �������� ������� �������
*/
class sumAreaF : public std::binary_function< float, Marker, float >
{
public:
	sumAreaF() {}

	float operator () (float area, const Marker& m)
	{		
		return area + 3.14 * m.get_d() * m.get_d() / 4;
	}

};

/**
������� ��� �������� �������� �� �������� � ������������� �������� �� �������
*/
class areaValidMarkerF : public std::unary_function< Marker, bool >
{
public:
	areaValidMarkerF(float meanArea, float Q) : m_meanArea(meanArea), m_Q(Q) {}

	bool operator () (const Marker &m) const
	{
		float area = 3.14 * m.get_d() * m.get_d() / 4;
		float minArea = m_meanArea / ( 1 + m_Q );
		float maxArea = m_meanArea * ( 1 + m_Q );
		//ATLTRACE("%g < %g < %g\n", minArea, area, maxArea);
		return ( area < minArea ) || ( area > maxArea );
	}

private:
	float m_meanArea;
	float m_Q;

};

#endif // !defined(AFX_MARKERSEARCHER_H__22E323BC_8FEC_4551_B633_1F13962EE173__INCLUDED_)