#pragma once

#include <float.h>
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>

#include "..\Data\ImageMap\IImageMap.h"
#include "..\Data\Coordinates\Coordinates.h"
#include "CalibrationTarget.h"
#include "Marker.h"
#include "MisalignmentCorr.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

typedef std::vector< Marker > MarkerColumn;
typedef std::vector< Marker >::iterator MarkerColumnIt;
typedef std::vector< Marker >::reverse_iterator MarkerColumnRevIt;

/**
���������������� ������:
	- ������������� �������� ������������� �������
	- ������������� ��������������� �������� ������� ��������� (������� �������, ������ �������)
	- ����� ��������� ������� �� ����������� � ������������� ���������
	- ����� ������������� �������������� ������������ ��������
	- ������������� ��������� �������� �� ����������� ��������� ����������� ������������� �������
	- ���������� ����������� ���� ������������� �������� �� ��������
	2005-12-10
	- �������������� ������� �������� �������� � ������������ � z-�������� ������������� ������
*/
class MarkerSearcher
{
public:
	MarkerSearcher() : m_pCT(0), MAX_INTENS_VAL(255), SEGMENTATION_INTENS_VAL(128)
	{
		//////////////////////////////////////////////////////////////////////////
		m_Z = 0.0f;
		m_xOrient = xRight;
		m_yOrient = yUp;
		m_bFlipXY = false;
		m_multilevelBackSideEnable = false;
		m_Q = 5.0f;
		m_minDotD = 1;
		m_maxDotD = 4000;
		m_minDotCount = 25;
		m_OriginAreaThreshold = 1.1f;
		m_boundaryPartX = 0.0f;
		m_boundaryPartY = 0.0f;
		//
		m_grayThrshld = false;
		m_grayValueManual = 128;
		m_grayValueAuto = 0;
		m_grayValueCurrent = 0;
		m_manualCenterMarker = false;
		m_manualCenterPos = Position(500, 500);
		m_windowSz = 100;
		//
		m_templCorrToogle = true;
		m_viewToggle = 0;
		m_bSaveSegmentation = false;
		m_radialDistCoeff = 0;
		m_radialDistCoeffAuto = 0;
		//////////////////////////////////////////////////////////////////////////
		m_id = s_id++;
		//////////////////////////////////////////////////////////////////////////
	}

	MarkerSearcher( 
		CalibrationTarget& CT,
		float Z,
		xImgCoordOrient xOrient,
		yImgCoordOrient yOrient,
		bool flipXYToogle,
		bool multilevelBackSideEnable,
		float Q,
		float minDotD,
		float maxDotD,
		int minDotCount,
		float OriginAreaThreshold,
		float boundaryPartX,
		float boundaryPartY,
		bool grayThrshld,
		int grayValue,
		bool manualCenterMarker,
		Position manualCenterPos,
		int windowSz,
		bool templCorrToogle,
		float rx,
		float ry,
		float tz,
		float ty,
		byte viewToggle,
		bool saveSegmentation,
		float radialDistCoeff = 0
	) : m_pCT(&CT), MAX_INTENS_VAL(255), SEGMENTATION_INTENS_VAL(128), m_mCorr(rx, ry, 0.0f, 0.0f, ty, tz)
	{
		//////////////////////////////////////////////////////////////////////////	
		m_Z = Z;
		m_xOrient = xOrient;
		m_yOrient = yOrient;
		m_bFlipXY = flipXYToogle;
		m_multilevelBackSideEnable = multilevelBackSideEnable;
		m_Q = Q;
		m_minDotD = minDotD;
		m_maxDotD = maxDotD;
		m_minDotCount = minDotCount;
		m_OriginAreaThreshold = OriginAreaThreshold;
		m_boundaryPartX = boundaryPartX;
		m_boundaryPartY = boundaryPartY;
		//
		m_grayThrshld = grayThrshld;
		m_grayValueManual = grayValue;
		m_grayValueAuto = 0;
		m_grayValueCurrent = 0;
		m_manualCenterMarker = manualCenterMarker;
		m_manualCenterPos = manualCenterPos;
		m_windowSz = windowSz;
		//
		m_templCorrToogle = templCorrToogle;
		m_viewToggle = viewToggle;
		m_bSaveSegmentation = saveSegmentation;
		m_radialDistCoeff = radialDistCoeff;
		m_radialDistCoeffAuto = 0;
		//////////////////////////////////////////////////////////////////////////
		m_id = s_id++;
	}

	virtual ~MarkerSearcher()
	{
		--s_id;
	}

	/**����� �������� �� �����������*/
	/**� ���������� ����������� ����� ��������!!!*/
	int Search(IImageMap* pI);

	float get_Z() const { return m_Z; };

	const unsigned short get_centerIndex()
	{
		unsigned short centerIndex = 0;

		MarkerColumnIt centerMarkerIt = std::find(m_markerSet.begin(), m_markerSet.end(), m_OriginMarker);

		if (centerMarkerIt != m_markerSet.end())
		{
			centerIndex = std::distance(m_markerSet.begin(), centerMarkerIt);
		}
		else
		{
			// �� ������ ������ ������������ �������
			//ATLASSERT(0);
			return 0;
		}

		return centerIndex;
	}

	const Marker& get_OriginMarker() const { return m_OriginMarker; }

	std::vector< CoordPair >& get_CoordPairSet() { return m_cpSet; }

	MarkerColumn& get_MarkerSet() { return m_markerSet; }

	const xImgCoordOrient& get_xOrient() const { return m_xOrient; }

	const yImgCoordOrient& get_yOrient() const { return m_yOrient; }

	/**����������� ����� ������������� ������� � ������� ����������� ��� ���� �������������� ��������*/
	static ModelRegion CommonModelRegionQuadrangle(const Quadrangle&, const Quadrangle&);

	static Quadrangle GetQuadrangleRegion(int coordToggle, bool bFlipXY, ModelCoord& c1, ModelCoord& c2, ModelCoord& c3, ModelCoord& c4);

	static void ClearID() { s_id = 0; }

	// == ���������� ��� ��������������
	void SortCoordPairSet();

	void SaveCoordPairSet(tstring);

	const float get_MeanD() const { return m_meanD; }
	const ushort get_GrayValueAuto() const { return m_grayValueAuto; }
	float get_RadialDistCoeffAuto() const { return m_radialDistCoeffAuto; }

private:
	/**�������� �������� ����������� � �����������*/
	int InitTransform(IImageMap* pI);

	/**������������ ����������� � ������� ��������*/
	int MarkerScan(IImageMap* pI);

	/**���������� ������������ ����� ������������ �������� �� ����������� � � ������� ������� ��������� ������������*/
	int BindImageMarkersToModel();

	int BindImageMarkersToModel2();

	/**�������������� ������� �������� �������� � ������������ � z-�������� ������������� ������*/
	int RetrieveColumnLevels();

	/**��������� ��������� ������ ������� ����� �����-���������� � ��� ��������*/
	int TemplateCorrelation(IImageMap* pI, std::vector<Marker>& ms);

	int ShiftToDiscrete(float x);

	int PreModelFit(std::vector< CoordPair >&, bool = false);

	static const short s_stateLevelTable[4];

	std::vector< CoordPair > m_cpSet;
	MarkerColumn m_markerSet;

	CalibrationTarget* m_pCT;
	Position m_Origin;
	Marker m_OriginMarker;
	float m_Z;

	//////////////////////////////////////////////////////////////////////////
	// ���������� ������ ���������� �� �������� ������

	// ���������� ������������� �������� ��� ������� �������
	float m_Q;


	// ����������� ���������� �������� � ������� �������
	float m_minDotD;
	// ������������ ���������� �������� � ������� �������
	float m_maxDotD;
	// ������� ������� �������� �� ����������� ����� ������ m_minDotD � m_maxDotD
	float m_meanD;

	// ����� ����������� (����������) ����������� �������� ����������� (����� �� ������ ����� � ������ �����������)
	// �� x ����������
	float m_boundaryPartX;

	// ����� ����������� (����������) ����������� �������� ����������� (����� �� ������ ����� � ������ �����������)
	// �� y ����������
	float m_boundaryPartY;

	// ����������� ����� ��������, ��� ������� ������ ����� ��������
	int m_minDotCount;

	// ����������� ������ �� ������ 
	float m_OriginAreaThreshold;

	//////////////////////////////////////////////////////////////////////////

	// ���� ������������� ������� ������ �����������
	bool m_grayThrshld;

	// �������� ������������� ��� ������� ������ �����������
	ushort m_grayValueManual;
	// �������� ������������� ��� ������� ������ �����������
	ushort m_grayValueAuto;
	// ������� ������������� ����� �����������
	ushort m_grayValueCurrent;

	// ���� ������������� ������� ������ ������ ������������ �������
	bool m_manualCenterMarker;

	// ��������������� ��������� ������������ �������
	Position m_manualCenterPos;

	// ������ ���� ������ ������������ ������� � px
	int m_windowSz;

	// ���� ������������� ��������� ���������� � ��������
	bool m_templCorrToogle;

	MisalignmentCorr m_mCorr;

	unsigned short m_centerIndex;

	byte m_viewToggle;

	bool m_bSaveSegmentation;

	float m_radialDistCoeff;
	float m_radialDistCoeffAuto;

	//////////////////////////////////////////////////////////////////////////
	// ����������� ���������� ���� ������� ���������

	// ��� x
	xImgCoordOrient m_xOrient;

	// ��� y
	yImgCoordOrient m_yOrient;

	bool m_bFlipXY;
	bool m_multilevelBackSideEnable;

	//////////////////////////////////////////////////////////////////////////

	// ������ � ������ ����� ��� �������
	CString _debug_f_name_;

	// id ����� �������� ������ ������ �� ����� � ���������
	static int s_id;
	int m_id;
	//

	// 2006_11_22 ��������� 16-�� ������ ����������
	int MAX_INTENS_VAL;
	int SEGMENTATION_INTENS_VAL;

};

/**
������� ��� ���������� ������ ���� �������
*/
class positionAverageF
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
�������� ��� ��������� ��������� ��������
*/

inline auto dCompareMarkerF = [](const Marker& m1, const Marker& m2)->bool
{
	return m1.get_d() < m2.get_d();
};
//class dCompareMarkerF : public std::binary_function< Marker, Marker, bool >
//{
//public:
//	dCompareMarkerF() {}
//	bool operator () (const Marker &m1, const Marker &m2) const
//	{
//		return m1.get_d() < m2.get_d();
//	}
//
//};

/**
�������� ���������� �������� � ���������� vector �� �������� ������������ � ��������� �������
*/

//1137 compareByDistToMarkerF //EK

//;class compareByDistToMarkerF : public std::binary_function< Marker, Marker, bool >
//{
//public:
//	compareByDistToMarkerF(const Marker &m) : m_refMarker(m) {}
//
//	bool operator () (const Marker &m, const Marker &m2)
//	{
//		return _hypot(m.get_x() - m_refMarker.get_x(), m.get_y() - m_refMarker.get_y())
//			 < _hypot(m2.get_x() - m_refMarker.get_x(), m2.get_y() - m_refMarker.get_y());
//	}
//
//private:
//	const Marker &m_refMarker
//
//};

/**
�������� ��� ���������� �������� �� ��������������� ��������� (� ������� ������� ���������) �� �����������
*/

inline auto directCompareMarkerXWCF = [](const Marker& m1, const Marker& m2)->bool
{
	return m1.get_X() < m2.get_X();
};

//class directCompareMarkerXWCF : public std::binary_function< Marker, Marker, bool >
//{
//public:
//	directCompareMarkerXWCF() {}
//
//	bool operator () (const Marker &m1, const Marker &m2) const
//	{
//		return m1.get_X() < m2.get_X();
//	}
//
//};

/**
�������� ��� ���������� �������� �� ��������������� ��������� (� ������� ������� ���������) �� ��������
*/

inline auto reverseCompareMarkerXWCF = [](const Marker& m1, const Marker& m2)->bool
{
	return m1.get_X() > m2.get_X();
};
//class reverseCompareMarkerXWCF : public std::binary_function< Marker, Marker, bool >
//	{
//	public:
//		reverseCompareMarkerXWCF() {}
//
//	bool operator () (const Marker &m1, const Marker &m2) const
//		{
//		return m1.get_X() > m2.get_X();
//		}
//
//	};

/**
�������� ��� ���������� �������� �� ��������������� ���������
*/

inline auto compareMarkerXF = [](const Marker& m1, const Marker& m2)-> bool
{
	return m1.get_x() < m2.get_x();
};
//class compareMarkerXF : public std::binary_function< Marker, Marker, bool >
//{
//public:
//	compareMarkerXF() {}
//
//	bool operator () (const Marker &m1, const Marker &m2) const
//	{
//		return m1.get_x() < m2.get_x();
//	}
//
//};

/**
�������� ��� ���������� ��� ��������� ������� �� ����������� � � ������� �.�. �� ����������� � ������� �.�.
*/

inline auto compareCoordPairWCSF = [](const CoordPair& cp1, const CoordPair& cp2)->bool
{
	if (cp1.m.Z != cp2.m.Z) return cp1.m.Z < cp2.m.Z;
	if (cp1.m.Y != cp2.m.Y) return cp1.m.Y < cp2.m.Y;
	if (cp1.m.X != cp2.m.X) return cp1.m.X < cp2.m.X;
	return false;
};
//class compareCoordPairWCSF : public std::binary_function< CoordPair, CoordPair, bool >
//{
//public:
//	compareCoordPairWCSF() {}
//
//	bool operator () (const CoordPair &cp1, const CoordPair &cp2) const
//	{
//		if ( cp1.m.Z != cp2.m.Z ) return cp1.m.Z < cp2.m.Z;
//		if ( cp1.m.Y != cp2.m.Y ) return cp1.m.Y < cp2.m.Y;
//		if ( cp1.m.X != cp2.m.X ) return cp1.m.X < cp2.m.X;
//		return false;
//	}
//
//};

/**
�������� ��� ���������� ��� ��������� ������� �� ����������� � � ������� �.�. �� ����������� � ������� �.�.
*/

inline auto compareMarkerWCSF = [](const Marker& m1, const Marker& m2)->bool
{
	float X1 = m1.get_X();
	float X2 = m2.get_X();
	float Y1 = m1.get_Y();
	float Y2 = m2.get_Y();
	float Z1 = m1.get_Z();
	float Z2 = m2.get_Z();
	if (Z1 != Z2) return Z1 < Z2;
	if (Y1 != Y2) return Y1 < Y2;
	if (X1 != X2) return X1 < X2;
	return false;
};
//class compareMarkerWCSF : public std::binary_function< Marker, Marker, bool >
//{
//public:
//	compareMarkerWCSF() {}
//
//	bool operator () (const Marker &m1, const Marker &m2) const
//	{
//		float X1 = m1.get_X();
//		float X2 = m2.get_X();
//		float Y1 = m1.get_Y();
//		float Y2 = m2.get_Y();
//		float Z1 = m1.get_Z();
//		float Z2 = m2.get_Z();
//		if ( Z1 != Z2 ) return Z1 < Z2;
//		if ( Y1 != Y2 ) return Y1 < Y2;
//		if ( X1 != X2 ) return X1 < X2;
//		return false;
//	}
//
//};

/**
�������� ��� ���������� �������� �� ������������� ���������
*/
inline auto compareMarkerYF = [](const Marker& m1, const Marker& m2)->bool
{
	return m1.get_y() < m2.get_y();
};
//class compareMarkerYF : public std::binary_function< Marker, Marker, bool >
//{
//public:
//	compareMarkerYF() {}
//
//	bool operator () (const Marker &m1, const Marker &m2) const
//	{
//		return m1.get_y() < m2.get_y();
//	}
//
//};

/**
������� ��� ������� ������� �������� ������� �������
*/

inline auto sumDiamsF = [](float D, const Marker& m)
{
	return D + m.get_d();
};

//class sumDiamsF : public std::binary_function< float, Marker, float >
//{
//public:
//	sumDiamsF() {}
//
//	float operator () (float D, const Marker& m)
//	{		
//		return D + m.get_d();
//	}
//
//};

/**
�������� ��� �������� �������� �� �������� � ������������� �������� �� �������
*/
class diameterValidMarkerF
{
public:
	diameterValidMarkerF(float meanD, float Q)
	{
		m_minD = meanD / Q;
		m_maxD = meanD * Q;
	}

	bool operator () (const Marker& m) const { return (m.get_d() < m_minD) || (m.get_d() > m_maxD); }

private:
	float m_minD;
	float m_maxD;

};

/**
�������� ��� �������� �������� ��������� �� ������� �������� �������
*/
class boundaryValidMarkerF
{
public:
	boundaryValidMarkerF(float boundaryPartX, float boundaryPartY, unsigned int w, unsigned int h) : m_boundaryPartX(boundaryPartX), m_boundaryPartY(boundaryPartY), m_w(w), m_h(h) {}

	bool operator () (const Marker& m) const
	{
		float x = m.get_x();
		float y = m.get_y();
		float d = m.get_d();

		return 	(x - d < m_w* m_boundaryPartX) ||
			(y - d < m_h* m_boundaryPartY) ||
			(x + d > m_w * (1 - m_boundaryPartX)) ||
			(y + d > m_h * (1 - m_boundaryPartY));
	}

private:
	float m_boundaryPartX;
	float m_boundaryPartY;
	unsigned int m_w;
	unsigned int m_h;

};

/**
�������� ��� �������� �������� ��������� �� ������� ���� ��� ������ ������������ �������
*/
class markersOutOfCenterWindowF
{
public:
	markersOutOfCenterWindowF(Position centerPosition, Position centerWindow)
	{
		int _w = centerWindow.m_i / 2;
		int _h = centerWindow.m_j / 2;
		m_xMin = centerPosition.m_i - _w;
		m_yMin = centerPosition.m_j - _h;
		m_xMax = centerPosition.m_i + _w;
		m_yMax = centerPosition.m_j + _h;
		//char buf[256];
		//sprintf(buf, "%d %d %d %d", m_xMin, m_yMin, m_xMax, m_yMax);
		//MessageBox(0, buf, "", 0);
	}

	bool operator () (const Marker& m) const
	{
		float x = m.get_x();
		float y = m.get_y();
		return x < m_xMin || y < m_yMin || x > m_xMax || y > m_yMax;
	}

private:
	int m_xMin;
	int m_yMin;
	int m_xMax;
	int m_yMax;

};