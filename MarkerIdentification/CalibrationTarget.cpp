// CalibrationTarget.cpp: implementation of the CalibrationTarget class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "CalibrationTarget.h"

#ifdef _DEBUG
#undef THIS_FILE
static TCHAR THIS_FILE[]=_T(__FILE__);
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CalibrationTarget::CalibrationTarget()
{
	/**Инициализация параметрами калибровочной матрицы по-умолчанию*/
	m_name = _T("");
	m_dotSpacing = 10.0f;//2.5f;//10.0f;
	m_zeroMarkerD = 2.0f;
	m_axisMarkerD = 1.0f;
	m_mainMarkerD = 1.5f,
	m_multiLevel = false;
	m_levelDist = 0.0f;
	m_bgType = false;
	
	m_hash = getHash();
}

CalibrationTarget::CalibrationTarget(tstring name, float dotSpacing, float zeroMarkerD, float axisMarkerD, float mainMarkerD, bool multiLevel, float levelDist, bool bgType)
{
	m_name = name;
	m_dotSpacing = dotSpacing;
	m_zeroMarkerD = zeroMarkerD;
	m_axisMarkerD = axisMarkerD;
	m_mainMarkerD = mainMarkerD;
	m_multiLevel = multiLevel;
	m_levelDist = levelDist;
	m_bgType = bgType;
	
	m_hash = getHash();
}


CalibrationTarget::~CalibrationTarget()
{
}

unsigned long CalibrationTarget::getHash()
{
	unsigned int multiLevelByte = m_multiLevel ? (byte)2 : (byte)0;
	unsigned int bgTypeByte = m_bgType ? (byte)1 : (byte)0;

	return ( (byte)(m_dotSpacing  * 10) ) << 24 |
		   ( (byte)(m_zeroMarkerD * 10) ) << 16 |
		   ( (byte)(m_axisMarkerD * 10) ) <<  8 |
		   ( (byte)(m_mainMarkerD * 10) )       |
		   ( (byte)(m_levelDist   * 10) ) <<  4 |
		           multiLevelByte             |
		           bgTypeByte;
}