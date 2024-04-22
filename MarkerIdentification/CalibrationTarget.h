// CalibrationTarget.h: interface for the CalibrationTarget class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CALIBRATIONTARGET_H__969189D5_89E1_4F21_9B6C_CBF566DAC0C9__INCLUDED_)
#define AFX_CALIBRATIONTARGET_H__969189D5_89E1_4F21_9B6C_CBF566DAC0C9__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "../../../ProcUtils/serializer/Serializer.h"
#include "../../../ProcUtils/serializer/SerializerLayout.h"
#include "../../../ProcUtils/serializer/SerializerPrim.h"
#include "../../../ProcUtils/serializer/SerializerXML.h"
#include "../../../ProcUtils/serializer/SerializerUtil.h"

namespace Serializer = ProcUtils::Serializer;

class CalibrationTarget
{
public:
	CalibrationTarget();
	CalibrationTarget(tstring name, float dotSpacing, float zeroMarkerD, float axisMarkerD, float mainMarkerD, bool multiLevel, float levelDist, bool bgType);
	virtual ~CalibrationTarget();

	const unsigned long get_hash() const {return m_hash;}
	void set_hash() {m_hash = getHash();}

	const tstring get_name() const {return m_name;}
	void set_name(tstring name) {m_name = name;}

	const float get_dotSpacing() const {return m_dotSpacing;}
	void set_dotSpacing(float dotSpacing) {m_dotSpacing = dotSpacing;}

	const float get_zeroMarkerD() const {return m_zeroMarkerD;}
	void set_zeroMarkerD(float zeroMarkerD) {m_zeroMarkerD = zeroMarkerD;}

	const float get_axisMarkerD() const {return m_axisMarkerD;}
	void set_axizMarkerD(float axisMarkerD) {m_axisMarkerD = axisMarkerD;}

	const float get_mainMarkerD() const {return m_mainMarkerD;}
	void set_mainMarkerD(float mainMarkerD) {m_mainMarkerD = mainMarkerD;}

	const bool get_multiLevel() const {return m_multiLevel;}
	void set_multiLevel(bool multiLevel) {m_multiLevel = multiLevel;}

	const float get_levelDist() const {return m_levelDist;}
	void set_levelDist(float levelDist) {m_levelDist = levelDist;}

	const bool get_bgType() const {return m_bgType;}
	void set_bgType(bool bgType) {m_bgType = bgType;}

	const float get_multiLevelZerroThickness() const {return m_multiLevelZerroThickness;}
	void set_multiLevelZerroThickness(float multiLevelZerroThickness) {m_multiLevelZerroThickness = multiLevelZerroThickness;}

private:
	unsigned long getHash();

public:
	unsigned long m_hash;
	tstring m_name;
	float m_dotSpacing;
	float m_mainMarkerD;
	float m_axisMarkerD;
	float m_zeroMarkerD;
	bool m_multiLevel;
	float m_levelDist;
	bool m_bgType;
	// позже сделать толщину по нулевому уровню многоуровневой двухсторонней мишени считываемой из xml файла!!!
	float m_multiLevelZerroThickness;


public:
	struct LayoutDefault : public Serializer::FieldAttributes::Layout<CalibrationTarget>
	{
		LayoutDefault()
		{
			Simple(L"hash", &CalibrationTarget::m_hash, unsigned long(0));
			Simple(L"name", &CalibrationTarget::m_name, std::wstring(L""));
			Simple(L"dotSpacing", &CalibrationTarget::m_dotSpacing, float(0.0f));
			Simple(L"mainMarkerD", &CalibrationTarget::m_mainMarkerD, float(0.0f));
			Simple(L"axisMarkerD", &CalibrationTarget::m_axisMarkerD, float(0.0f));
			Simple(L"zeroMarkerD", &CalibrationTarget::m_zeroMarkerD, float(0.0f));
			Simple(L"multiLevel", &CalibrationTarget::m_multiLevel, bool(false));
			Simple(L"levelDist", &CalibrationTarget::m_levelDist, float(0.0f));
			Simple(L"bgType", &CalibrationTarget::m_bgType, bool(false));
		}
	};

};

#endif // !defined(AFX_CALIBRATIONTARGET_H__969189D5_89E1_4F21_9B6C_CBF566DAC0C9__INCLUDED_)