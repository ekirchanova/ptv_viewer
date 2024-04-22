// Marker.cpp: implementation of the Marker class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Marker.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Marker::Marker() : m_x(0.0f), m_y(0.0f), m_d(0.0f), m_level(0), m_X(0.0f), m_Y(0.0f), m_Z(0.0f), m_used(false)
{	
}

Marker::~Marker()
{
}

Marker::Marker(const Marker& m) : m_x(m.m_x), m_y(m.m_y), m_d(m.m_d), m_level(m.m_level), m_X(m.m_X), m_Y(m.m_Y), m_Z(m.m_Z), m_used(m.m_used)
{	
}

void Marker::operator = (const Marker& m)
{
	m_x = m.m_x;
	m_y = m.m_y;
	m_d = m.m_d;
	m_level = m.m_level;
	m_X = m.m_X;
	m_Y = m.m_Y;
	m_Z = m.m_Z;
	m_used = m.m_used;
}

bool operator == (const Marker &m, const Marker &m2)
{
	return ( m.m_x == m2.m_x ) && ( m.m_y == m2.m_y );
}