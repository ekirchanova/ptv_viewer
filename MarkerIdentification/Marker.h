// Marker.h: interface for the Marker class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MARKER_H__EC97EFDA_6138_49E2_9AD3_A020C00A0B85__INCLUDED_)
#define AFX_MARKER_H__EC97EFDA_6138_49E2_9AD3_A020C00A0B85__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class Marker  
{
	friend bool operator == (const Marker &m, const Marker &m2);

public:
	Marker();
	Marker(const Marker&);
	Marker(float x, float y, float d) : m_x(x), m_y(y), m_d(d), m_level(0), m_X(0.0f), m_Y(0.0f), m_Z(0.0f), m_used(false) {};
	virtual ~Marker();

	void operator = (const Marker &);

	const float get_x() const {return m_x;}
	const float get_y() const {return m_y;}
	void set_x(float x) {m_x = x;}
	void set_y(float y) {m_y = y;}

	const float get_d() const {return m_d;}

	const float get_X() const {return m_X;}
	void set_X(float X) {m_X = X;}

	const float get_Y() const {return m_Y;}
	void set_Y(float Y) {m_Y = Y;}

	const float get_Z() const {return m_Z;}
	void set_Z(float Z) {m_Z = Z;}

	const short get_level() const {return m_level;}
	void set_level(short level) {m_level = level;}

	const bool get_used() const {return m_used;}
	void set_used(bool used) {m_used = used;}

private:
	/**положение маркра по оси x с подпиксельной точностью*/
	float m_x;

	/**положение маркра по оси y с подпиксельной точностью*/
	float m_y;
	
	/**характерный диаметр маркера*/
	float m_d;

	/**уровень плоскости, на которой находится маркер*/
	short m_level;

	/**преобразованная координата положения маркра по оси X с подпиксельной точностью*/
	float m_X;

	/**преобразованная координата положения маркра по оси Y с подпиксельной точностью*/
	float m_Y;

	/**преобразованная координата положения маркра по оси Z с подпиксельной точностью*/
	float m_Z;

	/**флаг означающий, что для маркера была найдена точка в мировой системе координат*/
	bool m_used;

};

bool operator == (const Marker &m, const Marker &m2);

#endif // !defined(AFX_MARKER_H__EC97EFDA_6138_49E2_9AD3_A020C00A0B85__INCLUDED_)