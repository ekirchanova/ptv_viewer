#pragma once

#include "..\Data\Coordinates\Coordinates.h"
#include "..\Math\tnt_cmat.h"

class MisalignmentCorr
{
public:
	MisalignmentCorr() : m_rx(3, 3), m_ry(3, 3), m_rz(3, 3), m_t(3), m_tmp(3), m_doCalc(false)
	{		 
		TNT::Vector<float>::AllocExtraBuf(3);
		Init(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
	};

	MisalignmentCorr(float rx, float ry, float rz, float tx, float ty, float tz) : m_rx(3, 3), m_ry(3, 3), m_rz(3, 3), m_t(3), m_tmp(3), m_doCalc(false)
	{		
		TNT::Vector<float>::AllocExtraBuf(3);
		Init(rx, ry, rz, tx, ty, tz);
	};

	~MisalignmentCorr(void)
	{
		TNT::Vector<float>::DeallocateExtraBuf();
	};
	
	ModelCoord& operator()(ModelCoord &mc);

	// �����������
	//const float dxx() const {return m_dxx;}
	//const float dxy() const {return m_dxy;}
	//const float dxz() const {return m_dxz;}

	//const float dyx() const {return m_dyx;}
	//const float dyy() const {return m_dyy;}
	//const float dyz() const {return m_dyz;}

	//const float dzx() const {return m_dzx;}
	//const float dzy() const {return m_dzy;}
	//const float dzz() const {return m_dzz;}

private:
	/**������������� ������ ��������������*/
	int Init(float rx, float ry, float rz, float tx, float ty, float tz);

private:
	/**������ ��������*/
	TNT::Vector<float> m_t;

	/**������� �������� Rx*/
	TNT::Matrix<float> m_rx;

	/**������� �������� Ry*/
	TNT::Matrix<float> m_ry;

	/**������� �������� Rz*/
	TNT::Matrix<float> m_rz;

	/**��������� �������*/
	TNT::Vector<float> m_tmp;

	/**����������� ����� �� ������������ ����� ��� ������� ���������*/
	bool m_doCalc;

	///**�����������*/
	//float m_dxx;
	//float m_dxy;
	//float m_dxz;

	//float m_dyx;
	//float m_dyy;
	//float m_dyz;

	//float m_dzx;
	//float m_dzy;
	//float m_dzz;

};