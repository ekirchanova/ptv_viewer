#include "stdafx.h"
#include ".\misalignmentcorr.h"
#include <cmath>

/** Производные вычислял в Математике 5.0 (ниже код для вычисления производных)

RX := {{1, 0 , 0}, {0, Cos[rx], Sin[rx]}, {0, -Sin[rx], Cos[rx]}}
RY := {{Cos[ry], 0, Sin[ry]}, {0, 1, 0}, {-Sin[ry], 0, Cos[ry]}}
RZ := {{Cos[rz], Sin[rz], 0}, {-Sin[rz], Cos[rz], 0}, {0, 0, 1}}
tZ := {0, 0, tz}
v := {x, y, z}
M := RX.RY.RZ.v + tZ
MatrixForm[M]

D[M[[1]], x]
D[M[[1]], y]
D[M[[1]], z]
D[M[[2]], x]
D[M[[2]], y]
D[M[[2]], z]
D[M[[3]], x]
D[M[[3]], y]
D[M[[3]], z]

*/
 
int MisalignmentCorr::Init(float rx, float ry, float rz, float tx, float ty, float tz)
{
	//rz = -0.2f;//ZALIPUHA!!!
	//tx = -1.0f;
	//ty = -1.0f;

	if ( rx == 0 && ry == 0 && tz == 0 && rz == 0 && tx == 0 && ty == 0 )
	{
		m_doCalc = false;
	}
	else
	{
		m_doCalc = true;
	}

	m_t[0] = tx;
	m_t[1] = ty;
	m_t[2] = tz;

	// преобразуем из градусов в радианы
	rx = -rx / 180 * PI;// знак минус, чтобы было положительное направление вращения вокруг оси x
	ry = ry / 180 * PI;
	rz = rz / 180 * PI;

	//вокруг X
	m_rx[0][0] = 1.0f;
	m_rx[0][1] = 0.0f;
	m_rx[0][2] = 0.0f;

	m_rx[1][0] = 0.0f;
	m_rx[1][1] = cos(rx);
	m_rx[1][2] = sin(rx);

	m_rx[2][0] = 0.0f;
	m_rx[2][1] = -sin(rx);
	m_rx[2][2] = cos(rx);

	//вокруг Y
	m_ry[0][0] = cos(ry);
	m_ry[0][1] = 0.0f;
	m_ry[0][2] = sin(ry);

	m_ry[1][0] = 0.0f;
	m_ry[1][1] = 1.0f;
	m_ry[1][2] = 0.0f;

	m_ry[2][0] = -sin(ry);
	m_ry[2][1] = 0.0f;
	m_ry[2][2] = cos(ry);

	//вокруг Z
	m_rz[0][0] = cos(rz);
	m_rz[0][1] = sin(rz);
	m_rz[0][2] = 0.0f;

	m_rz[1][0] = -sin(rz);
	m_rz[1][1] = cos(rz);
	m_rz[1][2] = 0.0f;

	m_rz[2][0] = 0.0f;
	m_rz[2][1] = 0.0f;
	m_rz[2][2] = 1.0f;

	//// Расчёт производных преобразования RX*RY*RZ*X+T
	//m_dxx = cos(ry) * cos(rz);
	//m_dxy = cos(ry) * sin(rz);
	//m_dxz = sin(ry);

	//m_dyx = -cos(rz) * sin(rx) * sin(ry) - cos(rx) * sin(rz);
	//m_dyy = cos(rx) * cos(rz) - sin(rx) * sin(ry) * sin(rz);
	//m_dyz = cos(ry) * sin(rx);

	//m_dzx = -cos(rx) * cos(rz) * sin(ry) + sin(rx) * sin(rz);
	//m_dzy = -cos(rz) * sin(rx) - cos(rx) * sin(ry) * sin(rz);
	//m_dzz = cos(rx) * cos(ry);

	return 1;
}

ModelCoord& MisalignmentCorr::operator()(ModelCoord &mc)
{
	// если преобразование "пустое", то пропустить расчёт
	if ( m_doCalc == false ) return mc;

	m_tmp[0] = mc.X;
	m_tmp[1] = mc.Y;
	m_tmp[2] = mc.Z;

	// повернуть и переместить (результат операции сразу записывается в вектор m_tmp)
	m_rx * m_tmp;
	m_ry * m_tmp;
	m_rz * m_tmp;
	m_tmp + m_t;

	mc.X = m_tmp[0];
	mc.Y = m_tmp[1];
	mc.Z = m_tmp[2];
	
	return mc;
}