// MarkerSearcher.cpp: implementation of the MarkerSearcher class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MarkerSearcher.h"

#include <ctime>
#include <numeric>
#include <map>
#include <climits>
#include <algorithm>
#include <fstream>

#include "..\ImagingModelFit\DltRevModelFit.h"
#include "preciseTargPointsLoader.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//#define _DEBUG_
#ifdef _DEBUG_
#include <fstream.h>
#endif

//#define READ_POINTS_FROM_FILE
#ifdef READ_POINTS_FROM_FILE
#include <fstream>
#endif

//18.01.17
template<typename T>
T inline sqr(T x)
{
	return x*x;
}
//18.01.17

// Эмпирические константы:

#define PRE_FIT2_DOMAIN_SIZE 1 // размер области из маркеров 11х11=(5+1+5)x(5+1+5) (центральный маркер в центре) для построения второй предварительной модели
							   // первая предварительная модель рассчитывается по кресту из маркеров (см. MarkerSearcher::MarkerScan)

#define MIN_DISTANCE_BETWEEN_ADJACENT_COLUMNS 0.4f // минимальное расстояние между столбцами из маркеров для различения
												   // столбцов с различными уровнями по Z для 3-х уровневой мишени (см. MarkerSearcher::RetrieveColumnLevels)

int MarkerSearcher::ShiftToDiscrete(float x)
{
	//здесь порог округления до целого равен 0.5
	//т.е. если дробная часть больше 0.5 берётся целая часть на 1 больше

	int ix = static_cast<int>(x);
	float dx = x - ix;
	if ( dx > 0.5 ) ++ix;
	else if ( dx < -0.5 ) --ix;

	return ix;
} 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

int MarkerSearcher::s_id = 0;

int MarkerSearcher::Search(IImageMap *pI)
{
#ifdef READ_POINTS_FROM_FILE

	std::string fname("C:\\work\\points.dat");
	std::string mes = "Reading marker points from " + fname;
	MessageBox(0, mes.c_str(), "3DLib.dll", 0);
	std::ifstream in; in.open(fname);
	int mySize;
	in >> mySize;
	std::stringstream sstr; sstr << "Input file " << fname << " has " << mySize << " points";
	MessageBox(0, sstr.str().c_str(), "3DLib.dll", 0);

	m_cpSet.clear();
	m_cpSet.resize(mySize);

	float X, Y, Z;
	ModelCoord mc;
	
	for (int i = 0; i < mySize; ++i)
	{
		// для вычисления модели
		in >> m_cpSet[i].i.x;
		in >> m_cpSet[i].i.y;
		in >> X;
		in >> Y;
		in >> Z;

		//std::stringstream sstr;sstr << my[i].i.x << " " << my[i].i.y << " " << my[i].m.X << " " << my[i].m.Y << " " << my[i].m.Z;
		//::MessageBox(0, sstr.str().c_str(), "", 0);

		mc = m_mCorr(ModelCoord(X, Y, Z));

		m_cpSet[i].m.X = mc.X;
		m_cpSet[i].m.Y = mc.Y;
		m_cpSet[i].m.Z = mc.Z;

		// для попадания маркеров в итоговую таблицу
		Marker m(m_cpSet[i].i.x, m_cpSet[i].i.y, 20);

		m.set_X(mc.X);
		m.set_Y(mc.Y);
		m.set_Z(mc.Z);
		m.set_used(true);
		m_markerSet.push_back(m);
	}

	return 1;
#endif // READ_POINTS_FROM_FILE

	MAX_INTENS_VAL = 255;
	SEGMENTATION_INTENS_VAL = ( MAX_INTENS_VAL + 1 ) / 2; // в функции InitTransform(pI) ниже значение будет переопределено на текущий порог бинаризации!!!

	clock_t dt = clock();

	int result = InitTransform(pI);
	// 1. найденно меньше двух максимумов на гистограмме интенсивности
	//	-
	//ATLASSERT(result);
	if (!result) return 0;

	ATLTRACE(_T("InitTransformTime = %d\n"), clock() - dt);
	dt = clock();

	result = MarkerScan(pI);
	//ATLASSERT(result);
	if (!result) return 0;

	ATLTRACE(_T("MarkerScanTime = %d\n"), clock() - dt);
	dt = clock();

	if ( m_pCT->get_multiLevel() )
	{
		result = RetrieveColumnLevels();
		//ATLASSERT(result);
		if (!result) return 0;

		ATLTRACE(_T("RetrieveColumnLevels = %d\n"), clock() - dt);
		dt = clock();
	}

	//result = ClipRectangle();
	//ATLASSERT(result);
	//if (!result) return 0;

	//ATLTRACE("ClipRectangleTime = %d\n", clock() - dt);
	//dt = clock();

	result = BindImageMarkersToModel();
	//ATLASSERT(result);
	if (!result) return 0;

	ATLTRACE(_T("BindImageMarkersToModelTime = %d\n"), clock() - dt);

	return 1;
}

const short MarkerSearcher::s_stateLevelTable[4] = {0, -1, 0, 1};

/**Восстановление уровней столбцов маркеров в соответствии с z-уровнями калибровочной мишени*/
int MarkerSearcher::RetrieveColumnLevels()
{
	MarkerColumnIt it, it_end;

	short state, level;
	float X, Y, X_last, Y_last;
	int _X = 0, _Y = 0;
	bool firstTime;

	//////////////////////////////////////////////////////////////////////////

	state = 3;
	level = s_stateLevelTable[state];

	X = 0.f;
	Y = 0.f;
	X_last = 0.f;
	Y_last = 0.f;
	firstTime = true;

	//FILE *out;
	//out = _tfopen(_T("D:\\debug2.txt"), _T("wb"));
	//unsigned int counter = 0;

	std::sort(m_markerSet.begin(), m_markerSet.end(), directCompareMarkerXWCF);

	it = m_markerSet.begin();
	it_end = m_markerSet.end();
	
	//int __X = 0;

	// рассмотрим X >= 0 (с учетом округления)
	for (; it != it_end; ++it)
	{
		if ( ShiftToDiscrete(it->get_X()) < 0 ) continue;
		
		X = it->get_X();
		Y = it->get_Y();
		
		if ( firstTime )
		{
			X_last = X;
			Y_last = Y;
			firstTime = false;
		}

		//ATLTRACE(_T("1: fabs(X_last - X) = %g, X = %g\tfabs(Y_last - Y) = %g, Y = %g\n"), fabs(X_last - X), X, fabs(Y_last - Y), Y);
		
		if ( X_last - X < -MIN_DISTANCE_BETWEEN_ADJACENT_COLUMNS )
		{
			++state;
			//++__X;
			if ( state == 4 ) state = 0;
			level = s_stateLevelTable[state];
		}

		it->set_level(level);
		
		//////////////////////////////////////////////////////////////////////////
		// Устанавливаем уровень 0 для трех маркеров в центральном столбце

		_X = ShiftToDiscrete(X);
		_Y = ShiftToDiscrete(Y);

		if ( _X == 0.f && ( _Y == -1.f || _Y == 0.f || _Y == 1.f ) )
		{
			it->set_level(0);
		}

		//////////////////////////////////////////////////////////////////////////

		//_ftprintf(out, _T("%u\t%g\t%g\t%g\t%g\t%d\n"), counter++, X, Y, X_last - X, Y_last - Y, level);
		
		X_last = X;
		Y_last = Y;
		
		//it->set_X(__X);
	}

	//fclose(out);
	//out = _tfopen(_T("D:\\debug3.txt"), _T("wb"));
	//counter = 0;

	//////////////////////////////////////////////////////////////////////////

	state = 0;
	level = s_stateLevelTable[state];

	X = 0.f;
	Y = 0.f;
	X_last = 0.f;
	Y_last = 0.f;
	firstTime = true;

	std::sort(m_markerSet.begin(), m_markerSet.end(), reverseCompareMarkerXWCF);

	it = m_markerSet.begin();
	it_end = m_markerSet.end();
	
	//__X = -1;

	// рассмотрим X < 0 (с учетом округления)
	for (; it != it_end; ++it)
	{
		if ( ShiftToDiscrete(it->get_X()) >= 0 ) continue;

		X = it->get_X();
		Y = it->get_Y();

		if ( firstTime )
		{
			X_last = X;
			Y_last = Y;
			firstTime = false;
		}

		//ATLTRACE(_T("2: fabs(X_last - X) = %g, X = %g\tfabs(Y_last - Y) = %g, Y = %g\n"), fabs(X_last - X), X, fabs(Y_last - Y), Y);			

		if ( X_last - X > MIN_DISTANCE_BETWEEN_ADJACENT_COLUMNS )
		{
			++state;
			//--__X;
			if ( state == 4 ) state = 0;
			level = s_stateLevelTable[state];
		}

		it->set_level(level);

		//_ftprintf(out, _T("%u\t%g\t%g\t%g\t%g\t%d\n"), counter++, X, Y, X_last - X, Y_last - Y, level);
		
		X_last = X;
		Y_last = Y;
		
		//it->set_X(__X);
	}

	//fclose(out);

	return 1;
}

int MarkerSearcher::InitTransform(IImageMap *pI)
{
	int i, j;
	int w = pI->GetWidth();
	int h = pI->GetHeight();

	ushort* pFFrame = pI->GetCurrentFramePtr();

//{
//
//	int k;
//	std::vector<int> intenceGraph;
//	intenceGraph.reserve(MAX_INTENS_VAL + 1);
//	for (k = 0; k <= MAX_INTENS_VAL; ++k) intenceGraph[k] = 0;
//	
//	for (j = 0; j < h; j++)
//		for (i = 0; i < w; i++)
//			intenceGraph[pFFrame[j * w + i]]++;
//
//	//////////////////////////////////////////////////////////////////////////
//	
//	/*
//	 *	обработка гuстограммы уровня серого	
//	 */
//
////#ifdef _DEBUG_
////
////{
////
////	fstream f;
////
////	f.open("C:\\gray_before.txt", ios::out);
////
////	for (k = 0; k < 256; ++k) {
////
////		f << k << "\t" << intenceGraph[k] <<  endl << flush;
////
////	}	
////
////	f.close();
////
////}
////	
////#endif
//
//	int I1, I2, I3, I4, I5, I6, I7;
//
//	// сглаживание
//	for (k = 0; k <= MAX_INTENS_VAL; ++k)
//	{
//		I1 = intenceGraph[k];
//		I2 = ( k + 1 ) <= MAX_INTENS_VAL ? intenceGraph[k + 1] : 0;
//		I3 = ( k + 2 ) <= MAX_INTENS_VAL ? intenceGraph[k + 2] : 0;
//		I4 = ( k + 3 ) <= MAX_INTENS_VAL ? intenceGraph[k + 3] : 0;
//		I5 = ( k + 4 ) <= MAX_INTENS_VAL ? intenceGraph[k + 4] : 0;
//		I6 = ( k + 5 ) <= MAX_INTENS_VAL ? intenceGraph[k + 5] : 0;
//		I7 = ( k + 6 ) <= MAX_INTENS_VAL ? intenceGraph[k + 6] : 0;
//		
//		intenceGraph[k] = int( float( I1 + I2 + I3 + I4 + I5 + I6 + I7 ) / 5 );
//	}
//	//
//
//	//////////////////////////////////////////////////////////////////////////
//
//	/*
//	 *	поиск локальных максимумов	
//	 */
//
//	std::vector<int> locMaxs;
//
////	for (k = 0; k <= MAX_INTENS_VAL - 2; ++k) {
////
////		if ( intenceGraph[k] < intenceGraph[k + 1] && intenceGraph[k + 2] < intenceGraph[k + 1] ) {
////			
////			locMaxs.push_back(k + 1);
////			//intenceGraph[k + 1] = 50000;
////
////		}
////
////	}
//
//	// процедура поиска с учётом наличия "горизонтальных плато" в распределениях
//	k = 0;
//	int k0, k1, k2;
//
//	while ( k <= MAX_INTENS_VAL - 2 )
//	{
//		k0 = k;
//		k1 = k + 1;
//		k2 = k + 2;
//
//		// если первый максимум в нуле и обрезан (есть только убывающая часть)
//		if ( k == 0 && intenceGraph[k] != 0 && intenceGraph[k] >= intenceGraph[k + 1] )
//		{
//			// идём по плато вправо
//			while ( intenceGraph[k] == intenceGraph[k + 1] ) ++k;
//			//
//
//			// если после плато пошёл опять рост - гоним дальше
//			if ( intenceGraph[k] < intenceGraph[k + 1] ) {++k;continue;}
//
//			locMaxs.push_back(k);
//			//intenceGraph[k1] = 50000;
//
//		// обработка "нормальных максимумов"
//		}
//		else if ( intenceGraph[k0] < intenceGraph[k1] && intenceGraph[k2] <= intenceGraph[k1] )
//		{
//			// идём по плато вправо
//			while ( intenceGraph[k + 1] == intenceGraph[k + 2] && k < ( MAX_INTENS_VAL - 2 ) ) ++k;
//			//
//
//			k1 = k + 1;
//			k2 = k + 2;
//
//			// если после плато пошёл опять рост - гоним дальше
//			if ( intenceGraph[k2] > intenceGraph[k1] ) {++k;continue;}
//
//			locMaxs.push_back(k1);
//			//intenceGraph[k1] = 50000;
//		}
//
//		++k;
//	}
//	//
//
//	//////////////////////////////////////////////////////////////////////////
//
//	// распечатка всех максимумов
//	std::vector<int>::iterator it, it_end = locMaxs.end();
//	ATLTRACE(_T("//////////////////////////////////////////////////////////////////////////\n"));
//	ATLTRACE(_T("local maximums:\n"));
//	for (it = locMaxs.begin(); it != locMaxs.end(); ++it) ATLTRACE("%d\n", *it);
//	ATLTRACE(_T("//////////////////////////////////////////////////////////////////////////\n"));
//	//
//
//	if ( locMaxs.size() < 2 ) return 0;
//
//	iThreshold = ( locMaxs[0] + *locMaxs.rbegin() ) / 2;
//
//}

m_grayValueAuto = 0;
const int numLinImageZones = 15;
for (int s = 0; s < numLinImageZones; ++s)
for (int t = 0; t < numLinImageZones; ++t)
{
	int w_sz = MAX_INTENS_VAL + 1;

	std::vector<float> g(w_sz, 0.f);

	for (j = h / numLinImageZones * s; j < h / numLinImageZones * ( s + 1 ); j++)
		for (i = w / numLinImageZones * t; i < w / numLinImageZones * ( t + 1 ); i++)
			g[pFFrame[j * w + i]]++;

	//copy(g.begin(), g.end(), ostream_iterator<float>(cout, " "));
	//cout << endl;

	//getch();

	// нормировка распределения + вычислить mu(max)
	float w_norm = 0.f;
	float mu_max = 0.f;
	for (i = 0; i < w_sz; ++i)
	{
		w_norm += g[i];
		mu_max += i * g[i];
	}
	//cout << "w_norm = " << w_norm << " mu_max = " << mu_max << endl;
	for (i = 0; i < w_sz; ++i) g[i] /= w_norm;
	mu_max /= w_norm;

	//copy(g.begin(), g.end(), ostream_iterator<float>(cout, " "));
	//cout << endl;

	//cout << "sum = " << accumulate(g.begin(), g.end(), 0.) << endl;

	// вычислить p0(x0), mu(x0)		
	float p0_x0 = g[0], mu_x0 = 0.f;
	float numerator, denominator, sigma, sigma_max = 0.f;
	int imax_elem = 0;
	for (i = 1; i < w_sz - 1; ++i)
	{
		p0_x0 += g[i];
		mu_x0 += i * g[i];

		numerator = ( mu_max * p0_x0 - mu_x0 );
		denominator = ( p0_x0 * ( 1 - p0_x0 ) );

		sigma = ( 0 == denominator ) ? 0 : numerator * numerator / denominator;
		if ( sigma > sigma_max ) {sigma_max = sigma; imax_elem = i;}

		//cout << i << " " << p0_x0 << " " << mu_x0 << " " << sigma << endl; 
		//	<< numerator << " " << denominator << endl;
	}

	m_grayValueAuto += imax_elem;

	if ( !m_grayThrshld )
	{		
		m_grayValueCurrent = imax_elem;
	}
	else
	{
		m_grayValueCurrent = m_grayValueManual;
	}

	// ЗАЛИПУХА!!!
	//if ( m_id == 0 || m_id == 1 || m_id == 2 ) iThreshold = 105;
	//ATLTRACE("MarkerSearcherID = %d\n", m_id);
	//

	ATLTRACE(_T("m_grayValueCurrent = %d\n"), m_grayValueCurrent);
	ATLTRACE("//////////////////////////////////////////////////////////////////////////\n");

	//{
	//	CString val, id;
	//	val.Format("%d", iThreshold);
	//	id.Format("%d", m_id + 1);
	//	//if ( 10 == m_id + 1 )
	//		MessageBox(0, (LPCSTR)val, id, 0);		
	//}

	//////////////////////////////////////////////////////////////////////////
	// Определяем тип мишени (черные метки или белые метки)
	const bool m_dotsType = m_pCT->get_bgType();

	if ( !m_dotsType )
	{
		//for (j = 0; j < h; j++)
			//for (i = 0; i < w; i++)
		for (j = h / numLinImageZones * s; j < h / numLinImageZones * ( s + 1 ); j++)
			for (i = w / numLinImageZones * t; i < w / numLinImageZones * ( t + 1 ); i++)
				pFFrame[j * w + i] = ( pFFrame[j * w + i] < m_grayValueCurrent ) ? MAX_INTENS_VAL : 0;
	}
	else
	{
		//for (j = 0; j < h; j++)
			//for (i = 0; i < w; i++)
		for (j = h / numLinImageZones * s; j < h / numLinImageZones * ( s + 1 ); j++)
			for (i = w / numLinImageZones * t; i < w / numLinImageZones * ( t + 1 ); i++)
				pFFrame[j * w + i] = ( pFFrame[j * w + i] < m_grayValueCurrent ) ? 0 : MAX_INTENS_VAL;
	}

}

	m_grayValueAuto /= numLinImageZones * numLinImageZones;

	if ( !m_grayThrshld )
		SEGMENTATION_INTENS_VAL = MAX_INTENS_VAL - 1;//m_grayValueAuto;
	else
		SEGMENTATION_INTENS_VAL = m_grayValueManual;

	//{
	//	CString val, id;
	//	val.Format("%d", m_grayValueAuto);
	//	id.Format("%d", SEGMENTATION_INTENS_VAL);
	//	MessageBox(0, (LPCSTR)val, id, 0);		
	//}

#ifdef _DEBUG_

{
	//////////////////////////////////////////////////////////////////////////

	intenceGraph[iThreshold] = 30000;

	//////////////////////////////////////////////////////////////////////////

	fstream f;

	f.open(m_id == 0 ? _T("C:\\Histogram_A.txt" : "C:\\Histogram_B.txt"), ios::out);

	for (k = 0; k <= MAX_INTENS_VAL; ++k)
	{
		f << k << _T("\t") << intenceGraph[k] <<  endl << flush;
	}

	f.close();

}

#endif

	return 1;
}

int MarkerSearcher::TemplateCorrelation(IImageMap *pI, std::vector<Marker> &ms)
{
	int w = pI->GetWidth();
	int h = pI->GetHeight();

	//FILE *gout;
	//gout = fopen("D:\\1\\templ_corr_res2.dat", "wb");

	//////////////////////////////////////////////////////////////////////////
	// перечитываем файл и инвертируем, так чтобы интенсивность маркера была больше
	// интенсивности фона

	if ( !pI->ReReadBmpFile() ) return 0;

	if ( !m_pCT->get_bgType() )
	{
		if ( !pI->Invert() ) return 0;
	}

	//////////////////////////////////////////////////////////////////////////

	byte corrSz = 5;
	byte half_corrSz = corrSz / 2;
	int corrSz2 = corrSz * corrSz;
	//std::vector< double > corr(corrSz2);
	double *corr = (double*)malloc(corrSz2*sizeof(double));

	int msSz = ms.size();
	for (int kk = 0; kk < msSz; ++kk)
	{
		float markDiam = ms[kk].get_d();
		int templSz = (int)markDiam;

		//делаем размер нечетным
		templSz = ( templSz % 2 == 0 ) ? templSz + 1 : templSz;

		int templSz2 = templSz * templSz;
		int half_templSz = templSz / 2;

		//////////////////////////////////////////////////////////////////////////
		// Создание шаблона маркера

		//float xx = ms[kk].get_x();
		//float yy = ms[kk].get_y();

		// начальное значение центра маркера
		float __x = ms[kk].get_x();
		float __y = ms[kk].get_y();
		ushort x = ms[kk].get_x();
		ushort y = ms[kk].get_y();

		if ( __x - x > 0.5 ) x += 1;
		if ( __y - y > 0.5 ) y += 1;

		//bool myfl = false;
		//if ( abs(x - 378) < 3 && abs(y - 619) < 3 )
		//{
		//	printf("");
		//	myfl = true;
		//}

		ATLTRACE("curr. mark: (i, x, y, d) = (%d, %d, %d, %d)\n", kk, x, y, templSz);

		ushort *imgData = pI->GetCurrentFramePtr();

		ushort markerCenterIntensity = imgData[ y * w + x ];

		std::vector<UINT> templ(templSz2, 0);

		//FILE *out;
		//if ( myfl ) out = fopen("D:\\template.dat", "wb");

		for (int j = 0; j < templSz; ++j)
		{
			for (int i = 0; i < templSz; ++i)
			{
				templ[j * templSz + i] = 
					( ( i - templSz / 2 ) * ( i - templSz / 2 ) +
					( j - templSz / 2 ) * ( j - templSz / 2 ) < templSz2 / 4 ) ?
						markerCenterIntensity : 0;

				//if ( myfl ) fprintf(out, "%d %d %d\n", i, j, templ[j * templSz + i]);
				//if ( myfl ) fprintf(out, "%d ", templ[j * templSz + i]);
			}
			//if ( myfl ) fprintf(out, "\n");
		}

		//if ( myfl ) fclose(out);

		//////////////////////////////////////////////////////////////////////////
		// Расчет кросс-корреляции изображения маркера с шаблоном

		for ( int i = 0; i < corrSz2; ++i ) corr[i] = 0;

		int ii, jj, s, t;

		UINT second;

		//FILE *out;
		//char _buf_[PATH_LENGTH];
		//sprintf(_buf_, "D:\\1\\corr_%d.dat", kk);

		//if ( kk == 29 || kk == 45 )
		//	out = fopen(_buf_, "wb");
		//if ( myfl ) out = fopen("D:\\corrfunc.dat", "wb");

		int _ii, _jj;
		float fx = 0.0f, fy = 0.0f;
		double corr_max = 0, corr_val;

		for (t = -half_corrSz; t <= half_corrSz; ++t)
		{
			for (s = -half_corrSz; s <= half_corrSz; ++s)
			{
				for (jj = 0; jj < templSz; ++jj)
				{
					for (ii = 0; ii < templSz; ++ii)
					{
						_ii = s + ii + x - half_templSz;
						_jj = t + jj + y - half_templSz;

						if ( _ii < 0 || _ii > w - 1 || _jj < 0 || _jj > h - 1 )
						{
							second = 0;
						}
						else
						{
							second = imgData[ _jj * w + _ii ];
						}

						corr[( t + half_corrSz ) * corrSz + ( s + half_corrSz )] += templ[jj * templSz + ii] * second;
					}
				}
				//if ( kk == 29 || kk == 45 )
				//	fprintf(out, "%d %d %d\n", s, t, corr[ ( t + half_corrSz ) * corrSz + ( s + half_corrSz ) ]);
				//if ( myfl) fprintf(out, "%g ", corr[ ( t + half_corrSz ) * corrSz + ( s + half_corrSz ) ]);

				corr_val = corr[ ( t + half_corrSz ) * corrSz + ( s + half_corrSz ) ];

				if ( corr_val > corr_max )
				{
					corr_max = corr_val;
					fx = s;
					fy = t;
				}

			}
			//if ( myfl ) fprintf(out, "\n");
		}

		//if ( kk == 29 || kk == 45 )
		//	fclose(out);
		//if ( myfl )	fclose(out);

		//////////////////////////////////////////////////////////////////////////
		// Подпиксельная интерполяция положения максимума кросс-корреляционной функции

		if ( fx < half_corrSz && fx > -half_corrSz && fy < half_corrSz && fy > -half_corrSz )
		{
			double log_A, log_B, log_C, notZero;

			log_A = log( (double)corr[(int)( ( fy + half_corrSz ) * corrSz + fx + half_corrSz - 1 )] );
			log_B = log( (double)corr[(int)( ( fy + half_corrSz ) * corrSz + fx + half_corrSz )] );
			log_C = log( (double)corr[(int)( ( fy + half_corrSz ) * corrSz + fx + half_corrSz + 1 )] );

			notZero = 2 * log_B - log_A - log_C;

			float tmp_x = fx;

			if ( notZero != 0 ) fx = x + fx + 0.5 * ( log_C -log_A ) / notZero;

			log_A = log( (double)corr[(int)( ( fy + half_corrSz - 1 ) * corrSz + tmp_x + half_corrSz )] );
			log_C = log( (double)corr[(int)( ( fy + half_corrSz + 1 ) * corrSz + tmp_x + half_corrSz )] );

			notZero = 2 * log_B - log_A - log_C;

			if ( notZero != 0 ) fy = y + fy + 0.5 * ( log_C - log_A ) / notZero;

			//ATLTRACE("cm: %g, %g; int: %d, %d; corr: %g, %g\n", 
			//	m_markerSet[kk].get_x(),
			//	m_markerSet[kk].get_y(),
			//	x,
			//	y,
			//	fx,
			//	fy
			//);

			//fprintf(gout, "%g %g %g %g\n", ms[kk].get_x(), ms[kk].get_y(), fx, fy);

			ms[kk] = Marker(fx, fy, markDiam);
		}
		else
		{
			ATLTRACE("out of corrSz");
		}

		//////////////////////////////////////////////////////////////////////////

	}

	free(corr);

	//fclose(gout);
	return 1;
}

int MarkerSearcher::MarkerScan(IImageMap *pI)
{
	int i, j;
	int w = pI->GetWidth();
	int h = pI->GetHeight();

	ushort* pFFrame = pI->GetCurrentFramePtr();
	
	// Пример сохранения изображения (после бинаризации)
	//char buf[256];
	//strcpy(buf, "D:\\!\\img_a");
	//pI->WriteFile(buf, NULL);
	
	int index, index2;
	int _i, _j;
	unsigned int dotArea;
	float realMinD = FLT_MAX;
	float realMaxD = 0;
	//unsigned int realMinDotArea = UINT_MAX;
	//unsigned int realMaxDotArea = 0;
	float x, y, D;
	uint curPos = 0;

	// очищаем поле маркеров
	m_markerSet.clear();
	
	// множество для пикселов области связности
	std::vector< Position > pixSet;
	pixSet.reserve(3000);//Примерное количество пикселей в маркере на изображении ~2010 (взято с запасом)
	//

	//////////////////////////////////////////////////////////////////////////

	for (i = 0; i < w; i++) 
	{	
		for (j = 0; j < h; j++)
		{
			index = j * w + i;
			
			//поиск затравки
			//т.е. начала связной области
			if ( pFFrame[index] == MAX_INTENS_VAL )
			{
				curPos = 0;

				pixSet.push_back(Position(i, j));

				//////////////////////////////////////////////////////////////////////////

				while ( curPos < pixSet.size() )
				{
					_i = pixSet[curPos].m_i;
					_j = pixSet[curPos].m_j;
					
					//вверх
					if ( _j + 1 < h )
					{
						index2 = (_j + 1) * w + _i;

						if ( pFFrame[index2] == MAX_INTENS_VAL )
						{
							pixSet.push_back(Position(_i, _j + 1));
							pFFrame[index2] = SEGMENTATION_INTENS_VAL;
						}
					}
					//

					//вниз
					if ( _j - 1 >= 0 )
					{
						index2 = (_j - 1) * w + _i;

						if ( pFFrame[index2] == MAX_INTENS_VAL )
						{
							pixSet.push_back(Position(_i, _j - 1));
							pFFrame[index2] = SEGMENTATION_INTENS_VAL;
						}
					}
					//

					//вправо
					if ( _i + 1 < w )
					{
						index2 = _j * w + _i + 1;

						if ( pFFrame[index2] == MAX_INTENS_VAL )
						{
							pixSet.push_back(Position(_i + 1, _j));
							pFFrame[index2] = SEGMENTATION_INTENS_VAL;
						}
					}
					//

					//влево
					if ( _i - 1 >= 0 )
					{
						index2 = _j * w + _i - 1;

						if ( pFFrame[index2] == MAX_INTENS_VAL )
						{
							pixSet.push_back(Position(_i - 1, _j));
							pFFrame[index2] = SEGMENTATION_INTENS_VAL;
						}
					}

					curPos++;
				}

				//////////////////////////////////////////////////////////////////////////
				
				//вычисляем центры масс маркеров, их характерный диаметр и сохраняем
				//информацию о маркерах в контейнер

				//Вычисляет координаты центра масс маркера
				ImageCoord ic(std::for_each(pixSet.begin(), pixSet.end(), positionAverageF()).result());
				x = ic.x;
				y = ic.y;

				dotArea = pixSet.size();
				//ATLTRACE("dotArea = %d\n", dotArea);
				//Вычисляет характерный диаметр маркера
				//Обратное действие от S = PI * D^2 / 4 будет D = sqrt(4 * S / PI)
				D = (float)sqrt(4 * dotArea / PI);
				//ATLTRACE("D = %g\n", D);

				realMinD = min(realMinD, D);
				realMaxD = max(realMaxD, D);

				// диаметр маркера должен лежать в указанном пользователем диапазоне
				if ( D > m_minDotD && D < m_maxDotD )
				{
					//на изображении 035nofz5.img из-за несовершенства алгоритма бинаризации и
					//алгоритма поиска связных областей (вверх, вниз, вправо, влево) к
					//большой связной области примыкает один пиксел по-диаганали
					//(см. файл CalibrationProc\bugs!\bug1.bmp),
					//который и определяется рядом с большой связной областью
					//всё это приводит к тому, что расстояни м/у ними оказывается равным 0
					//и алгоритм ломается
					//Вообще, надо переделывать "Бинаризацию" или вводить в "Поиск связных областей"
					//шаги по-диаганали (что увеличит время расчёта)

					m_markerSet.push_back(Marker(x, y, D));
				}
				//
				pixSet.clear();
			}
		}
	}

	// записывает сегментированное изображение на неиспользованый кадр двухкадрового изображения (для однокадрового - ничего не делает)
	if ( m_bSaveSegmentation ) pI->WriteFile();

	//////////////////////////////////////////////////////////////////////////

	/*
	 *	Валидация маркеров
	 */

	int numMarkersInField = (int)m_markerSet.size();
	m_meanD = std::accumulate(m_markerSet.begin(), m_markerSet.end(), 0.0f, sumDiamsF) / numMarkersInField;

	ATLTRACE(_T("meanD = %g\n"), m_meanD);

	ATLTRACE(_T("Marker diameter filter...\n"));
	m_markerSet.erase(std::remove_if(m_markerSet.begin(), m_markerSet.end(), 
		diameterValidMarkerF(m_meanD, m_Q)), m_markerSet.end());

	ATLTRACE(_T("Marker boundary position filter...\n"));
	m_markerSet.erase(std::remove_if(m_markerSet.begin(), m_markerSet.end(), 
		boundaryValidMarkerF(m_boundaryPartX, m_boundaryPartY, w, h)), m_markerSet.end());

	//////////////////////////////////////////////////////////////////////////
	
	/*
	 * Уточнение положения центра маркера путем кросс-корреляции с его шаблоном
	 */

	if ( m_templCorrToogle && !TemplateCorrelation(pI, m_markerSet) ) return 0;

	//////////////////////////////////////////////////////////////////////////

	// По-умолчанию поиск центра ведется путем выбора маркера максимального диаметра 
	// в окне расположенном в центре изображения размером указанном ниже
	Position centerPosition(w/2, h/2);
	//Position centerWindow(w/2, h/2);
	Position centerWindow(w, h);

	// Поиск центрального маркера выполняется либо автоматически, либо вблизи заданой позиции
	if ( m_manualCenterMarker )
	{
		centerPosition = m_manualCenterPos;
		centerWindow = Position(m_windowSz, m_windowSz);
	}

	//char __buf[255];
	//sprintf(__buf, "centerPosition: %d %d; centerWindow : %d %d", centerPosition.m_i, centerPosition.m_j, centerWindow.m_i, centerWindow.m_j);
	//MessageBox(0, __buf, "", 0);

	MarkerColumn _tmpMarkerSet(m_markerSet);
	_tmpMarkerSet.erase(std::remove_if(_tmpMarkerSet.begin(), _tmpMarkerSet.end(), 
		markersOutOfCenterWindowF(centerPosition, centerWindow)), _tmpMarkerSet.end());
	ATLTRACE(_T("Find marker with max diameter...\n"));

	////char __buf[255];
	//sprintf(__buf, "size: %d", _tmpMarkerSet.size());
	//MessageBox(0, __buf, "", 0);

	//for (int iii = 0; iii < _tmpMarkerSet.size(); ++iii)
	//{
	//	sprintf(__buf, "markers: %g %g", _tmpMarkerSet[iii].get_x(), _tmpMarkerSet[iii].get_y());
	//	MessageBox(0, __buf, "", 0);
	//}

	MarkerColumnIt centerMarkerIt = std::max_element(_tmpMarkerSet.begin(), _tmpMarkerSet.end(), dCompareMarkerF);
	if ( centerMarkerIt != _tmpMarkerSet.end() ) m_OriginMarker = *centerMarkerIt;

	//sprintf(__buf, "center: %g %g %g %g", m_OriginMarker.get_x(), m_OriginMarker.get_y(), m_OriginMarker.get_d(), m_meanD);
	//MessageBox(0, __buf, "", 0);

	// Старое
	//ATLTRACE(_T("Sort by distance to specified marker...\n"));
	//std::sort(_tmpMarkerSet.begin(), _tmpMarkerSet.end(), compareByDistToMarkerF(Marker(m_manualCenterPos.m_i, m_manualCenterPos.m_j, 0)));		
	//
	
	//////////////////////////////////////////////////////////////////////////
	
	/*
	 * ВЫВОД СООБЩЕНИЯ О НЕДОСТАТОЧНОМ КОЛИЧЕСТВЕ МАРКЕРОВ
	*/

	// минимальное количество маркеров 5 (центральный маркер и 4 ближайших к нему маркера)
	if ( numMarkersInField < m_minDotCount )
	{
		ATLTRACE(_T("m_markerSet.size() < m_minDotCount failed...\n"));
		
		TCHAR caption[125];
		_sntprintf(caption, 125, _T("Calibration Info (plane %d)"), m_id + 1);
		
		// Минимальное ограничивающее пороговое значение фильтрующее маркеры
		float minQ = min(realMaxD / m_meanD, m_meanD / realMinD);

		tstring strGrayLevelToDo;
		strGrayLevelToDo = m_pCT->get_bgType() ? _T("Decrease") : _T("Increase");
		
		TCHAR text[2048];
		_sntprintf(text, 2048, 
_T("\"Min dot count\" failed\n \
Recommendations (in brackets last user input):\n \
Decrease \"Min dot count\" = %d (%d), \"Boundary part\" (%g, %g) or \"Min marker D\" = %.1f (%.1f)\n \
Increase \"Diff. from mean D\" = %.1f (%.1f), \"Max marker D\" = %.1f (%.1f)\n \
%s binarization threshold is %d (%d) \n \
Mean marker diameter after min/max validation is %.1f pix \n \
Also you can adjust manually the grayscale threshold binarization level\n \
You may have an inverted calibration target type (dots color)"), 
			numMarkersInField, 
			m_minDotCount,
			m_boundaryPartX, 
			m_boundaryPartY, 
			realMinD, 
			m_minDotD, 
			minQ,
			m_Q,
			realMaxD,
			m_maxDotD,
			strGrayLevelToDo.c_str(),
			m_grayValueAuto,
			m_grayValueCurrent,
			m_meanD
		);
		
		::MessageBox(NULL, text, caption, MB_OK | MB_ICONINFORMATION);
		
		return 0; 
	}

	//////////////////////////////////////////////////////////////////////////
	
	/*
	 * ВЫВОД СООБЩЕНИЯ, О ТОМ ЧТО НАЙДЕННЫЙ ЦЕНТРАЛЬНЫЙ МАРКЕР НЕ ПРОШЕЛ КРИТЕРИЙ ОТБОРА
	 */

	// проверяем удовлетворяет ли центральный маркер
	// критерию по величине площади
	
	if ( m_OriginMarker.get_d() < m_meanD * m_OriginAreaThreshold )
	{
		ATLTRACE(_T("%g > %g\n"), m_OriginMarker.get_d(), m_meanD * m_OriginAreaThreshold);

		TCHAR caption[125];
		_sntprintf(caption, 125, _T("Calibration Info (plane %d)"), m_id + 1);

		tstring strGrayLevelToDo;
		strGrayLevelToDo = !m_pCT->get_bgType() ? _T("Decrease") : _T("Increase");

		TCHAR text[2048];
		_sntprintf(text, 2048,
_T("\"Center marker D threshold\" failed\n \
Recommendations (in brackets last user input):\n \
Decrease \"Center marker D threshold\" = %.1f (%.1f)\n \
%s binarization threshold is %d (%d) \n \
Mean marker diameter after min/max validation is %.1f pix \n \
Also you can adjust manually the center marker position"),
			m_OriginMarker.get_d() / m_meanD, 
			m_OriginAreaThreshold,
			strGrayLevelToDo.c_str(),
			m_grayValueAuto,
			m_grayValueCurrent,
			m_meanD
		);

		::MessageBox(NULL, text, caption, MB_OK | MB_ICONINFORMATION);

		return 0;
	}


	ATLTRACE(_T("ORIGIN_AREA_THRESHOLD %g > %g\n"), m_OriginMarker.get_d(), m_meanD * m_OriginAreaThreshold);
	//ATLTRACE(_T("CENTER_MARKER_POSITION (x, y) = (%d, %d)\n"), m_Origin.m_i, m_Origin.m_j);

	/////////////////////////////////////////////////////////////////////////

	{
		int j;

		// Центральный маркер
		Marker center = m_OriginMarker;
		//

		ATLTRACE(_T("CENTER_MARKER: (x, y, d) = (%g, %g, %g)\n"), center.get_x(), center.get_y(), center.get_d());

		MarkerColumn tmpMarkerSet2(m_markerSet);

		ATLTRACE(_T("Sort by distance from center marker...\n"));

		auto compareByDistToMarkerF = [m = center](const Marker& m1, const Marker& m2) ->bool
		{
			return _hypot(m1.get_x() - m.get_x(), m1.get_y() - m.get_y())
				< _hypot(m2.get_x() - m.get_x(), m2.get_y() - m.get_y());
		};

		std::sort(tmpMarkerSet2.begin(), tmpMarkerSet2.end(), compareByDistToMarkerF);

		//////////////////////////////////////////////////////////////////////////
		// СТОП!!!: ОСНОВНАЯ ОШИБКА ПОДХОДА ПО ПОИСКУ 4 БЛИЖАЙШИХ МАРКЕРОВ, ТО ЧТО
		// В КАЧЕСТВЕ БЛИЖАЙШЕГО МОЖЕТ ОПРЕДЕЛИТЬСЯ ДИАГОНАЛЬНЫЙ МАРКЕР (ЗА СЧЕТ РАСТЯЖЕНИЯ ИЗОБРАЖЕНИЯ)
		// ТОГДА ПРЕДВАРИТЕЛЬНАЯ МОДЕЛЬ БУДЕТ НАЙДЕНА НЕПРАВИЛЬНО
		// (ПРИМЕР ОШИБКИ ЛЕЖИТ ЗДЕСЬ: "\\big-piv\USERS\!exchange\misha\2007_11_26\image4_a.bmp", при указании 
		// положения центрального маркера вручную (349, 330) калибровка камеры рассчитывается неверно
		//////////////////////////////////////////////////////////////////////////

		//for (j = 0; j < (int)tmpMarkerSet2.size(); ++j)
		//{
		//	ATLTRACE(_T("MARKER: (x, y, d) = (%g, %g, %g)\n"), tmpMarkerSet2[j].get_x(), tmpMarkerSet2[j].get_y(), tmpMarkerSet2[j].get_d());
		//}

		//03.11.2016
		MarkerColumn tmpMarkerSet3(4);
		std::copy(tmpMarkerSet2.begin() + 1, tmpMarkerSet2.begin() + 5, tmpMarkerSet3.begin());
		//replaced by
		//MarkerColumn tmpMarkerSet3;
		//tmpMarkerSet3.push_back(tmpMarkerSet2[1]);
		////const double max_cos_value = 0.8660254;//cos(30)
		//const double max_cos_value = 0.966;//cos(15)
		//int i = 2;
		//do {
		//	ATLASSERT(i < tmpMarkerSet2.size());
		//	float vx = tmpMarkerSet2[i].get_x() - center.get_x();
		//	float vy = tmpMarkerSet2[i].get_y() - center.get_y();
		//	bool flag = true;
		//	for(int l = 0; l < tmpMarkerSet3.size(); ++l)
		//	{
		//		float vx2 = tmpMarkerSet3[l].get_x() - center.get_x();
		//		float vy2 = tmpMarkerSet3[l].get_y() - center.get_y();
		//		if((vx*vx2+vy*vy2)/sqrt((vx*vx+vy*vy)*(vx2*vx2+vy2*vy2)) > max_cos_value)
		//		{
		//			flag = false;
		//			break;
		//		}
		//	}
		//	if(flag)
		//	{
		//		tmpMarkerSet3.push_back(tmpMarkerSet2[i]);
		//	}
		//	++i;
		//}
		//while(tmpMarkerSet3.size() < 4);
		//replace end

		ATLTRACE(_T("Sort by x...\n"));
		std::sort(tmpMarkerSet3.begin(), tmpMarkerSet3.end(), compareMarkerXF);

		// Левый и правый маркеры по оси X
		Marker left = *tmpMarkerSet3.begin();
		Marker right = *tmpMarkerSet3.rbegin();
		//

		//for (j = 0; j < (int)tmpMarkerSet3.size(); ++j)
		//{
		//	ATLTRACE(_T("MARKER: (x, y, d) = (%g, %g, %g)\n"), tmpMarkerSet3[j].get_x(), tmpMarkerSet3[j].get_y(), tmpMarkerSet3[j].get_d());
		//}

		ATLTRACE(_T("Sort by y...\n"));
		std::sort(tmpMarkerSet3.begin(), tmpMarkerSet3.end(), compareMarkerYF);

		// Нижний и верхний маркеры по оси Y
		Marker bottom = *tmpMarkerSet3.begin();
		Marker top = *tmpMarkerSet3.rbegin();
		//

		//for (j = 0; j < (int)tmpMarkerSet3.size(); ++j)
		//{
		//	ATLTRACE(_T("MARKER: (x, y, d) = (%g, %g, %g)\n"), tmpMarkerSet3[j].get_x(), tmpMarkerSet3[j].get_y(), tmpMarkerSet3[j].get_d());
		//}

		//////////////////////////////////////////////////////////////////////////

		MarkerColumnIt it, it_end;
		ModelCoord c;

		std::vector< CoordPair > cpSet;

		//////////////////////////////////////////////////////////////////////////

//18.01.17
		//
		//float dx = top.get_x() - bottom.get_x();
		//top.set_x(center.get_x() - dx / 2);
		//bottom.set_x(center.get_x() + dx / 2);
		//
		const int precalc_size = 0;
		cpSet.reserve(5);
		cpSet.push_back(CoordPair(ImageCoord(center.get_x(), center.get_y()), ModelCoord(0.0f, 0.0f, 0.0f)));
		cpSet.push_back(CoordPair(ImageCoord(left.get_x(), left.get_y()), ModelCoord(-1.0f, 0.0f, 0.0f)));
		cpSet.push_back(CoordPair(ImageCoord(right.get_x(), right.get_y()), ModelCoord(1.0f, 0.0f, 0.0f)));
		cpSet.push_back(CoordPair(ImageCoord(bottom.get_x(), bottom.get_y()), ModelCoord(0.0f, -1.0f, 0.0f)));
		cpSet.push_back(CoordPair(ImageCoord(top.get_x(), top.get_y()), ModelCoord(0.0f, 1.0f, 0.0f)));

		//const int precalc_size = 2;
		//double v0_x = (right.get_x() - left.get_x())/2;
		//double v0_y = (right.get_y() - left.get_y())/2;
		//double v1_x = (top.get_x() - bottom.get_x())/2;
		//double v1_y = (top.get_y() - bottom.get_y())/2;

		//std::vector<std::vector<int>> d;
		//d.reserve( sqr(2*precalc_size +1) );
		//for(int y = -precalc_size; y <= precalc_size; ++y)
		//{
		//	for(int x = -precalc_size; x <= precalc_size; ++x)
		//	{
		//		std::vector<int> tmp(2);
		//		tmp[0] = x;
		//		tmp[1] = y;
		//		d.push_back(tmp);
		//	}
		//}
		//double err = min(sqr(v0_x)+sqr(v0_y), sqr(v1_x)+sqr(v1_y));
		//for(int k = 0; k < d.size(); ++k)
		//{
		//	Marker m(center.get_x()+v0_x*d[k][0]+v1_x*d[k][1],center.get_y()+v0_y*d[k][0]+v1_y*d[k][1],0);
		//	std::sort(tmpMarkerSet2.begin(), tmpMarkerSet2.end(), compareByDistToMarkerF(m));	//можно заменить сортировку на поиск ближайшего.
		//	ImageCoord ic(tmpMarkerSet2[0].get_x(), tmpMarkerSet2[0].get_y());
		//	if(sqr(m.get_x()-ic.x) + sqr(m.get_y()-ic.y) <= err)
		//	{
		//		ModelCoord mc(d[k][0],d[k][1],0.0f);
		//		cpSet.push_back(CoordPair(ic, mc));
		//	}
		//}
//18.01.17 replace end

		for (auto &cpit : cpSet) ATLTRACE("center seed: %g %g %g %g\n", cpit.i.x, cpit.i.y, cpit.m.X, cpit.m.Y);// 03.11.2016

		if ( !PreModelFit(cpSet) ) return 0;
		
		//////////////////////////////////////////////////////////////////////////
		
		TCHAR buf_[25];
		float X, Y;

		//unsigned char N = 1;
		unsigned char upperBoundN = sqrt((float)m_markerSet.size()) / 2;

		for (unsigned char N = precalc_size + 1; N <= upperBoundN; N+=2)	//9.11.16	//18.01.17
		{
			cpSet.clear();
			cpSet.reserve((2 * N + 1) * (2 * N + 1));

			int _c_ = 0;

			//////////////////////////////////////////////////////////////////////////
			
			it = m_markerSet.begin();
			it_end = m_markerSet.end();

			for (; it != it_end; ++it)
			{
				X = ShiftToDiscrete(it->get_X());
				Y = ShiftToDiscrete(it->get_Y());
				
				if ( X <= N && X >= -N && Y <= N && Y >= -N )
				{
					cpSet.push_back(CoordPair(ImageCoord(it->get_x(), it->get_y()), ModelCoord(X, Y, 0.f)));
					++_c_;
					ATLTRACE("%g %g %g %g %d = %g %g %d\n", it->get_x(), it->get_y(), X, Y, N, it->get_X(), it->get_Y(), _c_);// 03.11.2016
					int curSize = (2 * N + 1) * (2 * N + 1);
#ifndef NDEBUG
//					assert( _c_ <= curSize );
#endif
					if ( _c_ > curSize ) 
					{
						std::stringstream sstr;sstr << "fitted=" << _c_ << " > curSize=" << curSize;
						ATLTRACE(sstr.str().c_str());
						//::MessageBox(0, sstr.str().c_str(), "", 0);
						break;
					}
				}
			}
			
			//_stprintf(buf_, _T("%d"), _c_);
			//MessageBox(NULL, buf_, NULL, MB_OK);

			//{
			//	CString val;
			//	val.Format("%d %d", N, size_t(sqrt((float)m_markerSet.size()) / 2));
			//	MessageBox(0, (LPCSTR)val, "PreModel", 0);
			//}

			if ( !PreModelFit(cpSet, false) ) return 0;
		}
		
		//////////////////////////////////////////////////////////////////////////

		//_c_ = 0;
		//N = 3;
		//cpSet.clear();
		//cpSet.reserve((2 * N + 1) * (2 * N + 1));

		////////////////////////////////////////////////////////////////////////////

		//it = m_markerSet.begin();
		//it_end = m_markerSet.end();

		//for (; it != it_end; ++it)
		//{
		//	X = ShiftToDiscrete(it->get_X());
		//	Y = ShiftToDiscrete(it->get_Y());

		//	if ( X <= N && X >= -N && Y <= N && Y >= -N )
		//	{
		//		cpSet.push_back(CoordPair(ImageCoord(it->get_x(), it->get_y()), ModelCoord(X, Y, 0.f)));
		//		++_c_;
		//	}
		//}
		//
		////_stprintf(buf_, _T("%d"), _c_);
		////MessageBox(NULL, buf_, NULL, MB_OK);

		//if ( !PreModelFit(cpSet) ) return 0;

		//////////////////////////////////////////////////////////////////////////
		// вывод данных для отладки

		//unsigned int counter = 0;

		//float X_last = it->get_X(), Y_last = it->get_Y();

		//FILE *out;
		//TCHAR _buf_[256];
		//_stprintf(_buf_, _T("D:\\!\\%d.txt"), m_id + 1);
		//out = _tfopen(_buf_, _T("wb"));
		//
		//std::sort(m_markerSet.begin(), m_markerSet.end(), directCompareMarkerXWCF());
		//
		//it = m_markerSet.begin();
		//it_end = m_markerSet.end();

		//// Устанавливаем маркерам их положение в мировой системе координат по предварительной модели отображения
		//for (; it != it_end; ++it)
		//{
		//	_ftprintf(out, "%u\t%g\t%g\t%g\t%g\n", counter++, it->get_X(), it->get_Y(), it->get_x(), it->get_y()/*X_last - it->get_X(), Y_last - it->get_Y()*/);
		//	X_last = it->get_X();
		//	Y_last = it->get_Y();
		//}

		//fclose(out);

		//////////////////////////////////////////////////////////////////////////
	}

	//////////////////////////////////////////////////////////////////////////	

	return 1;
}

int MarkerSearcher::BindImageMarkersToModel()
{
	// проверка и выход если выполнение невозможно
	if ( m_pCT == 0 ) return 0;

	const float dotSpacing = m_pCT->get_dotSpacing();

	float x, y, X, Y, Z;
	ModelCoord mc;

	// Поиск ограничивающей области
	MarkerColumnIt it = m_markerSet.begin(), it_end = m_markerSet.end();

	//////////////////////////////////////////////////////////////////////////

	// 2008_05_05 Коррекция смещения столбцов по горизонтали для трехуровневых мишеней
	
	float X_corr = m_pCT->get_levelDist() / m_pCT->get_dotSpacing() * 0.84;// tag(40deg)~0.84

	// различаем правое и левое направление обзора
	if ( m_viewToggle == 0 ) X_corr = 0;				// Unknown view
	else if ( m_viewToggle == 1 ) X_corr = -X_corr;		// Left view
	/*else if ( m_viewToggle == 2 ) X_corr = X_corr;*/	// Right view

	// flipXY
    float tmp;

	// **************
	//ModelCoord tmp_mc;
	//PreciseTargPointsLoader ptpLoader("E:\\qqqq\\scanProc\\out2.txt");
	//PreciseTargPointsLoader ptpLoader("E:\\qqqq\\OpenCVcalib\\open_cv_test2.dat");
	//std::ofstream f_f_f("C:\\f_f_f.txt");
	// **************

	for (; it != it_end; ++it)
	{
		X = it->get_X();
		Y = it->get_Y();

		// 2008_05_05 Коррекция смещения столбцов по горизонтали для трехуровневых мишеней
		if ( m_pCT->get_multiLevel() )
			if ( it->get_level() == -1 )
				X -= X_corr;
			if ( it->get_level() == 1 )
				X += X_corr;

		X = ShiftToDiscrete(X) * m_xOrient * dotSpacing;
		Y = ShiftToDiscrete(Y) * m_yOrient * dotSpacing;

		// **************
 		//tmp_mc = ptpLoader.get(X / dotSpacing, Y / dotSpacing);		
		//f_f_f << it->get_x() << " " << it->get_y() << " " << X << " " << Y << std::endl;
		//f_f_f << it->get_x() << " " << it->get_y() << " " << tmp_mc.X << " " << tmp_mc.Y << std::endl;
		//f_f_f << X << " " << Y << " " << tmp_mc.X << " " << tmp_mc.Y << std::endl;
		//X = tmp_mc.X; Y = tmp_mc.Y;
		// **************

		//flipXY
		if ( m_bFlipXY )
		{
			tmp = X;
			X = Y;
			Y = tmp;
		}

		Z = m_pCT->get_multiLevel() ? it->get_level() * m_pCT->get_levelDist() : m_Z;

		// Преобразование для обратной стороны мишени
		if ( m_pCT->get_multiLevel() && m_multilevelBackSideEnable )
		{
			Z = -( Z + m_pCT->get_multiLevelZerroThickness() );
		}

		// Внимание!!!
		// Выполняется поворот плоскости мишени и ее смещения для коррекции рассогласования
		mc = m_mCorr(ModelCoord(X, Y, Z));

		// Для отображения на экране
		it->set_X(mc.X);
		it->set_Y(mc.Y);
		it->set_Z(mc.Z);
		it->set_used(true);
	}

	std::sort(m_markerSet.begin(), m_markerSet.end(), compareMarkerWCSF);

	// скопировать опрорные точки в виде пар (x,y)<->(X,Y,Z) в массив m_cpSet
	m_cpSet.clear();
	it = m_markerSet.begin(), it_end = m_markerSet.end();
	for (; it != it_end; ++it) m_cpSet.push_back(
		CoordPair(ImageCoord(it->get_x(), it->get_y()), ModelCoord(it->get_X(), it->get_Y(), it->get_Z()))
		);

	return 1;
}

ModelRegion MarkerSearcher::CommonModelRegionQuadrangle(const Quadrangle &lq, const Quadrangle &rq)
{
	ModelRegion r;

	float r_max_x = rq.lb.X > rq.lt.X ? rq.lb.X : rq.lt.X;
	float l_max_x = lq.lb.X > lq.lt.X ? lq.lb.X : lq.lt.X;

	r.lb.X = r_max_x > l_max_x ? r_max_x : l_max_x;

	float r_max_y = rq.lb.Y > rq.rb.Y ? rq.lb.Y : rq.rb.Y;
	float l_max_y = lq.lb.Y > lq.rb.Y ? lq.lb.Y : lq.rb.Y;

	r.lb.Y = r_max_y > l_max_y ? r_max_y : l_max_y;

	float r_min_x = rq.rt.X < rq.rb.X ? rq.rt.X : rq.rb.X;
	float l_min_x = lq.rt.X < lq.rb.X ? lq.rt.X : lq.rb.X;

	r.rt.X = r_min_x < l_min_x ? r_min_x : l_min_x;

	float r_min_y = rq.lt.Y < rq.rt.Y ? rq.lt.Y : rq.rt.Y;
	float l_min_y = lq.lt.Y < lq.rt.Y ? lq.lt.Y : lq.rt.Y;

	r.rt.Y = r_min_y < l_min_y ? r_min_y : l_min_y;

	return r;
}

int MarkerSearcher::PreModelFit(std::vector< CoordPair > &cpSet, bool bCalcDist)
{
	MarkerColumnIt it, it_end;
	ModelCoord c;

	DLTRevModelFit mfr(cpSet);

	if ( !mfr.DoFit() )
	{
		// параметры предварительной модели не найдены
		//ATLASSERT(0);
		return 0;
	}

	it = m_markerSet.begin();
	it_end = m_markerSet.end();

	//////////////////////////////////////////////////////////////////////////
	// Оценка коэффициента радиальной дисторсии по отклонению маркеров от истинного положения
	if ( bCalcDist )
	{
		float X, Y, kx = 0.f, kx_cur, r2;
		int count = 0;

		for (; it != it_end; ++it)
		{
			c = mfr.get_ModelCoord(ImageCoord(it->get_x(), it->get_y()));

			X = ShiftToDiscrete(it->get_X());
			Y = ShiftToDiscrete(it->get_Y());

			r2 = c.X * c.X + c.Y * c.Y;
			kx_cur = ( 1 - fabs(X) / fabs(c.X) ) / r2;
			if ( kx_cur > 0 && kx_cur < 0.0002 ) {kx += kx_cur;count++;}
			//ky += ( Y / c.Y - 1 ) / r2;

			//char buf_[256];
			////sprintf(buf_, "%g %g %g %g", c.X, X, c.Y, Y);
			//sprintf(buf_, "%g", ( fabs(X) / fabs(c.X) - 1 ) / r2);
			//::MessageBox(0, buf_, "one point", 0);
		}

		kx /= count;
		//ky /= m_markerSet.size();
		//k = ( kx + ky ) / 2;

		//char buf[256];
		////sprintf(buf, "%g %g %g", kx, ky, k);
		//sprintf(buf, "%g %d", kx, count);
		//::MessageBox(0, buf, "avearage", 0);

		m_radialDistCoeffAuto = kx;
	}
	//////////////////////////////////////////////////////////////////////////

	//char buf[256];
	////sprintf(buf, "%g %g %g", kx, ky, k);
	//sprintf(buf, "%g", m_radialDistCoeff);
	//::MessageBox(0, buf, "avearage", 0);

	it = m_markerSet.begin();

	// Устанавливаем маркерам их положение в мировой системе координат по предварительной модели отображения
	for (; it != it_end; ++it)
	{
		c = mfr.get_ModelCoord(ImageCoord(it->get_x(), it->get_y()), 0.f/*m_radialDistCoeff*/); // AAAAAAAAA!!! ZALIPUHA
		it->set_X(c.X);
		it->set_Y(c.Y);
	}

	printf("Current calibError=%g\n", mfr.get_erMCS());

	return 1;
}

Quadrangle MarkerSearcher::GetQuadrangleRegion(int coordToggle, bool bFlipXY, ModelCoord &c1, ModelCoord &c2, ModelCoord &c3, ModelCoord &c4)
{
	Quadrangle q;
	// выбираем правильное положение углов четырёхугольника
	switch (coordToggle)
	{
	case 0:
		//m_xOrient = xRight;
		//m_yOrient = yUp;

		if ( false == bFlipXY )
		{
			q.lb = c1;
			q.lt = c2;
			q.rt = c3;
			q.rb = c4;
		}
		else
		{
			q.lb = c1;
			q.lt = c4;
			q.rt = c3;
			q.rb = c2;
		}

		break;
	case 1:
		//m_xOrient = xLeft;
		//m_yOrient = yUp;

		if ( false == bFlipXY )
		{
			q.lb = c4;
			q.lt = c3;
			q.rt = c2;
			q.rb = c1;
		}
		else
		{
			q.lb = c4;
			q.lt = c1;
			q.rt = c2;
			q.rb = c3;
		}

		break;
	case 2:
		//m_xOrient = xRight;
		//m_yOrient = yDown;

		if ( false == bFlipXY )
		{
			q.lb = c2;
			q.lt = c1;
			q.rt = c4;
			q.rb = c3;
		}
		else
		{
			q.lb = c2;
			q.lt = c3;
			q.rt = c4;
			q.rb = c1;
		}

		break;
	case 3:
		//m_xOrient = xLeft;
		//m_yOrient = yDown;

		if ( false == bFlipXY )
		{
			q.lb = c3;
			q.lt = c4;
			q.rt = c1;
			q.rb = c2;
		}
		else
		{
			q.lb = c3;
			q.lt = c2;
			q.rt = c1;
			q.rb = c4;
		}

		break;
	default:
		break;
	}
	return q;
}

void MarkerSearcher::SortCoordPairSet()
{
	std::sort(m_cpSet.begin(), m_cpSet.end(), compareCoordPairWCSF);	
}

void MarkerSearcher::SaveCoordPairSet(tstring fName)
{
	std::ofstream out(fName.c_str());
	std::vector<CoordPair>::iterator it = m_cpSet.begin(), it_end = m_cpSet.end();
	for (; it != it_end; ++it) out << it->m.X << " " << it->m.Y << " " << it->m.Z << std::endl;

	//MarkerColumnIt it = m_markerSet.begin(), it_end = m_markerSet.end();
	//for (; it != it_end; ++it) out << it->get_x() << " " << it->get_y() << " " << it->get_d() << std::endl;
}
