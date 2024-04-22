// MarkerSearcher.cpp: implementation of the MarkerSearcher class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MarkerSearcher.h"
#include <cmath>
#include <ctime>
#include <numeric>
#include <map>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define MAX_INTENS_VAL 255
#define MARK_INTENS_VAL 128

//#define _DEBUG_

#ifdef _DEBUG_

#include <fstream.h>

#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

int MarkerSearcher::s_id = 0;

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
	float boundaryPart) : m_pCT(&CT)		
{
	m_Z = Z;
	m_xOrient = xOrient;
	m_yOrient = yOrient;
	m_Q = Q;
	m_minDotArea = minDotArea;
	m_maxDotArea = maxDotArea;
	m_minDotCount = minDotCount;
	m_OriginAreaThreshold = OriginAreaThreshold;
	m_boundaryPart = boundaryPart;
	
	//////////////////////////////////////////////////////////////////////////
	m_id = s_id++;
}

MarkerSearcher::~MarkerSearcher()
{

}

int MarkerSearcher::Search(IImageMap *pI)
{	 
	clock_t dt = clock();

	int result = InitTransform(pI);
	ATLASSERT(result);
	if (!result) return 0;

	ATLTRACE("InitTransformTime = %d\n", clock() - dt);	
	dt = clock();

	result = MarkerScan(pI);
	ATLASSERT(result);
	if (!result) return 0;

	ATLTRACE("MarkerScanTime = %d\n", clock() - dt);
	dt = clock();

	result = ClipRectangle();
	ATLASSERT(result);
	if (!result) return 0;
	
	ATLTRACE("ClipRectangleTime = %d\n", clock() - dt);
	dt = clock();

	result = BindImageMarkersToModel();
	ATLASSERT(result);
	if (!result) return 0;

	ATLTRACE("BindImageMarkersToModelTime = %d\n", clock() - dt);
	
	return 1;
}

int MarkerSearcher::InitTransform(IImageMap *pI)
{
	int i, j;
	int w = pI->GetWidth();
	int h = pI->GetHeight();
		
	short* pFFrame = pI->GetFirstFrame();

	int k;
	std::vector<int> intenceGraph;
	intenceGraph.reserve(MAX_INTENS_VAL + 1);
	for (k = 0; k <= MAX_INTENS_VAL; ++k) intenceGraph[k] = 0;

	for (j = 0; j < h; j++)
		for (i = 0; i < w; i++)
			intenceGraph[pFFrame[j * w + i]]++;

	//////////////////////////////////////////////////////////////////////////
	
	/*
	 *	обработка гuстограммы уровня серого	
	 */

//#ifdef _DEBUG_
//
//{
//
//	fstream f;
//
//	f.open("C:\\gray_before.txt", ios::out);
//
//	for (k = 0; k < 256; ++k) {
//
//		f << k << "\t" << intenceGraph[k] <<  endl << flush;
//
//	}	
//
//	f.close();
//
//}
//	
//#endif

	int I1, I2, I3, I4, I5, I6, I7;

	// сглаживание
	for (k = 0; k <= MAX_INTENS_VAL; ++k) {

		I1 = intenceGraph[k];
		I2 = ( k + 1 ) <= MAX_INTENS_VAL ? intenceGraph[k + 1] : 0;
		I3 = ( k + 2 ) <= MAX_INTENS_VAL ? intenceGraph[k + 2] : 0;
		I4 = ( k + 3 ) <= MAX_INTENS_VAL ? intenceGraph[k + 3] : 0;
		I5 = ( k + 4 ) <= MAX_INTENS_VAL ? intenceGraph[k + 4] : 0;
		I6 = ( k + 5 ) <= MAX_INTENS_VAL ? intenceGraph[k + 5] : 0;
		I7 = ( k + 6 ) <= MAX_INTENS_VAL ? intenceGraph[k + 6] : 0;
		
		intenceGraph[k] = int( float( I1 + I2 + I3 + I4 + I5 + I6 + I7 ) / 5 );

	}
	//

	//////////////////////////////////////////////////////////////////////////	

	/*
	 *	поиск локальных максимумов	
	 */	

	std::vector<int> locMaxs;

//	for (k = 0; k <= MAX_INTENS_VAL - 2; ++k) {
//
//		if ( intenceGraph[k] < intenceGraph[k + 1] && intenceGraph[k + 2] < intenceGraph[k + 1] ) {
//			
//			locMaxs.push_back(k + 1);
//			//intenceGraph[k + 1] = 50000;
//
//		}
//
//	}

	// процедура поиска с учётом наличия "горизонтальных плато" в распределениях
	k = 0;
	int k0, k1, k2;

	while ( k <= MAX_INTENS_VAL - 2 ) {

		k0 = k;
		k1 = k + 1;
		k2 = k + 2;

		// если первый максимум в нуле и обрезан (есть только убывающая часть)
		if ( k == 0 && intenceGraph[k] != 0 && intenceGraph[k] >= intenceGraph[k + 1] ) {

			// идём по плато вправо
			while ( intenceGraph[k] == intenceGraph[k + 1] ) ++k;
			//

			// если после плато пошёл опять рост - гоним дальше
			if ( intenceGraph[k] < intenceGraph[k + 1] ) {++k;continue;}

			locMaxs.push_back(k);
			//intenceGraph[k1] = 50000;

		// обработка "нормальных максимумов"
		} else if ( intenceGraph[k0] < intenceGraph[k1] && intenceGraph[k2] <= intenceGraph[k1] ) {

			// идём по плато вправо
			while ( intenceGraph[k + 1] == intenceGraph[k + 2] && k < ( MAX_INTENS_VAL - 2 ) ) ++k;
			//

			k1 = k + 1;
			k2 = k + 2;

			// если после плато пошёл опять рост - гоним дальше
			if ( intenceGraph[k2] > intenceGraph[k1] ) {++k;continue;}

			locMaxs.push_back(k1);
			//intenceGraph[k1] = 50000;

		}

		++k;

	}
	//

	//////////////////////////////////////////////////////////////////////////	

	// распечатка всех максимумов
    std::vector<int>::iterator it, it_end = locMaxs.end();
	ATLTRACE("//////////////////////////////////////////////////////////////////////////\n"); 
	ATLTRACE("local maximums:\n");
	for (it = locMaxs.begin(); it != locMaxs.end(); ++it) ATLTRACE("%d\n", *it);
 	ATLTRACE("//////////////////////////////////////////////////////////////////////////\n");
	//

	if ( locMaxs.size() < 2 ) return 0;

	int iThreshold = ( locMaxs[0] + *locMaxs.rbegin() ) / 2;

	// ЗАЛИПУХА!!!
	//if ( m_id == 0 || m_id == 1 || m_id == 2 ) iThreshold = 105;
	//ATLTRACE("MarkerSearcherID = %d\n", m_id);
	//

	ATLTRACE("iThreshold = %d\n", iThreshold);
	ATLTRACE("//////////////////////////////////////////////////////////////////////////\n");

	//////////////////////////////////////////////////////////////////////////
	// Определяем тип мишени (черные метки или белые метки)
	const bool m_dotsType = m_pCT->get_bgType();

	if ( !m_dotsType ) {

		for (j = 0; j < h; j++)
			for (i = 0; i < w; i++)
				pFFrame[j * w + i] = ( pFFrame[j * w + i] < iThreshold ) ? MAX_INTENS_VAL : 0;

	} else {

		for (j = 0; j < h; j++)
			for (i = 0; i < w; i++)
				pFFrame[j * w + i] = ( pFFrame[j * w + i] < iThreshold ) ? 0 : MAX_INTENS_VAL;

	}

	//////////////////////////////////////////////////////////////////////////

	intenceGraph[iThreshold] = 30000;

	//////////////////////////////////////////////////////////////////////////	

#ifdef _DEBUG_

{

	fstream f;

	f.open(m_id == 0 ? "C:\\Histogram_A.txt" : "C:\\Histogram_B.txt", ios::out);

	for (k = 0; k <= MAX_INTENS_VAL; ++k) {

		f << k << "\t" << intenceGraph[k] <<  endl << flush;

	}	

	f.close();

}

#endif

	return 1;	
}

int MarkerSearcher::MarkerScan(IImageMap *pI)
{
	int i, j;
	int w = pI->GetWidth();
	int h = pI->GetHeight();

	//будем считать, что нам дают однокадровые изображения
	//изображения маркеров, которые находятся на первом кадре
	//второй кадр при этом пуст, т.е. память под него не выделена

	short* pFFrame = pI->GetFirstFrame();
	int index, index2;
	int _i, _j;
	int dotArea;
	float x, y, D;
	short curPos = 0;

	MarkerColumn tmpMarkerSet;
	
	// множество для пикселов области связности
	std::vector< Position > pixSet;
	pixSet.reserve(3000);//Примерное количество пикселей в маркере на изображении ~2010 (взято с запасом)
	//

	// очищаем поле маркеров
	MarkerFieldIt fit, fit_end = m_markerSet.end();
	for (fit = m_markerSet.begin(); fit != fit_end; ++fit) fit->clear();
	m_markerSet.clear();	
	//	

	//////////////////////////////////////////////////////////////////////////	

	for (i = 0; i < w; i++) {
	
		for (j = 0; j < h; j++) {

			index = j * w + i;
			
			//поиск затравки
			//т.е. начала связной области
			if ( pFFrame[index] == MAX_INTENS_VAL ) {				
				
				curPos = 0;

				pixSet.push_back(Position(i, j));

				//////////////////////////////////////////////////////////////////////////				

				while ( curPos < pixSet.size() ) {

					_i = pixSet[curPos].m_i;
					_j = pixSet[curPos].m_j;
					
					//вверх
					if ( _j + 1 < h ) {
						
						index2 = (_j + 1) * w + _i;

						if ( pFFrame[index2] == MAX_INTENS_VAL ) {
							
							pixSet.push_back(Position(_i, _j + 1));
							pFFrame[index2] = MARK_INTENS_VAL;
							
						}

					}
					//

					//вниз
					if ( _j - 1 >= 0 ) {
						
						index2 = (_j - 1) * w + _i;

						if ( pFFrame[index2] == MAX_INTENS_VAL ) {

							pixSet.push_back(Position(_i, _j - 1));
							pFFrame[index2] = MARK_INTENS_VAL;

						}

					}
					//

					//вправо
					if ( _i + 1 < w ) {

						index2 = _j * w + _i + 1;

						if ( pFFrame[index2] == MAX_INTENS_VAL ) {
							
							pixSet.push_back(Position(_i + 1, _j));
							pFFrame[index2] = MARK_INTENS_VAL;

						}

					}
					//

					//влево
					if ( _i - 1 >= 0 ) {

						index2 = _j * w + _i - 1;

						if ( pFFrame[index2] == MAX_INTENS_VAL ) {
							
							pixSet.push_back(Position(_i - 1, _j));
							pFFrame[index2] = MARK_INTENS_VAL;
						
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
				
				// площадь маркера должна быть больше минимально возможной
				if ( dotArea > m_minDotArea && dotArea < m_maxDotArea ) {
				
					//Вычисляет характерный диаметр маркера
					//Обратное действие от S = PI * D^2 / 4 будет D = sqrt(4 * S / PI)
					D = (float)sqrt(4 * dotArea / 3.14);
					//ATLTRACE("D = %g\n", D);

					//на изображении 035nofz5.img из-за несовершенства алгоритма бинаризации и
					//алгоритма поиска связных областей (вверх, вниз, вправо, влево) к
					//большой связной области примыкает один пиксел по-диаганали
					//(см. файл CalibrationProc\bugs!\bug1.bmp),
					//который и определяется рядом с большой связной областью
					//всё это приводит к тому, что расстояни м/у ними оказывается равным 0
					//и алгоритм ломается
					//Вообще, надо переделывать "Бинаризацию" или вводить в "Поиск связных областей"
					//шаги по-диаганали (что увеличит время расчёта)
					
					tmpMarkerSet.push_back(Marker(x, y, D));
					
				}
				//
				
				pixSet.clear();
			
			}

		}

	}

	//////////////////////////////////////////////////////////////////////////

	/*
	 *	Валидация маркеров
	 */

	float meanArea = std::accumulate(tmpMarkerSet.begin(), tmpMarkerSet.end(), 0.0f, sumAreaF()) / tmpMarkerSet.size();

	ATLTRACE("meanArea = %g\n", meanArea);

	tmpMarkerSet.erase(std::remove_if(tmpMarkerSet.begin(), tmpMarkerSet.end(), 
		areaValidMarkerF(meanArea, m_Q)), tmpMarkerSet.end());

	//ATLTRACE("ok\n");

	// количество маркеров больше 25 = 5 * 5 (для возможности построения отображения)
	if ( tmpMarkerSet.size() < m_minDotCount ) {ATLTRACE("tmpMarkerSet.size() < m_minDotCount failed...\n");return 0;}
	
	//////////////////////////////////////////////////////////////////////////	

	// сортировка маркеров
	std::sort(tmpMarkerSet.begin(), tmpMarkerSet.end(), compareMarkerPositionF(h));

	//////////////////////////////////////////////////////////////////////////
	
	/*
	 *	Отсеиваем маркеры расположенные близко от границы и находим центральный маркер начала координат
	 */

	MarkerColumnIt cit, cit_end = tmpMarkerSet.end();
	MarkerColumn tmp;
	float R, prev_x;
	bool isFirst = true;

	for(cit = tmpMarkerSet.begin(); cit != cit_end; ++cit) {

		Marker &m = *cit;

		x = m.get_x();
		y = m.get_y();
		R = m.get_d() / 2;//Здесь R - критерий отсева маркеров на границе

		// первый раз prev_x делаем равным x
		if (isFirst) {
			prev_x = x;
			isFirst = false;
		}		

		//вставляем вертикальную линию в выходной массив маркеров
		if ( x - prev_x > R && !tmp.empty() ) {
			
			std::sort(tmp.begin(), tmp.end(), compareMarkerYF());
			m_markerSet.push_back(tmp);
			tmp.clear();
		
			//ATLTRACE("line added...\n");

		}			
		//	
		
		//Составляем вертикальную линию с отсевом маркеров на границе, т.к. они
		//могут быть обрезаны и их  положение может быть неточным
		if ( x - R > w * m_boundaryPart &&
			 y - R > h * m_boundaryPart &&
			 x + R < w * ( 1 - m_boundaryPart ) &&
			 y + R < h * ( 1 - m_boundaryPart ) ) {

			tmp.push_back(m);
			
			//ATLTRACE("element added...\n");

		}		

		prev_x = x;

	}
	
	if ( !tmp.empty() ) {

		std::sort(tmp.begin(), tmp.end(), compareMarkerYF());
		m_markerSet.push_back(tmp);

	}

	/////////////////////////////////////////////////////////////////////////	
	
	// ищем центральный маркер
	findMaxDMarkerF find;
	find = std::for_each(m_markerSet.begin(), m_markerSet.end(), find);
	m_OriginMarker = find.get_OriginMarker();
	m_Origin = find.get_Origin();

	// проверяем удовлетворяет ли центральный маркер
	// критерию по величине площади
	if ( ( 3.14 * m_OriginMarker.get_d() * m_OriginMarker.get_d() / 4 ) < meanArea * m_OriginAreaThreshold ) {
		ATLTRACE("( 3.14 * m_OriginMarker.get_d() * m_OriginMarker.get_d() / 4 ) < meanArea * m_OriginAreaThreshold failed...\n");
		return 0;
	}

	ATLTRACE("Center marker: x = %d\ty = %d\n", m_Origin.m_i, m_Origin.m_j);

	//////////////////////////////////////////////////////////////////////////	

	return 1;
}

int MarkerSearcher::ClipRectangle()
{
	// находим уравнения прямых для ниней и верхней норизонтальной линии из маркеров
	float x1, x2, y1, y2, k_down, b_down, k_up, b_up;
	int h_size;	
	Marker m;
	MarkerFieldIt minVertLineIt, nearMinVertLineIt;
	MarkerColumnIt eqElemIt;
	//	
	
	// поле из маркеров должно иметь не менее 4 вертикальных линий
	// это необходимо для правильной работы алгоритма клипирования
	// ( см. ++m_markerSet.begin() и --m_markerSet.end() )
	if ( m_markerSet.size() < 4 ) return 0;

	MarkerFieldIt afterFirstIt = m_markerSet.begin();afterFirstIt++;
	MarkerFieldIt beforeEndIt = m_markerSet.end();beforeEndIt--;
	
	// обрезка снизу
	minVertLineIt = std::min_element(afterFirstIt, beforeEndIt, 
		sizeEqualVertMarkerLineF<yLessMarkerF>(yLessMarkerF(m_OriginMarker.get_y() - m_OriginMarker.get_d() / 2)));

	ATLTRACE("horIndexDown: %d\n", std::distance(m_markerSet.begin(), minVertLineIt));
	//ATLTRACE("SizeDown: %d\n", minVertLineIt->size());	

	// шаг вперёд к следующему вертикальному столбику
	nearMinVertLineIt = minVertLineIt + 1;

	eqElemIt = std::find_if(nearMinVertLineIt->begin(), nearMinVertLineIt->end(), std::bind2nd(yEqualMarkerF(), (*minVertLineIt)[0]));

	if ( eqElemIt == nearMinVertLineIt->end() ) {
		// эквивалентный маркер не найден

		// шаг назад к предыдущему вертикальному столбику
		nearMinVertLineIt = minVertLineIt - 1;

		eqElemIt = std::find_if(nearMinVertLineIt->begin(), nearMinVertLineIt->end(), std::bind2nd(yEqualMarkerF(), (*minVertLineIt)[0]));

		if ( eqElemIt == nearMinVertLineIt->end() ) {
			// эквивалентный маркер не найден
			// продолжение невозможно
			return 0;
		}

	}
	
	ATLTRACE("vertIndexDown: %d\n", std::distance(nearMinVertLineIt->begin(), eqElemIt));
	
	//
	m = *eqElemIt;
	x1 = m.get_x();
	y1 = m.get_y();
	m = (*minVertLineIt)[0];
	x2 = m.get_x();
	y2 = m.get_y();
	k_down = ( y2 - y1 ) / ( x2 - x1 );
	b_down = ( y1 * x2  - y2 * x1 ) / ( x2 - x1 );
	//

	//Обрезка сверху
	minVertLineIt = std::min_element(afterFirstIt, beforeEndIt,
		sizeEqualVertMarkerLineF<yGreaterMarkerF>(yGreaterMarkerF(m_OriginMarker.get_y() + m_OriginMarker.get_d() / 2)));

	ATLTRACE("horIndexUp: %d\n", std::distance(m_markerSet.begin(), minVertLineIt));
	//ATLTRACE("SizeUp: %d\n", minVertLineIt->size());	

	// шаг вперёд к следующему вертикальному столбику
	nearMinVertLineIt = minVertLineIt + 1;
	h_size = minVertLineIt->size();

	eqElemIt = std::find_if(nearMinVertLineIt->begin(), nearMinVertLineIt->end(), std::bind2nd(yEqualMarkerF(), (*minVertLineIt)[h_size - 1]));

	if ( eqElemIt == nearMinVertLineIt->end() ) {
		// эквивалентный маркер не найден

		// шаг назад к предыдущему вертикальному столбику
		nearMinVertLineIt = minVertLineIt - 1;

		eqElemIt = std::find_if(nearMinVertLineIt->begin(), nearMinVertLineIt->end(), std::bind2nd(yEqualMarkerF(), (*minVertLineIt)[h_size - 1]));

		if ( eqElemIt == nearMinVertLineIt->end() ) {
			// эквивалентный маркер не найден	
			// продолжение невозможно
			return 0;
		}

	}

	ATLTRACE("vertIndexUp: %d\n", std::distance(nearMinVertLineIt->begin(), eqElemIt));

	//
	m = *eqElemIt;
	x1 = m.get_x();
	y1 = m.get_y();
	m = (*minVertLineIt)[h_size - 1];
	x2 = m.get_x();
	y2 = m.get_y();
	k_up = ( y2 - y1 ) / ( x2 - x1 );
	b_up = ( y1 * x2  - y2 * x1 ) / ( x2 - x1 );
	//

//	int i, j;
//
//	//рисуем ограничивающие прямые
//	for (i = 0; i < w; i++) {
//
//		j = k_down * i + b_down;
//		if ( j > 0 && j < h) pFFrame[w * j + i] = 255;
//
//		j = k_up * i + b_up;
//		if ( j > 0 && j < h) pFFrame[w * j + i] = 255;
//
//	}	

	// удаляем маркеры в столбцах по ограничивающим прямым
	std::for_each(m_markerSet.begin(), m_markerSet.end(),
		condDeleteMarkersF<rangeDeleteMarkersF>(rangeDeleteMarkersF(k_down, b_down, k_up, b_up)));

//	return 1;
	
	// удаляем вертикальные линии, в которых маркеров меньше, чем максимальное количество маркеров среди вертикальных столбцов
	m_markerSet.erase(
		std::remove_if(m_markerSet.begin(), m_markerSet.end(), 
			std::not1(std::bind2nd(equalSizeVertLineF(),
				*std::max_element(m_markerSet.begin(), m_markerSet.end(), compareSizeVertMarkerLineF())))),
		m_markerSet.end());

	// заново ищем центральный маркер
	findMaxDMarkerF find;
	find = std::for_each(m_markerSet.begin(), m_markerSet.end(), find);
	m_OriginMarker = find.get_OriginMarker();
	m_Origin = find.get_Origin();

	ATLTRACE("After clipping -> Center marker: x = %d\ty = %d\n", m_Origin.m_i, m_Origin.m_j);

	return 1;
}

int MarkerSearcher::BindImageMarkersToModel()
{
	// проверка и выход если выполнение невозможно
	if ( m_pCT == 0 || m_Origin.m_i < 1 || m_Origin.m_j < 1 )  {
		// m_Origin.m_i < 1 || m_Origin.m_j < 1 - нужно чтобы считать
		// типичое расстояние между маркерами на изображении
		return 0;
	}

	//
	m_cpSet.clear();
	m_r = ModelRegion();	

	const float dotSpacing = m_pCT->get_dotSpacing();
	const float dotSpacingHalf = dotSpacing / 2;

	//ATLTRACE("%g %g\n", dotSpacing, dotSpacingHalf);

	float x, y, X, Y, Z;
	//

	//
	MarkerFieldIt leftFromOriginLineIt;
	MarkerColumnIt leftFromOriginMarkerIt;
	MarkerColumnIt equalOriginMarkerIt;
	//
	
	//
	int i, j, j2;
	MarkerFieldIt hit, hit_begin = m_markerSet.begin(), hit_end = m_markerSet.end();		
	//

	//int counter = 1;

	ATLTRACE("xOrient =  %d, yOrient = %d\n", m_xOrient, m_yOrient);
	
	for (hit = hit_begin; hit != hit_end; ++hit) {

		i = std::distance(hit_begin, hit);

		MarkerColumn &vertLine = *hit;
		MarkerColumnIt vit, vit_begin = vertLine.begin(), vit_end = vertLine.end();

		for (vit = vit_begin; vit != vit_end; ++vit) {

			j = std::distance(vit_begin, vit);

			Marker &m = *vit;

			x = m.get_x();
			y = m.get_y();
			
			equalOriginMarkerIt = std::find_if(hit->begin(), hit->end(), 
				std::bind2nd(yEqualMarkerF(), m_OriginMarker));

			if ( equalOriginMarkerIt == hit->end() ) return 0;

 			j2 = std::distance(hit->begin(), equalOriginMarkerIt);

			X = ( i - m_Origin.m_i ) * dotSpacing * m_xOrient;
			Y = ( j - j2 ) * dotSpacing * m_yOrient;
			Z = m_Z;			

			//ATLTRACE("%d: (%d, %d)\n", counter, i - m_Origin.m_i, j - j2);
			//counter++;

			m_cpSet.push_back(CoordPair(ImageCoord(x, y), ModelCoord(X, Y, Z)));
			
			//////////////////////////////////////////////////////////////////////////
			//получаем прямоугольник (или объем т.к. есть Z-координата), описывающий область из маркеров в модельной системе координат			

			// ищем левую нижнюю			
			if ( X < m_r.lb.X ) m_r.lb.X = X;
			if ( Y < m_r.lb.Y ) m_r.lb.Y = Y;			
			
			// ищем правую верхнюю
			if ( X > m_r.rt.X ) m_r.rt.X = X;
			if ( Y > m_r.rt.Y ) m_r.rt.Y = Y;

		}

	}

	m_r.lb.X -= dotSpacingHalf;
	m_r.lb.Y -= dotSpacingHalf;
	m_r.rt.X += dotSpacingHalf;
	m_r.rt.Y += dotSpacingHalf;
	m_r.lb.Z = Z;
	m_r.rt.Z = Z;

	return 1;
}

ModelRegion MarkerSearcher::CommonModelRegion(const ModelRegion &mr1, const ModelRegion &mr2)
{
	return
		ModelRegion(ModelCoord(
			mr1.lb.X < mr2.lb.X ? mr2.lb.X : mr1.lb.X,
			mr1.lb.Y < mr2.lb.Y ? mr2.lb.Y : mr1.lb.Y,
			mr1.lb.Z < mr2.lb.Z ? mr2.lb.Z : mr1.lb.Z),
		ModelCoord(
			mr1.rt.X < mr2.rt.X ? mr1.rt.X : mr2.rt.X,
			mr1.rt.Y < mr2.rt.Y ? mr1.rt.Y : mr2.rt.Y,
			mr1.rt.Z < mr2.rt.Z ? mr1.rt.Z : mr2.rt.Z));
}

int MarkerSearcher::Intersection(MarkerSearcher &leftMS, MarkerSearcher &rightMS)
{
	MarkerField &rms = rightMS.get_MarkerSet();
	MarkerField &lms = leftMS.get_MarkerSet();

	int l_x_index = leftMS.get_Origin().m_i;
	int l_y_index = leftMS.get_Origin().m_j;

	int r_x_index = rightMS.get_Origin().m_i;
	int r_y_index = rightMS.get_Origin().m_j;

	int l_h_size = lms.size();
	int l_v_size = lms[0].size();

	int r_h_size = rms.size();
	int r_v_size = rms[0].size();
	
	int i;

	//выравнивавем по-вертикали снизу
	if ( l_y_index < r_y_index ) {

		//подрезаем правую

		//ATLTRACE("vert_down_right = %d\n", r_y_size - l_y_size);

		for (i = 0; i < r_h_size; i++) {

			rms[i].erase(rms[i].begin(), rms[i].begin() + r_y_index - l_y_index);

		}

		r_y_index = l_y_index;
		
	} else if ( l_y_index > r_y_index ) {

		//подрезаем левую

		//ATLTRACE("vert_down_left = %d\n", l_y_size - r_y_size);

		for (i = 0; i < l_h_size; i++) {			

			lms[i].erase(lms[i].begin(), lms[i].begin() + l_y_index - r_y_index);

		}

		l_y_index = r_y_index;

	}

	l_v_size = lms[0].size();
	r_v_size = rms[0].size();

	//выравнивавем по-вертикали сверху
	if ( l_v_size < r_v_size ) {

		//подрезаем правую

		ATLTRACE("vert_up_right = %d\n", r_v_size - l_v_size);

		for (i = 0; i < r_h_size; i++) {

			rms[i].erase(rms[i].begin() + l_v_size, rms[i].end());

		}

	} else if ( l_v_size > r_v_size ) {

		//подрезаем левую

		ATLTRACE("vert_up_left = %d\n", l_v_size - r_v_size);
		
		for (i = 0; i < l_h_size; i++) {

			lms[i].erase(lms[i].begin() + r_v_size, lms[i].end());

		}

	}

	//выравниваем по-горизонтали снизу
	if ( l_x_index < r_x_index ) {

		//подрезаем правую

		ATLTRACE("hor_down_right = %d\n", r_x_index - l_x_index);

		rms.erase(rms.begin(), rms.begin() + r_x_index - l_x_index);

		r_x_index = l_x_index;

	} else if ( l_x_index > r_x_index ) {

		//подрезаем левую

		ATLTRACE("hor_down_left = %d\n", r_x_index - r_x_index);

		lms.erase(lms.begin(), lms.begin() + l_x_index - r_x_index);

		l_x_index = r_x_index;

	}

	l_h_size = lms.size();
	r_h_size = rms.size();

	//выравниваем по-горизонтали сверху
	if ( l_h_size < r_h_size ) {

		//подрезаем правую

		ATLTRACE("hor_up_right = %d\n", r_h_size - l_h_size);

		rms.erase(rms.begin() + l_h_size, rms.end());

	} else if ( l_h_size > r_h_size ) {

		//подрезаем левую

		ATLTRACE("hor_up_left = %d\n", l_h_size - r_h_size);

		lms.erase(lms.begin() + r_h_size, lms.end());

	}

	leftMS.m_Origin = Position(l_x_index, l_y_index);
	rightMS.m_Origin = Position(r_x_index, r_y_index);

	return 1;
}
