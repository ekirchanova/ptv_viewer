#include "stdafx.h"
#include "PlanarContainer.h"

bool operator < (const Particle& p1, const Particle& p2)
{
	if ( ( p1.x < p2.x ) && (p1.y < p2.y) )
	{
		return true;
	}
	return false;
}