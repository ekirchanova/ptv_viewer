#include <vector>

struct Field
{
	Field(int w, int h): width(w), height(h) 
	{
		data.resize(w*h);
	}
	std::vector<BYTE> data;
	int width;
	int height;
};

enum 
{
	_psExist = 1,
	_psAnalyzed = 2
};

void ObtainParticle( Field & f, int indx, std::vector<int>& particleIndxs);
