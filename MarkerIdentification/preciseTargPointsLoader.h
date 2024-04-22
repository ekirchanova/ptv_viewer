class PreciseTargPointsLoader
{
public:
	PreciseTargPointsLoader(std::string fName) 
	{
		float X, Y;
		std::ifstream in(fName.c_str());
		in >> w; in >> h;		
		for (int i = 0; i < w * h; ++i)
		{
			in >> X; in >> Y;
			uc.push_back(ModelCoord(X, Y, 0.f));
			//std::cout << X << " " << Y << "\n";
		}
	}

	ModelCoord get(int i, int j)
	{
		int ii = i + int(w / 2);
		int jj = j + int(h / 2);
		return uc[w * jj + ii];
	}

private:
	std::vector<ModelCoord> uc;
	int w, h;

};