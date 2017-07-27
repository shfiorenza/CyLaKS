#include <vector>

class motor; 

class microtubule
{
public:
	
	microtubule();

	int polarity;
	int plus_end;
	int minus_end;
	int delta_x;
	int mt_index_adj;

	int length;
	int n_bound;
	int coord;
	
	std::vector<motor*> track;

private:
};
