#ifndef _NEIGHBOR_ENTRY_H
#define _NEIGHBOR_ENTRY_H

class Tubulin;
class AssociatedProtein;
class Kinesin;

struct neighbor_entry{

	int speciesID_;							// ID as follows:
	Tubulin *site_ = nullptr;				// 0 for tubulin
	AssociatedProtein *xlink_ = nullptr;		// 1 for xlink
	Kinesin *motor_ = nullptr;				// 2 for motor
	int distance_ = 0;						// distance in sites
	double weight_ = 0;						// statistical weight

}; 
#endif
