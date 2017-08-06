#ifndef _TUBULIN_H
#define _TUBULIN_H

class motor;

class tubulin
{
public:

	tubulin();

	bool mutant;

	int mt_index;
	int site_coord;

	motor *occupant; 
};
#endif
