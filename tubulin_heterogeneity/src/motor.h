#ifndef _MOTOR_H
#define _MOTOR_H

class motor
{
public:
	
	motor();

	int ID;				// unique identifier of motor 

	int mt_index;		// Index of microtubule motor is on (when bound) 
	int site_coord;	    // Coordinate of motor on microtubule (when bound)

	int motor_index;	// Index of this motor in motor_list

	bool bound;			// Bound or not
	
};
#endif
