class motor
{
public:
	
	motor();

	int ID;				// unique identifier of motor 

	int mt_index;		// Index of microtubule motor is on (when bound) 
	int site_coord;	    // Coordinate of motor on microtubule (when bound)
	int global_coord;	// Encodes both mt and site information

	int motor_entry;	// Entry of motor_list that corresponds to motor (when bound)

	bool bound;			// Bound or not
	
};
