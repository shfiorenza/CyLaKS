#ifndef _KINESIN_H
#define _KINESIN_H

#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPERTIES_H
typedef struct system_properties system_properties;
#endif

class Microtubule;
class Tubulin;
class AssociatedProtein;

class Kinesin{
	private:

	public:
		int ID_;						// Unique ID of this kinesin
		int speciesID_ = 2;				// Unique ID of this species (kinesin)

		bool bound_ = false; 			// Is this motor bound to a microtubule? (i.e. are
										// the heads of this motor bound to tubulin sites?) 
		bool tethered_ = false; 		// Is the tail of this motor bound to a prc1 protein?

		Microtubule *mt_ = nullptr; 	// Microtubule this kinesin is bound to
		Tubulin *site_ = nullptr; 		// Site this kinesin is bound to
		AssociatedProtein *anchor_ = nullptr; 	// PRC1 protein this kinesin is tethered to

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		Kinesin();
		void Initialize(system_parameters *parameters, 
			system_properties *properties, int ID);
};
#endif
