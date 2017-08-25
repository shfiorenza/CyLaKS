#ifndef _ASSOCIATED_PROTEIN_H
#define _ASSOCIATED_PROTEIN_H

#ifndef _PARAMETERS_H
typedef struct system_parameters system_parameters;
#endif
#ifndef _SYSTEM_PROPETIES_H
typedef struct system_properties system_properties;
#endif

class Microtubule;
class Tubulin;
class Kinesin;

class AssociatedProtein{
	private:

	public:
		int ID_;

		bool bound_ = false;
		bool tethered_ = false;

		Microtubule *mt_ = nullptr;
		Tubulin *site_ = nullptr;
		Kinesin *motor_ = nullptr;

		system_parameters *parameters_ = nullptr;
		system_properties *properties_ = nullptr;
	private:

	public:
		AssociatedProtein();
};
#endif
