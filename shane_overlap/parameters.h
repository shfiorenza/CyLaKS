/* The <system_parameters> structure contains thermodynamic, control, and
   configuration parameters for the simulation. */

typedef struct {

   int n_steps;				/* number of simulation runs */

   int pickup_time;	 		/* how frequent to output data */

   int n_protofilaments;    /* number of individual microtubules */

   long seed;				/* random number seed */

   double v_motor_g;		/* motor velocity um/min */

   double f_turning;		/* motor turning frequency 1/min */

   double c_motor;          /* bulk motor concentration nM */

   double kon;              /* motor attachment frequency. 1/(nM*s) */

   double koff;				/* motor dissociation frequency 1/min */

} system_parameters;
