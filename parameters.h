/* The <system_parameters> structure contains thermodynamic, control, and
   configuration parameters for the simulation. */

typedef struct {

   int n_steps;			/* number of simulation runs */

   int pickup_time;		/* how frequent to output data */

   int fc_flag;			/* the dependence of catastrophe frequency, 1 means density-controlled,
                                   2 means flus-controlled, others mean constant frequency. */

   long seed;			/* random number seed */

   double vg;			/* growing velocity um/min */

   double shrink;		/* shrinking velocity um/min */

   double v_motor_g;		/* motor velocity um/min */

   double f_turning;		/* motor turning frequency 1/min */

   double c_motor;              /* bulk motor concentration nM */

   double kon;                  /* motor attachment frequency. 1/(nM*s) */

   double c_kon;		/* motor attachment frequency times bulk motor concentration */

   double koff;			/* motor dissociation frequency 1/min */

   double koff_end;		/* motor dissociation frequency at the end site 1/min */

   double right_conc_right;	/* right motor conc. at right end. */

   double right_conc_left;	/* left motor conc. at right end. */

   double left_conc_right;	/* right motor conc. at left end. */

   double left_conc_left;	/* left motor conc. at left end. */

   int n_protofilaments;

   int flux_flag;		/* 1 means flux control, and 0 means density control. */

   int *pseudo_end;		/* pseudo_end[0] left, pseudo_end[1] right. 0 means no motor and 1 means a motor. */

   int *temp_end;

} system_parameters;
