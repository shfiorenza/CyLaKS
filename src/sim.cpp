#include "master_header.h"

int main(int argc, char *argv[]){

    int threads_provided, world_rank, world_size;
	system_parameters parameters;
	system_properties properties;

	// Initialize MPI with capacity for multiple threads per node
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &threads_provided);
	// Check that MPI node initiated its threads properly
	properties.wallace.CheckMPI(threads_provided);
	// Get world rank and size
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Check that input has the correct number of arguments
	properties.wallace.CheckArguments(argv[0], argc);
	// Parse parameters from yaml file into parameter structure
	properties.wallace.ParseParameters(&parameters, argv[1]);
	// Initialize sim objects (MTs, kinesin, MAPs, etc.)
	properties.wallace.InitializeSimulation(&properties);
	// Generate data files (on root node only)
	if(world_rank == 0) properties.wallace.GenerateDataFiles(argv[2]);
	// Synchronize MPI nodes (if necessary)
	if(world_size > 1) MPI_Barrier(MPI_COMM_WORLD);

	// Main KMC loop
	for(int i_step = 0; i_step < parameters.n_steps; i_step++){
		// Synchronize MPI nodes (if necessary)
		if(world_size > 1 ) MPI_Barrier(MPI_COMM_WORLD);
		properties.wallace.UpdateTimestep(i_step);
		properties.kinesin4.RunKMC();
		properties.prc1.RunKMC();
		properties.prc1.RunDiffusion();
		properties.microtubules.RunDiffusion();
	}

	// Cleanup stuff
	if(world_rank == 0) properties.wallace.CloseDataFiles();
	properties.wallace.OutputSimDuration();
	properties.gsl.CleanUp();
	MPI_Finalize(); 
	return 0;
}
