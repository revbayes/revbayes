## name
getNumProcesses
## title
Get the number of MPI processes
## description
Returns the number of processes used by the current RevBayes session.
## details
In an MPI-enabled build (i.e., `rb-mpi`, launched as `mpirun -np N ./rb-mpi`),
the function returns `N` (i.e., the size of the `MPI_COMM_WORLD` communicator).
In a non-MPI build, it always returns 1. The value is the number of MPI ranks,
not necessarily the number of hardware CPU cores.
## authors
David Černý
## see_also
## example
    # For `mpirun -np 16 ./rb-mpi`, this prints 16:
    getNumProcesses()
    
    # Set the number of MCMC replicates to the number of available processes,
    # but at least to 4:
    N_RUNS = max( [getNumProcesses(), 4] )
    
    # Create a simple model (unclamped), and run the number of replicates
    # specified above:
    x ~ dnExp(10)
    mymcmc = mcmc( model(x), [mvSlide(x, delta=0.1)], [mnScreen()], nruns=N_RUNS )
    mymcmc.run(100)
## references
