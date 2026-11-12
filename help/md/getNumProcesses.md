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
    
    # Record the process count for reproducibility bookkeeping
    print("Running on " + getNumProcesses() + " process(es).")
    
    # Guard-check the launch configuration
    if (getNumProcesses() != 32) {
        stop("Expected 32 processes, got " + getNumProcesses() + ".")
    }
## references
