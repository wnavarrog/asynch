ASYNCH is an efficient numerical solver for solving large systems of ordinary differential equations (ODEs) with a tree structure. Source files are written in the C programming language.

Features include:

-Designed for a distributed memory architecture for fast and efficient calculations. This is done using the Message Passing Interface (MPI).

-Communication between processes is performed asynchronously.

-The numerical integrators are asynchronous, i.e. different integrators may be used at each ODE with different step sizes.