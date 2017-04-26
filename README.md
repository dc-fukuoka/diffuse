*diffuse* - 3D diffusion equation solver by using CG method(no precondition) with MPI.
===============	  
the result was compared to serial version(cg.F90)
hybrid version is still under testing.
===============	  
input file: fort.11

    $ cat fort.11 
    &param0
    dt        = 1.0d-6  ! stride of time
    dx        = 1.0d-2  ! stride of x, y, z direction
    tol       = 1.0d-40 ! convergence tolerance of CG method
    diff_coef = 1.0d2   ! diffusion coefficient
    /
    &param1
    imax = 64 ! global size for x direction
    jmax = 64 ! global size for y direction
    kmax = 64 ! global size for z direction
    /
    &param2
    idiv = 1 ! # of division for x direction
    jdiv = 1 ! # of division for y direction
    kdiv = 4 ! # of division for z direction
    ! idiv*jdiv*kdiv must equal to total # of processes(np)
    ! imax/idiv, jmax/jdiv and kmax/kdiv must be an integer
    /
    
    &param3
    iter_max   = 262144 ! maximum CG method iteration
    tstep_max  = 100    ! maximum timestep
    freq_write = 100    ! frequency of writing the result, write each tstep_max/freq_write time
    /


how to run:
    
    $ make  
    $ vi fort.11  
    $ ./create_input  
    $ mpirun -np $NP ./diffuse  
    $ ./create_output  
    $ ./create_anime.sh  

cpu: Intel(R) Xeon(R) CPU E5-2450 0 @ 2.10GHz  
problem size: 64x64x64  
maximum timestep: 100  
with 1   core: 38.3777060508728 s  
with 32 cores:  2.2367570400238 s  

the animation shows that how a gaussian wave decays at k = kmax/2
![Alt text](./diffuse.gif?raw=true "diffuse.gif")
