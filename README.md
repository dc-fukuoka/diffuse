*diffuse* - 3D diffusion equation solver by using CG method(no precondition) + Crank-Nicolson method with MPI.
===============	  
the result was compared to serial version(cg.F90)  
  
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
    
    $ make # if you have intel compiler and intelmpi, try "make -f makefile.impi.mk" and if you want to disable Crank-Nicolson method, try to remove -D_CN from the makefile
    $ vi fort.11 # adjust the parameters  
    $ ./create_input  
    $ mpirun -np $NP ./diffuse # or mpirun -np $NP ./diffuse_hyb where $NP = idiv*jdiv*kdiv  
    $ ./create_output  
    $ ./create_anime.sh  
  
performance comparison
===============	  
cpu: Intel(R) Xeon(R) CPU E5-2450 0 @ 2.10GHz  
problem size: 128x128x128  
tstep_max   : 100  
freq_write  : 100  
  
with 1   core           : 239.6 s  
with 32 cores(flat MPI) : 17.6  s(division:2x4x4)  
with 32 cores(hybrid)   : 14.4  s(division:1x2x2)  

compare the results between implicit method(CG method) and explict method
===============	  

compare the results between implicit method and explicit method.
explicit method program is expl.F90.  
the animation shows that how a gaussian wave decays at k = kmax/2
# 1.

    dt        = 1.0d-6
    dx        = 1.0d-2
    tol       = 1.0d-40
    diff_coef = 1.0d2
    tstep_max  = 100
    freq_write = 100
  
stability condition of explicit method:  
3\*diff_coef\*dt/dx/dx <= 1/2  
= 3.0  
this is instable for explicit method.  
    
# implicit method
`$ eog gifs/diffuse.impl.1.gif`  
![Alt text](gifs/diffuse.impl.1.gif?raw=true "implicit method 1")
# explicit method
the result diverges due to the instability of explicit method...  
`$ eog gifs/diffuse.expl.1.gif`  
![Alt text](gifs/diffuse.expl.1.gif?raw=true "explicit method 1")

# 2.

    dt        = 1.0d-7
    dx        = 1.0d-2
    tol       = 1.0d-40
    diff_coef = 1.0d2
    tstep_max  = 1000
    freq_write = 100
  
stability condition of explicit method:  
3\*diff_coef\*dt/dx/dx <= 1/2  
= 0.3  
this is stable for explicit method.  
    
# implicit method
`$ eog gifs/diffuse.impl.2.gif`  
![Alt text](gifs/diffuse.impl.2.gif?raw=true "implicit method 2")
# explicit method
this is stable condition, in this case, explicit method converges.  
`$ eog gifs/diffuse.expl.2.gif`  
![Alt text](gifs/diffuse.expl.2.gif?raw=true "explicit method 2")
  
note: calculated diffusion speed is a little bit different between implicit method and explicit method, because implicit method case uses Crank-Nicolson method as well, that is time direction second orcer accuracy so implicit method case is more accurete.  
if Crank-Nicolson method is disabled, (remove -D_CN from the makefile) the calculated diffusion speed will be almost the same.
