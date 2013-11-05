mesh = protein.vol

shared = /home/parallels/git/mfes/nglib/pointcharges
shared = /home/parallels/git/mfes/nglib/energydiff
shared = /home/parallels/git/mfes/nglib/writePotatAscii

define constant eps0 = 8.8541878e-22
define constant q0 = 1.60217646e-19

define constant heapsize = 2000000

define coefficient epsilon_solv
(80*eps0),(4*eps0)

define coefficient epsilon_ref
(4*eps0),(4*eps0)

define fespace v -order=4 -type=h1ho -dirichlet=[1]

define gridfunction u_solv -fespace=v
define gridfunction u_ref -fespace=v

define bilinearform a_solv -fespace=v -symmetric
laplace epsilon_solv

define bilinearform a_ref -fespace=v -symmetric
laplace epsilon_ref

define linearform f -fespace=v

numproc pointcharges ps1 -linearform=f -pqrfile=born.pqr -interpolate

define preconditioner c -type=multigrid -bilinearform=a_solv -inverse=mumps
define preconditioner c0 -type=multigrid -bilinearform=a_ref -inverse=mumps

numproc bvp np1 -gridfunction=u_solv -bilinearform=a_solv -linearform=f -preconditioner=c  -maxsteps=10
numproc bvp np10 -gridfunction=u_ref -bilinearform=a_ref -linearform=f -preconditioner=c0  -maxsteps=10

numproc energydiff npeval -gridfunction=u_solv -gridfunction0=u_ref -pqrfile=born.pqr
