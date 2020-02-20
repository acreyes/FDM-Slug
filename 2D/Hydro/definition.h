
! slope limiters
#define MINMOD  1
#define MC      2
#define VANLEER 3

! Riemann solvers
#define HLL     1
#define HLLC    2
#define ROE     3

! MultiD vars
#define XDIM 1
#define YDIM 2
#define NDIM 2

! primitive vars
#define DENS_VAR 1
#define VELX_VAR 2
#define VELY_VAR 3
#define PRES_VAR 4
#define BDRY_VAR 5
#define EINT_VAR 6
#define GAMC_VAR 7
#define GAME_VAR 8
#define NUMB_VAR 8
#define NSYS_VAR 4 /* total number of equations of the conservative system */
/*#define BDRY_VAR 8 /* an extra boundary var to simulate an internal solid boundary structure */

! conservative vars
#define MOMX_VAR 2
#define MOMY_VAR 3
#define ENER_VAR 4

! waves
#define SHOCKLEFT 1
#define SLOWWLEFT 2
#define CTENTROPY 3
#define SHOCKRGHT 4
#define NUMB_WAVE 4

! setup parameters
#define MAX_STRING_LENGTH 80

! BC
#define PERIODIC -1
#define OUTFLOW -2
#define DBLMCH -3
#define OUTFLOW45 -4
#define INFLOW -5
#define REFLECT -6

#define PI 4.*ATAN(1.)
