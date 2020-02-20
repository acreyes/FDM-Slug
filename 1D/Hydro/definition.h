
! slope limiters
#define MINMOD  1
#define MC      2
#define VANLEER 3

! Riemann solvers
#define HLL     1
#define ROE     2

! primitive vars
#define DENS_VAR 1
#define VELX_VAR 2
#define PRES_VAR 3
#define EINT_VAR 4
#define GAMC_VAR 5
#define GAME_VAR 6
#define NUMB_VAR 6
#define NSYS_VAR 3 /* total number of equations of the conservative system */

! conservative vars
#define MOMX_VAR 2
#define ENER_VAR 3

! waves
#define SHOCKLEFT 1
#define CTENTROPY 2
#define SHOCKRGHT 3
#define NUMB_WAVE 3

! setup parameters
#define MAX_STRING_LENGTH 80

! BC
#define OUTFLOW  1
#define PERIODIC 2
#define REFLECT  3
#define USER     4

#define PI 4.*ATAN(1.)
