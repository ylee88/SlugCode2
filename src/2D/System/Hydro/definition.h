#define PI 4.*ATAN(1.)

! MultiD vars
#define XDIM 1
#define YDIM 2
#define NDIM 2

! primitive vars
#define DENS_VAR 1
#define VELX_VAR 2
#define VELY_VAR 3
#define PRES_VAR 4
#define BDRY_VAR 5 /* an extra boundary var to simulate an internal solid boundary structure */
#define EINT_VAR 6
#define GAMC_VAR 7
#define GAME_VAR 8
#define NUMB_VAR 8

! conservative vars
#define MOMX_VAR 2
#define MOMY_VAR 3
#define ENER_VAR 4
#define NSYS_VAR 4 /* total number of equations of the conservative system */

! waves
#define SHOCKLEFT 1
#define SLOWWLEFT 2
#define CTENTROPY 3
#define SHOCKRGHT 4
#define NUMB_WAVE 4

! setup parameters
#define MAX_STRING_LENGTH 100
