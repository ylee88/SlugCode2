#define PI 4.*ATAN(1.)

! MultiD vars
#define XDIM 1
#define YDIM 2
#define ZDIM 3
#define NDIM 3

! primitive vars
#define DENS_VAR 1
#define VELX_VAR 2
#define VELY_VAR 3
#define VELZ_VAR 4
#define PRES_VAR 5
#define BDRY_VAR 6 /* an extra boundary var to simulate an internal solid boundary structure */
#define EINT_VAR 7
#define GAMC_VAR 8
#define GAME_VAR 9
#define NUMB_VAR 9

! conservative vars
#define MOMX_VAR 2
#define MOMY_VAR 3
#define MOMZ_VAR 4
#define ENER_VAR 5
#define NSYS_VAR 5 /* total number of equations of the conservative system */

! waves
#define SHOCKLEFT 1
#define SLOWWLEFT 2
#define CTENTROPY 3
#define SHOCKRGHT 4
#define SLOWWRGHT 5
#define NUMB_WAVE 5

! setup parameters
#define MAX_STRING_LENGTH 100
