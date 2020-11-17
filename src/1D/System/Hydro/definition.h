#define PI 4.*ATAN(1.)

! MultiD vars
#define XDIM 1
#define NDIM 1

! primitive vars
#define DENS_VAR 1
#define VELX_VAR 2
#define PRES_VAR 3
#define BDRY_VAR 4 /* an extra boundary var to simulate an internal solid boundary structure */
#define EINT_VAR 5
#define GAMC_VAR 6
#define GAME_VAR 7
#define NUMB_VAR 7

! conservative vars
#define MOMX_VAR 2
#define ENER_VAR 3
#define NSYS_VAR 3 /* total number of equations of the conservative system */

! waves
#define SHOCKLEFT 1
#define CTENTROPY 2
#define SHOCKRGHT 3
#define NUMB_WAVE 3

! setup parameters
#define MAX_STRING_LENGTH 100
