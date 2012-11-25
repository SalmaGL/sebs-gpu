#define pi 3.1415926536f
#define deg2rad 0.01745329252f

#define kel2deg -273.15f  /* Kelvin -> Celsius */

#define k 0.4f  /* von Karman constant */
#define Sigma_SB 5.678E-8f  /* Stefan-Boltzmann's constant (W/m2/K4) */
#define T0 273.15f  

#define Rso 1366.0f  /* Solar Constant(W/m2) */
#define g 9.81f  /* Gravity accelaration (kg s-2) */

#define Rmax 6378137.0f  /* the earth's equatorial radius (m) */
#define Rmin 6356752.0f  /* the earth's polar radius (m) */
 
#define Rd 287.04f  /* Gas Constant for Dry air, from table 2.1 P25 of Brutsaert 2005 (J kg-1 K-1) */
#define Rv 461.5f  /* Gas Constant for Water vapor, from table 2.1 P25 of Brutsaert 2005 (J kg-1 K-1) */

#define Cpw 1846.0f  /* specific heat coefficient for water vapor, J Kg-1 K-1 */
#define Cpd 1005.0f  /* specific heat coefficient for dry air, J Kg-1 K-1 */

#define Cd 0.2f  /* Foliage drag coefficient */
#define Ct 0.01f  /* Heat transfer coefficient */

#define gammaConst 67.0f  /* psychometric constant (Pa K-1) */


#define Pr 0.7f  /* Prandtl Prandtl number */
#define Pr_u 1.0f  /* Turbulent Prandtl number for unstable case */
#define Pr_s 0.95f  /* Turbulent Prandtl number for stable case */

#define ri_i 60.0f  /* surface resistance of standard crop, s m-1 */

/* The latent heat of vaporization at 30C from Brutsaert, 1982, p.41,tab. 3.4,
   more exact values can be obtained from eqn(3.22, 3.24a,b) */
#define L_e 2.430f   /* MJ Kg-1 */
#define rho_w 0.998f  /* density of water [kg/(m2 mm)] */

#define PSI0 1.3656120718024247f 

#define DT 86400.0f

#define MaxAllowedError 0.01f

/* MODIS Tile Sizes */
#define dimx 120
#define dimy 120

#define SelectedDevice 0

#define nThreadsPerBlock 64

