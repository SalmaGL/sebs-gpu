#define pi 3.1415926536
#define deg2rad 0.01745329252

#define kel2deg -273.15 /* Kelvin -> Celsius */

#define k 0.4 /* von Karman constant */
#define Sigma_SB 5.678E-8 /* Stefan-Boltzmann's constant (W/m2/K4) */
#define T0 273.15 

#define Rso 1366.0 /* Solar Constant(W/m2) */
#define g 9.81 /* Gravity accelaration (kg s-2) */

#define Rmax 6378137.0 /* the earth's equatorial radius (m) */
#define Rmin 6356752.0 /* the earth's polar radius (m) */
 
#define Rd 287.04 /* Gas Constant for Dry air, from table 2.1 P25 of Brutsaert 2005 (J kg-1 K-1) */
#define Rv 461.5 /* Gas Constant for Water vapor, from table 2.1 P25 of Brutsaert 2005 (J kg-1 K-1) */

#define Cpw 1846.0 /* specific heat coefficient for water vapor, J Kg-1 K-1 */
#define Cpd 1005.0 /* specific heat coefficient for dry air, J Kg-1 K-1 */

#define Cd 0.2 /* Foliage drag coefficient */
#define Ct 0.01 /* Heat transfer coefficient */

#define gammaConst 67.0 /* psychometric constant (Pa K-1) */


#define Pr 0.7 /* Prandtl Prandtl number */
#define Pr_u 1.0 /* Turbulent Prandtl number for unstable case */
#define Pr_s 0.95 /* Turbulent Prandtl number for stable case */

#define ri_i 60.0 /* surface resistance of standard crop, s m-1 */

/* The latent heat of vaporization at 30C from Brutsaert, 1982, p.41,tab. 3.4,
   more exact values can be obtained from eqn(3.22, 3.24a,b) */
#define L_e 2.430  /* MJ Kg-1 */
#define rho_w 0.998 /* density of water [kg/(m2 mm)] */

#define PSI0 1.3656120718024247

#define dimx 120
#define dimy 120
