// Authors: Benjamin Fenker 2013
// Copyright Benjamin Fenker 2013

#ifndef INCLUDE_UNITS_H_
#define INCLUDE_UNITS_H_

/* Length [m] */
#define _m (1e0)                /* meters */
#define _cm (1e-2)              /* centimeters */
#define _nm (1e-9)              /* nanometers */
#define _fm (1e-15)             /* femtometers */

/* Mass [kg] */
#define _kg (1e0)               /* kilograms */

/* Time [s]*/
#define _s (1e0)                /* seconds */
#define _us (1e-6)              /* microseconds */
#define _ns (1e-9)              /* nanoseconds */

/* Current [A] */
#define _A (1e0)                /* amperes */

/* Temperature  [K] */
#define _K (1e0)                /* degrees kelvin */

/* Luminous Intensity [cd] */
#define _cd (1e0)               /* candelas */

/* Amount [mol] */
#define _mol (1e0)              /* moles */

/* Area [m^2] */
#define _m2 (1e0)               /* square-meters */
#define _cm2 (1e-4)             /* square-centimeters */

/* Frequency [s^-1] */
#define _Hz (1e0)               /* hertz */
#define _MHz (1e6)              /* megahertz */

/* Energy [J] */
#define _J (1e0)                /* joules */
#define _eV (1.602177e-19)      /* electronvolts */
#define _MeV (1.602177e-13)     /* megaelectronvolts */

/* Power [W = kg*m^2/s^3] */
#define _W (1e0)                /* watts */
#define _mW (1e-3)              /* milliwatts */
#define _uW (1e-6)              /* microwatts */

/* Magnetic Field [T = Tesla] */
#define _T (1e0)                /* tesla */
#define _G (1e-4)               /* gauss */

/* Electric charge [C = couloumb] */
#define _C (1e0)                /* couloumbs */

/* Electrical voltage [V = volt] */
#define _V (1e0)                         /* volts */

/* Electric dipole moment [C*m = couloumb-meter] */
#define _Cm (1e0)                        /* couloumb-meters */

/* Physical constants */
#define _speed_of_light (2.99792458e8) /* m/s */
#define _planck_h (6.62606896e-34)     /* m^3 / kg s^2 */
#define _planck_hbar (1.05457162825e-34) /* m^3 / kg s^2 */
#define _epsilon_0 (8.854e-12)           /* s^4 A^2 / kg m^3 */
#define _elementary_charge (1.60217646e-19) /* C */
#endif  // INCLUDE_UNITS_H_
