//
//  Defining input parameters for quadroptRoot.C

//  Defining macros

#define TRUE         1
#define FALSE        0

#define MAXORD       200    // Maximal possible value of polynomial order ( > 1)

// Precision of cl_F corresponding to MAXORD.
//
// A good guess is:                          DIGIT = 70+2.8*MAXORD
//
// but one has to check this by two runs with increasing precision.

// 700
#define DIGIT        700    // Precision of cl_F


//  Define constants to the desired precision

cl_F  ONE =    "1.0e+0_700";      // Precise 1
cl_F  TWO =    "2.0e+0_700";      // Precise 2
cl_F  ZERO =   "0.0e+0_700";      // Precise 0
cl_F  HALF =   "0.5e+0_700";      // Precise 0.5
cl_F  HUND =   "100.e+0_700";

// Define basic parameters to the desired precision

int   MAXPOW = 48;

/* cl_F  ALPHA =    "-0.500e+0_700", */
/*   EPSILON =  "0.1e+0_700", */
/*   LAMBDA =   "1.00e+0_700"; */

cl_F  ALPHA =    "0.500e+0_700",
  EPSILON =  "0.0043e+0_700",
  LAMBDA =   "1.e+0_700";

//  Define output format and files

char Format[] = "C";
char Filename[] = "recur_A25_8_002.cff";
char Filenamr[] = "roots_A25_8_002.cff";

/******************************************************************************/




