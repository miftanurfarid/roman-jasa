
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* auditory filterbank constants */

#define MAX_CHANNEL         128                         /* maxmimum number of filters */
#define BW_CORRECTION       1.019                       /* ERB bandwidth correction 4th order */
#define SAMPLING_FREQUENCY  44100                       /* Hz */
#define MAX_DELAY           400                         /* corresponds to 80 Hz */
#define MAX_WINDOW          880                 /* use a window of 20 ms */
#define MAX_DISP            44
#define MAX_BUFFER_SIZE     2*MAX_DISP+MAX_WINDOW       /* buffer size */
#define OFFSET              440                         /* compute acg every 10 ms */
#define HWIN                MAX_WINDOW/2.0

/* hair cell constants from Meddis 1988 paper */


#define MED_Y 5.05
#define MED_G 2000.0
#define MED_L 2500.0
#define MED_R 6580.0
#define MED_X 66.31
#define MED_A 3.0
#define MED_B 300.0
#define MED_H 48000.0
#define MED_M 1.0

/* outer/middle ear */

#define MIDDLE_EAR_SIZE 29
#define DB 60.0

/* frequency scale definitions from Moore and Glasberg 1990 */

#define erb(f) (24.7*(4.37e-3*(f)+1.0))
#define hzToERBrate(f) (21.4*log10(4.37e-3*(f)+1.0))
#define ERBrateToHz(f) ((pow(10.0,((f)/21.4))-1.0)/4.37e-3)
#define sqr(x) ((x)*(x))

#define LEFT 0
#define RIGHT 1


#define MAX_FRAME 700
#define MAX_SIGNAL 100000

/* constants of name files */
#define TARGET_DIR "target"

#define INPUTFILE_LEFT "Left"
#define INPUTFILE_RIGHT "Right"

#define INPUTFILE_TARGET "Left1" // target signal at better ear  
#define INPUTFILE_NOISE "Left2"  // interference signal at better ear
#define INPUTFILE_MIXT INPUTFILE_LEFT  // mixture signal at better ear

#define INPUTFILE_MASK "mask"
#define INPUTFILE_ITD "itd"
#define INPUTFILE_IID "iid"
#define INPUTFILE_R "ratio"
#define INPUTFILE_RESYNTH "resynth"

#define TRAINFILE_ITD "ITD"
#define TRAINFILE_IID "IID"
#define TRAINFILE_R "RATIO"

#define TEST_DIR "test"

/* constants for training */

#define TRAIN 1
#define COMPUTE_R 1


