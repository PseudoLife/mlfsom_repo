
 ***** COLSPOT ***** (VERSION Jan 31, 2020  BUILT=20200131)  29-Feb-2020

 INPUT PARAMETER VALUES
 ----------------------
 MAXIMUM_NUMBER_OF_PROCESSORS=   0  (0: automatic choice of OpenMP threads)
 MAXIMUM_NUMBER_OF_JOBS=   0  (0:automatic choice of forked main tasks)
 number of forked COLSPOT tasks   1
 CLUSTER_NODES=
 SECONDS=     0
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/single/lm_h97_c01_rt_10_???.cbf CBF
 DATA_RANGE=       1       1
 BACKGROUND_PIXEL=    6.00
 STRONG_PIXEL=    3.00
 MAXIMUM_NUMBER_OF_STRONG_PIXELS=   1500000
 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=    6 SPOT_MAXIMUM-CENTROID=    3.00
 SPOT_RANGE=       1       1
 ------------------------------------------------------------------------

 Task # 1 scans a subset of       1 images in the range      1 ...     1
 ------------------------------------------------------------------------

 RESULTS FOR FORKED COLSPOT TASK # 1
 USING    1 OUT OF    1 PROCESSORS
 PROCESS #  1 scans     1  images in the range       1 ...      1

 *** DEFINITION OF SYMBOLS ***
 NBKG      = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NSTRONG   = NUMBER OF STRONG PIXELS ON DATA IMAGE
 I/O-FLAG  = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT

 FRAME #      NBKG     NSTRONG   I/O-FLAG
       1   5222611        3858       0

 NUMBER OF STRONG PIXELS EXTRACTED FROM IMAGES        3858

 NUMBER OF DIFFRACTION SPOTS LOCATED                   261
 IGNORED BECAUSE OF SPOT CLOSE TO UNTRUSTED REGION      25
 WEAK SPOTS OMITTED                                      0
 NUMBER OF DIFFRACTION SPOTS ACCEPTED                  236

 total elapsed wall-clock time for COLSPOT       0.4 sec

 TASK     cpu time (sec)    elapsed wall-clock time (sec)
   1           0.3                     0.3
