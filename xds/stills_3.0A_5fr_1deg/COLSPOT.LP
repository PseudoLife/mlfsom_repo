
 ***** COLSPOT ***** (VERSION Jan 31, 2020  BUILT=20200131)  13-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 MAXIMUM_NUMBER_OF_PROCESSORS=   4  (0: automatic choice of OpenMP threads)
 MAXIMUM_NUMBER_OF_JOBS=  16  (0:automatic choice of forked main tasks)
 number of forked COLSPOT tasks   1
 CLUSTER_NODES=
 SECONDS=     0
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_3.0A_5fr_1deg/stills_3.0A_5fr_1deg_???_001.img SMV
 DATA_RANGE=       5       5
 BACKGROUND_PIXEL=    6.00
 STRONG_PIXEL=    3.00
 MAXIMUM_NUMBER_OF_STRONG_PIXELS=   1500000
 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=    6 SPOT_MAXIMUM-CENTROID=    3.00
 SPOT_RANGE=       5       5
 ------------------------------------------------------------------------

 Task # 1 scans a subset of       1 images in the range      5 ...     5
 ------------------------------------------------------------------------

 RESULTS FOR FORKED COLSPOT TASK # 1
 USING    1 OUT OF    1 PROCESSORS
 PROCESS #  1 scans     1  images in the range       5 ...      5

 *** DEFINITION OF SYMBOLS ***
 NBKG      = NUMBER OF BACKGROUND PIXELS ON DATA IMAGE
 NSTRONG   = NUMBER OF STRONG PIXELS ON DATA IMAGE
 I/O-FLAG  = ERROR CODE AFTER ACCESSING DATA IMAGE
              0: NO ERROR
             -1: CANNOT OPEN OR READ IMAGE FILE
             -3: WRONG DATA FORMAT

 FRAME #      NBKG     NSTRONG   I/O-FLAG
       5   5224391        2727       0

 NUMBER OF STRONG PIXELS EXTRACTED FROM IMAGES        2727

 NUMBER OF DIFFRACTION SPOTS LOCATED                   203
 IGNORED BECAUSE OF SPOT CLOSE TO UNTRUSTED REGION      18
 WEAK SPOTS OMITTED                                      0
 NUMBER OF DIFFRACTION SPOTS ACCEPTED                  185

 total elapsed wall-clock time for COLSPOT       0.4 sec

 TASK     cpu time (sec)    elapsed wall-clock time (sec)
   1           0.4                     0.4
