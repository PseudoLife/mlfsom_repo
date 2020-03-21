
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.285    77.285    38.675  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=       2       2
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       2
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  1.000000 -0.000556  0.000711
 FRACTION_OF_POLARIZATION=   0.900
 POLARIZATION_PLANE_NORMAL=  0.000000  1.000000  0.000000
 AIR=   0.00034
 DETECTOR=PILATUS
 SILICON=  3.942868 SENSOR_THICKNESS=  0.320000
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463  NY=  2527    QX=  0.172000  QY=  0.172000
 ORGX=   1263.00  ORGY=   1231.00
 DETECTOR_DISTANCE=   250.400
 DIRECTION_OF_DETECTOR_X-AXIS=   0.00000   0.00000   1.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000  -1.00000   0.00000
 BEAM_DIVERGENCE_E.S.D.=   0.039
 REFLECTING_RANGE_E.S.D.=   0.069
 MINIMUM_ZETA=     0.050
 MAXIMUM_ERROR_OF_SPOT_POSITION=   3.0
 MAXIMUM_ERROR_OF_SPINDLE_POSITION=   2.0
 MINIMUM_I/SIGMA=   3.0
 REFLECTIONS/CORRECTION_FACTOR=      50
 SNRC=  0.10 I/sigma for use in scaling factor determination
 STRICT_ABSORPTION_CORRECTION=FALSE 
 CORRECTIONS= DECAY MODULATION ABSORPTION                      
 REFERENCE_DATA_SET=                                                  
 FIT_B-FACTOR_TO_REFERENCE_DATA_SET=FALSE 

    1558 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       2 ...       2
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1558 REFLECTIONS ACCEPTED



 ******************************************************************************
                      AUTOMATIC SPACE GROUP ASSIGNMENT
 ******************************************************************************

 XDS adopts the following approach.
 (1) it looks for possible symmetries of the crystal lattice
 (2) it computes a redundancy independent R-factor for all enantiomorphous
     point groups compatible with the observed lattice symmetry.
 (3) it selects the group which explains the intensity data at an acceptable,
     redundancy-independent R-factor (Rmeas, Rrim) using a minimum number of
     unique reflections.

 This approach does not test for the presence of screw axes. Consequently,
 orthorhombic cell axes will be specified in increasing length (following
 conventions), despite the possibility that different assignments for the
 cell axes could become necessary for space groups P222(1) and P2(1)2(1)2
 containing one or two screw axes, respectively.

 The user can always override the automatic decisions by specifying the
 correct space group number and unit cell constants in XDS.INP and repeating
 the CORRECT step of XDS. This provides a simple way to rename orthorhombic
 cell constants if screw axes are present.

 In addition, the user has the option to specify in XDS.INP
 (a) a reference data set or
 (b) a reindexing transformation or
 (c) the three basis vectors (if known from processing a previous data set
 taken at the same crystal orientation in a multi-wavelength experiment).
 These features of XDS are useful for resolving the issue of alternative
 settings of polar or rhombohedral cells (like P4, P6, R3).

 =====> Specifications provided by the user in XDS.INP:
 SPACE_GROUP_NUMBER=   96 (0: unknown symmetry and unit cell)
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_B-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_C-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 REIDX=  0  0  0  0  0  0  0  0  0  0  0  0 (all 0 : not specified)
 TEST_RESOLUTION_RANGE=    8.00    4.50 (Angstrom) for space group determination
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               48
 NUMBER OF ACCEPTED UNIQUE REFLECTIONS FROM REFERENCE DATA SET       0


 *********** DETERMINATION OF LATTICE CHARACTER AND BRAVAIS LATTICE ***********

 The CHARACTER OF A LATTICE is defined by the metrical parameters of its
 reduced cell as described in the INTERNATIONAL TABLES FOR CRYSTALLOGRAPHY
 Volume A, p. 746 (KLUWER ACADEMIC PUBLISHERS, DORDRECHT/BOSTON/LONDON, 1989).
 Note that more than one lattice character may have the same BRAVAIS LATTICE.

 A lattice character is marked "*" to indicate a lattice consistent with the
 observed locations of the diffraction spots. These marked lattices must have
 low values for the QUALITY OF FIT and their implicated UNIT CELL CONSTANTS
 should not violate the ideal values by more than
 MAX_CELL_AXIS_ERROR=  0.03
 MAX_CELL_ANGLE_ERROR=           2.0 (Degrees)

 The REINDEXING TRANSFORMATION REIDX() consists of 12 integers that relate
 the original indices H,K,L from file INTEGRATE.HKL to the indices H',K',L'
 with respect to the new cell.
    H'=(REIDX(1)*H+REIDX( 2)*K+REIDX( 3)*L)/IDXV+REIDX( 4)
    K'=(REIDX(5)*H+REIDX( 6)*K+REIDX( 7)*L)/IDXV+REIDX( 8)
    L'=(REIDX(9)*H+REIDX(10)*K+REIDX(11)*L)/IDXV+REIDX(12)
 The value of the integer IDXV depends on the lattice type used for
 specifying reflections on file INTEGRATE.HKL; IDXV is
 1 for a primitive,
 2 for a face or body centred,
 3 for a rhombohedral,
 4 for a lattice centred on all faces.
 IDXV is set by XDS and cannot be input by the user.

  LATTICE-  BRAVAIS-   QUALITY  UNIT CELL CONSTANTS (ANGSTROEM & DEGREES)    REINDEXING TRANSFORMATION
 CHARACTER  LATTICE     OF FIT      a      b      c   alpha  beta gamma

 *  44        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0 -1  0  0  0  0 -1  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  25        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  20        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.3   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.3   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    28        mC        250.0      38.7  159.3   77.3  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    29        mC        250.0      38.7  159.3   77.3  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    36        oC        250.0      38.7  159.3   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.7  159.3   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    27        mC        500.0     159.3   38.7  109.3  90.0 133.3  76.0   -2  0 -1  0  0  0 -1  0  1 -1  0  0
    19        oI        500.0      38.7  109.3  115.9  90.0  70.5  90.0    0  0 -1  0  1 -1  0  0 -1 -1 -1  0
    26        oF        625.0      38.7  159.3  159.3  86.6 104.0 104.0    0  0  1  0 -2  0 -1  0  0 -2 -1  0
    18        tI        625.0     109.3  115.9   38.7  70.5  90.0  90.0    1 -1  0  0  1  1  1  0  0  0  1  0
    30        mC        998.2      77.3  172.8   38.7  90.0  90.0  63.4   -1  0  0  0 -1  2  0  0  0  0 -1  0
    41        mC        998.2     172.8   77.3   38.7  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    40        oC        998.2      77.3  172.8   38.7  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.2      77.3   77.3   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     115.9  115.9  115.9  96.4  96.4 141.0    1 -1  1  0 -1  1  1  0 -1 -1 -1  0
     2        hR        999.0      86.4   86.4  115.9 116.5  63.5 101.6    1  0  1  0  0 -1 -1  0  1  1 -1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.4   86.4  115.9 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.3   86.4   86.4  78.4  50.8  50.8   -1 -1  0  0  0 -1  1  0 -1  0  1  0
     7        tI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    0 -1  1  0 -1  0  1  0 -1 -1  0  0
     8        oI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    1  0 -1  0  0  1 -1  0  1  1  0  0
     9        hR        999.0      38.7   86.4  247.4 102.1  99.0 116.6    0  0  1  0 -1  0 -1  0  1 -3 -1  0
    10        mC        999.0      86.4   86.4   77.3  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.7   77.3  177.1  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.4   86.4  177.1 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.4   86.4   86.4  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.1  109.3   38.7  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.7   77.3  177.1 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.4  177.1   77.3 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

 ****** LATTICE SYMMETRY IMPLICATED BY SPACE GROUP SYMMETRY ******

BRAVAIS-            POSSIBLE SPACE-GROUPS FOR PROTEIN CRYSTALS
  TYPE                     [SPACE GROUP NUMBER,SYMBOL]
  aP      [1,P1]
  mP      [3,P2] [4,P2(1)]
 mC,mI    [5,C2]
  oP      [16,P222] [17,P222(1)] [18,P2(1)2(1)2] [19,P2(1)2(1)2(1)]
  oC      [21,C222] [20,C222(1)]
  oF      [22,F222]
  oI      [23,I222] [24,I2(1)2(1)2(1)]
  tP      [75,P4] [76,P4(1)] [77,P4(2)] [78,P4(3)] [89,P422] [90,P42(1)2]
          [91,P4(1)22] [92,P4(1)2(1)2] [93,P4(2)22] [94,P4(2)2(1)2]
          [95,P4(3)22] [96,P4(3)2(1)2]
  tI      [79,I4] [80,I4(1)] [97,I422] [98,I4(1)22]
  hP      [143,P3] [144,P3(1)] [145,P3(2)] [149,P312] [150,P321] [151,P3(1)12]
          [152,P3(1)21] [153,P3(2)12] [154,P3(2)21] [168,P6] [169,P6(1)]
          [170,P6(5)] [171,P6(2)] [172,P6(4)] [173,P6(3)] [177,P622]
          [178,P6(1)22] [179,P6(5)22] [180,P6(2)22] [181,P6(4)22] [182,P6(3)22]
  hR      [146,R3] [155,R32]
  cP      [195,P23] [198,P2(1)3] [207,P432] [208,P4(2)32] [212,P4(3)32]
          [213,P4(1)32]
  cF      [196,F23] [209,F432] [210,F4(1)32]
  cI      [197,I23] [199,I2(1)3] [211,I432] [214,I4(1)32]


 Several solutions exist for the given space group that explain the data. The
 selected solution (marked "*") correlates best with the given reference data set.
 The independent settings of the cell can be obtained by application of the
 REINDEXING TRANSFORMATION to the original indices H,K,L from file INTEGRATE.HKL

 CORRELATION = correlation factor with the reference data set
 NPAIR       = number of unique reflection pairs correlated 
 Rmeas=Rrim  = redundancy independent R-factor (intensities)
               For definition see:
               Rmeas : Diederichs & Karplus,
                       Nature Struct. Biol. 4, 269-275 (1997);
               Rrim  : Weiss & Hilgenfeld,
                       J.Appl.Cryst. 30,203-205 (1997).
 COMPARED    = number of reflections used for calculating Rmeas
 ESD         = Agreement with given unit cell geometry

 CORRELATION  NPAIR  Rmeas  COMPARED  ESD       REINDEX TRANSFORMATION

  *   0.00        0    5.4      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.657   -23.401     0.042
 UNIT_CELL_B-AXIS=    23.401    73.657     0.012
 UNIT_CELL_C-AXIS=    -0.022     0.001    38.676
 REFLECTIONS H,K,L FROM "INTEGRATE.HKL" WILL BE REINDEXED BY
 H' = (  1*H +  0*K +  0*L)/ 1 +  0
 K' = (  0*H +  1*K +  0*L)/ 1 +  0
 L' = (  0*H +  0*K +  1*L)/ 1 +  0

 *****  8 EQUIVALENT POSITIONS IN SPACE GROUP # 96 *****

    If x',y',z' is an equivalent position to x,y,z, then
        x'=x*ML(1)+y*ML( 2)+z*ML( 3)+ML( 4)/12.0
        y'=x*ML(5)+y*ML( 6)+z*ML( 7)+ML( 8)/12.0
        z'=x*ML(9)+y*ML(10)+z*ML(11)+ML(12)/12.0

    #    1  2  3  4    5  6  7  8    9 10 11 12 
    1    1  0  0  0    0  1  0  0    0  0  1  0
    2   -1  0  0  0    0 -1  0  0    0  0  1  6
    3    0 -1  0  6    1  0  0  6    0  0  1  9
    4    0  1  0  6   -1  0  0  6    0  0  1  3
    5    0  1  0  0    1  0  0  0    0  0 -1  0
    6    0 -1  0  0   -1  0  0  0    0  0 -1  6
    7   -1  0  0  6    0  1  0  6    0  0 -1  9
    8    1  0  0  6    0 -1  0  6    0  0 -1  3

 ******************************************************************************
      MEAN DISCREPANCIES BETWEEN OBSERVED AND CALCULATED SPOT LOCATIONS
 ******************************************************************************

 The discrepancies in X- and Y-coordinates of the spots are depicted in the
 two images DX-CORRECTIONS.cbf and DY-CORRECTIONS.cbf for inspection with
 the XDS-Viewer.



 ******************************************************************************
  REFINEMENT OF DIFFRACTION PARAMETERS USING ALL IMAGES
 ******************************************************************************


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       352 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.16
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.06
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.267    77.267    38.669  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.3E-01 1.3E-01 6.3E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012942  0.012942  0.025861  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.644   -23.386     0.030
 COORDINATES OF UNIT CELL B-AXIS    23.386    73.644     0.008
 COORDINATES OF UNIT CELL C-AXIS    -0.016     0.001    38.669
 CRYSTAL MOSAICITY (DEGREES)     0.069
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999976  0.006496  0.002198
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.85   1232.90
 DETECTOR ORIGIN (PIXELS) AT                     1260.65   1242.35
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.23
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.27    77.27    38.67  90.000  90.000  90.000


 ******************************************************************************
        MEAN INTENSITY AS FUNCTION OF SPINDLE POSITION WITHIN DATA IMAGE
 ******************************************************************************

 This statistics could serve as a diagnostic tool for detecting
 shutter problems (suggested by Kay Diederichs).
 Data are corrected for this effect if PATCH_SHUTTER_PROBLEM=TRUE in XDS.INP.
 ===> Selected : PATCH_SHUTTER_PROBLEM=FALSE

 INTERVAL = Angular interval in units of the oscillation range
            covered by a data image.
 NUMBER   = Number of reflections in each interval
 INTENSITY= Mean intensity of the reflections in each interval
 FACTOR   = Correction factor applied to intensities;
            INTENSITY(overall)/INTENSITY(interval)

   INTERVAL     NUMBER      INTENSITY      FACTOR
  0.05  0.15       143        337.195       1.019
  0.15  0.25       186        338.481       1.015
  0.25  0.35       183        295.823       1.161
  0.35  0.45       163        220.404       1.559
  0.45  0.55       181        422.410       0.813
  0.55  0.65       192        303.193       1.133
  0.65  0.75       188        248.757       1.381
  0.75  0.85       165        334.389       1.027
  0.85  0.95       153        636.649       0.540
 
 -0.05  0.95      1558        343.541       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               243.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.985
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     1.0 XMAX=     2.0 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59736 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        528


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               243.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.984
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    11.5 XMAX=  2454.6 NXBIN=    3
 YMIN=     6.0 YMAX=  2520.3 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        528


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               243.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.979
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     1.0 XMAX=     2.0 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1263
 DETECTOR_SURFACE_POSITION=    1648    1690
 DETECTOR_SURFACE_POSITION=     818    1690
 DETECTOR_SURFACE_POSITION=     818     836
 DETECTOR_SURFACE_POSITION=    1648     836
 DETECTOR_SURFACE_POSITION=    2172    1663
 DETECTOR_SURFACE_POSITION=    1622    2230
 DETECTOR_SURFACE_POSITION=     844    2230
 DETECTOR_SURFACE_POSITION=     294    1663
 DETECTOR_SURFACE_POSITION=     294     863
 DETECTOR_SURFACE_POSITION=     844     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2172     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        528


 ******************************************************************************
    CORRECTION PARAMETERS FOR THE STANDARD ERROR OF REFLECTION INTENSITIES
 ******************************************************************************

 The variance v0(I) of the intensity I obtained from counting statistics is
 replaced by v(I)=a*(v0(I)+b*I^2). The model parameters a, b are chosen to
 minimize the discrepancies between v(I) and the variance estimated from
 sample statistics of symmetry related reflections. This model implicates
 an asymptotic limit ISa=1/SQRT(a*b) for the highest I/Sigma(I) that the
 experimental setup can produce (Diederichs (2010) Acta Cryst D66, 733-740).

     a        b          ISa
 1.243E+00  3.021E-03   16.32


 ******************************************************************************
      STANDARD ERROR OF REFLECTION INTENSITIES AS FUNCTION OF RESOLUTION
      FOR DATA SET  XDS_ASCII.HKL                                     
 ******************************************************************************

  I/Sigma  = mean intensity/Sigma of a reflection in shell
  Chi^2    = goodness of fit between sample variances of
             symmetry-related intensities and their errors
             (Chi^2 = 1 for perfect agreement)
  R-FACTOR
  observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
  expected = expected R-FACTOR derived from Sigma(I)

   NUMBER  = number of reflections in resolution shell
             used for calculation of R-FACTOR
  ACCEPTED = number of accepted reflections
  REJECTED = number of rejected reflections (MISFITS),
             recognized by comparison with symmetry-related
             reflections.

 RESOLUTION RANGE  I/Sigma  Chi^2  R-FACTOR  R-FACTOR  NUMBER ACCEPTED REJECTED
                                   observed  expected

   24.434   3.395    14.81   1.13      2.40      3.52      66      95      20
    3.395   2.412    12.47   1.60      5.14      3.84     132     175      12
    2.412   1.973     9.51   1.04      5.02      4.53     180     223      10
    1.973   1.710     0.05   0.67    591.12    722.39     202     255       0
    1.710   1.530     0.03   0.71   1088.53   1274.04     218     296       0
    1.530   1.397    -0.00   0.63   4222.71   5325.19     246     319       2
    1.397   1.294     0.07   0.59    324.29    412.91     100     148       2
   --------------------------------------------------------------------------
   24.434   1.294     3.80   0.86      5.12      5.39    1144    1511      46
 


 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR VARIOUS SUBSETS OF INCLUDED DATA IMAGES
 ******************************************************************************

 Data set statistics is reported below several times, each with a different
 upper limit on the number of images included. This provides the user with
 the information for deciding which data images should be excluded from the
 final data set because of radiation damage or other defects. If the user
 decides for a subset of "good" images that differs from the specification
 as given by the input parameter (see DATA_RANGE= in XDS.INP) the CORRECT-
 step of XDS must be repeated with the new parameter values for DATA_RANGE=.


 ******************************************************************************
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       2       2)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.680   0.2868E+01  0.1204E+02         0.24           1 
   47    0    0     1.644  -0.2383E+02  0.1212E+02        -1.97           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -6.6% OF MEAN INTENSITY


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.84          72      47      1285        3.7%       2.2%      3.5%       50   18.14      3.2%    99.9*     0    0.000       0
     2.73         112      71      2098        3.4%       4.7%      3.7%       82   16.48      6.6%    99.2*     0    0.000       0
     2.23         166     102      2665        3.8%       4.8%      4.0%      128   14.99      6.8%    98.9*     0    0.000       0
     1.94         170     101      3101        3.3%       5.4%      5.3%      138    8.69      7.6%    99.1*     0    0.000       0
     1.73         196     118      3479        3.4%     583.2%    698.2%      156    0.02    824.7%    -8.7      0    0.000       0
     1.58         227     140      3850        3.6%     620.4%    699.2%      174    0.10    877.3%   -18.5      0    0.000       0
     1.47         255     160      4143        3.9%    3858.8%   5211.1%      190   -0.03   5457.2%    16.1      0    0.000       0
     1.37         224     137      4428        3.1%    4289.1%   5210.7%      174   -0.01   6065.7%   -17.9      0    0.000       0
     1.29          89      63      4744        1.3%     186.4%    227.7%       52    0.15    263.6%    34.9      0    0.000       0
    total        1511     939     29793        3.2%       5.1%      5.4%     1144    4.74      7.2%    99.7*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1558
 NUMBER OF REJECTED MISFITS                                46
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1511
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    939



 ******************************************************************************
    WILSON STATISTICS OF DATA SET  "XDS_ASCII.HKL"
 ******************************************************************************

 Data is divided into resolution shells and a straight line 
 A - 2*B*SS is fitted to log<I>, where
   RES    = mean resolution (Angstrom) in shell
   SS     = mean of (sin(THETA)/LAMBDA)**2 in shell
   <I>    = mean reflection intensity in shell
   BO     = (A - log<I>)/(2*SS)
    #     = number of reflections in resolution shell

   WILSON LINE (using all data) : A=   8.687 B=  34.786 CORRELATION=  0.90
      #      RES      SS        <I>       log(<I>)       BO
      17     7.549   0.004  3.2279E+03       8.080      69.2
      25     4.655   0.012  3.0407E+03       8.020      28.9
      32     3.508   0.020  2.2179E+03       7.704      24.2
      34     2.993   0.028  7.6009E+02       6.633      36.8
      41     2.644   0.036  9.5636E+02       6.863      25.5
      55     2.394   0.044  7.0802E+02       6.562      24.3
      44     2.205   0.051  4.8982E+02       6.194      24.2
      48     2.054   0.059  2.6249E+02       5.570      26.3
      53     1.930   0.067 -9.9019E-01    -999.999   -1000.0
      64     1.826   0.075 -5.2618E-02    -999.999   -1000.0
      56     1.736   0.083  1.7798E+00       0.576      48.9
      61     1.663   0.090  1.1763E+00       0.162      47.1
      66     1.595   0.098 -3.4732E-01    -999.999   -1000.0
      84     1.534   0.106  2.0598E-01      -1.580      48.3
      81     1.478   0.114 -7.8680E-01    -999.999   -1000.0
      62     1.430   0.122 -7.0741E-01    -999.999   -1000.0
      55     1.390   0.129 -6.9672E-02    -999.999   -1000.0
      41     1.350   0.137  1.5399E+00       0.432      30.1
      20     1.313   0.145  2.5546E+00       0.938      26.7


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.549        0.916        0.423        0.155
       3     4.655        1.240        0.739        0.329
       4     3.508        4.596        6.515        6.767
       4     2.993        1.355        0.669        0.250
       3     2.644        0.245        0.059        0.011
       5     2.394        0.667        0.281        0.089
       3     2.205        1.112        0.416        0.112
       8     2.054        0.213        0.066        0.016
       7     1.930        9.036       10.268       16.595
      11     1.736        5.229        1.368        4.189
      11     1.663       17.179       10.633       78.737
       4     1.595       31.417       80.466      153.002
       3     1.534      170.370      191.822     2626.037
       2     1.350       92.111      432.336     1451.048
      94   overall       10.458       17.112      124.715


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.549        1.241        1.781        2.104
      22     4.655        1.110        1.291        1.310
      28     3.508        0.574        0.366        0.206
      30     2.993        0.762        0.539        0.315
      38     2.644        1.165        1.258        1.210
      50     2.394        0.954        0.785        0.541
      41     2.205        0.785        0.564        0.348
      40     2.054        1.292        1.405        1.393
      53     1.826       22.107        0.012      226.082
      45     1.736        4.800        3.229       12.587
      50     1.663       23.264       33.105      238.937
      81     1.534      733.500     2279.645   291754.984
      54     1.390       27.266       -1.323      457.881
      39     1.350       25.269      101.674      864.494
      20     1.313       19.688       85.598      553.919
     845   overall       73.877      227.659    28035.116

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.549  0.00  0.00  0.00  0.00  0.00  0.36  0.33  0.64  0.61  0.59
       3     4.655  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     3.508  1.01  0.72  0.60  0.53  0.48  0.45  0.42  0.40  0.38  0.73
       4     2.993  0.00  0.00  0.00  0.00  0.00  0.00  0.42  0.40  0.38  0.37
       3     2.644  1.34  0.96  0.80  0.70  0.64  0.59  1.12  1.06  1.01  0.98
       5     2.394  0.00  0.00  0.00  0.42  1.15  1.07  1.00  0.95  0.91  0.88
       3     2.205  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.054  2.01  1.45  1.50  1.32  1.20  1.11  1.26  1.19  1.14  1.10
       7     1.930  2.30  1.65  1.37  1.51  1.37  1.27  1.20  1.14  1.09  1.05
      11     1.736  1.46  1.05  0.87  0.77  0.70  0.65  0.61  0.58  0.55  0.53
      11     1.663  1.10  0.79  0.65  0.58  0.52  0.49  0.46  0.43  0.42  0.53
       4     1.595  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       3     1.534  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     1.350  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.76  0.73
      94   overall  0.99  0.71  0.61  0.58  0.57  0.55  0.57  0.56  0.55  0.56


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.549  1.75  2.30  1.93  2.27  1.91  1.66  1.49  1.36  1.26  1.19
      22     4.655  1.91  1.76  1.23  1.10  1.16  1.21  1.17  1.07  1.07  1.08
      28     3.508  1.13  0.79  0.83  0.87  1.09  1.11  1.06  1.23  1.14  1.13
      30     2.993  1.75  1.65  1.29  1.31  1.10  1.03  1.06  1.09  1.12  1.05
      38     2.644  1.11  1.16  1.22  1.20  1.07  0.93  1.05  1.10  1.06  1.04
      50     2.394  1.47  1.32  1.08  1.15  1.07  1.11  1.11  1.09  1.04  0.98
      41     2.205  1.54  1.21  0.94  0.81  1.05  0.92  0.92  0.97  1.03  1.04
      40     2.054  2.10  1.24  1.16  1.06  0.95  0.89  0.84  0.95  0.97  0.91
      53     1.826  5.16  2.71  1.89  1.49  1.25  1.09  0.97  0.89  0.83  0.78
      45     1.736  4.20  2.21  1.54  1.35  1.19  1.03  0.93  0.85  0.79  0.77
      50     1.663  4.62  2.43  1.70  1.40  1.22  1.06  0.95  0.87  0.81  0.76
      81     1.534  5.45  2.86  2.00  1.57  1.32  1.15  1.03  0.94  0.89  0.84
      54     1.390  5.06  2.66  1.86  1.46  1.22  1.07  0.96  0.87  0.81  0.76
      39     1.350  5.39  2.97  2.08  1.71  1.43  1.25  1.12  1.02  0.95  0.89
      20     1.313  4.73  2.48  1.74  1.36  1.27  1.11  0.99  0.91  0.93  0.95
     845   overall  2.51  1.48  1.10  0.95  0.85  0.77  0.72  0.70  0.67  0.65


 List of    157 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   31   22   10    1.80  329.95 -0.1736E+02  0.1287E+02 "alien"
   33   32   16    1.38  314.03 -0.2188E+02  0.1413E+02 "alien"
   34   22    4    1.87  259.63 -0.1366E+02  0.9456E+01 "alien"
   33   23    7    1.81  206.18 -0.1085E+02  0.1326E+02 "alien"
   43    1    2    1.79  179.18 -0.9428E+01  0.8334E+01 "alien"
   28   21   22    1.38  176.55 -0.1230E+02  0.7754E+01 "alien"
   31   30   18    1.38  171.96 -0.1198E+02  0.9370E+01 "alien"
   36   32   14    1.39  161.49 -0.1125E+02  0.1578E+02 "alien"
   50    5   11    1.41  133.61 -0.9309E+01  0.1365E+02 "alien"
   42    1    3    1.82  132.11 -0.6951E+01  0.8678E+01 "alien"
   41    1    3    1.86  128.65 -0.6770E+01  0.8960E+01 "alien"
   39   33   10    1.41  123.30 -0.8590E+01  0.1030E+02 "alien"
   33   22    6    1.86  121.17 -0.6376E+01  0.8686E+01 "alien"
   15    7   20    1.79  111.59 -0.5872E+01  0.4536E+01 "alien"
   34    1   12    1.86  111.11 -0.5847E+01  0.7492E+01 "alien"
   30   21   10    1.85  110.29 -0.5803E+01  0.1193E+02 "alien"
   29   21   11    1.84  110.15 -0.5796E+01  0.8025E+01 "alien"
   54    4    6    1.39  109.76 -0.7647E+01  0.9595E+01 "alien"
   28   12   23    1.40  106.75 -0.7438E+01  0.9466E+01 "alien"
   51    4   10    1.41  104.31 -0.7268E+01  0.1374E+02 "alien"
   27   20   22    1.40  100.25 -0.6985E+01  0.7682E+01 "alien"
   42    1    4    1.81   96.66 -0.5086E+01  0.8615E+01 "alien"
   46    7   16    1.37   95.23 -0.6635E+01  0.8947E+01 "alien"
   51    2    1    1.51   94.64  0.1949E+02  0.1124E+02 "alien"
   47    2    9    1.53   88.83  0.1830E+02  0.1607E+02 "alien"
   17   15   18    1.82   87.97 -0.4629E+01  0.5150E+01 "alien"
   43    1    1    1.79   83.83 -0.4411E+01  0.8422E+01 "alien"
   37   29    9    1.54   82.16  0.1692E+02  0.9680E+01 "alien"
   41    5   15    1.51   81.32  0.1675E+02  0.1232E+02 "alien"
   49    2    7    1.52   80.98  0.1668E+02  0.9673E+01 "alien"
   30    3   15    1.82   78.96 -0.4155E+01  0.9159E+01 "alien"
   34   31   15    1.41   75.07 -0.5231E+01  0.9255E+01 "alien"
   54    4    7    1.38   72.92 -0.5080E+01  0.9689E+01 "alien"
   40   33    9    1.41   70.69 -0.4925E+01  0.9965E+01 "alien"
   46    1    6    1.63   69.12 -0.2401E+02  0.1251E+02 "alien"
   38    9   20    1.38   67.02 -0.4669E+01  0.1164E+02 "alien"
   14    5   20    1.81   65.87 -0.3466E+01  0.3942E+01 "alien"
   27   20   12    1.87   62.29 -0.3278E+01  0.7695E+01 "alien"
   42    4   14    1.53   61.98  0.1277E+02  0.8481E+01 "alien"
   39   30    7    1.51   59.25  0.1220E+02  0.9724E+01 "alien"
   42    1    1    1.84   51.50 -0.2710E+01  0.8841E+01 "alien"
   38   29    8    1.53   48.77  0.1005E+02  0.9747E+01 "alien"
   29   27   16    1.52   47.27  0.9738E+01  0.1243E+02 "alien"
   24   13   24    1.40   47.14 -0.3285E+01  0.7653E+01 "alien"
   28    8   20    1.56   45.90  0.9454E+01  0.9848E+01 "alien"
   50    1    2    1.54   44.43  0.9151E+01  0.9464E+01 "alien"
   47    1    6    1.59   42.22 -0.1466E+02  0.1312E+02 "alien"
   31    2   14    1.85   41.01 -0.2158E+01  0.7008E+01 "alien"
   28   21   12    1.82   40.18 -0.2114E+01  0.7828E+01 "alien"
   47    3   10    1.51   40.17  0.8275E+01  0.8979E+01 "alien"
   33   27   12    1.58   40.05 -0.1391E+02  0.1593E+02 "alien"
   32   28   14    1.52   39.43  0.8122E+01  0.1504E+02 "alien"
   35   24    4    1.79   38.25 -0.2012E+01  0.1075E+02 "alien"
   19    6   23    1.54   38.21  0.7871E+01  0.5134E+01 "alien"
   38   30    8    1.52   38.08  0.7844E+01  0.1634E+02 "alien"
   48    6   14    1.38   36.98 -0.2576E+01  0.1287E+02 "alien"
   48    1    2    1.60   36.51 -0.1268E+02  0.1636E+02 "alien"
   40    1    6    1.85   35.44 -0.1865E+01  0.8400E+01 "alien"
   21   19   16    1.81   34.37 -0.1809E+01  0.6789E+01 "alien"
   32   27   14    1.53   34.25  0.7056E+01  0.1251E+02 "alien"
   42    3   12    1.59   34.19 -0.1188E+02  0.1193E+02 "alien"
   35    5   17    1.58   33.06 -0.1148E+02  0.1127E+02 "alien"
   41    1    5    1.83   32.08 -0.1688E+01  0.8846E+01 "alien"
   21   17   20    1.60   30.67 -0.1065E+02  0.8241E+01 "alien"
   41    7   17    1.44   29.08 -0.2057E+02  0.1560E+02 "alien"
   54    5    8    1.37   28.99  0.4464E+02  0.1052E+02 "alien"
   43    1    0    1.80   28.22 -0.1264E+02  0.1213E+02 "alien"
   51    3    6    1.47   27.99 -0.2202E+02  0.1554E+02 "alien"
   38    0    9    1.84   27.69 -0.1181E+02  0.7807E+01 "alien"
   30   26   15    1.55   27.42  0.5647E+01  0.1161E+02 "alien"
   33   28   13    1.53   26.81  0.5522E+01  0.1332E+02 "alien"
   48    1    4    1.59   26.24 -0.9112E+01  0.8932E+01 "alien"
   37   30   10    1.50   25.91 -0.2038E+02  0.9533E+01 "alien"
   32    1   12    1.93   25.50 -0.2525E+02  0.1119E+02 "alien"
   40    1    4    1.89   25.38 -0.2513E+02  0.1412E+02 "alien"
   36   32   13    1.41   24.97 -0.1766E+02  0.1398E+02 "alien"
   36   28    9    1.58   24.61 -0.8547E+01  0.1488E+02 "alien"
   26   13   24    1.38   23.86 -0.1662E+01  0.9009E+01 "alien"
   32   11   22    1.39   23.45 -0.1633E+01  0.1016E+02 "alien"
   51    2    2    1.51   23.15  0.4768E+01  0.9275E+01 "alien"
   18   18   17    1.82   22.90 -0.6397E+01  0.6952E+01 "alien"
   53    7   13    1.30   22.31  0.5700E+02  0.1437E+02 "alien"
   27    9   21    1.52   22.26  0.4586E+01  0.6485E+01 "alien"
   37   31   12    1.43   22.16 -0.1568E+02  0.1433E+02 "alien"
   44    3   12    1.54   21.87  0.4504E+01  0.8676E+01 "alien"
   42    1    0    1.84   20.95 -0.4967E+01  0.1279E+02 "alien"
   39   34   11    1.37   20.61 -0.1436E+01  0.1030E+02 "alien"
   30   26   20    1.37   20.60 -0.1435E+01  0.1305E+02 "alien"
   39    0    8    1.83   20.13 -0.4461E+01  0.8262E+01 "alien"
   53    4    8    1.39   19.78 -0.1378E+01  0.9398E+01 "alien"
   46    3   11    1.51   19.69  0.4056E+01  0.8850E+01 "alien"
   36   33   14    1.37   19.61 -0.1366E+01  0.9758E+01 "alien"
   28   27   17    1.50   19.32 -0.1520E+02  0.8204E+01 "alien"
   21   10   22    1.55   19.30  0.3975E+01  0.7188E+01 "alien"
   34   27   11    1.59   18.46 -0.6410E+01  0.8803E+01 "alien"
   27    2   14    1.98   17.81 -0.1763E+02  0.1077E+02 "alien"
   31   31   17    1.39   17.69 -0.4276E+01  0.1580E+02 "alien"
   42    8   18    1.38   17.59 -0.1226E+01  0.8691E+01 "alien"
   38    6   17    1.51   17.53 -0.1379E+02  0.7765E+01 "alien"
   27   24   16    1.60   17.45 -0.6062E+01  0.1231E+02 "alien"
   40   32    8    1.44   16.91 -0.1196E+02  0.9905E+01 "alien"
   43    2   10    1.63   16.78  0.1974E+02  0.9444E+01 "alien"
   40   31    5    1.50   16.62 -0.1308E+02  0.1070E+02 "alien"
   38   33   12    1.39   16.39 -0.1142E+01  0.9910E+01 "alien"
   24    9   21    1.57   16.02 -0.5566E+01  0.7231E+01 "alien"
   25   22   20    1.49   15.49 -0.1219E+02  0.1303E+02 "alien"
   39    2    3    1.96   15.12 -0.1497E+02  0.1051E+02 "alien"
   44    4   13    1.51   15.11  0.3113E+01  0.8438E+01 "alien"
   34   22    2    1.90   14.73 -0.1458E+02  0.9163E+01 "alien"
   30    7   19    1.58   14.50 -0.5036E+01  0.9572E+01 "alien"
   37   26    5    1.67   14.46  0.1701E+02  0.9262E+01 "alien"
   50    3    7    1.49   14.41 -0.1134E+02  0.1754E+02 "alien"
   26   24   19    1.49   14.17 -0.1115E+02  0.8837E+01 "alien"
   22    9   21    1.60   14.07 -0.4886E+01  0.6775E+01 "alien"
   34    5   17    1.60   14.07 -0.4885E+01  0.7295E+01 "alien"
   41   31    1    1.50   13.99 -0.1101E+02  0.9628E+01 "alien"
   52    3    6    1.45   13.95 -0.9870E+01  0.9158E+01 "alien"
   33    7   19    1.52   13.86  0.2855E+01  0.7300E+01 "alien"
   22   19   15    1.85   13.72 -0.7220E+00  0.1094E+02 "alien"
   20   17   19    1.67   13.69  0.1611E+02  0.9369E+01 "alien"
   52    2    2    1.48   13.59 -0.1069E+02  0.9182E+01 "alien"
   28   24   15    1.63   13.49 -0.4686E+01  0.7916E+01 "alien"
   50    5   12    1.39   13.15 -0.9162E+00  0.9603E+01 "alien"
   49    7   15    1.34   12.95  0.1994E+02  0.1078E+02 "alien"
   49    4   11    1.43   12.89 -0.9122E+01  0.8861E+01 "alien"
   39   31    9    1.46   12.89 -0.1014E+02  0.1045E+02 "alien"
   52    4    8    1.42   12.75 -0.9018E+01  0.9413E+01 "alien"
   49    6   13    1.39   12.72 -0.8861E+00  0.1397E+02 "alien"
   53    2    0    1.46   12.44 -0.2316E+02  0.1509E+02 "alien"
   15   12   23    1.55   12.43  0.2561E+01  0.5145E+01 "alien"
   37    1    7    1.95   12.32 -0.1220E+02  0.8385E+01 "alien"
   51    2    3    1.50   11.89 -0.9351E+01  0.9063E+01 "alien"
   32   29   15    1.47   11.81 -0.9290E+01  0.9116E+01 "alien"
   36   30   12    1.47   11.80 -0.9284E+01  0.9330E+01 "alien"
   31    5   17    1.67   11.70  0.1377E+02  0.1074E+02 "alien"
   30   12   23    1.38   11.44 -0.7971E+00  0.6888E+01 "alien"
   52    3    5    1.46   11.42 -0.8982E+01  0.9341E+01 "alien"
   50    1    0    1.55   11.36  0.5150E+01  0.1266E+02 "alien"
   50    6   13    1.36   11.24  0.1731E+02  0.9841E+01 "alien"
   23   17   21    1.52   11.22  0.2312E+01  0.8980E+01 "alien"
   18    7   22    1.61   11.12 -0.3862E+01  0.5740E+01 "alien"
   45    5   14    1.45   11.09 -0.7847E+01  0.8593E+01 "alien"
   38    4   15    1.59   11.03 -0.3831E+01  0.7818E+01 "alien"
   40   31    7    1.47   11.00 -0.8652E+01  0.9937E+01 "alien"
   31   28   15    1.50   10.99 -0.8645E+01  0.1306E+02 "alien"
   33   33   17    1.34   10.96  0.3618E+02  0.1490E+02 "alien"
   49    5   12    1.41   10.93 -0.7733E+01  0.9247E+01 "alien"
   53    3    4    1.44   10.92 -0.7726E+01  0.9536E+01 "alien"
   42   32    1    1.46   10.80 -0.8499E+01  0.1069E+02 "alien"
   33   25   10    1.68   10.72  0.1261E+02  0.9752E+01 "alien"
   40   31    6    1.49   10.66 -0.8388E+01  0.9755E+01 "alien"
   26   10   22    1.48   10.61 -0.8352E+01  0.8564E+01 "alien"
   22   16   21    1.55   10.57  0.2177E+01  0.6057E+01 "alien"
   25   18   12    1.98   10.42 -0.1032E+02  0.1101E+02 "alien"
   35   28   11    1.55   10.19  0.2099E+01  0.8888E+01 "alien"
   43    3   12    1.57   10.17 -0.3531E+01  0.7943E+01 "alien"
   46    2    8    1.59   10.02 -0.3480E+01  0.1633E+02 "alien"


 List of    157 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   22.31  0.5700E+02  0.1437E+02
   49    7   15    1.34   12.95  0.1994E+02  0.1078E+02
   33   33   17    1.34   10.96  0.3618E+02  0.1490E+02
   50    6   13    1.36   11.24  0.1731E+02  0.9841E+01
   54    5    8    1.37   28.99  0.4464E+02  0.1052E+02
   46    7   16    1.37   95.23 -0.6635E+01  0.8947E+01
   30   26   20    1.37   20.60 -0.1435E+01  0.1305E+02
   36   33   14    1.37   19.61 -0.1366E+01  0.9758E+01
   39   34   11    1.37   20.61 -0.1436E+01  0.1030E+02
   28   21   22    1.38  176.55 -0.1230E+02  0.7754E+01
   30   12   23    1.38   11.44 -0.7971E+00  0.6888E+01
   31   30   18    1.38  171.96 -0.1198E+02  0.9370E+01
   26   13   24    1.38   23.86 -0.1662E+01  0.9009E+01
   33   32   16    1.38  314.03 -0.2188E+02  0.1413E+02
   54    4    7    1.38   72.92 -0.5080E+01  0.9689E+01
   48    6   14    1.38   36.98 -0.2576E+01  0.1287E+02
   38    9   20    1.38   67.02 -0.4669E+01  0.1164E+02
   42    8   18    1.38   17.59 -0.1226E+01  0.8691E+01
   49    6   13    1.39   12.72 -0.8861E+00  0.1397E+02
   38   33   12    1.39   16.39 -0.1142E+01  0.9910E+01
   36   32   14    1.39  161.49 -0.1125E+02  0.1578E+02
   50    5   12    1.39   13.15 -0.9162E+00  0.9603E+01
   53    4    8    1.39   19.78 -0.1378E+01  0.9398E+01
   32   11   22    1.39   23.45 -0.1633E+01  0.1016E+02
   31   31   17    1.39   17.69 -0.4276E+01  0.1580E+02
   54    4    6    1.39  109.76 -0.7647E+01  0.9595E+01
   27   20   22    1.40  100.25 -0.6985E+01  0.7682E+01
   24   13   24    1.40   47.14 -0.3285E+01  0.7653E+01
   28   12   23    1.40  106.75 -0.7438E+01  0.9466E+01
   51    4   10    1.41  104.31 -0.7268E+01  0.1374E+02
   34   31   15    1.41   75.07 -0.5231E+01  0.9255E+01
   40   33    9    1.41   70.69 -0.4925E+01  0.9965E+01
   39   33   10    1.41  123.30 -0.8590E+01  0.1030E+02
   50    5   11    1.41  133.61 -0.9309E+01  0.1365E+02
   49    5   12    1.41   10.93 -0.7733E+01  0.9247E+01
   36   32   13    1.41   24.97 -0.1766E+02  0.1398E+02
   52    4    8    1.42   12.75 -0.9018E+01  0.9413E+01
   37   31   12    1.43   22.16 -0.1568E+02  0.1433E+02
   49    4   11    1.43   12.89 -0.9122E+01  0.8861E+01
   41    7   17    1.44   29.08 -0.2057E+02  0.1560E+02
   53    3    4    1.44   10.92 -0.7726E+01  0.9536E+01
   40   32    8    1.44   16.91 -0.1196E+02  0.9905E+01
   52    3    6    1.45   13.95 -0.9870E+01  0.9158E+01
   45    5   14    1.45   11.09 -0.7847E+01  0.8593E+01
   53    2    0    1.46   12.44 -0.2316E+02  0.1509E+02
   52    3    5    1.46   11.42 -0.8982E+01  0.9341E+01
   39   31    9    1.46   12.89 -0.1014E+02  0.1045E+02
   42   32    1    1.46   10.80 -0.8499E+01  0.1069E+02
   36   30   12    1.47   11.80 -0.9284E+01  0.9330E+01
   32   29   15    1.47   11.81 -0.9290E+01  0.9116E+01
   40   31    7    1.47   11.00 -0.8652E+01  0.9937E+01
   51    3    6    1.47   27.99 -0.2202E+02  0.1554E+02
   52    2    2    1.48   13.59 -0.1069E+02  0.9182E+01
   26   10   22    1.48   10.61 -0.8352E+01  0.8564E+01
   25   22   20    1.49   15.49 -0.1219E+02  0.1303E+02
   40   31    6    1.49   10.66 -0.8388E+01  0.9755E+01
   50    3    7    1.49   14.41 -0.1134E+02  0.1754E+02
   26   24   19    1.49   14.17 -0.1115E+02  0.8837E+01
   37   30   10    1.50   25.91 -0.2038E+02  0.9533E+01
   28   27   17    1.50   19.32 -0.1520E+02  0.8204E+01
   40   31    5    1.50   16.62 -0.1308E+02  0.1070E+02
   41   31    1    1.50   13.99 -0.1101E+02  0.9628E+01
   31   28   15    1.50   10.99 -0.8645E+01  0.1306E+02
   51    2    3    1.50   11.89 -0.9351E+01  0.9063E+01
   38    6   17    1.51   17.53 -0.1379E+02  0.7765E+01
   44    4   13    1.51   15.11  0.3113E+01  0.8438E+01
   51    2    2    1.51   23.15  0.4768E+01  0.9275E+01
   47    3   10    1.51   40.17  0.8275E+01  0.8979E+01
   39   30    7    1.51   59.25  0.1220E+02  0.9724E+01
   51    2    1    1.51   94.64  0.1949E+02  0.1124E+02
   46    3   11    1.51   19.69  0.4056E+01  0.8850E+01
   41    5   15    1.51   81.32  0.1675E+02  0.1232E+02
   49    2    7    1.52   80.98  0.1668E+02  0.9673E+01
   38   30    8    1.52   38.08  0.7844E+01  0.1634E+02
   29   27   16    1.52   47.27  0.9738E+01  0.1243E+02
   32   28   14    1.52   39.43  0.8122E+01  0.1504E+02
   33    7   19    1.52   13.86  0.2855E+01  0.7300E+01
   23   17   21    1.52   11.22  0.2312E+01  0.8980E+01
   27    9   21    1.52   22.26  0.4586E+01  0.6485E+01
   42    4   14    1.53   61.98  0.1277E+02  0.8481E+01
   33   28   13    1.53   26.81  0.5522E+01  0.1332E+02
   38   29    8    1.53   48.77  0.1005E+02  0.9747E+01
   47    2    9    1.53   88.83  0.1830E+02  0.1607E+02
   32   27   14    1.53   34.25  0.7056E+01  0.1251E+02
   37   29    9    1.54   82.16  0.1692E+02  0.9680E+01
   44    3   12    1.54   21.87  0.4504E+01  0.8676E+01
   50    1    2    1.54   44.43  0.9151E+01  0.9464E+01
   19    6   23    1.54   38.21  0.7871E+01  0.5134E+01
   50    1    0    1.55   11.36  0.5150E+01  0.1266E+02
   22   16   21    1.55   10.57  0.2177E+01  0.6057E+01
   35   28   11    1.55   10.19  0.2099E+01  0.8888E+01
   15   12   23    1.55   12.43  0.2561E+01  0.5145E+01
   30   26   15    1.55   27.42  0.5647E+01  0.1161E+02
   21   10   22    1.55   19.30  0.3975E+01  0.7188E+01
   28    8   20    1.56   45.90  0.9454E+01  0.9848E+01
   43    3   12    1.57   10.17 -0.3531E+01  0.7943E+01
   24    9   21    1.57   16.02 -0.5566E+01  0.7231E+01
   35    5   17    1.58   33.06 -0.1148E+02  0.1127E+02
   36   28    9    1.58   24.61 -0.8547E+01  0.1488E+02
   33   27   12    1.58   40.05 -0.1391E+02  0.1593E+02
   30    7   19    1.58   14.50 -0.5036E+01  0.9572E+01
   46    2    8    1.59   10.02 -0.3480E+01  0.1633E+02
   48    1    4    1.59   26.24 -0.9112E+01  0.8932E+01
   34   27   11    1.59   18.46 -0.6410E+01  0.8803E+01
   38    4   15    1.59   11.03 -0.3831E+01  0.7818E+01
   47    1    6    1.59   42.22 -0.1466E+02  0.1312E+02
   42    3   12    1.59   34.19 -0.1188E+02  0.1193E+02
   34    5   17    1.60   14.07 -0.4885E+01  0.7295E+01
   27   24   16    1.60   17.45 -0.6062E+01  0.1231E+02
   21   17   20    1.60   30.67 -0.1065E+02  0.8241E+01
   22    9   21    1.60   14.07 -0.4886E+01  0.6775E+01
   48    1    2    1.60   36.51 -0.1268E+02  0.1636E+02
   18    7   22    1.61   11.12 -0.3862E+01  0.5740E+01
   46    1    6    1.63   69.12 -0.2401E+02  0.1251E+02
   28   24   15    1.63   13.49 -0.4686E+01  0.7916E+01
   43    2   10    1.63   16.78  0.1974E+02  0.9444E+01
   37   26    5    1.67   14.46  0.1701E+02  0.9262E+01
   31    5   17    1.67   11.70  0.1377E+02  0.1074E+02
   20   17   19    1.67   13.69  0.1611E+02  0.9369E+01
   33   25   10    1.68   10.72  0.1261E+02  0.9752E+01
   15    7   20    1.79  111.59 -0.5872E+01  0.4536E+01
   43    1    2    1.79  179.18 -0.9428E+01  0.8334E+01
   35   24    4    1.79   38.25 -0.2012E+01  0.1075E+02
   43    1    1    1.79   83.83 -0.4411E+01  0.8422E+01
   43    1    0    1.80   28.22 -0.1264E+02  0.1213E+02
   31   22   10    1.80  329.95 -0.1736E+02  0.1287E+02
   42    1    4    1.81   96.66 -0.5086E+01  0.8615E+01
   21   19   16    1.81   34.37 -0.1809E+01  0.6789E+01
   14    5   20    1.81   65.87 -0.3466E+01  0.3942E+01
   33   23    7    1.81  206.18 -0.1085E+02  0.1326E+02
   17   15   18    1.82   87.97 -0.4629E+01  0.5150E+01
   30    3   15    1.82   78.96 -0.4155E+01  0.9159E+01
   18   18   17    1.82   22.90 -0.6397E+01  0.6952E+01
   42    1    3    1.82  132.11 -0.6951E+01  0.8678E+01
   28   21   12    1.82   40.18 -0.2114E+01  0.7828E+01
   41    1    5    1.83   32.08 -0.1688E+01  0.8846E+01
   39    0    8    1.83   20.13 -0.4461E+01  0.8262E+01
   42    1    1    1.84   51.50 -0.2710E+01  0.8841E+01
   38    0    9    1.84   27.69 -0.1181E+02  0.7807E+01
   29   21   11    1.84  110.15 -0.5796E+01  0.8025E+01
   42    1    0    1.84   20.95 -0.4967E+01  0.1279E+02
   31    2   14    1.85   41.01 -0.2158E+01  0.7008E+01
   40    1    6    1.85   35.44 -0.1865E+01  0.8400E+01
   22   19   15    1.85   13.72 -0.7220E+00  0.1094E+02
   30   21   10    1.85  110.29 -0.5803E+01  0.1193E+02
   34    1   12    1.86  111.11 -0.5847E+01  0.7492E+01
   41    1    3    1.86  128.65 -0.6770E+01  0.8960E+01
   33   22    6    1.86  121.17 -0.6376E+01  0.8686E+01
   27   20   12    1.87   62.29 -0.3278E+01  0.7695E+01
   34   22    4    1.87  259.63 -0.1366E+02  0.9456E+01
   40    1    4    1.89   25.38 -0.2513E+02  0.1412E+02
   34   22    2    1.90   14.73 -0.1458E+02  0.9163E+01
   32    1   12    1.93   25.50 -0.2525E+02  0.1119E+02
   37    1    7    1.95   12.32 -0.1220E+02  0.8385E+01
   39    2    3    1.96   15.12 -0.1497E+02  0.1051E+02
   25   18   12    1.98   10.42 -0.1032E+02  0.1101E+02
   27    2   14    1.98   17.81 -0.1763E+02  0.1077E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        89
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1558
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1378
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      180


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
