
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)   7-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.307    77.307    38.685  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/demo_undamaged/myfake_???.img SMV
 DATA_RANGE=       1      50
 ROTATION_AXIS=  0.000969 -0.000005  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       1
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000805  0.000651
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
 BEAM_DIVERGENCE_E.S.D.=   0.040
 REFLECTING_RANGE_E.S.D.=   0.090
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

   91392 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       1 ...      50
       0 OVERLOADED REFLECTIONS (IGNORED)
      81 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
   91311 REFLECTIONS ACCEPTED



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
 UNIT_CELL_CONSTANTS=    77.00    77.00    38.00  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_B-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_C-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 REIDX=  0  0  0  0  0  0  0  0  0  0  0  0 (all 0 : not specified)
 TEST_RESOLUTION_RANGE=    8.00    4.50 (Angstrom) for space group determination
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas             2177
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

 *  44        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  20        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    27        mC        500.0     159.4   38.7  109.3  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.7  109.3  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.3  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    30        mC        998.4      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.4      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    41        mC        998.4     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    22        hP        998.4      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.4   86.4  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.3   86.4   86.4  78.4  50.8  50.8    1  1  0  0  1  0 -1  0  0  1 -1  0
     7        tI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    1  0 -1  0  0  1 -1  0  1  1  0  0
     8        oI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    0 -1  1  0 -1  0  1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.4  247.5 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.4   86.4   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.1  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.4   86.4  177.1 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.4   86.4   86.4  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.1  109.3   38.7  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.7   77.3  177.1 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.4  177.1   77.3 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    5.0    1647    0.02  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.00    77.00    38.00  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.670   -23.433     0.022
 UNIT_CELL_B-AXIS=    23.433    73.670    -0.034
 UNIT_CELL_C-AXIS=    -0.005     0.020    38.685
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM     13500 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM AXIS ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.12
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.197    77.197    38.634  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.7E-02 1.7E-02 8.2E-03 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012954  0.012954  0.025884  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.567   -23.394     0.028
 COORDINATES OF UNIT CELL B-AXIS    23.394    73.567     0.009
 COORDINATES OF UNIT CELL C-AXIS    -0.015     0.000    38.634
 CRYSTAL MOSAICITY (DEGREES)     0.090
 LAB COORDINATES OF ROTATION AXIS -0.000361  0.000012  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006291  0.002034
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1260.91   1242.03
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.01
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.20    77.20    38.63  90.000  90.000  90.000


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
 -0.05  0.05      8900        524.896       0.907
  0.05  0.15      9159        499.581       0.953
  0.15  0.25      9194        454.389       1.048
  0.25  0.35      9118        449.225       1.060
  0.35  0.45      9106        462.643       1.029
  0.45  0.55      9263        475.769       1.001
  0.55  0.65      9023        463.643       1.027
  0.65  0.75      9168        473.212       1.006
  0.75  0.85      9270        463.212       1.028
  0.85  0.95      9110        495.161       0.961
 
 -0.05  0.95     91311        476.038       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED      200
 DEGREES OF FREEDOM OF CHI^2 FIT             24123.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.644
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     0.1 XMAX=    49.9 NXBIN=   10
 YMIN= 0.00067 YMAX= 0.60134 NYBIN=   20
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS      38758


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED      756
 DEGREES OF FREEDOM OF CHI^2 FIT             24094.8
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.641
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     5.9 XMAX=  2457.3 NXBIN=   27
 YMIN=     5.7 YMAX=  2521.3 NYBIN=   28
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS      38758


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED      130
 DEGREES OF FREEDOM OF CHI^2 FIT             24123.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.639
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     0.1 XMAX=    49.9 NXBIN=   10
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1648    1691
 DETECTOR_SURFACE_POSITION=     815    1691
 DETECTOR_SURFACE_POSITION=     815     836
 DETECTOR_SURFACE_POSITION=    1648     836
 DETECTOR_SURFACE_POSITION=    2174    1664
 DETECTOR_SURFACE_POSITION=    1622    2230
 DETECTOR_SURFACE_POSITION=     841    2230
 DETECTOR_SURFACE_POSITION=     289    1664
 DETECTOR_SURFACE_POSITION=     289     863
 DETECTOR_SURFACE_POSITION=     841     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2174     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS      38758


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
 8.195E-01  1.787E-03   26.13


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

   38.598  12.237    25.92   1.25      3.02      2.54      90      93      29
   12.237   8.879    24.71   1.51      3.20      2.67     194     196      45
    8.879   7.314    24.55   1.58      3.12      2.58     215     221      36
    7.314   6.363    24.19   1.26      2.92      2.66     269     271      52
    6.363   5.707    24.04   1.40      3.04      2.65     295     296      51
    5.707   5.219    24.23   1.15      2.39      2.70     330     334      56
    5.219   4.838    24.13   1.11      2.61      2.70     364     369      60
    4.838   4.530    24.30   1.11      2.48      2.72     410     414      67
    4.530   4.274    24.80   0.92      2.09      2.72     476     481      49
    4.274   4.058    24.36   1.03      2.47      2.70     488     492      59
    4.058   3.871    23.98   1.02      2.25      2.70     498     503      63
    3.871   3.707    23.38   1.17      2.79      2.66     463     476      56
    3.707   3.563    24.40   0.95      2.13      2.66     560     565      43
    3.563   3.435    22.37   0.98      2.30      2.68     499     505      65
    3.435   3.319    23.84   0.99      2.38      2.70     528     535      63
    3.319   3.214    23.44   1.11      2.34      2.66     478     485      60
    3.214   3.119    21.63   1.19      2.65      2.76     612     628      57
    3.119   3.032    22.04   1.25      2.75      2.76     595     605     102
    3.032   2.951    22.19   1.11      2.37      2.72     553     563      84
    2.951   2.877    20.05   1.17      2.77      2.85     637     648      85
    2.877   2.808    21.41   1.10      2.71      2.84     679     686      97
    2.808   2.744    20.78   1.26      3.00      2.87     627     635     122
    2.744   2.684    20.62   1.30      2.85      2.87     715     726      78
    2.684   2.628    20.42   1.25      3.00      2.91     736     746      78
    2.628   2.575    19.68   1.15      2.84      2.94     696     704     120
    2.575   2.525    19.90   1.31      3.17      2.88     662     676     118
    2.525   2.478    10.08   1.21      3.84      3.51     784     788      61
    2.478   2.434    -0.00   0.88   5646.15   6037.12     913     917       1
    2.434   2.391     0.01   0.96   5164.33   5165.57     837     843       1
    2.391   2.351    -0.04   0.94    -99.90    -99.90     937     942       2
    2.351   2.313     0.02   0.96   2182.67   2243.80     832     836       2
    2.313   2.277    -0.03   0.99    -99.90    -99.90     840     842       1
    2.277   2.242    -0.03   0.97    -99.90    -99.90     967     972       2
    2.242   2.209    -0.01   0.93    -99.90    -99.90     857     864       7
    2.209   2.177     0.01   0.96   2225.90   2288.57    1036    1040       3
    2.177   2.147     0.10   0.83    674.51    741.40     964     969       2
    2.147   2.118     0.05   0.89   1280.15   1321.02     966     978       5
    2.118   2.090    -0.01   0.89    -99.90    -99.90     995    1005       1
    2.090   2.063     0.05   0.91   1114.86   1161.51    1012    1017       2
    2.063   2.037     0.03   0.86   2457.30   2643.17    1030    1038       5
    2.037   2.012     0.00   0.97  13867.04  13938.65     962     980       5
    2.012   1.988    -0.01   0.90    -99.90    -99.90    1030    1038       0
    1.988   1.965    -0.04   0.88    -99.90    -99.90     919     940       5
    1.965   1.943     0.00   0.94   2337.55   2400.59    1017    1031       4
    1.943   1.921    -0.04   0.95    -99.90    -99.90    1074    1089       3
    1.921   1.900    -0.02   0.87    -99.90    -99.90    1133    1147       2
    1.900   1.880    -0.00   0.99  12283.74  12031.65    1061    1068       5
    1.880   1.860    -0.05   0.84    -99.90    -99.90    1113    1120       5
    1.860   1.841     0.02   0.91   2487.28   2591.48     944     955       5
    1.841   1.823    -0.04   0.97    -99.90    -99.90    1158    1160       5
    1.823   1.805    -0.00   0.86   9025.76   9919.71    1101    1108       2
    1.805   1.787    -0.03   0.94    -99.90    -99.90    1136    1138       0
    1.787   1.770    -0.02   0.95    -99.90    -99.90    1186    1189       2
    1.770   1.754     0.01   0.87   3730.56   3969.67    1123    1130       2
    1.754   1.738     0.01   0.84   4469.93   4948.30    1289    1292       4
    1.738   1.722    -0.00   0.97  11840.79  12092.51    1300    1303       4
    1.722   1.707    -0.05   0.85    -99.90    -99.90    1112    1119       1
    1.707   1.693    -0.06   0.95    -99.90    -99.90    1176    1185       2
    1.693   1.678     0.00   0.88  25308.28  26344.98    1245    1247       3
    1.678   1.664    -0.06   0.98    -99.90    -99.90    1320    1331       3
    1.664   1.651    -0.03   0.92    -99.90    -99.90    1192    1205       4
    1.651   1.637    -0.06   0.87    -99.90    -99.90    1153    1161       7
    1.637   1.624    -0.03   1.02    -99.90    -99.90    1219    1222       5
    1.624   1.611     0.06   0.92    990.33   1017.36    1294    1312       1
    1.611   1.599     0.05   0.88    974.57   1018.93    1266    1280       2
    1.599   1.587    -0.01   0.87    -99.90    -99.90    1182    1187       4
    1.587   1.575    -0.02   0.99    -99.90    -99.90    1229    1236       1
    1.575   1.563    -0.01   0.96  26965.10  27186.84    1202    1210       2
    1.563   1.552    -0.00   0.99  27696.00  27420.03    1291    1294       1
    1.552   1.541    -0.04   0.88    -99.90    -99.90    1236    1240       0
    1.541   1.530    -0.01   0.99    -99.90    -99.90    1382    1382       1
    1.530   1.519    -0.01   0.92    -99.90    -99.90    1343    1346       2
    1.519   1.509    -0.01   0.98    -99.90    -99.90    1417    1420       4
    1.509   1.499     0.01   0.91   5564.90   5751.74    1293    1307       5
    1.499   1.489     0.01   0.89   3405.71   3616.20    1311    1317       3
    1.489   1.479     0.02   1.03   2908.86   2901.20    1457    1461       4
    1.479   1.469    -0.00   0.93   5391.93   5567.99    1467    1469       0
    1.469   1.460    -0.02   0.90    -99.90    -99.90    1373    1373       3
    1.460   1.451    -0.05   0.96    -99.90    -99.90    1340    1352       1
    1.451   1.442    -0.01   0.93  21339.38  22018.11    1306    1328       1
    1.442   1.433    -0.03   1.01    -99.90    -99.90    1238    1257       4
    1.433   1.424    -0.04   0.91    -99.90    -99.90    1228    1245       5
    1.424   1.415    -0.01   0.96  24274.31  24592.94    1180    1199       1
    1.415   1.407    -0.01   0.98    -99.90    -99.90    1228    1251       7
    1.407   1.399    -0.02   0.98   1082.77   1086.30    1044    1064      14
    1.399   1.390     0.08   0.96    313.56    318.38     909     914       6
    1.390   1.382    -0.07   0.95    -99.90    -99.90     980     981      13
    1.382   1.375    -0.01   0.97   5399.62   5429.84     952     959       9
    1.375   1.367    -0.00   0.90   1916.04   1971.40     805     815       8
    1.367   1.359    -0.02   1.00    -99.90    -99.90     775     779      14
    1.359   1.352     0.04   0.90   1165.82   1224.72     720     726      12
    1.352   1.344    -0.03   0.87    -99.90    -99.90     653     678       4
    1.344   1.337    -0.01   0.88   2500.08   2693.90     537     571      10
    1.337   1.330    -0.01   0.93   2866.68   2993.74     494     535       3
    1.330   1.323     0.09   0.89    555.71    585.39     420     462      20
    1.323   1.316    -0.00   0.81    -99.90    -99.90     373     388       5
    1.316   1.309     0.12   0.93    591.74    580.96     231     260      11
    1.309   1.303     0.01   0.69    -99.90    -99.90     136     183       0
    1.303   1.296     0.03   0.91    -99.90    -99.90     100     150       4
    1.296   1.290     0.01   0.84    -99.90    -99.90      68     117       4
   --------------------------------------------------------------------------
   38.598   1.290     3.31   0.96      4.14      4.31   88062   89110    2157
 


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
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1      11
 ******************************************************************************


 AVERAGE INTENSITY FOR      13 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -0.3% OF MEAN INTENSITY


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

     3.85         832     489      1260       38.8%       2.4%      2.2%      674   31.11      3.3%    99.8*    49    1.194       4
     2.73        1417     788      2111       37.3%       2.5%      2.3%     1185   28.60      3.4%    99.7*    75*   1.153      32
     2.23        2038     974      2693       36.2%       4.2%      3.6%     1925    9.76      5.6%    99.8*    24    0.789      95
     1.93        2439    1214      3113       39.0%    1706.9%   1783.4%     2305    0.03   2341.4%    -2.3     -1    0.448      69
     1.73        2761    1402      3521       39.8%    2752.7%   2908.7%     2629    0.01   3840.4%     3.5    -31    0.406      41
     1.58        3014    1567      3881       40.4%     -99.9%    -99.9%     2831    0.00    -99.9%    -3.1    -19    0.426      27
     1.46        3270    1677      4179       40.1%    5607.8%   5766.0%     3095   -0.00   7846.8%    -4.4    -19    0.502      42
     1.37        2762    1398      4492       31.1%     -99.9%    -99.9%     2642    0.00    -99.9%    -1.4    -10    0.577      40
     1.29        1068     671      4733       14.2%     -99.9%    -99.9%      782    0.00    -99.9%    18.7*     4    0.563       6
    total       19601   10180     29983       34.0%       3.8%      3.7%    18068    4.64      5.3%    99.9*     4    0.627     356


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     20049
 NUMBER OF REJECTED MISFITS                               411
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                   13
 NUMBER OF ACCEPTED OBSERVATIONS                        19625
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                  10195



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1      21
 ******************************************************************************


 AVERAGE INTENSITY FOR      23 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -0.2% OF MEAN INTENSITY


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

     3.85        1583     658      1260       52.2%       2.4%      2.5%     1457   36.17      3.0%    99.8*    80*   1.438     111
     2.73        2705    1214      2111       57.5%       2.5%      2.5%     2441   31.79      3.3%    99.8*    82*   1.416     217
     2.23        3896    1665      2693       61.8%       4.2%      3.7%     3743   10.27      5.4%    99.8*    22*   0.799     329
     1.93        4633    2026      3113       65.1%   32992.6%  34632.2%     4452   -0.01  43385.6%     0.2     -2    0.470     357
     1.73        5276    2313      3521       65.7%     -99.9%    -99.9%     5111    0.00    -99.9%     5.6      3    0.497     368
     1.58        5776    2547      3881       65.6%     -99.9%    -99.9%     5552    0.00    -99.9%    -2.4     -4    0.489     408
     1.46        6250    2717      4179       65.0%   49469.7%  49896.4%     6043   -0.01  65160.4%    -3.5     -5    0.523     450
     1.37        5269    2307      4492       51.4%     -99.9%    -99.9%     5112    0.00    -99.9%    -2.3     -3    0.531     358
     1.29        2065    1147      4733       24.2%    1344.6%   1357.5%     1646    0.07   1800.2%    10.6    -17    0.468      53
    total       37453   16594     29983       55.3%       3.9%      3.9%    35557    4.78      5.0%    99.9*     8    0.653    2651


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     38387
 NUMBER OF REJECTED MISFITS                               853
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                   23
 NUMBER OF ACCEPTED OBSERVATIONS                        37511
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                  16624



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1      30
 ******************************************************************************


 AVERAGE INTENSITY FOR      35 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -0.1% OF MEAN INTENSITY


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

     3.85        2273     684      1260       54.3%       2.5%      2.6%     2209   42.45      3.0%    99.8*    81*   1.520     215
     2.73        3868    1282      2111       60.7%       2.5%      2.7%     3734   36.97      3.0%    99.9*    83*   1.396     474
     2.23        5571    1759      2693       65.3%       4.3%      4.0%     5486   11.94      5.1%    99.9*    21*   0.804     897
     1.93        6618    2116      3113       68.0%    4612.9%   4827.3%     6535   -0.00   5582.8%     2.1     -2    0.501    1143
     1.73        7555    2419      3521       68.7%   42778.2%  44571.2%     7473   -0.01  51824.2%     5.5      0    0.500    1270
     1.58        8248    2659      3881       68.5%     -99.9%    -99.9%     8116    0.00    -99.9%    -0.9      0    0.488    1372
     1.46        8931    2834      4179       67.8%    6040.7%   6199.6%     8820   -0.01   7309.5%     0.5      0    0.527    1479
     1.37        7550    2553      4492       56.8%   38737.5%  39369.2%     7415   -0.03  47528.4%    -2.3     -1    0.511    1118
     1.29        2954    1344      4733       28.4%    1541.3%   1591.4%     2548    0.08   1948.0%    12.5*     2    0.474     207
    total       53568   17650     29983       58.9%       4.0%      4.2%    52336    5.51      4.8%    99.9*     7    0.616    8175


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     54882
 NUMBER OF REJECTED MISFITS                              1212
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                   35
 NUMBER OF ACCEPTED OBSERVATIONS                        53635
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                  17681



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1      40
 ******************************************************************************


 AVERAGE INTENSITY FOR      44 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     0.1% OF MEAN INTENSITY


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

     3.85        3019     781      1260       62.0%       2.6%      2.7%     2984   46.04      2.9%    99.9*    79*   1.551     320
     2.73        5159    1452      2111       68.8%       2.5%      2.7%     5096   40.25      2.9%    99.9*    82*   1.496     660
     2.23        7402    2019      2693       75.0%       4.2%      4.1%     7350   12.83      5.0%    99.9*    21*   0.820    1206
     1.93        8834    2460      3113       79.0%    4686.9%   4918.3%     8767    0.01   5514.5%     4.8      3    0.512    1571
     1.73       10081    2824      3521       80.2%     -99.9%    -99.9%    10017    0.00    -99.9%     4.9      0    0.493    1810
     1.58       11027    3136      3881       80.8%     -99.9%    -99.9%    10926    0.00    -99.9%     0.1      0    0.493    2034
     1.46       11957    3396      4179       81.3%   22615.3%  23138.0%    11885   -0.02  26796.9%     0.6     -2    0.518    2194
     1.37       10050    3124      4492       69.5%    3468.8%   3525.9%     9884   -0.01   4168.9%    72.6*     2    0.511    1625
     1.29        3959    1585      4733       33.5%    2766.8%   2872.2%     3594    0.05   3430.7%     9.6      2    0.487     351
    total       71488   20777     29983       69.3%       4.1%      4.3%    70503    5.78      4.8%    99.9*     8    0.621   11771


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     73275
 NUMBER OF REJECTED MISFITS                              1653
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                   44
 NUMBER OF ACCEPTED OBSERVATIONS                        71578
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                  20814



 ******************************************************************************
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       1      50)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

    2    0    0    38.598   0.2304E+04  0.8847E+02        26.04           1 
    3    0    0    25.732   0.4242E+00  0.1020E+01         0.42           2*
    4    0    0    19.299   0.1170E+05  0.3173E+03        36.87           2 
    5    0    0    15.439   0.4030E+00  0.1680E+01         0.24           2*
    6    0    0    12.866   0.3431E+05  0.9311E+03        36.85           2 
    7    0    0    11.028   0.3012E+01  0.1759E+01         1.71           2*
    8    0    0     9.650   0.4361E+03  0.1301E+02        33.53           2 
    9    0    0     8.577  -0.1964E+01  0.2304E+01        -0.85           2*
   10    0    0     7.720   0.2408E+04  0.9374E+02        25.69           1 
   11    0    0     7.018  -0.1009E+01  0.4038E+01        -0.25           2*
   12    0    0     6.433   0.1074E+03  0.5299E+01        20.27           2 
   13    0    0     5.938  -0.1304E+02  0.6860E+01        -1.90           1*
   14    0    0     5.514   0.3890E+04  0.1069E+03        36.40           2 
   15    0    0     5.146  -0.2013E+01  0.4716E+01        -0.43           2*
   16    0    0     4.825   0.3393E+04  0.9387E+02        36.14           2 
   17    0    0     4.541   0.3329E+01  0.5478E+01         0.61           2*
   18    0    0     4.289   0.1928E+05  0.5237E+03        36.82           2 
   19    0    0     4.063   0.5510E+01  0.7962E+01         0.69           2*
   20    0    0     3.860   0.1032E+05  0.3982E+03        25.93           1 
   21    0    0     3.676  -0.3527E+01  0.7289E+01        -0.48           2*
   22    0    0     3.509   0.2733E+04  0.7664E+02        35.65           2 
   23    0    0     3.356   0.8654E+00  0.8922E+01         0.10           1*
   24    0    0     3.217   0.2091E+04  0.5991E+02        34.90           2 
   25    0    0     3.088   0.2983E+01  0.6761E+01         0.44           2*
   26    0    0     2.969   0.5925E+04  0.1632E+03        36.32           2 
   27    0    0     2.859  -0.8872E+00  0.9900E+01        -0.09           2*
   28    0    0     2.757   0.1340E+04  0.3987E+02        33.60           2 
   29    0    0     2.662  -0.8360E+01  0.9222E+01        -0.91           2*
   30    0    0     2.573   0.1420E+02  0.1351E+02         1.05           1 
   31    0    0     2.490   0.1156E+02  0.7500E+01         1.54           2*
   32    0    0     2.412   0.8101E+01  0.1060E+02         0.76           2 
   33    0    0     2.339   0.4553E+01  0.1086E+02         0.42           1*
   34    0    0     2.270  -0.1053E+02  0.7795E+01        -1.35           2 
   35    0    0     2.206  -0.1456E+02  0.1127E+02        -1.29           2*
   36    0    0     2.144  -0.1581E+02  0.8123E+01        -1.95           2 
   37    0    0     2.086   0.1066E+01  0.1146E+02         0.09           2*
   38    0    0     2.031  -0.1256E+02  0.9859E+01        -1.27           2 
   39    0    0     1.979  -0.5599E+01  0.1141E+02        -0.49           1*
   40    0    0     1.930   0.1747E+02  0.1134E+02         1.54           2 
   41    0    0     1.883   0.6981E+00  0.7791E+01         0.09           2*
   42    0    0     1.838   0.2156E+00  0.1021E+02         0.02           2 
   43    0    0     1.795   0.1324E+02  0.1058E+02         1.25           2*
   44    0    0     1.754  -0.7110E+01  0.7705E+01        -0.92           2 
   45    0    0     1.715  -0.6409E+01  0.1082E+02        -0.59           2*
   46    0    0     1.678  -0.7450E+01  0.7781E+01        -0.96           2 
   47    0    0     1.642   0.3224E+00  0.7934E+01         0.04           2*
   48    0    0     1.608   0.5122E+01  0.1581E+02         0.32           1 
   49    0    0     1.575   0.1359E+02  0.1166E+02         1.17           1*
   50    0    0     1.544   0.1444E+02  0.1612E+02         0.90           1 
   51    0    0     1.514   0.1233E+02  0.1565E+02         0.79           1*
   52    0    0     1.485  -0.7448E+00  0.1182E+02        -0.06           1 

 AVERAGE INTENSITY FOR      44 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     0.1% OF MEAN INTENSITY


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

     3.85        3733     971      1260       77.1%       2.6%      2.7%     3690   45.86      3.0%    99.8*    79*   1.620     366
     2.73        6403    1769      2111       83.8%       2.5%      2.7%     6304   40.60      2.9%    99.9*    83*   1.590     714
     2.23        9174    2415      2693       89.7%       4.2%      4.1%     9100   13.04      4.9%    99.9*    21*   0.829    1323
     1.93       11050    2922      3113       93.9%    3545.1%   3698.9%    10930    0.02   4126.0%     5.1      1    0.508    1761
     1.73       12566    3224      3521       91.6%     -99.9%    -99.9%    12498    0.00    -99.9%     5.5     -1    0.495    2142
     1.58       13723    3523      3881       90.8%     -99.9%    -99.9%    13623    0.00    -99.9%    -0.2      1    0.499    2425
     1.46       14869    3749      4179       89.7%   41938.5%  42853.7%    14822   -0.02  48580.5%     0.4     -1    0.513    2686
     1.37       12529    3472      4492       77.3%    3293.2%   3343.3%    12379   -0.02   3865.1%    72.8*     2    0.515    2063
     1.29        4927    1748      4733       36.9%    2914.3%   3045.3%     4589    0.04   3519.8%    18.0*     6    0.501     580
    total       88974   23793     29983       79.4%       4.1%      4.3%    87935    6.20      4.8%    99.9*     8    0.620   14060


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     91311
 NUMBER OF REJECTED MISFITS                              2157
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                   44
 NUMBER OF ACCEPTED OBSERVATIONS                        89110
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                  23838



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

   WILSON LINE (using all data) : A=   7.566 B=  35.150 CORRELATION=  0.77
      #      RES      SS        <I>       log(<I>)       BO
     323     7.295   0.005  5.7132E+03       8.651    -115.4
     572     4.516   0.012  4.7279E+03       8.461     -36.5
     763     3.529   0.020  3.6281E+03       8.196     -15.7
     879     2.989   0.028  2.0652E+03       7.633      -1.2
    1020     2.641   0.036  1.3096E+03       7.177       5.4
    1192     2.392   0.044  3.4629E+01       3.545      46.0
    1279     2.200   0.052  2.4218E-01      -1.418      87.0
    1392     2.050   0.059  2.4922E-01      -1.389      75.3
    1474     1.927   0.067 -9.2303E-02    -999.999   -1000.0
    1528     1.822   0.075 -7.8564E-02    -999.999   -1000.0
    1591     1.734   0.083 -1.4534E-01    -999.999   -1000.0
    1688     1.657   0.091 -1.5505E-01    -999.999   -1000.0
    1709     1.589   0.099  3.1997E-01      -1.140      43.9
    1747     1.529   0.107 -1.2272E-01    -999.999   -1000.0
    1822     1.476   0.115  1.6933E-02      -4.079      50.7
    1851     1.428   0.123 -4.2735E-03    -999.999   -1000.0
    1332     1.385   0.130  1.2788E+00       0.246      28.1
    1021     1.345   0.138  6.0601E-01      -0.501      29.2
     609     1.309   0.146  2.3081E+00       0.836      23.1


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

     119     7.295        0.963        0.803        0.652
     141     4.516        0.539        0.292        0.131
     135     3.529        1.386        2.001        2.882
     142     2.989        0.724        0.623        0.500
     152     2.641        0.776        0.738        0.647
     152     2.392        2.366       12.184       48.634
     160     2.050      424.358     1506.465    52321.399
     155     1.822       23.238        0.415      177.454
     162     1.657       26.339       19.087      166.242
     146     1.589      289.598      541.509    30447.802
     152     1.476    86746.943  2736829.146*************
     138     1.428       39.028       59.637      443.417
      46     1.385       13.842       -4.156       39.541
      15     1.309       15.984       12.588       55.088
    2486   overall     5377.245   167404.083135459186.815


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

     204     7.295        0.963        0.975        0.958
     431     4.516        1.136        1.159        1.103
     628     3.529        0.896        0.750        0.580
     737     2.989        1.125        1.381        1.970
     868     2.641        1.056        1.128        1.310
    1040     2.392       35.114     1163.005    36698.060
    1115     2.200      357.326      965.380    85022.382
    1232     2.050      358.668      706.641    89815.598
    1563     1.589      223.188       16.191    34077.337
    1286     1.385      274.592    46841.783  6602950.403
     984     1.345      333.952    42272.153  6554793.974
     594     1.309      101.927     8481.762   644256.137
   21306   overall     6634.425   131799.923353014090.652

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

     119     7.295  1.12  0.95  0.85  0.82  0.84  0.87  0.86  0.87  0.90  0.91
     141     4.516  1.11  1.17  1.06  1.12  1.12  1.09  1.09  1.06  1.07  1.07
     135     3.529  1.01  1.03  0.96  1.03  0.97  0.95  0.94  0.95  0.98  0.99
     142     2.989  1.13  1.08  1.03  1.07  1.06  1.07  1.07  1.06  1.06  1.05
     152     2.641  1.06  1.08  1.11  1.14  1.11  1.11  1.11  1.10  1.10  1.11
     152     2.392  2.97  2.45  2.21  2.03  1.84  1.72  1.63  1.55  1.48  1.43
     160     2.050  2.04  1.46  1.22  1.07  0.97  0.90  0.85  0.81  0.77  0.74
     155     1.822  2.18  1.57  1.30  1.15  1.04  0.98  0.92  0.87  0.83  0.80
     162     1.657  1.94  1.39  1.16  1.02  0.92  0.86  0.81  0.77  0.73  0.71
     146     1.589  1.74  1.27  1.05  0.93  0.86  0.79  0.75  0.71  0.68  0.65
     152     1.476  2.04  1.47  1.22  1.07  0.97  0.90  0.85  0.81  0.77  0.74
     138     1.428  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
      46     1.385  1.58  1.13  0.94  0.83  0.75  0.70  0.66  0.62  0.60  0.57
      15     1.309  2.15  1.54  1.28  1.13  1.02  0.95  1.00  0.95  0.91  0.88
    2486   overall  1.26  1.00  0.88  0.82  0.77  0.73  0.71  0.68  0.67  0.66


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

     204     7.295  1.13  1.03  1.12  1.14  1.10  1.11  1.09  1.05  1.04  1.02
     431     4.516  0.95  0.87  1.01  0.95  0.97  0.96  0.98  0.97  0.97  0.98
     628     3.529  1.10  1.09  1.01  0.98  0.96  0.97  0.94  0.98  0.99  0.99
     737     2.989  1.34  1.11  1.03  1.01  0.99  0.97  0.96  0.96  0.97  0.96
     868     2.641  0.99  0.95  0.95  0.99  0.99  0.95  0.97  0.95  0.96  0.97
    1040     2.392  7.44  4.69  3.56  2.85  2.42  2.12  1.91  1.75  1.62  1.52
    1115     2.200  5.20  2.74  1.93  1.53  1.28  1.12  1.00  0.92  0.85  0.80
    1232     2.050  5.08  2.68  1.88  1.48  1.25  1.09  0.98  0.90  0.84  0.79
    1563     1.589  5.20  2.75  1.93  1.52  1.28  1.12  1.01  0.93  0.86  0.81
    1286     1.385  5.41  2.87  2.04  1.62  1.37  1.21  1.10  1.02  0.96  0.91
     984     1.345  5.22  2.75  1.94  1.53  1.30  1.14  1.03  0.94  0.88  0.83
     594     1.309  5.38  2.88  2.05  1.64  1.40  1.25  1.14  1.05  0.99  0.95
   21306   overall  2.18  1.25  0.93  0.76  0.67  0.60  0.55  0.52  0.49  0.47


 List of    500 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   45   29    5    1.42 6292.82 -0.2689E+02  0.1172E+02 "alien"
   36    8   20    1.42 6068.98 -0.2594E+02  0.9002E+01 "alien"
   54    5    1    1.42 5965.16 -0.2549E+02  0.8726E+01 "alien"
   44   11   15    1.42 5822.19 -0.2488E+02  0.1113E+02 "alien"
   53    9    3    1.43 5757.79 -0.2461E+02  0.1195E+02 "alien"
   36   11   20    1.41 5755.51 -0.2460E+02  0.9287E+01 "alien"
   50   20    3    1.42 5115.49 -0.2186E+02  0.1213E+02 "alien"
   46   14   12    1.44 5092.47 -0.2176E+02  0.1153E+02 "alien"
   53    7    2    1.44 5078.17 -0.2170E+02  0.1114E+02 "alien"
   53    3    6    1.42 4867.45 -0.2080E+02  0.1368E+02 "alien"
   35    5   21    1.41 4807.54 -0.2055E+02  0.8749E+01 "alien"
   51   19    1    1.42 4743.20 -0.2027E+02  0.1149E+02 "alien"
   33    8   21    1.43 4701.47 -0.2009E+02  0.8736E+01 "alien"
   45   13   14    1.41 4665.44 -0.1994E+02  0.9166E+01 "alien"
   46   27    5    1.42 4665.32 -0.1994E+02  0.8948E+01 "alien"
   47   27    4    1.41 4537.63 -0.1939E+02  0.1284E+02 "alien"
   41   35    1    1.43 4442.77 -0.1899E+02  0.8488E+01 "alien"
   53   11    3    1.42 4409.57 -0.1884E+02  0.1193E+02 "alien"
   48   22    5    1.44 4408.36 -0.1884E+02  0.1233E+02 "alien"
   50    4   11    1.41 4397.74 -0.1879E+02  0.7327E+01 "alien"
   51   18    2    1.42 4320.76 -0.1846E+02  0.1165E+02 "alien"
   51    1   10    1.41 4213.29 -0.1801E+02  0.6881E+01 "alien"
   30    9   22    1.43 4126.17 -0.1763E+02  0.8296E+01 "alien"
   48   24    2    1.43 4098.94 -0.1752E+02  0.1199E+02 "alien"
   50   13    8    1.43 4018.82 -0.1717E+02  0.8487E+01 "alien"
   51    2   10    1.41 3996.84 -0.1708E+02  0.7222E+01 "alien"
   50   22    1    1.41 3883.74 -0.1660E+02  0.1263E+02 "alien"
   45   28    7    1.41 3878.69 -0.1658E+02  0.8841E+01 "alien"
   46   27    2    1.44 3814.92 -0.1630E+02  0.1190E+02 "alien"
   42   23   13    1.42 3758.58 -0.1606E+02  0.1158E+02 "alien"
   48   18    8    1.44 3673.31 -0.1570E+02  0.8021E+01 "alien"
   42    5   17    1.42 3671.98 -0.1569E+02  0.1058E+02 "alien"
   34   28   15    1.45 3616.27 -0.1545E+02  0.6187E+01 "alien"
   47    4   14    1.41 3594.78 -0.1536E+02  0.9641E+01 "alien"
   48   26    3    1.41 3589.40 -0.1534E+02  0.8969E+01 "alien"
   35    2   21    1.41 3564.43 -0.1523E+02  0.8608E+01 "alien"
   41   22   13    1.45 3537.88 -0.1512E+02  0.8141E+01 "alien"
   54    4    3    1.42 3536.66 -0.1511E+02  0.8765E+01 "alien"
   47    3   14    1.41 3504.11 -0.1497E+02  0.9773E+01 "alien"
   49   20    4    1.44 3487.74 -0.1491E+02  0.9107E+01 "alien"
   51   16    4    1.43 3463.89 -0.1480E+02  0.1030E+02 "alien"
   46   24    7    1.44 3447.65 -0.1473E+02  0.6142E+01 "alien"
   40   25   13    1.43 3437.08 -0.1469E+02  0.6824E+01 "alien"
   52   16    3    1.41 3430.83 -0.1466E+02  0.1139E+02 "alien"
   51   16    5    1.42 3347.42 -0.1431E+02  0.8631E+01 "alien"
   47    2   13    1.44 3328.19 -0.1422E+02  0.7811E+01 "alien"
   34   32   14    1.42 3318.34 -0.1418E+02  0.7082E+01 "alien"
   53    3    2    1.45 3306.76 -0.1413E+02  0.1175E+02 "alien"
   42   32    6    1.43 3299.93 -0.1410E+02  0.6860E+01 "alien"
   36   30   13    1.44 3263.93 -0.1395E+02  0.5294E+01 "alien"
   30   21   20    1.42 3243.81 -0.1386E+02  0.5316E+01 "alien"
   44   27    7    1.44 3217.04 -0.1375E+02  0.6175E+01 "alien"
   45   16   13    1.42 3199.07 -0.1367E+02  0.1162E+02 "alien"
   54    8    1    1.41 3191.29 -0.1364E+02  0.8760E+01 "alien"
   41    5   17    1.44 3170.32 -0.1355E+02  0.8686E+01 "alien"
   53    7    4    1.43 3163.97 -0.1352E+02  0.1143E+02 "alien"
   39   37    2    1.43 3154.23 -0.1348E+02  0.8460E+01 "alien"
   37   26   15    1.42 3142.63 -0.1343E+02  0.7729E+01 "alien"
   37   36    7    1.44 3139.97 -0.1342E+02  0.7975E+01 "alien"
   40   35    3    1.44 3133.44 -0.1339E+02  0.7167E+01 "alien"
   37    5   19    1.45 3124.23 -0.1335E+02  0.9873E+01 "alien"
   35   18   19    1.41 3123.13 -0.1335E+02  0.6976E+01 "alien"
   42   11   16    1.43 3117.47 -0.1332E+02  0.7730E+01 "alien"
   46   28    2    1.43 3072.71 -0.1313E+02  0.1171E+02 "alien"
   46   28    3    1.42 3060.49 -0.1308E+02  0.8985E+01 "alien"
   52    7    8    1.41 3047.43 -0.1302E+02  0.1247E+02 "alien"
   37    3   20    1.42 3036.61 -0.1298E+02  0.7948E+01 "alien"
   51   14    5    1.43 3030.29 -0.1295E+02  0.8761E+01 "alien"
   44   32    2    1.42 3029.33 -0.1295E+02  0.8702E+01 "alien"
   43   31    4    1.44 3027.57 -0.1294E+02  0.1148E+02 "alien"
   34   13   20    1.43 3025.49 -0.1293E+02  0.6361E+01 "alien"
   39   10   18    1.43 3019.84 -0.1291E+02  0.8380E+01 "alien"
   42    6   17    1.42 3017.85 -0.1290E+02  0.8691E+01 "alien"
   42   17   14    1.45 3010.48 -0.1287E+02  0.1096E+02 "alien"
   44   19   12    1.44 3007.22 -0.1285E+02  0.8244E+01 "alien"
   50   19    1    1.44 3004.45 -0.1284E+02  0.8765E+01 "alien"
   28    7   23    1.42 2977.11 -0.1272E+02  0.7455E+01 "alien"
   48    2   12    1.44 2958.87 -0.1264E+02  0.6788E+01 "alien"
   33    7   21    1.43 2951.08 -0.1261E+02  0.8604E+01 "alien"
   37    1   20    1.42 2947.63 -0.1260E+02  0.8867E+01 "alien"
   47   26    5    1.41 2935.15 -0.1254E+02  0.1265E+02 "alien"
   53   13    1    1.41 2910.26 -0.1244E+02  0.1227E+02 "alien"
   44   30    6    1.41 2906.32 -0.1242E+02  0.1217E+02 "alien"
   41    6   17    1.44 2900.28 -0.1239E+02  0.8902E+01 "alien"
   35    6   20    1.44 2896.74 -0.1238E+02  0.8029E+01 "alien"
   32   23   18    1.45 2881.85 -0.1232E+02  0.5391E+01 "alien"
   53   11    4    1.41 2880.50 -0.1231E+02  0.1116E+02 "alien"
   48   13   10    1.44 2863.03 -0.1224E+02  0.7669E+01 "alien"
   38   31   11    1.44 2852.16 -0.1219E+02  0.5903E+01 "alien"
   49   20    5    1.43 2835.94 -0.1212E+02  0.8951E+01 "alien"
   52    8    6    1.43 2835.91 -0.1212E+02  0.1042E+02 "alien"
   49   22    1    1.44 2832.71 -0.1211E+02  0.1178E+02 "alien"
   37   13   18    1.45 2824.04 -0.1207E+02  0.5838E+01 "alien"
   49    5   12    1.41 2823.20 -0.1207E+02  0.5533E+01 "alien"
   34    2   21    1.43 2821.97 -0.1206E+02  0.8307E+01 "alien"
   51    2    9    1.43 2819.97 -0.1205E+02  0.7396E+01 "alien"
   50   15    8    1.41 2816.71 -0.1204E+02  0.8848E+01 "alien"
   40   34    5    1.44 2798.23 -0.1196E+02  0.8125E+01 "alien"
   44    6   15    1.44 2793.97 -0.1194E+02  0.7698E+01 "alien"
   50   19    2    1.44 2780.94 -0.1188E+02  0.8468E+01 "alien"
   16   12   25    1.43 2768.49 -0.1183E+02  0.5510E+01 "alien"
   49   14    9    1.43 2766.34 -0.1182E+02  0.8018E+01 "alien"
   40   26   12    1.45 2748.28 -0.1174E+02  0.9612E+01 "alien"
   41   28   11    1.42 2747.22 -0.1174E+02  0.7668E+01 "alien"
   38   23   16    1.41 2719.26 -0.1162E+02  0.6392E+01 "alien"
   45    8   14    1.44 2718.89 -0.1162E+02  0.7837E+01 "alien"
   48   19    9    1.41 2714.80 -0.1160E+02  0.6121E+01 "alien"
   42   34    2    1.42 2691.57 -0.1150E+02  0.5900E+01 "alien"
   49   16    7    1.45 2683.55 -0.1147E+02  0.8441E+01 "alien"
   51   11    8    1.41 2640.10 -0.1128E+02  0.6101E+01 "alien"
   43   12   15    1.44 2618.17 -0.1119E+02  0.7476E+01 "alien"
   36   33   11    1.44 2611.09 -0.1116E+02  0.7696E+01 "alien"
   18    9   25    1.43 2603.06 -0.1112E+02  0.5362E+01 "alien"
   45   20   11    1.43 2594.73 -0.1109E+02  0.5877E+01 "alien"
   36   28   14    1.44 2594.13 -0.1109E+02  0.7685E+01 "alien"
   24    6   24    1.43 2591.20 -0.1107E+02  0.8688E+01 "alien"
   39   36    3    1.45 2586.47 -0.1105E+02  0.7363E+01 "alien"
   53   12    2    1.42 2544.14 -0.1087E+02  0.1172E+02 "alien"
   38   29   12    1.44 2519.29 -0.1077E+02  0.8307E+01 "alien"
   39   35    6    1.44 2516.95 -0.1076E+02  0.8085E+01 "alien"
   50   16    8    1.41 2501.07 -0.1069E+02  0.1242E+02 "alien"
   45   22   11    1.41 2488.50 -0.1063E+02  0.6803E+01 "alien"
   41   35    3    1.42 2486.16 -0.1062E+02  0.8343E+01 "alien"
   41   13   16    1.44 2474.67 -0.1058E+02  0.7597E+01 "alien"
   35   33   13    1.41 2466.48 -0.1054E+02  0.7997E+01 "alien"
   47   15   11    1.43 2454.84 -0.1049E+02  0.7751E+01 "alien"
   49   22    3    1.43 2446.42 -0.1045E+02  0.8818E+01 "alien"
   43   32    4    1.42 2445.06 -0.1045E+02  0.1262E+02 "alien"
   39    6   19    1.41 2432.40 -0.1039E+02  0.8743E+01 "alien"
   48    3   13    1.41 2407.48 -0.1029E+02  0.7322E+01 "alien"
   49   14    8    1.45 2404.38 -0.1028E+02  0.6761E+01 "alien"
   39    4   19    1.41 2393.01 -0.1023E+02  0.8168E+01 "alien"
   37   21   16    1.45 2382.02 -0.1018E+02  0.7681E+01 "alien"
   53    6    2    1.44 2380.19 -0.1017E+02  0.1137E+02 "alien"
   39   23   14    1.45 2372.18 -0.1014E+02  0.9240E+01 "alien"
   43   26   11    1.41 2370.89 -0.1013E+02  0.6112E+01 "alien"
   39   14   17    1.44 2364.46 -0.1010E+02  0.7380E+01 "alien"
   53    7    3    1.44 2353.80 -0.1006E+02  0.1201E+02 "alien"
   33   18   19    1.44 2338.80 -0.9995E+01  0.7489E+01 "alien"
   50   18    3    1.44 2337.72 -0.9990E+01  0.8934E+01 "alien"
   35   20   18    1.43 2334.82 -0.9978E+01  0.9829E+01 "alien"
   34   30   14    1.45 2328.29 -0.9950E+01  0.6097E+01 "alien"
   25    6   24    1.42 2321.35 -0.9920E+01  0.9104E+01 "alien"
   41   25   12    1.44 2298.87 -0.9824E+01  0.8284E+01 "alien"
   46   20   11    1.41 2297.35 -0.9818E+01  0.7191E+01 "alien"
   41   19   15    1.42 2283.85 -0.9760E+01  0.1123E+02 "alien"
   42    8   16    1.45 2280.21 -0.9745E+01  0.7176E+01 "alien"
   31   12   21    1.44 2279.86 -0.9743E+01  0.5894E+01 "alien"
   27   26   19    1.45 2278.13 -0.9736E+01  0.6000E+01 "alien"
   36   24   16    1.43 2268.81 -0.9696E+01  0.5986E+01 "alien"
   51   13    7    1.42 2253.51 -0.9630E+01  0.1191E+02 "alien"
   36   35    9    1.45 2251.25 -0.9621E+01  0.6882E+01 "alien"
   31   16   21    1.41 2249.08 -0.9612E+01  0.5019E+01 "alien"
   50    7    9    1.44 2234.90 -0.9551E+01  0.4813E+01 "alien"
   41   29   10    1.43 2234.64 -0.9550E+01  0.8086E+01 "alien"
   38   23   15    1.44 2219.68 -0.9486E+01  0.5589E+01 "alien"
   42   34    1    1.43 2196.97 -0.9389E+01  0.6013E+01 "alien"
   47   16   10    1.44 2196.65 -0.9387E+01  0.8331E+01 "alien"
   41    4   18    1.41 2192.43 -0.9369E+01  0.8728E+01 "alien"
   40   34    7    1.42 2185.13 -0.9338E+01  0.7317E+01 "alien"
   53    3    5    1.43 2175.16 -0.9296E+01  0.8835E+01 "alien"
   53   12    3    1.41 2174.43 -0.9293E+01  0.1207E+02 "alien"
   51   12    5    1.45 2160.96 -0.9235E+01  0.6040E+01 "alien"
   33   31   15    1.42 2160.55 -0.9233E+01  0.7223E+01 "alien"
   52    5    7    1.43 2160.05 -0.9231E+01  0.5745E+01 "alien"
   47   24    7    1.41 2156.98 -0.9218E+01  0.8735E+01 "alien"
   43   13   15    1.43 2156.45 -0.9216E+01  0.1104E+02 "alien"
   40   11   17    1.44 2155.08 -0.9210E+01  0.7051E+01 "alien"
   32   31   16    1.41 2153.73 -0.9204E+01  0.7707E+01 "alien"
   47    9   11    1.47 2140.21  0.3624E+02  0.1126E+02 "alien"
   41   32    7    1.43 2131.55 -0.9109E+01  0.6792E+01 "alien"
   47    9   12    1.44 2128.93 -0.9098E+01  0.1143E+02 "alien"
   52    4    8    1.42 2118.74 -0.9055E+01  0.6340E+01 "alien"
   49   11   11    1.41 2117.01 -0.9047E+01  0.8244E+01 "alien"
   37   25   16    1.41 2105.99 -0.9000E+01  0.6443E+01 "alien"
   51    1    8    1.44 2105.70 -0.8999E+01  0.7207E+01 "alien"
   51    5    8    1.44 2104.20 -0.8992E+01  0.6037E+01 "alien"
   45    9   15    1.41 2102.95 -0.8987E+01  0.7519E+01 "alien"
   45   19   11    1.44 2101.78 -0.8982E+01  0.5842E+01 "alien"
   23    9   24    1.43 2100.49 -0.8977E+01  0.8812E+01 "alien"
   36   28   15    1.41 2094.77 -0.8952E+01  0.6469E+01 "alien"
   39    3   19    1.42 2093.98 -0.8949E+01  0.1024E+02 "alien"
   33   15   20    1.43 2077.73 -0.8879E+01  0.6578E+01 "alien"
   41   16   15    1.45 2076.96 -0.8876E+01  0.6271E+01 "alien"
   36   31   12    1.45 2063.70 -0.8819E+01  0.7001E+01 "alien"
   51   17    1    1.43 2061.12 -0.8808E+01  0.1241E+02 "alien"
   37   34   10    1.43 2057.90 -0.8795E+01  0.6294E+01 "alien"
   35   31   14    1.42 2057.35 -0.8792E+01  0.6183E+01 "alien"
   47   13   12    1.42 2045.93 -0.8743E+01  0.5878E+01 "alien"
   50   21    3    1.41 2044.86 -0.8739E+01  0.1297E+02 "alien"
   51   17    2    1.43 2036.56 -0.8703E+01  0.1224E+02 "alien"
   50    3   11    1.41 2028.85 -0.8670E+01  0.7342E+01 "alien"
   44   26   10    1.41 2028.38 -0.8668E+01  0.8177E+01 "alien"
   27   26   20    1.41 2025.78 -0.8657E+01  0.6804E+01 "alien"
   46   12   12    1.45 2025.20 -0.8655E+01  0.8058E+01 "alien"
   41   29   11    1.41 2018.21 -0.8625E+01  0.5597E+01 "alien"
   22    8   24    1.45 2015.45 -0.8613E+01  0.3990E+01 "alien"
   39   36    5    1.43 2012.60 -0.8601E+01  0.7611E+01 "alien"
   37    2   20    1.42 2010.92 -0.8594E+01  0.9259E+01 "alien"
   40   30   11    1.41 2001.92 -0.8555E+01  0.6014E+01 "alien"
   24   19   22    1.44 2001.10 -0.8552E+01  0.5749E+01 "alien"
   50    9    8    1.45 1985.94 -0.8487E+01  0.6227E+01 "alien"
   39    9   18    1.43 1980.33 -0.8463E+01  0.1027E+02 "alien"
   39    5   18    1.45 1971.77 -0.8426E+01  0.8101E+01 "alien"
   34   23   18    1.41 1971.75 -0.8426E+01  0.6097E+01 "alien"
   45    5   14    1.45 1962.68 -0.8388E+01  0.7866E+01 "alien"
   52    9    5    1.44 1959.22 -0.8373E+01  0.1174E+02 "alien"
   39    1   19    1.42 1945.43 -0.8314E+01  0.7011E+01 "alien"
   42   33    6    1.41 1945.11 -0.8313E+01  0.7222E+01 "alien"
   52   17    1    1.41 1944.09 -0.8308E+01  0.8910E+01 "alien"
   47   28    2    1.41 1943.33 -0.8305E+01  0.1241E+02 "alien"
   54    5    2    1.42 1942.93 -0.8303E+01  0.8729E+01 "alien"
   37   35   10    1.41 1940.51 -0.8293E+01  0.7569E+01 "alien"
   28   10   23    1.41 1936.08 -0.8274E+01  0.5386E+01 "alien"
   28   13   22    1.44 1932.56 -0.8259E+01  0.5753E+01 "alien"
   36   23   16    1.45 1931.25 -0.8253E+01  0.7275E+01 "alien"
   52   15    4    1.41 1929.71 -0.8247E+01  0.1233E+02 "alien"
   53    4    6    1.42 1902.33 -0.8130E+01  0.1176E+02 "alien"
   44   20   13    1.41 1900.56 -0.8122E+01  0.6698E+01 "alien"
   37   23   16    1.43 1894.96 -0.8098E+01  0.7554E+01 "alien"
   31   24   18    1.45 1891.00 -0.8081E+01  0.5588E+01 "alien"
   45   22    9    1.45 1883.70 -0.8050E+01  0.5974E+01 "alien"
   34   24   17    1.44 1880.07 -0.8035E+01  0.5375E+01 "alien"
   34   31   15    1.41 1873.92 -0.8008E+01  0.5281E+01 "alien"
   53    6    5    1.42 1872.94 -0.8004E+01  0.8727E+01 "alien"
   51   15    5    1.43 1872.47 -0.8002E+01  0.8578E+01 "alien"
   38   11   18    1.44 1868.37 -0.7985E+01  0.6754E+01 "alien"
   42   14   16    1.41 1854.98 -0.7927E+01  0.6418E+01 "alien"
   46    5   14    1.43 1849.10 -0.7902E+01  0.6355E+01 "alien"
   48   10   12    1.41 1848.02 -0.7898E+01  0.5835E+01 "alien"
   46   26    4    1.44 1847.97 -0.7897E+01  0.1279E+02 "alien"
   32   20   19    1.44 1846.02 -0.7889E+01  0.5621E+01 "alien"
   53   13    2    1.41 1842.21 -0.7873E+01  0.1205E+02 "alien"
   40   36    3    1.43 1838.36 -0.7856E+01  0.6163E+01 "alien"
   42    9   16    1.44 1838.00 -0.7855E+01  0.5511E+01 "alien"
   43   24   12    1.41 1832.79 -0.7832E+01  0.8415E+01 "alien"
   53   10    5    1.41 1831.95 -0.7829E+01  0.8444E+01 "alien"
   46   22    9    1.43 1830.83 -0.7824E+01  0.7937E+01 "alien"
   42   28    9    1.44 1826.54 -0.7806E+01  0.8552E+01 "alien"
   38   29   13    1.42 1825.81 -0.7803E+01  0.9903E+01 "alien"
   39   16   17    1.43 1825.70 -0.7802E+01  0.6173E+01 "alien"
   52   15    1    1.43 1825.23 -0.7800E+01  0.8522E+01 "alien"
   52    2    6    1.45 1817.40 -0.7767E+01  0.7624E+01 "alien"
   48   16   10    1.42 1817.04 -0.7765E+01  0.7834E+01 "alien"
   35    5   20    1.45 1815.68 -0.7759E+01  0.9052E+01 "alien"
   31   30   17    1.41 1802.61 -0.7704E+01  0.6127E+01 "alien"
   38   18   17    1.43 1801.73 -0.7700E+01  0.6264E+01 "alien"
   48   23    6    1.41 1800.97 -0.7697E+01  0.8743E+01 "alien"
   47   26    2    1.43 1800.70 -0.7695E+01  0.8747E+01 "alien"
   43   32    1    1.44 1800.61 -0.7695E+01  0.8720E+01 "alien"
   51   18    3    1.42 1782.85 -0.7619E+01  0.8846E+01 "alien"
   49   16    8    1.43 1781.38 -0.7613E+01  0.6022E+01 "alien"
   50   10    8    1.44 1775.92 -0.7589E+01  0.8261E+01 "alien"
   41   24   13    1.43 1774.46 -0.7583E+01  0.1112E+02 "alien"
   46   15   13    1.41 1773.79 -0.7580E+01  0.6692E+01 "alien"
   41    2   18    1.41 1769.72 -0.7563E+01  0.6986E+01 "alien"
   36   32   13    1.41 1768.89 -0.7559E+01  0.6644E+01 "alien"
   46   13   13    1.42 1765.30 -0.7544E+01  0.1109E+02 "alien"
   49    4   12    1.41 1745.88 -0.7461E+01  0.6924E+01 "alien"
   41   33    5    1.44 1739.62 -0.7434E+01  0.7141E+01 "alien"
   25   20   22    1.42 1735.44 -0.7416E+01  0.4287E+01 "alien"
   49   19    6    1.43 1734.97 -0.7414E+01  0.1096E+02 "alien"
   41   11   17    1.42 1734.82 -0.7414E+01  0.7140E+01 "alien"
   54    7    1    1.42 1732.47 -0.7404E+01  0.1079E+02 "alien"
   43   12   16    1.41 1730.40 -0.7395E+01  0.6552E+01 "alien"
   34    1   21    1.43 1722.89 -0.7363E+01  0.1147E+02 "alien"
   44    6   16    1.41 1722.03 -0.7359E+01  0.1077E+02 "alien"
   29   26   19    1.42 1712.33 -0.7318E+01  0.5054E+01 "alien"
   47   21    7    1.45 1711.55 -0.7314E+01  0.8511E+01 "alien"
   32    7   21    1.45 1707.41 -0.7297E+01  0.8232E+01 "alien"
   49   19    7    1.42 1705.23 -0.7287E+01  0.1191E+02 "alien"
   40    7   18    1.42 1694.47 -0.7241E+01  0.1042E+02 "alien"
   42   30    8    1.43 1692.29 -0.7232E+01  0.8839E+01 "alien"
   27   20   21    1.44 1690.93 -0.7226E+01  0.4892E+01 "alien"
   39   22   15    1.43 1682.60 -0.7191E+01  0.6493E+01 "alien"
   29   17   21    1.44 1675.82 -0.7162E+01  0.4848E+01 "alien"
   38   25   15    1.42 1667.38 -0.7126E+01  0.5676E+01 "alien"
   37   30   13    1.42 1665.93 -0.7119E+01  0.8365E+01 "alien"
   50    9   10    1.41 1661.06 -0.7099E+01  0.8396E+01 "alien"
   39   38    3    1.41 1658.58 -0.7088E+01  0.7464E+01 "alien"
   46   10   14    1.41 1658.49 -0.7088E+01  0.9425E+01 "alien"
   19    3   25    1.44 1644.91 -0.7030E+01  0.4336E+01 "alien"
   53    8    1    1.44 1644.27 -0.7027E+01  0.1081E+02 "alien"
   54    8    2    1.41 1643.35 -0.7023E+01  0.8405E+01 "alien"
   53    5    3    1.44 1640.09 -0.7009E+01  0.1227E+02 "alien"
   49   22    4    1.42 1639.77 -0.7008E+01  0.9055E+01 "alien"
   22    1   25    1.41 1639.29 -0.7006E+01  0.6976E+01 "alien"
   50    5   10    1.43 1625.22 -0.6945E+01  0.5429E+01 "alien"
   41   23   13    1.44 1620.95 -0.6927E+01  0.8132E+01 "alien"
   49   24    2    1.41 1617.88 -0.6914E+01  0.8725E+01 "alien"
   35   32   13    1.43 1612.48 -0.6891E+01  0.8382E+01 "alien"
   45   19   12    1.42 1605.39 -0.6861E+01  0.7803E+01 "alien"
   23   10   24    1.43 1601.29 -0.6843E+01  0.5169E+01 "alien"
   43   19   13    1.44 1597.47 -0.6827E+01  0.5758E+01 "alien"
   44   31    4    1.42 1595.24 -0.6817E+01  0.9067E+01 "alien"
   44   29    6    1.43 1590.95 -0.6799E+01  0.8668E+01 "alien"
   39   35    8    1.41 1590.88 -0.6799E+01  0.7388E+01 "alien"
   43   27    9    1.43 1589.75 -0.6794E+01  0.6196E+01 "alien"
   51   17    5    1.41 1583.26 -0.6766E+01  0.8930E+01 "alien"
   43    1   17    1.41 1583.11 -0.6765E+01  0.6466E+01 "alien"
   28    6   23    1.43 1580.90 -0.6756E+01  0.7429E+01 "alien"
   47   23    8    1.41 1575.18 -0.6732E+01  0.8965E+01 "alien"
   45   30    2    1.42 1574.47 -0.6729E+01  0.1215E+02 "alien"
   30   25   19    1.42 1574.02 -0.6727E+01  0.5199E+01 "alien"
   38   36    6    1.44 1569.52 -0.6707E+01  0.5969E+01 "alien"
   52   11    2    1.45 1555.75 -0.6649E+01  0.8233E+01 "alien"
   38   37    5    1.43 1554.89 -0.6645E+01  0.7517E+01 "alien"
   54    9    2    1.41 1553.79 -0.6640E+01  0.8353E+01 "alien"
   26   15   22    1.45 1552.52 -0.6635E+01  0.4559E+01 "alien"
   29   18   21    1.43 1551.00 -0.6628E+01  0.4970E+01 "alien"
   39   20   16    1.42 1535.98 -0.6564E+01  0.6365E+01 "alien"
   47   27    3    1.42 1535.93 -0.6564E+01  0.1277E+02 "alien"
   48   21    5    1.45 1533.92 -0.6555E+01  0.8687E+01 "alien"
   48   15   10    1.43 1530.69 -0.6541E+01  0.7732E+01 "alien"
   42   10   17    1.41 1528.87 -0.6534E+01  0.5499E+01 "alien"
   31    6   22    1.43 1521.34 -0.6502E+01  0.8095E+01 "alien"
   37   31   12    1.43 1516.92 -0.6483E+01  0.6817E+01 "alien"
   40   21   14    1.45 1508.57  0.2554E+02  0.1132E+02 "alien"
   37   25   15    1.44 1505.06 -0.6432E+01  0.5645E+01 "alien"
   54    4    1    1.42 1504.75 -0.6431E+01  0.8624E+01 "alien"
   49    6   11    1.43 1487.80 -0.6358E+01  0.4856E+01 "alien"
   34   27   16    1.43 1481.58 -0.6332E+01  0.6153E+01 "alien"
   29   12   22    1.43 1475.52 -0.6306E+01  0.4471E+01 "alien"
   41   33    4    1.45 1473.61 -0.6298E+01  0.7073E+01 "alien"
   41    5   18    1.41 1470.30 -0.6283E+01  0.1044E+02 "alien"
   23   21   22    1.43 1458.41 -0.6233E+01  0.4086E+01 "alien"
   23   11   24    1.42 1457.14 -0.6227E+01  0.5436E+01 "alien"
   36   20   18    1.41 1454.14 -0.6214E+01  0.7543E+01 "alien"
   46   24    6    1.45 1452.47 -0.6207E+01  0.8351E+01 "alien"
   43   25   10    1.44 1451.24 -0.6202E+01  0.7786E+01 "alien"
   39   26   13    1.44 1443.21 -0.6168E+01  0.6654E+01 "alien"
   49   21    2    1.44 1432.86 -0.6123E+01  0.1240E+02 "alien"
   39   24   14    1.44 1432.32 -0.6121E+01  0.6156E+01 "alien"
   30   13   22    1.41 1431.08 -0.6116E+01  0.8609E+01 "alien"
   23   20   22    1.44 1422.41 -0.6079E+01  0.5559E+01 "alien"
   50    5    9    1.45 1422.16 -0.6078E+01  0.6425E+01 "alien"
   22    5   25    1.41 1421.64 -0.6075E+01  0.7356E+01 "alien"
   45   23    6    1.49 1418.16  0.2401E+02  0.8496E+01 "alien"
   29   15   22    1.41 1416.13 -0.6052E+01  0.8617E+01 "alien"
   53    5    5    1.43 1409.88 -0.6025E+01  0.8584E+01 "alien"
   35   34   12    1.42 1408.01 -0.6017E+01  0.7850E+01 "alien"
   35   30   14    1.43 1401.70 -0.5990E+01  0.4909E+01 "alien"
   42   26   12    1.41 1399.78 -0.5982E+01  0.6089E+01 "alien"
   46   29    2    1.42 1386.71 -0.5926E+01  0.8645E+01 "alien"
   44   19   11    1.46 1383.88  0.2343E+02  0.9727E+01 "alien"
   52   15    2    1.42 1382.19 -0.5907E+01  0.8344E+01 "alien"
   52   13    2    1.44 1379.38 -0.5895E+01  0.8379E+01 "alien"
   52    7    5    1.45 1377.93 -0.5889E+01  0.6483E+01 "alien"
   36    5   20    1.43 1367.99 -0.5846E+01  0.7925E+01 "alien"
   42   20   13    1.45 1364.20 -0.5830E+01  0.7696E+01 "alien"
   27   11   23    1.42 1358.95 -0.5808E+01  0.5272E+01 "alien"
   50   13    9    1.41 1355.67 -0.5794E+01  0.6093E+01 "alien"
   50   15    6    1.44 1344.39 -0.5745E+01  0.8458E+01 "alien"
   47   20    9    1.43 1336.15 -0.5710E+01  0.7012E+01 "alien"
   29   20   21    1.41 1335.00 -0.5705E+01  0.5350E+01 "alien"
   15   12   25    1.44 1331.74 -0.5691E+01  0.4429E+01 "alien"
   24    8   24    1.42 1330.40 -0.5686E+01  0.5785E+01 "alien"
   47    6   13    1.43 1327.26 -0.5672E+01  0.6181E+01 "alien"
   41   31    7    1.45 1321.82 -0.5649E+01  0.8181E+01 "alien"
   37   11   19    1.43 1317.46 -0.5630E+01  0.7004E+01 "alien"
   51   14    3    1.45 1307.35 -0.5587E+01  0.1204E+02 "alien"
   30   19   20    1.44 1303.68 -0.5571E+01  0.5324E+01 "alien"
   33   12   20    1.45 1301.40 -0.5562E+01  0.5938E+01 "alien"
   21   18   23    1.44 1299.34 -0.5553E+01  0.3533E+01 "alien"
   38    7   19    1.43 1298.06 -0.5547E+01  0.1002E+02 "alien"
   52   14    1    1.43 1296.32 -0.5540E+01  0.8596E+01 "alien"
   24    7   24    1.43 1295.12 -0.5535E+01  0.8161E+01 "alien"
   37   14   19    1.41 1290.76 -0.5516E+01  0.6788E+01 "alien"
   30   22   20    1.41 1289.73 -0.5512E+01  0.5531E+01 "alien"
   48   26    2    1.41 1286.40 -0.5497E+01  0.8729E+01 "alien"
   46    2   15    1.41 1285.42 -0.5493E+01  0.6052E+01 "alien"
   44   27    9    1.41 1285.11 -0.5492E+01  0.6200E+01 "alien"
   31   11   22    1.41 1283.91 -0.5487E+01  0.8235E+01 "alien"
   25   19   22    1.43 1281.40 -0.5476E+01  0.4004E+01 "alien"
   36    3   20    1.43 1277.36 -0.5459E+01  0.8417E+01 "alien"
   29   27   19    1.41 1276.55 -0.5455E+01  0.4441E+01 "alien"
   48   11   10    1.45 1273.93  0.2157E+02  0.8272E+01 "alien"
   53    4    2    1.45 1272.97 -0.5440E+01  0.1184E+02 "alien"
   49   23    1    1.43 1271.86 -0.5435E+01  0.8782E+01 "alien"
   52    8    5    1.44 1270.42 -0.5429E+01  0.6330E+01 "alien"
   35   31   11    1.49 1268.73  0.2148E+02  0.9072E+01 "alien"
   18   15   24    1.45 1262.17 -0.5394E+01  0.2829E+01 "alien"
   41   35    5    1.41 1260.72 -0.5388E+01  0.6087E+01 "alien"
   32   30   16    1.42 1246.01 -0.5325E+01  0.5609E+01 "alien"
   32   11   21    1.43 1241.69 -0.5306E+01  0.5988E+01 "alien"
   34   22   16    1.50 1240.62  0.2101E+02  0.1355E+02 "alien"
   42   30    7    1.44 1237.87 -0.5290E+01  0.6159E+01 "alien"
   42   29    8    1.44 1236.20 -0.5283E+01  0.6175E+01 "alien"
   49   16    3    1.49 1235.41  0.2092E+02  0.6312E+01 "alien"
   46   18   12    1.41 1233.98 -0.5273E+01  0.8352E+01 "alien"
   34   28   16    1.42 1233.15 -0.5270E+01  0.5525E+01 "alien"
   41   31    9    1.42 1231.50 -0.5263E+01  0.6803E+01 "alien"
   33   22   18    1.44 1229.71 -0.5255E+01  0.5433E+01 "alien"
   35    4   20    1.45 1224.60 -0.5233E+01  0.9156E+01 "alien"
   28   17   21    1.45 1224.00 -0.5231E+01  0.8207E+01 "alien"
   41   17   15    1.44 1221.06 -0.5218E+01  0.7531E+01 "alien"
   32   27   16    1.47 1215.72  0.2059E+02  0.5543E+01 "alien"
   21   20   23    1.42 1213.39 -0.5185E+01  0.3930E+01 "alien"
   37   30   12    1.45 1213.31 -0.5185E+01  0.7067E+01 "alien"
   36   35   10    1.43 1212.83 -0.5183E+01  0.7581E+01 "alien"
   34    5   21    1.42 1209.31 -0.5168E+01  0.8512E+01 "alien"
   28   14   22    1.43 1196.68 -0.5114E+01  0.4447E+01 "alien"
   47   12   12    1.43 1193.10 -0.5099E+01  0.6693E+01 "alien"
   44   29    4    1.45 1192.19 -0.5095E+01  0.9262E+01 "alien"
   44    1   16    1.42 1191.08 -0.5090E+01  0.6666E+01 "alien"
   44   31    2    1.43 1189.68 -0.5084E+01  0.8648E+01 "alien"
   38    5   19    1.43 1187.70 -0.5076E+01  0.1018E+02 "alien"
   25   21   22    1.41 1186.68 -0.5071E+01  0.5755E+01 "alien"
   29   24   20    1.41 1179.26 -0.5040E+01  0.6248E+01 "alien"
   50    9    9    1.43 1174.40 -0.5019E+01  0.5452E+01 "alien"
   47   18   11    1.41 1172.59 -0.5011E+01  0.8158E+01 "alien"
   37   28   14    1.42 1167.50 -0.4989E+01  0.6399E+01 "alien"
   36   13   19    1.43 1167.43 -0.4989E+01  0.5752E+01 "alien"
   47    2   14    1.41 1164.49 -0.4976E+01  0.1066E+02 "alien"
   51    2    7    1.46 1160.93  0.1966E+02  0.7606E+01 "alien"
   47   25    3    1.44 1157.22 -0.4945E+01  0.8809E+01 "alien"
   25   23   21    1.43 1156.88 -0.4944E+01  0.4415E+01 "alien"
   45   29    1    1.44 1150.11 -0.4915E+01  0.8790E+01 "alien"
   42   22   12    1.45 1150.05  0.1947E+02  0.8330E+01 "alien"
   51   17    3    1.43 1141.39 -0.4878E+01  0.1190E+02 "alien"
   43   28    5    1.48 1135.72  0.1923E+02  0.8735E+01 "alien"
   29    8   22    1.45 1133.70 -0.4845E+01  0.7763E+01 "alien"
   45   17   12    1.44 1128.85 -0.4824E+01  0.1105E+02 "alien"
   46   19   10    1.44 1128.67 -0.4823E+01  0.5995E+01 "alien"
   32    5   22    1.41 1127.69 -0.4819E+01  0.8647E+01 "alien"
   42   31    8    1.41 1127.06 -0.4817E+01  0.6585E+01 "alien"
   41   30   10    1.41 1126.01 -0.4812E+01  0.6387E+01 "alien"
   49   10   11    1.41 1125.74 -0.4811E+01  0.7879E+01 "alien"
   49   24    1    1.41 1125.45 -0.4810E+01  0.8956E+01 "alien"
   29    5   23    1.41 1124.76 -0.4807E+01  0.7293E+01 "alien"
   21    9   25    1.41 1122.36 -0.4796E+01  0.5255E+01 "alien"
   40   30   10    1.43 1118.63 -0.4781E+01  0.6709E+01 "alien"
   43   21   13    1.42 1116.50 -0.4771E+01  0.8225E+01 "alien"
   35   16   19    1.43 1115.33 -0.4766E+01  0.6921E+01 "alien"
   34   33   11    1.48 1111.50  0.1882E+02  0.7497E+01 "alien"
   24   12   23    1.45 1101.78 -0.4708E+01  0.3806E+01 "alien"
   45   13   13    1.44 1098.54 -0.4695E+01  0.7566E+01 "alien"
   18   10   25    1.43 1096.33 -0.4685E+01  0.4147E+01 "alien"
   45   24   10    1.41 1095.71 -0.4683E+01  0.6949E+01 "alien"
   54    6    4    1.41 1095.30 -0.4681E+01  0.9065E+01 "alien"
   36   33   12    1.42 1094.95 -0.4679E+01  0.5867E+01 "alien"
   44   18   12    1.45 1092.22 -0.4668E+01  0.5864E+01 "alien"
   40   33    6    1.45 1090.81 -0.4662E+01  0.7962E+01 "alien"
   45   28    2    1.45 1088.06  0.1842E+02  0.1029E+02 "alien"
   54    4    2    1.42 1081.14 -0.4620E+01  0.8440E+01 "alien"
   42   20   14    1.42 1069.55 -0.4571E+01  0.6798E+01 "alien"
   38   20   17    1.41 1069.48 -0.4570E+01  0.7701E+01 "alien"
   39   18   16    1.44 1058.66 -0.4524E+01  0.7375E+01 "alien"
   27    7   23    1.44 1057.91 -0.4521E+01  0.7376E+01 "alien"
   48    9   11    1.44 1049.85 -0.4487E+01  0.6558E+01 "alien"
   44   11   13    1.48 1049.55  0.1777E+02  0.1091E+02 "alien"
   46   22    8    1.44 1047.84 -0.4478E+01  0.6100E+01 "alien"
   29   19   21    1.42 1045.99 -0.4470E+01  0.5033E+01 "alien"
   43   30    2    1.47 1044.50  0.1769E+02  0.8674E+01 "alien"
   35   29   15    1.42 1038.43 -0.4438E+01  0.6116E+01 "alien"
   44   23   10    1.44 1037.10 -0.4432E+01  0.6768E+01 "alien"
   44   10   15    1.43 1035.54 -0.4425E+01  0.1114E+02 "alien"
   35   32   12    1.45 1033.33  0.1750E+02  0.7930E+01 "alien"
   44   22   11    1.43 1028.49 -0.4395E+01  0.8312E+01 "alien"
   44   30    5    1.42 1025.20 -0.4381E+01  0.8814E+01 "alien"
   19    4   25    1.44 1021.85 -0.4367E+01  0.4501E+01 "alien"
   46   22    4    1.50 1018.90  0.1725E+02  0.8030E+01 "alien"
   26   16   22    1.44 1018.81 -0.4354E+01  0.5549E+01 "alien"
   28    9   23    1.41 1018.26 -0.4352E+01  0.5080E+01 "alien"
   23   16   23    1.43 1017.46 -0.4348E+01  0.3763E+01 "alien"
   42   19   14    1.43 1017.30 -0.4347E+01  0.6756E+01 "alien"
   17   11   25    1.43 1013.36 -0.4331E+01  0.3954E+01 "alien"
   43   33    3    1.42 1012.55 -0.4327E+01  0.6606E+01 "alien"
   40    9   17    1.45 1007.03 -0.4304E+01  0.7185E+01 "alien"
   47    3   13    1.44 1001.16 -0.4279E+01  0.7647E+01 "alien"
   48   12   10    1.45 1000.33 -0.4275E+01  0.7663E+01 "alien"
   21   19   23    1.43  995.99 -0.4256E+01  0.3658E+01 "alien"
   43   27    8    1.45  993.67 -0.4246E+01  0.8746E+01 "alien"
   49   10   10    1.43  992.59 -0.4242E+01  0.5954E+01 "alien"
   32   10   21    1.44  991.28 -0.4236E+01  0.5722E+01 "alien"
   38   20   16    1.44  991.06 -0.4235E+01  0.7394E+01 "alien"
   45   26    1    1.48  983.37  0.1665E+02  0.6381E+01 "alien"
   44   15   13    1.45  981.00 -0.4192E+01  0.5693E+01 "alien"
   39   37    4    1.42  978.61 -0.4182E+01  0.6316E+01 "alien"
   31    5   22    1.43  977.59 -0.4178E+01  0.7874E+01 "alien"
   21    7   25    1.41  975.08 -0.4167E+01  0.5177E+01 "alien"
   45   29    4    1.43  974.25 -0.4164E+01  0.1320E+02 "alien"
   29   22   20    1.43  971.80 -0.4153E+01  0.5395E+01 "alien"
   40   19   15    1.44  967.30 -0.4134E+01  0.6536E+01 "alien"
   49    2   11    1.44  965.14 -0.4125E+01  0.7448E+01 "alien"
   40    1   18    1.43  964.13 -0.4120E+01  0.5977E+01 "alien"
   38   26   15    1.41  963.17 -0.4116E+01  0.5217E+01 "alien"
   39   25   12    1.48  962.79  0.1630E+02  0.6631E+01 "alien"
   48    6   13    1.41  960.75 -0.4106E+01  0.5855E+01 "alien"
   40   27   11    1.46  959.34  0.1624E+02  0.1181E+02 "alien"
   46   18   10    1.45  959.28 -0.4100E+01  0.5904E+01 "alien"
   51    9    7    1.44  956.07 -0.4086E+01  0.6881E+01 "alien"
   48   24    1    1.44  955.56 -0.4084E+01  0.1134E+02 "alien"
   43   27    7    1.47  953.55  0.1615E+02  0.6286E+01 "alien"
   29    6   23    1.41  953.23 -0.4074E+01  0.7378E+01 "alien"
   23    8   24    1.44  952.86 -0.4072E+01  0.5945E+01 "alien"
   35   12   18    1.50  951.94  0.1612E+02  0.6542E+01 "alien"
   40   27   12    1.43  948.44 -0.4053E+01  0.8434E+01 "alien"
   49   21    5    1.42  946.59 -0.4045E+01  0.8791E+01 "alien"
   37   20   15    1.49  945.41  0.1601E+02  0.7288E+01 "alien"


 List of    500 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   46    2   15    1.41 1285.42 -0.5493E+01  0.6052E+01
   38   26   15    1.41  963.17 -0.4116E+01  0.5217E+01
   42   10   17    1.41 1528.87 -0.6534E+01  0.5499E+01
   21    9   25    1.41 1122.36 -0.4796E+01  0.5255E+01
   47   18   11    1.41 1172.59 -0.5011E+01  0.8158E+01
   37   25   16    1.41 2105.99 -0.9000E+01  0.6443E+01
   48   26    3    1.41 3589.40 -0.1534E+02  0.8969E+01
   54    6    4    1.41 1095.30 -0.4681E+01  0.9065E+01
   46   15   13    1.41 1773.79 -0.7580E+01  0.6692E+01
   34   31   15    1.41 1873.92 -0.8008E+01  0.5281E+01
   43   12   16    1.41 1730.40 -0.7395E+01  0.6552E+01
   42   26   12    1.41 1399.78 -0.5982E+01  0.6089E+01
   46   18   12    1.41 1233.98 -0.5273E+01  0.8352E+01
   31   30   17    1.41 1802.61 -0.7704E+01  0.6127E+01
   48    6   13    1.41  960.75 -0.4106E+01  0.5855E+01
   31   11   22    1.41 1283.91 -0.5487E+01  0.8235E+01
   36   11   20    1.41 5755.51 -0.2460E+02  0.9287E+01
   29   24   20    1.41 1179.26 -0.5040E+01  0.6248E+01
   54    9    2    1.41 1553.79 -0.6640E+01  0.8353E+01
   50   16    8    1.41 2501.07 -0.1069E+02  0.1242E+02
   29   27   19    1.41 1276.55 -0.5455E+01  0.4441E+01
   44   26   10    1.41 2028.38 -0.8668E+01  0.8177E+01
   44   20   13    1.41 1900.56 -0.8122E+01  0.6698E+01
   35    5   21    1.41 4807.54 -0.2055E+02  0.8749E+01
   47   28    2    1.41 1943.33 -0.8305E+01  0.1241E+02
   53   10    5    1.41 1831.95 -0.7829E+01  0.8444E+01
   52    7    8    1.41 3047.43 -0.1302E+02  0.1247E+02
   43   26   11    1.41 2370.89 -0.1013E+02  0.6112E+01
   47    4   14    1.41 3594.78 -0.1536E+02  0.9641E+01
   32   31   16    1.41 2153.73 -0.9204E+01  0.7707E+01
   37   14   19    1.41 1290.76 -0.5516E+01  0.6788E+01
   41   35    5    1.41 1260.72 -0.5388E+01  0.6087E+01
   41   29   11    1.41 2018.21 -0.8625E+01  0.5597E+01
   49   11   11    1.41 2117.01 -0.9047E+01  0.8244E+01
   45   28    7    1.41 3878.69 -0.1658E+02  0.8841E+01
   22    5   25    1.41 1421.64 -0.6075E+01  0.7356E+01
   45    9   15    1.41 2102.95 -0.8987E+01  0.7519E+01
   51    2   10    1.41 3996.84 -0.1708E+02  0.7222E+01
   43    1   17    1.41 1583.11 -0.6765E+01  0.6466E+01
   27   26   20    1.41 2025.78 -0.8657E+01  0.6804E+01
   47   27    4    1.41 4537.63 -0.1939E+02  0.1284E+02
   29   20   21    1.41 1335.00 -0.5705E+01  0.5350E+01
   39   35    8    1.41 1590.88 -0.6799E+01  0.7388E+01
   30   13   22    1.41 1431.08 -0.6116E+01  0.8609E+01
   51    1   10    1.41 4213.29 -0.1801E+02  0.6881E+01
   39   38    3    1.41 1658.58 -0.7088E+01  0.7464E+01
   49    5   12    1.41 2823.20 -0.1207E+02  0.5533E+01
   47    3   14    1.41 3504.11 -0.1497E+02  0.9773E+01
   45   24   10    1.41 1095.71 -0.4683E+01  0.6949E+01
   43   24   12    1.41 1832.79 -0.7832E+01  0.8415E+01
   41    5   18    1.41 1470.30 -0.6283E+01  0.1044E+02
   50    4   11    1.41 4397.74 -0.1879E+02  0.7327E+01
   46   20   11    1.41 2297.35 -0.9818E+01  0.7191E+01
   46   10   14    1.41 1658.49 -0.7088E+01  0.9425E+01
   29   15   22    1.41 1416.13 -0.6052E+01  0.8617E+01
   25   21   22    1.41 1186.68 -0.5071E+01  0.5755E+01
   39    6   19    1.41 2432.40 -0.1039E+02  0.8743E+01
   38   20   17    1.41 1069.48 -0.4570E+01  0.7701E+01
   52   17    1    1.41 1944.09 -0.8308E+01  0.8910E+01
   42   33    6    1.41 1945.11 -0.8313E+01  0.7222E+01
   28   10   23    1.41 1936.08 -0.8274E+01  0.5386E+01
   54    8    2    1.41 1643.35 -0.7023E+01  0.8405E+01
   48   26    2    1.41 1286.40 -0.5497E+01  0.8729E+01
   52   16    3    1.41 3430.83 -0.1466E+02  0.1139E+02
   47    2   14    1.41 1164.49 -0.4976E+01  0.1066E+02
   38   23   16    1.41 2719.26 -0.1162E+02  0.6392E+01
   36   32   13    1.41 1768.89 -0.7559E+01  0.6644E+01
   44    6   16    1.41 1722.03 -0.7359E+01  0.1077E+02
   53   13    2    1.41 1842.21 -0.7873E+01  0.1205E+02
   53   11    4    1.41 2880.50 -0.1231E+02  0.1116E+02
   47   23    8    1.41 1575.18 -0.6732E+01  0.8965E+01
   37   35   10    1.41 1940.51 -0.8293E+01  0.7569E+01
   49   24    2    1.41 1617.88 -0.6914E+01  0.8725E+01
   52   15    4    1.41 1929.71 -0.8247E+01  0.1233E+02
   50   13    9    1.41 1355.67 -0.5794E+01  0.6093E+01
   45   22   11    1.41 2488.50 -0.1063E+02  0.6803E+01
   50    3   11    1.41 2028.85 -0.8670E+01  0.7342E+01
   49    4   12    1.41 1745.88 -0.7461E+01  0.6924E+01
   41    4   18    1.41 2192.43 -0.9369E+01  0.8728E+01
   35   18   19    1.41 3123.13 -0.1335E+02  0.6976E+01
   51   17    5    1.41 1583.26 -0.6766E+01  0.8930E+01
   35    2   21    1.41 3564.43 -0.1523E+02  0.8608E+01
   36   20   18    1.41 1454.14 -0.6214E+01  0.7543E+01
   29    6   23    1.41  953.23 -0.4074E+01  0.7378E+01
   53   12    3    1.41 2174.43 -0.9293E+01  0.1207E+02
   35   33   13    1.41 2466.48 -0.1054E+02  0.7997E+01
   48   19    9    1.41 2714.80 -0.1160E+02  0.6121E+01
   44   27    9    1.41 1285.11 -0.5492E+01  0.6200E+01
   50   22    1    1.41 3883.74 -0.1660E+02  0.1263E+02
   48    3   13    1.41 2407.48 -0.1029E+02  0.7322E+01
   21    7   25    1.41  975.08 -0.4167E+01  0.5177E+01
   47   26    5    1.41 2935.15 -0.1254E+02  0.1265E+02
   49   10   11    1.41 1125.74 -0.4811E+01  0.7879E+01
   54    8    1    1.41 3191.29 -0.1364E+02  0.8760E+01
   40   30   11    1.41 2001.92 -0.8555E+01  0.6014E+01
   42   14   16    1.41 1854.98 -0.7927E+01  0.6418E+01
   53   13    1    1.41 2910.26 -0.1244E+02  0.1227E+02
   32    5   22    1.41 1127.69 -0.4819E+01  0.8647E+01
   30   22   20    1.41 1289.73 -0.5512E+01  0.5531E+01
   49   24    1    1.41 1125.45 -0.4810E+01  0.8956E+01
   47   24    7    1.41 2156.98 -0.9218E+01  0.8735E+01
   42   31    8    1.41 1127.06 -0.4817E+01  0.6585E+01
   50   15    8    1.41 2816.71 -0.1204E+02  0.8848E+01
   22    1   25    1.41 1639.29 -0.7006E+01  0.6976E+01
   41   30   10    1.41 1126.01 -0.4812E+01  0.6387E+01
   50    9   10    1.41 1661.06 -0.7099E+01  0.8396E+01
   44   30    6    1.41 2906.32 -0.1242E+02  0.1217E+02
   48   10   12    1.41 1848.02 -0.7898E+01  0.5835E+01
   34   23   18    1.41 1971.75 -0.8426E+01  0.6097E+01
   41    2   18    1.41 1769.72 -0.7563E+01  0.6986E+01
   36   28   15    1.41 2094.77 -0.8952E+01  0.6469E+01
   39    4   19    1.41 2393.01 -0.1023E+02  0.8168E+01
   29    5   23    1.41 1124.76 -0.4807E+01  0.7293E+01
   31   16   21    1.41 2249.08 -0.9612E+01  0.5019E+01
   51   11    8    1.41 2640.10 -0.1128E+02  0.6101E+01
   28    9   23    1.41 1018.26 -0.4352E+01  0.5080E+01
   50   21    3    1.41 2044.86 -0.8739E+01  0.1297E+02
   48   23    6    1.41 1800.97 -0.7697E+01  0.8743E+01
   45   13   14    1.41 4665.44 -0.1994E+02  0.9166E+01
   44   32    2    1.42 3029.33 -0.1295E+02  0.8702E+01
   52    4    8    1.42 2118.74 -0.9055E+01  0.6340E+01
   37    3   20    1.42 3036.61 -0.1298E+02  0.7948E+01
   43   33    3    1.42 1012.55 -0.4327E+01  0.6606E+01
   47   27    3    1.42 1535.93 -0.6564E+01  0.1277E+02
   46   29    2    1.42 1386.71 -0.5926E+01  0.8645E+01
   39    3   19    1.42 2093.98 -0.8949E+01  0.1024E+02
   37    2   20    1.42 2010.92 -0.8594E+01  0.9259E+01
   54    7    1    1.42 1732.47 -0.7404E+01  0.1079E+02
   53   12    2    1.42 2544.14 -0.1087E+02  0.1172E+02
   53    4    6    1.42 1902.33 -0.8130E+01  0.1176E+02
   35   31   14    1.42 2057.35 -0.8792E+01  0.6183E+01
   54    4    3    1.42 3536.66 -0.1511E+02  0.8765E+01
   42   23   13    1.42 3758.58 -0.1606E+02  0.1158E+02
   38   25   15    1.42 1667.38 -0.7126E+01  0.5676E+01
   37    1   20    1.42 2947.63 -0.1260E+02  0.8867E+01
   30   25   19    1.42 1574.02 -0.6727E+01  0.5199E+01
   51   19    1    1.42 4743.20 -0.2027E+02  0.1149E+02
   53   11    3    1.42 4409.57 -0.1884E+02  0.1193E+02
   45   29    5    1.42 6292.82 -0.2689E+02  0.1172E+02
   51   13    7    1.42 2253.51 -0.9630E+01  0.1191E+02
   41   31    9    1.42 1231.50 -0.5263E+01  0.6803E+01
   43   21   13    1.42 1116.50 -0.4771E+01  0.8225E+01
   35   29   15    1.42 1038.43 -0.4438E+01  0.6116E+01
   39    1   19    1.42 1945.43 -0.8314E+01  0.7011E+01
   29   19   21    1.42 1045.99 -0.4470E+01  0.5033E+01
   34   32   14    1.42 3318.34 -0.1418E+02  0.7082E+01
   34   28   16    1.42 1233.15 -0.5270E+01  0.5525E+01
   27   11   23    1.42 1358.95 -0.5808E+01  0.5272E+01
   53    3    6    1.42 4867.45 -0.2080E+02  0.1368E+02
   51   18    3    1.42 1782.85 -0.7619E+01  0.8846E+01
   45   19   12    1.42 1605.39 -0.6861E+01  0.7803E+01
   44   31    4    1.42 1595.24 -0.6817E+01  0.9067E+01
   25    6   24    1.42 2321.35 -0.9920E+01  0.9104E+01
   36   33   12    1.42 1094.95 -0.4679E+01  0.5867E+01
   46   13   13    1.42 1765.30 -0.7544E+01  0.1109E+02
   38   29   13    1.42 1825.81 -0.7803E+01  0.9903E+01
   48   16   10    1.42 1817.04 -0.7765E+01  0.7834E+01
   44    1   16    1.42 1191.08 -0.5090E+01  0.6666E+01
   29   26   19    1.42 1712.33 -0.7318E+01  0.5054E+01
   49   19    7    1.42 1705.23 -0.7287E+01  0.1191E+02
   25   20   22    1.42 1735.44 -0.7416E+01  0.4287E+01
   36    8   20    1.42 6068.98 -0.2594E+02  0.9002E+01
   54    5    2    1.42 1942.93 -0.8303E+01  0.8729E+01
   51   16    5    1.42 3347.42 -0.1431E+02  0.8631E+01
   35   34   12    1.42 1408.01 -0.6017E+01  0.7850E+01
   41   11   17    1.42 1734.82 -0.7414E+01  0.7140E+01
   45   16   13    1.42 3199.07 -0.1367E+02  0.1162E+02
   44   11   15    1.42 5822.19 -0.2488E+02  0.1113E+02
   39   37    4    1.42  978.61 -0.4182E+01  0.6316E+01
   42    6   17    1.42 3017.85 -0.1290E+02  0.8691E+01
   21   20   23    1.42 1213.39 -0.5185E+01  0.3930E+01
   47   13   12    1.42 2045.93 -0.8743E+01  0.5878E+01
   40   34    7    1.42 2185.13 -0.9338E+01  0.7317E+01
   23   11   24    1.42 1457.14 -0.6227E+01  0.5436E+01
   49   22    4    1.42 1639.77 -0.7008E+01  0.9055E+01
   33   31   15    1.42 2160.55 -0.9233E+01  0.7223E+01
   41   28   11    1.42 2747.22 -0.1174E+02  0.7668E+01
   54    4    2    1.42 1081.14 -0.4620E+01  0.8440E+01
   42   20   14    1.42 1069.55 -0.4571E+01  0.6798E+01
   32   30   16    1.42 1246.01 -0.5325E+01  0.5609E+01
   28    7   23    1.42 2977.11 -0.1272E+02  0.7455E+01
   54    5    1    1.42 5965.16 -0.2549E+02  0.8726E+01
   52   15    2    1.42 1382.19 -0.5907E+01  0.8344E+01
   46   27    5    1.42 4665.32 -0.1994E+02  0.8948E+01
   53    6    5    1.42 1872.94 -0.8004E+01  0.8727E+01
   37   30   13    1.42 1665.93 -0.7119E+01  0.8365E+01
   37   26   15    1.42 3142.63 -0.1343E+02  0.7729E+01
   39   20   16    1.42 1535.98 -0.6564E+01  0.6365E+01
   42    5   17    1.42 3671.98 -0.1569E+02  0.1058E+02
   40    7   18    1.42 1694.47 -0.7241E+01  0.1042E+02
   41   35    3    1.42 2486.16 -0.1062E+02  0.8343E+01
   49   21    5    1.42  946.59 -0.4045E+01  0.8791E+01
   34    5   21    1.42 1209.31 -0.5168E+01  0.8512E+01
   51   18    2    1.42 4320.76 -0.1846E+02  0.1165E+02
   45   30    2    1.42 1574.47 -0.6729E+01  0.1215E+02
   41   19   15    1.42 2283.85 -0.9760E+01  0.1123E+02
   24    8   24    1.42 1330.40 -0.5686E+01  0.5785E+01
   30   21   20    1.42 3243.81 -0.1386E+02  0.5316E+01
   43   32    4    1.42 2445.06 -0.1045E+02  0.1262E+02
   54    4    1    1.42 1504.75 -0.6431E+01  0.8624E+01
   42   34    2    1.42 2691.57 -0.1150E+02  0.5900E+01
   46   28    3    1.42 3060.49 -0.1308E+02  0.8985E+01
   50   20    3    1.42 5115.49 -0.2186E+02  0.1213E+02
   44   30    5    1.42 1025.20 -0.4381E+01  0.8814E+01
   37   28   14    1.42 1167.50 -0.4989E+01  0.6399E+01
   44   10   15    1.43 1035.54 -0.4425E+01  0.1114E+02
   38    7   19    1.43 1298.06 -0.5547E+01  0.1002E+02
   49   23    1    1.43 1271.86 -0.5435E+01  0.8782E+01
   53    5    5    1.43 1409.88 -0.6025E+01  0.8584E+01
   52   15    1    1.43 1825.23 -0.7800E+01  0.8522E+01
   47   20    9    1.43 1336.15 -0.5710E+01  0.7012E+01
   28    6   23    1.43 1580.90 -0.6756E+01  0.7429E+01
   40   36    3    1.43 1838.36 -0.7856E+01  0.6163E+01
   41   24   13    1.43 1774.46 -0.7583E+01  0.1112E+02
   42   32    6    1.43 3299.93 -0.1410E+02  0.6860E+01
   37   11   19    1.43 1317.46 -0.5630E+01  0.7004E+01
   39   16   17    1.43 1825.70 -0.7802E+01  0.6173E+01
   45   29    4    1.43  974.25 -0.4164E+01  0.1320E+02
   31    6   22    1.43 1521.34 -0.6502E+01  0.8095E+01
   23   10   24    1.43 1601.29 -0.6843E+01  0.5169E+01
   51    2    9    1.43 2819.97 -0.1205E+02  0.7396E+01
   48   15   10    1.43 1530.69 -0.6541E+01  0.7732E+01
   47   12   12    1.43 1193.10 -0.5099E+01  0.6693E+01
   53    9    3    1.43 5757.79 -0.2461E+02  0.1195E+02
   51   17    3    1.43 1141.39 -0.4878E+01  0.1190E+02
   51   15    5    1.43 1872.47 -0.8002E+01  0.8578E+01
   29   18   21    1.43 1551.00 -0.6628E+01  0.4970E+01
   24    7   24    1.43 1295.12 -0.5535E+01  0.8161E+01
   52    5    7    1.43 2160.05 -0.9231E+01  0.5745E+01
   50   13    8    1.43 4018.82 -0.1717E+02  0.8487E+01
   37   34   10    1.43 2057.90 -0.8795E+01  0.6294E+01
   50    5   10    1.43 1625.22 -0.6945E+01  0.5429E+01
   42   34    1    1.43 2196.97 -0.9389E+01  0.6013E+01
   35   32   13    1.43 1612.48 -0.6891E+01  0.8382E+01
   46    5   14    1.43 1849.10 -0.7902E+01  0.6355E+01
   46   22    9    1.43 1830.83 -0.7824E+01  0.7937E+01
   35   16   19    1.43 1115.33 -0.4766E+01  0.6921E+01
   29   22   20    1.43  971.80 -0.4153E+01  0.5395E+01
   34   13   20    1.43 3025.49 -0.1293E+02  0.6361E+01
   53    7    4    1.43 3163.97 -0.1352E+02  0.1143E+02
   38   18   17    1.43 1801.73 -0.7700E+01  0.6264E+01
   41   29   10    1.43 2234.64 -0.9550E+01  0.8086E+01
   49   22    3    1.43 2446.42 -0.1045E+02  0.8818E+01
   51   16    4    1.43 3463.89 -0.1480E+02  0.1030E+02
   44   29    6    1.43 1590.95 -0.6799E+01  0.8668E+01
   34    2   21    1.43 2821.97 -0.1206E+02  0.8307E+01
   49   14    9    1.43 2766.34 -0.1182E+02  0.8018E+01
   36   35   10    1.43 1212.83 -0.5183E+01  0.7581E+01
   49    6   11    1.43 1487.80 -0.6358E+01  0.4856E+01
   37   23   16    1.43 1894.96 -0.8098E+01  0.7554E+01
   47    6   13    1.43 1327.26 -0.5672E+01  0.6181E+01
   42   30    8    1.43 1692.29 -0.7232E+01  0.8839E+01
   18   10   25    1.43 1096.33 -0.4685E+01  0.4147E+01
   35   20   18    1.43 2334.82 -0.9978E+01  0.9829E+01
   25   19   22    1.43 1281.40 -0.5476E+01  0.4004E+01
   31    5   22    1.43  977.59 -0.4178E+01  0.7874E+01
   36    5   20    1.43 1367.99 -0.5846E+01  0.7925E+01
   53    3    5    1.43 2175.16 -0.9296E+01  0.8835E+01
   34    1   21    1.43 1722.89 -0.7363E+01  0.1147E+02
   29   12   22    1.43 1475.52 -0.6306E+01  0.4471E+01
   47   15   11    1.43 2454.84 -0.1049E+02  0.7751E+01
   39   36    5    1.43 2012.60 -0.8601E+01  0.7611E+01
   43   13   15    1.43 2156.45 -0.9216E+01  0.1104E+02
   46   28    2    1.43 3072.71 -0.1313E+02  0.1171E+02
   25   23   21    1.43 1156.88 -0.4944E+01  0.4415E+01
   39   10   18    1.43 3019.84 -0.1291E+02  0.8380E+01
   21   19   23    1.43  995.99 -0.4256E+01  0.3658E+01
   33    8   21    1.43 4701.47 -0.2009E+02  0.8736E+01
   30    9   22    1.43 4126.17 -0.1763E+02  0.8296E+01
   44   31    2    1.43 1189.68 -0.5084E+01  0.8648E+01
   38   37    5    1.43 1554.89 -0.6645E+01  0.7517E+01
   49   16    8    1.43 1781.38 -0.7613E+01  0.6022E+01
   28   14   22    1.43 1196.68 -0.5114E+01  0.4447E+01
   52    8    6    1.43 2835.91 -0.1212E+02  0.1042E+02
   24    6   24    1.43 2591.20 -0.1107E+02  0.8688E+01
   33   15   20    1.43 2077.73 -0.8879E+01  0.6578E+01
   38    5   19    1.43 1187.70 -0.5076E+01  0.1018E+02
   41   35    1    1.43 4442.77 -0.1899E+02  0.8488E+01
   23    9   24    1.43 2100.49 -0.8977E+01  0.8812E+01
   45   20   11    1.43 2594.73 -0.1109E+02  0.5877E+01
   35   30   14    1.43 1401.70 -0.5990E+01  0.4909E+01
   42   19   14    1.43 1017.30 -0.4347E+01  0.6756E+01
   42   11   16    1.43 3117.47 -0.1332E+02  0.7730E+01
   34   27   16    1.43 1481.58 -0.6332E+01  0.6153E+01
   36   13   19    1.43 1167.43 -0.4989E+01  0.5752E+01
   51   17    2    1.43 2036.56 -0.8703E+01  0.1224E+02
   39   37    2    1.43 3154.23 -0.1348E+02  0.8460E+01
   49   19    6    1.43 1734.97 -0.7414E+01  0.1096E+02
   32   11   21    1.43 1241.69 -0.5306E+01  0.5988E+01
   17   11   25    1.43 1013.36 -0.4331E+01  0.3954E+01
   37   31   12    1.43 1516.92 -0.6483E+01  0.6817E+01
   50    9    9    1.43 1174.40 -0.5019E+01  0.5452E+01
   52   14    1    1.43 1296.32 -0.5540E+01  0.8596E+01
   40   27   12    1.43  948.44 -0.4053E+01  0.8434E+01
   39   22   15    1.43 1682.60 -0.7191E+01  0.6493E+01
   44   22   11    1.43 1028.49 -0.4395E+01  0.8312E+01
   23   21   22    1.43 1458.41 -0.6233E+01  0.4086E+01
   36    3   20    1.43 1277.36 -0.5459E+01  0.8417E+01
   43   27    9    1.43 1589.75 -0.6794E+01  0.6196E+01
   47   26    2    1.43 1800.70 -0.7695E+01  0.8747E+01
   49   20    5    1.43 2835.94 -0.1212E+02  0.8951E+01
   41   32    7    1.43 2131.55 -0.9109E+01  0.6792E+01
   18    9   25    1.43 2603.06 -0.1112E+02  0.5362E+01
   49   10   10    1.43  992.59 -0.4242E+01  0.5954E+01
   40   25   13    1.43 3437.08 -0.1469E+02  0.6824E+01
   40   30   10    1.43 1118.63 -0.4781E+01  0.6709E+01
   33    7   21    1.43 2951.08 -0.1261E+02  0.8604E+01
   23   16   23    1.43 1017.46 -0.4348E+01  0.3763E+01
   51   14    5    1.43 3030.29 -0.1295E+02  0.8761E+01
   48   24    2    1.43 4098.94 -0.1752E+02  0.1199E+02
   39    9   18    1.43 1980.33 -0.8463E+01  0.1027E+02
   16   12   25    1.43 2768.49 -0.1183E+02  0.5510E+01
   40    1   18    1.43  964.13 -0.4120E+01  0.5977E+01
   36   24   16    1.43 2268.81 -0.9696E+01  0.5986E+01
   51   17    1    1.43 2061.12 -0.8808E+01  0.1241E+02
   53    7    3    1.44 2353.80 -0.1006E+02  0.1201E+02
   47    3   13    1.44 1001.16 -0.4279E+01  0.7647E+01
   23    8   24    1.44  952.86 -0.4072E+01  0.5945E+01
   37   25   15    1.44 1505.06 -0.6432E+01  0.5645E+01
   43   12   15    1.44 2618.17 -0.1119E+02  0.7476E+01
   29   17   21    1.44 1675.82 -0.7162E+01  0.4848E+01
   27    7   23    1.44 1057.91 -0.4521E+01  0.7376E+01
   39   35    6    1.44 2516.95 -0.1076E+02  0.8085E+01
   27   20   21    1.44 1690.93 -0.7226E+01  0.4892E+01
   49   22    1    1.44 2832.71 -0.1211E+02  0.1178E+02
   52   13    2    1.44 1379.38 -0.5895E+01  0.8379E+01
   45   17   12    1.44 1128.85 -0.4824E+01  0.1105E+02
   38   31   11    1.44 2852.16 -0.1219E+02  0.5903E+01
   49    2   11    1.44  965.14 -0.4125E+01  0.7448E+01
   48   22    5    1.44 4408.36 -0.1884E+02  0.1233E+02
   47    2   13    1.44 3328.19 -0.1422E+02  0.7811E+01
   46   24    7    1.44 3447.65 -0.1473E+02  0.6142E+01
   46   14   12    1.44 5092.47 -0.2176E+02  0.1153E+02
   34   24   17    1.44 1880.07 -0.8035E+01  0.5375E+01
   28   13   22    1.44 1932.56 -0.8259E+01  0.5753E+01
   52    9    5    1.44 1959.22 -0.8373E+01  0.1174E+02
   43   19   13    1.44 1597.47 -0.6827E+01  0.5758E+01
   41   23   13    1.44 1620.95 -0.6927E+01  0.8132E+01
   32   10   21    1.44  991.28 -0.4236E+01  0.5722E+01
   48   24    1    1.44  955.56 -0.4084E+01  0.1134E+02
   38   36    6    1.44 1569.52 -0.6707E+01  0.5969E+01
   48   18    8    1.44 3673.31 -0.1570E+02  0.8021E+01
   48    2   12    1.44 2958.87 -0.1264E+02  0.6788E+01
   51    5    8    1.44 2104.20 -0.8992E+01  0.6037E+01
   41   25   12    1.44 2298.87 -0.9824E+01  0.8284E+01
   39   24   14    1.44 1432.32 -0.6121E+01  0.6156E+01
   51    9    7    1.44  956.07 -0.4086E+01  0.6881E+01
   21   18   23    1.44 1299.34 -0.5553E+01  0.3533E+01
   43   32    1    1.44 1800.61 -0.7695E+01  0.8720E+01
   53    8    1    1.44 1644.27 -0.7027E+01  0.1081E+02
   50   19    2    1.44 2780.94 -0.1188E+02  0.8468E+01
   46   19   10    1.44 1128.67 -0.4823E+01  0.5995E+01
   40   11   17    1.44 2155.08 -0.9210E+01  0.7051E+01
   53    7    2    1.44 5078.17 -0.2170E+02  0.1114E+02
   43   31    4    1.44 3027.57 -0.1294E+02  0.1148E+02
   43   25   10    1.44 1451.24 -0.6202E+01  0.7786E+01
   50    7    9    1.44 2234.90 -0.9551E+01  0.4813E+01
   19    4   25    1.44 1021.85 -0.4367E+01  0.4501E+01
   48   13   10    1.44 2863.03 -0.1224E+02  0.7669E+01
   41   13   16    1.44 2474.67 -0.1058E+02  0.7597E+01
   44   19   12    1.44 3007.22 -0.1285E+02  0.8244E+01
   39   26   13    1.44 1443.21 -0.6168E+01  0.6654E+01
   45    8   14    1.44 2718.89 -0.1162E+02  0.7837E+01
   42   28    9    1.44 1826.54 -0.7806E+01  0.8552E+01
   38   23   15    1.44 2219.68 -0.9486E+01  0.5589E+01
   39   14   17    1.44 2364.46 -0.1010E+02  0.7380E+01
   41    6   17    1.44 2900.28 -0.1239E+02  0.8902E+01
   36   30   13    1.44 3263.93 -0.1395E+02  0.5294E+01
   44    6   15    1.44 2793.97 -0.1194E+02  0.7698E+01
   45   29    1    1.44 1150.11 -0.4915E+01  0.8790E+01
   47   25    3    1.44 1157.22 -0.4945E+01  0.8809E+01
   53    5    3    1.44 1640.09 -0.7009E+01  0.1227E+02
   41   33    5    1.44 1739.62 -0.7434E+01  0.7141E+01
   24   19   22    1.44 2001.10 -0.8552E+01  0.5749E+01
   45   19   11    1.44 2101.78 -0.8982E+01  0.5842E+01
   45   13   13    1.44 1098.54 -0.4695E+01  0.7566E+01
   50   15    6    1.44 1344.39 -0.5745E+01  0.8458E+01
   41   17   15    1.44 1221.06 -0.5218E+01  0.7531E+01
   48    9   11    1.44 1049.85 -0.4487E+01  0.6558E+01
   36   33   11    1.44 2611.09 -0.1116E+02  0.7696E+01
   52    8    5    1.44 1270.42 -0.5429E+01  0.6330E+01
   39   18   16    1.44 1058.66 -0.4524E+01  0.7375E+01
   42    9   16    1.44 1838.00 -0.7855E+01  0.5511E+01
   33   22   18    1.44 1229.71 -0.5255E+01  0.5433E+01
   38   20   16    1.44  991.06 -0.4235E+01  0.7394E+01
   31   12   21    1.44 2279.86 -0.9743E+01  0.5894E+01
   19    3   25    1.44 1644.91 -0.7030E+01  0.4336E+01
   32   20   19    1.44 1846.02 -0.7889E+01  0.5621E+01
   50   19    1    1.44 3004.45 -0.1284E+02  0.8765E+01
   47    9   12    1.44 2128.93 -0.9098E+01  0.1143E+02
   49   20    4    1.44 3487.74 -0.1491E+02  0.9107E+01
   26   16   22    1.44 1018.81 -0.4354E+01  0.5549E+01
   15   12   25    1.44 1331.74 -0.5691E+01  0.4429E+01
   44   23   10    1.44 1037.10 -0.4432E+01  0.6768E+01
   47   16   10    1.44 2196.65 -0.9387E+01  0.8331E+01
   36   28   14    1.44 2594.13 -0.1109E+02  0.7685E+01
   23   20   22    1.44 1422.41 -0.6079E+01  0.5559E+01
   46   27    2    1.44 3814.92 -0.1630E+02  0.1190E+02
   53    6    2    1.44 2380.19 -0.1017E+02  0.1137E+02
   40   35    3    1.44 3133.44 -0.1339E+02  0.7167E+01
   44   27    7    1.44 3217.04 -0.1375E+02  0.6175E+01
   37   36    7    1.44 3139.97 -0.1342E+02  0.7975E+01
   42   29    8    1.44 1236.20 -0.5283E+01  0.6175E+01
   50   18    3    1.44 2337.72 -0.9990E+01  0.8934E+01
   38   29   12    1.44 2519.29 -0.1077E+02  0.8307E+01
   41    5   17    1.44 3170.32 -0.1355E+02  0.8686E+01
   42   30    7    1.44 1237.87 -0.5290E+01  0.6159E+01
   40   19   15    1.44  967.30 -0.4134E+01  0.6536E+01
   38   11   18    1.44 1868.37 -0.7985E+01  0.6754E+01
   30   19   20    1.44 1303.68 -0.5571E+01  0.5324E+01
   35    6   20    1.44 2896.74 -0.1238E+02  0.8029E+01
   49   21    2    1.44 1432.86 -0.6123E+01  0.1240E+02
   51    1    8    1.44 2105.70 -0.8999E+01  0.7207E+01
   46   26    4    1.44 1847.97 -0.7897E+01  0.1279E+02
   40   34    5    1.44 2798.23 -0.1196E+02  0.8125E+01
   46   22    8    1.44 1047.84 -0.4478E+01  0.6100E+01
   50   10    8    1.44 1775.92 -0.7589E+01  0.8261E+01
   33   18   19    1.44 2338.80 -0.9995E+01  0.7489E+01
   39   36    3    1.45 2586.47 -0.1105E+02  0.7363E+01
   52    7    5    1.45 1377.93 -0.5889E+01  0.6483E+01
   49   16    7    1.45 2683.55 -0.1147E+02  0.8441E+01
   49   14    8    1.45 2404.38 -0.1028E+02  0.6761E+01
   52    2    6    1.45 1817.40 -0.7767E+01  0.7624E+01
   40   26   12    1.45 2748.28 -0.1174E+02  0.9612E+01
   42    8   16    1.45 2280.21 -0.9745E+01  0.7176E+01
   18   15   24    1.45 1262.17 -0.5394E+01  0.2829E+01
   50    5    9    1.45 1422.16 -0.6078E+01  0.6425E+01
   22    8   24    1.45 2015.45 -0.8613E+01  0.3990E+01
   48   12   10    1.45 1000.33 -0.4275E+01  0.7663E+01
   36   23   16    1.45 1931.25 -0.8253E+01  0.7275E+01
   35    5   20    1.45 1815.68 -0.7759E+01  0.9052E+01
   32   23   18    1.45 2881.85 -0.1232E+02  0.5391E+01
   27   26   19    1.45 2278.13 -0.9736E+01  0.6000E+01
   47   21    7    1.45 1711.55 -0.7314E+01  0.8511E+01
   51   12    5    1.45 2160.96 -0.9235E+01  0.6040E+01
   48   21    5    1.45 1533.92 -0.6555E+01  0.8687E+01
   36   35    9    1.45 2251.25 -0.9621E+01  0.6882E+01
   37   30   12    1.45 1213.31 -0.5185E+01  0.7067E+01
   52   11    2    1.45 1555.75 -0.6649E+01  0.8233E+01
   53    4    2    1.45 1272.97 -0.5440E+01  0.1184E+02
   44   29    4    1.45 1192.19 -0.5095E+01  0.9262E+01
   41   22   13    1.45 3537.88 -0.1512E+02  0.8141E+01
   39    5   18    1.45 1971.77 -0.8426E+01  0.8101E+01
   46   18   10    1.45  959.28 -0.4100E+01  0.5904E+01
   42   20   13    1.45 1364.20 -0.5830E+01  0.7696E+01
   34   30   14    1.45 2328.29 -0.9950E+01  0.6097E+01
   34   28   15    1.45 3616.27 -0.1545E+02  0.6187E+01
   35    4   20    1.45 1224.60 -0.5233E+01  0.9156E+01
   41   31    7    1.45 1321.82 -0.5649E+01  0.8181E+01
   29    8   22    1.45 1133.70 -0.4845E+01  0.7763E+01
   50    9    8    1.45 1985.94 -0.8487E+01  0.6227E+01
   44   15   13    1.45  981.00 -0.4192E+01  0.5693E+01
   46   24    6    1.45 1452.47 -0.6207E+01  0.8351E+01
   42   17   14    1.45 3010.48 -0.1287E+02  0.1096E+02
   41   16   15    1.45 2076.96 -0.8876E+01  0.6271E+01
   37    5   19    1.45 3124.23 -0.1335E+02  0.9873E+01
   46   12   12    1.45 2025.20 -0.8655E+01  0.8058E+01
   44   18   12    1.45 1092.22 -0.4668E+01  0.5864E+01
   40    9   17    1.45 1007.03 -0.4304E+01  0.7185E+01
   53    3    2    1.45 3306.76 -0.1413E+02  0.1175E+02
   41   33    4    1.45 1473.61 -0.6298E+01  0.7073E+01
   28   17   21    1.45 1224.00 -0.5231E+01  0.8207E+01
   32    7   21    1.45 1707.41 -0.7297E+01  0.8232E+01
   43   27    8    1.45  993.67 -0.4246E+01  0.8746E+01
   26   15   22    1.45 1552.52 -0.6635E+01  0.4559E+01
   51   14    3    1.45 1307.35 -0.5587E+01  0.1204E+02
   40   33    6    1.45 1090.81 -0.4662E+01  0.7962E+01
   45    5   14    1.45 1962.68 -0.8388E+01  0.7866E+01
   39   23   14    1.45 2372.18 -0.1014E+02  0.9240E+01
   45   22    9    1.45 1883.70 -0.8050E+01  0.5974E+01
   37   21   16    1.45 2382.02 -0.1018E+02  0.7681E+01
   24   12   23    1.45 1101.78 -0.4708E+01  0.3806E+01
   36   31   12    1.45 2063.70 -0.8819E+01  0.7001E+01
   37   13   18    1.45 2824.04 -0.1207E+02  0.5838E+01
   31   24   18    1.45 1891.00 -0.8081E+01  0.5588E+01
   33   12   20    1.45 1301.40 -0.5562E+01  0.5938E+01
   45   28    2    1.45 1088.06  0.1842E+02  0.1029E+02
   48   11   10    1.45 1273.93  0.2157E+02  0.8272E+01
   35   32   12    1.45 1033.33  0.1750E+02  0.7930E+01
   40   21   14    1.45 1508.57  0.2554E+02  0.1132E+02
   42   22   12    1.45 1150.05  0.1947E+02  0.8330E+01
   40   27   11    1.46  959.34  0.1624E+02  0.1181E+02
   51    2    7    1.46 1160.93  0.1966E+02  0.7606E+01
   44   19   11    1.46 1383.88  0.2343E+02  0.9727E+01
   32   27   16    1.47 1215.72  0.2059E+02  0.5543E+01
   43   27    7    1.47  953.55  0.1615E+02  0.6286E+01
   47    9   11    1.47 2140.21  0.3624E+02  0.1126E+02
   43   30    2    1.47 1044.50  0.1769E+02  0.8674E+01
   43   28    5    1.48 1135.72  0.1923E+02  0.8735E+01
   44   11   13    1.48 1049.55  0.1777E+02  0.1091E+02
   34   33   11    1.48 1111.50  0.1882E+02  0.7497E+01
   39   25   12    1.48  962.79  0.1630E+02  0.6631E+01
   45   26    1    1.48  983.37  0.1665E+02  0.6381E+01
   45   23    6    1.49 1418.16  0.2401E+02  0.8496E+01
   49   16    3    1.49 1235.41  0.2092E+02  0.6312E+01
   35   31   11    1.49 1268.73  0.2148E+02  0.9072E+01
   37   20   15    1.49  945.41  0.1601E+02  0.7288E+01
   46   22    4    1.50 1018.90  0.1725E+02  0.8030E+01
   35   12   18    1.50  951.94  0.1612E+02  0.6542E+01
   34   22   16    1.50 1240.62  0.2101E+02  0.1355E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT      5147
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"       91311
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)   69927
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))    21384


 cpu time used                  5.4 sec
 elapsed wall-clock time        2.4 sec
