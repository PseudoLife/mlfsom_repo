
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.539    77.539    38.796  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      30      30
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      30
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000940  0.000681
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
 BEAM_DIVERGENCE_E.S.D.=   0.046
 REFLECTING_RANGE_E.S.D.=   0.150
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

    1275 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      30 ...      30
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1275 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               41
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

 *  44        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  31        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  23        oC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  35        mP          0.0      77.5   38.8   77.5  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  20        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  33        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  21        tP          0.0      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
 *  32        oP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    37        mC        250.0     159.9   38.8   77.5  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.9   38.8   77.5  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  159.9   77.5  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.9   77.5  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  159.9   77.5  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.8  159.9   77.5  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    19        oI        500.0      38.8  109.7  116.3  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    27        mC        500.0     159.9   38.8  109.7  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    18        tI        625.0     109.7  116.3   38.8  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    26        oF        625.0      38.8  159.9  159.9  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    22        hP        998.6      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    30        mC        998.6      77.5  173.4   38.8  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.6      77.5  173.4   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    41        mC        998.6     173.4   77.5   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.3  116.3  116.3  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    1  0  1  0  0 -1 -1  0  1  1 -1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.7   86.7  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.7   86.7   86.7  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.7   86.7  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.7   86.7  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.7  248.2 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.7   86.7   77.5  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.5  177.7  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.7   86.7  177.7 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.7   86.7   86.7  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.7  109.7   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.5  177.7 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.7  177.7   77.5 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    4.6      17    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.897   -23.483     0.054
 UNIT_CELL_B-AXIS=    23.483    73.897     0.023
 UNIT_CELL_C-AXIS=    -0.029    -0.003    38.797
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       167 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.23
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.437    77.437    38.749  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  7.4E-02 7.4E-02 3.1E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012914  0.012914  0.025807  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.808   -23.428     0.048
 COORDINATES OF UNIT CELL B-AXIS    23.428    73.808     0.016
 COORDINATES OF UNIT CELL C-AXIS    -0.025    -0.001    38.749
 CRYSTAL MOSAICITY (DEGREES)     0.150
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999975  0.006815  0.001842
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.93   1232.92
 DETECTOR ORIGIN (PIXELS) AT                     1261.25   1242.83
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.05
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.44    77.44    38.75  90.000  90.000  90.000


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
  0.05  0.15        46         71.018       1.742
  0.15  0.25       166        137.043       0.903
  0.25  0.35       175        124.281       0.996
  0.35  0.45       159         61.111       2.025
  0.45  0.55       166        218.825       0.565
  0.55  0.65       181         88.421       1.399
  0.65  0.75       198        134.305       0.921
  0.75  0.85       139        123.656       1.001
  0.85  0.95        45         92.884       1.332
 
 -0.05  0.95      1275        123.742       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               209.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.927
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    29.1 XMAX=    29.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59480 NYBIN=    9
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        452


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               209.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.923
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    41.8 XMAX=  2453.4 NXBIN=    3
 YMIN=     8.5 YMAX=  2519.2 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        452


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               209.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.918
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    29.1 XMAX=    29.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1248    1264
 DETECTOR_SURFACE_POSITION=    1657    1690
 DETECTOR_SURFACE_POSITION=     838    1690
 DETECTOR_SURFACE_POSITION=     838     837
 DETECTOR_SURFACE_POSITION=    1657     837
 DETECTOR_SURFACE_POSITION=    2175    1664
 DETECTOR_SURFACE_POSITION=    1632    2229
 DETECTOR_SURFACE_POSITION=     864    2229
 DETECTOR_SURFACE_POSITION=     321    1664
 DETECTOR_SURFACE_POSITION=     321     864
 DETECTOR_SURFACE_POSITION=     864     299
 DETECTOR_SURFACE_POSITION=    1632     299
 DETECTOR_SURFACE_POSITION=    2175     864
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        452


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
 8.529E-01  1.245E-03   30.69


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

   24.488   3.154    21.92   1.30      2.14      2.04      64      87      20
    3.154   2.240     7.33   1.18      7.75      6.56     174     199      12
    2.240   1.831     0.87   0.98     57.88     57.01     198     221       0
    1.831   1.587     0.09   0.82    510.74    548.96     236     263       2
    1.587   1.420     0.03   0.96    -99.90    -99.90     252     301       2
    1.420   1.297    -0.01   0.88    -99.90    -99.90     120     163       4
   --------------------------------------------------------------------------
   24.488   1.297     2.91   0.98      8.49      8.36    1044    1234      40
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      30      30)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.683  -0.1641E+02  0.1218E+02        -1.35           1 
   47    0    0     1.648  -0.2306E+00  0.1198E+02        -0.02           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
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

     3.85          56      35      1285        2.7%       2.0%      2.0%       42   29.35      2.8%    99.9*     0    0.000       0
     2.74          79      47      2098        2.2%       4.4%      3.5%       64   18.28      6.2%    99.4*     0    0.000       0
     2.24         151      85      2665        3.2%       9.1%      7.9%      132    8.55     12.9%    98.4*     0    0.000       0
     1.94         142      78      3101        2.5%      40.6%     37.7%      128    1.83     57.4%    81.1*     0    0.000       0
     1.74         157      86      3480        2.5%     -99.9%    -99.9%      142    0.00    -99.9%     4.3      0    0.000       0
     1.59         185     103      3850        2.7%     275.9%    292.6%      164    0.21    390.1%     2.5      0    0.000       0
     1.47         191     112      4143        2.7%     -99.9%    -99.9%      158    0.00    -99.9%     2.7      0    0.000       0
     1.38         194     113      4428        2.6%     719.1%    667.5%      162    0.11   1016.9%   -12.0      0    0.000       0
     1.30          79      53      4742        1.1%     -99.9%    -99.9%       52    0.00    -99.9%     9.0      0    0.000       0
    total        1234     712     29792        2.4%       8.5%      8.4%     1044    3.90     12.0%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1275
 NUMBER OF REJECTED MISFITS                                40
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1234
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    712



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

   WILSON LINE (using all data) : A=   7.002 B=  29.620 CORRELATION=  0.90
      #      RES      SS        <I>       log(<I>)       BO
      14     7.394   0.005  1.9450E+03       7.573     -62.4
      17     4.702   0.011  1.9225E+03       7.561     -24.7
      22     3.548   0.020  7.4031E+02       6.607       9.9
      24     2.994   0.028  2.0567E+02       5.326      30.0
      32     2.633   0.036  1.6517E+02       5.107      26.3
      47     2.402   0.043  8.1515E+01       4.401      30.0
      33     2.205   0.051  3.5729E+01       3.576      33.3
      39     2.055   0.059  1.5107E+01       2.715      36.2
      38     1.934   0.067 -1.0066E+00    -999.999   -1000.0
      47     1.829   0.075 -1.1956E+00    -999.999   -1000.0
      42     1.739   0.083  1.6143E+00       0.479      39.5
      43     1.668   0.090  2.3185E+00       0.841      34.3
      48     1.598   0.098  9.6980E-01      -0.031      35.9
      53     1.537   0.106  2.9564E-01      -1.219      38.9
      66     1.482   0.114  2.2452E-01      -1.494      37.3
      51     1.434   0.122  2.3686E+00       0.862      25.2
      44     1.392   0.129 -3.2728E-01    -999.999   -1000.0
      36     1.352   0.137  1.6923E+00       0.526      23.7
      16     1.315   0.145  3.6378E+00       1.291      19.8


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.394        1.916        1.379        0.776
       3     4.702        1.424        0.929        0.450
       2     3.548        0.411        0.091        0.015
       3     2.994        2.163        1.166        0.464
       2     2.633        0.468        0.155        0.037
       5     2.402        1.053        0.597        0.266
       1     2.205        1.693        0.763        0.246
       6     2.055        0.362        0.132        0.040
      10     1.829       20.018       20.225       79.154
       7     1.739        7.866        1.075        9.532
       3     1.598       83.054      409.277     1511.352
       3     1.537      200.350      -64.552     4508.027
       5     1.482      156.743      732.628     2531.804
       1     1.392       35.845       74.342      110.132
       1     1.352       33.822       68.137       98.049
      76   overall       24.009       62.517      407.392


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.394        0.535        0.391        0.225
      14     4.702        1.276        1.623        1.634
      20     3.548        0.958        0.915        0.800
      21     2.994        0.760        0.650        0.464
      30     2.633        0.970        0.805        0.569
      42     2.402        0.947        0.852        0.650
      32     2.205        0.847        0.685        0.512
      33     2.055        1.480        1.901        2.395
      32     1.934       39.278       73.699      880.510
      35     1.739       11.772       16.691       69.209
      32     1.668       10.083       13.355       48.415
      45     1.598       48.453      124.596     1483.012
      50     1.537      590.274     2271.635   193856.647
      48     1.434        8.139       17.797       55.240
      35     1.352       25.351       90.241      492.535
      16     1.315        7.625       11.062       28.886
     636   overall      136.661       64.600    40583.560

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.394  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.73
       3     4.702  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.548  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.994  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.633  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.402  0.00  0.00  0.00  0.00  0.38  0.71  0.67  0.64  0.91  0.88
       1     2.205  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.055  0.67  0.96  0.80  0.70  0.96  0.89  0.84  1.33  1.27  1.22
      10     1.829  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  1.03
       7     1.739  1.15  0.83  0.69  0.91  0.82  0.76  0.72  0.68  0.65  0.63
       3     1.598  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     1.537  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       5     1.482  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     1.392  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       1     1.352  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      76   overall  0.69  0.53  0.44  0.42  0.43  0.42  0.40  0.42  0.44  0.48


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.394  2.10  3.31  2.31  1.82  1.78  1.55  1.59  1.45  1.35  1.27
      14     4.702  3.00  1.97  1.65  1.30  1.27  1.58  1.42  1.30  1.44  1.36
      20     3.548  1.05  0.83  0.77  0.61  0.76  0.89  1.09  1.18  1.10  1.11
      21     2.994  2.00  2.10  1.65  1.44  1.45  1.37  1.23  1.30  1.28  1.21
      30     2.633  1.05  1.29  1.03  1.01  0.93  0.96  0.93  0.91  0.84  1.00
      42     2.402  0.75  0.66  1.01  1.08  1.15  1.11  1.23  1.25  1.20  1.13
      32     2.205  1.64  1.21  1.09  0.85  0.71  0.76  0.81  0.85  0.84  0.84
      33     2.055  1.27  1.00  1.05  1.10  1.08  1.28  1.20  1.10  1.02  0.96
      32     1.934  4.27  2.59  2.05  1.80  1.51  1.32  1.18  1.08  1.05  1.04
      35     1.739  4.50  2.52  1.76  1.56  1.38  1.20  1.14  1.04  0.96  0.90
      32     1.668  3.94  2.07  1.45  1.23  1.03  0.97  0.87  0.79  0.74  0.74
      45     1.598  5.84  3.06  2.14  1.69  1.41  1.23  1.10  1.01  0.94  0.88
      50     1.537  5.04  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      48     1.434  4.16  2.30  1.61  1.26  1.16  1.02  0.91  0.87  0.84  0.79
      35     1.352  5.70  2.99  2.20  1.73  1.45  1.33  1.19  1.09  1.01  0.99
      16     1.315  5.25  2.76  1.93  1.52  1.27  1.25  1.12  1.13  1.05  0.99
     636   overall  2.68  1.62  1.24  1.05  0.94  0.89  0.84  0.81  0.77  0.74


 List of     91 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   33   27   13    1.55   84.79  0.2507E+02  0.1480E+02 "alien"
   31   27   15    1.52   79.66  0.2355E+02  0.1588E+02 "alien"
   51    2    4    1.50   73.72  0.1655E+02  0.9250E+01 "alien"
   51    3    7    1.46   71.39  0.1603E+02  0.9364E+01 "alien"
   35   29   12    1.51   70.17  0.1576E+02  0.1402E+02 "alien"
   32    8   20    1.49   70.15  0.1575E+02  0.1190E+02 "alien"
   34   28   12    1.54   68.57  0.2027E+02  0.1469E+02 "alien"
   26   10   22    1.49   63.91  0.1435E+02  0.1073E+02 "alien"
   39   34   11    1.38   63.78 -0.2087E+02  0.1697E+02 "alien"
   40   30    4    1.53   58.28  0.1723E+02  0.9197E+01 "alien"
   38   31   10    1.46   57.57  0.1293E+02  0.9202E+01 "alien"
   42   32    2    1.46   56.50  0.1269E+02  0.9750E+01 "alien"
   23   13   23    1.46   53.83  0.1209E+02  0.9930E+01 "alien"
   36   29   11    1.51   53.13  0.1571E+02  0.1407E+02 "alien"
   45    3   11    1.54   49.09  0.1451E+02  0.9177E+01 "alien"
   25   10   22    1.50   46.02  0.1033E+02  0.9807E+01 "alien"
   26   23   20    1.46   44.44  0.9977E+01  0.1284E+02 "alien"
   41   31    1    1.51   43.03  0.9662E+01  0.9507E+01 "alien"
   29    9   21    1.49   41.69  0.9360E+01  0.7101E+01 "alien"
   44    4   13    1.51   40.83  0.1207E+02  0.8967E+01 "alien"
   36   30   12    1.47   38.92  0.8739E+01  0.9072E+01 "alien"
   49    3    9    1.48   37.80  0.8487E+01  0.9011E+01 "alien"
   37    9   20    1.40   36.61 -0.1198E+02  0.8189E+01 "alien"
   45    7   16    1.39   34.74 -0.1137E+02  0.1124E+02 "alien"
   33   11   22    1.38   32.19 -0.1054E+02  0.1160E+02 "alien"
   37    4   15    1.62   32.17  0.3120E+02  0.1332E+02 "alien"
   40   33    9    1.41   31.80 -0.1041E+02  0.9197E+01 "alien"
   36   33   14    1.38   31.54 -0.1032E+02  0.1056E+02 "alien"
   25   13   24    1.39   31.17 -0.1020E+02  0.6991E+01 "alien"
   36   10   21    1.38   31.14 -0.1019E+02  0.8233E+01 "alien"
   41   31    4    1.49   30.97  0.6954E+01  0.9458E+01 "alien"
   30   28   19    1.39   30.33 -0.9927E+01  0.8559E+01 "alien"
   28   23   21    1.40   29.36 -0.9608E+01  0.7871E+01 "alien"
   52    4    9    1.40   29.12 -0.9530E+01  0.9718E+01 "alien"
   51    5   11    1.39   28.57 -0.9352E+01  0.9461E+01 "alien"
   29   28   17    1.47   28.30  0.6353E+01  0.8383E+01 "alien"
   35   32   15    1.38   27.55 -0.9015E+01  0.1734E+02 "alien"
   42    4   14    1.53   26.28  0.7769E+01  0.8471E+01 "alien"
   40   30    3    1.54   26.12  0.7721E+01  0.9519E+01 "alien"
   39   31    9    1.46   25.73  0.5778E+01  0.1051E+02 "alien"
   48    6   14    1.39   25.38 -0.8307E+01  0.1314E+02 "alien"
   27   18   23    1.38   25.36 -0.8298E+01  0.8013E+01 "alien"
   39   31    8    1.48   24.68  0.5541E+01  0.9500E+01 "alien"
   54    4    6    1.40   24.02 -0.7863E+01  0.9644E+01 "alien"
   30   26   15    1.56   23.17  0.6849E+01  0.1244E+02 "alien"
   51    2    3    1.51   22.74  0.5105E+01  0.9300E+01 "alien"
   47    2    8    1.56   22.23  0.6573E+01  0.8672E+01 "alien"
   30   25   21    1.35   22.04  0.3729E+02  0.1288E+02 "alien"
   31    7   19    1.56   21.09  0.6235E+01  0.1158E+02 "alien"
   38    9   20    1.39   18.91 -0.6189E+01  0.9137E+01 "alien"
   25   20   21    1.47   18.83  0.4227E+01  0.1113E+02 "alien"
   51    4   10    1.41   18.60 -0.6088E+01  0.1561E+02 "alien"
   26   19   12    1.93   18.42 -0.1854E+02  0.7619E+01 "alien"
   27   26   17    1.53   18.32  0.5415E+01  0.1276E+02 "alien"
   41   32    5    1.46   18.23  0.4092E+01  0.9298E+01 "alien"
   30   28   16    1.49   18.20  0.4087E+01  0.1420E+02 "alien"
   27   25   19    1.46   18.00  0.4041E+01  0.1351E+02 "alien"
   25   23   19    1.52   17.74  0.5245E+01  0.1242E+02 "alien"
   51    6   13    1.35   16.76  0.2836E+02  0.1706E+02 "alien"
   39   28    4    1.59   16.41  0.1592E+02  0.1099E+02 "alien"
   32   27   14    1.54   15.65  0.4626E+01  0.1007E+02 "alien"
   47    1    5    1.61   15.28  0.1482E+02  0.8926E+01 "alien"
   28   21   22    1.38   15.23 -0.4986E+01  0.7888E+01 "alien"
   37    8   19    1.44   14.87  0.3521E+02  0.1468E+02 "alien"
   38    5   16    1.55   14.55  0.4301E+01  0.8235E+01 "alien"
   52    6   12    1.34   14.37  0.2432E+02  0.1006E+02 "alien"
   38    6   17    1.51   14.29  0.3208E+01  0.8960E+01 "alien"
   32   23    9    1.79   14.25 -0.1704E+02  0.9516E+01 "alien"
   41   31    2    1.50   13.79  0.3096E+01  0.9090E+01 "alien"
   33   26   12    1.60   13.32  0.1292E+02  0.9977E+01 "alien"
   46    1    6    1.63   13.21  0.1281E+02  0.1282E+02 "alien"
   30   25   14    1.61   13.00  0.1261E+02  0.7988E+01 "alien"
   45    8   17    1.36   12.71  0.2151E+02  0.1550E+02 "alien"
   26   25   17    1.56   12.66  0.3743E+01  0.1247E+02 "alien"
   27   19   11    1.95   12.26 -0.1234E+02  0.8663E+01 "alien"
   29   23   13    1.71   12.19  0.1967E+02  0.7809E+01 "alien"
   25   25   18    1.54   12.11  0.8284E+01  0.7932E+01 "alien"
   52    2    0    1.49   11.97  0.6157E+01  0.1272E+02 "alien"
   38    4   15    1.59   11.95  0.1159E+02  0.7778E+01 "alien"
   34   22    0    1.91   11.81 -0.2697E+02  0.1237E+02 "alien"
   40   29    2    1.56   11.80  0.3488E+01  0.9316E+01 "alien"
   39   28    0    1.61   11.67  0.2544E+02  0.1292E+02 "alien"
   32   21    6    1.93   11.53 -0.1161E+02  0.8366E+01 "alien"
   22   14   22    1.51   11.37  0.3361E+01  0.9422E+01 "alien"
   26   26   18    1.51   11.30  0.5567E+01  0.7973E+01 "alien"
   36    3   14    1.69   11.19  0.2593E+02  0.1317E+02 "alien"
   53    2    0    1.46   10.84  0.5181E+01  0.1528E+02 "alien"
   36    1    8    1.97   10.77 -0.1084E+02  0.8667E+01 "alien"
   40    3   13    1.62   10.46  0.1014E+02  0.8838E+01 "alien"
   35   24    5    1.78   10.29  0.1661E+02  0.8516E+01 "alien"
   35   23    3    1.83   10.08 -0.1205E+02  0.8541E+01 "alien"


 List of     91 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   52    6   12    1.34   14.37  0.2432E+02  0.1006E+02
   51    6   13    1.35   16.76  0.2836E+02  0.1706E+02
   30   25   21    1.35   22.04  0.3729E+02  0.1288E+02
   45    8   17    1.36   12.71  0.2151E+02  0.1550E+02
   36   33   14    1.38   31.54 -0.1032E+02  0.1056E+02
   27   18   23    1.38   25.36 -0.8298E+01  0.8013E+01
   39   34   11    1.38   63.78 -0.2087E+02  0.1697E+02
   28   21   22    1.38   15.23 -0.4986E+01  0.7888E+01
   36   10   21    1.38   31.14 -0.1019E+02  0.8233E+01
   35   32   15    1.38   27.55 -0.9015E+01  0.1734E+02
   33   11   22    1.38   32.19 -0.1054E+02  0.1160E+02
   30   28   19    1.39   30.33 -0.9927E+01  0.8559E+01
   48    6   14    1.39   25.38 -0.8307E+01  0.1314E+02
   38    9   20    1.39   18.91 -0.6189E+01  0.9137E+01
   51    5   11    1.39   28.57 -0.9352E+01  0.9461E+01
   45    7   16    1.39   34.74 -0.1137E+02  0.1124E+02
   25   13   24    1.39   31.17 -0.1020E+02  0.6991E+01
   54    4    6    1.40   24.02 -0.7863E+01  0.9644E+01
   28   23   21    1.40   29.36 -0.9608E+01  0.7871E+01
   37    9   20    1.40   36.61 -0.1198E+02  0.8189E+01
   52    4    9    1.40   29.12 -0.9530E+01  0.9718E+01
   51    4   10    1.41   18.60 -0.6088E+01  0.1561E+02
   40   33    9    1.41   31.80 -0.1041E+02  0.9197E+01
   37    8   19    1.44   14.87  0.3521E+02  0.1468E+02
   53    2    0    1.46   10.84  0.5181E+01  0.1528E+02
   23   13   23    1.46   53.83  0.1209E+02  0.9930E+01
   51    3    7    1.46   71.39  0.1603E+02  0.9364E+01
   39   31    9    1.46   25.73  0.5778E+01  0.1051E+02
   41   32    5    1.46   18.23  0.4092E+01  0.9298E+01
   38   31   10    1.46   57.57  0.1293E+02  0.9202E+01
   42   32    2    1.46   56.50  0.1269E+02  0.9750E+01
   26   23   20    1.46   44.44  0.9977E+01  0.1284E+02
   27   25   19    1.46   18.00  0.4041E+01  0.1351E+02
   25   20   21    1.47   18.83  0.4227E+01  0.1113E+02
   29   28   17    1.47   28.30  0.6353E+01  0.8383E+01
   36   30   12    1.47   38.92  0.8739E+01  0.9072E+01
   39   31    8    1.48   24.68  0.5541E+01  0.9500E+01
   49    3    9    1.48   37.80  0.8487E+01  0.9011E+01
   26   10   22    1.49   63.91  0.1435E+02  0.1073E+02
   52    2    0    1.49   11.97  0.6157E+01  0.1272E+02
   30   28   16    1.49   18.20  0.4087E+01  0.1420E+02
   41   31    4    1.49   30.97  0.6954E+01  0.9458E+01
   32    8   20    1.49   70.15  0.1575E+02  0.1190E+02
   29    9   21    1.49   41.69  0.9360E+01  0.7101E+01
   51    2    4    1.50   73.72  0.1655E+02  0.9250E+01
   41   31    2    1.50   13.79  0.3096E+01  0.9090E+01
   25   10   22    1.50   46.02  0.1033E+02  0.9807E+01
   41   31    1    1.51   43.03  0.9662E+01  0.9507E+01
   26   26   18    1.51   11.30  0.5567E+01  0.7973E+01
   35   29   12    1.51   70.17  0.1576E+02  0.1402E+02
   51    2    3    1.51   22.74  0.5105E+01  0.9300E+01
   38    6   17    1.51   14.29  0.3208E+01  0.8960E+01
   44    4   13    1.51   40.83  0.1207E+02  0.8967E+01
   36   29   11    1.51   53.13  0.1571E+02  0.1407E+02
   22   14   22    1.51   11.37  0.3361E+01  0.9422E+01
   25   23   19    1.52   17.74  0.5245E+01  0.1242E+02
   31   27   15    1.52   79.66  0.2355E+02  0.1588E+02
   40   30    4    1.53   58.28  0.1723E+02  0.9197E+01
   42    4   14    1.53   26.28  0.7769E+01  0.8471E+01
   27   26   17    1.53   18.32  0.5415E+01  0.1276E+02
   25   25   18    1.54   12.11  0.8284E+01  0.7932E+01
   40   30    3    1.54   26.12  0.7721E+01  0.9519E+01
   32   27   14    1.54   15.65  0.4626E+01  0.1007E+02
   45    3   11    1.54   49.09  0.1451E+02  0.9177E+01
   34   28   12    1.54   68.57  0.2027E+02  0.1469E+02
   33   27   13    1.55   84.79  0.2507E+02  0.1480E+02
   38    5   16    1.55   14.55  0.4301E+01  0.8235E+01
   30   26   15    1.56   23.17  0.6849E+01  0.1244E+02
   47    2    8    1.56   22.23  0.6573E+01  0.8672E+01
   40   29    2    1.56   11.80  0.3488E+01  0.9316E+01
   26   25   17    1.56   12.66  0.3743E+01  0.1247E+02
   31    7   19    1.56   21.09  0.6235E+01  0.1158E+02
   39   28    4    1.59   16.41  0.1592E+02  0.1099E+02
   38    4   15    1.59   11.95  0.1159E+02  0.7778E+01
   33   26   12    1.60   13.32  0.1292E+02  0.9977E+01
   47    1    5    1.61   15.28  0.1482E+02  0.8926E+01
   30   25   14    1.61   13.00  0.1261E+02  0.7988E+01
   39   28    0    1.61   11.67  0.2544E+02  0.1292E+02
   40    3   13    1.62   10.46  0.1014E+02  0.8838E+01
   37    4   15    1.62   32.17  0.3120E+02  0.1332E+02
   46    1    6    1.63   13.21  0.1281E+02  0.1282E+02
   36    3   14    1.69   11.19  0.2593E+02  0.1317E+02
   29   23   13    1.71   12.19  0.1967E+02  0.7809E+01
   35   24    5    1.78   10.29  0.1661E+02  0.8516E+01
   32   23    9    1.79   14.25 -0.1704E+02  0.9516E+01
   35   23    3    1.83   10.08 -0.1205E+02  0.8541E+01
   34   22    0    1.91   11.81 -0.2697E+02  0.1237E+02
   26   19   12    1.93   18.42 -0.1854E+02  0.7619E+01
   32   21    6    1.93   11.53 -0.1161E+02  0.8366E+01
   27   19   11    1.95   12.26 -0.1234E+02  0.8663E+01
   36    1    8    1.97   10.77 -0.1084E+02  0.8667E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        49
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1275
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1155
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      120


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.3 sec
