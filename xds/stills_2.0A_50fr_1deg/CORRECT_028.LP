
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.517    77.517    38.791  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      28      28
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      28
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000914  0.000686
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
 REFLECTING_RANGE_E.S.D.=   0.134
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

    1329 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      28 ...      28
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1329 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               42
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

 *  44        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  31        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  35        mP          0.0      77.5   38.8   77.5  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  33        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.6  109.6   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.5   77.5   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.8   38.8   77.5  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.8   38.8   77.5  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.8  159.8   77.5  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.8  159.8   77.5  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.8  159.8   77.5  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.8  159.8   77.5  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    27        mC        500.0     159.8   38.8  109.6  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.8  109.6  116.3  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.8  159.8  159.8  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.6  116.3   38.8  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        998.4     173.3   77.5   38.8  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.4      77.5  173.3   38.8  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.4      77.5  173.3   38.8  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.4      77.5   77.5   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.3  116.3  116.3  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.7   86.7  109.6  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.6   86.7   86.7  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.7   86.7  109.6  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.7   86.7  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.8   86.7  248.2 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.7   86.7   77.5  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.8   77.5  177.6  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.7   86.7  177.6 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.7   86.7   86.7  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.6  109.6   38.8  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.8   77.5  177.6 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.7  177.6   77.5 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    4.2      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.867   -23.509     0.071
 UNIT_CELL_B-AXIS=    23.509    73.867     0.020
 UNIT_CELL_C-AXIS=    -0.037     0.001    38.790
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       189 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.24
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.454    77.454    38.766  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  7.4E-02 7.4E-02 3.8E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012911  0.012911  0.025796  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.823   -23.438     0.041
 COORDINATES OF UNIT CELL B-AXIS    23.438    73.823     0.008
 COORDINATES OF UNIT CELL C-AXIS    -0.021     0.003    38.766
 CRYSTAL MOSAICITY (DEGREES)     0.134
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999975  0.006780  0.001900
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.88   1232.91
 DETECTOR ORIGIN (PIXELS) AT                     1261.12   1242.77
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.18
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.45    77.45    38.77  90.000  90.000  90.000


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
  0.05  0.15        60        109.087       1.182
  0.15  0.25       166        135.225       0.953
  0.25  0.35       186        116.866       1.103
  0.35  0.45       154         79.928       1.613
  0.45  0.55       173        218.829       0.589
  0.55  0.65       183         82.213       1.568
  0.65  0.75       201        152.887       0.843
  0.75  0.85       143        128.265       1.005
  0.85  0.95        63        100.459       1.283
 
 -0.05  0.95      1329        128.926       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        8
 DEGREES OF FREEDOM OF CHI^2 FIT               207.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.887
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    27.1 XMAX=    27.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59443 NYBIN=    8
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        444


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               207.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.884
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    11.0 XMAX=  2453.7 NXBIN=    3
 YMIN=     8.3 YMAX=  2520.0 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        444


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               207.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.873
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    27.1 XMAX=    27.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1647    1691
 DETECTOR_SURFACE_POSITION=     817    1691
 DETECTOR_SURFACE_POSITION=     817     838
 DETECTOR_SURFACE_POSITION=    1647     838
 DETECTOR_SURFACE_POSITION=    2171    1664
 DETECTOR_SURFACE_POSITION=    1621    2230
 DETECTOR_SURFACE_POSITION=     843    2230
 DETECTOR_SURFACE_POSITION=     293    1664
 DETECTOR_SURFACE_POSITION=     293     864
 DETECTOR_SURFACE_POSITION=     843     299
 DETECTOR_SURFACE_POSITION=    1621     299
 DETECTOR_SURFACE_POSITION=    2171     864
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        444


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
 8.676E-01  1.005E-03   33.87


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

   24.493   3.155    24.12   1.35      2.10      1.88      66      90      24
    3.155   2.240     8.44   1.10      6.08      5.41     182     210       8
    2.240   1.832     1.11   0.92     42.56     43.95     200     228       0
    1.832   1.587     0.09   1.05    479.28    447.37     240     273       2
    1.587   1.420    -0.04   0.90    -99.90    -99.90     260     315       4
    1.420   1.297    -0.15   0.71    -99.90    -99.90     126     172       2
   --------------------------------------------------------------------------
   24.493   1.297     3.25   0.98      7.89      7.77    1074    1288      40
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      28      28)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.684  -0.1492E+02  0.1201E+02        -1.24           1 
   47    0    0     1.648   0.3811E+01  0.1198E+02         0.32           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     2.9% OF MEAN INTENSITY


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

     3.85          55      34      1285        2.6%       2.0%      1.8%       42   32.69      2.8%    99.9*     0    0.000       0
     2.74          88      52      2098        2.5%       4.3%      3.1%       72   20.62      6.0%    99.4*     0    0.000       0
     2.24         157      90      2665        3.4%       6.6%      6.9%      134    9.63      9.3%    99.0*     0    0.000       0
     1.94         149      84      3101        2.7%      30.3%     30.4%      130    2.16     42.9%    90.0*     0    0.000       0
     1.74         156      86      3481        2.5%     347.0%    324.4%      140    0.15    490.7%   -36.2      0    0.000       0
     1.59         194     110      3847        2.9%     517.8%    510.0%      168    0.08    732.3%    -5.7      0    0.000       0
     1.47         203     120      4142        2.9%     -99.9%    -99.9%      166    0.00    -99.9%    16.6      0    0.000       0
     1.38         197     115      4433        2.6%     -99.9%    -99.9%      164    0.00    -99.9%   -11.1      0    0.000       0
     1.30          82      56      4737        1.2%    3132.4%   4138.9%       52    0.10   4429.8%    14.4      0    0.000       0
    total        1281     747     29789        2.5%       7.8%      7.7%     1068    4.34     11.1%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1329
 NUMBER OF REJECTED MISFITS                                40
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1288
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    751



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

   WILSON LINE (using all data) : A=   7.242 B=  35.705 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      14     7.396   0.005  1.9413E+03       7.571     -36.0
      16     4.702   0.011  2.0998E+03       7.650     -18.0
      24     3.516   0.020  7.7487E+02       6.653      14.6
      25     2.998   0.028  2.2428E+02       5.413      32.9
      36     2.640   0.036  2.0839E+02       5.339      26.5
      48     2.404   0.043  9.9180E+01       4.597      30.6
      36     2.209   0.051  4.4856E+01       3.803      33.6
      43     2.058   0.059  1.6378E+01       2.796      37.7
      38     1.935   0.067  5.6669E-01      -0.568      58.5
      47     1.830   0.075  1.5305E+00       0.426      45.6
      43     1.740   0.083  1.5497E+00       0.438      41.2
      49     1.669   0.090  1.7960E+00       0.586      37.1
      48     1.599   0.098 -8.5873E-01    -999.999   -1000.0
      59     1.538   0.106  1.3417E-01      -2.009      43.8
      70     1.482   0.114 -1.3482E-01    -999.999   -1000.0
      52     1.434   0.122  1.2100E-01      -2.112      38.5
      44     1.392   0.129 -1.8374E+00    -999.999   -1000.0
      38     1.353   0.137  4.5323E-02      -3.094      37.8
      17     1.316   0.144  4.8965E+00       1.589      19.6


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.396        1.911        1.376        0.775
       3     4.702        1.287        0.797        0.367
       2     3.516        0.400        0.088        0.014
       3     2.998        2.177        1.200        0.491
       2     2.640        0.399        0.123        0.027
       5     2.404        0.891        0.469        0.194
       1     2.209        0.747        0.223        0.048
       6     2.058        0.387        0.191        0.070
      10     1.830        8.468       12.545       20.683
       7     1.740       10.899       11.546       30.345
       3     1.538      958.614    12235.515   113099.917
       6     1.482       26.016       72.292      154.445
       3     1.434     1531.943    25962.718   319534.566
      76   overall      175.588     -202.343    47138.494


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.396        0.537        0.392        0.225
      13     4.702        1.220        1.465        1.393
      22     3.516        0.910        0.828        0.707
      22     2.998        0.727        0.572        0.372
      34     2.640        1.215        1.526        1.799
      43     2.404        0.963        0.853        0.638
      35     2.209        0.946        0.800        0.590
      37     2.058        1.693        2.461        3.321
      32     1.935      119.552     -121.444    10756.158
      37     1.830        9.652       10.633       64.343
      36     1.740        8.994        2.324       42.283
      38     1.669       11.001       18.735       56.332
      37     1.353    21664.025   736274.618176900116.901
      17     1.316        4.819       14.824       42.315
     671   overall     1641.029    44893.812 10657186.802

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.396  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.38  0.73
       3     4.702  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.516  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.998  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.640  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.404  0.00  0.00  0.00  0.42  0.38  0.71  0.67  0.64  0.91  0.88
       1     2.209  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.058  2.01  1.45  1.20  1.76  1.60  1.48  1.40  1.33  1.27  1.22
      10     1.830  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       7     1.740  1.73  1.24  1.03  0.91  0.82  0.76  0.72  0.68  0.65  0.63
       3     1.538  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     1.482  3.36  2.41  2.00  1.76  1.60  1.48  1.40  1.33  1.27  1.22
       3     1.434  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      76   overall  0.95  0.69  0.57  0.58  0.53  0.52  0.48  0.48  0.48  0.50


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.396  3.15  3.31  2.31  1.82  1.78  1.55  1.59  1.45  1.35  1.27
      13     4.702  1.62  2.12  1.48  1.17  1.37  1.53  1.38  1.26  1.43  1.34
      22     3.516  0.96  0.75  0.70  0.55  0.69  0.81  1.08  1.16  1.07  1.08
      22     2.998  1.43  2.26  1.58  1.52  1.39  1.31  1.35  1.24  1.23  1.15
      34     2.640  0.62  1.30  1.36  1.16  1.05  0.91  0.99  1.07  1.09  1.12
      43     2.404  1.22  0.90  1.17  1.13  1.12  1.13  1.20  1.10  1.10  1.07
      35     2.209  1.80  1.42  0.99  0.87  0.87  0.82  0.79  0.83  0.96  0.99
      37     2.058  2.84  1.79  1.36  1.23  1.24  1.14  1.02  0.98  0.91  0.98
      32     1.935  4.93  2.59  1.81  1.42  1.19  1.04  0.93  0.85  0.79  0.74
      37     1.830  4.54  2.39  1.77  1.39  1.17  1.02  0.91  0.83  0.77  0.73
      36     1.740  3.79  1.99  1.39  1.10  0.92  0.86  0.77  0.76  0.75  0.75
      38     1.669  3.87  2.03  1.52  1.28  1.07  0.93  0.84  0.81  0.75  0.71
      37     1.353  5.96  3.13  2.19  1.72  1.44  1.26  1.13  1.03  0.96  0.90
      17     1.316  3.09  1.62  1.13  1.07  1.20  1.17  1.05  1.18  1.29  1.21
     671   overall  1.83  1.18  0.90  0.76  0.70  0.65  0.63  0.61  0.60  0.58


 List of    123 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   53    6   11    1.34  446.86  0.2025E+02  0.1036E+02 "alien"
   52    5   11    1.37  428.30  0.1941E+02  0.1706E+02 "alien"
   45    8   17    1.36  362.02  0.1641E+02  0.1487E+02 "alien"
   41   10   20    1.33  282.92  0.1282E+02  0.9156E+01 "alien"
   37    8   19    1.44  282.33  0.3416E+02  0.1435E+02 "alien"
   32   30   19    1.34  274.03  0.1242E+02  0.1434E+02 "alien"
   37   10   21    1.36  261.02  0.1183E+02  0.1293E+02 "alien"
   50    7   14    1.34  256.75  0.1164E+02  0.1071E+02 "alien"
   38   34   13    1.35  234.78  0.1064E+02  0.1003E+02 "alien"
   26   25   17    1.56  213.51  0.2865E+02  0.1294E+02 "alien"
   39   31    9    1.46  192.90 -0.2601E+02  0.1046E+02 "alien"
   27   13   24    1.37  177.65  0.8052E+01  0.1028E+02 "alien"
   46    8   17    1.34  156.73  0.7103E+01  0.9312E+01 "alien"
   49    4   11    1.44  150.49  0.1821E+02  0.9597E+01 "alien"
   35   33   15    1.37  139.86  0.6339E+01  0.1505E+02 "alien"
   34   29   13    1.50  120.85 -0.1629E+02  0.9399E+01 "alien"
   26   24   19    1.49  118.57 -0.1599E+02  0.1310E+02 "alien"
   45    6   15    1.42  117.76  0.1425E+02  0.9548E+01 "alien"
   52    2    3    1.48  113.10 -0.1525E+02  0.9189E+01 "alien"
   40   30    5    1.52  110.61  0.1484E+02  0.9385E+01 "alien"
   52    4    8    1.42  109.18  0.1321E+02  0.1036E+02 "alien"
   39   29    5    1.56  106.09  0.1423E+02  0.9159E+01 "alien"
   33   30   15    1.44  103.62  0.1254E+02  0.9017E+01 "alien"
   51    2    5    1.49  103.58 -0.1397E+02  0.9306E+01 "alien"
   37   31   11    1.46  102.92 -0.1387E+02  0.1559E+02 "alien"
   44    3   12    1.54   95.94  0.1287E+02  0.9054E+01 "alien"
   42    7   17    1.42   94.15  0.1139E+02  0.9672E+01 "alien"
   33   29   14    1.49   92.73 -0.1250E+02  0.1048E+02 "alien"
   39   29    6    1.55   92.28  0.1238E+02  0.9815E+01 "alien"
   37   34   14    1.35   91.53  0.4149E+01  0.9799E+01 "alien"
   45    3   11    1.54   89.92  0.1206E+02  0.9130E+01 "alien"
   29    9   21    1.50   89.13 -0.1202E+02  0.6734E+01 "alien"
   38   32   11    1.43   88.17  0.1067E+02  0.9725E+01 "alien"
   43    5   15    1.47   87.07 -0.1174E+02  0.9459E+01 "alien"
   39   31    8    1.48   86.22 -0.1162E+02  0.9388E+01 "alien"
   42    9   19    1.35   84.15  0.3814E+01  0.8893E+01 "alien"
   37   30   11    1.48   78.04 -0.1052E+02  0.1089E+02 "alien"
   38    7   18    1.47   77.80 -0.1049E+02  0.8555E+01 "alien"
   46    4   12    1.49   74.30 -0.1002E+02  0.9151E+01 "alien"
   31    7   19    1.56   71.57  0.9602E+01  0.1170E+02 "alien"
   40   30    3    1.54   71.05  0.9533E+01  0.9591E+01 "alien"
   37   30   10    1.50   70.75 -0.9539E+01  0.1085E+02 "alien"
   42   32    3    1.46   68.26  0.8259E+01  0.1181E+02 "alien"
   25   20   21    1.47   63.81 -0.8603E+01  0.1074E+02 "alien"
   41   31    2    1.50   62.20 -0.8385E+01  0.9190E+01 "alien"
   38   29    7    1.56   61.44  0.8243E+01  0.9628E+01 "alien"
   31    8   20    1.51   59.19  0.7942E+01  0.7191E+01 "alien"
   51    3    8    1.45   58.19  0.7041E+01  0.1099E+02 "alien"
   47    2    9    1.54   57.33  0.7692E+01  0.1617E+02 "alien"
   47    2    8    1.56   56.23  0.7545E+01  0.8699E+01 "alien"
   40   29    1    1.57   55.69  0.7471E+01  0.9246E+01 "alien"
   42    4   14    1.53   52.43  0.7035E+01  0.8441E+01 "alien"
   25   23   19    1.52   50.91  0.6830E+01  0.1252E+02 "alien"
   27    9   21    1.53   45.87  0.6155E+01  0.9382E+01 "alien"
   34   11   22    1.37   43.83  0.1986E+01  0.7906E+01 "alien"
   40   29    2    1.56   43.61  0.5851E+01  0.9309E+01 "alien"
   38   31   10    1.46   42.57 -0.5739E+01  0.9153E+01 "alien"
   38    5   16    1.55   40.63  0.5451E+01  0.8179E+01 "alien"
   30    9   21    1.48   40.04 -0.5399E+01  0.7378E+01 "alien"
   41   32    4    1.47   38.01 -0.5124E+01  0.1154E+02 "alien"
   38    8   19    1.43   37.12  0.4492E+01  0.8213E+01 "alien"
   37   29    9    1.54   36.84  0.4942E+01  0.9171E+01 "alien"
   48    5   13    1.41   34.89  0.4222E+01  0.1001E+02 "alien"
   37   29   10    1.52   34.03  0.4566E+01  0.1508E+02 "alien"
   28   20   11    1.90   33.82  0.1917E+02  0.8169E+01 "alien"
   36   30   12    1.47   33.44 -0.4509E+01  0.9038E+01 "alien"
   42    6   16    1.46   33.04  0.3998E+01  0.9404E+01 "alien"
   53    3    4    1.44   31.88  0.3858E+01  0.9703E+01 "alien"
   42    5   15    1.49   31.19 -0.4205E+01  0.8479E+01 "alien"
   35   30   13    1.46   30.74 -0.4145E+01  0.9215E+01 "alien"
   26    9   21    1.54   29.81  0.4000E+01  0.1037E+02 "alien"
   53    3    6    1.42   29.76  0.3601E+01  0.1118E+02 "alien"
   34    6   18    1.55   28.70  0.3851E+01  0.7688E+01 "alien"
   40   31    6    1.49   28.39 -0.3828E+01  0.9235E+01 "alien"
   36    1    8    1.97   28.36  0.1607E+02  0.8933E+01 "alien"
   39    2    2    1.97   27.73  0.1571E+02  0.9161E+01 "alien"
   36   34   15    1.34   27.68  0.1255E+01  0.1356E+02 "alien"
   31    9   21    1.46   26.74 -0.3606E+01  0.1167E+02 "alien"
   30   28   16    1.49   24.08 -0.3246E+01  0.1455E+02 "alien"
   51    6   13    1.35   23.78  0.1078E+01  0.1626E+02 "alien"
   48    2    8    1.53   23.03  0.3090E+01  0.9035E+01 "alien"
   53    2    1    1.46   22.82 -0.3076E+01  0.1150E+02 "alien"
   28   28   18    1.45   20.82  0.1123E+02  0.9929E+01 "alien"
   33   26   12    1.60   20.75 -0.1782E+02  0.9511E+01 "alien"
   21   10   22    1.56   20.13  0.2701E+01  0.8540E+01 "alien"
   34    7   19    1.51   19.98 -0.2694E+01  0.1217E+02 "alien"
   46    3   11    1.52   19.97  0.2679E+01  0.9091E+01 "alien"
   37   27    7    1.62   19.95 -0.1713E+02  0.8648E+01 "alien"
   21   17   14    1.99   19.56  0.1108E+02  0.1112E+02 "alien"
   50    2    6    1.51   19.32 -0.2604E+01  0.8930E+01 "alien"
   40   30    4    1.53   18.94  0.2542E+01  0.8968E+01 "alien"
   29   29   18    1.42   18.85  0.8671E+01  0.1513E+02 "alien"
   47    1    6    1.60   18.43 -0.1583E+02  0.1233E+02 "alien"
   38   28    5    1.61   18.36 -0.1577E+02  0.9937E+01 "alien"
   25   25   18    1.54   18.13  0.8732E+01  0.7981E+01 "alien"
   24   24   18    1.57   18.10  0.8696E+01  0.7955E+01 "alien"
   41   32    5    1.46   18.02 -0.2429E+01  0.9152E+01 "alien"
   42   32    1    1.47   16.65 -0.2245E+01  0.9465E+01 "alien"
   25   10   22    1.50   15.89 -0.2142E+01  0.9022E+01 "alien"
   39    6   17    1.49   15.85 -0.2137E+01  0.8085E+01 "alien"
   37    1    7    1.96   15.35  0.8699E+01  0.8445E+01 "alien"
   39   29    4    1.57   14.80 -0.1271E+02  0.1042E+02 "alien"
   44    5   14    1.48   14.74 -0.1987E+01  0.1562E+02 "alien"
   39   33   10    1.41   14.60 -0.2683E+02  0.1773E+02 "alien"
   45    2   10    1.57   14.15 -0.1215E+02  0.9589E+01 "alien"
   40    1    5    1.88   14.05  0.2150E+02  0.1204E+02 "alien"
   31   27   15    1.52   13.83  0.1855E+01  0.1434E+02 "alien"
   45    5   14    1.46   13.32  0.1612E+01  0.8685E+01 "alien"
   33    6   18    1.58   12.97 -0.1114E+02  0.7040E+01 "alien"
   54    3    3    1.42   12.86  0.1556E+01  0.9326E+01 "alien"
   30   20    8    1.96   12.45  0.7057E+01  0.8383E+01 "alien"
   52    2    2    1.48   12.23 -0.1648E+01  0.9213E+01 "alien"
   34   26   10    1.64   12.19  0.2189E+02  0.9177E+01 "alien"
   32   21    6    1.93   11.82  0.6701E+01  0.8544E+01 "alien"
   52    2    1    1.49   11.64 -0.1569E+01  0.9067E+01 "alien"
   53    7   13    1.30   11.42  0.5590E+02  0.1476E+02 "alien"
   25   18   12    1.98   11.01  0.6238E+01  0.1256E+02 "alien"
   39   28    0    1.61   10.63 -0.1917E+02  0.1218E+02 "alien"
   34   22    1    1.91   10.57  0.5992E+01  0.9034E+01 "alien"
   36   25    4    1.74   10.36  0.1606E+02  0.9519E+01 "alien"
   29    1   13    1.99   10.32  0.5850E+01  0.8019E+01 "alien"
   32   30   16    1.43   10.05  0.1216E+01  0.8753E+01 "alien"
   28   24   15    1.63   10.00 -0.8589E+01  0.7775E+01 "alien"


 List of    123 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   11.42  0.5590E+02  0.1476E+02
   41   10   20    1.33  282.92  0.1282E+02  0.9156E+01
   32   30   19    1.34  274.03  0.1242E+02  0.1434E+02
   36   34   15    1.34   27.68  0.1255E+01  0.1356E+02
   46    8   17    1.34  156.73  0.7103E+01  0.9312E+01
   50    7   14    1.34  256.75  0.1164E+02  0.1071E+02
   53    6   11    1.34  446.86  0.2025E+02  0.1036E+02
   51    6   13    1.35   23.78  0.1078E+01  0.1626E+02
   37   34   14    1.35   91.53  0.4149E+01  0.9799E+01
   42    9   19    1.35   84.15  0.3814E+01  0.8893E+01
   38   34   13    1.35  234.78  0.1064E+02  0.1003E+02
   45    8   17    1.36  362.02  0.1641E+02  0.1487E+02
   37   10   21    1.36  261.02  0.1183E+02  0.1293E+02
   35   33   15    1.37  139.86  0.6339E+01  0.1505E+02
   52    5   11    1.37  428.30  0.1941E+02  0.1706E+02
   34   11   22    1.37   43.83  0.1986E+01  0.7906E+01
   27   13   24    1.37  177.65  0.8052E+01  0.1028E+02
   39   33   10    1.41   14.60 -0.2683E+02  0.1773E+02
   48    5   13    1.41   34.89  0.4222E+01  0.1001E+02
   52    4    8    1.42  109.18  0.1321E+02  0.1036E+02
   29   29   18    1.42   18.85  0.8671E+01  0.1513E+02
   42    7   17    1.42   94.15  0.1139E+02  0.9672E+01
   53    3    6    1.42   29.76  0.3601E+01  0.1118E+02
   54    3    3    1.42   12.86  0.1556E+01  0.9326E+01
   45    6   15    1.42  117.76  0.1425E+02  0.9548E+01
   38   32   11    1.43   88.17  0.1067E+02  0.9725E+01
   38    8   19    1.43   37.12  0.4492E+01  0.8213E+01
   32   30   16    1.43   10.05  0.1216E+01  0.8753E+01
   49    4   11    1.44  150.49  0.1821E+02  0.9597E+01
   33   30   15    1.44  103.62  0.1254E+02  0.9017E+01
   53    3    4    1.44   31.88  0.3858E+01  0.9703E+01
   37    8   19    1.44  282.33  0.3416E+02  0.1435E+02
   51    3    8    1.45   58.19  0.7041E+01  0.1099E+02
   28   28   18    1.45   20.82  0.1123E+02  0.9929E+01
   45    5   14    1.46   13.32  0.1612E+01  0.8685E+01
   42   32    3    1.46   68.26  0.8259E+01  0.1181E+02
   42    6   16    1.46   33.04  0.3998E+01  0.9404E+01
   53    2    1    1.46   22.82 -0.3076E+01  0.1150E+02
   37   31   11    1.46  102.92 -0.1387E+02  0.1559E+02
   39   31    9    1.46  192.90 -0.2601E+02  0.1046E+02
   41   32    5    1.46   18.02 -0.2429E+01  0.9152E+01
   38   31   10    1.46   42.57 -0.5739E+01  0.9153E+01
   31    9   21    1.46   26.74 -0.3606E+01  0.1167E+02
   35   30   13    1.46   30.74 -0.4145E+01  0.9215E+01
   42   32    1    1.47   16.65 -0.2245E+01  0.9465E+01
   38    7   18    1.47   77.80 -0.1049E+02  0.8555E+01
   25   20   21    1.47   63.81 -0.8603E+01  0.1074E+02
   43    5   15    1.47   87.07 -0.1174E+02  0.9459E+01
   36   30   12    1.47   33.44 -0.4509E+01  0.9038E+01
   41   32    4    1.47   38.01 -0.5124E+01  0.1154E+02
   37   30   11    1.48   78.04 -0.1052E+02  0.1089E+02
   52    2    3    1.48  113.10 -0.1525E+02  0.9189E+01
   44    5   14    1.48   14.74 -0.1987E+01  0.1562E+02
   30    9   21    1.48   40.04 -0.5399E+01  0.7378E+01
   39   31    8    1.48   86.22 -0.1162E+02  0.9388E+01
   52    2    2    1.48   12.23 -0.1648E+01  0.9213E+01
   33   29   14    1.49   92.73 -0.1250E+02  0.1048E+02
   52    2    1    1.49   11.64 -0.1569E+01  0.9067E+01
   39    6   17    1.49   15.85 -0.2137E+01  0.8085E+01
   46    4   12    1.49   74.30 -0.1002E+02  0.9151E+01
   30   28   16    1.49   24.08 -0.3246E+01  0.1455E+02
   51    2    5    1.49  103.58 -0.1397E+02  0.9306E+01
   40   31    6    1.49   28.39 -0.3828E+01  0.9235E+01
   26   24   19    1.49  118.57 -0.1599E+02  0.1310E+02
   42    5   15    1.49   31.19 -0.4205E+01  0.8479E+01
   29    9   21    1.50   89.13 -0.1202E+02  0.6734E+01
   34   29   13    1.50  120.85 -0.1629E+02  0.9399E+01
   37   30   10    1.50   70.75 -0.9539E+01  0.1085E+02
   41   31    2    1.50   62.20 -0.8385E+01  0.9190E+01
   25   10   22    1.50   15.89 -0.2142E+01  0.9022E+01
   50    2    6    1.51   19.32 -0.2604E+01  0.8930E+01
   34    7   19    1.51   19.98 -0.2694E+01  0.1217E+02
   31    8   20    1.51   59.19  0.7942E+01  0.7191E+01
   37   29   10    1.52   34.03  0.4566E+01  0.1508E+02
   46    3   11    1.52   19.97  0.2679E+01  0.9091E+01
   40   30    5    1.52  110.61  0.1484E+02  0.9385E+01
   25   23   19    1.52   50.91  0.6830E+01  0.1252E+02
   31   27   15    1.52   13.83  0.1855E+01  0.1434E+02
   27    9   21    1.53   45.87  0.6155E+01  0.9382E+01
   40   30    4    1.53   18.94  0.2542E+01  0.8968E+01
   48    2    8    1.53   23.03  0.3090E+01  0.9035E+01
   42    4   14    1.53   52.43  0.7035E+01  0.8441E+01
   25   25   18    1.54   18.13  0.8732E+01  0.7981E+01
   47    2    9    1.54   57.33  0.7692E+01  0.1617E+02
   40   30    3    1.54   71.05  0.9533E+01  0.9591E+01
   37   29    9    1.54   36.84  0.4942E+01  0.9171E+01
   44    3   12    1.54   95.94  0.1287E+02  0.9054E+01
   26    9   21    1.54   29.81  0.4000E+01  0.1037E+02
   45    3   11    1.54   89.92  0.1206E+02  0.9130E+01
   39   29    6    1.55   92.28  0.1238E+02  0.9815E+01
   38    5   16    1.55   40.63  0.5451E+01  0.8179E+01
   34    6   18    1.55   28.70  0.3851E+01  0.7688E+01
   38   29    7    1.56   61.44  0.8243E+01  0.9628E+01
   21   10   22    1.56   20.13  0.2701E+01  0.8540E+01
   47    2    8    1.56   56.23  0.7545E+01  0.8699E+01
   39   29    5    1.56  106.09  0.1423E+02  0.9159E+01
   40   29    2    1.56   43.61  0.5851E+01  0.9309E+01
   26   25   17    1.56  213.51  0.2865E+02  0.1294E+02
   31    7   19    1.56   71.57  0.9602E+01  0.1170E+02
   24   24   18    1.57   18.10  0.8696E+01  0.7955E+01
   40   29    1    1.57   55.69  0.7471E+01  0.9246E+01
   45    2   10    1.57   14.15 -0.1215E+02  0.9589E+01
   39   29    4    1.57   14.80 -0.1271E+02  0.1042E+02
   33    6   18    1.58   12.97 -0.1114E+02  0.7040E+01
   47    1    6    1.60   18.43 -0.1583E+02  0.1233E+02
   33   26   12    1.60   20.75 -0.1782E+02  0.9511E+01
   38   28    5    1.61   18.36 -0.1577E+02  0.9937E+01
   39   28    0    1.61   10.63 -0.1917E+02  0.1218E+02
   37   27    7    1.62   19.95 -0.1713E+02  0.8648E+01
   28   24   15    1.63   10.00 -0.8589E+01  0.7775E+01
   34   26   10    1.64   12.19  0.2189E+02  0.9177E+01
   36   25    4    1.74   10.36  0.1606E+02  0.9519E+01
   40    1    5    1.88   14.05  0.2150E+02  0.1204E+02
   28   20   11    1.90   33.82  0.1917E+02  0.8169E+01
   34   22    1    1.91   10.57  0.5992E+01  0.9034E+01
   32   21    6    1.93   11.82  0.6701E+01  0.8544E+01
   37    1    7    1.96   15.35  0.8699E+01  0.8445E+01
   30   20    8    1.96   12.45  0.7057E+01  0.8383E+01
   36    1    8    1.97   28.36  0.1607E+02  0.8933E+01
   39    2    2    1.97   27.73  0.1571E+02  0.9161E+01
   25   18   12    1.98   11.01  0.6238E+01  0.1256E+02
   29    1   13    1.99   10.32  0.5850E+01  0.8019E+01
   21   17   14    1.99   19.56  0.1108E+02  0.1112E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        94
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1329
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1128
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      201


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
