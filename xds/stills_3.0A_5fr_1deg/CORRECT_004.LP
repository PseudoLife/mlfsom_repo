
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.335    77.335    38.714  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_3.0A_5fr_1deg/stills_3.0A_5fr_1deg_???_001.img SMV
 DATA_RANGE=       4       4
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       4
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001017  0.000626
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
 REFLECTING_RANGE_E.S.D.=   0.168
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

    1211 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       4 ...       4
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1211 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               37
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
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0  1  0  0  1  0  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  75.9    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  75.9    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    27        mC        500.0     159.4   38.7  109.4  90.0 133.3  75.9    0  2  1  0  0  0  1  0  1 -1  0  0
    19        oI        500.0      38.7  109.4  116.0  90.0  70.5  90.0    0  0  1  0  1 -1  0  0  1  1  1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.1 104.1    0  0 -1  0  0  2  1  0  2  0  1  0
    18        tI        625.0     109.4  116.0   38.7  70.5  90.0  90.0    1 -1  0  0 -1 -1 -1  0  0  0 -1  0
    30        mC        997.6      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0 -2  1  0  0  0  0  1  0
    41        mC        997.6     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    40        oC        997.6      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        997.6      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0    1 -1 -1  0 -1  1 -1  0  1  1  1  0
     2        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.4   86.5   86.5  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.7   86.5  247.6 102.1  99.0 116.6    0  0 -1  0  0  1  1  0  3 -1  1  0
    10        mC        999.0      86.5   86.5   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.2  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.5   86.5  177.2 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.2  109.4   38.7  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.7   77.3  177.2 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.5  177.2   77.3 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    4.6      13    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.569   -23.837    -0.170
 UNIT_CELL_B-AXIS=    23.837    73.569    -0.053
 UNIT_CELL_C-AXIS=     0.089    -0.001    38.714
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       104 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.16
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.378    77.378    38.722  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.5E-01 2.5E-01 1.4E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012924  0.012924  0.025825  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.751   -23.413     0.014
 COORDINATES OF UNIT CELL B-AXIS    23.413    73.751     0.005
 COORDINATES OF UNIT CELL C-AXIS    -0.007     0.000    38.722
 CRYSTAL MOSAICITY (DEGREES)     0.168
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999976  0.006560  0.002220
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.85   1232.90
 DETECTOR ORIGIN (PIXELS) AT                     1260.62   1242.45
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.47
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.38    77.38    38.72  90.000  90.000  90.000


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
  0.05  0.15        46         59.652       1.782
  0.15  0.25       150         85.020       1.250
  0.25  0.35       162        107.163       0.992
  0.35  0.45       152         49.274       2.157
  0.45  0.55       160        219.133       0.485
  0.55  0.65       190         70.010       1.518
  0.65  0.75       175        144.383       0.736
  0.75  0.85       152         87.932       1.209
  0.85  0.95        24         57.415       1.851
 
 -0.05  0.95      1211        106.293       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT                88.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.693
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     3.1 XMAX=     3.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59568 NYBIN=    4
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        202


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT                88.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.687
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    33.1 XMAX=  2439.4 NXBIN=    2
 YMIN=     9.1 YMAX=  2519.9 NYBIN=    2
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        202


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT                81.5
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.702
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     3.1 XMAX=     3.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1236    1265
 DETECTOR_SURFACE_POSITION=    1645    1691
 DETECTOR_SURFACE_POSITION=     828    1691
 DETECTOR_SURFACE_POSITION=     828     838
 DETECTOR_SURFACE_POSITION=    1645     838
 DETECTOR_SURFACE_POSITION=    2161    1664
 DETECTOR_SURFACE_POSITION=    1619    2230
 DETECTOR_SURFACE_POSITION=     853    2230
 DETECTOR_SURFACE_POSITION=     311    1664
 DETECTOR_SURFACE_POSITION=     311     865
 DETECTOR_SURFACE_POSITION=     853     299
 DETECTOR_SURFACE_POSITION=    1619     299
 DETECTOR_SURFACE_POSITION=    2161     865
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        202


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
 6.873E-01  9.472E-04   39.19


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

   24.469   2.881    21.40   1.30      2.13      1.77      70     107      26
    2.881   2.044     0.04   0.98   1150.03   1099.74     190     240       6
    2.044   1.671    -0.03   1.09    -99.90    -99.90     222     282       4
    1.671   1.448    -0.08   0.87    -99.90    -99.90     244     325       6
    1.448   1.296    -0.12   0.91    -99.90    -99.90     148     212       2
   --------------------------------------------------------------------------
   24.469   1.296     1.92   0.99      8.92      8.59     874    1166      44
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       4       4)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.682  -0.1326E+02  0.1074E+02        -1.24           1 
   47    0    0     1.646  -0.7809E+01  0.1060E+02        -0.74           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -7.1% OF MEAN INTENSITY


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

     3.84          45      32      1285        2.5%       1.8%      1.5%       26   35.63      2.5%    99.9*     0    0.000       0
     2.74          82      53      2098        2.5%       5.6%      4.3%       58   14.27      8.0%    99.5*     0    0.000       0
     2.24         141      87      2665        3.3%     627.6%    625.3%      108    0.08    887.6%    19.7      0    0.000       0
     1.94         133      79      3101        2.5%   17663.5%  15256.8%      108   -0.08  24979.9%   -19.7      0    0.000       0
     1.74         151      90      3479        2.6%    1263.0%   1278.8%      122    0.01   1786.1%     0.1      0    0.000       0
     1.59         172     105      3851        2.7%     -99.9%    -99.9%      134    0.00    -99.9%     7.3      0    0.000       0
     1.47         186     118      4143        2.8%     -99.9%    -99.9%      136    0.00    -99.9%    13.0      0    0.000       0
     1.37         175     107      4427        2.4%     -99.9%    -99.9%      136    0.00    -99.9%    -1.6      0    0.000       0
     1.30          69      51      4743        1.1%     -99.9%    -99.9%       36    0.00    -99.9%   -34.9      0    0.000       0
    total        1154     722     29792        2.4%       8.8%      8.5%      864    2.61     12.5%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1211
 NUMBER OF REJECTED MISFITS                                44
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1166
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    729



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

   WILSON LINE (using all data) : A=   6.543 B=  33.405 CORRELATION=  0.84
      #      RES      SS        <I>       log(<I>)       BO
      13     7.539   0.004  2.2603E+03       7.723    -134.2
      16     4.759   0.011  2.0688E+03       7.635     -49.5
      21     3.526   0.020  6.7747E+02       6.518       0.6
      28     3.001   0.028  1.2870E+02       4.857      30.3
      33     2.638   0.036 -1.8246E+00    -999.999   -1000.0
      48     2.398   0.043  1.3231E+00       0.280      72.0
      34     2.209   0.051  2.4346E+00       0.890      55.1
      40     2.056   0.059 -6.2176E-01    -999.999   -1000.0
      39     1.931   0.067 -2.0588E+00    -999.999   -1000.0
      50     1.829   0.075  5.5931E-01      -0.581      47.7
      43     1.739   0.083  1.3525E+00       0.302      37.7
      46     1.667   0.090 -7.1857E-01    -999.999   -1000.0
      46     1.598   0.098 -8.6458E-01    -999.999   -1000.0
      61     1.536   0.106  1.3682E+00       0.313      29.4
      63     1.481   0.114  1.1434E+00       0.134      28.1
      48     1.433   0.122  9.9620E-01      -0.004      26.9
      43     1.391   0.129 -1.4213E+00    -999.999   -1000.0
      36     1.352   0.137 -7.3647E-01    -999.999   -1000.0
      14     1.316   0.144 -6.2454E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.539        1.842        1.356        0.780
       3     4.759        1.104        0.618        0.259
       2     3.526        0.574        0.151        0.028
       3     3.001        1.770        1.411        0.804
       5     2.398       11.179        7.658       22.347
      10     1.829      129.764      276.067     3730.448
       7     1.739       17.886       50.372      125.781
      10     1.667       22.969       23.782       95.432
       3     1.536       15.959       33.530       55.246
       5     1.481       45.849       34.738      237.704
       1     1.352      138.293      563.364     1639.272
      74   overall       25.513       45.625      543.697


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

       9     7.539        0.473        0.306        0.155
      13     4.759        1.244        1.490        1.416
      19     3.526        1.019        1.072        1.020
      25     3.001        1.743        2.501        2.910
      43     2.398       19.980       37.660      169.949
      33     2.209        8.139       11.313       31.746
      34     2.056       31.865       27.128      465.516
      40     1.829      124.426      670.380    13386.643
      36     1.739       11.853        2.809       53.140
      58     1.536       22.134       32.599      303.797
      58     1.481       57.571      215.211     3335.480
      46     1.433       44.950       80.168     1400.623
     648   overall       20.533       72.646     1267.520

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.539  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.38  0.73
       3     4.759  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.526  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     3.001  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       5     2.398  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
      10     1.829  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
       7     1.739  1.15  0.83  0.69  0.60  0.55  0.51  0.48  0.45  0.43  0.42
      10     1.667  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
       3     1.536  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       5     1.481  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.59
       1     1.352  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      74   overall  1.31  0.94  0.78  0.69  0.62  0.58  0.54  0.54  0.51  0.53


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       9     7.539  2.34  3.06  2.14  2.02  1.98  1.72  1.55  1.41  1.31  1.23
      13     4.759  1.62  2.12  1.48  1.17  1.37  1.53  1.38  1.26  1.30  1.34
      19     3.526  0.55  1.16  0.81  0.80  0.94  1.17  1.25  1.24  1.15  1.08
      25     3.001  5.04  3.31  2.31  1.82  1.63  1.42  1.43  1.31  1.28  1.20
      43     2.398  5.38  2.82  1.97  1.55  1.30  1.13  1.02  0.93  0.86  0.85
      33     2.209  5.41  2.84  1.99  1.65  1.46  1.34  1.20  1.10  1.02  0.96
      34     2.056  4.64  2.60  1.82  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      40     1.829  5.52  2.90  2.03  1.59  1.33  1.16  1.04  0.95  0.88  0.83
      36     1.739  5.25  2.76  2.04  1.60  1.34  1.17  1.05  0.96  0.89  0.88
      58     1.536  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.85
      58     1.481  5.44  2.85  2.00  1.57  1.31  1.15  1.03  0.94  0.87  0.85
      46     1.433  4.80  2.52  1.76  1.38  1.16  1.06  0.95  0.87  0.81  0.76
     648   overall  3.08  1.73  1.21  0.97  0.84  0.76  0.69  0.64  0.60  0.57


 List of     68 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   34   22    4    1.87   52.49  0.2936E+02  0.1378E+02 "alien"
   36   26    6    1.68   46.51 -0.3342E+02  0.1362E+02 "alien"
   49    4   10    1.46   42.75  0.4889E+02  0.1581E+02 "alien"
   26   20   13    1.85   39.47  0.2208E+02  0.6952E+01 "alien"
   28   19   23    1.36   32.45 -0.2390E+02  0.1016E+02 "alien"
   41   33    6    1.43   30.65  0.3054E+02  0.1465E+02 "alien"
   35    5   17    1.58   28.11 -0.2430E+02  0.1149E+02 "alien"
   37   10   21    1.36   27.13 -0.1998E+02  0.1179E+02 "alien"
   32   20    2    2.04   27.00 -0.1679E+02  0.8319E+01 "alien"
   42    2   11    1.63   24.48 -0.1759E+02  0.7213E+01 "alien"
   27   25   19    1.46   24.03  0.2747E+02  0.1126E+02 "alien"
   36   24    1    1.79   23.66  0.1323E+02  0.7659E+01 "alien"
   39   30    6    1.53   21.70  0.2969E+02  0.1469E+02 "alien"
   54    5    8    1.37   21.08 -0.1552E+02  0.1336E+02 "alien"
   36    1    7    2.00   20.87 -0.1298E+02  0.1260E+02 "alien"
   28   19    9    2.02   20.71 -0.1288E+02  0.1221E+02 "alien"
   44    1    8    1.65   20.12 -0.1446E+02  0.7333E+01 "alien"
   51    4   10    1.41   18.99 -0.2699E+02  0.1511E+02 "alien"
   35   11   22    1.35   18.98 -0.1398E+02  0.1058E+02 "alien"
   23   15   22    1.49   18.90  0.2161E+02  0.8996E+01 "alien"
   40   10   20    1.35   18.83 -0.1387E+02  0.1352E+02 "alien"
   35   26    9    1.64   18.79 -0.1350E+02  0.7540E+01 "alien"
   43    1    1    1.80   18.66  0.1043E+02  0.7361E+01 "alien"
   54    3    3    1.42   17.95  0.1788E+02  0.8625E+01 "alien"
   34   26   10    1.64   17.28 -0.1241E+02  0.7613E+01 "alien"
   30   22   11    1.79   16.63  0.9303E+01  0.7275E+01 "alien"
   23   21   19    1.58   16.02 -0.1385E+02  0.6509E+01 "alien"
   20   18   16    1.85   15.81  0.8840E+01  0.8435E+01 "alien"
   45    4   13    1.48   15.79  0.1805E+02  0.1159E+02 "alien"
   54    4    7    1.38   15.63 -0.2221E+02  0.8609E+01 "alien"
   28   16    3    2.36   15.42  0.2041E+02  0.7876E+01 "alien"
   36    9   20    1.42   14.91  0.1486E+02  0.1138E+02 "alien"
   43    1    0    1.80   14.36  0.2171E+02  0.1097E+02 "alien"
   31   27   20    1.35   14.32 -0.1055E+02  0.1232E+02 "alien"
   34    1    8    2.06   14.18 -0.8819E+01  0.8120E+01 "alien"
   44    2   10    1.60   13.98 -0.1209E+02  0.7158E+01 "alien"
   37   26    3    1.70   13.62 -0.9790E+01  0.8745E+01 "alien"
   50    2    6    1.50   13.51  0.1545E+02  0.8440E+01 "alien"
   38   28    6    1.59   13.20 -0.1141E+02  0.1114E+02 "alien"
   36   26    7    1.66   13.05 -0.9374E+01  0.7386E+01 "alien"
   54    3    1    1.43   12.99  0.1294E+02  0.8583E+01 "alien"
   24   15    9    2.31   12.97  0.1716E+02  0.1167E+02 "alien"
   37    2    4    2.04   12.89 -0.8013E+01  0.7742E+01 "alien"
   38    2    0    2.03   12.84 -0.1940E+02  0.1151E+02 "alien"
   43    5   15    1.47   12.78  0.1461E+02  0.1215E+02 "alien"
   25   20   14    1.82   12.44  0.6958E+01  0.7525E+01 "alien"
   29    2    7    2.40   12.23  0.1619E+02  0.7321E+01 "alien"
   43    1    2    1.79   12.19  0.6820E+01  0.7621E+01 "alien"
   27    1   10    2.30   12.09  0.1600E+02  0.9752E+01 "alien"
   30   21    9    1.90   12.07 -0.2485E+02  0.1285E+02 "alien"
   26   14    2    2.60   11.91 -0.2173E+02  0.1303E+02 "alien"
   38    0    9    1.84   11.85  0.1507E+02  0.7471E+01 "alien"
   40    2   12    1.66   11.83 -0.8500E+01  0.7838E+01 "alien"
   40   29    0    1.57   11.74 -0.2291E+02  0.1116E+02 "alien"
   53    3    6    1.42   11.54  0.1150E+02  0.1573E+02 "alien"
   27    3    5    2.67   11.34 -0.2069E+02  0.1057E+02 "alien"
   40   30    2    1.54   11.26  0.1541E+02  0.1522E+02 "alien"
   53    5   10    1.36   11.10 -0.8172E+01  0.9216E+01 "alien"
   45    3   11    1.54   11.03  0.1509E+02  0.7953E+01 "alien"
   26   24   19    1.49   10.82  0.1238E+02  0.1053E+02 "alien"
   35   28   11    1.55   10.82  0.1480E+02  0.7826E+01 "alien"
   27   18   23    1.38   10.74 -0.1526E+02  0.9917E+01 "alien"
   53    6   11    1.34   10.62 -0.7825E+01  0.9069E+01 "alien"
   32   27   14    1.54   10.50  0.1437E+02  0.1254E+02 "alien"
   31   30   18    1.38   10.49 -0.1491E+02  0.9361E+01 "alien"
   40   30    3    1.54   10.45  0.1430E+02  0.8470E+01 "alien"
   29   22   22    1.36   10.03 -0.7386E+01  0.7816E+01 "alien"
   47    4   12    1.46   10.03  0.1147E+02  0.8581E+01 "alien"


 List of     68 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    6   11    1.34   10.62 -0.7825E+01  0.9069E+01
   40   10   20    1.35   18.83 -0.1387E+02  0.1352E+02
   31   27   20    1.35   14.32 -0.1055E+02  0.1232E+02
   35   11   22    1.35   18.98 -0.1398E+02  0.1058E+02
   29   22   22    1.36   10.03 -0.7386E+01  0.7816E+01
   28   19   23    1.36   32.45 -0.2390E+02  0.1016E+02
   53    5   10    1.36   11.10 -0.8172E+01  0.9216E+01
   37   10   21    1.36   27.13 -0.1998E+02  0.1179E+02
   54    5    8    1.37   21.08 -0.1552E+02  0.1336E+02
   27   18   23    1.38   10.74 -0.1526E+02  0.9917E+01
   31   30   18    1.38   10.49 -0.1491E+02  0.9361E+01
   54    4    7    1.38   15.63 -0.2221E+02  0.8609E+01
   51    4   10    1.41   18.99 -0.2699E+02  0.1511E+02
   36    9   20    1.42   14.91  0.1486E+02  0.1138E+02
   53    3    6    1.42   11.54  0.1150E+02  0.1573E+02
   54    3    3    1.42   17.95  0.1788E+02  0.8625E+01
   54    3    1    1.43   12.99  0.1294E+02  0.8583E+01
   41   33    6    1.43   30.65  0.3054E+02  0.1465E+02
   49    4   10    1.46   42.75  0.4889E+02  0.1581E+02
   47    4   12    1.46   10.03  0.1147E+02  0.8581E+01
   27   25   19    1.46   24.03  0.2747E+02  0.1126E+02
   43    5   15    1.47   12.78  0.1461E+02  0.1215E+02
   45    4   13    1.48   15.79  0.1805E+02  0.1159E+02
   26   24   19    1.49   10.82  0.1238E+02  0.1053E+02
   23   15   22    1.49   18.90  0.2161E+02  0.8996E+01
   50    2    6    1.50   13.51  0.1545E+02  0.8440E+01
   39   30    6    1.53   21.70  0.2969E+02  0.1469E+02
   40   30    3    1.54   10.45  0.1430E+02  0.8470E+01
   32   27   14    1.54   10.50  0.1437E+02  0.1254E+02
   45    3   11    1.54   11.03  0.1509E+02  0.7953E+01
   40   30    2    1.54   11.26  0.1541E+02  0.1522E+02
   35   28   11    1.55   10.82  0.1480E+02  0.7826E+01
   40   29    0    1.57   11.74 -0.2291E+02  0.1116E+02
   23   21   19    1.58   16.02 -0.1385E+02  0.6509E+01
   35    5   17    1.58   28.11 -0.2430E+02  0.1149E+02
   38   28    6    1.59   13.20 -0.1141E+02  0.1114E+02
   44    2   10    1.60   13.98 -0.1209E+02  0.7158E+01
   42    2   11    1.63   24.48 -0.1759E+02  0.7213E+01
   34   26   10    1.64   17.28 -0.1241E+02  0.7613E+01
   35   26    9    1.64   18.79 -0.1350E+02  0.7540E+01
   44    1    8    1.65   20.12 -0.1446E+02  0.7333E+01
   40    2   12    1.66   11.83 -0.8500E+01  0.7838E+01
   36   26    7    1.66   13.05 -0.9374E+01  0.7386E+01
   36   26    6    1.68   46.51 -0.3342E+02  0.1362E+02
   37   26    3    1.70   13.62 -0.9790E+01  0.8745E+01
   36   24    1    1.79   23.66  0.1323E+02  0.7659E+01
   30   22   11    1.79   16.63  0.9303E+01  0.7275E+01
   43    1    2    1.79   12.19  0.6820E+01  0.7621E+01
   43    1    1    1.80   18.66  0.1043E+02  0.7361E+01
   43    1    0    1.80   14.36  0.2171E+02  0.1097E+02
   25   20   14    1.82   12.44  0.6958E+01  0.7525E+01
   38    0    9    1.84   11.85  0.1507E+02  0.7471E+01
   26   20   13    1.85   39.47  0.2208E+02  0.6952E+01
   20   18   16    1.85   15.81  0.8840E+01  0.8435E+01
   34   22    4    1.87   52.49  0.2936E+02  0.1378E+02
   30   21    9    1.90   12.07 -0.2485E+02  0.1285E+02
   36    1    7    2.00   20.87 -0.1298E+02  0.1260E+02
   28   19    9    2.02   20.71 -0.1288E+02  0.1221E+02
   38    2    0    2.03   12.84 -0.1940E+02  0.1151E+02
   32   20    2    2.04   27.00 -0.1679E+02  0.8319E+01
   37    2    4    2.04   12.89 -0.8013E+01  0.7742E+01
   34    1    8    2.06   14.18 -0.8819E+01  0.8120E+01
   27    1   10    2.30   12.09  0.1600E+02  0.9752E+01
   24   15    9    2.31   12.97  0.1716E+02  0.1167E+02
   28   16    3    2.36   15.42  0.2041E+02  0.7876E+01
   29    2    7    2.40   12.23  0.1619E+02  0.7321E+01
   26   14    2    2.60   11.91 -0.2173E+02  0.1303E+02
   27    3    5    2.67   11.34 -0.2069E+02  0.1057E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        28
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1211
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1127
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       84


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
