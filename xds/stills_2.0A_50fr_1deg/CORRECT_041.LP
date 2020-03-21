
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.603    77.603    38.844  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      41      41
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      41
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001008  0.000645
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
 BEAM_DIVERGENCE_E.S.D.=   0.048
 REFLECTING_RANGE_E.S.D.=   0.187
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

    1139 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      41 ...      41
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1139 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               38
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

 *  31        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0  1  0  0  1  0  0  0
 *  44        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  35        mP          0.0      77.6   38.8   77.6  90.0  90.0  90.0    0  1  0  0  0  0 -1  0 -1  0  0  0
 *  33        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  34        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  20        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  25        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1 -1  0  0  1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0  1  0  0  0  0  0  1  0
    28        mC        250.0      38.8  160.0   77.6  90.0  90.0  75.9    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.8  160.0   77.6  90.0  90.0  75.9    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    37        mC        250.0     160.0   38.8   77.6  90.0  90.0  75.9    2  0  1  0  0  0  1  0  0 -1  0  0
    39        mC        250.0     160.0   38.8   77.6  90.0  90.0  75.9    0 -2 -1  0  0  0 -1  0  1  0  0  0
    36        oC        250.0      38.8  160.0   77.6  90.0  90.0 104.1    0  0  1  0 -2  0 -1  0  0 -1  0  0
    38        oC        250.0      38.8  160.0   77.6  90.0  90.0 104.1    0  0 -1  0  0  2  1  0  1  0  0  0
    27        mC        500.0     160.0   38.8  109.7  90.0 133.3  75.9    0  2  1  0  0  0  1  0  1 -1  0  0
    19        oI        500.0      38.8  109.7  116.4  90.0  70.5  90.0    0  0  1  0  1 -1  0  0  1  1  1  0
    26        oF        625.0      38.8  160.0  160.0  86.6 104.1 104.1    0  0 -1  0  0  2  1  0  2  0  1  0
    18        tI        625.0     109.7  116.4   38.8  70.5  90.0  90.0    1 -1  0  0 -1 -1 -1  0  0  0 -1  0
    30        mC        997.9      77.6  173.5   38.8  90.0  90.0  63.4    0  1  0  0 -2  1  0  0  0  0  1  0
    41        mC        997.9     173.5   77.6   38.8  90.0  90.0  63.4   -2  1  0  0  0  1  0  0  0  0 -1  0
    40        oC        997.9      77.6  173.5   38.8  90.0  90.0 116.6    0  1  0  0  2 -1  0  0  0  0 -1  0
    22        hP        997.9      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0  1  0  0  0  0  0  1  0
     1        cF        999.0     116.4  116.4  116.4  96.4  96.4 141.0    1 -1 -1  0 -1  1 -1  0  1  1  1  0
     2        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
     5        cI        999.0      86.8   86.8  109.7  50.8  50.8  78.4    1  0  1  0  0 -1  1  0  1 -1  0  0
     4        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6    0  1  1  0  1  0 -1  0 -1  1 -1  0
     6        tI        999.0     109.7   86.8   86.8  78.4  50.8  50.8    1  1  0  0  1  0 -1  0  0  1 -1  0
     7        tI        999.0      86.8   86.8  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.8   86.8  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.8   86.8  248.5 102.1  99.0 116.6    0  0 -1  0  0  1  1  0  3 -1  1  0
    10        mC        999.0      86.8   86.8   77.6  90.0  90.0 126.8    0  1 -1  0  0 -1 -1  0 -1  0  0  0
    11        tP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.8   77.6  177.8  64.1  77.4  90.0    0  0  1  0  0 -1  0  0  2 -1  1  0
    16        oF        999.0      86.8   86.8  177.8 107.0 119.2  53.2    0 -1  1  0  0 -1 -1  0  2  1 -1  0
    14        mC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    0 -1  1  0  0 -1 -1  0  1  0  0  0
    17        mC        999.0      86.8   86.8   86.8  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.8  109.7   38.8  90.0  77.4 108.0    1  2 -1  0  1 -1  0  0  0  0 -1  0
    42        oI        999.0      38.8   77.6  177.8 115.9 102.6  90.0    0  0 -1  0  0  1  0  0  2 -1  1  0
    43        mI        999.0      86.8  177.8   77.6 115.9 153.4  60.8    0  1 -1  0  2  1 -1  0  0 -1  0  0

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

  *   0.00        0    4.5      14    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.922   -23.622     0.100
 UNIT_CELL_B-AXIS=    23.622    73.922     0.041
 UNIT_CELL_C-AXIS=    -0.054    -0.004    38.842
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       118 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.31
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.633    77.633    38.846  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.6E-01 1.6E-01 7.8E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012881  0.012881  0.025742  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.985   -23.517     0.037
 COORDINATES OF UNIT CELL B-AXIS    23.517    73.985     0.015
 COORDINATES OF UNIT CELL C-AXIS    -0.020    -0.001    38.846
 CRYSTAL MOSAICITY (DEGREES)     0.187
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006347  0.002001
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.90   1232.87
 DETECTOR ORIGIN (PIXELS) AT                     1260.98   1242.11
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.42
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.63    77.63    38.85  90.000  90.000  90.000


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
  0.05  0.15        15          3.936      25.855
  0.15  0.25       141        104.432       0.975
  0.25  0.35       160        105.183       0.968
  0.35  0.45       158         49.347       2.062
  0.45  0.55       167        131.183       0.776
  0.55  0.65       169        102.840       0.990
  0.65  0.75       195        126.358       0.805
  0.75  0.85       123         99.943       1.018
  0.85  0.95        11         26.733       3.807
 
 -0.05  0.95      1139        101.778       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               149.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.866
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    40.1 XMAX=    40.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59181 NYBIN=    6
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        318


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               149.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.862
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    55.9 XMAX=  2451.2 NXBIN=    2
 YMIN=     6.9 YMAX=  2519.1 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        318


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               145.3
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.861
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    40.1 XMAX=    40.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1254    1263
 DETECTOR_SURFACE_POSITION=    1660    1690
 DETECTOR_SURFACE_POSITION=     847    1690
 DETECTOR_SURFACE_POSITION=     847     836
 DETECTOR_SURFACE_POSITION=    1660     836
 DETECTOR_SURFACE_POSITION=    2174    1663
 DETECTOR_SURFACE_POSITION=    1635    2229
 DETECTOR_SURFACE_POSITION=     872    2229
 DETECTOR_SURFACE_POSITION=     333    1663
 DETECTOR_SURFACE_POSITION=     333     863
 DETECTOR_SURFACE_POSITION=     872     297
 DETECTOR_SURFACE_POSITION=    1635     297
 DETECTOR_SURFACE_POSITION=    2174     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        318


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
 7.649E-01  1.378E-03   30.80


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

   24.550   2.890    17.47   1.17      2.20      2.20      80     108      18
    2.890   2.051     2.84   1.19     24.83     21.79     180     221       4
    2.051   1.677     0.06   1.05    730.95    695.88     218     267       2
    1.677   1.453     0.00   0.74  27445.24  31948.11     238     300       4
    1.453   1.300    -0.03   0.79   1387.24   1589.64     156     208       6
   --------------------------------------------------------------------------
   24.550   1.300     2.29   0.96      9.19      9.38     872    1104      34
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      41      41)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.688  -0.4417E+00  0.1298E+02        -0.03           1 
   47    0    0     1.652  -0.1303E+02  0.1200E+02        -1.09           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS   -12.2% OF MEAN INTENSITY


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

     3.86          56      36      1285        2.8%       2.0%      2.0%       40   27.54      2.8%    99.9*     0    0.000       0
     2.74          71      45      2098        2.1%       4.6%      4.7%       52   13.99      6.5%    99.5*     0    0.000       0
     2.25         138      79      2665        3.0%      20.6%     17.1%      118    4.18     29.2%    92.0*     0    0.000       0
     1.95         120      72      3100        2.3%     140.2%    131.7%       96    0.57    198.3%    52.9*     0    0.000       0
     1.74         138      82      3481        2.4%    1371.6%   1327.3%      112    0.04   1939.7%     3.8      0    0.000       0
     1.59         175     104      3850        2.7%    2273.9%   2832.4%      142    0.06   3215.7%     8.5      0    0.000       0
     1.47         162      99      4144        2.4%     -99.9%    -99.9%      126    0.00    -99.9%     9.8      0    0.000       0
     1.38         168      99      4427        2.2%    1942.5%   2283.0%      138   -0.03   2747.1%     3.6      0    0.000       0
     1.30          75      51      4741        1.1%     296.4%    368.5%       48   -0.01    419.2%    13.7      0    0.000       0
    total        1103     667     29791        2.2%       9.2%      9.4%      872    3.00     13.0%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1139
 NUMBER OF REJECTED MISFITS                                34
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1104
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    668



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

   WILSON LINE (using all data) : A=   7.082 B=  35.195 CORRELATION=  0.96
      #      RES      SS        <I>       log(<I>)       BO
      13     7.329   0.005  2.0289E+03       7.615     -57.3
      19     4.764   0.011  1.4409E+03       7.273      -8.7
      18     3.588   0.019  5.5496E+02       6.319      19.7
      25     2.990   0.028  1.1768E+02       4.768      41.4
      30     2.642   0.036  8.7966E+01       4.477      36.4
      43     2.411   0.043  3.3677E+01       3.517      41.5
      33     2.213   0.051  1.2815E+01       2.551      44.4
      34     2.062   0.059  4.6158E+00       1.529      47.2
      36     1.940   0.066 -1.5461E+00    -999.999   -1000.0
      44     1.833   0.074 -4.9057E-01    -999.999   -1000.0
      43     1.743   0.082  3.0320E+00       1.109      36.3
      42     1.672   0.089  1.6957E+00       0.528      36.7
      50     1.602   0.097  7.1719E-01      -0.332      38.1
      46     1.543   0.105 -8.9977E-01    -999.999   -1000.0
      57     1.486   0.113 -1.1426E+00    -999.999   -1000.0
      44     1.437   0.121  8.9452E-01      -0.111      29.7
      40     1.396   0.128  2.4955E-01      -1.388      33.0
      34     1.355   0.136 -6.0652E-01    -999.999   -1000.0
      15     1.319   0.144 -5.1089E-01    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.329        1.691        1.181        0.644
       3     4.764        1.581        1.085        0.553
       2     3.588        0.338        0.069        0.010
       3     2.990        2.367        1.378        0.598
       3     2.411        0.197        0.037        0.005
       1     2.213        0.269        0.048        0.006
       6     2.062        2.555        1.803        0.990
       8     1.833       64.922      382.613     1989.916
       7     1.743        3.061        2.094        1.103
      11     1.672        5.561        4.229        6.450
       3     1.602      116.231      346.386     2062.137
       3     1.437       22.525       39.413       85.950
      70   overall        7.633       17.740       93.356


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

       9     7.329        0.557        0.393        0.218
      16     4.764        1.332        1.822        1.987
      16     3.588        1.086        1.230        1.264
      22     2.990        0.653        0.502        0.331
      29     2.642        1.447        2.405        3.785
      40     2.411        1.038        0.969        0.793
      32     2.213        1.791        2.106        2.364
      28     2.062        3.622        2.587        4.757
      36     1.743        7.741        9.844       33.821
      31     1.672       10.587       -1.366       87.997
      47     1.602       58.529      -14.434     1420.770
      41     1.437       73.753       95.032     3153.636
      40     1.396      784.525    -6277.050   328149.183
     596   overall       63.864     -414.833    22359.731

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.329  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.76  0.73
       3     4.764  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     3.588  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.73
       3     2.990  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.411  0.00  0.00  0.00  0.00  0.64  1.19  1.12  1.06  1.01  0.98
       1     2.213  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  1.52  1.47
       6     2.062  0.67  0.48  0.40  0.35  0.64  0.59  0.56  0.53  0.51  0.49
       8     1.833  3.02  2.17  1.80  1.58  1.44  1.34  1.26  1.19  1.14  1.10
       7     1.743  0.58  0.41  0.34  0.30  0.27  0.25  0.24  0.23  0.22  0.42
      11     1.672  1.10  1.05  0.87  0.77  0.70  0.65  0.61  0.58  0.55  0.53
       3     1.602  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.437  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      70   overall  0.92  0.70  0.58  0.51  0.52  0.51  0.48  0.48  0.50  0.52


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       9     7.329  3.50  3.06  2.14  2.02  1.69  1.72  1.55  1.41  1.31  1.23
      16     4.764  3.28  1.72  1.69  1.33  1.27  1.25  1.49  1.36  1.37  1.38
      16     3.588  1.31  0.69  0.72  0.76  0.79  1.11  1.12  1.13  1.16  1.09
      22     2.990  1.43  1.00  1.40  1.10  1.39  1.41  1.44  1.32  1.30  1.22
      29     2.642  0.72  0.95  0.93  1.15  1.05  1.15  1.10  1.06  1.05  1.15
      40     2.411  1.05  0.97  0.87  0.83  0.95  0.94  0.94  1.00  0.97  1.03
      32     2.213  2.96  2.07  1.45  1.33  1.27  1.18  1.06  0.96  0.90  0.99
      28     2.062  4.13  2.17  1.65  1.41  1.27  1.11  0.99  0.97  0.90  0.90
      36     1.743  5.25  2.76  1.93  1.52  1.27  1.11  1.05  1.01  0.94  0.92
      31     1.672  4.41  2.31  1.62  1.27  1.15  1.00  0.96  0.88  0.87  0.82
      47     1.602  4.70  2.46  1.72  1.36  1.19  1.04  0.93  0.85  0.79  0.74
      41     1.437  4.87  2.56  1.79  1.41  1.18  1.03  0.97  0.89  0.82  0.77
      40     1.396  4.73  2.48  1.74  1.36  1.14  1.00  0.89  0.82  0.76  0.71
     596   overall  2.26  1.30  0.98  0.83  0.76  0.71  0.68  0.64  0.61  0.61


 List of     95 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   32   31   17    1.39   72.39  0.1806E+02  0.1623E+02 "alien"
   29   26   20    1.39   64.26  0.1604E+02  0.1462E+02 "alien"
   41   11   21    1.30   50.41 -0.2576E+02  0.1393E+02 "alien"
   48    6   14    1.39   47.84  0.1194E+02  0.1388E+02 "alien"
   52    5   10    1.39   46.52  0.1161E+02  0.1134E+02 "alien"
   53    4    7    1.41   46.51  0.1161E+02  0.1166E+02 "alien"
   45    7   16    1.40   44.85  0.1119E+02  0.1527E+02 "alien"
   35   32   15    1.38   42.31  0.1056E+02  0.1618E+02 "alien"
   52    5   11    1.37   38.93 -0.2361E+02  0.1749E+02 "alien"
   28   21   22    1.38   38.29  0.9555E+01  0.8832E+01 "alien"
   35   11   22    1.36   37.08 -0.2249E+02  0.1267E+02 "alien"
   36   34   15    1.34   34.83 -0.2113E+02  0.1326E+02 "alien"
   37    1   11    1.80   33.81 -0.1659E+02  0.8454E+01 "alien"
   31   30   17    1.41   31.79  0.7932E+01  0.9761E+01 "alien"
   37   31   12    1.44   30.54  0.2732E+02  0.1522E+02 "alien"
   28   23   21    1.40   28.65  0.7150E+01  0.8503E+01 "alien"
   38   32   11    1.43   28.52  0.2551E+02  0.9644E+01 "alien"
   50    7   14    1.34   28.19 -0.1710E+02  0.1133E+02 "alien"
   31   26   21    1.33   27.37 -0.1398E+02  0.1327E+02 "alien"
   38   29    7    1.56   27.28 -0.2455E+02  0.1039E+02 "alien"
   48    9   17    1.30   27.28 -0.1393E+02  0.1775E+02 "alien"
   40    1    5    1.88   26.25 -0.1288E+02  0.1245E+02 "alien"
   40    1    6    1.86   25.48 -0.1250E+02  0.9078E+01 "alien"
   51    6   13    1.35   24.74 -0.1501E+02  0.1695E+02 "alien"
   37   33   13    1.39   24.13  0.6021E+01  0.1024E+02 "alien"
   34   29   13    1.50   24.04 -0.2747E+02  0.1402E+02 "alien"
   28   21   12    1.83   23.87 -0.1171E+02  0.7525E+01 "alien"
   38   33   12    1.39   23.50  0.5864E+01  0.9321E+01 "alien"
   21   19   16    1.82   23.33 -0.1145E+02  0.1020E+02 "alien"
   52    4    9    1.41   22.99  0.5738E+01  0.9858E+01 "alien"
   35   23    1    1.85   22.48 -0.1103E+02  0.8183E+01 "alien"
   31   22   10    1.81   22.41 -0.1099E+02  0.1237E+02 "alien"
   36   10   21    1.38   21.79  0.5437E+01  0.8768E+01 "alien"
   35    1   12    1.83   21.75 -0.1067E+02  0.8217E+01 "alien"
   37   29   10    1.52   21.64 -0.1947E+02  0.1527E+02 "alien"
   35   23    3    1.83   21.16 -0.1038E+02  0.8236E+01 "alien"
   33    7   19    1.53   20.81 -0.1872E+02  0.7681E+01 "alien"
   22   18   20    1.58   20.64  0.1481E+02  0.1027E+02 "alien"
   17   15   18    1.83   20.62 -0.1012E+02  0.8540E+01 "alien"
   51    7   14    1.33   20.02 -0.1023E+02  0.1405E+02 "alien"
   35   23    2    1.85   19.75 -0.9689E+01  0.8197E+01 "alien"
   39   33   11    1.40   18.95  0.4729E+01  0.1032E+02 "alien"
   27   24   20    1.44   18.19  0.1627E+02  0.1267E+02 "alien"
   31    9   21    1.47   18.08 -0.2066E+02  0.1213E+02 "alien"
   33    6   18    1.58   18.07  0.1296E+02  0.7850E+01 "alien"
   27   21   13    1.81   17.34 -0.8506E+01  0.7326E+01 "alien"
   25   24   17    1.60   17.32  0.1242E+02  0.1360E+02 "alien"
   50    2    6    1.51   16.94 -0.1935E+02  0.9139E+01 "alien"
   46    8   17    1.34   16.61 -0.1008E+02  0.9712E+01 "alien"
   35   31   14    1.42   16.29  0.1457E+02  0.9296E+01 "alien"
   40   29    2    1.57   16.26 -0.1463E+02  0.9307E+01 "alien"
   31   26   14    1.58   16.24  0.1165E+02  0.8681E+01 "alien"
   42   32    3    1.46   15.80  0.1413E+02  0.1143E+02 "alien"
   37   32   12    1.42   15.51  0.1387E+02  0.1644E+02 "alien"
   42    1    3    1.83   15.29 -0.7502E+01  0.1266E+02 "alien"
   28   24   15    1.63   15.07  0.1081E+02  0.7911E+01 "alien"
   51    2    4    1.50   14.66 -0.1675E+02  0.8876E+01 "alien"
   38   28    5    1.61   14.47  0.1038E+02  0.1593E+02 "alien"
   48    7   15    1.36   14.37 -0.8718E+01  0.1032E+02 "alien"
   54    4    6    1.40   14.01  0.3496E+01  0.1057E+02 "alien"
   36   33   14    1.38   13.98  0.3489E+01  0.1561E+02 "alien"
   29    1   13    1.99   13.88 -0.2146E+02  0.8194E+01 "alien"
   25   20   14    1.83   13.85 -0.6793E+01  0.7399E+01 "alien"
   31   29   16    1.46   13.37  0.1196E+02  0.9754E+01 "alien"
   32   30   19    1.34   13.20 -0.8008E+01  0.1343E+02 "alien"
   31   30   18    1.38   12.96  0.3233E+01  0.1493E+02 "alien"
   43    3   12    1.57   12.87  0.9228E+01  0.8616E+01 "alien"
   35    5   17    1.58   12.85  0.9217E+01  0.1312E+02 "alien"
   39   29    6    1.55   12.57 -0.1131E+02  0.9686E+01 "alien"
   53    2    0    1.46   12.36 -0.3322E+02  0.1548E+02 "alien"
   45    6   15    1.43   12.35  0.1105E+02  0.1046E+02 "alien"
   43    2   11    1.61   12.27  0.8799E+01  0.1358E+02 "alien"
   55    3    1    1.41   12.15  0.3033E+01  0.9012E+01 "alien"
   40   29    0    1.57   11.96  0.1966E+02  0.1242E+02 "alien"
   25   20   21    1.47   11.84 -0.1353E+02  0.1093E+02 "alien"
   37   29    9    1.54   11.81 -0.1063E+02  0.9857E+01 "alien"
   33   33   17    1.34   11.77 -0.1615E+02  0.1389E+02 "alien"
   35    6   18    1.54   11.55 -0.1039E+02  0.1329E+02 "alien"
   24   22   19    1.55   11.38 -0.1024E+02  0.1180E+02 "alien"
   42    8   18    1.39   11.38  0.2839E+01  0.9856E+01 "alien"
   34   23    6    1.82   11.29 -0.5539E+01  0.8093E+01 "alien"
   34    3   15    1.71   11.21  0.3400E+02  0.1418E+02 "alien"
   31    8   20    1.52   11.14 -0.1002E+02  0.7871E+01 "alien"
   52    2    4    1.47   11.13 -0.1272E+02  0.1683E+02 "alien"
   34   22    4    1.88   11.07 -0.5430E+01  0.9847E+01 "alien"
   28   23   14    1.70   11.01  0.1867E+02  0.7586E+01 "alien"
   50    1    0    1.55   10.93 -0.2106E+02  0.1280E+02 "alien"
   30   25   21    1.35   10.78 -0.6538E+01  0.9479E+01 "alien"
   24   23   17    1.63   10.74  0.7706E+01  0.1207E+02 "alien"
   35   34   16    1.33   10.62 -0.5425E+01  0.1360E+02 "alien"
   21   17   20    1.61   10.33  0.7409E+01  0.9418E+01 "alien"
   36    0   10    1.89   10.29 -0.1037E+02  0.8250E+01 "alien"
   41    3   13    1.60   10.26  0.7361E+01  0.1209E+02 "alien"
   51    2    5    1.49   10.21 -0.1167E+02  0.9796E+01 "alien"
   50    3    8    1.48   10.12 -0.1156E+02  0.8833E+01 "alien"


 List of     95 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   41   11   21    1.30   50.41 -0.2576E+02  0.1393E+02
   48    9   17    1.30   27.28 -0.1393E+02  0.1775E+02
   51    7   14    1.33   20.02 -0.1023E+02  0.1405E+02
   35   34   16    1.33   10.62 -0.5425E+01  0.1360E+02
   31   26   21    1.33   27.37 -0.1398E+02  0.1327E+02
   32   30   19    1.34   13.20 -0.8008E+01  0.1343E+02
   36   34   15    1.34   34.83 -0.2113E+02  0.1326E+02
   46    8   17    1.34   16.61 -0.1008E+02  0.9712E+01
   33   33   17    1.34   11.77 -0.1615E+02  0.1389E+02
   50    7   14    1.34   28.19 -0.1710E+02  0.1133E+02
   51    6   13    1.35   24.74 -0.1501E+02  0.1695E+02
   30   25   21    1.35   10.78 -0.6538E+01  0.9479E+01
   35   11   22    1.36   37.08 -0.2249E+02  0.1267E+02
   48    7   15    1.36   14.37 -0.8718E+01  0.1032E+02
   52    5   11    1.37   38.93 -0.2361E+02  0.1749E+02
   36   33   14    1.38   13.98  0.3489E+01  0.1561E+02
   28   21   22    1.38   38.29  0.9555E+01  0.8832E+01
   36   10   21    1.38   21.79  0.5437E+01  0.8768E+01
   31   30   18    1.38   12.96  0.3233E+01  0.1493E+02
   35   32   15    1.38   42.31  0.1056E+02  0.1618E+02
   32   31   17    1.39   72.39  0.1806E+02  0.1623E+02
   37   33   13    1.39   24.13  0.6021E+01  0.1024E+02
   52    5   10    1.39   46.52  0.1161E+02  0.1134E+02
   48    6   14    1.39   47.84  0.1194E+02  0.1388E+02
   42    8   18    1.39   11.38  0.2839E+01  0.9856E+01
   29   26   20    1.39   64.26  0.1604E+02  0.1462E+02
   38   33   12    1.39   23.50  0.5864E+01  0.9321E+01
   45    7   16    1.40   44.85  0.1119E+02  0.1527E+02
   39   33   11    1.40   18.95  0.4729E+01  0.1032E+02
   54    4    6    1.40   14.01  0.3496E+01  0.1057E+02
   28   23   21    1.40   28.65  0.7150E+01  0.8503E+01
   52    4    9    1.41   22.99  0.5738E+01  0.9858E+01
   55    3    1    1.41   12.15  0.3033E+01  0.9012E+01
   53    4    7    1.41   46.51  0.1161E+02  0.1166E+02
   31   30   17    1.41   31.79  0.7932E+01  0.9761E+01
   35   31   14    1.42   16.29  0.1457E+02  0.9296E+01
   37   32   12    1.42   15.51  0.1387E+02  0.1644E+02
   45    6   15    1.43   12.35  0.1105E+02  0.1046E+02
   38   32   11    1.43   28.52  0.2551E+02  0.9644E+01
   37   31   12    1.44   30.54  0.2732E+02  0.1522E+02
   27   24   20    1.44   18.19  0.1627E+02  0.1267E+02
   31   29   16    1.46   13.37  0.1196E+02  0.9754E+01
   42   32    3    1.46   15.80  0.1413E+02  0.1143E+02
   53    2    0    1.46   12.36 -0.3322E+02  0.1548E+02
   31    9   21    1.47   18.08 -0.2066E+02  0.1213E+02
   25   20   21    1.47   11.84 -0.1353E+02  0.1093E+02
   52    2    4    1.47   11.13 -0.1272E+02  0.1683E+02
   50    3    8    1.48   10.12 -0.1156E+02  0.8833E+01
   51    2    5    1.49   10.21 -0.1167E+02  0.9796E+01
   34   29   13    1.50   24.04 -0.2747E+02  0.1402E+02
   51    2    4    1.50   14.66 -0.1675E+02  0.8876E+01
   50    2    6    1.51   16.94 -0.1935E+02  0.9139E+01
   31    8   20    1.52   11.14 -0.1002E+02  0.7871E+01
   37   29   10    1.52   21.64 -0.1947E+02  0.1527E+02
   33    7   19    1.53   20.81 -0.1872E+02  0.7681E+01
   35    6   18    1.54   11.55 -0.1039E+02  0.1329E+02
   37   29    9    1.54   11.81 -0.1063E+02  0.9857E+01
   39   29    6    1.55   12.57 -0.1131E+02  0.9686E+01
   24   22   19    1.55   11.38 -0.1024E+02  0.1180E+02
   50    1    0    1.55   10.93 -0.2106E+02  0.1280E+02
   38   29    7    1.56   27.28 -0.2455E+02  0.1039E+02
   40   29    2    1.57   16.26 -0.1463E+02  0.9307E+01
   40   29    0    1.57   11.96  0.1966E+02  0.1242E+02
   43    3   12    1.57   12.87  0.9228E+01  0.8616E+01
   31   26   14    1.58   16.24  0.1165E+02  0.8681E+01
   33    6   18    1.58   18.07  0.1296E+02  0.7850E+01
   22   18   20    1.58   20.64  0.1481E+02  0.1027E+02
   35    5   17    1.58   12.85  0.9217E+01  0.1312E+02
   41    3   13    1.60   10.26  0.7361E+01  0.1209E+02
   25   24   17    1.60   17.32  0.1242E+02  0.1360E+02
   43    2   11    1.61   12.27  0.8799E+01  0.1358E+02
   38   28    5    1.61   14.47  0.1038E+02  0.1593E+02
   21   17   20    1.61   10.33  0.7409E+01  0.9418E+01
   24   23   17    1.63   10.74  0.7706E+01  0.1207E+02
   28   24   15    1.63   15.07  0.1081E+02  0.7911E+01
   28   23   14    1.70   11.01  0.1867E+02  0.7586E+01
   34    3   15    1.71   11.21  0.3400E+02  0.1418E+02
   37    1   11    1.80   33.81 -0.1659E+02  0.8454E+01
   27   21   13    1.81   17.34 -0.8506E+01  0.7326E+01
   31   22   10    1.81   22.41 -0.1099E+02  0.1237E+02
   34   23    6    1.82   11.29 -0.5539E+01  0.8093E+01
   21   19   16    1.82   23.33 -0.1145E+02  0.1020E+02
   17   15   18    1.83   20.62 -0.1012E+02  0.8540E+01
   25   20   14    1.83   13.85 -0.6793E+01  0.7399E+01
   35    1   12    1.83   21.75 -0.1067E+02  0.8217E+01
   42    1    3    1.83   15.29 -0.7502E+01  0.1266E+02
   28   21   12    1.83   23.87 -0.1171E+02  0.7525E+01
   35   23    3    1.83   21.16 -0.1038E+02  0.8236E+01
   35   23    2    1.85   19.75 -0.9689E+01  0.8197E+01
   35   23    1    1.85   22.48 -0.1103E+02  0.8183E+01
   40    1    6    1.86   25.48 -0.1250E+02  0.9078E+01
   34   22    4    1.88   11.07 -0.5430E+01  0.9847E+01
   40    1    5    1.88   26.25 -0.1288E+02  0.1245E+02
   36    0   10    1.89   10.29 -0.1037E+02  0.8250E+01
   29    1   13    1.99   13.88 -0.2146E+02  0.8194E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        28
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1139
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1065
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       74


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
