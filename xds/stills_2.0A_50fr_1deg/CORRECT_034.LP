
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.561    77.561    38.817  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      34      34
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      34
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000974  0.000670
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
 BEAM_DIVERGENCE_E.S.D.=   0.047
 REFLECTING_RANGE_E.S.D.=   0.153
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

    1266 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      34 ...      34
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1266 REFLECTIONS ACCEPTED



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

 *  44        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.6   38.8   77.6  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  25        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  33        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  32        oP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  23        oC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1  1  0  0 -1 -1  0  0  0  0  1  0
 *  21        tP          0.0      77.6   77.6   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.9   38.8   77.6  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.9   38.8   77.6  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  159.9   77.6  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.9   77.6  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  159.9   77.6  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  159.9   77.6  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     159.9   38.8  109.7  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.7  116.4  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  159.9  159.9  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.7  116.4   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        998.1     173.4   77.6   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.1      77.6  173.4   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.1      77.6  173.4   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.1      77.6   77.6   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.4  116.4  116.4  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.7   86.7  116.4 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.7   86.7  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.7   86.7  116.4 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.7   86.7   86.7  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.7   86.7  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.7   86.7  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.7  248.3 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.7   86.7   77.6  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.7   86.7   77.6  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.6  177.7  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.7   86.7  177.7 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.7   86.7   77.6  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.7   86.7   86.7  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.7  109.7   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.6  177.7 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.7  177.7   77.6 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    4.2      17    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.905   -23.532     0.049
 UNIT_CELL_B-AXIS=    23.532    73.905     0.018
 UNIT_CELL_C-AXIS=    -0.026    -0.001    38.817
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       158 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.32
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.476    77.476    38.776  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.3E-01 2.3E-01 1.0E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012907  0.012907  0.025789  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.831   -23.484     0.067
 COORDINATES OF UNIT CELL B-AXIS    23.484    73.831     0.012
 COORDINATES OF UNIT CELL C-AXIS    -0.034     0.004    38.776
 CRYSTAL MOSAICITY (DEGREES)     0.153
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006177  0.002367
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.89   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1260.45   1241.86
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.10
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.48    77.48    38.78  90.000  90.000  90.000


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
  0.05  0.15        48        123.751       0.903
  0.15  0.25       161        111.158       1.005
  0.25  0.35       178        110.332       1.012
  0.35  0.45       166         54.443       2.052
  0.45  0.55       158        203.707       0.548
  0.55  0.65       185         78.257       1.427
  0.65  0.75       194        123.985       0.901
  0.75  0.85       140        111.974       0.998
  0.85  0.95        36         69.575       1.605
 
 -0.05  0.95      1266        111.699       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        8
 DEGREES OF FREEDOM OF CHI^2 FIT               197.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.980
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    33.1 XMAX=    33.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59411 NYBIN=    8
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        420


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               197.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.979
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    12.6 XMAX=  2452.2 NXBIN=    3
 YMIN=     5.9 YMAX=  2519.3 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        420


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               195.2
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.980
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    33.1 XMAX=    33.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1263
 DETECTOR_SURFACE_POSITION=    1647    1689
 DETECTOR_SURFACE_POSITION=     818    1689
 DETECTOR_SURFACE_POSITION=     818     836
 DETECTOR_SURFACE_POSITION=    1647     836
 DETECTOR_SURFACE_POSITION=    2170    1663
 DETECTOR_SURFACE_POSITION=    1621    2229
 DETECTOR_SURFACE_POSITION=     844    2229
 DETECTOR_SURFACE_POSITION=     295    1663
 DETECTOR_SURFACE_POSITION=     295     862
 DETECTOR_SURFACE_POSITION=     844     297
 DETECTOR_SURFACE_POSITION=    1621     297
 DETECTOR_SURFACE_POSITION=    2170     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        420


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
 8.071E-01  9.958E-04   35.27


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

   24.500   3.156    23.47   1.45      1.94      1.83      64      87      20
    3.156   2.241     6.42   1.41      9.26      7.36     172     201       8
    2.241   1.832     0.59   0.96     88.43     85.10     190     214       0
    1.832   1.588    -0.08   0.76    -99.90    -99.90     238     268       0
    1.588   1.421     0.02   0.79    847.68    963.82     238     292       6
    1.421   1.297     0.00   0.93    867.97    852.52     112     163       6
   --------------------------------------------------------------------------
   24.500   1.297     2.81   0.98      8.59      8.52    1014    1225      40
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      34      34)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.684   0.1800E+00  0.1214E+02         0.01           1 
   47    0    0     1.648   0.2088E+02  0.1268E+02         1.65           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    17.9% OF MEAN INTENSITY


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

     3.85          56      35      1285        2.7%       1.8%      1.7%       42   32.48      2.5%   100.0*     0    0.000       0
     2.74          84      50      2098        2.4%       4.9%      3.6%       68   17.88      7.0%    99.4*     0    0.000       0
     2.24         148      85      2665        3.2%      12.4%     10.5%      126    6.61     17.6%    97.0*     0    0.000       0
     1.94         137      76      3101        2.5%      54.8%     54.0%      122    1.30     77.5%    71.3*     0    0.000       0
     1.74         156      86      3481        2.5%     -99.9%    -99.9%      140    0.00    -99.9%    18.7      0    0.000       0
     1.59         187     105      3847        2.7%   22966.5%  23852.2%      164   -0.02  32479.5%     6.6      0    0.000       0
     1.47         186     111      4142        2.7%     887.6%    940.9%      150   -0.01   1255.3%   -10.8      0    0.000       0
     1.38         186     112      4431        2.5%     599.6%    675.0%      148    0.06    848.0%    22.6      0    0.000       0
     1.30          83      57      4743        1.2%   41911.7%  42792.9%       52   -0.11  59272.0%   -11.8      0    0.000       0
    total        1223     717     29793        2.4%       8.6%      8.5%     1012    3.73     12.1%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1266
 NUMBER OF REJECTED MISFITS                                40
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1225
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    718



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

   WILSON LINE (using all data) : A=   6.861 B=  30.651 CORRELATION=  0.89
      #      RES      SS        <I>       log(<I>)       BO
      14     7.177   0.005  1.9410E+03       7.571     -73.2
      17     4.705   0.011  1.7691E+03       7.478     -27.3
      21     3.578   0.020  6.4739E+02       6.473       9.9
      27     2.999   0.028  1.7626E+02       5.172      30.4
      29     2.643   0.036  1.3699E+02       4.920      27.1
      49     2.403   0.043  6.1476E+01       4.119      31.7
      34     2.210   0.051  2.7005E+01       3.296      34.8
      38     2.057   0.059  1.1303E+01       2.425      37.5
      38     1.935   0.067  3.5974E-01      -1.022      59.0
      46     1.830   0.075 -2.3089E+00    -999.999   -1000.0
      44     1.739   0.083 -6.3908E-01    -999.999   -1000.0
      42     1.669   0.090 -6.6894E-01    -999.999   -1000.0
      50     1.599   0.098  4.2991E-01      -0.844      39.4
      55     1.539   0.106 -5.2718E-01    -999.999   -1000.0
      63     1.483   0.114 -6.9276E-01    -999.999   -1000.0
      49     1.435   0.121  1.1968E+00       0.180      27.5
      45     1.392   0.129  2.3430E+00       0.851      23.3
      38     1.353   0.137 -2.4572E+00    -999.999   -1000.0
      18     1.316   0.144  2.9810E-01      -1.210      28.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.177        1.837        1.298        0.717
       3     4.705        1.369        0.869        0.410
       2     3.578        0.412        0.092        0.015
       3     2.999        1.960        0.997        0.372
       5     2.403        0.787        0.408        0.168
       1     2.210        0.353        0.073        0.011
       6     2.057        0.680        0.245        0.099
       6     1.935      597.067     7563.022    87578.590
       3     1.599       81.756      215.929     1027.597
       3     1.539       42.087      157.728      430.633
       3     1.435       20.078       16.199       43.822
       1     1.353      101.104      352.164      876.178
      75   overall       52.195      614.521     7049.233


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.177        0.539        0.389        0.223
      14     4.705        1.281        1.630        1.642
      19     3.578        1.038        1.074        1.006
      24     2.999        0.766        0.609        0.399
      28     2.643        1.353        1.888        2.441
      44     2.403        0.976        0.848        0.623
      33     2.210        0.892        0.674        0.452
      32     2.057        1.833        1.871        1.881
      47     1.599      244.468     -165.434    28960.971
      46     1.435       46.609       85.795     1055.977
      44     1.392       11.170        7.226       52.906
      18     1.316     1462.167    49209.975  2311236.394
     642   overall       83.193     1311.504    70112.956

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.177  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.38  0.73
       3     4.705  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.578  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.999  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       5     2.403  0.00  0.00  0.00  0.42  0.38  0.36  0.67  0.64  0.61  0.88
       1     2.210  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.057  2.01  1.45  1.20  1.06  0.96  0.89  0.84  1.06  1.01  0.98
       6     1.935  0.67  0.48  0.40  0.35  0.32  0.30  0.28  0.27  0.25  0.24
       3     1.599  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.539  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.59  1.52  1.47
       3     1.435  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       1     1.353  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      75   overall  0.59  0.42  0.35  0.34  0.31  0.29  0.29  0.34  0.32  0.37


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.177  3.15  2.76  1.93  2.12  1.78  1.55  1.59  1.45  1.35  1.27
      14     4.705  3.00  1.58  1.65  1.30  1.27  1.58  1.42  1.30  1.44  1.36
      19     3.578  1.11  1.16  1.02  0.80  0.94  1.05  1.25  1.24  1.15  1.08
      24     2.999  1.31  1.84  1.45  1.39  1.48  1.48  1.32  1.29  1.19  1.12
      28     2.643  1.13  1.18  1.24  1.08  1.00  0.95  1.14  1.10  1.08  1.13
      44     2.403  1.19  1.00  1.32  1.10  1.04  0.96  1.08  1.11  1.15  1.08
      33     2.210  1.91  1.17  0.94  0.74  0.69  0.74  0.84  0.77  0.71  0.86
      32     2.057  2.96  1.55  1.21  0.95  0.87  0.76  0.68  0.74  0.84  0.79
      47     1.599  5.14  2.70  1.89  1.48  1.24  1.08  0.97  0.93  0.86  0.81
      46     1.435  5.25  2.88  2.01  1.58  1.33  1.16  1.04  0.95  0.88  0.83
      44     1.392  4.06  2.13  1.49  1.17  0.98  0.86  0.77  0.70  0.65  0.61
      18     1.316  7.01  3.68  2.57  2.02  1.69  1.48  1.32  1.21  1.12  1.05
     642   overall  1.80  1.09  0.86  0.70  0.63  0.59  0.58  0.55  0.54  0.52


 List of     92 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30  175.71  0.5238E+02  0.1513E+02 "alien"
   40   29    3    1.56   71.17 -0.3752E+02  0.1659E+02 "alien"
   45    9   18    1.33   54.24  0.1617E+02  0.1075E+02 "alien"
   30   21   23    1.32   53.74  0.1602E+02  0.1316E+02 "alien"
   31   27   15    1.52   51.35 -0.2707E+02  0.1409E+02 "alien"
   34   22    1    1.91   49.28  0.1773E+02  0.9258E+01 "alien"
   47    1    6    1.60   47.35  0.2036E+02  0.1205E+02 "alien"
   51    8   15    1.30   44.98  0.1341E+02  0.1435E+02 "alien"
   39    2    3    1.96   43.53  0.1566E+02  0.1147E+02 "alien"
   25   18   12    1.98   41.86  0.1506E+02  0.1251E+02 "alien"
   37    4   15    1.62   39.11  0.1681E+02  0.1263E+02 "alien"
   39   28    2    1.61   39.05  0.1679E+02  0.8728E+01 "alien"
   31   26   14    1.58   33.56  0.1443E+02  0.8573E+01 "alien"
   40    2    1    1.93   33.30  0.1198E+02  0.1084E+02 "alien"
   38    1    7    1.91   31.02  0.1116E+02  0.8466E+01 "alien"
   35   27   10    1.60   30.64  0.1317E+02  0.8490E+01 "alien"
   51    2    2    1.51   29.99 -0.1581E+02  0.1616E+02 "alien"
   43    3   12    1.57   29.81  0.1282E+02  0.8337E+01 "alien"
   32   23   10    1.75   29.75 -0.1901E+02  0.1391E+02 "alien"
   36    3   14    1.70   28.40 -0.1900E+02  0.1268E+02 "alien"
   29    2   14    1.92   28.36  0.1020E+02  0.8180E+01 "alien"
   42   32    2    1.46   26.35 -0.1825E+02  0.9371E+01 "alien"
   28   23   14    1.69   25.98 -0.1738E+02  0.7427E+01 "alien"
   45    2   10    1.57   25.78  0.1108E+02  0.9842E+01 "alien"
   41   33    7    1.42   25.41  0.3041E+02  0.1441E+02 "alien"
   33   32   18    1.33   23.45  0.6992E+01  0.1313E+02 "alien"
   28   20   11    1.90   23.02  0.8281E+01  0.7857E+01 "alien"
   17   16   16    1.96   22.96  0.8261E+01  0.6426E+01 "alien"
   41   32    5    1.46   22.86 -0.1584E+02  0.9179E+01 "alien"
   15   12   18    1.90   22.53  0.8105E+01  0.5715E+01 "alien"
   35   28   11    1.55   21.86 -0.1152E+02  0.8669E+01 "alien"
   40   29    2    1.56   21.79 -0.1149E+02  0.9167E+01 "alien"
   53    3    4    1.44   21.62  0.2588E+02  0.1001E+02 "alien"
   36    2   13    1.74   21.42 -0.1369E+02  0.7359E+01 "alien"
   30   28   16    1.49   20.64 -0.1430E+02  0.1381E+02 "alien"
   19   18   18    1.74   20.16 -0.1288E+02  0.6211E+01 "alien"
   48    3    9    1.51   20.13 -0.1395E+02  0.1514E+02 "alien"
   40    2    0    1.93   20.06  0.3022E+02  0.1472E+02 "alien"
   47    3   10    1.51   19.98 -0.1053E+02  0.8813E+01 "alien"
   33    6   18    1.58   19.88  0.8546E+01  0.7214E+01 "alien"
   38    4   15    1.60   19.74  0.8487E+01  0.7740E+01 "alien"
   48    1    3    1.60   19.43  0.8351E+01  0.1039E+02 "alien"
   27   10   22    1.47   19.18 -0.1328E+02  0.1107E+02 "alien"
   52    7   13    1.32   18.71  0.5579E+01  0.1031E+02 "alien"
   52    2    3    1.48   18.46 -0.1278E+02  0.9255E+01 "alien"
   39   28    1    1.61   18.30  0.7867E+01  0.8765E+01 "alien"
   39   32    9    1.45   17.62  0.2108E+02  0.1682E+02 "alien"
   35   25    7    1.71   17.38 -0.1111E+02  0.8266E+01 "alien"
   51    3    7    1.46   16.92 -0.1172E+02  0.8953E+01 "alien"
   30    9   21    1.48   16.89 -0.1170E+02  0.7532E+01 "alien"
   22   14   22    1.52   16.51 -0.8702E+01  0.8992E+01 "alien"
   40   30    3    1.54   16.28 -0.8583E+01  0.9885E+01 "alien"
   48    2    7    1.55   16.07 -0.8474E+01  0.8872E+01 "alien"
   34   22    0    1.91   15.58  0.1657E+02  0.1362E+02 "alien"
   36   26    7    1.66   15.23 -0.1019E+02  0.7942E+01 "alien"
   38   28    6    1.59   15.01  0.6451E+01  0.1192E+02 "alien"
   25   10   22    1.50   14.97 -0.1037E+02  0.9072E+01 "alien"
   38    6   17    1.51   14.74 -0.1021E+02  0.9150E+01 "alien"
   34   26   10    1.64   14.68 -0.9819E+01  0.8868E+01 "alien"
   33    3   15    1.73   14.62 -0.9345E+01  0.7467E+01 "alien"
   53    3    6    1.42   14.54  0.1740E+02  0.1175E+02 "alien"
   21   20   17    1.73   14.43 -0.9222E+01  0.7187E+01 "alien"
   34   28   12    1.54   14.34 -0.7561E+01  0.1504E+02 "alien"
   44    1    1    1.76   14.28 -0.9124E+01  0.9758E+01 "alien"
   27   22   14    1.73   14.25 -0.9110E+01  0.7291E+01 "alien"
   36   27    8    1.62   14.10  0.6061E+01  0.9653E+01 "alien"
   41   31    0    1.51   14.00 -0.2555E+02  0.1329E+02 "alien"
   52    3    7    1.44   13.97  0.1671E+02  0.1079E+02 "alien"
   24   10   22    1.52   12.86 -0.6781E+01  0.8741E+01 "alien"
   43    2   11    1.60   12.84  0.5521E+01  0.1005E+02 "alien"
   37   27    6    1.64   12.64 -0.8452E+01  0.1207E+02 "alien"
   41   31    2    1.50   12.55 -0.8695E+01  0.8866E+01 "alien"
   46    3   11    1.52   12.33 -0.6503E+01  0.8800E+01 "alien"
   39   29    5    1.56   11.94 -0.6293E+01  0.8481E+01 "alien"
   34    3   15    1.71   11.79 -0.7537E+01  0.1278E+02 "alien"
   53    2    0    1.46   11.59 -0.1794E+02  0.1527E+02 "alien"
   34   29   13    1.50   11.45 -0.7931E+01  0.1028E+02 "alien"
   52    2    2    1.48   11.05 -0.7658E+01  0.8933E+01 "alien"
   27   24   16    1.61   10.99  0.4725E+01  0.1315E+02 "alien"
   43    0    7    1.71   10.88 -0.1484E+02  0.1018E+02 "alien"
   23   23   18    1.60   10.86  0.9950E+01  0.1091E+02 "alien"
   48    1    4    1.59   10.63  0.4571E+01  0.8576E+01 "alien"
   54    4    7    1.39   10.60  0.2484E+02  0.9790E+01 "alien"
   38   28    5    1.61   10.57  0.4546E+01  0.1075E+02 "alien"
   42    2   11    1.63   10.56 -0.7067E+01  0.8127E+01 "alien"
   46    0    4    1.66   10.51 -0.1464E+02  0.9280E+01 "alien"
   39    7   18    1.45   10.31  0.1234E+02  0.8642E+01 "alien"
   47    1    5    1.61   10.26  0.4412E+01  0.8855E+01 "alien"
   40    1   10    1.73   10.26 -0.6558E+01  0.7859E+01 "alien"
   23   13   23    1.46   10.22 -0.7079E+01  0.9218E+01 "alien"
   38   26    0    1.68   10.10 -0.1371E+02  0.1465E+02 "alien"
   24   14   23    1.44   10.03  0.1200E+02  0.9942E+01 "alien"


 List of     92 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   51    8   15    1.30   44.98  0.1341E+02  0.1435E+02
   53    7   13    1.30  175.71  0.5238E+02  0.1513E+02
   30   21   23    1.32   53.74  0.1602E+02  0.1316E+02
   52    7   13    1.32   18.71  0.5579E+01  0.1031E+02
   33   32   18    1.33   23.45  0.6992E+01  0.1313E+02
   45    9   18    1.33   54.24  0.1617E+02  0.1075E+02
   54    4    7    1.39   10.60  0.2484E+02  0.9790E+01
   41   33    7    1.42   25.41  0.3041E+02  0.1441E+02
   53    3    6    1.42   14.54  0.1740E+02  0.1175E+02
   52    3    7    1.44   13.97  0.1671E+02  0.1079E+02
   24   14   23    1.44   10.03  0.1200E+02  0.9942E+01
   53    3    4    1.44   21.62  0.2588E+02  0.1001E+02
   39   32    9    1.45   17.62  0.2108E+02  0.1682E+02
   39    7   18    1.45   10.31  0.1234E+02  0.8642E+01
   53    2    0    1.46   11.59 -0.1794E+02  0.1527E+02
   23   13   23    1.46   10.22 -0.7079E+01  0.9218E+01
   51    3    7    1.46   16.92 -0.1172E+02  0.8953E+01
   41   32    5    1.46   22.86 -0.1584E+02  0.9179E+01
   42   32    2    1.46   26.35 -0.1825E+02  0.9371E+01
   27   10   22    1.47   19.18 -0.1328E+02  0.1107E+02
   52    2    3    1.48   18.46 -0.1278E+02  0.9255E+01
   30    9   21    1.48   16.89 -0.1170E+02  0.7532E+01
   52    2    2    1.48   11.05 -0.7658E+01  0.8933E+01
   30   28   16    1.49   20.64 -0.1430E+02  0.1381E+02
   34   29   13    1.50   11.45 -0.7931E+01  0.1028E+02
   41   31    2    1.50   12.55 -0.8695E+01  0.8866E+01
   25   10   22    1.50   14.97 -0.1037E+02  0.9072E+01
   41   31    0    1.51   14.00 -0.2555E+02  0.1329E+02
   48    3    9    1.51   20.13 -0.1395E+02  0.1514E+02
   38    6   17    1.51   14.74 -0.1021E+02  0.9150E+01
   51    2    2    1.51   29.99 -0.1581E+02  0.1616E+02
   47    3   10    1.51   19.98 -0.1053E+02  0.8813E+01
   22   14   22    1.52   16.51 -0.8702E+01  0.8992E+01
   24   10   22    1.52   12.86 -0.6781E+01  0.8741E+01
   46    3   11    1.52   12.33 -0.6503E+01  0.8800E+01
   31   27   15    1.52   51.35 -0.2707E+02  0.1409E+02
   40   30    3    1.54   16.28 -0.8583E+01  0.9885E+01
   34   28   12    1.54   14.34 -0.7561E+01  0.1504E+02
   48    2    7    1.55   16.07 -0.8474E+01  0.8872E+01
   35   28   11    1.55   21.86 -0.1152E+02  0.8669E+01
   40   29    3    1.56   71.17 -0.3752E+02  0.1659E+02
   39   29    5    1.56   11.94 -0.6293E+01  0.8481E+01
   40   29    2    1.56   21.79 -0.1149E+02  0.9167E+01
   43    3   12    1.57   29.81  0.1282E+02  0.8337E+01
   45    2   10    1.57   25.78  0.1108E+02  0.9842E+01
   31   26   14    1.58   33.56  0.1443E+02  0.8573E+01
   33    6   18    1.58   19.88  0.8546E+01  0.7214E+01
   38   28    6    1.59   15.01  0.6451E+01  0.1192E+02
   48    1    4    1.59   10.63  0.4571E+01  0.8576E+01
   38    4   15    1.60   19.74  0.8487E+01  0.7740E+01
   47    1    6    1.60   47.35  0.2036E+02  0.1205E+02
   35   27   10    1.60   30.64  0.1317E+02  0.8490E+01
   23   23   18    1.60   10.86  0.9950E+01  0.1091E+02
   48    1    3    1.60   19.43  0.8351E+01  0.1039E+02
   43    2   11    1.60   12.84  0.5521E+01  0.1005E+02
   38   28    5    1.61   10.57  0.4546E+01  0.1075E+02
   27   24   16    1.61   10.99  0.4725E+01  0.1315E+02
   39   28    2    1.61   39.05  0.1679E+02  0.8728E+01
   47    1    5    1.61   10.26  0.4412E+01  0.8855E+01
   39   28    1    1.61   18.30  0.7867E+01  0.8765E+01
   37    4   15    1.62   39.11  0.1681E+02  0.1263E+02
   36   27    8    1.62   14.10  0.6061E+01  0.9653E+01
   42    2   11    1.63   10.56 -0.7067E+01  0.8127E+01
   37   27    6    1.64   12.64 -0.8452E+01  0.1207E+02
   34   26   10    1.64   14.68 -0.9819E+01  0.8868E+01
   46    0    4    1.66   10.51 -0.1464E+02  0.9280E+01
   36   26    7    1.66   15.23 -0.1019E+02  0.7942E+01
   38   26    0    1.68   10.10 -0.1371E+02  0.1465E+02
   28   23   14    1.69   25.98 -0.1738E+02  0.7427E+01
   36    3   14    1.70   28.40 -0.1900E+02  0.1268E+02
   34    3   15    1.71   11.79 -0.7537E+01  0.1278E+02
   35   25    7    1.71   17.38 -0.1111E+02  0.8266E+01
   43    0    7    1.71   10.88 -0.1484E+02  0.1018E+02
   40    1   10    1.73   10.26 -0.6558E+01  0.7859E+01
   33    3   15    1.73   14.62 -0.9345E+01  0.7467E+01
   27   22   14    1.73   14.25 -0.9110E+01  0.7291E+01
   21   20   17    1.73   14.43 -0.9222E+01  0.7187E+01
   19   18   18    1.74   20.16 -0.1288E+02  0.6211E+01
   36    2   13    1.74   21.42 -0.1369E+02  0.7359E+01
   32   23   10    1.75   29.75 -0.1901E+02  0.1391E+02
   44    1    1    1.76   14.28 -0.9124E+01  0.9758E+01
   28   20   11    1.90   23.02  0.8281E+01  0.7857E+01
   15   12   18    1.90   22.53  0.8105E+01  0.5715E+01
   34   22    1    1.91   49.28  0.1773E+02  0.9258E+01
   38    1    7    1.91   31.02  0.1116E+02  0.8466E+01
   34   22    0    1.91   15.58  0.1657E+02  0.1362E+02
   29    2   14    1.92   28.36  0.1020E+02  0.8180E+01
   40    2    1    1.93   33.30  0.1198E+02  0.1084E+02
   40    2    0    1.93   20.06  0.3022E+02  0.1472E+02
   17   16   16    1.96   22.96  0.8261E+01  0.6426E+01
   39    2    3    1.96   43.53  0.1566E+02  0.1147E+02
   25   18   12    1.98   41.86  0.1506E+02  0.1251E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        37
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1266
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1167
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       99


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
