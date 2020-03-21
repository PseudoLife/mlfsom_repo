
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.620    77.620    38.855  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      43      43
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      43
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001021  0.000689
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
 REFLECTING_RANGE_E.S.D.=   0.201
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

    1108 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      43 ...      43
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1108 REFLECTIONS ACCEPTED



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

 *  44        aP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.6   38.9   77.6  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  33        mP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  32        oP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  23        oC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.6   77.6   38.9  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     160.0   38.9   77.6  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     160.0   38.9   77.6  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.9  160.0   77.6  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.9  160.0   77.6  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.9  160.0   77.6  90.0  90.0  75.9    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.9  160.0   77.6  90.0  90.0  75.9    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     160.0   38.9  109.8  90.0 133.3  75.9   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.9  109.8  116.4  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.9  160.0  160.0  86.6 104.1 104.1    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.8  116.4   38.9  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        997.8     173.6   77.6   38.9  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        997.8      77.6  173.6   38.9  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        997.8      77.6  173.6   38.9  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        997.8      77.6   77.6   38.9  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.4  116.4  116.4  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.8   86.8   86.8  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.8   86.8  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.9   86.8  248.5 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.8   86.8   77.6  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.9   77.6  177.9  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.8   86.8  177.9 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.8   86.8   86.8  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.9  109.8   38.9  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.9   77.6  177.9 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.8  177.9   77.6 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    4.5      13    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.981   -23.492     0.033
 UNIT_CELL_B-AXIS=    23.492    73.981     0.017
 UNIT_CELL_C-AXIS=    -0.018    -0.003    38.854
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       109 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.28
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.627    77.627    38.850  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  3.0E-01 3.0E-01 1.2E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012882  0.012882  0.025740  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.984   -23.501     0.022
 COORDINATES OF UNIT CELL B-AXIS    23.501    73.984     0.000
 COORDINATES OF UNIT CELL C-AXIS    -0.011     0.003    38.850
 CRYSTAL MOSAICITY (DEGREES)     0.201
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006369  0.002239
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.87
 DETECTOR ORIGIN (PIXELS) AT                     1260.61   1242.14
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.38
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
  0.15  0.25       134         99.302       1.012
  0.25  0.35       156        110.811       0.907
  0.35  0.45       162         40.885       2.457
  0.45  0.55       167        150.699       0.667
  0.55  0.65       168         73.543       1.366
  0.65  0.75       182        134.586       0.746
  0.75  0.85       118        100.561       0.999
  0.85  0.95         9         25.185       3.989
 
 -0.05  0.95      1108        100.458       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               149.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.878
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    42.1 XMAX=    42.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59182 NYBIN=    6
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        314


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               149.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.876
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    86.5 XMAX=  2451.3 NXBIN=    2
 YMIN=     5.8 YMAX=  2517.8 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        314


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               147.1
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.878
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    42.1 XMAX=    42.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1269    1262
 DETECTOR_SURFACE_POSITION=    1671    1688
 DETECTOR_SURFACE_POSITION=     867    1688
 DETECTOR_SURFACE_POSITION=     867     835
 DETECTOR_SURFACE_POSITION=    1671     835
 DETECTOR_SURFACE_POSITION=    2178    1662
 DETECTOR_SURFACE_POSITION=    1645    2227
 DETECTOR_SURFACE_POSITION=     892    2227
 DETECTOR_SURFACE_POSITION=     360    1662
 DETECTOR_SURFACE_POSITION=     360     862
 DETECTOR_SURFACE_POSITION=     892     296
 DETECTOR_SURFACE_POSITION=    1645     296
 DETECTOR_SURFACE_POSITION=    2178     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        314


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
 8.626E-01  8.975E-04   35.94


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

   24.548   2.890    18.45   0.99      2.25      2.04      80     110      12
    2.890   2.051     2.44   1.24     30.40     26.00     178     222       0
    2.051   1.677    -0.01   0.85   3102.95   3375.80     208     262       2
    1.677   1.453    -0.03   0.88   2414.26   2536.85     230     294       0
    1.453   1.300     0.07   1.12    606.79    559.46     142     201       4
   --------------------------------------------------------------------------
   24.548   1.300     2.36   1.00     10.23      9.87     838    1089      18
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      43      43)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.688  -0.7013E+01  0.1358E+02        -0.52           1 
   47    0    0     1.652   0.3618E+01  0.1290E+02         0.28           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     3.5% OF MEAN INTENSITY


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

     3.86          56      37      1285        2.9%       1.9%      1.7%       38   29.43      2.7%    99.9*     0    0.000       0
     2.74          73      46      2098        2.2%       6.1%      4.8%       54   13.97      8.6%    99.1*     0    0.000       0
     2.25         137      80      2665        3.0%      22.7%     20.3%      114    3.59     32.1%    90.5*     0    0.000       0
     1.95         121      73      3101        2.4%     132.9%    144.3%       96    0.51    188.0%    34.7      0    0.000       0
     1.74         136      82      3481        2.4%    1728.3%   1805.9%      108    0.02   2444.2%    -6.9      0    0.000       0
     1.59         169     102      3848        2.7%    2721.7%   2609.0%      134   -0.09   3849.1%    -5.6      0    0.000       0
     1.47         160      98      4143        2.4%    4445.2%   4526.3%      124   -0.02   6286.5%     3.8      0    0.000       0
     1.38         164     101      4428        2.3%     491.1%    478.7%      126    0.08    694.5%   -17.0      0    0.000       0
     1.30          71      50      4744        1.1%     -99.9%    -99.9%       42    0.00    -99.9%     7.2      0    0.000       0
    total        1087     669     29793        2.2%      10.2%      9.8%      836    3.07     14.4%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1108
 NUMBER OF REJECTED MISFITS                                18
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1089
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    670



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

   WILSON LINE (using all data) : A=   6.424 B=  30.835 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      15     7.471   0.004  1.7693E+03       7.478    -117.7
      18     4.755   0.011  1.4268E+03       7.263     -38.0
      20     3.571   0.020  5.0852E+02       6.232       4.9
      24     2.984   0.028  1.0606E+02       4.664      31.3
      31     2.641   0.036  8.5245E+01       4.446      27.6
      42     2.411   0.043  2.7865E+01       3.327      36.0
      34     2.214   0.051  1.1829E+01       2.471      38.8
      35     2.062   0.059  3.3363E+00       1.205      44.4
      37     1.939   0.066  1.3071E+00       0.268      46.3
      43     1.833   0.074  8.1772E-01      -0.201      44.5
      43     1.743   0.082 -1.3737E+00    -999.999   -1000.0
      42     1.672   0.089 -1.3658E+00    -999.999   -1000.0
      47     1.603   0.097  4.5326E-01      -0.791      37.1
      47     1.543   0.105 -1.2785E+00    -999.999   -1000.0
      56     1.486   0.113 -5.9949E-02    -999.999   -1000.0
      45     1.437   0.121  2.6305E+00       0.967      22.5
      41     1.395   0.128  1.5884E+00       0.463      23.2
      34     1.355   0.136  1.5677E-01      -1.853      30.4
      15     1.318   0.144 -1.5665E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.471        2.204        1.767        1.105
       3     4.755        1.486        0.989        0.489
       2     3.571        0.364        0.078        0.012
       3     2.984        2.558        1.538        0.690
       3     2.411        0.324        0.082        0.016
       1     2.214        0.021        0.001        0.000
       6     2.062       13.199       15.798       26.685
       6     1.939       24.522       68.389      142.030
       7     1.743       30.299      -21.266      172.814
       3     1.543       10.294       15.635       18.099
       1     1.395      410.033     2876.195    14410.869
       1     1.355     1469.781    19519.540   185164.829
      71   overall       32.355      321.343     2833.925


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.471        0.590        0.473        0.296
      15     4.755        1.346        1.809        1.929
      18     3.571        1.088        1.292        1.408
      21     2.984        0.652        0.523        0.362
      30     2.641        1.346        1.961        2.717
      39     2.411        1.014        0.919        0.756
      33     2.214        1.047        0.810        0.543
      29     2.062        5.741        6.885       14.669
      31     1.939       37.980      -22.473      485.342
      35     1.833       54.897       18.769     1797.648
      44     1.603      203.874      292.360    18091.487
      51     1.486       55.003       33.788     1467.925
      42     1.437        6.789        4.840       20.951
      40     1.395       17.804        2.569      203.200
     598   overall      221.233     2096.128   326104.011

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.471  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.73
       3     4.755  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.571  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.76  0.73
       3     2.984  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.411  0.00  0.00  0.00  0.00  0.00  0.59  1.12  1.06  1.01  0.98
       1     2.214  0.00  0.00  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       6     2.062  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       6     1.939  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.80  0.76  0.73
       7     1.743  2.88  2.07  1.71  1.51  1.37  1.27  1.20  1.14  1.09  1.05
       3     1.543  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       1     1.395  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     1.355  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      71   overall  0.74  0.53  0.47  0.42  0.38  0.38  0.38  0.38  0.41  0.43


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.471  2.87  3.51  2.46  2.21  1.85  1.61  1.63  1.49  1.38  1.29
      15     4.755  3.50  1.84  1.80  1.42  1.19  1.18  1.46  1.33  1.35  1.37
      18     3.571  1.17  0.61  0.64  0.67  0.71  0.74  1.21  1.21  1.22  1.14
      21     2.984  1.50  1.31  1.47  1.16  1.45  1.27  1.32  1.38  1.28  1.21
      30     2.641  0.70  1.29  0.90  0.91  1.02  1.11  1.06  1.03  1.12  1.11
      39     2.411  1.35  0.85  0.99  0.93  1.04  0.97  0.92  0.88  0.91  0.97
      33     2.214  1.91  1.17  1.05  0.92  0.85  0.74  0.90  0.88  0.87  0.81
      29     2.062  4.71  2.47  1.86  1.46  1.23  1.07  0.96  1.00  0.99  0.93
      31     1.939  4.75  2.49  1.74  1.37  1.15  1.00  0.90  0.88  0.82  0.82
      35     1.833  5.10  2.68  1.87  1.47  1.23  1.14  1.02  0.93  0.87  0.81
      44     1.603  5.02  2.63  1.84  1.45  1.21  1.06  0.99  0.91  0.84  0.79
      51     1.486  5.98  3.14  2.19  1.72  1.45  1.26  1.13  1.03  0.96  0.90
      42     1.437  3.00  1.71  1.19  1.01  0.85  0.74  0.71  0.65  0.60  0.60
      40     1.395  5.52  2.90  2.03  1.59  1.33  1.16  1.04  1.00  0.93  0.87
     598   overall  2.69  1.54  1.16  0.95  0.85  0.77  0.75  0.72  0.69  0.67


 List of     79 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   52    2    3    1.48  416.97 -0.2500E+02  0.9873E+01 "alien"
   38    7   18    1.47  363.64 -0.2180E+02  0.1278E+02 "alien"
   41   32    5    1.47  302.70 -0.1815E+02  0.9901E+01 "alien"
   42   32    2    1.47  292.99 -0.1756E+02  0.1024E+02 "alien"
   51    6   13    1.35  193.53  0.3034E+02  0.1851E+02 "alien"
   51    2    4    1.50  189.85 -0.1138E+02  0.9263E+01 "alien"
   34   33   16    1.36  176.89  0.2773E+02  0.1708E+02 "alien"
   25   20   21    1.47  171.09 -0.1026E+02  0.1163E+02 "alien"
   34    7   19    1.51  138.85 -0.8324E+01  0.1231E+02 "alien"
   41   10   20    1.34  129.28  0.2027E+02  0.1089E+02 "alien"
   45    8   17    1.36  119.95  0.1880E+02  0.1684E+02 "alien"
   41   31    2    1.51  118.45 -0.7101E+01  0.9620E+01 "alien"
   41   31    3    1.50  116.77 -0.7000E+01  0.9693E+01 "alien"
   48    3   10    1.49  113.65 -0.6813E+01  0.1527E+02 "alien"
   51    3    7    1.47  103.91 -0.6230E+01  0.9631E+01 "alien"
   52    2    2    1.49  102.64 -0.6153E+01  0.9407E+01 "alien"
   43    5   15    1.47   99.03 -0.5937E+01  0.1609E+02 "alien"
   38   34   13    1.36   98.25  0.1540E+02  0.1056E+02 "alien"
   32   12   23    1.36   86.03  0.1349E+02  0.1351E+02 "alien"
   28   19   23    1.36   72.05  0.1130E+02  0.1385E+02 "alien"
   41   31    1    1.51   66.43 -0.3983E+01  0.1048E+02 "alien"
   38   31   10    1.47   65.49 -0.3926E+01  0.9872E+01 "alien"
   42    5   15    1.50   63.18 -0.3788E+01  0.9347E+01 "alien"
   29   22   22    1.36   61.18  0.9592E+01  0.1370E+02 "alien"
   54    5    9    1.36   59.87  0.9386E+01  0.1049E+02 "alien"
   34   11   22    1.37   58.45  0.9163E+01  0.9610E+01 "alien"
   47    1    5    1.62   44.40  0.2012E+02  0.1013E+02 "alien"
   39   28    2    1.61   41.96  0.1902E+02  0.9223E+01 "alien"
   50    2    6    1.51   40.97 -0.2456E+01  0.9720E+01 "alien"
   39    6   17    1.49   36.04 -0.2161E+01  0.8954E+01 "alien"
   41   31    0    1.51   34.75 -0.3285E+02  0.1512E+02 "alien"
   44    2   10    1.61   31.65  0.1435E+02  0.9039E+01 "alien"
   53    2    0    1.46   31.29 -0.2123E+02  0.1705E+02 "alien"
   50    3    8    1.48   30.79 -0.1846E+01  0.9644E+01 "alien"
   35    5   17    1.58   29.47  0.1336E+02  0.1368E+02 "alien"
   26   26   18    1.51   29.13 -0.1616E+02  0.9073E+01 "alien"
   43    3   12    1.57   27.42  0.1243E+02  0.9527E+01 "alien"
   52    2    0    1.49   26.17 -0.1109E+02  0.1366E+02 "alien"
   28   24   15    1.63   25.72  0.1166E+02  0.8598E+01 "alien"
   38    4   15    1.60   25.06  0.1136E+02  0.8789E+01 "alien"
   31   29   19    1.36   24.98  0.3917E+01  0.1003E+02 "alien"
   21   19   16    1.82   24.91  0.2037E+02  0.1228E+02 "alien"
   45   10   19    1.30   23.97 -0.3755E+02  0.1617E+02 "alien"
   35   27   10    1.60   22.92  0.1039E+02  0.9182E+01 "alien"
   46    1    6    1.63   20.06  0.9091E+01  0.1355E+02 "alien"
   53    6   11    1.35   19.46  0.3051E+01  0.1166E+02 "alien"
   31   27   15    1.53   19.10 -0.2442E+02  0.1436E+02 "alien"
   39    4   15    1.57   18.29  0.8289E+01  0.1330E+02 "alien"
   33   33   17    1.34   18.28  0.1041E+02  0.1588E+02 "alien"
   53    5   10    1.37   18.20  0.2854E+01  0.1055E+02 "alien"
   39   28    1    1.62   18.03  0.8173E+01  0.9484E+01 "alien"
   35   23    2    1.85   17.91  0.1465E+02  0.9269E+01 "alien"
   41    1    5    1.84   17.73  0.1449E+02  0.1313E+02 "alien"
   33    6   18    1.58   17.10  0.7752E+01  0.8717E+01 "alien"
   39   33   11    1.40   16.06  0.2551E+02  0.1111E+02 "alien"
   36   27    9    1.60   15.77  0.7149E+01  0.1031E+02 "alien"
   33   25   11    1.66   15.60 -0.2131E+02  0.1363E+02 "alien"
   21   17   14    2.00   14.90  0.1947E+02  0.1319E+02 "alien"
   29   21   11    1.85   14.52  0.1187E+02  0.8396E+01 "alien"
   55    3    0    1.41   13.64  0.5571E+02  0.1454E+02 "alien"
   31    8   20    1.52   13.31 -0.1701E+02  0.1150E+02 "alien"
   34   24    7    1.77   13.21 -0.1815E+02  0.1589E+02 "alien"
   28   21   12    1.83   13.07  0.1069E+02  0.8215E+01 "alien"
   41    1    4    1.86   12.32  0.1008E+02  0.9262E+01 "alien"
   34   22    1    1.91   12.18  0.1592E+02  0.1019E+02 "alien"
   50    2    5    1.52   11.99 -0.1533E+02  0.9727E+01 "alien"
   39    2    1    1.99   11.94  0.1560E+02  0.9282E+01 "alien"
   25   19   13    1.90   11.93  0.1560E+02  0.8778E+01 "alien"
   21   19   19    1.64   11.83 -0.1615E+02  0.1077E+02 "alien"
   38    6   17    1.51   11.47 -0.6878E+00  0.1296E+02 "alien"
   40   29    2    1.57   11.28 -0.1442E+02  0.9959E+01 "alien"
   35    2   13    1.78   11.12 -0.1528E+02  0.1205E+02 "alien"
   32   23    9    1.79   11.08  0.9059E+01  0.1563E+02 "alien"
   32   23   10    1.76   11.07 -0.1521E+02  0.1451E+02 "alien"
   36   32   14    1.39   10.78  0.1713E+02  0.1689E+02 "alien"
   39    2   12    1.69   10.63 -0.1452E+02  0.8766E+01 "alien"
   37   26    3    1.70   10.29 -0.1405E+02  0.1626E+02 "alien"
   46    2    9    1.57   10.22 -0.1306E+02  0.8832E+01 "alien"
   39    2    2    1.98   10.08  0.1318E+02  0.9579E+01 "alien"


 List of     79 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   45   10   19    1.30   23.97 -0.3755E+02  0.1617E+02
   41   10   20    1.34  129.28  0.2027E+02  0.1089E+02
   33   33   17    1.34   18.28  0.1041E+02  0.1588E+02
   53    6   11    1.35   19.46  0.3051E+01  0.1166E+02
   51    6   13    1.35  193.53  0.3034E+02  0.1851E+02
   32   12   23    1.36   86.03  0.1349E+02  0.1351E+02
   38   34   13    1.36   98.25  0.1540E+02  0.1056E+02
   34   33   16    1.36  176.89  0.2773E+02  0.1708E+02
   54    5    9    1.36   59.87  0.9386E+01  0.1049E+02
   29   22   22    1.36   61.18  0.9592E+01  0.1370E+02
   28   19   23    1.36   72.05  0.1130E+02  0.1385E+02
   31   29   19    1.36   24.98  0.3917E+01  0.1003E+02
   45    8   17    1.36  119.95  0.1880E+02  0.1684E+02
   53    5   10    1.37   18.20  0.2854E+01  0.1055E+02
   34   11   22    1.37   58.45  0.9163E+01  0.9610E+01
   36   32   14    1.39   10.78  0.1713E+02  0.1689E+02
   39   33   11    1.40   16.06  0.2551E+02  0.1111E+02
   55    3    0    1.41   13.64  0.5571E+02  0.1454E+02
   53    2    0    1.46   31.29 -0.2123E+02  0.1705E+02
   51    3    7    1.47  103.91 -0.6230E+01  0.9631E+01
   41   32    5    1.47  302.70 -0.1815E+02  0.9901E+01
   38   31   10    1.47   65.49 -0.3926E+01  0.9872E+01
   42   32    2    1.47  292.99 -0.1756E+02  0.1024E+02
   38    7   18    1.47  363.64 -0.2180E+02  0.1278E+02
   25   20   21    1.47  171.09 -0.1026E+02  0.1163E+02
   43    5   15    1.47   99.03 -0.5937E+01  0.1609E+02
   50    3    8    1.48   30.79 -0.1846E+01  0.9644E+01
   52    2    3    1.48  416.97 -0.2500E+02  0.9873E+01
   52    2    2    1.49  102.64 -0.6153E+01  0.9407E+01
   48    3   10    1.49  113.65 -0.6813E+01  0.1527E+02
   39    6   17    1.49   36.04 -0.2161E+01  0.8954E+01
   52    2    0    1.49   26.17 -0.1109E+02  0.1366E+02
   42    5   15    1.50   63.18 -0.3788E+01  0.9347E+01
   41   31    3    1.50  116.77 -0.7000E+01  0.9693E+01
   51    2    4    1.50  189.85 -0.1138E+02  0.9263E+01
   41   31    2    1.51  118.45 -0.7101E+01  0.9620E+01
   50    2    6    1.51   40.97 -0.2456E+01  0.9720E+01
   34    7   19    1.51  138.85 -0.8324E+01  0.1231E+02
   41   31    1    1.51   66.43 -0.3983E+01  0.1048E+02
   26   26   18    1.51   29.13 -0.1616E+02  0.9073E+01
   41   31    0    1.51   34.75 -0.3285E+02  0.1512E+02
   38    6   17    1.51   11.47 -0.6878E+00  0.1296E+02
   31    8   20    1.52   13.31 -0.1701E+02  0.1150E+02
   50    2    5    1.52   11.99 -0.1533E+02  0.9727E+01
   31   27   15    1.53   19.10 -0.2442E+02  0.1436E+02
   40   29    2    1.57   11.28 -0.1442E+02  0.9959E+01
   46    2    9    1.57   10.22 -0.1306E+02  0.8832E+01
   39    4   15    1.57   18.29  0.8289E+01  0.1330E+02
   43    3   12    1.57   27.42  0.1243E+02  0.9527E+01
   33    6   18    1.58   17.10  0.7752E+01  0.8717E+01
   35    5   17    1.58   29.47  0.1336E+02  0.1368E+02
   38    4   15    1.60   25.06  0.1136E+02  0.8789E+01
   35   27   10    1.60   22.92  0.1039E+02  0.9182E+01
   36   27    9    1.60   15.77  0.7149E+01  0.1031E+02
   44    2   10    1.61   31.65  0.1435E+02  0.9039E+01
   39   28    2    1.61   41.96  0.1902E+02  0.9223E+01
   47    1    5    1.62   44.40  0.2012E+02  0.1013E+02
   39   28    1    1.62   18.03  0.8173E+01  0.9484E+01
   46    1    6    1.63   20.06  0.9091E+01  0.1355E+02
   28   24   15    1.63   25.72  0.1166E+02  0.8598E+01
   21   19   19    1.64   11.83 -0.1615E+02  0.1077E+02
   33   25   11    1.66   15.60 -0.2131E+02  0.1363E+02
   39    2   12    1.69   10.63 -0.1452E+02  0.8766E+01
   37   26    3    1.70   10.29 -0.1405E+02  0.1626E+02
   32   23   10    1.76   11.07 -0.1521E+02  0.1451E+02
   34   24    7    1.77   13.21 -0.1815E+02  0.1589E+02
   35    2   13    1.78   11.12 -0.1528E+02  0.1205E+02
   32   23    9    1.79   11.08  0.9059E+01  0.1563E+02
   21   19   16    1.82   24.91  0.2037E+02  0.1228E+02
   28   21   12    1.83   13.07  0.1069E+02  0.8215E+01
   41    1    5    1.84   17.73  0.1449E+02  0.1313E+02
   35   23    2    1.85   17.91  0.1465E+02  0.9269E+01
   29   21   11    1.85   14.52  0.1187E+02  0.8396E+01
   41    1    4    1.86   12.32  0.1008E+02  0.9262E+01
   25   19   13    1.90   11.93  0.1560E+02  0.8778E+01
   34   22    1    1.91   12.18  0.1592E+02  0.1019E+02
   39    2    2    1.98   10.08  0.1318E+02  0.9579E+01
   39    2    1    1.99   11.94  0.1560E+02  0.9282E+01
   21   17   14    2.00   14.90  0.1947E+02  0.1319E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        45
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1108
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1018
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       90


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
