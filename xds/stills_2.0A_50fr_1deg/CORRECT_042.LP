
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.612    77.612    38.847  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      42      42
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      42
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001005  0.000646
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
 REFLECTING_RANGE_E.S.D.=   0.192
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

    1126 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      42 ...      42
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1126 REFLECTIONS ACCEPTED



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

 *  44        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  1  0  0  0  0 -1  0  0
 *  35        mP          0.0      77.6   38.8   77.6  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  25        mC          0.0     109.8  109.8   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  34        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  33        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  23        oC          0.0     109.8  109.8   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  32        oP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  20        mC          0.0     109.8  109.8   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.6   77.6   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     160.0   38.8   77.6  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     160.0   38.8   77.6  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  160.0   77.6  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  160.0   77.6  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  160.0   77.6  90.0  90.0  75.9    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  160.0   77.6  90.0  90.0  75.9    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     160.0   38.8  109.8  90.0 133.3  75.9   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.8  116.4  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  160.0  160.0  86.6 104.1 104.1    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.8  116.4   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        997.9     173.5   77.6   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        997.9      77.6  173.5   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        997.9      77.6  173.5   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        997.9      77.6   77.6   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.4  116.4  116.4  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.8   86.8   86.8  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.8   86.8  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.8  248.5 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.8   86.8   77.6  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.6  177.8  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.8   86.8  177.8 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.8   86.8   86.8  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.8  109.8   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.6  177.8 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.8  177.8   77.6 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    4.5      15    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.975   -23.479     0.063
 UNIT_CELL_B-AXIS=    23.479    73.975     0.022
 UNIT_CELL_C-AXIS=    -0.033    -0.001    38.847
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       114 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.32
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.608    77.608    38.839  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.5E-01 2.5E-01 1.4E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012885  0.012885  0.025747  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.969   -23.486     0.135
 COORDINATES OF UNIT CELL B-AXIS    23.486    73.969     0.047
 COORDINATES OF UNIT CELL C-AXIS    -0.071    -0.002    38.839
 CRYSTAL MOSAICITY (DEGREES)     0.192
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999976  0.006492  0.002486
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.89   1232.86
 DETECTOR ORIGIN (PIXELS) AT                     1260.27   1242.31
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.33
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.61    77.61    38.84  90.000  90.000  90.000


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
  0.05  0.15        13         29.080       3.367
  0.15  0.25       137         88.483       1.106
  0.25  0.35       166        102.633       0.954
  0.35  0.45       161         35.833       2.732
  0.45  0.55       154        180.759       0.542
  0.55  0.65       174         70.239       1.394
  0.65  0.75       187        125.421       0.781
  0.75  0.85       127         87.139       1.123
  0.85  0.95         7         49.580       1.975
 
 -0.05  0.95      1126         97.897       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               150.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.016
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    41.1 XMAX=    41.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59213 NYBIN=    6
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        326


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               150.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.013
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    54.2 XMAX=  2450.7 NXBIN=    2
 YMIN=     5.9 YMAX=  2518.2 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        326


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               145.6
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.015
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    41.1 XMAX=    41.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1252    1262
 DETECTOR_SURFACE_POSITION=    1659    1689
 DETECTOR_SURFACE_POSITION=     845    1689
 DETECTOR_SURFACE_POSITION=     845     835
 DETECTOR_SURFACE_POSITION=    1659     835
 DETECTOR_SURFACE_POSITION=    2174    1662
 DETECTOR_SURFACE_POSITION=    1634    2228
 DETECTOR_SURFACE_POSITION=     871    2228
 DETECTOR_SURFACE_POSITION=     331    1662
 DETECTOR_SURFACE_POSITION=     331     862
 DETECTOR_SURFACE_POSITION=     871     296
 DETECTOR_SURFACE_POSITION=    1634     296
 DETECTOR_SURFACE_POSITION=    2174     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        326


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
 7.803E-01  9.953E-04   35.88


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

   24.542   2.890    18.07   1.79      2.52      2.01      80     100      22
    2.890   2.050     2.88   1.29     20.15     17.27     194     222       4
    2.050   1.676    -0.08   0.78    -99.90    -99.90     228     266       0
    1.676   1.452    -0.01   0.84    -99.90    -99.90     248     292       2
    1.452   1.300    -0.05   0.72    -99.90    -99.90     164     213       4
   --------------------------------------------------------------------------
   24.542   1.300     2.21   0.98     10.21     10.00     914    1093      32
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      42      42)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.687  -0.1199E+01  0.1321E+02        -0.09           1 
   47    0    0     1.651   0.1381E+02  0.1257E+02         1.10           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    13.7% OF MEAN INTENSITY


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

     3.86          55      34      1285        2.6%       2.2%      1.7%       42   30.59      3.1%    99.9*     0    0.000       0
     2.74          66      38      2098        1.8%       6.5%      4.7%       56   13.20      9.2%    99.3*     0    0.000       0
     2.24         136      75      2665        2.8%      19.2%     17.1%      122    4.02     27.2%    93.0*     0    0.000       0
     1.95         119      69      3101        2.2%     144.9%    135.8%      100    0.59    204.9%    33.8      0    0.000       0
     1.74         142      82      3479        2.4%     -99.9%    -99.9%      120    0.00    -99.9%    -6.5      0    0.000       0
     1.59         167      95      3850        2.5%    1724.8%   1875.7%      144    0.01   2439.2%    13.2      0    0.000       0
     1.47         159      91      4143        2.2%   19786.3%  22675.8%      136    0.02  27982.0%     2.0      0    0.000       0
     1.38         169      98      4428        2.2%     -99.9%    -99.9%      142    0.00    -99.9%    24.1      0    0.000       0
     1.30          76      52      4740        1.1%     -99.9%    -99.9%       48    0.00    -99.9%   -19.8      0    0.000       0
    total        1089     634     29789        2.1%      10.2%      9.9%      910    2.95     14.4%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1126
 NUMBER OF REJECTED MISFITS                                32
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1093
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    636



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

   WILSON LINE (using all data) : A=   6.952 B=  35.023 CORRELATION=  0.91
      #      RES      SS        <I>       log(<I>)       BO
      12     7.272   0.005  2.0710E+03       7.636     -72.3
      18     4.754   0.011  1.4461E+03       7.277     -14.7
      12     3.650   0.019  5.5102E+02       6.312      17.1
      25     2.989   0.028  1.1063E+02       4.706      40.1
      29     2.634   0.036  8.9419E+01       4.493      34.1
      39     2.408   0.043  2.4412E+01       3.195      43.6
      31     2.216   0.051  1.2984E+01       2.564      43.1
      35     2.062   0.059  3.6454E+00       1.293      48.1
      34     1.937   0.067 -1.4497E+00    -999.999   -1000.0
      44     1.832   0.074 -5.5986E-02    -999.999   -1000.0
      40     1.745   0.082  6.0341E-01      -0.505      45.4
      40     1.671   0.090 -1.3187E-01    -999.999   -1000.0
      44     1.603   0.097  3.5112E-01      -1.047      41.1
      44     1.543   0.105  2.1501E+00       0.766      29.4
      52     1.486   0.113 -1.2672E+00    -999.999   -1000.0
      45     1.437   0.121  2.0483E+00       0.717      25.7
      39     1.396   0.128 -1.1562E+00    -999.999   -1000.0
      34     1.355   0.136 -3.4800E-01    -999.999   -1000.0
      17     1.317   0.144 -3.2911E-01    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.272        1.583        1.066        0.561
       3     4.754        1.460        0.958        0.468
       1     3.650        0.417        0.093        0.015
       3     2.989        2.529        1.472        0.627
       1     2.634        0.017        0.001        0.000
       3     2.408        0.147        0.022        0.003
       1     2.216        0.389        0.084        0.013
       6     1.937      190.751     1319.199     7902.921
       9     1.832       32.580      102.222      410.122
      10     1.671       28.462      -26.481      135.355
       3     1.603      271.658      996.384     8498.336
       3     1.543       13.420        3.110       23.353
      70   overall       18.709       34.904      419.258


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

       8     7.272        0.558        0.391        0.210
      15     4.754        1.346        1.806        1.924
      11     3.650        1.188        1.508        1.639
      22     2.989        0.626        0.475        0.313
      28     2.634        1.258        1.688        2.161
      36     2.408        0.912        0.651        0.397
      30     2.216        1.440        1.590        1.512
      29     2.062        5.841        9.204       21.678
      33     1.745       75.212      233.307     3042.201
      41     1.603      434.995     4836.411   177184.259
      41     1.543        7.820       14.517       52.252
      42     1.437       16.661       10.459      107.757
     564   overall       38.421      367.853    13071.666

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.272  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.76  0.73
       3     4.754  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       1     3.650  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.989  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     2.634  0.00  0.00  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       3     2.408  0.00  0.00  0.00  0.70  0.64  0.59  1.12  1.06  1.52  1.47
       1     2.216  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     1.937  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       9     1.832  2.24  1.61  1.33  1.17  1.07  0.99  0.93  0.88  0.85  0.81
      10     1.671  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
       3     1.603  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       3     1.543  1.34  0.96  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      70   overall  0.98  0.70  0.65  0.60  0.55  0.51  0.50  0.50  0.52  0.50


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       8     7.272  3.94  3.45  2.41  2.27  1.91  1.66  1.49  1.36  1.26  1.19
      15     4.754  3.50  1.84  1.80  1.42  1.19  1.18  1.46  1.33  1.35  1.37
      11     3.650  0.96  0.50  0.35  1.10  1.16  1.01  1.08  1.16  1.23  1.15
      22     2.989  0.48  1.25  1.40  1.38  1.39  1.41  1.35  1.24  1.15  1.15
      28     2.634  0.75  1.18  0.96  0.97  1.00  1.11  1.14  1.04  1.02  1.07
      36     2.408  1.46  0.77  0.86  0.76  0.78  0.80  0.72  0.81  0.89  0.92
      30     2.216  3.50  1.84  1.41  1.42  1.19  1.18  1.13  1.09  1.01  0.95
      29     2.062  3.62  2.28  1.73  1.36  1.23  1.07  0.96  0.88  0.81  0.82
      33     1.745  5.09  2.84  1.99  1.56  1.31  1.14  1.02  0.94  0.87  0.81
      41     1.603  5.38  2.96  2.07  1.63  1.36  1.24  1.11  1.02  0.95  0.93
      41     1.543  3.84  2.29  1.69  1.41  1.24  1.14  1.11  1.02  0.95  0.89
      42     1.437  3.75  1.97  1.38  1.08  0.97  0.84  0.80  0.78  0.72  0.68
     564   overall  1.94  1.17  0.90  0.78  0.70  0.66  0.63  0.59  0.57  0.56


 List of     90 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   37    1   11    1.80  353.90 -0.1981E+02  0.8325E+01 "alien"
   43    1    1    1.80  343.70 -0.1924E+02  0.8071E+01 "alien"
   34   23    5    1.84  255.09 -0.1428E+02  0.1005E+02 "alien"
   33   23    8    1.79  220.40 -0.1234E+02  0.8019E+01 "alien"
   33   22    6    1.87  189.39 -0.1060E+02  0.8074E+01 "alien"
   42    1    4    1.81  182.91 -0.1024E+02  0.9076E+01 "alien"
   34   23    6    1.81  170.71 -0.9557E+01  0.8027E+01 "alien"
   37   27    6    1.64  142.93 -0.1885E+02  0.1309E+02 "alien"
   34   22    4    1.88  122.10 -0.6836E+01  0.9714E+01 "alien"
   35   23    1    1.85  117.13 -0.6558E+01  0.8351E+01 "alien"
   37   26    3    1.70  107.97 -0.1424E+02  0.1535E+02 "alien"
   37   26    5    1.68  107.25 -0.1414E+02  0.8559E+01 "alien"
   35    5   17    1.58  103.87  0.3647E+02  0.1410E+02 "alien"
   42    1    9    1.70  102.79 -0.1355E+02  0.8562E+01 "alien"
   28   21   12    1.83   98.24 -0.5500E+01  0.7607E+01 "alien"
   42    2   11    1.64   97.45 -0.1285E+02  0.8140E+01 "alien"
   20   18   16    1.86   97.23 -0.5443E+01  0.1010E+02 "alien"
   36   24    1    1.79   92.68 -0.5189E+01  0.8194E+01 "alien"
   40    1    6    1.86   77.40 -0.4333E+01  0.9155E+01 "alien"
   35   23    2    1.84   74.17 -0.4153E+01  0.8325E+01 "alien"
   36   27    8    1.63   73.92  0.2595E+02  0.1083E+02 "alien"
   46    9   18    1.31   72.90 -0.2399E+02  0.1470E+02 "alien"
   43    1    2    1.80   71.79 -0.4019E+01  0.8267E+01 "alien"
   32   25   12    1.65   68.63 -0.9050E+01  0.8332E+01 "alien"
   22   21   17    1.70   65.74 -0.8669E+01  0.9789E+01 "alien"
   31    3   15    1.80   65.07 -0.3643E+01  0.1084E+02 "alien"
   32   24   11    1.70   64.18 -0.8464E+01  0.8558E+01 "alien"
   31   26   21    1.33   61.99 -0.2040E+02  0.1350E+02 "alien"
   32   30   19    1.34   60.71 -0.2113E+02  0.1434E+02 "alien"
   26   20   13    1.85   59.95 -0.3356E+01  0.7392E+01 "alien"
   34   25    9    1.69   57.80 -0.7622E+01  0.7969E+01 "alien"
   45    8   17    1.36   54.00 -0.1879E+02  0.1596E+02 "alien"
   36   26    7    1.67   52.21 -0.6885E+01  0.8207E+01 "alien"
   35   25    8    1.69   48.33 -0.6373E+01  0.9905E+01 "alien"
   41   10   20    1.34   44.54 -0.1550E+02  0.9857E+01 "alien"
   21   19   16    1.82   44.35 -0.2483E+01  0.1076E+02 "alien"
   35   11   22    1.36   41.72 -0.1452E+02  0.1252E+02 "alien"
   33    4   16    1.68   38.23 -0.5041E+01  0.7392E+01 "alien"
   29   20   23    1.34   38.20 -0.1329E+02  0.1349E+02 "alien"
   45    2    9    1.60   36.51  0.1282E+02  0.1520E+02 "alien"
   48    7   15    1.36   36.04 -0.1254E+02  0.1008E+02 "alien"
   21   17   20    1.61   35.08  0.1232E+02  0.9843E+01 "alien"
   28   19   23    1.36   34.44 -0.1198E+02  0.8877E+01 "alien"
   39   28    2    1.61   33.94  0.1192E+02  0.8605E+01 "alien"
   37   25    2    1.73   33.21  0.2004E+02  0.9829E+01 "alien"
   40   11   21    1.32   32.54 -0.1071E+02  0.1414E+02 "alien"
   32   26   13    1.59   28.43  0.9983E+01  0.8720E+01 "alien"
   30   22   11    1.80   28.35 -0.1587E+01  0.7894E+01 "alien"
   43    1    0    1.80   27.24 -0.1187E+02  0.1162E+02 "alien"
   34   24    7    1.77   26.14  0.1577E+02  0.1541E+02 "alien"
   35   34   16    1.33   25.67 -0.8447E+01  0.1379E+02 "alien"
   47    0    1    1.65   23.34 -0.1697E+02  0.8609E+01 "alien"
   39   34   12    1.36   23.19 -0.8071E+01  0.1135E+02 "alien"
   35   23    3    1.83   23.03 -0.1289E+01  0.8446E+01 "alien"
   50    5   12    1.39   21.68 -0.2507E+02  0.1044E+02 "alien"
   29   22   22    1.36   21.39 -0.7444E+01  0.9084E+01 "alien"
   24   23   17    1.63   21.39  0.7509E+01  0.1236E+02 "alien"
   38   26    0    1.69   20.36 -0.1152E+02  0.1502E+02 "alien"
   33    1   12    1.90   20.28 -0.2939E+02  0.1206E+02 "alien"
   45    9   18    1.33   19.44 -0.6396E+01  0.1512E+02 "alien"
   36    0   10    1.88   19.38 -0.4300E+01  0.8333E+01 "alien"
   44    0    6    1.70   18.85 -0.9444E+01  0.9096E+01 "alien"
   36   25    4    1.74   18.05  0.1089E+02  0.1074E+02 "alien"
   52    5   10    1.39   17.38 -0.2009E+02  0.1095E+02 "alien"
   34    5   17    1.61   16.03  0.5627E+01  0.7970E+01 "alien"
   38    4   15    1.60   15.90  0.5583E+01  0.8017E+01 "alien"
   34   24    8    1.74   14.30  0.8629E+01  0.8628E+01 "alien"
   40   29    0    1.57   14.25  0.1343E+02  0.1274E+02 "alien"
   31   27   15    1.53   14.09  0.3029E+02  0.1570E+02 "alien"
   52    2    2    1.49   13.90 -0.1761E+02  0.8885E+01 "alien"
   33   26   12    1.60   13.83  0.4856E+01  0.9686E+01 "alien"
   20   19   18    1.71   13.78  0.8318E+01  0.1072E+02 "alien"
   55    3    0    1.41   13.32 -0.3868E+02  0.1398E+02 "alien"
   27   20   12    1.88   13.20 -0.7389E+00  0.7561E+01 "alien"
   47    1    6    1.60   13.13  0.4611E+01  0.1244E+02 "alien"
   53    3    4    1.45   13.00  0.2663E+02  0.1053E+02 "alien"
   43    2   11    1.61   12.33  0.4331E+01  0.1402E+02 "alien"
   37    2   13    1.72   12.13  0.7319E+01  0.9312E+01 "alien"
   35   24    6    1.76   11.87  0.7160E+01  0.8650E+01 "alien"
   27   21   13    1.81   11.80 -0.6608E+00  0.7511E+01 "alien"
   19   16   19    1.71   11.66  0.7036E+01  0.9223E+01 "alien"
   39    8   19    1.41   11.57 -0.1337E+02  0.1384E+02 "alien"
   40    2   12    1.66   11.41 -0.1504E+01  0.8384E+01 "alien"
   39    2    2    1.98   11.23 -0.1628E+02  0.8935E+01 "alien"
   30    1   13    1.96   11.22 -0.1627E+02  0.1290E+02 "alien"
   40   31    7    1.48   11.00 -0.1394E+02  0.1018E+02 "alien"
   35   32   15    1.38   10.95 -0.1266E+02  0.1659E+02 "alien"
   53    2    1    1.46   10.75 -0.1362E+02  0.1687E+02 "alien"
   54    4    6    1.40   10.74 -0.1242E+02  0.1017E+02 "alien"
   42   10   20    1.32   10.06 -0.3311E+01  0.1562E+02 "alien"


 List of     90 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   46    9   18    1.31   72.90 -0.2399E+02  0.1470E+02
   40   11   21    1.32   32.54 -0.1071E+02  0.1414E+02
   42   10   20    1.32   10.06 -0.3311E+01  0.1562E+02
   35   34   16    1.33   25.67 -0.8447E+01  0.1379E+02
   45    9   18    1.33   19.44 -0.6396E+01  0.1512E+02
   31   26   21    1.33   61.99 -0.2040E+02  0.1350E+02
   41   10   20    1.34   44.54 -0.1550E+02  0.9857E+01
   32   30   19    1.34   60.71 -0.2113E+02  0.1434E+02
   29   20   23    1.34   38.20 -0.1329E+02  0.1349E+02
   35   11   22    1.36   41.72 -0.1452E+02  0.1252E+02
   29   22   22    1.36   21.39 -0.7444E+01  0.9084E+01
   28   19   23    1.36   34.44 -0.1198E+02  0.8877E+01
   39   34   12    1.36   23.19 -0.8071E+01  0.1135E+02
   48    7   15    1.36   36.04 -0.1254E+02  0.1008E+02
   45    8   17    1.36   54.00 -0.1879E+02  0.1596E+02
   35   32   15    1.38   10.95 -0.1266E+02  0.1659E+02
   52    5   10    1.39   17.38 -0.2009E+02  0.1095E+02
   50    5   12    1.39   21.68 -0.2507E+02  0.1044E+02
   54    4    6    1.40   10.74 -0.1242E+02  0.1017E+02
   55    3    0    1.41   13.32 -0.3868E+02  0.1398E+02
   39    8   19    1.41   11.57 -0.1337E+02  0.1384E+02
   53    3    4    1.45   13.00  0.2663E+02  0.1053E+02
   53    2    1    1.46   10.75 -0.1362E+02  0.1687E+02
   40   31    7    1.48   11.00 -0.1394E+02  0.1018E+02
   52    2    2    1.49   13.90 -0.1761E+02  0.8885E+01
   31   27   15    1.53   14.09  0.3029E+02  0.1570E+02
   40   29    0    1.57   14.25  0.1343E+02  0.1274E+02
   35    5   17    1.58  103.87  0.3647E+02  0.1410E+02
   32   26   13    1.59   28.43  0.9983E+01  0.8720E+01
   38    4   15    1.60   15.90  0.5583E+01  0.8017E+01
   47    1    6    1.60   13.13  0.4611E+01  0.1244E+02
   45    2    9    1.60   36.51  0.1282E+02  0.1520E+02
   33   26   12    1.60   13.83  0.4856E+01  0.9686E+01
   43    2   11    1.61   12.33  0.4331E+01  0.1402E+02
   34    5   17    1.61   16.03  0.5627E+01  0.7970E+01
   21   17   20    1.61   35.08  0.1232E+02  0.9843E+01
   39   28    2    1.61   33.94  0.1192E+02  0.8605E+01
   36   27    8    1.63   73.92  0.2595E+02  0.1083E+02
   24   23   17    1.63   21.39  0.7509E+01  0.1236E+02
   42    2   11    1.64   97.45 -0.1285E+02  0.8140E+01
   37   27    6    1.64  142.93 -0.1885E+02  0.1309E+02
   32   25   12    1.65   68.63 -0.9050E+01  0.8332E+01
   47    0    1    1.65   23.34 -0.1697E+02  0.8609E+01
   40    2   12    1.66   11.41 -0.1504E+01  0.8384E+01
   36   26    7    1.67   52.21 -0.6885E+01  0.8207E+01
   37   26    5    1.68  107.25 -0.1414E+02  0.8559E+01
   33    4   16    1.68   38.23 -0.5041E+01  0.7392E+01
   38   26    0    1.69   20.36 -0.1152E+02  0.1502E+02
   35   25    8    1.69   48.33 -0.6373E+01  0.9905E+01
   34   25    9    1.69   57.80 -0.7622E+01  0.7969E+01
   42    1    9    1.70  102.79 -0.1355E+02  0.8562E+01
   32   24   11    1.70   64.18 -0.8464E+01  0.8558E+01
   37   26    3    1.70  107.97 -0.1424E+02  0.1535E+02
   44    0    6    1.70   18.85 -0.9444E+01  0.9096E+01
   22   21   17    1.70   65.74 -0.8669E+01  0.9789E+01
   19   16   19    1.71   11.66  0.7036E+01  0.9223E+01
   20   19   18    1.71   13.78  0.8318E+01  0.1072E+02
   37    2   13    1.72   12.13  0.7319E+01  0.9312E+01
   37   25    2    1.73   33.21  0.2004E+02  0.9829E+01
   34   24    8    1.74   14.30  0.8629E+01  0.8628E+01
   36   25    4    1.74   18.05  0.1089E+02  0.1074E+02
   35   24    6    1.76   11.87  0.7160E+01  0.8650E+01
   34   24    7    1.77   26.14  0.1577E+02  0.1541E+02
   36   24    1    1.79   92.68 -0.5189E+01  0.8194E+01
   33   23    8    1.79  220.40 -0.1234E+02  0.8019E+01
   31    3   15    1.80   65.07 -0.3643E+01  0.1084E+02
   30   22   11    1.80   28.35 -0.1587E+01  0.7894E+01
   43    1    2    1.80   71.79 -0.4019E+01  0.8267E+01
   43    1    1    1.80  343.70 -0.1924E+02  0.8071E+01
   37    1   11    1.80  353.90 -0.1981E+02  0.8325E+01
   43    1    0    1.80   27.24 -0.1187E+02  0.1162E+02
   27   21   13    1.81   11.80 -0.6608E+00  0.7511E+01
   42    1    4    1.81  182.91 -0.1024E+02  0.9076E+01
   34   23    6    1.81  170.71 -0.9557E+01  0.8027E+01
   21   19   16    1.82   44.35 -0.2483E+01  0.1076E+02
   28   21   12    1.83   98.24 -0.5500E+01  0.7607E+01
   35   23    3    1.83   23.03 -0.1289E+01  0.8446E+01
   34   23    5    1.84  255.09 -0.1428E+02  0.1005E+02
   35   23    2    1.84   74.17 -0.4153E+01  0.8325E+01
   35   23    1    1.85  117.13 -0.6558E+01  0.8351E+01
   26   20   13    1.85   59.95 -0.3356E+01  0.7392E+01
   20   18   16    1.86   97.23 -0.5443E+01  0.1010E+02
   40    1    6    1.86   77.40 -0.4333E+01  0.9155E+01
   33   22    6    1.87  189.39 -0.1060E+02  0.8074E+01
   27   20   12    1.88   13.20 -0.7389E+00  0.7561E+01
   34   22    4    1.88  122.10 -0.6836E+01  0.9714E+01
   36    0   10    1.88   19.38 -0.4300E+01  0.8333E+01
   33    1   12    1.90   20.28 -0.2939E+02  0.1206E+02
   30    1   13    1.96   11.22 -0.1627E+02  0.1290E+02
   39    2    2    1.98   11.23 -0.1628E+02  0.8935E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        46
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1126
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1021
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      105


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
