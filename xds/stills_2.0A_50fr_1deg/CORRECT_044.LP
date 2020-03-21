
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.630    77.630    38.856  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      44      44
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      44
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001031  0.000644
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
 BEAM_DIVERGENCE_E.S.D.=   0.049
 REFLECTING_RANGE_E.S.D.=   0.200
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

    1100 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      44 ...      44
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1100 REFLECTIONS ACCEPTED



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
 *  32        oP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.6   77.6   38.9  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     160.0   38.9   77.6  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     160.0   38.9   77.6  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.9  160.0   77.6  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.9  160.0   77.6  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.9  160.0   77.6  90.0  90.0  75.9    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.9  160.0   77.6  90.0  90.0  75.9    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     160.0   38.9  109.8  90.0 133.3  75.9   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.9  109.8  116.5  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.9  160.0  160.0  86.6 104.1 104.1    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.8  116.5   38.9  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    30        mC        997.9      77.6  173.6   38.9  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        997.9      77.6  173.6   38.9  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    41        mC        997.9     173.6   77.6   38.9  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    22        hP        997.9      77.6   77.6   38.9  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.5  116.5  116.5  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.8   86.8  116.5 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.8   86.8  116.5 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
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

  *   0.00        0    4.9      13    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.988   -23.494     0.078
 UNIT_CELL_B-AXIS=    23.494    73.989     0.031
 UNIT_CELL_C-AXIS=    -0.042    -0.003    38.856
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       103 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.27
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.594    77.594    38.827  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.4E-01 2.4E-01 1.2E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012888  0.012888  0.025755  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.957   -23.476     0.019
 COORDINATES OF UNIT CELL B-AXIS    23.476    73.957     0.008
 COORDINATES OF UNIT CELL C-AXIS    -0.010    -0.001    38.827
 CRYSTAL MOSAICITY (DEGREES)     0.200
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999976  0.006495  0.002326
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.82   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1260.44   1242.34
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.23
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.59    77.59    38.83  90.000  90.000  90.000


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
  0.05  0.15        12          4.556      21.228
  0.15  0.25       129         93.199       1.038
  0.25  0.35       157        100.283       0.964
  0.35  0.45       157         42.479       2.277
  0.45  0.55       165        125.101       0.773
  0.55  0.65       170        101.543       0.952
  0.65  0.75       185        122.451       0.790
  0.75  0.85       118         94.115       1.028
  0.85  0.95         7         33.337       2.901
 
 -0.05  0.95      1100         96.716       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               146.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.740
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    43.1 XMAX=    43.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59759 NYBIN=    6
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        306


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               146.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.737
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    82.8 XMAX=  2451.1 NXBIN=    2
 YMIN=     6.2 YMAX=  2517.8 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        306


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               143.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.735
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    43.1 XMAX=    43.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1267    1262
 DETECTOR_SURFACE_POSITION=    1669    1689
 DETECTOR_SURFACE_POSITION=     865    1689
 DETECTOR_SURFACE_POSITION=     865     835
 DETECTOR_SURFACE_POSITION=    1669     835
 DETECTOR_SURFACE_POSITION=    2177    1662
 DETECTOR_SURFACE_POSITION=    1644    2227
 DETECTOR_SURFACE_POSITION=     890    2227
 DETECTOR_SURFACE_POSITION=     357    1662
 DETECTOR_SURFACE_POSITION=     357     862
 DETECTOR_SURFACE_POSITION=     890     297
 DETECTOR_SURFACE_POSITION=    1644     297
 DETECTOR_SURFACE_POSITION=    2177     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        306


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
 7.973E-01  1.943E-03   25.40


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

   24.537   2.877    15.16   1.41      2.81      2.60      78     106      14
    2.877   2.041     2.29   0.88     26.30     27.65     188     228       0
    2.041   1.668    -0.04   0.89    -99.90    -99.90     208     262       0
    1.668   1.446     0.00   0.93   5059.51   5290.69     238     293       0
    1.446   1.294     0.10   0.89    481.72    495.17     138     194       2
   --------------------------------------------------------------------------
   24.537   1.294     1.98   0.95     10.23     10.33     850    1083      16
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      44      44)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.687  -0.3265E+02  0.1311E+02        -2.49           1 
   47    0    0     1.651   0.9666E+01  0.1241E+02         0.78           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     9.5% OF MEAN INTENSITY


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

     3.84          56      37      1297        2.9%       2.5%      2.3%       38   23.53      3.6%    99.8*     0    0.000       0
     2.73          71      44      2129        2.1%       5.6%      5.4%       54   11.52      7.9%    99.4*     0    0.000       0
     2.23         134      78      2698        2.9%      20.5%     20.3%      112    3.57     28.9%    90.8*     0    0.000       0
     1.94         121      72      3140        2.3%     134.2%    162.1%       98    0.45    189.7%    38.2      0    0.000       0
     1.73         152      91      3535        2.6%     -99.9%    -99.9%      122    0.00    -99.9%    -4.7      0    0.000       0
     1.58         151      90      3881        2.3%     624.6%    651.3%      122    0.03    883.3%    -8.3      0    0.000       0
     1.47         168     101      4207        2.4%     -99.9%    -99.9%      134    0.00    -99.9%     6.2      0    0.000       0
     1.37         160      95      4496        2.1%     638.9%    647.3%      130    0.15    903.6%    24.0      0    0.000       0
     1.29          70      50      4759        1.1%     630.4%    672.1%       40    0.04    891.5%    12.0      0    0.000       0
    total        1083     658     30142        2.2%      10.2%     10.3%      850    2.57     14.5%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1100
 NUMBER OF REJECTED MISFITS                                16
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1083
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    658



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

   WILSON LINE (using all data) : A=   6.222 B=  24.441 CORRELATION=  0.83
      #      RES      SS        <I>       log(<I>)       BO
      15     7.256   0.005  1.7833E+03       7.486    -133.1
      18     4.763   0.011  1.4272E+03       7.263     -47.3
      17     3.569   0.020  4.7025E+02       6.153       1.8
      24     2.983   0.028  1.0273E+02       4.632      28.3
      33     2.632   0.036  7.7540E+01       4.351      25.9
      40     2.404   0.043  2.8664E+01       3.356      33.1
      33     2.209   0.051  1.0119E+01       2.314      38.1
      36     2.052   0.059  4.2658E+00       1.451      40.2
      38     1.928   0.067 -7.8897E-01    -999.999   -1000.0
      42     1.825   0.075 -7.4245E-01    -999.999   -1000.0
      45     1.737   0.083 -1.2096E+00    -999.999   -1000.0
      42     1.663   0.090  1.0432E+00       0.042      34.2
      45     1.595   0.098  8.5611E-01      -0.155      32.5
      49     1.533   0.106 -2.2726E-01    -999.999   -1000.0
      49     1.479   0.114 -9.6297E-01    -999.999   -1000.0
      47     1.431   0.122 -3.1098E-01    -999.999   -1000.0
      41     1.389   0.130  4.8757E+00       1.584      17.9
      33     1.348   0.137  4.5539E+00       1.516      17.1
      11     1.310   0.146 -4.6990E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.256        2.146        1.702        1.053
       3     4.763        1.382        0.881        0.418
       1     3.569        0.505        0.124        0.022
       3     2.983        2.520        1.510        0.676
       2     2.632        0.020        0.001        0.000
       2     2.404        0.401        0.089        0.014
       5     2.052        2.986        3.343        2.868
       8     1.737       35.385       49.456      198.229
       4     1.595       52.585       60.618      407.826
       3     1.431       64.120      235.414      812.632
       1     1.389      306.668     1860.346     8061.034
      71   overall       10.973       29.899      152.825


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.256        0.573        0.447        0.275
      15     4.763        1.282        1.664        1.723
      16     3.569        1.140        1.480        1.734
      21     2.983        0.618        0.450        0.279
      31     2.632        1.424        2.204        3.258
      38     2.404        1.017        0.943        0.772
      31     2.209        1.826        1.561        1.461
      31     2.052        2.450        1.854        3.137
      31     1.928       54.607      -65.158     1842.169
      32     1.663       29.622       37.323      311.517
      41     1.595       70.116      -56.441     1697.213
      40     1.389        5.284        6.239       32.079
      32     1.348        4.485        7.556       16.000
     587   overall        7.580       -0.605      139.156

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.256  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.73
       3     4.763  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       1     3.569  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.983  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.632  2.01  1.45  1.20  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       2     2.404  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.73
       5     2.052  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.59
       8     1.737  1.51  1.08  0.90  0.79  0.72  0.67  0.63  0.60  0.57  0.55
       4     1.595  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       3     1.431  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       1     1.389  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      71   overall  0.57  0.41  0.34  0.33  0.30  0.28  0.26  0.25  0.26  0.29


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.256  2.87  3.01  2.10  2.21  1.85  1.61  1.63  1.49  1.38  1.29
      15     4.763  2.80  1.47  1.54  1.21  1.19  1.18  1.46  1.33  1.46  1.37
      16     3.569  1.31  0.69  0.48  0.76  0.95  0.83  1.12  1.25  1.16  1.19
      21     2.983  1.50  1.58  1.29  1.30  1.45  1.37  1.32  1.30  1.28  1.21
      31     2.632  1.02  0.89  1.00  0.98  1.15  1.07  0.96  0.88  1.03  1.02
      38     2.404  1.38  1.31  1.12  1.04  1.00  0.99  1.05  1.05  1.02  0.96
      31     2.209  3.39  1.78  1.37  1.08  0.90  0.86  0.77  0.76  0.71  0.71
      31     2.052  3.73  1.96  1.49  1.37  1.15  1.00  1.03  1.00  0.92  0.87
      31     1.928  5.76  3.03  2.12  1.66  1.39  1.22  1.09  1.00  0.92  0.87
      32     1.663  4.27  2.24  1.57  1.23  1.03  0.90  0.81  0.74  0.68  0.69
      41     1.595  4.36  2.42  1.69  1.33  1.12  0.97  0.87  0.80  0.74  0.69
      40     1.389  4.99  2.62  1.83  1.44  1.33  1.16  1.04  0.95  0.93  0.87
      32     1.348  3.94  2.24  1.57  1.23  1.11  1.04  0.93  0.91  0.95  0.94
     587   overall  2.13  1.25  0.94  0.80  0.74  0.67  0.64  0.61  0.60  0.58


 List of     82 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   24   22   19    1.55   98.69 -0.2243E+02  0.1244E+02 "alien"
   27   24   20    1.44   89.79 -0.2792E+02  0.1240E+02 "alien"
   31   27   15    1.53   79.02 -0.1796E+02  0.1424E+02 "alien"
   42    7   17    1.42   74.23 -0.2308E+02  0.1405E+02 "alien"
   41   33    7    1.42   71.99 -0.2239E+02  0.9824E+01 "alien"
   42    4   14    1.53   69.50 -0.1580E+02  0.8813E+01 "alien"
   40   30    4    1.53   66.43 -0.1510E+02  0.8797E+01 "alien"
   37   29   10    1.52   66.24 -0.1505E+02  0.1431E+02 "alien"
   41   31    1    1.51   64.90 -0.1475E+02  0.9918E+01 "alien"
   51    4   10    1.41   62.07 -0.1930E+02  0.1726E+02 "alien"
   45    3   11    1.55   59.38 -0.1349E+02  0.9676E+01 "alien"
   31   30   17    1.41   57.76 -0.1796E+02  0.1005E+02 "alien"
   35   29   12    1.51   57.59 -0.1309E+02  0.9503E+01 "alien"
   49    5   12    1.42   51.10 -0.1589E+02  0.1014E+02 "alien"
   29   27   19    1.41   47.08 -0.1464E+02  0.1386E+02 "alien"
   38   32   11    1.43   42.43 -0.1320E+02  0.9741E+01 "alien"
   48    2    8    1.53   38.31 -0.8706E+01  0.9785E+01 "alien"
   40    1    5    1.88   37.11 -0.2928E+02  0.1277E+02 "alien"
   39   30    7    1.52   36.39 -0.8271E+01  0.8943E+01 "alien"
   42    1    4    1.81   33.33 -0.2475E+02  0.1303E+02 "alien"
   30   26   15    1.56   32.99 -0.7497E+01  0.1219E+02 "alien"
   50    2    6    1.51   32.88 -0.7471E+01  0.9401E+01 "alien"
   41    4   14    1.56   29.55 -0.6715E+01  0.9768E+01 "alien"
   47    3   10    1.52   29.03 -0.6597E+01  0.9372E+01 "alien"
   50    2    5    1.52   28.36 -0.6445E+01  0.9501E+01 "alien"
   26   21   21    1.45   23.08 -0.7176E+01  0.8783E+01 "alien"
   33    1   12    1.90   22.92 -0.1808E+02  0.1219E+02 "alien"
   44    3   12    1.55   21.93 -0.4983E+01  0.1019E+02 "alien"
   48    1    4    1.59   21.50  0.1841E+02  0.9113E+01 "alien"
   42   32    3    1.46   20.91 -0.2013E+02  0.1173E+02 "alien"
   48    3    9    1.51   20.32 -0.4619E+01  0.1690E+02 "alien"
   36    1    8    1.97   19.70 -0.1554E+02  0.9380E+01 "alien"
   34   24    8    1.74   18.99 -0.2297E+02  0.8901E+01 "alien"
   27   21   13    1.81   18.59 -0.1380E+02  0.7726E+01 "alien"
   30   25   14    1.62   18.27  0.1564E+02  0.8743E+01 "alien"
   54    3    1    1.43   18.06 -0.5617E+01  0.1031E+02 "alien"
   32   30   16    1.43   17.66 -0.5492E+01  0.9475E+01 "alien"
   41    1    5    1.84   17.58 -0.1306E+02  0.1253E+02 "alien"
   22   16   21    1.55   16.89 -0.3838E+01  0.1073E+02 "alien"
   40   31    7    1.48   16.61 -0.1599E+02  0.1033E+02 "alien"
   43    5   15    1.47   16.38 -0.1577E+02  0.1544E+02 "alien"
   31    8   20    1.52   16.28 -0.3699E+01  0.8292E+01 "alien"
   51    2    5    1.49   15.94 -0.1535E+02  0.1001E+02 "alien"
   37   27    6    1.64   15.86  0.1654E+02  0.1357E+02 "alien"
   41   32    6    1.45   15.76 -0.4902E+01  0.9966E+01 "alien"
   40    2   12    1.66   15.45  0.1612E+02  0.8955E+01 "alien"
   27   26   17    1.53   15.35 -0.3488E+01  0.1477E+02 "alien"
   35    2   13    1.78   15.32 -0.1853E+02  0.1179E+02 "alien"
   38   30    9    1.50   15.18 -0.1462E+02  0.9066E+01 "alien"
   20   18   16    1.86   15.17 -0.1126E+02  0.1001E+02 "alien"
   34   23    6    1.81   14.91 -0.1107E+02  0.8271E+01 "alien"
   28   25   16    1.57   14.69  0.1257E+02  0.1412E+02 "alien"
   30   21   10    1.86   14.60 -0.1084E+02  0.1287E+02 "alien"
   31   26   14    1.58   14.30  0.1224E+02  0.8973E+01 "alien"
   35    5   17    1.58   13.88  0.1188E+02  0.1319E+02 "alien"
   37    1   11    1.80   13.77 -0.1022E+02  0.8661E+01 "alien"
   25   25   18    1.54   13.69 -0.8023E+01  0.8209E+01 "alien"
   25   20   21    1.47   13.22 -0.1273E+02  0.1113E+02 "alien"
   43    2   11    1.61   13.16  0.1127E+02  0.1464E+02 "alien"
   36   27    8    1.62   13.13  0.1124E+02  0.1506E+02 "alien"
   45    1    8    1.62   13.12  0.1123E+02  0.1060E+02 "alien"
   26   19   12    1.93   12.87 -0.1015E+02  0.8336E+01 "alien"
   55    3    0    1.41   12.65  0.1479E+03  0.1842E+02 "alien"
   34   23    5    1.84   12.64 -0.9387E+01  0.1058E+02 "alien"
   55    3    2    1.41   12.44  0.6063E+02  0.1228E+02 "alien"
   34    1   12    1.86   12.07 -0.8963E+01  0.8026E+01 "alien"
   35   26    9    1.65   12.05  0.1257E+02  0.8690E+01 "alien"
   54    3    0    1.43   11.91 -0.8456E+01  0.1524E+02 "alien"
   22   21   17    1.70   11.82 -0.1430E+02  0.1053E+02 "alien"
   52    2    1    1.49   11.77 -0.1134E+02  0.9175E+01 "alien"
   51    3    7    1.46   11.67 -0.1124E+02  0.9337E+01 "alien"
   48    1    5    1.58   11.66  0.9985E+01  0.8948E+01 "alien"
   50    1    1    1.55   11.52 -0.2619E+01  0.9372E+01 "alien"
   33    6   18    1.58   11.34  0.9706E+01  0.8435E+01 "alien"
   38   31   10    1.47   11.02 -0.1061E+02  0.9453E+01 "alien"
   42   32    0    1.47   10.72 -0.2178E+02  0.1306E+02 "alien"
   30    1   13    1.95   10.68 -0.8422E+01  0.1282E+02 "alien"
   28    9   21    1.51   10.65 -0.2421E+01  0.1096E+02 "alien"
   41   31    4    1.49   10.42 -0.1003E+02  0.9847E+01 "alien"
   46    1    7    1.61   10.41  0.8912E+01  0.8914E+01 "alien"
   54    3    4    1.42   10.33 -0.3211E+01  0.1007E+02 "alien"
   26   26   18    1.51   10.03 -0.4607E+01  0.8606E+01 "alien"


 List of     82 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   55    3    2    1.41   12.44  0.6063E+02  0.1228E+02
   55    3    0    1.41   12.65  0.1479E+03  0.1842E+02
   51    4   10    1.41   62.07 -0.1930E+02  0.1726E+02
   31   30   17    1.41   57.76 -0.1796E+02  0.1005E+02
   29   27   19    1.41   47.08 -0.1464E+02  0.1386E+02
   49    5   12    1.42   51.10 -0.1589E+02  0.1014E+02
   54    3    4    1.42   10.33 -0.3211E+01  0.1007E+02
   42    7   17    1.42   74.23 -0.2308E+02  0.1405E+02
   41   33    7    1.42   71.99 -0.2239E+02  0.9824E+01
   38   32   11    1.43   42.43 -0.1320E+02  0.9741E+01
   32   30   16    1.43   17.66 -0.5492E+01  0.9475E+01
   54    3    1    1.43   18.06 -0.5617E+01  0.1031E+02
   54    3    0    1.43   11.91 -0.8456E+01  0.1524E+02
   27   24   20    1.44   89.79 -0.2792E+02  0.1240E+02
   26   21   21    1.45   23.08 -0.7176E+01  0.8783E+01
   41   32    6    1.45   15.76 -0.4902E+01  0.9966E+01
   42   32    3    1.46   20.91 -0.2013E+02  0.1173E+02
   51    3    7    1.46   11.67 -0.1124E+02  0.9337E+01
   38   31   10    1.47   11.02 -0.1061E+02  0.9453E+01
   42   32    0    1.47   10.72 -0.2178E+02  0.1306E+02
   25   20   21    1.47   13.22 -0.1273E+02  0.1113E+02
   43    5   15    1.47   16.38 -0.1577E+02  0.1544E+02
   40   31    7    1.48   16.61 -0.1599E+02  0.1033E+02
   52    2    1    1.49   11.77 -0.1134E+02  0.9175E+01
   41   31    4    1.49   10.42 -0.1003E+02  0.9847E+01
   51    2    5    1.49   15.94 -0.1535E+02  0.1001E+02
   38   30    9    1.50   15.18 -0.1462E+02  0.9066E+01
   50    2    6    1.51   32.88 -0.7471E+01  0.9401E+01
   41   31    1    1.51   64.90 -0.1475E+02  0.9918E+01
   26   26   18    1.51   10.03 -0.4607E+01  0.8606E+01
   35   29   12    1.51   57.59 -0.1309E+02  0.9503E+01
   48    3    9    1.51   20.32 -0.4619E+01  0.1690E+02
   28    9   21    1.51   10.65 -0.2421E+01  0.1096E+02
   31    8   20    1.52   16.28 -0.3699E+01  0.8292E+01
   47    3   10    1.52   29.03 -0.6597E+01  0.9372E+01
   39   30    7    1.52   36.39 -0.8271E+01  0.8943E+01
   37   29   10    1.52   66.24 -0.1505E+02  0.1431E+02
   50    2    5    1.52   28.36 -0.6445E+01  0.9501E+01
   31   27   15    1.53   79.02 -0.1796E+02  0.1424E+02
   40   30    4    1.53   66.43 -0.1510E+02  0.8797E+01
   48    2    8    1.53   38.31 -0.8706E+01  0.9785E+01
   42    4   14    1.53   69.50 -0.1580E+02  0.8813E+01
   27   26   17    1.53   15.35 -0.3488E+01  0.1477E+02
   25   25   18    1.54   13.69 -0.8023E+01  0.8209E+01
   44    3   12    1.55   21.93 -0.4983E+01  0.1019E+02
   45    3   11    1.55   59.38 -0.1349E+02  0.9676E+01
   50    1    1    1.55   11.52 -0.2619E+01  0.9372E+01
   24   22   19    1.55   98.69 -0.2243E+02  0.1244E+02
   22   16   21    1.55   16.89 -0.3838E+01  0.1073E+02
   41    4   14    1.56   29.55 -0.6715E+01  0.9768E+01
   30   26   15    1.56   32.99 -0.7497E+01  0.1219E+02
   28   25   16    1.57   14.69  0.1257E+02  0.1412E+02
   31   26   14    1.58   14.30  0.1224E+02  0.8973E+01
   33    6   18    1.58   11.34  0.9706E+01  0.8435E+01
   48    1    5    1.58   11.66  0.9985E+01  0.8948E+01
   35    5   17    1.58   13.88  0.1188E+02  0.1319E+02
   48    1    4    1.59   21.50  0.1841E+02  0.9113E+01
   43    2   11    1.61   13.16  0.1127E+02  0.1464E+02
   46    1    7    1.61   10.41  0.8912E+01  0.8914E+01
   30   25   14    1.62   18.27  0.1564E+02  0.8743E+01
   45    1    8    1.62   13.12  0.1123E+02  0.1060E+02
   36   27    8    1.62   13.13  0.1124E+02  0.1506E+02
   37   27    6    1.64   15.86  0.1654E+02  0.1357E+02
   35   26    9    1.65   12.05  0.1257E+02  0.8690E+01
   40    2   12    1.66   15.45  0.1612E+02  0.8955E+01
   22   21   17    1.70   11.82 -0.1430E+02  0.1053E+02
   34   24    8    1.74   18.99 -0.2297E+02  0.8901E+01
   35    2   13    1.78   15.32 -0.1853E+02  0.1179E+02
   37    1   11    1.80   13.77 -0.1022E+02  0.8661E+01
   27   21   13    1.81   18.59 -0.1380E+02  0.7726E+01
   42    1    4    1.81   33.33 -0.2475E+02  0.1303E+02
   34   23    6    1.81   14.91 -0.1107E+02  0.8271E+01
   34   23    5    1.84   12.64 -0.9387E+01  0.1058E+02
   41    1    5    1.84   17.58 -0.1306E+02  0.1253E+02
   20   18   16    1.86   15.17 -0.1126E+02  0.1001E+02
   30   21   10    1.86   14.60 -0.1084E+02  0.1287E+02
   34    1   12    1.86   12.07 -0.8963E+01  0.8026E+01
   40    1    5    1.88   37.11 -0.2928E+02  0.1277E+02
   33    1   12    1.90   22.92 -0.1808E+02  0.1219E+02
   26   19   12    1.93   12.87 -0.1015E+02  0.8336E+01
   30    1   13    1.95   10.68 -0.8422E+01  0.1282E+02
   36    1    8    1.97   19.70 -0.1554E+02  0.9380E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        31
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1100
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1034
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       66


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
