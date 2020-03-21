
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.625    77.625    38.864  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      45      45
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      45
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001036  0.000635
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
 REFLECTING_RANGE_E.S.D.=   0.198
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

    1096 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      45 ...      45
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1096 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               36
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
 *  31        aP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  1  0  0  0  0 -1  0  0
 *  35        mP          0.0      77.6   38.9   77.6  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  25        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  34        mP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  33        mP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
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
    41        mC        997.4     173.6   77.6   38.9  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        997.4      77.6  173.6   38.9  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        997.4      77.6  173.6   38.9  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        997.4      77.6   77.6   38.9  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
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

  *   0.00        0    4.0      13    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.976   -23.524     0.108
 UNIT_CELL_B-AXIS=    23.524    73.976     0.039
 UNIT_CELL_C-AXIS=    -0.058    -0.002    38.864
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       100 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.28
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.616    77.616    38.840  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.4E-01 2.4E-01 1.1E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012884  0.012884  0.025746  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.986   -23.460     0.044
 COORDINATES OF UNIT CELL B-AXIS    23.460    73.986     0.023
 COORDINATES OF UNIT CELL C-AXIS    -0.024    -0.005    38.840
 CRYSTAL MOSAICITY (DEGREES)     0.198
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999975  0.006627  0.002286
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.89   1232.91
 DETECTOR ORIGIN (PIXELS) AT                     1260.56   1242.55
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.27
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.62    77.62    38.84  90.000  90.000  90.000


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
  0.05  0.15        11        160.852       0.579
  0.15  0.25       124         70.307       1.326
  0.25  0.35       160        103.671       0.899
  0.35  0.45       161         31.109       2.996
  0.45  0.55       156        164.465       0.567
  0.55  0.65       172         66.934       1.393
  0.65  0.75       181        118.218       0.788
  0.75  0.85       127         89.284       1.044
  0.85  0.95         4         41.252       2.259
 
 -0.05  0.95      1096         93.206       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               151.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.776
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    44.1 XMAX=    44.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59721 NYBIN=    6
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        316


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               151.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.775
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    83.3 XMAX=  2450.7 NXBIN=    2
 YMIN=     6.2 YMAX=  2518.3 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        316


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               151.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.772
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    44.1 XMAX=    44.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1267    1262
 DETECTOR_SURFACE_POSITION=    1669    1689
 DETECTOR_SURFACE_POSITION=     865    1689
 DETECTOR_SURFACE_POSITION=     865     836
 DETECTOR_SURFACE_POSITION=    1669     836
 DETECTOR_SURFACE_POSITION=    2177    1662
 DETECTOR_SURFACE_POSITION=    1644    2228
 DETECTOR_SURFACE_POSITION=     890    2228
 DETECTOR_SURFACE_POSITION=     357    1662
 DETECTOR_SURFACE_POSITION=     357     862
 DETECTOR_SURFACE_POSITION=     890     297
 DETECTOR_SURFACE_POSITION=    1644     297
 DETECTOR_SURFACE_POSITION=    2177     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        316


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
 8.268E-01  2.184E-03   23.53


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

   24.544   2.878    14.09   1.40      2.91      2.78      78     102      14
    2.878   2.042     2.20   0.93     23.64     24.65     190     226       0
    2.042   1.669    -0.01   0.85   1392.05   1520.92     216     258       4
    1.669   1.446    -0.02   0.80    -99.90    -99.90     240     293       0
    1.446   1.294    -0.03   0.97    940.49    899.64     144     196       2
   --------------------------------------------------------------------------
   24.544   1.294     1.79   0.92     10.59     10.89     868    1075      20
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      45      45)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.687  -0.3237E+00  0.1348E+02        -0.02           1 
   47    0    0     1.651  -0.6800E+01  0.1214E+02        -0.56           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -6.9% OF MEAN INTENSITY


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

     3.84          56      36      1297        2.8%       2.6%      2.5%       40   22.18      3.7%    99.8*     0    0.000       0
     2.73          68      40      2129        1.9%       5.8%      5.7%       56   10.66      8.3%    99.4*     0    0.000       0
     2.24         131      76      2696        2.8%      19.9%     21.6%      110    3.33     28.2%    92.8*     0    0.000       0
     1.94         119      70      3142        2.2%     144.9%    128.8%       98    0.48    204.9%    13.3      0    0.000       0
     1.73         151      87      3535        2.5%    2271.2%   2643.6%      128   -0.05   3211.9%    10.1      0    0.000       0
     1.58         150      88      3881        2.3%     -99.9%    -99.9%      124    0.00    -99.9%    15.4      0    0.000       0
     1.47         169     100      4205        2.4%     -99.9%    -99.9%      138    0.00    -99.9%     5.9      0    0.000       0
     1.37         159      94      4498        2.1%    1765.4%   1960.3%      130   -0.02   2496.6%    16.8      0    0.000       0
     1.29          72      50      4757        1.1%     663.2%    490.6%       44   -0.13    937.9%   -13.3      0    0.000       0
    total        1075     641     30140        2.1%      10.6%     10.9%      868    2.33     15.0%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1096
 NUMBER OF REJECTED MISFITS                                20
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1075
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    641



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

   WILSON LINE (using all data) : A=   6.350 B=  30.474 CORRELATION=  0.85
      #      RES      SS        <I>       log(<I>)       BO
      14     7.340   0.005  1.8424E+03       7.519    -125.9
      18     4.765   0.011  1.3891E+03       7.236     -40.2
      14     3.613   0.019  4.6583E+02       6.144       5.4
      24     2.984   0.028  9.4149E+01       4.545      32.1
      32     2.628   0.036  7.0807E+01       4.260      28.9
      38     2.405   0.043  2.5223E+01       3.228      36.1
      33     2.209   0.051  1.2702E+01       2.542      37.2
      35     2.054   0.059  2.6809E+00       0.986      45.3
      35     1.929   0.067 -3.7696E+00    -999.999   -1000.0
      41     1.825   0.075  5.7833E-01      -0.548      45.9
      42     1.738   0.083  2.3825E+00       0.868      33.1
      42     1.664   0.090 -1.4074E+00    -999.999   -1000.0
      45     1.596   0.098  4.7428E-01      -0.746      36.1
      48     1.534   0.106 -1.0422E+00    -999.999   -1000.0
      49     1.479   0.114  2.2957E-01      -1.472      34.2
      48     1.431   0.122 -2.2060E-01    -999.999   -1000.0
      39     1.389   0.130 -1.7487E+00    -999.999   -1000.0
      33     1.349   0.137 -2.0492E+00    -999.999   -1000.0
      11     1.310   0.146  3.1350E+00       1.143      17.9


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.340        1.963        1.493        0.886
       3     4.765        1.372        0.875        0.415
       1     3.613        0.438        0.100        0.016
       3     2.984        2.651        1.609        0.725
       2     2.628        0.009        0.000        0.000
       2     2.405        0.157        0.030        0.004
       2     2.209        0.976        0.431        0.147
       5     2.054        1.936        1.694        1.260
       8     1.738        7.801        8.148        7.999
       4     1.596      198.483     1620.313    10753.544
       4     1.534       49.437      190.840      580.093
       3     1.479      492.720     4547.707    30554.123
      71   overall       57.231       38.752     4175.034


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.340        0.568        0.429        0.254
      15     4.765        1.290        1.680        1.744
      13     3.613        1.209        1.639        1.930
      21     2.984        0.618        0.447        0.270
      30     2.628        1.461        2.220        3.163
      36     2.405        0.968        0.738        0.488
      31     2.209        1.200        1.021        0.728
      30     2.054        7.184       11.397       31.748
      33     1.825      124.704      139.050     8278.074
      34     1.738        7.345        3.980       20.310
      11     1.310        5.424        1.526       10.266
     570   overall      115.491      146.598    46895.149

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.340  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.38  0.73
       3     4.765  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       1     3.613  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.984  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.628  2.01  1.45  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       2     2.405  0.00  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  1.47
       2     2.209  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.054  2.42  1.74  1.44  1.27  1.15  1.07  1.00  1.27  1.22  1.17
       8     1.738  1.51  1.08  0.90  0.79  0.72  0.67  0.63  0.60  0.57  0.55
       4     1.596  1.01  0.72  0.60  0.53  0.48  0.45  0.42  0.40  0.38  0.37
       4     1.534  3.02  2.17  1.80  1.58  1.44  1.34  1.26  1.19  1.14  1.10
       3     1.479  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      71   overall  0.74  0.57  0.51  0.45  0.41  0.38  0.35  0.38  0.36  0.41


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.340  3.15  3.31  2.31  2.12  1.78  1.55  1.59  1.45  1.35  1.27
      15     4.765  2.80  1.47  1.54  1.21  1.19  1.18  1.46  1.45  1.46  1.37
      13     3.613  0.81  0.42  0.59  0.70  0.98  0.85  1.22  1.26  1.17  1.22
      21     2.984  1.00  1.31  1.65  1.59  1.45  1.37  1.23  1.30  1.20  1.21
      30     2.628  1.05  0.92  1.03  1.11  1.02  1.18  1.13  1.21  1.12  1.11
      36     2.405  1.46  0.77  0.75  0.67  0.78  0.86  0.88  0.91  0.98  0.92
      31     2.209  2.37  1.78  1.74  1.37  1.15  1.07  1.03  1.00  0.92  0.87
      30     2.054  3.85  2.02  1.41  1.21  1.10  0.96  0.93  0.91  0.84  0.84
      33     1.825  4.78  2.51  1.75  1.38  1.16  1.01  0.90  0.83  0.77  0.72
      34     1.738  4.02  2.11  1.48  1.16  0.97  0.85  0.76  0.69  0.64  0.70
      11     1.310  1.91  1.00  0.70  0.55  0.46  0.60  0.54  0.50  0.61  0.58
     570   overall  1.22  0.75  0.63  0.54  0.49  0.47  0.47  0.46  0.44  0.43


 List of     76 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   42    6   16    1.46  113.86  0.2614E+02  0.1382E+02 "alien"
   49    4   11    1.44  103.66 -0.2287E+02  0.1029E+02 "alien"
   39    8   19    1.41  102.66 -0.2265E+02  0.1360E+02 "alien"
   46    5   14    1.44  101.04 -0.2229E+02  0.1558E+02 "alien"
   37   32   12    1.42   92.40 -0.2038E+02  0.1769E+02 "alien"
   52    2    3    1.48   88.13  0.2023E+02  0.9856E+01 "alien"
   41   32    5    1.47   83.17  0.1909E+02  0.1032E+02 "alien"
   43    6   16    1.44   81.90 -0.1807E+02  0.1292E+02 "alien"
   34   31   15    1.41   68.34 -0.1508E+02  0.9544E+01 "alien"
   40   31    7    1.48   62.13  0.1426E+02  0.1022E+02 "alien"
   34   29   13    1.50   60.02  0.1378E+02  0.1559E+02 "alien"
   47    4   12    1.47   57.21  0.1313E+02  0.1067E+02 "alien"
   35   30   13    1.47   55.00  0.1263E+02  0.9586E+01 "alien"
   45    1    8    1.62   53.11  0.2519E+02  0.1098E+02 "alien"
   33    8   20    1.48   51.05  0.1172E+02  0.1323E+02 "alien"
   39   28    3    1.60   47.04  0.2231E+02  0.9406E+01 "alien"
   39   31    8    1.48   45.17  0.1037E+02  0.1100E+02 "alien"
   54    3    4    1.42   40.46 -0.8925E+01  0.1010E+02 "alien"
   25   20   21    1.47   39.56  0.9082E+01  0.1204E+02 "alien"
   39    7   18    1.45   39.46 -0.8705E+01  0.9170E+01 "alien"
   40   33    9    1.41   39.10 -0.8626E+01  0.9565E+01 "alien"
   42    1    4    1.81   36.61  0.2117E+02  0.9697E+01 "alien"
   36   30   12    1.47   35.94  0.8251E+01  0.9633E+01 "alien"
   52    3    7    1.44   33.70 -0.7434E+01  0.1110E+02 "alien"
   29    9   21    1.50   33.23  0.7630E+01  0.8345E+01 "alien"
   26   21   21    1.45   32.80 -0.7237E+01  0.8624E+01 "alien"
   32   30   16    1.43   32.14 -0.7090E+01  0.9314E+01 "alien"
   40   29    1    1.57   29.50  0.1399E+02  0.9198E+01 "alien"
   43    2   11    1.61   29.11  0.1381E+02  0.1445E+02 "alien"
   51    4    9    1.43   26.29 -0.5800E+01  0.1031E+02 "alien"
   31   29   16    1.46   25.61  0.5880E+01  0.1050E+02 "alien"
   27   24   16    1.61   25.58  0.1213E+02  0.1358E+02 "alien"
   35   31   14    1.42   25.21 -0.5563E+01  0.9404E+01 "alien"
   43    1    3    1.79   23.80  0.1376E+02  0.9746E+01 "alien"
   29   28   17    1.47   23.52  0.5400E+01  0.9741E+01 "alien"
   30   22   11    1.80   23.28  0.1346E+02  0.8342E+01 "alien"
   35   23    3    1.83   22.59  0.1306E+02  0.8937E+01 "alien"
   29   29   18    1.42   22.47 -0.2538E+02  0.1459E+02 "alien"
   40    1    6    1.86   20.71  0.1198E+02  0.1002E+02 "alien"
   45    4   13    1.49   20.06  0.4605E+01  0.1032E+02 "alien"
   28    9   21    1.51   19.89 -0.2073E+02  0.1111E+02 "alien"
   43    1    2    1.80   19.82  0.1147E+02  0.8597E+01 "alien"
   48    1    5    1.58   19.82  0.9400E+01  0.9061E+01 "alien"
   31    3   15    1.80   18.31  0.1059E+02  0.1167E+02 "alien"
   45    2   10    1.58   17.80  0.8441E+01  0.1096E+02 "alien"
   46    1    7    1.61   17.71  0.8398E+01  0.8999E+01 "alien"
   47    2    8    1.56   16.38 -0.1707E+02  0.8970E+01 "alien"
   42   32    0    1.47   16.16  0.1145E+02  0.1376E+02 "alien"
   54    3    1    1.43   16.12 -0.3556E+01  0.1027E+02 "alien"
   39   28    0    1.62   15.58  0.2183E+02  0.1318E+02 "alien"
   52    2    0    1.49   15.09  0.9872E+01  0.1332E+02 "alien"
   36   24    2    1.79   15.03  0.8691E+01  0.8812E+01 "alien"
   20   18   16    1.86   14.34  0.8293E+01  0.1028E+02 "alien"
   42    5   15    1.50   14.14  0.3246E+01  0.9176E+01 "alien"
   30   29   17    1.44   14.00 -0.3089E+01  0.9339E+01 "alien"
   21   17   20    1.61   14.00  0.6638E+01  0.1029E+02 "alien"
   42   32    3    1.46   13.98  0.3210E+01  0.1167E+02 "alien"
   46    4   12    1.49   13.51  0.3102E+01  0.1057E+02 "alien"
   41    1    4    1.86   13.45  0.7777E+01  0.8958E+01 "alien"
   31   24   12    1.69   13.28 -0.1869E+02  0.8115E+01 "alien"
   45    8   17    1.36   13.14 -0.2692E+02  0.1697E+02 "alien"
   35   34   16    1.33   12.71 -0.2604E+02  0.1438E+02 "alien"
   48    2    8    1.53   12.58 -0.1311E+02  0.9689E+01 "alien"
   24   23   17    1.63   12.49 -0.1757E+02  0.1320E+02 "alien"
   40   29    2    1.57   12.48  0.5917E+01  0.9661E+01 "alien"
   44    3   12    1.55   12.46 -0.1298E+02  0.1435E+02 "alien"
   28   23   14    1.70   11.86 -0.1670E+02  0.8065E+01 "alien"
   33   25   11    1.66   11.10 -0.1562E+02  0.1337E+02 "alien"
   36    2    4    2.10   11.06  0.2966E+02  0.9986E+01 "alien"
   39   29    6    1.55   11.01 -0.1147E+02  0.1023E+02 "alien"
   51    2    5    1.49   10.71  0.2459E+01  0.1022E+02 "alien"
   34    6   18    1.56   10.50 -0.1094E+02  0.8263E+01 "alien"
   50    2    5    1.52   10.30 -0.1074E+02  0.9427E+01 "alien"
   31   22   10    1.81   10.25  0.5928E+01  0.1373E+02 "alien"
   33   32   18    1.33   10.15 -0.2080E+02  0.1361E+02 "alien"
   37   27    7    1.62   10.01  0.4749E+01  0.8926E+01 "alien"


 List of     76 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   33   32   18    1.33   10.15 -0.2080E+02  0.1361E+02
   35   34   16    1.33   12.71 -0.2604E+02  0.1438E+02
   45    8   17    1.36   13.14 -0.2692E+02  0.1697E+02
   39    8   19    1.41  102.66 -0.2265E+02  0.1360E+02
   34   31   15    1.41   68.34 -0.1508E+02  0.9544E+01
   40   33    9    1.41   39.10 -0.8626E+01  0.9565E+01
   54    3    4    1.42   40.46 -0.8925E+01  0.1010E+02
   29   29   18    1.42   22.47 -0.2538E+02  0.1459E+02
   35   31   14    1.42   25.21 -0.5563E+01  0.9404E+01
   37   32   12    1.42   92.40 -0.2038E+02  0.1769E+02
   32   30   16    1.43   32.14 -0.7090E+01  0.9314E+01
   51    4    9    1.43   26.29 -0.5800E+01  0.1031E+02
   54    3    1    1.43   16.12 -0.3556E+01  0.1027E+02
   46    5   14    1.44  101.04 -0.2229E+02  0.1558E+02
   52    3    7    1.44   33.70 -0.7434E+01  0.1110E+02
   43    6   16    1.44   81.90 -0.1807E+02  0.1292E+02
   49    4   11    1.44  103.66 -0.2287E+02  0.1029E+02
   30   29   17    1.44   14.00 -0.3089E+01  0.9339E+01
   26   21   21    1.45   32.80 -0.7237E+01  0.8624E+01
   39    7   18    1.45   39.46 -0.8705E+01  0.9170E+01
   31   29   16    1.46   25.61  0.5880E+01  0.1050E+02
   42   32    3    1.46   13.98  0.3210E+01  0.1167E+02
   42    6   16    1.46  113.86  0.2614E+02  0.1382E+02
   41   32    5    1.47   83.17  0.1909E+02  0.1032E+02
   47    4   12    1.47   57.21  0.1313E+02  0.1067E+02
   35   30   13    1.47   55.00  0.1263E+02  0.9586E+01
   42   32    0    1.47   16.16  0.1145E+02  0.1376E+02
   25   20   21    1.47   39.56  0.9082E+01  0.1204E+02
   29   28   17    1.47   23.52  0.5400E+01  0.9741E+01
   36   30   12    1.47   35.94  0.8251E+01  0.9633E+01
   40   31    7    1.48   62.13  0.1426E+02  0.1022E+02
   33    8   20    1.48   51.05  0.1172E+02  0.1323E+02
   52    2    3    1.48   88.13  0.2023E+02  0.9856E+01
   39   31    8    1.48   45.17  0.1037E+02  0.1100E+02
   45    4   13    1.49   20.06  0.4605E+01  0.1032E+02
   52    2    0    1.49   15.09  0.9872E+01  0.1332E+02
   46    4   12    1.49   13.51  0.3102E+01  0.1057E+02
   51    2    5    1.49   10.71  0.2459E+01  0.1022E+02
   42    5   15    1.50   14.14  0.3246E+01  0.9176E+01
   29    9   21    1.50   33.23  0.7630E+01  0.8345E+01
   34   29   13    1.50   60.02  0.1378E+02  0.1559E+02
   28    9   21    1.51   19.89 -0.2073E+02  0.1111E+02
   50    2    5    1.52   10.30 -0.1074E+02  0.9427E+01
   48    2    8    1.53   12.58 -0.1311E+02  0.9689E+01
   44    3   12    1.55   12.46 -0.1298E+02  0.1435E+02
   39   29    6    1.55   11.01 -0.1147E+02  0.1023E+02
   34    6   18    1.56   10.50 -0.1094E+02  0.8263E+01
   47    2    8    1.56   16.38 -0.1707E+02  0.8970E+01
   40   29    2    1.57   12.48  0.5917E+01  0.9661E+01
   40   29    1    1.57   29.50  0.1399E+02  0.9198E+01
   45    2   10    1.58   17.80  0.8441E+01  0.1096E+02
   48    1    5    1.58   19.82  0.9400E+01  0.9061E+01
   39   28    3    1.60   47.04  0.2231E+02  0.9406E+01
   43    2   11    1.61   29.11  0.1381E+02  0.1445E+02
   21   17   20    1.61   14.00  0.6638E+01  0.1029E+02
   27   24   16    1.61   25.58  0.1213E+02  0.1358E+02
   46    1    7    1.61   17.71  0.8398E+01  0.8999E+01
   39   28    0    1.62   15.58  0.2183E+02  0.1318E+02
   37   27    7    1.62   10.01  0.4749E+01  0.8926E+01
   45    1    8    1.62   53.11  0.2519E+02  0.1098E+02
   24   23   17    1.63   12.49 -0.1757E+02  0.1320E+02
   33   25   11    1.66   11.10 -0.1562E+02  0.1337E+02
   31   24   12    1.69   13.28 -0.1869E+02  0.8115E+01
   28   23   14    1.70   11.86 -0.1670E+02  0.8065E+01
   36   24    2    1.79   15.03  0.8691E+01  0.8812E+01
   43    1    3    1.79   23.80  0.1376E+02  0.9746E+01
   31    3   15    1.80   18.31  0.1059E+02  0.1167E+02
   30   22   11    1.80   23.28  0.1346E+02  0.8342E+01
   43    1    2    1.80   19.82  0.1147E+02  0.8597E+01
   31   22   10    1.81   10.25  0.5928E+01  0.1373E+02
   42    1    4    1.81   36.61  0.2117E+02  0.9697E+01
   35   23    3    1.83   22.59  0.1306E+02  0.8937E+01
   20   18   16    1.86   14.34  0.8293E+01  0.1028E+02
   41    1    4    1.86   13.45  0.7777E+01  0.8958E+01
   40    1    6    1.86   20.71  0.1198E+02  0.1002E+02
   36    2    4    2.10   11.06  0.2966E+02  0.9986E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        40
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1096
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1007
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       89


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
