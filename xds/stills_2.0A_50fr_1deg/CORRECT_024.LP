
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.491    77.491    38.774  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      24      24
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      24
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000901  0.000697
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
 BEAM_DIVERGENCE_E.S.D.=   0.045
 REFLECTING_RANGE_E.S.D.=   0.121
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

    1377 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      24 ...      24
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1377 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               43
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
 *  31        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0 -1  0  0  0  0 -1  0  0
 *  34        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  33        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.5   38.8   77.5  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  25        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  32        oP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  20        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.8   38.8   77.5  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.8   38.8   77.5  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    28        mC        250.0      38.8  159.8   77.5  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    29        mC        250.0      38.8  159.8   77.5  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    36        oC        250.0      38.8  159.8   77.5  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.8   77.5  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    27        mC        500.0     159.8   38.8  109.6  90.0 133.3  76.0   -2  0 -1  0  0  0 -1  0  1 -1  0  0
    19        oI        500.0      38.8  109.6  116.2  90.0  70.5  90.0    0  0 -1  0  1 -1  0  0 -1 -1 -1  0
    26        oF        625.0      38.8  159.8  159.8  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.6  116.2   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        998.5     173.3   77.5   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.5      77.5  173.3   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.5      77.5  173.3   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.5      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.2  116.2  116.2  96.4  96.4 141.0    1 -1  1  0 -1  1  1  0 -1 -1 -1  0
     2        hR        999.0      86.7   86.7  116.2 116.5  63.5 101.6    1  0  1  0  0 -1 -1  0  1  1 -1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.7   86.7  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.7   86.7  116.2 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.6   86.7   86.7  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.7   86.7  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.7   86.7  109.6  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.7  248.1 102.1  99.0 116.6    0  0  1  0 -1  0 -1  0  1 -3 -1  0
    10        mC        999.0      86.7   86.7   77.5  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.5  177.6  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.7   86.7  177.6 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.7   86.7   86.7  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.6  109.6   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.5  177.6 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.7  177.6   77.5 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    5.5      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.850   -23.477     0.044
 UNIT_CELL_B-AXIS=    23.477    73.850     0.017
 UNIT_CELL_C-AXIS=    -0.023    -0.001    38.775
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       221 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.29
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.396    77.396    38.735  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  7.4E-02 7.4E-02 3.0E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012921  0.012921  0.025816  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.759   -23.445     0.068
 COORDINATES OF UNIT CELL B-AXIS    23.445    73.759     0.020
 COORDINATES OF UNIT CELL C-AXIS    -0.035     0.001    38.735
 CRYSTAL MOSAICITY (DEGREES)     0.121
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006526  0.002076
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.90   1232.90
 DETECTOR ORIGIN (PIXELS) AT                     1260.88   1242.39
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.09
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.40    77.40    38.74  90.000  90.000  90.000


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
  0.05  0.15        77        114.959       1.283
  0.15  0.25       173        166.952       0.884
  0.25  0.35       190        129.868       1.136
  0.35  0.45       145         89.142       1.655
  0.45  0.55       180        183.047       0.806
  0.55  0.65       187        138.380       1.066
  0.65  0.75       201        178.353       0.827
  0.75  0.85       148        140.984       1.046
  0.85  0.95        76        161.610       0.913
 
 -0.05  0.95      1377        147.536       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               246.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.789
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    23.1 XMAX=    23.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59535 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               246.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.786
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     9.5 XMAX=  2454.8 NXBIN=    3
 YMIN=     8.3 YMAX=  2520.1 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               246.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.782
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    23.1 XMAX=    23.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1647    1691
 DETECTOR_SURFACE_POSITION=     817    1691
 DETECTOR_SURFACE_POSITION=     817     838
 DETECTOR_SURFACE_POSITION=    1647     838
 DETECTOR_SURFACE_POSITION=    2172    1664
 DETECTOR_SURFACE_POSITION=    1621    2230
 DETECTOR_SURFACE_POSITION=     843    2230
 DETECTOR_SURFACE_POSITION=     292    1664
 DETECTOR_SURFACE_POSITION=     292     864
 DETECTOR_SURFACE_POSITION=     843     299
 DETECTOR_SURFACE_POSITION=    1621     299
 DETECTOR_SURFACE_POSITION=    2172     864
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


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
 9.766E-01  4.673E-03   14.80


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

   24.475   3.153    12.73   1.50      3.55      3.95      86     108      14
    3.153   2.239     7.24   1.06      6.57      6.11     196     225       2
    2.239   1.830     1.32   0.55     25.48     34.17     206     239       0
    1.830   1.586     0.01   0.61  28044.26  34667.14     246     283       0
    1.586   1.419     0.03   0.63    738.15    924.40     276     324       4
    1.419   1.296    -0.02   0.78    -99.90    -99.90     132     175       2
   --------------------------------------------------------------------------
   24.475   1.296     2.46   0.77      7.59      8.57    1142    1354      22
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      24      24)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.683   0.3574E+02  0.1316E+02         2.72           1 
   47    0    0     1.647  -0.2252E+02  0.1239E+02        -1.82           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS   -14.3% OF MEAN INTENSITY


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

     3.85          67      40      1285        3.1%       3.0%      3.9%       54   16.68      4.2%    99.8*     0    0.000       0
     2.74         103      60      2098        2.9%       6.0%      4.5%       86   13.49      8.5%    98.0*     0    0.000       0
     2.24         163      92      2665        3.5%       7.1%      7.0%      142    8.73     10.1%    98.7*     0    0.000       0
     1.94         155      87      3101        2.8%      18.8%     23.4%      136    2.62     26.7%    95.3*     0    0.000       0
     1.74         165      94      3481        2.7%     -99.9%    -99.9%      142    0.00    -99.9%     3.3      0    0.000       0
     1.59         202     115      3848        3.0%    2033.0%   2629.1%      174    0.05   2875.1%     9.3      0    0.000       0
     1.47         208     121      4143        2.9%    4668.4%   5909.1%      174   -0.02   6602.1%    16.8      0    0.000       0
     1.37         207     116      4429        2.6%     500.2%    610.5%      182    0.04    707.4%     6.8      0    0.000       0
     1.30          84      58      4743        1.2%     -99.9%    -99.9%       52    0.00    -99.9%   -25.6      0    0.000       0
    total        1354     783     29793        2.6%       7.6%      8.6%     1142    3.21     10.7%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1377
 NUMBER OF REJECTED MISFITS                                22
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1354
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    783



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

   WILSON LINE (using all data) : A=   7.446 B=  36.366 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      17     7.326   0.005  1.9245E+03       7.562     -12.5
      19     4.699   0.011  2.0249E+03       7.613      -7.4
      26     3.518   0.020  9.0294E+02       6.806      15.9
      30     2.997   0.028  2.9330E+02       5.681      31.7
      37     2.642   0.036  2.6287E+02       5.572      26.2
      49     2.400   0.043  1.3728E+02       4.922      29.1
      38     2.209   0.051  6.8160E+01       4.222      31.5
      44     2.056   0.059  2.5280E+01       3.230      35.7
      42     1.931   0.067  9.2984E-01      -0.073      56.1
      50     1.830   0.075  1.2601E+00       0.231      48.3
      45     1.740   0.083 -2.2557E+00    -999.999   -1000.0
      51     1.667   0.090  1.1527E+00       0.142      40.6
      52     1.598   0.098  2.1833E-01      -1.522      45.8
      61     1.537   0.106  3.3174E-02      -3.406      51.3
      68     1.480   0.114  1.0588E+00       0.057      32.4
      52     1.434   0.122  7.7195E-01      -0.259      31.7
      45     1.392   0.129 -5.1593E-01    -999.999   -1000.0
      39     1.352   0.137 -2.0553E+00    -999.999   -1000.0
      18     1.314   0.145  1.7722E+00       0.572      23.7


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.326        1.875        1.300        0.722
       3     4.699        1.640        1.147        0.596
       2     3.518        0.397        0.087        0.014
       4     2.997        1.256        0.590        0.209
       2     2.642        0.262        0.066        0.012
       5     2.400        0.842        0.435        0.177
       2     2.209        1.267        0.513        0.151
       7     2.056        0.230        0.082        0.022
       6     1.931       27.545       48.982      114.890
      10     1.830       22.677        7.609       97.683
      11     1.667       15.990       20.552       64.396
       3     1.537    55706.812  7479757.304748954350.450
      81   overall     2092.076   276921.876 27739713.096


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.326        0.471        0.343        0.200
      16     4.699        1.240        1.643        1.743
      24     3.518        0.844        0.685        0.521
      26     2.997        0.862        0.753        0.566
      35     2.642        1.180        1.376        1.488
      44     2.400        0.947        0.823        0.607
      36     2.209        0.778        0.576        0.381
      37     2.056        1.661        2.544        3.635
      36     1.931       48.697      191.475     1833.492
      40     1.830       17.050        9.852      134.715
      40     1.667       18.757       19.340      139.664
      49     1.598      754.433     1502.630   327665.995
      62     1.480       33.830       -4.225      488.801
      49     1.434       68.756      126.054     2303.966
      44     1.392       29.601      -12.664      407.790
      18     1.314       39.690      213.513     1384.449
     702   overall     2430.350    50125.305 28578443.620

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.326  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.32  0.30  0.29
       3     4.699  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     3.518  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       4     2.997  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       2     2.642  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.400  0.00  0.00  0.00  0.42  0.38  0.71  0.67  0.95  0.91  0.88
       2     2.209  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       7     2.056  1.15  1.65  1.71  1.51  1.37  1.27  1.20  1.14  1.30  1.26
       6     1.931  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      10     1.830  1.21  0.87  0.72  0.63  0.58  0.53  0.50  0.48  0.46  0.44
      11     1.667  1.46  1.05  0.87  0.77  0.70  0.65  0.61  0.58  0.55  0.53
       3     1.537  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      81   overall  0.70  0.57  0.50  0.47  0.43  0.44  0.41  0.43  0.43  0.42


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.326  2.63  2.76  2.25  1.77  1.91  1.66  1.66  1.51  1.40  1.32
      16     4.699  1.97  1.38  1.21  1.14  1.11  1.66  1.49  1.36  1.37  1.38
      24     3.518  0.88  0.69  0.64  0.51  0.64  0.74  1.08  0.98  1.05  0.99
      26     2.997  1.62  1.91  1.34  1.28  1.37  1.28  1.15  1.26  1.23  1.16
      35     2.642  0.60  1.26  1.10  1.13  0.94  0.95  1.02  1.09  1.16  1.08
      44     2.400  0.96  1.00  0.96  1.03  1.10  1.11  1.08  1.11  1.11  1.04
      36     2.209  1.17  0.77  0.64  0.59  0.64  0.86  1.05  1.01  0.98  1.01
      37     2.056  1.99  1.34  1.04  0.82  0.82  1.02  1.02  1.03  0.96  0.98
      36     1.931  6.71  3.52  2.47  1.94  1.62  1.42  1.27  1.16  1.08  1.01
      40     1.830  4.99  2.76  1.93  1.52  1.33  1.16  1.04  0.95  0.88  0.83
      40     1.667  4.47  2.34  1.64  1.29  1.08  1.00  0.94  0.86  0.80  0.75
      49     1.598  5.15  2.70  1.89  1.49  1.24  1.09  0.97  0.89  0.83  0.77
      62     1.480  4.41  2.40  1.68  1.32  1.11  0.97  0.90  0.82  0.76  0.74
      49     1.434  4.29  2.36  1.65  1.30  1.09  0.95  0.85  0.78  0.76  0.71
      44     1.392  5.49  2.88  2.02  1.59  1.33  1.16  1.04  0.95  0.92  0.86
      18     1.314  6.42  3.37  2.57  2.02  1.69  1.48  1.32  1.21  1.12  1.05
     702   overall  2.87  1.71  1.26  1.04  0.93  0.89  0.85  0.81  0.78  0.75


 List of     93 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   32   27   14    1.54  642.53  0.2132E+02  0.1092E+02 "alien"
   37   29   10    1.52  419.34  0.1391E+02  0.1671E+02 "alien"
   31    7   19    1.56  393.51  0.1305E+02  0.1161E+02 "alien"
   39   29    5    1.56  389.01  0.1291E+02  0.9336E+01 "alien"
   40   30    4    1.53  361.42  0.1199E+02  0.9752E+01 "alien"
   51    2    2    1.51  352.44  0.1169E+02  0.1189E+02 "alien"
   38   29    8    1.54  340.62  0.1130E+02  0.9193E+01 "alien"
   25   23   19    1.52  302.53  0.1004E+02  0.1217E+02 "alien"
   27    9   21    1.53  295.64  0.9808E+01  0.9359E+01 "alien"
   39   29    6    1.55  227.65  0.7552E+01  0.1004E+02 "alien"
   41    4   14    1.55  220.73  0.7323E+01  0.8772E+01 "alien"
   36   28   10    1.55  198.70  0.6592E+01  0.9128E+01 "alien"
   46    3   11    1.52  180.77  0.5997E+01  0.9182E+01 "alien"
   22   16   21    1.55  148.41  0.4923E+01  0.6550E+01 "alien"
   48    2    7    1.55  146.55  0.4862E+01  0.9188E+01 "alien"
   48    2    8    1.53  134.40  0.4459E+01  0.9945E+01 "alien"
   28    9   21    1.51  125.60  0.4167E+01  0.6988E+01 "alien"
   22   10   22    1.54  118.68  0.3937E+01  0.8304E+01 "alien"
   24   22   19    1.55  112.83  0.3743E+01  0.8255E+01 "alien"
   38   28    6    1.59   89.16  0.1947E+02  0.1335E+02 "alien"
   43    2   11    1.60   85.16  0.1859E+02  0.1036E+02 "alien"
   32   26   13    1.59   71.83  0.1568E+02  0.8981E+01 "alien"
   35   28   11    1.55   61.81  0.2051E+01  0.8996E+01 "alien"
   39   28    3    1.60   60.80  0.1327E+02  0.9447E+01 "alien"
   50    1    1    1.55   57.54  0.1909E+01  0.9289E+01 "alien"
   41    3   13    1.59   56.75  0.1239E+02  0.8719E+01 "alien"
   33    5   17    1.63   56.46  0.1233E+02  0.7349E+01 "alien"
   38    5   16    1.55   54.30  0.1801E+01  0.8266E+01 "alien"
   33   26   12    1.60   49.96  0.1091E+02  0.9881E+01 "alien"
   45    2    9    1.60   48.58  0.1061E+02  0.1033E+02 "alien"
   51    4   10    1.41   39.66 -0.2046E+02  0.1214E+02 "alien"
   24   24   18    1.57   36.67  0.2312E+02  0.8451E+01 "alien"
   46    2    9    1.57   36.45  0.1209E+01  0.8982E+01 "alien"
   48    1    4    1.59   35.46  0.7742E+01  0.9350E+01 "alien"
   37   27    7    1.62   34.03  0.7429E+01  0.9031E+01 "alien"
   33   22    5    1.89   33.55  0.3119E+02  0.1725E+02 "alien"
   51    3    8    1.45   33.54  0.2589E+02  0.1166E+02 "alien"
   34    5   17    1.60   30.82  0.6729E+01  0.7625E+01 "alien"
   47    1    6    1.60   26.18  0.5716E+01  0.1315E+02 "alien"
   53    7   13    1.30   25.81  0.4575E+02  0.1490E+02 "alien"
   38    9   20    1.39   22.09 -0.1140E+02  0.9205E+01 "alien"
   42    8   18    1.39   22.07 -0.1139E+02  0.9483E+01 "alien"
   52    3    6    1.45   21.91  0.1691E+02  0.9679E+01 "alien"
   48    6   14    1.39   20.33 -0.1049E+02  0.1368E+02 "alien"
   49    5   12    1.41   19.67  0.1518E+02  0.9812E+01 "alien"
   45    7   16    1.39   19.55 -0.1009E+02  0.1068E+02 "alien"
   22   18   20    1.58   19.32  0.4217E+01  0.6811E+01 "alien"
   50    8   15    1.32   19.08  0.3381E+02  0.1877E+02 "alien"
   27   18   23    1.38   18.82 -0.9711E+01  0.8599E+01 "alien"
   37   33   13    1.38   18.82 -0.9708E+01  0.9834E+01 "alien"
   37   29    9    1.54   18.68  0.6197E+00  0.9181E+01 "alien"
   38    4   15    1.59   18.54  0.4049E+01  0.7919E+01 "alien"
   52    5   10    1.38   18.28 -0.9432E+01  0.1032E+02 "alien"
   20   13   21    1.60   17.17  0.3748E+01  0.9188E+01 "alien"
   55    3    0    1.41   16.89 -0.2841E+02  0.1226E+02 "alien"
   31   25   13    1.63   16.63  0.3631E+01  0.8141E+01 "alien"
   32   30   16    1.43   16.52  0.1275E+02  0.9442E+01 "alien"
   37   31   11    1.46   16.32  0.1728E+02  0.1192E+02 "alien"
   36    1    8    1.96   15.75  0.1465E+02  0.9006E+01 "alien"
   17    9   21    1.68   15.23  0.1756E+02  0.7010E+01 "alien"
   39    4   15    1.57   15.05  0.3287E+01  0.1163E+02 "alien"
   47    1    5    1.61   15.05  0.3285E+01  0.9038E+01 "alien"
   37   32   12    1.42   15.02  0.1160E+02  0.1081E+02 "alien"
   38    8   19    1.43   14.81  0.1143E+02  0.8395E+01 "alien"
   39    2    1    1.98   14.73  0.1369E+02  0.9354E+01 "alien"
   26   13   24    1.38   14.50 -0.7481E+01  0.9819E+01 "alien"
   36   31   13    1.43   14.34  0.1107E+02  0.1012E+02 "alien"
   41   31    1    1.50   14.33  0.1517E+02  0.9893E+01 "alien"
   30   28   16    1.49   14.19  0.1502E+02  0.1343E+02 "alien"
   42    4   14    1.53   14.05  0.4661E+00  0.8559E+01 "alien"
   50    1    0    1.55   14.01  0.1226E+01  0.1296E+02 "alien"
   41   31    4    1.49   13.97  0.1479E+02  0.1004E+02 "alien"
   40   32    7    1.46   12.66  0.9772E+01  0.1096E+02 "alien"
   44    5   14    1.48   12.48  0.1322E+02  0.1626E+02 "alien"
   52    2    2    1.48   11.70  0.1239E+02  0.9686E+01 "alien"
   31   22   10    1.80   11.66  0.1470E+02  0.1401E+02 "alien"
   30   21    9    1.90   11.66  0.1084E+02  0.1589E+02 "alien"
   42    6   16    1.46   11.62  0.8966E+01  0.9323E+01 "alien"
   38    7   18    1.47   11.59  0.1227E+02  0.8910E+01 "alien"
   32   21    6    1.93   11.54  0.1073E+02  0.8937E+01 "alien"
   32   30   19    1.33   11.40 -0.2342E+02  0.1515E+02 "alien"
   54    3    2    1.43   11.25  0.8687E+01  0.9981E+01 "alien"
   24   19   14    1.87   11.21  0.1413E+02  0.1267E+02 "alien"
   40    1    5    1.88   10.98  0.1384E+02  0.1254E+02 "alien"
   49    6   14    1.36   10.88 -0.2235E+02  0.1701E+02 "alien"
   37   27    6    1.63   10.64  0.1226E+02  0.1276E+02 "alien"
   38    1    7    1.91   10.58  0.9835E+01  0.8811E+01 "alien"
   39   29    4    1.57   10.48  0.2288E+01  0.1022E+02 "alien"
   45    1    8    1.62   10.40  0.2270E+01  0.9549E+01 "alien"
   39   32    9    1.45   10.30  0.7951E+01  0.1076E+02 "alien"
   33   25   11    1.65   10.22  0.1179E+02  0.1410E+02 "alien"
   34   31   15    1.41   10.21 -0.5266E+01  0.9262E+01 "alien"
   37    1    7    1.96   10.00  0.9300E+01  0.8915E+01 "alien"


 List of     93 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   25.81  0.4575E+02  0.1490E+02
   50    8   15    1.32   19.08  0.3381E+02  0.1877E+02
   32   30   19    1.33   11.40 -0.2342E+02  0.1515E+02
   49    6   14    1.36   10.88 -0.2235E+02  0.1701E+02
   27   18   23    1.38   18.82 -0.9711E+01  0.8599E+01
   26   13   24    1.38   14.50 -0.7481E+01  0.9819E+01
   37   33   13    1.38   18.82 -0.9708E+01  0.9834E+01
   52    5   10    1.38   18.28 -0.9432E+01  0.1032E+02
   48    6   14    1.39   20.33 -0.1049E+02  0.1368E+02
   38    9   20    1.39   22.09 -0.1140E+02  0.9205E+01
   42    8   18    1.39   22.07 -0.1139E+02  0.9483E+01
   45    7   16    1.39   19.55 -0.1009E+02  0.1068E+02
   55    3    0    1.41   16.89 -0.2841E+02  0.1226E+02
   51    4   10    1.41   39.66 -0.2046E+02  0.1214E+02
   34   31   15    1.41   10.21 -0.5266E+01  0.9262E+01
   49    5   12    1.41   19.67  0.1518E+02  0.9812E+01
   37   32   12    1.42   15.02  0.1160E+02  0.1081E+02
   38    8   19    1.43   14.81  0.1143E+02  0.8395E+01
   32   30   16    1.43   16.52  0.1275E+02  0.9442E+01
   54    3    2    1.43   11.25  0.8687E+01  0.9981E+01
   36   31   13    1.43   14.34  0.1107E+02  0.1012E+02
   39   32    9    1.45   10.30  0.7951E+01  0.1076E+02
   51    3    8    1.45   33.54  0.2589E+02  0.1166E+02
   52    3    6    1.45   21.91  0.1691E+02  0.9679E+01
   42    6   16    1.46   11.62  0.8966E+01  0.9323E+01
   40   32    7    1.46   12.66  0.9772E+01  0.1096E+02
   37   31   11    1.46   16.32  0.1728E+02  0.1192E+02
   38    7   18    1.47   11.59  0.1227E+02  0.8910E+01
   44    5   14    1.48   12.48  0.1322E+02  0.1626E+02
   52    2    2    1.48   11.70  0.1239E+02  0.9686E+01
   30   28   16    1.49   14.19  0.1502E+02  0.1343E+02
   41   31    4    1.49   13.97  0.1479E+02  0.1004E+02
   41   31    1    1.50   14.33  0.1517E+02  0.9893E+01
   28    9   21    1.51  125.60  0.4167E+01  0.6988E+01
   51    2    2    1.51  352.44  0.1169E+02  0.1189E+02
   37   29   10    1.52  419.34  0.1391E+02  0.1671E+02
   46    3   11    1.52  180.77  0.5997E+01  0.9182E+01
   25   23   19    1.52  302.53  0.1004E+02  0.1217E+02
   27    9   21    1.53  295.64  0.9808E+01  0.9359E+01
   40   30    4    1.53  361.42  0.1199E+02  0.9752E+01
   48    2    8    1.53  134.40  0.4459E+01  0.9945E+01
   42    4   14    1.53   14.05  0.4661E+00  0.8559E+01
   38   29    8    1.54  340.62  0.1130E+02  0.9193E+01
   32   27   14    1.54  642.53  0.2132E+02  0.1092E+02
   37   29    9    1.54   18.68  0.6197E+00  0.9181E+01
   22   10   22    1.54  118.68  0.3937E+01  0.8304E+01
   39   29    6    1.55  227.65  0.7552E+01  0.1004E+02
   50    1    1    1.55   57.54  0.1909E+01  0.9289E+01
   48    2    7    1.55  146.55  0.4862E+01  0.9188E+01
   24   22   19    1.55  112.83  0.3743E+01  0.8255E+01
   50    1    0    1.55   14.01  0.1226E+01  0.1296E+02
   22   16   21    1.55  148.41  0.4923E+01  0.6550E+01
   35   28   11    1.55   61.81  0.2051E+01  0.8996E+01
   38    5   16    1.55   54.30  0.1801E+01  0.8266E+01
   41    4   14    1.55  220.73  0.7323E+01  0.8772E+01
   36   28   10    1.55  198.70  0.6592E+01  0.9128E+01
   39   29    5    1.56  389.01  0.1291E+02  0.9336E+01
   31    7   19    1.56  393.51  0.1305E+02  0.1161E+02
   24   24   18    1.57   36.67  0.2312E+02  0.8451E+01
   46    2    9    1.57   36.45  0.1209E+01  0.8982E+01
   39    4   15    1.57   15.05  0.3287E+01  0.1163E+02
   39   29    4    1.57   10.48  0.2288E+01  0.1022E+02
   22   18   20    1.58   19.32  0.4217E+01  0.6811E+01
   32   26   13    1.59   71.83  0.1568E+02  0.8981E+01
   38   28    6    1.59   89.16  0.1947E+02  0.1335E+02
   48    1    4    1.59   35.46  0.7742E+01  0.9350E+01
   41    3   13    1.59   56.75  0.1239E+02  0.8719E+01
   38    4   15    1.59   18.54  0.4049E+01  0.7919E+01
   47    1    6    1.60   26.18  0.5716E+01  0.1315E+02
   45    2    9    1.60   48.58  0.1061E+02  0.1033E+02
   39   28    3    1.60   60.80  0.1327E+02  0.9447E+01
   33   26   12    1.60   49.96  0.1091E+02  0.9881E+01
   43    2   11    1.60   85.16  0.1859E+02  0.1036E+02
   34    5   17    1.60   30.82  0.6729E+01  0.7625E+01
   20   13   21    1.60   17.17  0.3748E+01  0.9188E+01
   47    1    5    1.61   15.05  0.3285E+01  0.9038E+01
   37   27    7    1.62   34.03  0.7429E+01  0.9031E+01
   45    1    8    1.62   10.40  0.2270E+01  0.9549E+01
   33    5   17    1.63   56.46  0.1233E+02  0.7349E+01
   31   25   13    1.63   16.63  0.3631E+01  0.8141E+01
   37   27    6    1.63   10.64  0.1226E+02  0.1276E+02
   33   25   11    1.65   10.22  0.1179E+02  0.1410E+02
   17    9   21    1.68   15.23  0.1756E+02  0.7010E+01
   31   22   10    1.80   11.66  0.1470E+02  0.1401E+02
   24   19   14    1.87   11.21  0.1413E+02  0.1267E+02
   40    1    5    1.88   10.98  0.1384E+02  0.1254E+02
   33   22    5    1.89   33.55  0.3119E+02  0.1725E+02
   30   21    9    1.90   11.66  0.1084E+02  0.1589E+02
   38    1    7    1.91   10.58  0.9835E+01  0.8811E+01
   32   21    6    1.93   11.54  0.1073E+02  0.8937E+01
   37    1    7    1.96   10.00  0.9300E+01  0.8915E+01
   36    1    8    1.96   15.75  0.1465E+02  0.9006E+01
   39    2    1    1.98   14.73  0.1369E+02  0.9354E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        44
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1377
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1277
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      100


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
