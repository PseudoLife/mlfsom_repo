
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.494    77.494    38.781  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      25      25
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      25
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000899  0.000686
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
 REFLECTING_RANGE_E.S.D.=   0.127
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

    1354 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      25 ...      25
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1354 REFLECTIONS ACCEPTED



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
 *  31        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  35        mP          0.0      77.5   38.8   77.5  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  25        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  34        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  23        oC          0.0     109.6  109.6   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  32        oP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  20        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
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
    41        mC        998.2     173.3   77.5   38.8  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.2      77.5  173.3   38.8  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.2      77.5  173.3   38.8  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.2      77.5   77.5   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.3  116.3  116.3  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.7   86.7  109.6  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.6   86.7   86.7  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.7   86.7  109.6  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.7   86.7  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.8   86.7  248.1 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
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

  *   0.00        0    4.4      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.850   -23.481     0.037
 UNIT_CELL_B-AXIS=    23.481    73.850     0.016
 UNIT_CELL_C-AXIS=    -0.020    -0.002    38.781
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       211 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.25
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.455    77.455    38.761  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  6.9E-02 6.9E-02 3.4E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012911  0.012911  0.025799  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.819   -23.453     0.031
 COORDINATES OF UNIT CELL B-AXIS    23.453    73.819     0.011
 COORDINATES OF UNIT CELL C-AXIS    -0.016    -0.001    38.761
 CRYSTAL MOSAICITY (DEGREES)     0.127
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006633  0.001785
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.92   1232.90
 DETECTOR ORIGIN (PIXELS) AT                     1261.32   1242.55
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.26
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.46    77.46    38.76  90.000  90.000  90.000


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
  0.05  0.15        62        113.307       1.238
  0.15  0.25       170        154.648       0.907
  0.25  0.35       189        126.388       1.110
  0.35  0.45       150         86.467       1.623
  0.45  0.55       170        193.512       0.725
  0.55  0.65       195        128.477       1.092
  0.65  0.75       194        155.665       0.901
  0.75  0.85       154        159.587       0.879
  0.85  0.95        70        101.229       1.386
 
 -0.05  0.95      1354        140.313       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               234.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.016
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    24.1 XMAX=    24.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59449 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        502


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               234.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.015
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     9.9 XMAX=  2454.2 NXBIN=    3
 YMIN=     7.6 YMAX=  2520.2 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        502


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               234.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.011
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    24.1 XMAX=    24.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1647    1691
 DETECTOR_SURFACE_POSITION=     817    1691
 DETECTOR_SURFACE_POSITION=     817     837
 DETECTOR_SURFACE_POSITION=    1647     837
 DETECTOR_SURFACE_POSITION=    2172    1664
 DETECTOR_SURFACE_POSITION=    1621    2230
 DETECTOR_SURFACE_POSITION=     843    2230
 DETECTOR_SURFACE_POSITION=     293    1664
 DETECTOR_SURFACE_POSITION=     293     864
 DETECTOR_SURFACE_POSITION=     843     298
 DETECTOR_SURFACE_POSITION=    1621     298
 DETECTOR_SURFACE_POSITION=    2172     864
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        502


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
 8.566E-01  1.143E-03   31.96


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

   24.493   3.155    24.23   1.79      2.39      1.95      66      90      28
    3.155   2.240     9.73   1.31      5.80      4.76     186     217       6
    2.240   1.832     1.45   0.85     29.76     32.02     196     230       2
    1.832   1.587    -0.00   0.95   3171.39   3263.02     240     277       2
    1.587   1.420     0.04   0.73   2519.57   2823.49     272     322       2
    1.420   1.297    -0.09   0.64   1625.97   2018.00     128     175       2
   --------------------------------------------------------------------------
   24.493   1.297     3.53   0.95      7.37      7.21    1088    1311      42
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      25      25)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.684   0.2265E+01  0.1201E+02         0.19           1 
   47    0    0     1.648  -0.7238E+00  0.1209E+02        -0.06           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -0.5% OF MEAN INTENSITY


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

     3.85          56      35      1285        2.7%       2.2%      1.9%       42   32.00      3.1%    99.9*     0    0.000       0
     2.74          96      57      2098        2.7%       4.1%      3.0%       78   21.64      5.8%    99.4*     0    0.000       0
     2.24         155      89      2665        3.3%       7.1%      6.0%      132   10.94     10.1%    99.0*     0    0.000       0
     1.94         149      85      3101        2.7%      20.9%     21.8%      128    2.88     29.5%    93.8*     0    0.000       0
     1.74         162      92      3479        2.6%    4380.8%   4980.6%      140   -0.03   6195.4%    16.8      0    0.000       0
     1.59         196     112      3851        2.9%    1959.8%   1950.2%      168    0.01   2771.6%     3.6      0    0.000       0
     1.47         208     121      4143        2.9%     544.9%    618.3%      174    0.12    770.6%    10.6      0    0.000       0
     1.38         203     117      4427        2.6%     -99.9%    -99.9%      172    0.00    -99.9%    15.6      0    0.000       0
     1.30          84      57      4743        1.2%     214.6%    241.1%       54    0.14    303.6%    60.7*     0    0.000       0
    total        1309     765     29792        2.6%       7.4%      7.2%     1088    4.68     10.4%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1354
 NUMBER OF REJECTED MISFITS                                42
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1311
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    767



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

   WILSON LINE (using all data) : A=   7.358 B=  32.809 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      14     7.652   0.004  1.9613E+03       7.581     -26.2
      17     4.653   0.012  2.1198E+03       7.659     -13.1
      23     3.511   0.020  8.6510E+02       6.763      14.7
      30     2.998   0.028  2.7617E+02       5.621      31.2
      35     2.645   0.036  2.4900E+02       5.517      25.8
      48     2.401   0.043  1.2297E+02       4.812      29.3
      36     2.211   0.051  6.5533E+01       4.183      31.1
      44     2.058   0.059  2.3483E+01       3.156      35.6
      40     1.935   0.067  4.5963E-01      -0.777      60.9
      50     1.831   0.075 -1.1439E+00    -999.999   -1000.0
      45     1.741   0.082  9.9278E-01      -0.007      44.6
      49     1.668   0.090 -3.9179E-01    -999.999   -1000.0
      51     1.599   0.098  2.5151E+00       0.922      32.9
      60     1.538   0.106  1.5902E+00       0.464      32.6
      69     1.482   0.114  1.0378E+00       0.037      32.1
      52     1.435   0.121 -2.4471E+00    -999.999   -1000.0
      46     1.392   0.129 -1.6553E+00    -999.999   -1000.0
      39     1.353   0.137  2.4111E+00       0.880      23.7
      17     1.316   0.144 -7.3841E-01    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.652        1.897        1.373        0.777
       3     4.653        1.406        0.911        0.438
       2     3.511        0.407        0.090        0.014
       3     2.998        1.585        0.742        0.258
       2     2.645        0.287        0.075        0.014
       5     2.401        0.881        0.481        0.208
       1     2.211        0.339        0.068        0.010
       7     2.058        0.520        0.210        0.065
       7     1.741       32.657       19.775      219.728
      11     1.668       31.906       23.741      200.016
       6     1.482       47.637      119.824      319.821
      77   overall       27.689      -75.769      846.047


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.652        0.534        0.389        0.224
      14     4.653        1.240        1.558        1.554
      21     3.511        0.869        0.766        0.633
      27     2.998        0.860        0.763        0.581
      33     2.645        1.208        1.415        1.509
      43     2.401        0.926        0.801        0.586
      35     2.211        0.880        0.685        0.473
      37     2.058        1.315        1.718        2.162
      34     1.935      194.131     1177.282    18786.378
      38     1.741       41.218      104.005     1085.429
      48     1.599        5.949        4.969       13.746
      57     1.538       12.871       17.916       80.786
      63     1.482       42.471      -54.957     1021.711
      38     1.353       21.799      132.198     1000.956
     688   overall       18.769       68.360     1145.169

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.652  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.38  0.73
       3     4.653  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.511  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.998  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.645  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.401  0.00  0.00  0.00  0.00  0.38  0.71  0.67  0.95  0.91  0.88
       1     2.211  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       7     2.058  1.73  1.24  1.37  1.21  1.10  1.02  0.96  0.91  0.87  0.84
       7     1.741  1.73  1.24  1.03  0.91  0.82  0.76  0.72  0.68  0.65  0.63
      11     1.668  2.93  2.10  1.75  1.54  1.40  1.30  1.22  1.16  1.11  1.07
       6     1.482  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      77   overall  0.94  0.68  0.59  0.52  0.50  0.49  0.46  0.48  0.45  0.48


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.652  2.10  3.31  2.31  1.82  1.78  1.55  1.59  1.45  1.35  1.27
      14     4.653  1.50  1.97  1.38  1.30  1.09  1.58  1.42  1.30  1.44  1.36
      21     3.511  0.50  0.53  0.55  0.43  0.61  0.74  1.04  1.04  1.04  1.05
      27     2.998  1.56  1.63  1.29  1.35  1.32  1.31  1.25  1.28  1.19  1.11
      33     2.645  0.64  1.50  1.17  1.19  1.00  1.01  1.02  1.16  1.17  1.10
      43     2.401  0.49  0.90  0.90  1.06  1.18  1.13  1.02  1.06  1.14  1.07
      35     2.211  1.50  0.79  0.66  0.87  0.87  0.76  0.85  0.93  0.96  0.90
      37     2.058  1.42  0.89  0.73  0.90  1.03  0.90  0.97  0.98  1.05  1.03
      34     1.935  4.64  2.43  1.70  1.34  1.12  0.98  0.88  0.80  0.74  0.70
      38     1.741  4.42  2.32  1.73  1.36  1.14  0.99  0.94  0.86  0.80  0.75
      48     1.599  4.16  2.18  1.53  1.39  1.16  1.02  0.95  0.87  0.81  0.76
      57     1.538  4.79  2.61  1.83  1.54  1.29  1.13  1.01  0.96  0.89  0.83
      63     1.482  4.84  2.54  1.84  1.44  1.21  1.06  0.95  0.89  0.86  0.83
      38     1.353  4.98  2.61  2.13  1.76  1.47  1.28  1.15  1.05  0.98  0.96
     688   overall  2.23  1.38  1.04  0.93  0.84  0.77  0.74  0.72  0.70  0.67


 List of     59 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   29   20   10    1.91   51.39  0.2362E+02  0.1320E+02 "alien"
   21   17   14    1.99   42.08  0.1934E+02  0.1084E+02 "alien"
   40    2   12    1.66   41.27 -0.1617E+02  0.7928E+01 "alien"
   50    8   15    1.32   38.61 -0.2851E+02  0.1631E+02 "alien"
   41   11   21    1.30   36.52 -0.2697E+02  0.1287E+02 "alien"
   32   25   12    1.64   35.58 -0.1394E+02  0.7998E+01 "alien"
   29    1   13    1.99   34.23  0.1573E+02  0.7965E+01 "alien"
   34   26   10    1.64   32.90 -0.1289E+02  0.8412E+01 "alien"
   36    3   14    1.70   32.71 -0.1281E+02  0.1188E+02 "alien"
   29    2   14    1.92   32.60  0.1498E+02  0.7719E+01 "alien"
   21   19   19    1.64   32.43 -0.1270E+02  0.9686E+01 "alien"
   54    5    9    1.36   30.50  0.7355E+02  0.1008E+02 "alien"
   32   23   10    1.75   28.36  0.2815E+02  0.1494E+02 "alien"
   36   32   14    1.39   27.46 -0.4545E+02  0.1629E+02 "alien"
   39    2   12    1.69   27.38 -0.1073E+02  0.7626E+01 "alien"
   36   26    7    1.66   23.11 -0.9054E+01  0.8045E+01 "alien"
   38    2    4    1.99   21.64  0.9946E+01  0.1032E+02 "alien"
   39    2    1    1.98   21.01  0.9656E+01  0.8814E+01 "alien"
   25   19   13    1.90   19.21  0.8831E+01  0.7644E+01 "alien"
   39    1   11    1.73   18.44  0.1831E+02  0.1320E+02 "alien"
   45    1    7    1.64   18.12 -0.7098E+01  0.8142E+01 "alien"
   52    2    1    1.49   17.44  0.1810E+02  0.9078E+01 "alien"
   23   22   17    1.66   17.42 -0.6824E+01  0.7106E+01 "alien"
   22   21   17    1.70   17.18 -0.6729E+01  0.9386E+01 "alien"
   34    4   16    1.65   16.46 -0.6447E+01  0.7556E+01 "alien"
   42   10   20    1.32   16.43 -0.1213E+02  0.1268E+02 "alien"
   25   20   21    1.47   16.18  0.1680E+02  0.1067E+02 "alien"
   34   22    1    1.91   16.12  0.7407E+01  0.8613E+01 "alien"
   27   19   11    1.95   15.33  0.7044E+01  0.8229E+01 "alien"
   25   18   12    1.98   15.31  0.7038E+01  0.1200E+02 "alien"
   33   29   14    1.49   15.31  0.1588E+02  0.9981E+01 "alien"
   51    8   15    1.30   15.13 -0.1117E+02  0.1363E+02 "alien"
   43    1    8    1.69   14.91 -0.5843E+01  0.1548E+02 "alien"
   38   26    0    1.68   14.83 -0.1626E+02  0.1502E+02 "alien"
   52    7   13    1.32   14.80 -0.1093E+02  0.9906E+01 "alien"
   26   24   19    1.49   14.48  0.1503E+02  0.1278E+02 "alien"
   44    3   12    1.54   14.19  0.2256E+02  0.8852E+01 "alien"
   37   30   11    1.48   14.00  0.1453E+02  0.1085E+02 "alien"
   34   22    4    1.88   13.28 -0.1519E+02  0.9558E+01 "alien"
   35    0   10    1.92   13.22  0.1517E+02  0.7992E+01 "alien"
   38    9   20    1.39   13.13 -0.2173E+02  0.8769E+01 "alien"
   27   21   13    1.80   12.91 -0.1477E+02  0.7073E+01 "alien"
   43    1    3    1.78   12.70  0.1261E+02  0.8586E+01 "alien"
   23   13   23    1.46   12.51  0.1299E+02  0.8919E+01 "alien"
   48    8   16    1.33   12.20 -0.9012E+01  0.9718E+01 "alien"
   51    6   13    1.35   12.14  0.2926E+02  0.1609E+02 "alien"
   38    1    6    1.94   11.96  0.5499E+01  0.9743E+01 "alien"
   26   23   20    1.46   11.46  0.1190E+02  0.1253E+02 "alien"
   26   19   12    1.93   11.33  0.5206E+01  0.7755E+01 "alien"
   46    0    2    1.68   10.78 -0.8955E+01  0.8400E+01 "alien"
   31    9   21    1.46   10.77  0.1118E+02  0.1191E+02 "alien"
   37   30   10    1.50   10.68  0.1109E+02  0.1040E+02 "alien"
   37   33   13    1.38   10.68 -0.1768E+02  0.9303E+01 "alien"
   37   27    6    1.64   10.65 -0.4174E+01  0.1173E+02 "alien"
   33   25   11    1.65   10.47 -0.4103E+01  0.1292E+02 "alien"
   42   32    1    1.47   10.37  0.1076E+02  0.9490E+01 "alien"
   41   31    3    1.50   10.27  0.1066E+02  0.9171E+01 "alien"
   39   29    5    1.56   10.17  0.1616E+02  0.8926E+01 "alien"
   45    9   18    1.33   10.15 -0.7493E+01  0.9992E+01 "alien"


 List of     59 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   51    8   15    1.30   15.13 -0.1117E+02  0.1363E+02
   41   11   21    1.30   36.52 -0.2697E+02  0.1287E+02
   50    8   15    1.32   38.61 -0.2851E+02  0.1631E+02
   42   10   20    1.32   16.43 -0.1213E+02  0.1268E+02
   52    7   13    1.32   14.80 -0.1093E+02  0.9906E+01
   45    9   18    1.33   10.15 -0.7493E+01  0.9992E+01
   48    8   16    1.33   12.20 -0.9012E+01  0.9718E+01
   51    6   13    1.35   12.14  0.2926E+02  0.1609E+02
   54    5    9    1.36   30.50  0.7355E+02  0.1008E+02
   37   33   13    1.38   10.68 -0.1768E+02  0.9303E+01
   38    9   20    1.39   13.13 -0.2173E+02  0.8769E+01
   36   32   14    1.39   27.46 -0.4545E+02  0.1629E+02
   23   13   23    1.46   12.51  0.1299E+02  0.8919E+01
   31    9   21    1.46   10.77  0.1118E+02  0.1191E+02
   26   23   20    1.46   11.46  0.1190E+02  0.1253E+02
   42   32    1    1.47   10.37  0.1076E+02  0.9490E+01
   25   20   21    1.47   16.18  0.1680E+02  0.1067E+02
   37   30   11    1.48   14.00  0.1453E+02  0.1085E+02
   33   29   14    1.49   15.31  0.1588E+02  0.9981E+01
   52    2    1    1.49   17.44  0.1810E+02  0.9078E+01
   26   24   19    1.49   14.48  0.1503E+02  0.1278E+02
   41   31    3    1.50   10.27  0.1066E+02  0.9171E+01
   37   30   10    1.50   10.68  0.1109E+02  0.1040E+02
   44    3   12    1.54   14.19  0.2256E+02  0.8852E+01
   39   29    5    1.56   10.17  0.1616E+02  0.8926E+01
   21   19   19    1.64   32.43 -0.1270E+02  0.9686E+01
   37   27    6    1.64   10.65 -0.4174E+01  0.1173E+02
   34   26   10    1.64   32.90 -0.1289E+02  0.8412E+01
   32   25   12    1.64   35.58 -0.1394E+02  0.7998E+01
   45    1    7    1.64   18.12 -0.7098E+01  0.8142E+01
   33   25   11    1.65   10.47 -0.4103E+01  0.1292E+02
   34    4   16    1.65   16.46 -0.6447E+01  0.7556E+01
   40    2   12    1.66   41.27 -0.1617E+02  0.7928E+01
   36   26    7    1.66   23.11 -0.9054E+01  0.8045E+01
   23   22   17    1.66   17.42 -0.6824E+01  0.7106E+01
   46    0    2    1.68   10.78 -0.8955E+01  0.8400E+01
   38   26    0    1.68   14.83 -0.1626E+02  0.1502E+02
   43    1    8    1.69   14.91 -0.5843E+01  0.1548E+02
   39    2   12    1.69   27.38 -0.1073E+02  0.7626E+01
   36    3   14    1.70   32.71 -0.1281E+02  0.1188E+02
   22   21   17    1.70   17.18 -0.6729E+01  0.9386E+01
   39    1   11    1.73   18.44  0.1831E+02  0.1320E+02
   32   23   10    1.75   28.36  0.2815E+02  0.1494E+02
   43    1    3    1.78   12.70  0.1261E+02  0.8586E+01
   27   21   13    1.80   12.91 -0.1477E+02  0.7073E+01
   34   22    4    1.88   13.28 -0.1519E+02  0.9558E+01
   25   19   13    1.90   19.21  0.8831E+01  0.7644E+01
   34   22    1    1.91   16.12  0.7407E+01  0.8613E+01
   29   20   10    1.91   51.39  0.2362E+02  0.1320E+02
   29    2   14    1.92   32.60  0.1498E+02  0.7719E+01
   35    0   10    1.92   13.22  0.1517E+02  0.7992E+01
   26   19   12    1.93   11.33  0.5206E+01  0.7755E+01
   38    1    6    1.94   11.96  0.5499E+01  0.9743E+01
   27   19   11    1.95   15.33  0.7044E+01  0.8229E+01
   39    2    1    1.98   21.01  0.9656E+01  0.8814E+01
   25   18   12    1.98   15.31  0.7038E+01  0.1200E+02
   29    1   13    1.99   34.23  0.1573E+02  0.7965E+01
   21   17   14    1.99   42.08  0.1934E+02  0.1084E+02
   38    2    4    1.99   21.64  0.9946E+01  0.1032E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        17
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1354
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1285
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       69


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
