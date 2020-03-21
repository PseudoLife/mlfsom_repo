
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.323    77.323    38.697  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=       7       7
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       7
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  1.000000 -0.000660  0.000714
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
 BEAM_DIVERGENCE_E.S.D.=   0.040
 REFLECTING_RANGE_E.S.D.=   0.087
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

    1496 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       7 ...       7
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1496 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               46
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
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    27        mC        500.0     159.4   38.7  109.4  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.7  109.4  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.4  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        998.2     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.2      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.2      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.2      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.4   86.5   86.5  78.4  50.8  50.8    1  1  0  0  1  0 -1  0  0  1 -1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    1  0 -1  0  0  1 -1  0  1  1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0 -1  1  0 -1  0  1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.5  247.6 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.5   86.5   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.2  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.5   86.5  177.2 107.0 119.2  53.2    0 -1  1  0  0 -1 -1  0  2  1 -1  0
    14        mC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2    0 -1 -1  0  0 -1  1  0 -1  0  1  0
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

  *   0.00        0    5.9      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.692   -23.416     0.027
 UNIT_CELL_B-AXIS=    23.416    73.692     0.010
 UNIT_CELL_C-AXIS=    -0.014    -0.000    38.696
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       331 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.21
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.285    77.285    38.679  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.9E-01 1.9E-01 8.8E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012939  0.012939  0.025854  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.654   -23.410     0.021
 COORDINATES OF UNIT CELL B-AXIS    23.410    73.654     0.007
 COORDINATES OF UNIT CELL C-AXIS    -0.011    -0.000    38.679
 CRYSTAL MOSAICITY (DEGREES)     0.087
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006357  0.001876
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.88   1232.90
 DETECTOR ORIGIN (PIXELS) AT                     1261.15   1242.15
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.17
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.28    77.28    38.68  90.000  90.000  90.000


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
  0.05  0.15       124        251.488       1.093
  0.15  0.25       173        243.839       1.128
  0.25  0.35       192        250.508       1.098
  0.35  0.45       151        162.372       1.693
  0.45  0.55       180        403.968       0.681
  0.55  0.65       197        227.406       1.209
  0.65  0.75       179        212.804       1.292
  0.75  0.85       169        235.704       1.167
  0.85  0.95       131        533.823       0.515
 
 -0.05  0.95      1496        274.969       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               233.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.018
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     6.1 XMAX=     6.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59707 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               233.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.015
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    12.2 XMAX=  2438.5 NXBIN=    3
 YMIN=     6.9 YMAX=  2519.3 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               233.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.007
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     6.1 XMAX=     6.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1225    1263
 DETECTOR_SURFACE_POSITION=    1637    1690
 DETECTOR_SURFACE_POSITION=     813    1690
 DETECTOR_SURFACE_POSITION=     813     836
 DETECTOR_SURFACE_POSITION=    1637     836
 DETECTOR_SURFACE_POSITION=    2158    1663
 DETECTOR_SURFACE_POSITION=    1612    2229
 DETECTOR_SURFACE_POSITION=     839    2229
 DETECTOR_SURFACE_POSITION=     293    1663
 DETECTOR_SURFACE_POSITION=     293     863
 DETECTOR_SURFACE_POSITION=     839     297
 DETECTOR_SURFACE_POSITION=    1612     297
 DETECTOR_SURFACE_POSITION=    2158     863
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
 1.371E+00  2.218E-03   18.13


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

   24.440   3.396    15.99   1.44      2.64      3.20      64      89      20
    3.396   2.413    12.84   1.56      4.91      3.71     130     169      14
    2.413   1.973     8.12   0.97      5.73      5.37     180     216       8
    1.973   1.710    -0.04   0.59    -99.90    -99.90     190     244       0
    1.710   1.530     0.00   0.74    -99.90    -99.90     214     278       0
    1.530   1.398    -0.03   0.52    -99.90    -99.90     244     310       0
    1.398   1.294     0.04   0.54   1613.59   2030.18     100     145       2
   --------------------------------------------------------------------------
   24.440   1.294     3.68   0.82      5.57      5.88    1122    1451      44
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       7       7)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.680  -0.7371E+01  0.1305E+02        -0.56           1 
   47    0    0     1.644   0.9950E+01  0.1306E+02         0.76           1*

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

     3.84          67      43      1285        3.3%       2.6%      3.2%       48   19.92      3.7%    99.8*     0    0.000       0
     2.73         106      67      2098        3.2%       4.4%      3.5%       78   17.72      6.2%    99.2*     0    0.000       0
     2.24         161      97      2665        3.6%       4.9%      4.1%      128   14.50      6.9%    99.0*     0    0.000       0
     1.94         165      95      3101        3.1%       6.4%      6.8%      140    7.27      9.1%    98.8*     0    0.000       0
     1.73         188     115      3481        3.3%     -99.9%    -99.9%      146    0.00    -99.9%     0.7      0    0.000       0
     1.58         216     131      3848        3.4%     -99.9%    -99.9%      170    0.00    -99.9%   -19.2      0    0.000       0
     1.47         241     150      4143        3.6%     -99.9%    -99.9%      182    0.00    -99.9%    -2.7      0    0.000       0
     1.37         220     131      4429        3.0%    6050.8%   8130.2%      178    0.02   8557.1%     6.2      0    0.000       0
     1.29          87      61      4743        1.3%     -99.9%    -99.9%       52    0.00    -99.9%   -27.5      0    0.000       0
    total        1451     890     29793        3.0%       5.6%      5.9%     1122    4.64      7.9%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1496
 NUMBER OF REJECTED MISFITS                                44
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1451
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    890



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

   WILSON LINE (using all data) : A=   8.727 B=  41.673 CORRELATION=  0.89
      #      RES      SS        <I>       log(<I>)       BO
      18     7.395   0.005  3.0400E+03       8.020      77.4
      21     4.619   0.012  2.9760E+03       7.998      31.1
      30     3.493   0.020  1.3877E+03       7.235      36.4
      31     2.994   0.028  7.1462E+02       6.572      38.6
      40     2.646   0.036  7.4619E+02       6.615      29.6
      52     2.394   0.044  4.3291E+02       6.071      30.5
      42     2.204   0.051  3.2032E+02       5.769      28.7
      45     2.055   0.059  1.5030E+02       5.013      31.4
      50     1.930   0.067 -6.8987E-01    -999.999   -1000.0
      62     1.827   0.075 -2.8420E-01    -999.999   -1000.0
      54     1.738   0.083  6.5450E-01      -0.424      55.3
      58     1.663   0.090  3.5861E-01      -1.026      54.0
      61     1.595   0.098  2.9634E-01      -1.216      50.6
      77     1.533   0.106 -9.0733E-02    -999.999   -1000.0
      79     1.479   0.114  8.1332E-02      -2.509      49.2
      57     1.431   0.122  3.9655E-03      -5.530      58.4
      54     1.390   0.129  4.7640E-01      -0.742      36.6
      40     1.351   0.137  3.3344E-01      -1.098      35.8
      19     1.313   0.145  2.8250E+00       1.039      26.5


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.395        0.974        0.463        0.174
       3     4.619        1.148        0.662        0.285
       3     3.493        0.357        0.091        0.016
       4     2.994        0.997        0.423        0.136
       3     2.646        0.249        0.059        0.010
       5     2.394        0.872        0.426        0.156
       3     2.204        1.199        0.462        0.128
       8     2.055        0.253        0.083        0.021
      11     1.827       15.738        4.289       47.110
      11     1.738       73.415      347.693     1281.738
       6     1.479     9320.363   445282.998 23209610.056
       3     1.431  1021023.300408849117.595*************
      91   overall    34300.561 13507890.841*************


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      13     7.395        1.406        2.310        3.103
      18     4.619        1.095        1.229        1.175
      27     3.493        0.852        0.647        0.450
      27     2.994        0.859        0.754        0.599
      37     2.646        1.129        1.190        1.129
      47     2.394        0.896        0.711        0.474
      39     2.204        0.774        0.561        0.349
      37     2.055        1.382        1.694        1.903
      47     1.663      190.731      934.852    16475.858
      58     1.595      398.246    -4514.653   180779.100
      74     1.533       31.812       35.529      838.292
      54     1.390      106.567      273.231     4533.815
      38     1.351      251.900     1178.619    34261.258
      19     1.313       14.395       41.427      175.354
     799   overall   117801.604 -6622384.912*************

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.395  0.00  0.00  0.00  0.00  0.00  0.36  0.33  0.64  0.61  0.59
       3     4.619  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.493  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.994  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.646  1.34  0.96  0.80  0.70  0.64  0.59  0.56  1.06  1.01  0.98
       5     2.394  0.00  0.00  0.00  0.42  0.77  1.07  1.00  0.95  0.91  0.88
       3     2.204  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.055  1.51  1.08  1.50  1.32  1.20  1.11  1.05  1.19  1.14  1.10
      11     1.827  2.93  2.10  1.75  1.54  1.40  1.30  1.22  1.16  1.11  1.07
      11     1.738  0.37  0.26  0.22  0.19  0.17  0.16  0.15  0.14  0.14  0.13
       6     1.479  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       3     1.431  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      91   overall  0.80  0.57  0.53  0.49  0.46  0.49  0.46  0.49  0.47  0.47


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      13     7.395  2.42  2.55  1.78  2.10  1.95  1.70  1.53  1.40  1.30  1.22
      18     4.619  1.17  1.84  1.29  1.01  1.27  1.23  1.10  1.01  1.12  1.14
      27     3.493  0.78  0.61  0.57  0.67  0.66  0.74  0.88  0.94  0.94  0.88
      27     2.994  1.56  1.43  1.00  1.12  1.13  1.15  1.18  1.14  1.06  1.05
      37     2.646  0.85  1.19  1.25  1.15  1.03  0.96  1.07  1.03  1.05  1.03
      47     2.394  1.34  1.29  0.99  1.10  1.08  1.04  1.06  1.08  1.11  1.04
      39     2.204  0.81  1.27  0.99  0.78  0.85  0.97  0.92  0.98  1.12  1.05
      37     2.055  1.70  1.49  1.04  1.23  1.03  0.90  0.91  0.98  1.00  0.98
      47     1.663  5.59  2.93  2.05  1.61  1.35  1.18  1.06  0.97  0.90  0.84
      58     1.595  4.71  2.47  1.73  1.36  1.14  0.99  0.89  0.85  0.78  0.74
      74     1.533  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      54     1.390  4.86  2.66  1.86  1.46  1.22  1.07  0.96  0.87  0.84  0.79
      38     1.351  4.98  2.61  1.83  1.44  1.20  1.05  0.94  0.86  0.80  0.75
      19     1.313  5.53  2.90  2.03  1.60  1.34  1.17  1.05  0.96  0.89  0.83
     799   overall  2.24  1.39  1.01  0.87  0.77  0.71  0.67  0.65  0.63  0.60


 List of    200 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   36    9   20    1.42 4244.36  0.1683E+02  0.1243E+02 "alien"
   39    7   18    1.44 4213.98  0.1671E+02  0.9422E+01 "alien"
   26   19   22    1.42 3321.43  0.1317E+02  0.1087E+02 "alien"
   46    5   14    1.43 3141.33  0.1246E+02  0.1488E+02 "alien"
   49    5   12    1.41 2257.17  0.8951E+01  0.1007E+02 "alien"
   37   31   12    1.43 2009.52  0.7969E+01  0.1687E+02 "alien"
   51    4    9    1.43 1931.86  0.7661E+01  0.1044E+02 "alien"
   30   29   17    1.44 1930.53  0.7655E+01  0.9644E+01 "alien"
   54    3    1    1.43 1770.82  0.7022E+01  0.1023E+02 "alien"
   38    8   19    1.42 1761.54  0.6985E+01  0.8807E+01 "alien"
   37    8   19    1.44 1494.89  0.5928E+01  0.9885E+01 "alien"
   48    5   13    1.41 1389.08  0.5508E+01  0.1057E+02 "alien"
   36   31   13    1.43 1101.52  0.4368E+01  0.1051E+02 "alien"
   32   30   16    1.42  851.81  0.3378E+01  0.9731E+01 "alien"
   28   26   19    1.43  809.63  0.3211E+01  0.9389E+01 "alien"
   53    3    6    1.42  780.90  0.3097E+01  0.1272E+02 "alien"
   54    3    2    1.43  514.52  0.2040E+01  0.1048E+02 "alien"
   24   14   23    1.44  454.92  0.1804E+01  0.7066E+01 "alien"
   50    4   10    1.43  424.03  0.1681E+01  0.1016E+02 "alien"
   21   13   24    1.43  421.72  0.1672E+01  0.7637E+01 "alien"
   43    6   16    1.43  385.12  0.1527E+01  0.9562E+01 "alien"
   27   24   20    1.43  367.04  0.1455E+01  0.1270E+02 "alien"
   37   32   12    1.42  337.29  0.1338E+01  0.1060E+02 "alien"
   53    2    1    1.46  307.18  0.2498E+02  0.1346E+02 "alien"
   40   30    2    1.54  264.77 -0.2402E+02  0.1133E+02 "alien"
   53    3    4    1.44  226.12  0.8967E+00  0.1043E+02 "alien"
   37   31   11    1.46  197.07  0.1603E+02  0.1109E+02 "alien"
   52    3    6    1.45  178.01  0.7059E+00  0.1011E+02 "alien"
   51    2    1    1.51  167.27 -0.1518E+02  0.1893E+02 "alien"
   35   30   13    1.46  149.78  0.1218E+02  0.1009E+02 "alien"
   33    8   20    1.47  127.84  0.1040E+02  0.1233E+02 "alien"
   27   25   19    1.46  119.51  0.9720E+01  0.1008E+02 "alien"
   47    3   10    1.51  119.51 -0.1084E+02  0.9956E+01 "alien"
   39    5   16    1.53  117.77 -0.1069E+02  0.8633E+01 "alien"
   44    4   13    1.51  113.82 -0.1033E+02  0.9363E+01 "alien"
   37   29    9    1.54  113.15 -0.1027E+02  0.9688E+01 "alien"
   51    2    2    1.51  112.47 -0.1020E+02  0.1031E+02 "alien"
   41   32    5    1.46  110.96  0.9025E+01  0.1094E+02 "alien"
   31   28   15    1.50  109.78  0.8929E+01  0.1484E+02 "alien"
   25   23   19    1.52  109.17 -0.9905E+01  0.1118E+02 "alien"
   24   16   22    1.47  105.77  0.8603E+01  0.1144E+02 "alien"
   46    4   12    1.49  104.30  0.8483E+01  0.9660E+01 "alien"
   48    3    9    1.51   98.61  0.8020E+01  0.1552E+02 "alien"
   29    9   21    1.49   94.05  0.7649E+01  0.7326E+01 "alien"
   34    6   18    1.55   91.28 -0.8282E+01  0.8275E+01 "alien"
   44    3   12    1.54   91.16 -0.8272E+01  0.9426E+01 "alien"
   42   32    1    1.46   91.13  0.7412E+01  0.1126E+02 "alien"
   39   29    5    1.56   89.85 -0.8152E+01  0.1000E+02 "alien"
   40   31    6    1.49   88.70  0.7214E+01  0.1043E+02 "alien"
   36   30   12    1.47   87.86  0.7146E+01  0.1023E+02 "alien"
   33   28   13    1.53   84.62 -0.7678E+01  0.1522E+02 "alien"
   45    3   11    1.54   84.01 -0.7622E+01  0.9425E+01 "alien"
   44    5   14    1.48   82.63  0.6720E+01  0.1401E+02 "alien"
   23   17   21    1.52   82.44 -0.7480E+01  0.9987E+01 "alien"
   33   29   14    1.48   82.00  0.6669E+01  0.1014E+02 "alien"
   50    2    5    1.51   81.78 -0.7420E+01  0.9692E+01 "alien"
   54    3    4    1.41   74.18  0.2942E+00  0.1044E+02 "alien"
   37   30   10    1.50   73.37  0.5967E+01  0.1065E+02 "alien"
   41   31    2    1.50   73.35  0.5966E+01  0.1081E+02 "alien"
   32   29   15    1.47   72.15  0.5868E+01  0.1001E+02 "alien"
   30    9   21    1.48   71.54  0.5818E+01  0.1011E+02 "alien"
   39   30    7    1.51   68.30 -0.6197E+01  0.1025E+02 "alien"
   30    8   20    1.53   66.92 -0.6072E+01  0.7471E+01 "alien"
   41   31    4    1.49   66.15  0.5380E+01  0.1070E+02 "alien"
   21   10   22    1.55   62.23 -0.5646E+01  0.7928E+01 "alien"
   21   12   22    1.54   61.78 -0.5606E+01  0.8312E+01 "alien"
   37    7   18    1.48   61.66  0.5015E+01  0.1387E+02 "alien"
   42    5   15    1.49   60.31  0.4905E+01  0.9049E+01 "alien"
   36   28   10    1.55   58.56 -0.5313E+01  0.1469E+02 "alien"
   52    2    4    1.47   57.03  0.4638E+01  0.1800E+02 "alien"
   41   32    4    1.47   55.85  0.4542E+01  0.1096E+02 "alien"
   33   12   23    1.34   54.39  0.1814E+02  0.1549E+02 "alien"
   38   29    7    1.55   53.39 -0.4845E+01  0.9979E+01 "alien"
   41   31    3    1.49   53.32  0.4337E+01  0.1055E+02 "alien"
   27   21   13    1.80   50.91 -0.1447E+02  0.7964E+01 "alien"
   51    6   13    1.34   50.01  0.1668E+02  0.1715E+02 "alien"
   47    4   12    1.46   49.18  0.4000E+01  0.1018E+02 "alien"
   33   22    6    1.87   48.47 -0.1378E+02  0.9430E+01 "alien"
   27   24   16    1.60   48.46  0.1436E+02  0.1426E+02 "alien"
   37   29   10    1.51   47.80 -0.4337E+01  0.1606E+02 "alien"
   43    1    2    1.79   47.03 -0.1337E+02  0.9179E+01 "alien"
   35   23    3    1.83   46.80 -0.1330E+02  0.9497E+01 "alien"
   29   29   18    1.42   45.93  0.8819E+01  0.1574E+02 "alien"
   37   34   14    1.34   45.58  0.1520E+02  0.1511E+02 "alien"
   39   32   10    1.42   45.54  0.1806E+00  0.1092E+02 "alien"
   54    3    0    1.43   44.95  0.7804E+01  0.1472E+02 "alien"
   31    6   18    1.61   44.32  0.1313E+02  0.1189E+02 "alien"
   35   25    8    1.68   44.21  0.1586E+02  0.1425E+02 "alien"
   39   29    4    1.57   43.81  0.1298E+02  0.1065E+02 "alien"
   31   26   14    1.57   43.67  0.1294E+02  0.9621E+01 "alien"
   33   25   11    1.65   42.10  0.1510E+02  0.1432E+02 "alien"
   49    6   14    1.36   41.97  0.1399E+02  0.1772E+02 "alien"
   26   23   16    1.64   41.40  0.1485E+02  0.1516E+02 "alien"
   37   27    7    1.61   40.38  0.1197E+02  0.1021E+02 "alien"
   26   25   17    1.56   39.17 -0.3554E+01  0.9050E+01 "alien"
   37    5   16    1.57   37.88  0.1123E+02  0.1390E+02 "alien"
   50    2    6    1.50   37.22  0.3027E+01  0.9853E+01 "alien"
   41    1    4    1.85   36.70 -0.1043E+02  0.9323E+01 "alien"
   43    4   13    1.53   36.07 -0.3273E+01  0.1576E+02 "alien"
   22    9   21    1.60   35.68  0.1057E+02  0.8974E+01 "alien"
   34    1   12    1.86   35.45 -0.1007E+02  0.8031E+01 "alien"
   36   26    7    1.66   34.80  0.1248E+02  0.9720E+01 "alien"
   30    3   15    1.82   33.93 -0.9643E+01  0.1074E+02 "alien"
   27   20   12    1.87   33.81 -0.9609E+01  0.8472E+01 "alien"
   28   21   22    1.38   33.68  0.1605E+02  0.8796E+01 "alien"
   36   26    6    1.68   33.31  0.1194E+02  0.1519E+02 "alien"
   43    3   12    1.57   32.49  0.9628E+01  0.9314E+01 "alien"
   50    5   11    1.41   32.17  0.1532E+02  0.1703E+02 "alien"
   32   23    9    1.78   32.10 -0.9124E+01  0.9053E+01 "alien"
   34    4   16    1.65   32.05  0.1149E+02  0.8691E+01 "alien"
   53    2    0    1.46   31.03  0.2786E+02  0.1794E+02 "alien"
   32   27   14    1.53   30.26 -0.2745E+01  0.1320E+02 "alien"
   25   10   22    1.50   30.14  0.2451E+01  0.8765E+01 "alien"
   20    5   24    1.48   30.05  0.2444E+01  0.6475E+01 "alien"
   22   19   15    1.85   28.88 -0.8207E+01  0.1319E+02 "alien"
   49    2    7    1.52   28.75 -0.2609E+01  0.9893E+01 "alien"
   42    1    9    1.69   26.63  0.9548E+01  0.9159E+01 "alien"
   21   20   18    1.67   26.56  0.9524E+01  0.1161E+02 "alien"
   39   34   11    1.37   26.34  0.1255E+02  0.1176E+02 "alien"
   42    1    1    1.84   26.18 -0.7441E+01  0.9771E+01 "alien"
   45    7   16    1.39   25.19  0.1200E+02  0.1038E+02 "alien"
   24   20   20    1.52   24.94 -0.2263E+01  0.1296E+02 "alien"
   26   10   22    1.49   24.90  0.2025E+01  0.9310E+01 "alien"
   46    8   17    1.34   24.25  0.8084E+01  0.1043E+02 "alien"
   38   33   12    1.39   24.09  0.1148E+02  0.1093E+02 "alien"
   47    1    6    1.59   24.08  0.7137E+01  0.1377E+02 "alien"
   42   32    0    1.46   23.66  0.1091E+02  0.1588E+02 "alien"
   39   28    3    1.60   23.31  0.6908E+01  0.1012E+02 "alien"
   41    8   18    1.40   23.09  0.1100E+02  0.1670E+02 "alien"
   31   29   19    1.36   22.66  0.7555E+01  0.1020E+02 "alien"
   33   11   22    1.38   22.29  0.1062E+02  0.8831E+01 "alien"
   24   24   18    1.56   22.28 -0.1018E+02  0.8920E+01 "alien"
   31    2   14    1.85   22.26 -0.6326E+01  0.7361E+01 "alien"
   30    1   13    1.95   21.91 -0.1512E+02  0.1204E+02 "alien"
   38   28    5    1.60   21.25  0.6298E+01  0.9889E+01 "alien"
   20   13   21    1.60   21.13  0.6263E+01  0.8284E+01 "alien"
   24   23   17    1.63   21.07  0.6244E+01  0.8159E+01 "alien"
   41   31    0    1.50   20.60  0.7340E+01  0.1487E+02 "alien"
   25   25   18    1.53   20.54 -0.8122E+01  0.8757E+01 "alien"
   42    4   14    1.53   20.26 -0.1838E+01  0.8952E+01 "alien"
   33    1   12    1.89   19.82 -0.1368E+02  0.8169E+01 "alien"
   31   24   12    1.68   19.29  0.6917E+01  0.9028E+01 "alien"
   37    1   11    1.80   18.79 -0.5341E+01  0.8352E+01 "alien"
   25   18   12    1.98   18.52 -0.1278E+02  0.1248E+02 "alien"
   37   25    1    1.73   18.19  0.1190E+02  0.9884E+01 "alien"
   36   24    3    1.77   18.18  0.1190E+02  0.1062E+02 "alien"
   33   26   12    1.60   17.95  0.5320E+01  0.1324E+02 "alien"
   41   10   20    1.33   17.62  0.5874E+01  0.9389E+01 "alien"
   39    2   12    1.69   17.55  0.6295E+01  0.8619E+01 "alien"
   31   22   10    1.80   17.53 -0.4983E+01  0.1352E+02 "alien"
   44    2   10    1.60   17.39  0.5153E+01  0.9448E+01 "alien"
   30   29   18    1.40   17.19  0.8189E+01  0.1480E+02 "alien"
   38    2    4    1.99   16.98 -0.1172E+02  0.1136E+02 "alien"
   35    1   12    1.82   16.72 -0.4752E+01  0.8560E+01 "alien"
   24   10   22    1.51   16.71 -0.1517E+01  0.8148E+01 "alien"
   33   22    5    1.89   16.62 -0.1147E+02  0.1057E+02 "alien"
   24   13   24    1.40   16.55  0.7885E+01  0.8800E+01 "alien"
   42    9   19    1.35   16.47  0.5491E+01  0.9618E+01 "alien"
   53    7   13    1.30   16.23  0.4584E+02  0.1475E+02 "alien"
   39    1    6    1.89   16.00 -0.1104E+02  0.9010E+01 "alien"
   29    2   14    1.92   16.00 -0.1104E+02  0.7699E+01 "alien"
   54    4    6    1.39   15.65  0.7458E+01  0.1093E+02 "alien"
   43    1    0    1.80   15.04 -0.1213E+02  0.1300E+02 "alien"
   33    6   18    1.57   14.99  0.4443E+01  0.7631E+01 "alien"
   34    5   17    1.60   14.72  0.4362E+01  0.8135E+01 "alien"
   34   23    6    1.81   14.50 -0.4122E+01  0.9303E+01 "alien"
   15    7   19    1.87   14.47 -0.4114E+01  0.5581E+01 "alien"
   47    1    5    1.61   14.38  0.4261E+01  0.9413E+01 "alien"
   28   19   23    1.35   14.05  0.4684E+01  0.1123E+02 "alien"
   37    2   13    1.71   13.97  0.9143E+01  0.9555E+01 "alien"
   20    8   23    1.52   13.78 -0.1251E+01  0.7793E+01 "alien"
   21   17   20    1.60   13.75  0.4076E+01  0.9146E+01 "alien"
   39    4   15    1.57   13.35  0.3955E+01  0.1214E+02 "alien"
   30   28   16    1.49   13.21  0.1074E+01  0.1019E+02 "alien"
   32   24   11    1.69   12.94  0.4639E+01  0.9583E+01 "alien"
   40   33    9    1.41   12.92  0.6157E+01  0.1086E+02 "alien"
   29   20    9    1.95   12.92 -0.8913E+01  0.1397E+02 "alien"
   32   23   10    1.75   12.63  0.8263E+01  0.1608E+02 "alien"
   50    6   13    1.36   12.62  0.4209E+01  0.1039E+02 "alien"
   35   23    2    1.84   12.57 -0.3573E+01  0.9911E+01 "alien"
   45    0    2    1.71   12.46  0.1931E+02  0.1756E+02 "alien"
   53    4    8    1.39   12.15  0.5789E+01  0.1069E+02 "alien"
   54    5    8    1.37   11.90  0.3968E+01  0.1219E+02 "alien"
   37    4   15    1.62   11.75  0.3482E+01  0.8940E+01 "alien"
   44    1    0    1.76   11.61  0.1700E+02  0.1405E+02 "alien"
   17   17   17    1.86   11.57 -0.7340E+01  0.8688E+01 "alien"
   40    2   12    1.66   11.38  0.4082E+01  0.8759E+01 "alien"
   17   15   25    1.41   11.38  0.5422E+01  0.7353E+01 "alien"
   39   28    2    1.60   11.22  0.3324E+01  0.9830E+01 "alien"
   27   23   15    1.66   11.16  0.4001E+01  0.8779E+01 "alien"
   21    9   21    1.62   11.11  0.3293E+01  0.7725E+01 "alien"
   29   26   20    1.38   10.99  0.5237E+01  0.1628E+02 "alien"
   31   21    8    1.90   10.80 -0.7452E+01  0.9144E+01 "alien"
   46    0    3    1.67   10.67  0.8059E+01  0.9541E+01 "alien"
   42    1    4    1.81   10.67 -0.3032E+01  0.9219E+01 "alien"
   36    3   14    1.69   10.55  0.3784E+01  0.9283E+01 "alien"
   35   24    6    1.75   10.52  0.6886E+01  0.9764E+01 "alien"
   43    0    7    1.71   10.51  0.1433E+02  0.1464E+02 "alien"
   43    2   11    1.60   10.47  0.3104E+01  0.1287E+02 "alien"
   39   28    1    1.61   10.42  0.3089E+01  0.1004E+02 "alien"


 List of    200 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   16.23  0.4584E+02  0.1475E+02
   41   10   20    1.33   17.62  0.5874E+01  0.9389E+01
   33   12   23    1.34   54.39  0.1814E+02  0.1549E+02
   46    8   17    1.34   24.25  0.8084E+01  0.1043E+02
   51    6   13    1.34   50.01  0.1668E+02  0.1715E+02
   37   34   14    1.34   45.58  0.1520E+02  0.1511E+02
   42    9   19    1.35   16.47  0.5491E+01  0.9618E+01
   28   19   23    1.35   14.05  0.4684E+01  0.1123E+02
   31   29   19    1.36   22.66  0.7555E+01  0.1020E+02
   49    6   14    1.36   41.97  0.1399E+02  0.1772E+02
   50    6   13    1.36   12.62  0.4209E+01  0.1039E+02
   54    5    8    1.37   11.90  0.3968E+01  0.1219E+02
   39   34   11    1.37   26.34  0.1255E+02  0.1176E+02
   28   21   22    1.38   33.68  0.1605E+02  0.8796E+01
   33   11   22    1.38   22.29  0.1062E+02  0.8831E+01
   29   26   20    1.38   10.99  0.5237E+01  0.1628E+02
   38   33   12    1.39   24.09  0.1148E+02  0.1093E+02
   45    7   16    1.39   25.19  0.1200E+02  0.1038E+02
   53    4    8    1.39   12.15  0.5789E+01  0.1069E+02
   54    4    6    1.39   15.65  0.7458E+01  0.1093E+02
   24   13   24    1.40   16.55  0.7885E+01  0.8800E+01
   41    8   18    1.40   23.09  0.1100E+02  0.1670E+02
   30   29   18    1.40   17.19  0.8189E+01  0.1480E+02
   40   33    9    1.41   12.92  0.6157E+01  0.1086E+02
   17   15   25    1.41   11.38  0.5422E+01  0.7353E+01
   50    5   11    1.41   32.17  0.1532E+02  0.1703E+02
   48    5   13    1.41 1389.08  0.5508E+01  0.1057E+02
   49    5   12    1.41 2257.17  0.8951E+01  0.1007E+02
   54    3    4    1.41   74.18  0.2942E+00  0.1044E+02
   29   29   18    1.42   45.93  0.8819E+01  0.1574E+02
   36    9   20    1.42 4244.36  0.1683E+02  0.1243E+02
   26   19   22    1.42 3321.43  0.1317E+02  0.1087E+02
   37   32   12    1.42  337.29  0.1338E+01  0.1060E+02
   53    3    6    1.42  780.90  0.3097E+01  0.1272E+02
   38    8   19    1.42 1761.54  0.6985E+01  0.8807E+01
   32   30   16    1.42  851.81  0.3378E+01  0.9731E+01
   39   32   10    1.42   45.54  0.1806E+00  0.1092E+02
   54    3    2    1.43  514.52  0.2040E+01  0.1048E+02
   51    4    9    1.43 1931.86  0.7661E+01  0.1044E+02
   36   31   13    1.43 1101.52  0.4368E+01  0.1051E+02
   54    3    1    1.43 1770.82  0.7022E+01  0.1023E+02
   54    3    0    1.43   44.95  0.7804E+01  0.1472E+02
   46    5   14    1.43 3141.33  0.1246E+02  0.1488E+02
   50    4   10    1.43  424.03  0.1681E+01  0.1016E+02
   21   13   24    1.43  421.72  0.1672E+01  0.7637E+01
   43    6   16    1.43  385.12  0.1527E+01  0.9562E+01
   37   31   12    1.43 2009.52  0.7969E+01  0.1687E+02
   27   24   20    1.43  367.04  0.1455E+01  0.1270E+02
   28   26   19    1.43  809.63  0.3211E+01  0.9389E+01
   30   29   17    1.44 1930.53  0.7655E+01  0.9644E+01
   24   14   23    1.44  454.92  0.1804E+01  0.7066E+01
   53    3    4    1.44  226.12  0.8967E+00  0.1043E+02
   37    8   19    1.44 1494.89  0.5928E+01  0.9885E+01
   39    7   18    1.44 4213.98  0.1671E+02  0.9422E+01
   52    3    6    1.45  178.01  0.7059E+00  0.1011E+02
   53    2    1    1.46  307.18  0.2498E+02  0.1346E+02
   37   31   11    1.46  197.07  0.1603E+02  0.1109E+02
   53    2    0    1.46   31.03  0.2786E+02  0.1794E+02
   41   32    5    1.46  110.96  0.9025E+01  0.1094E+02
   47    4   12    1.46   49.18  0.4000E+01  0.1018E+02
   35   30   13    1.46  149.78  0.1218E+02  0.1009E+02
   27   25   19    1.46  119.51  0.9720E+01  0.1008E+02
   42   32    1    1.46   91.13  0.7412E+01  0.1126E+02
   42   32    0    1.46   23.66  0.1091E+02  0.1588E+02
   52    2    4    1.47   57.03  0.4638E+01  0.1800E+02
   36   30   12    1.47   87.86  0.7146E+01  0.1023E+02
   41   32    4    1.47   55.85  0.4542E+01  0.1096E+02
   24   16   22    1.47  105.77  0.8603E+01  0.1144E+02
   32   29   15    1.47   72.15  0.5868E+01  0.1001E+02
   33    8   20    1.47  127.84  0.1040E+02  0.1233E+02
   44    5   14    1.48   82.63  0.6720E+01  0.1401E+02
   30    9   21    1.48   71.54  0.5818E+01  0.1011E+02
   20    5   24    1.48   30.05  0.2444E+01  0.6475E+01
   33   29   14    1.48   82.00  0.6669E+01  0.1014E+02
   37    7   18    1.48   61.66  0.5015E+01  0.1387E+02
   26   10   22    1.49   24.90  0.2025E+01  0.9310E+01
   46    4   12    1.49  104.30  0.8483E+01  0.9660E+01
   30   28   16    1.49   13.21  0.1074E+01  0.1019E+02
   41   31    4    1.49   66.15  0.5380E+01  0.1070E+02
   40   31    6    1.49   88.70  0.7214E+01  0.1043E+02
   42    5   15    1.49   60.31  0.4905E+01  0.9049E+01
   29    9   21    1.49   94.05  0.7649E+01  0.7326E+01
   41   31    3    1.49   53.32  0.4337E+01  0.1055E+02
   37   30   10    1.50   73.37  0.5967E+01  0.1065E+02
   41   31    2    1.50   73.35  0.5966E+01  0.1081E+02
   25   10   22    1.50   30.14  0.2451E+01  0.8765E+01
   50    2    6    1.50   37.22  0.3027E+01  0.9853E+01
   31   28   15    1.50  109.78  0.8929E+01  0.1484E+02
   41   31    0    1.50   20.60  0.7340E+01  0.1487E+02
   48    3    9    1.51   98.61  0.8020E+01  0.1552E+02
   44    4   13    1.51  113.82 -0.1033E+02  0.9363E+01
   51    2    2    1.51  112.47 -0.1020E+02  0.1031E+02
   47    3   10    1.51  119.51 -0.1084E+02  0.9956E+01
   39   30    7    1.51   68.30 -0.6197E+01  0.1025E+02
   37   29   10    1.51   47.80 -0.4337E+01  0.1606E+02
   51    2    1    1.51  167.27 -0.1518E+02  0.1893E+02
   24   10   22    1.51   16.71 -0.1517E+01  0.8148E+01
   50    2    5    1.51   81.78 -0.7420E+01  0.9692E+01
   49    2    7    1.52   28.75 -0.2609E+01  0.9893E+01
   25   23   19    1.52  109.17 -0.9905E+01  0.1118E+02
   23   17   21    1.52   82.44 -0.7480E+01  0.9987E+01
   20    8   23    1.52   13.78 -0.1251E+01  0.7793E+01
   24   20   20    1.52   24.94 -0.2263E+01  0.1296E+02
   39    5   16    1.53  117.77 -0.1069E+02  0.8633E+01
   42    4   14    1.53   20.26 -0.1838E+01  0.8952E+01
   30    8   20    1.53   66.92 -0.6072E+01  0.7471E+01
   33   28   13    1.53   84.62 -0.7678E+01  0.1522E+02
   25   25   18    1.53   20.54 -0.8122E+01  0.8757E+01
   43    4   13    1.53   36.07 -0.3273E+01  0.1576E+02
   32   27   14    1.53   30.26 -0.2745E+01  0.1320E+02
   37   29    9    1.54  113.15 -0.1027E+02  0.9688E+01
   44    3   12    1.54   91.16 -0.8272E+01  0.9426E+01
   21   12   22    1.54   61.78 -0.5606E+01  0.8312E+01
   45    3   11    1.54   84.01 -0.7622E+01  0.9425E+01
   40   30    2    1.54  264.77 -0.2402E+02  0.1133E+02
   34    6   18    1.55   91.28 -0.8282E+01  0.8275E+01
   38   29    7    1.55   53.39 -0.4845E+01  0.9979E+01
   36   28   10    1.55   58.56 -0.5313E+01  0.1469E+02
   21   10   22    1.55   62.23 -0.5646E+01  0.7928E+01
   39   29    5    1.56   89.85 -0.8152E+01  0.1000E+02
   26   25   17    1.56   39.17 -0.3554E+01  0.9050E+01
   24   24   18    1.56   22.28 -0.1018E+02  0.8920E+01
   39    4   15    1.57   13.35  0.3955E+01  0.1214E+02
   43    3   12    1.57   32.49  0.9628E+01  0.9314E+01
   39   29    4    1.57   43.81  0.1298E+02  0.1065E+02
   31   26   14    1.57   43.67  0.1294E+02  0.9621E+01
   33    6   18    1.57   14.99  0.4443E+01  0.7631E+01
   37    5   16    1.57   37.88  0.1123E+02  0.1390E+02
   47    1    6    1.59   24.08  0.7137E+01  0.1377E+02
   39   28    3    1.60   23.31  0.6908E+01  0.1012E+02
   33   26   12    1.60   17.95  0.5320E+01  0.1324E+02
   44    2   10    1.60   17.39  0.5153E+01  0.9448E+01
   43    2   11    1.60   10.47  0.3104E+01  0.1287E+02
   34    5   17    1.60   14.72  0.4362E+01  0.8135E+01
   20   13   21    1.60   21.13  0.6263E+01  0.8284E+01
   38   28    5    1.60   21.25  0.6298E+01  0.9889E+01
   27   24   16    1.60   48.46  0.1436E+02  0.1426E+02
   21   17   20    1.60   13.75  0.4076E+01  0.9146E+01
   22    9   21    1.60   35.68  0.1057E+02  0.8974E+01
   39   28    2    1.60   11.22  0.3324E+01  0.9830E+01
   47    1    5    1.61   14.38  0.4261E+01  0.9413E+01
   39   28    1    1.61   10.42  0.3089E+01  0.1004E+02
   37   27    7    1.61   40.38  0.1197E+02  0.1021E+02
   31    6   18    1.61   44.32  0.1313E+02  0.1189E+02
   37    4   15    1.62   11.75  0.3482E+01  0.8940E+01
   21    9   21    1.62   11.11  0.3293E+01  0.7725E+01
   24   23   17    1.63   21.07  0.6244E+01  0.8159E+01
   26   23   16    1.64   41.40  0.1485E+02  0.1516E+02
   33   25   11    1.65   42.10  0.1510E+02  0.1432E+02
   34    4   16    1.65   32.05  0.1149E+02  0.8691E+01
   40    2   12    1.66   11.38  0.4082E+01  0.8759E+01
   36   26    7    1.66   34.80  0.1248E+02  0.9720E+01
   27   23   15    1.66   11.16  0.4001E+01  0.8779E+01
   46    0    3    1.67   10.67  0.8059E+01  0.9541E+01
   21   20   18    1.67   26.56  0.9524E+01  0.1161E+02
   36   26    6    1.68   33.31  0.1194E+02  0.1519E+02
   31   24   12    1.68   19.29  0.6917E+01  0.9028E+01
   35   25    8    1.68   44.21  0.1586E+02  0.1425E+02
   39    2   12    1.69   17.55  0.6295E+01  0.8619E+01
   42    1    9    1.69   26.63  0.9548E+01  0.9159E+01
   36    3   14    1.69   10.55  0.3784E+01  0.9283E+01
   32   24   11    1.69   12.94  0.4639E+01  0.9583E+01
   37    2   13    1.71   13.97  0.9143E+01  0.9555E+01
   43    0    7    1.71   10.51  0.1433E+02  0.1464E+02
   45    0    2    1.71   12.46  0.1931E+02  0.1756E+02
   37   25    1    1.73   18.19  0.1190E+02  0.9884E+01
   32   23   10    1.75   12.63  0.8263E+01  0.1608E+02
   35   24    6    1.75   10.52  0.6886E+01  0.9764E+01
   44    1    0    1.76   11.61  0.1700E+02  0.1405E+02
   36   24    3    1.77   18.18  0.1190E+02  0.1062E+02
   32   23    9    1.78   32.10 -0.9124E+01  0.9053E+01
   43    1    2    1.79   47.03 -0.1337E+02  0.9179E+01
   37    1   11    1.80   18.79 -0.5341E+01  0.8352E+01
   43    1    0    1.80   15.04 -0.1213E+02  0.1300E+02
   27   21   13    1.80   50.91 -0.1447E+02  0.7964E+01
   31   22   10    1.80   17.53 -0.4983E+01  0.1352E+02
   42    1    4    1.81   10.67 -0.3032E+01  0.9219E+01
   34   23    6    1.81   14.50 -0.4122E+01  0.9303E+01
   30    3   15    1.82   33.93 -0.9643E+01  0.1074E+02
   35    1   12    1.82   16.72 -0.4752E+01  0.8560E+01
   35   23    3    1.83   46.80 -0.1330E+02  0.9497E+01
   35   23    2    1.84   12.57 -0.3573E+01  0.9911E+01
   42    1    1    1.84   26.18 -0.7441E+01  0.9771E+01
   31    2   14    1.85   22.26 -0.6326E+01  0.7361E+01
   41    1    4    1.85   36.70 -0.1043E+02  0.9323E+01
   22   19   15    1.85   28.88 -0.8207E+01  0.1319E+02
   34    1   12    1.86   35.45 -0.1007E+02  0.8031E+01
   17   17   17    1.86   11.57 -0.7340E+01  0.8688E+01
   33   22    6    1.87   48.47 -0.1378E+02  0.9430E+01
   15    7   19    1.87   14.47 -0.4114E+01  0.5581E+01
   27   20   12    1.87   33.81 -0.9609E+01  0.8472E+01
   33   22    5    1.89   16.62 -0.1147E+02  0.1057E+02
   39    1    6    1.89   16.00 -0.1104E+02  0.9010E+01
   33    1   12    1.89   19.82 -0.1368E+02  0.8169E+01
   31   21    8    1.90   10.80 -0.7452E+01  0.9144E+01
   29    2   14    1.92   16.00 -0.1104E+02  0.7699E+01
   30    1   13    1.95   21.91 -0.1512E+02  0.1204E+02
   29   20    9    1.95   12.92 -0.8913E+01  0.1397E+02
   25   18   12    1.98   18.52 -0.1278E+02  0.1248E+02
   38    2    4    1.99   16.98 -0.1172E+02  0.1136E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT       140
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1496
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1225
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      271


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.4 sec
