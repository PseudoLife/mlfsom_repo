
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.419    77.419    38.754  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      18      18
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      18
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000827  0.000700
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
 BEAM_DIVERGENCE_E.S.D.=   0.044
 REFLECTING_RANGE_E.S.D.=   0.107
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

    1417 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      18 ...      18
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1417 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               44
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

 *  44        aP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0 -1  0  0  0  0 -1  0  0
 *  34        mP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  33        mP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.4   38.8   77.4  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  25        mC          0.0     109.5  109.5   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  32        oP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  20        mC          0.0     109.5  109.5   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.5  109.5   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.4   77.4   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.6   38.8   77.4  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.6   38.8   77.4  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  159.6   77.4  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.6   77.4  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  159.6   77.4  90.0  90.0  75.9    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    29        mC        250.0      38.8  159.6   77.4  90.0  90.0  75.9    0  0  1  0  2  0  1  0  0  1  0  0
    27        mC        500.0     159.6   38.8  109.5  90.0 133.3  75.9   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.5  116.1  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  159.6  159.6  86.6 104.1 104.1    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.5  116.1   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        997.7     173.1   77.4   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        997.7      77.4  173.1   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        997.7      77.4  173.1   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        997.7      77.4   77.4   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.1  116.1  116.1  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.6   86.6  116.1 116.5  63.5 101.6    1  0  1  0  0 -1 -1  0  1  1 -1  0
     3        cP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.6   86.6  116.1 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.5   86.6   86.6  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.6   86.6  109.5  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.6  247.9 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.6   86.6   77.4  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.6   86.6   77.4  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.4  177.4  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.6   86.6  177.4 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.6   86.6   77.4  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.6   86.6   86.6  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.4  109.5   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.4  177.4 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.6  177.4   77.4 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    5.6      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.785   -23.440     0.020
 UNIT_CELL_B-AXIS=    23.440    73.785     0.014
 UNIT_CELL_C-AXIS=    -0.011    -0.003    38.754
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       267 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.28
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.336    77.336    38.705  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.1E-01 1.1E-01 5.6E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012931  0.012931  0.025837  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.698   -23.443     0.028
 COORDINATES OF UNIT CELL B-AXIS    23.443    73.698     0.010
 COORDINATES OF UNIT CELL C-AXIS    -0.015    -0.001    38.705
 CRYSTAL MOSAICITY (DEGREES)     0.107
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006343  0.001752
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.91   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1261.36   1242.11
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.04
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.34    77.34    38.70  90.000  90.000  90.000


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
  0.05  0.15        91        118.108       1.521
  0.15  0.25       171        211.530       0.849
  0.25  0.35       196        135.515       1.326
  0.35  0.45       145        111.116       1.617
  0.45  0.55       173        254.905       0.705
  0.55  0.65       191        160.384       1.120
  0.65  0.75       185        186.975       0.961
  0.75  0.85       175        154.939       1.159
  0.85  0.95        90        317.084       0.567
 
 -0.05  0.95      1417        179.651       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               236.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.973
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    17.1 XMAX=    17.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59627 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        502


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               236.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.971
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     8.5 XMAX=  2455.9 NXBIN=    3
 YMIN=     6.5 YMAX=  2516.8 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        502


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               236.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.968
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    17.1 XMAX=    17.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1262
 DETECTOR_SURFACE_POSITION=    1648    1688
 DETECTOR_SURFACE_POSITION=     817    1688
 DETECTOR_SURFACE_POSITION=     817     835
 DETECTOR_SURFACE_POSITION=    1648     835
 DETECTOR_SURFACE_POSITION=    2173    1661
 DETECTOR_SURFACE_POSITION=    1622    2227
 DETECTOR_SURFACE_POSITION=     843    2227
 DETECTOR_SURFACE_POSITION=     292    1661
 DETECTOR_SURFACE_POSITION=     292     862
 DETECTOR_SURFACE_POSITION=     843     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2173     862
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
 1.041E+00  2.117E-03   21.31


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

   24.456   3.150    17.99   1.55      2.90      2.79      76     104      22
    3.150   2.237    10.11   1.21      4.72      4.38     188     224       8
    2.237   1.829     2.16   0.77     17.10     19.43     214     249       0
    1.829   1.585     0.02   0.79   3270.67   3536.35     244     290       0
    1.585   1.418     0.00   0.70   1115.83   1339.91     280     343       0
    1.418   1.295     0.02   0.91   1259.91   1260.34     126     176       0
   --------------------------------------------------------------------------
   24.456   1.295     3.38   0.90      6.60      6.82    1128    1386      30
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      18      18)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681   0.3803E+01  0.1247E+02         0.31           1 
   47    0    0     1.645   0.1051E+02  0.1264E+02         0.83           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     5.6% OF MEAN INTENSITY


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

     3.84          64      41      1285        3.2%       2.6%      2.7%       46   22.78      3.6%    99.8*     0    0.000       0
     2.73         100      60      2098        2.9%       3.6%      3.3%       80   18.73      5.0%    99.3*     0    0.000       0
     2.24         164      95      2665        3.6%       6.0%      4.9%      138   12.00      8.5%    99.2*     0    0.000       0
     1.94         161      91      3101        2.9%      11.4%     13.4%      140    4.34     16.2%    97.8*     0    0.000       0
     1.74         172     100      3481        2.9%     -99.9%    -99.9%      144    0.00    -99.9%     7.0      0    0.000       0
     1.58         206     119      3848        3.1%    1414.2%   1549.8%      174    0.06   1999.9%   -10.7      0    0.000       0
     1.47         223     133      4143        3.2%    1937.4%   2101.5%      180   -0.02   2739.9%   -17.1      0    0.000       0
     1.37         211     122      4429        2.8%    1438.8%   1858.7%      178   -0.02   2034.8%     1.4      0    0.000       0
     1.30          85      61      4740        1.3%     451.3%    387.4%       48    0.08    638.2%    42.0      0    0.000       0
    total        1386     822     29790        2.8%       6.6%      6.8%     1128    4.38      9.3%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1417
 NUMBER OF REJECTED MISFITS                                30
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1386
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    822



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

   WILSON LINE (using all data) : A=   7.656 B=  32.387 CORRELATION=  0.90
      #      RES      SS        <I>       log(<I>)       BO
      17     7.253   0.005  2.0421E+03       7.622       3.6
      20     4.585   0.012  2.4797E+03       7.816      -6.7
      27     3.511   0.020  1.0799E+03       6.985      16.6
      29     2.991   0.028  3.9238E+02       5.972      30.1
      37     2.643   0.036  3.7994E+02       5.940      24.0
      52     2.397   0.044  1.9524E+02       5.274      27.4
      39     2.205   0.051  1.1806E+02       4.771      28.1
      45     2.056   0.059  4.8433E+01       3.880      31.9
      45     1.930   0.067 -1.2094E+00    -999.999   -1000.0
      53     1.828   0.075  4.6709E-01      -0.761      56.3
      50     1.739   0.083 -6.5271E-01    -999.999   -1000.0
      52     1.665   0.090  4.5368E-01      -0.790      46.8
      52     1.596   0.098  2.2408E+00       0.807      34.9
      68     1.536   0.106 -1.2642E+00    -999.999   -1000.0
      74     1.480   0.114  1.2575E+00       0.229      32.5
      57     1.432   0.122 -1.2151E+00    -999.999   -1000.0
      45     1.390   0.129  4.5574E-01      -0.786      32.6
      42     1.351   0.137  5.1821E-01      -0.657      30.3
      18     1.314   0.145  2.6550E+00       0.976      23.1


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.253        1.836        1.230        0.658
       3     4.585        1.276        0.785        0.359
       2     3.511        0.425        0.096        0.015
       4     2.991        1.182        0.541        0.188
       3     2.643        0.297        0.069        0.012
       5     2.397        1.032        0.604        0.274
       2     2.205        1.365        0.566        0.170
       8     2.056        0.297        0.114        0.034
      10     1.828      136.879      579.180     4276.001
       8     1.739       41.760       49.603      345.871
       3     1.536       19.893       41.101       79.115
       6     1.480       23.661       52.307      111.337
       3     1.432       17.810       19.606       47.427
      84   overall       37.519       -4.829     1026.480


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.253        0.454        0.317        0.178
      17     4.585        1.093        1.232        1.162
      25     3.511        0.810        0.627        0.459
      25     2.991        0.812        0.642        0.430
      34     2.643        1.133        1.250        1.284
      47     2.397        0.962        0.850        0.636
      37     2.205        0.781        0.586        0.388
      37     2.056        1.492        2.035        2.504
      43     1.828      127.691      -13.270     6353.887
      41     1.665      132.581      352.000    12288.946
      49     1.596        8.074       12.290       41.905
      68     1.480       25.392       85.949      859.231
      45     1.390      176.375       47.404    10834.133
      41     1.351      533.128    10408.519   388624.324
      18     1.314       14.128       23.201       92.145
     738   overall       58.714      609.545    23388.295

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.253  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.32  0.30  0.29
       3     4.585  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.511  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       4     2.991  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.643  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       5     2.397  0.00  0.00  0.00  0.42  0.38  1.07  1.00  0.95  0.91  0.88
       2     2.205  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.056  1.51  1.08  1.20  1.32  1.20  1.34  1.26  1.19  1.14  1.10
      10     1.828  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.73
       8     1.739  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       3     1.536  1.34  0.96  0.80  0.70  0.64  0.59  1.12  1.06  1.01  0.98
       6     1.480  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.98
       3     1.432  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      84   overall  0.91  0.65  0.57  0.55  0.50  0.55  0.54  0.53  0.51  0.54


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.253  2.63  2.76  1.93  1.52  1.91  1.66  1.66  1.51  1.40  1.32
      17     4.585  1.24  1.62  1.13  1.07  1.20  1.30  1.17  1.18  1.29  1.21
      25     3.511  0.84  0.88  0.62  0.49  0.61  0.71  0.95  0.94  0.88  1.08
      25     2.991  1.68  1.99  1.39  1.09  1.12  1.15  1.11  1.23  1.15  1.14
      34     2.643  0.62  1.30  1.36  1.16  0.97  0.98  1.05  1.07  1.04  1.07
      47     2.397  1.12  0.94  1.07  1.16  1.24  1.18  1.06  1.08  1.11  1.08
      37     2.205  0.85  0.75  0.83  0.74  0.89  0.90  0.97  1.03  1.00  0.98
      37     2.056  0.85  1.04  0.73  0.98  1.17  1.20  1.07  1.03  1.05  1.03
      43     1.828  5.38  2.82  1.97  1.55  1.30  1.13  1.02  0.93  0.86  0.81
      41     1.665  5.13  2.69  2.07  1.63  1.36  1.19  1.07  0.97  0.90  0.85
      49     1.596  4.29  2.25  1.57  1.30  1.14  1.04  0.93  0.85  0.79  0.77
      68     1.480  5.10  2.76  1.93  1.52  1.27  1.11  1.02  0.93  0.92  0.88
      45     1.390  5.60  2.94  2.06  1.62  1.36  1.18  1.06  0.97  0.90  0.84
      41     1.351  4.87  2.56  1.98  1.55  1.30  1.14  1.02  0.97  0.90  0.85
      18     1.314  5.84  3.06  2.14  1.69  1.41  1.23  1.10  1.11  1.03  0.97
     738   overall  2.45  1.50  1.13  0.95  0.88  0.81  0.76  0.74  0.71  0.69


 List of     87 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   54    5    9    1.35  123.85  0.6418E+02  0.1049E+02 "alien"
   54    5    8    1.37  104.18  0.5399E+02  0.1893E+02 "alien"
   36   26    6    1.68   46.16  0.2094E+02  0.1552E+02 "alien"
   46    6   15    1.40   36.42  0.1660E+02  0.9989E+01 "alien"
   31   28   15    1.50   32.65  0.4105E+02  0.1740E+02 "alien"
   33   23    7    1.82   32.39  0.1513E+02  0.1567E+02 "alien"
   49    7   15    1.34   32.36  0.1677E+02  0.1608E+02 "alien"
   48    6   14    1.38   30.31  0.1382E+02  0.1394E+02 "alien"
   38    2   12    1.72   30.17 -0.1969E+02  0.1436E+02 "alien"
   28   23   14    1.69   29.74  0.1349E+02  0.8032E+01 "alien"
   37   34   14    1.34   29.48  0.1528E+02  0.1046E+02 "alien"
   36   32   13    1.41   28.31 -0.3440E+02  0.1820E+02 "alien"
   31   30   18    1.38   28.27  0.1288E+02  0.9969E+01 "alien"
   19   14   20    1.67   27.83  0.1263E+02  0.5746E+01 "alien"
   20    9   21    1.63   27.70  0.1257E+02  0.8229E+01 "alien"
   33   11   22    1.38   27.62  0.1259E+02  0.1147E+02 "alien"
   53    3    3    1.45   27.42 -0.3332E+02  0.1783E+02 "alien"
   43    1    3    1.78   27.35 -0.1785E+02  0.8793E+01 "alien"
   37    1   11    1.80   27.28  0.1274E+02  0.8421E+01 "alien"
   28   13   24    1.36   25.87  0.1341E+02  0.1216E+02 "alien"
   28   23   21    1.39   25.84  0.1178E+02  0.8431E+01 "alien"
   51    4   10    1.41   25.31  0.1153E+02  0.1560E+02 "alien"
   41    1    3    1.87   25.12  0.1173E+02  0.1161E+02 "alien"
   20   18   16    1.85   23.78  0.1111E+02  0.9704E+01 "alien"
   43    1    2    1.79   21.98  0.1027E+02  0.8709E+01 "alien"
   29   26   20    1.39   21.96  0.1001E+02  0.1589E+02 "alien"
   24   19   14    1.87   21.77  0.1017E+02  0.1202E+02 "alien"
   50    5   12    1.39   21.60  0.9846E+01  0.1032E+02 "alien"
   27   20   12    1.87   21.12  0.9867E+01  0.8011E+01 "alien"
   36    2   13    1.74   21.06 -0.1374E+02  0.7821E+01 "alien"
   52    6   12    1.34   20.75  0.1075E+02  0.1037E+02 "alien"
   37    9   20    1.40   20.14  0.9180E+01  0.8513E+01 "alien"
   36   32   14    1.39   19.15  0.8729E+01  0.1827E+02 "alien"
   37   26    3    1.70   19.10  0.8665E+01  0.9269E+01 "alien"
   38    3   14    1.64   19.08  0.8658E+01  0.8431E+01 "alien"
   28   21   22    1.38   18.89  0.8611E+01  0.8083E+01 "alien"
   37   26    4    1.68   18.58  0.8431E+01  0.9262E+01 "alien"
   35   11   22    1.35   17.95  0.9300E+01  0.1207E+02 "alien"
   46    7   16    1.37   17.72  0.8074E+01  0.9565E+01 "alien"
   31    4   16    1.73   17.45 -0.1139E+02  0.1305E+02 "alien"
   34    4   16    1.65   17.14  0.7776E+01  0.8072E+01 "alien"
   30   23   12    1.73   16.97 -0.1108E+02  0.8117E+01 "alien"
   39   34   12    1.36   16.68  0.8646E+01  0.1665E+02 "alien"
   31   30   17    1.41   16.58  0.7558E+01  0.1009E+02 "alien"
   36   25    4    1.74   16.36 -0.1068E+02  0.8940E+01 "alien"
   47    2    9    1.54   15.91 -0.2011E+02  0.1620E+02 "alien"
   20   19   18    1.71   15.90 -0.1038E+02  0.6876E+01 "alien"
   18    9   21    1.66   15.79  0.7162E+01  0.5158E+01 "alien"
   27   23   15    1.67   15.61  0.7084E+01  0.7913E+01 "alien"
   34   11   22    1.37   15.61  0.8088E+01  0.8070E+01 "alien"
   39    0    9    1.80   14.94  0.1968E+02  0.9140E+01 "alien"
   35    1   12    1.82   14.88  0.6949E+01  0.8521E+01 "alien"
   23   17   21    1.52   14.66 -0.1854E+02  0.1069E+02 "alien"
   52    2    3    1.48   14.59  0.1834E+02  0.9821E+01 "alien"
   34   23    5    1.83   14.23  0.6647E+01  0.9223E+01 "alien"
   51    4    9    1.43   14.02 -0.1704E+02  0.9738E+01 "alien"
   49    4   11    1.44   14.02 -0.1704E+02  0.9872E+01 "alien"
   30   22   11    1.79   13.74  0.6416E+01  0.8187E+01 "alien"
   42    2   11    1.63   13.57  0.6157E+01  0.8665E+01 "alien"
   32   23   10    1.75   13.47 -0.8795E+01  0.1540E+02 "alien"
   39   33   10    1.41   12.96  0.5908E+01  0.1158E+02 "alien"
   43    1    0    1.80   12.54  0.1394E+02  0.1298E+02 "alien"
   54    4    7    1.38   12.38  0.5643E+01  0.1020E+02 "alien"
   28    2   14    1.95   12.33 -0.1491E+02  0.1178E+02 "alien"
   33   25   10    1.68   12.19  0.5532E+01  0.1386E+02 "alien"
   46    1    6    1.63   12.16  0.2726E+02  0.1389E+02 "alien"
   42   32    3    1.46   12.10 -0.1470E+02  0.1731E+02 "alien"
   37   31   12    1.43   12.02 -0.1461E+02  0.1598E+02 "alien"
   45    0    4    1.69   12.00  0.1250E+02  0.9149E+01 "alien"
   16    8   20    1.77   11.91 -0.7771E+01  0.5892E+01 "alien"
   45    0    3    1.70   11.67 -0.1711E+02  0.9079E+01 "alien"
   34    1   12    1.86   11.67  0.5450E+01  0.8043E+01 "alien"
   31   27   15    1.52   11.65 -0.1473E+02  0.1482E+02 "alien"
   32   28   20    1.33   11.59  0.3078E+02  0.1495E+02 "alien"
   34    7   19    1.50   11.54  0.1451E+02  0.1153E+02 "alien"
   44    1    1    1.76   11.49 -0.7502E+01  0.9420E+01 "alien"
   30    1   13    1.95   11.36 -0.1374E+02  0.1202E+02 "alien"
   53    5    9    1.38   10.90  0.4968E+01  0.1182E+02 "alien"
   42    1    2    1.83   10.78  0.5033E+01  0.9214E+01 "alien"
   37   27    6    1.63   10.74  0.4871E+01  0.1295E+02 "alien"
   41   31    3    1.49   10.58  0.1330E+02  0.9929E+01 "alien"
   36   24    3    1.77   10.58 -0.6903E+01  0.9679E+01 "alien"
   38    1    7    1.91   10.45 -0.1263E+02  0.8771E+01 "alien"
   46    4   12    1.49   10.38  0.1305E+02  0.9543E+01 "alien"
   22   14   22    1.51   10.36 -0.1310E+02  0.8232E+01 "alien"
   44    3   12    1.54   10.21 -0.1291E+02  0.9090E+01 "alien"
   25   13   24    1.39   10.00  0.4558E+01  0.6345E+01 "alien"


 List of     87 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   32   28   20    1.33   11.59  0.3078E+02  0.1495E+02
   49    7   15    1.34   32.36  0.1677E+02  0.1608E+02
   52    6   12    1.34   20.75  0.1075E+02  0.1037E+02
   37   34   14    1.34   29.48  0.1528E+02  0.1046E+02
   35   11   22    1.35   17.95  0.9300E+01  0.1207E+02
   54    5    9    1.35  123.85  0.6418E+02  0.1049E+02
   28   13   24    1.36   25.87  0.1341E+02  0.1216E+02
   39   34   12    1.36   16.68  0.8646E+01  0.1665E+02
   34   11   22    1.37   15.61  0.8088E+01  0.8070E+01
   54    5    8    1.37  104.18  0.5399E+02  0.1893E+02
   46    7   16    1.37   17.72  0.8074E+01  0.9565E+01
   53    5    9    1.38   10.90  0.4968E+01  0.1182E+02
   28   21   22    1.38   18.89  0.8611E+01  0.8083E+01
   31   30   18    1.38   28.27  0.1288E+02  0.9969E+01
   33   11   22    1.38   27.62  0.1259E+02  0.1147E+02
   54    4    7    1.38   12.38  0.5643E+01  0.1020E+02
   48    6   14    1.38   30.31  0.1382E+02  0.1394E+02
   29   26   20    1.39   21.96  0.1001E+02  0.1589E+02
   36   32   14    1.39   19.15  0.8729E+01  0.1827E+02
   50    5   12    1.39   21.60  0.9846E+01  0.1032E+02
   25   13   24    1.39   10.00  0.4558E+01  0.6345E+01
   28   23   21    1.39   25.84  0.1178E+02  0.8431E+01
   46    6   15    1.40   36.42  0.1660E+02  0.9989E+01
   37    9   20    1.40   20.14  0.9180E+01  0.8513E+01
   51    4   10    1.41   25.31  0.1153E+02  0.1560E+02
   31   30   17    1.41   16.58  0.7558E+01  0.1009E+02
   39   33   10    1.41   12.96  0.5908E+01  0.1158E+02
   36   32   13    1.41   28.31 -0.3440E+02  0.1820E+02
   51    4    9    1.43   14.02 -0.1704E+02  0.9738E+01
   37   31   12    1.43   12.02 -0.1461E+02  0.1598E+02
   49    4   11    1.44   14.02 -0.1704E+02  0.9872E+01
   53    3    3    1.45   27.42 -0.3332E+02  0.1783E+02
   42   32    3    1.46   12.10 -0.1470E+02  0.1731E+02
   52    2    3    1.48   14.59  0.1834E+02  0.9821E+01
   46    4   12    1.49   10.38  0.1305E+02  0.9543E+01
   41   31    3    1.49   10.58  0.1330E+02  0.9929E+01
   34    7   19    1.50   11.54  0.1451E+02  0.1153E+02
   31   28   15    1.50   32.65  0.4105E+02  0.1740E+02
   22   14   22    1.51   10.36 -0.1310E+02  0.8232E+01
   31   27   15    1.52   11.65 -0.1473E+02  0.1482E+02
   23   17   21    1.52   14.66 -0.1854E+02  0.1069E+02
   47    2    9    1.54   15.91 -0.2011E+02  0.1620E+02
   44    3   12    1.54   10.21 -0.1291E+02  0.9090E+01
   46    1    6    1.63   12.16  0.2726E+02  0.1389E+02
   42    2   11    1.63   13.57  0.6157E+01  0.8665E+01
   37   27    6    1.63   10.74  0.4871E+01  0.1295E+02
   20    9   21    1.63   27.70  0.1257E+02  0.8229E+01
   38    3   14    1.64   19.08  0.8658E+01  0.8431E+01
   34    4   16    1.65   17.14  0.7776E+01  0.8072E+01
   18    9   21    1.66   15.79  0.7162E+01  0.5158E+01
   27   23   15    1.67   15.61  0.7084E+01  0.7913E+01
   19   14   20    1.67   27.83  0.1263E+02  0.5746E+01
   36   26    6    1.68   46.16  0.2094E+02  0.1552E+02
   33   25   10    1.68   12.19  0.5532E+01  0.1386E+02
   37   26    4    1.68   18.58  0.8431E+01  0.9262E+01
   28   23   14    1.69   29.74  0.1349E+02  0.8032E+01
   45    0    4    1.69   12.00  0.1250E+02  0.9149E+01
   37   26    3    1.70   19.10  0.8665E+01  0.9269E+01
   45    0    3    1.70   11.67 -0.1711E+02  0.9079E+01
   20   19   18    1.71   15.90 -0.1038E+02  0.6876E+01
   38    2   12    1.72   30.17 -0.1969E+02  0.1436E+02
   30   23   12    1.73   16.97 -0.1108E+02  0.8117E+01
   31    4   16    1.73   17.45 -0.1139E+02  0.1305E+02
   36   25    4    1.74   16.36 -0.1068E+02  0.8940E+01
   36    2   13    1.74   21.06 -0.1374E+02  0.7821E+01
   32   23   10    1.75   13.47 -0.8795E+01  0.1540E+02
   44    1    1    1.76   11.49 -0.7502E+01  0.9420E+01
   16    8   20    1.77   11.91 -0.7771E+01  0.5892E+01
   36   24    3    1.77   10.58 -0.6903E+01  0.9679E+01
   43    1    3    1.78   27.35 -0.1785E+02  0.8793E+01
   30   22   11    1.79   13.74  0.6416E+01  0.8187E+01
   43    1    2    1.79   21.98  0.1027E+02  0.8709E+01
   37    1   11    1.80   27.28  0.1274E+02  0.8421E+01
   43    1    0    1.80   12.54  0.1394E+02  0.1298E+02
   39    0    9    1.80   14.94  0.1968E+02  0.9140E+01
   33   23    7    1.82   32.39  0.1513E+02  0.1567E+02
   35    1   12    1.82   14.88  0.6949E+01  0.8521E+01
   34   23    5    1.83   14.23  0.6647E+01  0.9223E+01
   42    1    2    1.83   10.78  0.5033E+01  0.9214E+01
   20   18   16    1.85   23.78  0.1111E+02  0.9704E+01
   34    1   12    1.86   11.67  0.5450E+01  0.8043E+01
   24   19   14    1.87   21.77  0.1017E+02  0.1202E+02
   41    1    3    1.87   25.12  0.1173E+02  0.1161E+02
   27   20   12    1.87   21.12  0.9867E+01  0.8011E+01
   38    1    7    1.91   10.45 -0.1263E+02  0.8771E+01
   30    1   13    1.95   11.36 -0.1374E+02  0.1202E+02
   28    2   14    1.95   12.33 -0.1491E+02  0.1178E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        32
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1417
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1339
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       78


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
