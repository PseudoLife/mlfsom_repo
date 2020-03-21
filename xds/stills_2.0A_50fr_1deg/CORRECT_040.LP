
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.604    77.604    38.843  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      40      40
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      40
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000999  0.000648
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
 REFLECTING_RANGE_E.S.D.=   0.183
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

    1166 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      40 ...      40
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1166 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               39
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

 *  44        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.8  109.8   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  31        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  23        oC          0.0     109.8  109.8   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  35        mP          0.0      77.6   38.8   77.6  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  20        mC          0.0     109.8  109.8   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  33        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  21        tP          0.0      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
 *  32        oP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    37        mC        250.0     160.0   38.8   77.6  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     160.0   38.8   77.6  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.8  160.0   77.6  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.8  160.0   77.6  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.8  160.0   77.6  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  160.0   77.6  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    19        oI        500.0      38.8  109.8  116.4  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    27        mC        500.0     160.0   38.8  109.8  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    18        tI        625.0     109.8  116.4   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    26        oF        625.0      38.8  160.0  160.0  86.6 104.0 104.1    0  0 -1  0 -2  0  1  0  0  2  1  0
    41        mC        998.0     173.5   77.6   38.8  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    22        hP        998.0      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.0      77.6  173.5   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.0      77.6  173.5   38.8  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
     1        cF        999.0     116.4  116.4  116.4  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.8   86.8  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.8   86.8   86.8  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.8   86.8  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.8   86.8  248.5 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.8   86.8   77.6  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.8   77.6  177.8  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.8   86.8  177.8 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.8   86.8   86.8  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.8  109.8   38.8  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.8   77.6  177.8 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.8  177.8   77.6 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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
 UNIT_CELL_A-AXIS=    73.966   -23.487     0.064
 UNIT_CELL_B-AXIS=    23.487    73.966     0.021
 UNIT_CELL_C-AXIS=    -0.034    -0.000    38.842
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       123 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.25
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.569    77.569    38.813  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.4E-01 1.4E-01 5.5E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012892  0.012892  0.025765  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.927   -23.489     0.013
 COORDINATES OF UNIT CELL B-AXIS    23.489    73.927    -0.002
 COORDINATES OF UNIT CELL C-AXIS    -0.006     0.003    38.813
 CRYSTAL MOSAICITY (DEGREES)     0.183
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006578  0.001735
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.89   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1261.37   1242.45
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.23
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.57    77.57    38.81  90.000  90.000  90.000


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
  0.05  0.15        22         24.297       4.207
  0.15  0.25       146         94.490       1.082
  0.25  0.35       157        111.186       0.919
  0.35  0.45       172         49.543       2.063
  0.45  0.55       159        139.954       0.730
  0.55  0.65       174        100.034       1.022
  0.65  0.75       191        127.674       0.801
  0.75  0.85       131        102.722       0.995
  0.85  0.95        14         98.400       1.039
 
 -0.05  0.95      1166        102.218       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        7
 DEGREES OF FREEDOM OF CHI^2 FIT               181.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.875
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    39.1 XMAX=    39.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59281 NYBIN=    7
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        376


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               181.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.872
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    53.6 XMAX=  2451.5 NXBIN=    3
 YMIN=     5.9 YMAX=  2518.2 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        376


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               181.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.867
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    39.1 XMAX=    39.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1253    1262
 DETECTOR_SURFACE_POSITION=    1660    1689
 DETECTOR_SURFACE_POSITION=     845    1689
 DETECTOR_SURFACE_POSITION=     845     835
 DETECTOR_SURFACE_POSITION=    1660     835
 DETECTOR_SURFACE_POSITION=    2174    1662
 DETECTOR_SURFACE_POSITION=    1634    2228
 DETECTOR_SURFACE_POSITION=     871    2228
 DETECTOR_SURFACE_POSITION=     331    1662
 DETECTOR_SURFACE_POSITION=     331     862
 DETECTOR_SURFACE_POSITION=     871     296
 DETECTOR_SURFACE_POSITION=    1634     296
 DETECTOR_SURFACE_POSITION=    2174     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        376


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
 7.300E-01  1.562E-03   29.62


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

   24.529   2.888    17.17   1.39      2.42      2.28      88     116      12
    2.888   2.049     3.32   1.02     15.33     14.44     198     231       2
    2.049   1.675    -0.03   0.90    -99.90    -99.90     230     271       0
    1.675   1.452     0.07   0.82    731.79    798.40     250     310       2
    1.452   1.299     0.13   0.93    368.46    384.77     166     215       6
   --------------------------------------------------------------------------
   24.529   1.299     2.45   0.96      9.09      9.18     932    1143      22
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      40      40)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.686   0.3653E+01  0.1235E+02         0.30           1 
   47    0    0     1.650   0.5518E+01  0.1167E+02         0.47           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     5.1% OF MEAN INTENSITY


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

     3.85          59      37      1285        2.9%       2.1%      2.0%       44   26.90      3.0%    99.9*     0    0.000       0
     2.74          80      48      2098        2.3%       5.1%      4.2%       64   15.05      7.3%    99.3*     0    0.000       0
     2.24         138      78      2665        2.9%      14.3%     13.8%      120    4.84     20.2%    96.4*     0    0.000       0
     1.94         128      74      3099        2.4%      92.0%     99.8%      108    0.73    130.1%    48.9*     0    0.000       0
     1.74         140      82      3482        2.4%     -99.9%    -99.9%      116    0.00    -99.9%    12.5      0    0.000       0
     1.59         176     102      3850        2.6%     -99.9%    -99.9%      148    0.00    -99.9%    -8.4      0    0.000       0
     1.47         171     103      4144        2.5%     552.6%    671.4%      136    0.16    781.5%    26.0      0    0.000       0
     1.38         173     101      4426        2.3%     386.6%    440.2%      144    0.15    546.7%    11.4      0    0.000       0
     1.30          78      52      4742        1.1%     178.4%    164.1%       52    0.29    252.2%   -21.4      0    0.000       0
    total        1143     677     29791        2.3%       9.1%      9.2%      932    3.21     12.9%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1166
 NUMBER OF REJECTED MISFITS                                22
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1143
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    677



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

   WILSON LINE (using all data) : A=   6.446 B=  23.258 CORRELATION=  0.90
      #      RES      SS        <I>       log(<I>)       BO
      15     7.253   0.005  1.7983E+03       7.495    -110.3
      18     4.751   0.011  1.4992E+03       7.313     -39.1
      20     3.568   0.020  5.4812E+02       6.306       3.6
      26     2.987   0.028  1.3183E+02       4.882      27.9
      30     2.639   0.036  1.0355E+02       4.640      25.2
      42     2.407   0.043  3.4899E+01       3.552      33.5
      33     2.211   0.051  1.4513E+01       2.675      36.9
      36     2.061   0.059  7.1613E+00       1.969      38.1
      37     1.938   0.067 -1.1292E+00    -999.999   -1000.0
      44     1.831   0.075 -8.8103E-01    -999.999   -1000.0
      41     1.742   0.082 -1.7673E+00    -999.999   -1000.0
      42     1.670   0.090 -7.8179E-01    -999.999   -1000.0
      49     1.601   0.098 -1.2238E-01    -999.999   -1000.0
      48     1.541   0.105  3.0072E+00       1.101      25.4
      59     1.485   0.113  8.1225E-01      -0.208      29.4
      45     1.436   0.121  3.1466E+00       1.146      21.9
      41     1.395   0.129 -9.5769E-02    -999.999   -1000.0
      34     1.354   0.136  3.2343E+00       1.174      19.3
      17     1.317   0.144  3.1710E+00       1.154      18.3


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.253        2.097        1.621        0.979
       3     4.751        1.548        1.053        0.532
       2     3.568        0.379        0.081        0.013
       3     2.987        2.167        1.176        0.473
       3     2.407        0.082        0.010        0.001
       1     2.211        1.278        0.501        0.140
       6     2.061        1.269        0.617        0.286
      11     1.670       27.402       84.648      212.616
       3     1.541        2.677        1.859        0.922
      72   overall        1.989       -0.299        3.538


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.253        0.572        0.444        0.272
      15     4.751        1.329        1.774        1.879
      18     3.568        1.022        1.110        1.110
      23     2.987        0.793        0.692        0.504
      29     2.639        1.303        1.781        2.293
      39     2.407        1.111        1.077        0.951
      32     2.211        1.144        1.047        0.841
      30     2.061        1.664        1.859        2.118
      31     1.938       35.780      -11.593      500.122
      46     1.601       42.066        9.083      836.507
      45     1.541        5.831        7.339       18.632
      54     1.485       57.592      -35.785     2593.413
      42     1.436        4.327        3.945        7.245
      40     1.395       62.683       76.806     1798.539
      33     1.354        5.275        6.084       15.750
      17     1.317        9.613       21.637       61.923
     605   overall        6.813       -1.021      236.391

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.253  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.73
       3     4.751  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     3.568  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.73
       3     2.987  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.407  0.00  0.96  0.80  0.70  1.28  1.19  1.12  1.59  1.52  1.47
       1     2.211  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.061  1.34  0.96  0.80  0.70  0.64  0.59  0.84  0.80  0.76  0.73
      11     1.670  3.66  2.63  2.18  1.92  1.75  1.62  1.67  1.59  1.52  1.47
       3     1.541  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      72   overall  0.78  0.60  0.50  0.44  0.43  0.40  0.42  0.42  0.42  0.45


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.253  3.82  3.01  2.10  1.93  1.85  1.61  1.63  1.49  1.38  1.29
      15     4.751  3.50  1.84  1.80  1.42  1.19  1.18  1.46  1.33  1.46  1.37
      18     3.568  1.17  0.61  0.43  0.67  0.71  0.86  1.21  1.21  1.22  1.14
      23     2.987  1.83  1.68  1.34  1.45  1.55  1.35  1.47  1.34  1.25  1.17
      29     2.639  0.72  0.95  1.20  1.05  1.05  1.07  1.03  1.00  1.05  1.09
      39     2.407  1.62  0.99  0.89  0.86  0.85  0.85  0.81  0.84  0.82  0.93
      32     2.211  2.30  1.55  1.09  0.95  1.03  0.97  1.12  1.02  0.95  0.89
      30     2.061  2.80  1.47  1.29  1.11  0.93  0.81  0.93  0.85  0.84  0.95
      31     1.938  6.10  3.20  2.24  1.76  1.48  1.29  1.15  1.05  1.03  0.97
      46     1.601  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      45     1.541  4.44  2.33  1.80  1.42  1.19  1.03  0.93  0.89  0.82  0.77
      54     1.485  4.48  2.35  1.64  1.29  1.08  0.99  0.88  0.84  0.78  0.79
      42     1.436  3.75  2.10  1.56  1.23  1.09  1.00  0.90  0.82  0.76  0.72
      40     1.395  5.52  2.90  2.03  1.59  1.33  1.16  1.04  0.95  0.88  0.83
      33     1.354  2.87  1.67  1.29  1.10  0.92  0.87  0.78  0.77  0.77  0.81
      17     1.317  5.56  2.92  2.04  1.61  1.35  1.17  1.05  0.96  0.89  0.84
     605   overall  3.04  1.71  1.29  1.08  0.96  0.88  0.85  0.80  0.77  0.75


 List of     78 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   55    3    1    1.41  223.89 -0.2144E+02  0.9549E+01 "alien"
   43    7   17    1.40  208.00 -0.1992E+02  0.9055E+01 "alien"
   27   18   23    1.38  184.46 -0.1767E+02  0.1108E+02 "alien"
   35   27   10    1.60  178.97 -0.2190E+02  0.8296E+01 "alien"
   48    6   14    1.39  161.28 -0.1545E+02  0.1257E+02 "alien"
   40   33    9    1.41  137.97 -0.1321E+02  0.8941E+01 "alien"
   45    2    9    1.60  132.46 -0.1621E+02  0.1494E+02 "alien"
   36   33   14    1.38  117.61 -0.1126E+02  0.1149E+02 "alien"
   36   27    8    1.62  114.24 -0.1398E+02  0.9762E+01 "alien"
   52    4    9    1.41  113.23 -0.1084E+02  0.9436E+01 "alien"
   39   29    4    1.57  110.92 -0.1358E+02  0.1554E+02 "alien"
   32   26   13    1.59   94.83 -0.1161E+02  0.8196E+01 "alien"
   39   28    2    1.61   94.55 -0.1157E+02  0.8210E+01 "alien"
   37   33   13    1.39   86.06 -0.8241E+01  0.9897E+01 "alien"
   37    9   20    1.41   86.03 -0.8239E+01  0.8586E+01 "alien"
   30   12   23    1.38   83.68 -0.8014E+01  0.1166E+02 "alien"
   47    1    6    1.60   80.26 -0.9823E+01  0.1182E+02 "alien"
   38   28    5    1.61   80.02 -0.9793E+01  0.1495E+02 "alien"
   28   21   22    1.38   71.05 -0.6804E+01  0.8606E+01 "alien"
   43    2   11    1.60   67.26 -0.8232E+01  0.9476E+01 "alien"
   31   30   17    1.41   65.02 -0.6227E+01  0.9380E+01 "alien"
   36   32   14    1.39   61.95 -0.5933E+01  0.1598E+02 "alien"
   38   33   12    1.39   61.49 -0.5889E+01  0.9086E+01 "alien"
   30   25   14    1.61   53.21 -0.6512E+01  0.7638E+01 "alien"
   33    6   18    1.58   50.38 -0.6166E+01  0.7334E+01 "alien"
   33   25   11    1.65   44.34 -0.3466E+02  0.1322E+02 "alien"
   38   28    6    1.59   43.40 -0.5311E+01  0.1188E+02 "alien"
   45    2   10    1.57   40.17 -0.4917E+01  0.9679E+01 "alien"
   25   24   17    1.60   36.61 -0.4481E+01  0.1284E+02 "alien"
   39    4   15    1.57   34.01 -0.4162E+01  0.1133E+02 "alien"
   54    4    7    1.39   33.07 -0.3167E+01  0.9420E+01 "alien"
   47    1    5    1.61   32.94 -0.4031E+01  0.9168E+01 "alien"
   28   23   21    1.40   32.37 -0.3100E+01  0.8352E+01 "alien"
   41   31    3    1.50   32.31  0.2624E+02  0.8988E+01 "alien"
   55    3    0    1.41   27.19 -0.2018E+02  0.1362E+02 "alien"
   29   27   19    1.41   24.04 -0.2302E+01  0.9638E+01 "alien"
   40    3   13    1.62   22.95 -0.2809E+01  0.8802E+01 "alien"
   39    2   12    1.69   22.41 -0.1752E+02  0.7825E+01 "alien"
   34   26   10    1.64   22.18 -0.1734E+02  0.9175E+01 "alien"
   42    2   11    1.63   21.35 -0.1669E+02  0.7982E+01 "alien"
   46    1    6    1.63   20.96 -0.2565E+01  0.1343E+02 "alien"
   36   27    9    1.60   20.40 -0.2497E+01  0.9358E+01 "alien"
   23   21   19    1.58   19.55 -0.2393E+01  0.7730E+01 "alien"
   40   29    0    1.57   18.85 -0.8764E+01  0.1178E+02 "alien"
   36   26    7    1.67   17.98 -0.1406E+02  0.7853E+01 "alien"
   37   26    4    1.69   17.46 -0.1365E+02  0.8071E+01 "alien"
   37   26    5    1.67   17.39 -0.1360E+02  0.8275E+01 "alien"
   26   20   13    1.85   17.39 -0.1532E+02  0.7008E+01 "alien"
   50    3    8    1.48   17.03  0.1383E+02  0.8891E+01 "alien"
   30   28   19    1.39   17.01 -0.1629E+01  0.8803E+01 "alien"
   29    1   13    1.99   16.10 -0.1818E+02  0.8006E+01 "alien"
   32   23    9    1.79   15.35 -0.1352E+02  0.1323E+02 "alien"
   39   31    9    1.46   15.04  0.1221E+02  0.1036E+02 "alien"
   42   32    2    1.46   14.69  0.1193E+02  0.9444E+01 "alien"
   44    2   10    1.60   14.09 -0.1724E+01  0.8180E+01 "alien"
   53    2    1    1.46   13.94  0.1133E+02  0.1154E+02 "alien"
   28   23   14    1.69   13.67 -0.1068E+02  0.7238E+01 "alien"
   31   22   10    1.81   13.43 -0.1183E+02  0.1336E+02 "alien"
   35    1   12    1.83   13.31 -0.1173E+02  0.7982E+01 "alien"
   27   19   11    1.96   13.07 -0.1476E+02  0.8351E+01 "alien"
   39   28    0    1.62   12.69 -0.3737E+01  0.1201E+02 "alien"
   35   23    0    1.85   11.92 -0.2400E+02  0.1149E+02 "alien"
   41    1    4    1.86   11.76 -0.1036E+02  0.8530E+01 "alien"
   34   22    1    1.91   11.73 -0.1324E+02  0.8930E+01 "alien"
   30   21   23    1.32   11.41  0.3619E+02  0.1475E+02 "alien"
   34    4   16    1.66   11.21 -0.8763E+01  0.7661E+01 "alien"
   42   32    1    1.47   11.02  0.8951E+01  0.8867E+01 "alien"
   43    1    1    1.80   11.00 -0.9689E+01  0.7971E+01 "alien"
   33   29   14    1.49   10.87  0.8831E+01  0.1522E+02 "alien"
   40    2    0    1.94   10.85 -0.2607E+02  0.1354E+02 "alien"
   51    2    3    1.51   10.77  0.8747E+01  0.9467E+01 "alien"
   32   24   11    1.70   10.77 -0.8417E+01  0.8295E+01 "alien"
   39   28    3    1.60   10.59 -0.1296E+01  0.8557E+01 "alien"
   52    2    1    1.49   10.46  0.8497E+01  0.8921E+01 "alien"
   34   25    9    1.69   10.33 -0.8074E+01  0.7821E+01 "alien"
   33   22    6    1.87   10.26 -0.9043E+01  0.7858E+01 "alien"
   19   14   20    1.67   10.15 -0.7938E+01  0.8469E+01 "alien"
   49    3    9    1.48   10.14  0.8239E+01  0.8916E+01 "alien"


 List of     78 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   30   21   23    1.32   11.41  0.3619E+02  0.1475E+02
   36   33   14    1.38  117.61 -0.1126E+02  0.1149E+02
   27   18   23    1.38  184.46 -0.1767E+02  0.1108E+02
   28   21   22    1.38   71.05 -0.6804E+01  0.8606E+01
   30   12   23    1.38   83.68 -0.8014E+01  0.1166E+02
   37   33   13    1.39   86.06 -0.8241E+01  0.9897E+01
   54    4    7    1.39   33.07 -0.3167E+01  0.9420E+01
   30   28   19    1.39   17.01 -0.1629E+01  0.8803E+01
   48    6   14    1.39  161.28 -0.1545E+02  0.1257E+02
   38   33   12    1.39   61.49 -0.5889E+01  0.9086E+01
   36   32   14    1.39   61.95 -0.5933E+01  0.1598E+02
   28   23   21    1.40   32.37 -0.3100E+01  0.8352E+01
   43    7   17    1.40  208.00 -0.1992E+02  0.9055E+01
   37    9   20    1.41   86.03 -0.8239E+01  0.8586E+01
   52    4    9    1.41  113.23 -0.1084E+02  0.9436E+01
   55    3    1    1.41  223.89 -0.2144E+02  0.9549E+01
   55    3    0    1.41   27.19 -0.2018E+02  0.1362E+02
   31   30   17    1.41   65.02 -0.6227E+01  0.9380E+01
   40   33    9    1.41  137.97 -0.1321E+02  0.8941E+01
   29   27   19    1.41   24.04 -0.2302E+01  0.9638E+01
   53    2    1    1.46   13.94  0.1133E+02  0.1154E+02
   39   31    9    1.46   15.04  0.1221E+02  0.1036E+02
   42   32    2    1.46   14.69  0.1193E+02  0.9444E+01
   42   32    1    1.47   11.02  0.8951E+01  0.8867E+01
   50    3    8    1.48   17.03  0.1383E+02  0.8891E+01
   49    3    9    1.48   10.14  0.8239E+01  0.8916E+01
   33   29   14    1.49   10.87  0.8831E+01  0.1522E+02
   52    2    1    1.49   10.46  0.8497E+01  0.8921E+01
   41   31    3    1.50   32.31  0.2624E+02  0.8988E+01
   51    2    3    1.51   10.77  0.8747E+01  0.9467E+01
   40   29    0    1.57   18.85 -0.8764E+01  0.1178E+02
   39    4   15    1.57   34.01 -0.4162E+01  0.1133E+02
   45    2   10    1.57   40.17 -0.4917E+01  0.9679E+01
   39   29    4    1.57  110.92 -0.1358E+02  0.1554E+02
   33    6   18    1.58   50.38 -0.6166E+01  0.7334E+01
   23   21   19    1.58   19.55 -0.2393E+01  0.7730E+01
   32   26   13    1.59   94.83 -0.1161E+02  0.8196E+01
   38   28    6    1.59   43.40 -0.5311E+01  0.1188E+02
   25   24   17    1.60   36.61 -0.4481E+01  0.1284E+02
   47    1    6    1.60   80.26 -0.9823E+01  0.1182E+02
   35   27   10    1.60  178.97 -0.2190E+02  0.8296E+01
   45    2    9    1.60  132.46 -0.1621E+02  0.1494E+02
   36   27    9    1.60   20.40 -0.2497E+01  0.9358E+01
   39   28    3    1.60   10.59 -0.1296E+01  0.8557E+01
   44    2   10    1.60   14.09 -0.1724E+01  0.8180E+01
   43    2   11    1.60   67.26 -0.8232E+01  0.9476E+01
   38   28    5    1.61   80.02 -0.9793E+01  0.1495E+02
   39   28    2    1.61   94.55 -0.1157E+02  0.8210E+01
   47    1    5    1.61   32.94 -0.4031E+01  0.9168E+01
   30   25   14    1.61   53.21 -0.6512E+01  0.7638E+01
   39   28    0    1.62   12.69 -0.3737E+01  0.1201E+02
   40    3   13    1.62   22.95 -0.2809E+01  0.8802E+01
   36   27    8    1.62  114.24 -0.1398E+02  0.9762E+01
   46    1    6    1.63   20.96 -0.2565E+01  0.1343E+02
   42    2   11    1.63   21.35 -0.1669E+02  0.7982E+01
   34   26   10    1.64   22.18 -0.1734E+02  0.9175E+01
   33   25   11    1.65   44.34 -0.3466E+02  0.1322E+02
   34    4   16    1.66   11.21 -0.8763E+01  0.7661E+01
   36   26    7    1.67   17.98 -0.1406E+02  0.7853E+01
   19   14   20    1.67   10.15 -0.7938E+01  0.8469E+01
   37   26    5    1.67   17.39 -0.1360E+02  0.8275E+01
   37   26    4    1.69   17.46 -0.1365E+02  0.8071E+01
   34   25    9    1.69   10.33 -0.8074E+01  0.7821E+01
   39    2   12    1.69   22.41 -0.1752E+02  0.7825E+01
   28   23   14    1.69   13.67 -0.1068E+02  0.7238E+01
   32   24   11    1.70   10.77 -0.8417E+01  0.8295E+01
   32   23    9    1.79   15.35 -0.1352E+02  0.1323E+02
   43    1    1    1.80   11.00 -0.9689E+01  0.7971E+01
   31   22   10    1.81   13.43 -0.1183E+02  0.1336E+02
   35    1   12    1.83   13.31 -0.1173E+02  0.7982E+01
   35   23    0    1.85   11.92 -0.2400E+02  0.1149E+02
   26   20   13    1.85   17.39 -0.1532E+02  0.7008E+01
   41    1    4    1.86   11.76 -0.1036E+02  0.8530E+01
   33   22    6    1.87   10.26 -0.9043E+01  0.7858E+01
   34   22    1    1.91   11.73 -0.1324E+02  0.8930E+01
   40    2    0    1.94   10.85 -0.2607E+02  0.1354E+02
   27   19   11    1.96   13.07 -0.1476E+02  0.8351E+01
   29    1   13    1.99   16.10 -0.1818E+02  0.8006E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        42
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1166
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1074
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       92


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
