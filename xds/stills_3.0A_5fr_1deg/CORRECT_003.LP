
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.324    77.324    38.701  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_3.0A_5fr_1deg/stills_3.0A_5fr_1deg_???_001.img SMV
 DATA_RANGE=       3       3
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       3
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001013  0.000645
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
 BEAM_DIVERGENCE_E.S.D.=   0.042
 REFLECTING_RANGE_E.S.D.=   0.134
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

    1335 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       3 ...       3
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1335 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               40
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

 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0  1  0  0  1  0  0  0
 *  44        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  25        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    27        mC        500.0     159.4   38.7  109.4  90.0 133.3  76.0    0  2  1  0  0  0  1  0  1 -1  0  0
    19        oI        500.0      38.7  109.4  116.0  90.0  70.5  90.0    0  0  1  0  1 -1  0  0  1  1  1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0 -1  0  0  2  1  0  2  0  1  0
    18        tI        625.0     109.4  116.0   38.7  70.5  90.0  90.0    1 -1  0  0 -1 -1 -1  0  0  0 -1  0
    30        mC        998.0      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0 -2  1  0  0  0  0  1  0
    41        mC        998.0     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    40        oC        998.0      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0    1 -1 -1  0 -1  1 -1  0  1  1  1  0
     2        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    1  0  1  0  0 -1  1  0  1 -1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.4   86.5   86.5  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.7   86.5  247.6 102.1  99.0 116.6    0  0 -1  0  0  1  1  0  3 -1  1  0
    10        mC        999.0      86.5   86.5   77.3  90.0  90.0 126.8    0  1 -1  0  0 -1 -1  0 -1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.2  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.5   86.5  177.2 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1  1  0  0 -1 -1  0  1  0  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
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

  *   0.00        0    4.9      15    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.608   -23.684     0.112
 UNIT_CELL_B-AXIS=    23.684    73.608     0.044
 UNIT_CELL_C-AXIS=    -0.060    -0.004    38.701
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       125 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.22
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     76.820    76.820    38.463  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  5.4E-01 5.4E-01 2.6E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.013017  0.013017  0.025999  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.183   -23.356     0.099
 COORDINATES OF UNIT CELL B-AXIS    23.356    73.183     0.043
 COORDINATES OF UNIT CELL C-AXIS    -0.054    -0.006    38.463
 CRYSTAL MOSAICITY (DEGREES)     0.134
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999984  0.005206  0.002320
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.89   1232.82
 DETECTOR ORIGIN (PIXELS) AT                     1260.54   1240.34
 CRYSTAL TO DETECTOR DISTANCE (mm)       248.66
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    76.82    76.82    38.46  90.000  90.000  90.000


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
  0.05  0.15        74         98.810       1.326
  0.15  0.25       168        115.144       1.138
  0.25  0.35       173        117.949       1.111
  0.35  0.45       146         45.203       2.898
  0.45  0.55       183        249.280       0.526
  0.55  0.65       187         89.942       1.456
  0.65  0.75       180        146.117       0.897
  0.75  0.85       159         63.489       2.063
  0.85  0.95        65        344.469       0.380
 
 -0.05  0.95      1335        131.000       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT                99.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.785
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     2.1 XMAX=     2.9 NXBIN=    1
 YMIN= 0.00169 YMAX= 0.60411 NYBIN=    4
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        226


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT                99.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.782
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    36.8 XMAX=  2441.4 NXBIN=    2
 YMIN=     7.1 YMAX=  2519.7 NYBIN=    2
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        226


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT                97.1
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.784
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     2.1 XMAX=     2.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1239    1263
 DETECTOR_SURFACE_POSITION=    1648    1690
 DETECTOR_SURFACE_POSITION=     831    1690
 DETECTOR_SURFACE_POSITION=     831     837
 DETECTOR_SURFACE_POSITION=    1648     837
 DETECTOR_SURFACE_POSITION=    2163    1663
 DETECTOR_SURFACE_POSITION=    1622    2229
 DETECTOR_SURFACE_POSITION=     856    2229
 DETECTOR_SURFACE_POSITION=     315    1663
 DETECTOR_SURFACE_POSITION=     315     863
 DETECTOR_SURFACE_POSITION=     856     298
 DETECTOR_SURFACE_POSITION=    1622     298
 DETECTOR_SURFACE_POSITION=    2163     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        226


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
 8.630E-01  8.520E-04   36.88


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

   24.293   3.130    27.82   1.91      2.11      1.68      50      85      30
    3.130   2.222     1.46   1.07     20.26     19.53     170     213       0
    2.222   1.817    -0.04   0.89    -99.90    -99.90     194     237       0
    1.817   1.575    -0.07   0.92    -99.90    -99.90     226     275       0
    1.575   1.409    -0.03   0.79   4327.29   4842.12     254     324       0
    1.409   1.287     0.08   1.02    -99.90    -99.90     118     168       2
   --------------------------------------------------------------------------
   24.293   1.287     2.04   0.97      8.18      7.95    1012    1302      32
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       3       3)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.670   0.2350E+02  0.1163E+02         2.02           1 
   47    0    0     1.634   0.4221E+01  0.1145E+02         0.37           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     3.0% OF MEAN INTENSITY


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

     3.82          50      35      1285        2.7%       1.9%      1.6%       30   34.21      2.7%   100.0*     0    0.000       0
     2.72          92      60      2099        2.9%       3.8%      2.9%       64   16.58      5.3%    99.7*     0    0.000       0
     2.22         156      93      2663        3.5%    1024.4%   1014.1%      126    0.10   1448.7%   -26.6      0    0.000       0
     1.93         153      89      3102        2.9%     -99.9%    -99.9%      128    0.00    -99.9%   -11.2      0    0.000       0
     1.72         164      97      3481        2.8%     -99.9%    -99.9%      134    0.00    -99.9%     8.2      0    0.000       0
     1.57         195     116      3844        3.0%     -99.9%    -99.9%      158    0.00    -99.9%   -20.0      0    0.000       0
     1.46         212     129      4146        3.1%     -99.9%    -99.9%      166    0.00    -99.9%    -5.9      0    0.000       0
     1.36         199     119      4434        2.7%     -99.9%    -99.9%      160    0.00    -99.9%    -8.8      0    0.000       0
     1.29          81      58      4730        1.2%     -99.9%    -99.9%       46    0.00    -99.9%     3.9      0    0.000       0
    total        1302     796     29784        2.7%       8.2%      7.9%     1012    2.75     11.6%   100.0*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1335
 NUMBER OF REJECTED MISFITS                                32
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1302
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    796



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

   WILSON LINE (using all data) : A=   5.725 B=  17.962 CORRELATION=  0.45
      #      RES      SS        <I>       log(<I>)       BO
      14     7.329   0.005  2.4859E+03       7.818    -224.9
      18     4.615   0.012  2.4699E+03       7.812     -88.9
      26     3.473   0.021  9.7827E+02       6.886     -28.0
      29     2.967   0.028  1.7818E+02       5.183       9.5
      38     2.622   0.036 -1.2245E+00    -999.999   -1000.0
      50     2.382   0.044  1.1441E+00       0.135      63.4
      39     2.188   0.052  8.7670E-01      -0.132      56.1
      42     2.042   0.060 -1.1269E+00    -999.999   -1000.0
      44     1.918   0.068  2.2278E+00       0.801      36.2
      54     1.816   0.076 -3.6100E-01    -999.999   -1000.0
      47     1.727   0.084 -1.2991E-01    -999.999   -1000.0
      50     1.654   0.091 -1.1868E+00    -999.999   -1000.0
      52     1.586   0.099 -1.5005E+00    -999.999   -1000.0
      66     1.526   0.107 -7.3622E-01    -999.999   -1000.0
      69     1.471   0.116 -3.2015E-01    -999.999   -1000.0
      55     1.422   0.124 -6.0114E-01    -999.999   -1000.0
      47     1.381   0.131 -2.2897E-01    -999.999   -1000.0
      39     1.343   0.139  5.7967E+01       4.060       6.0
      17     1.305   0.147 -3.5382E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.329        1.496        0.918        0.448
       3     4.615        1.253        0.762        0.344
       3     3.473        0.339        0.084        0.015
       4     2.967        1.407        1.154        0.677
       3     2.622        5.539        7.181        6.899
       5     2.382       34.647      102.979      259.893
       2     2.188      173.497     1083.431     4916.911
       8     2.042       22.153       27.009      103.540
       7     1.918        5.466        4.853        5.725
      11     1.816       35.657       82.750      345.373
       6     1.471       50.928      -20.430      587.948
       2     1.343        0.004        0.000        0.000
      86   overall        6.702       31.715      129.993


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

       9     7.329        0.415        0.252        0.119
      15     4.615        1.132        1.284        1.189
      23     3.473        0.874        0.694        0.510
      25     2.967        1.971        3.048        3.822
      45     2.382       15.813       17.300      188.684
      37     2.188       44.915        0.711      867.315
      37     1.918        8.541        9.444       56.684
      38     1.727       31.922      -38.280      623.121
      37     1.343       19.486      246.363     2337.895
     710   overall        4.930       14.624      182.123

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.329  0.00  0.00  0.00  0.00  0.00  0.00  0.33  0.32  0.30  0.59
       3     4.615  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.473  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.967  3.02  2.17  1.80  1.58  1.44  1.34  1.26  1.19  1.14  1.10
       3     2.622  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       5     2.382  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.59
       2     2.188  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.042  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       7     1.918  2.30  1.65  1.37  1.21  1.10  1.02  0.96  0.91  0.87  0.84
      11     1.816  2.20  1.58  1.31  1.15  1.05  0.97  0.91  0.87  0.83  0.80
       6     1.471  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       2     1.343  2.01  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      86   overall  1.31  0.98  0.81  0.71  0.65  0.60  0.58  0.55  0.53  0.55


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       9     7.329  3.50  3.06  2.14  2.02  1.98  1.72  1.55  1.41  1.31  1.23
      15     4.615  1.40  1.84  1.29  1.21  1.19  1.33  1.19  1.21  1.24  1.27
      23     3.473  0.46  0.72  0.50  0.53  0.66  0.87  0.86  0.95  0.88  1.10
      25     2.967  5.88  3.53  2.47  2.06  1.73  1.51  1.35  1.31  1.21  1.20
      45     2.382  5.37  2.82  1.97  1.62  1.41  1.28  1.19  1.09  1.05  1.02
      37     2.188  5.40  2.83  1.98  1.56  1.31  1.14  1.02  0.98  0.91  0.86
      37     1.918  3.98  2.09  1.46  1.23  1.10  0.96  0.86  0.79  0.73  0.68
      38     1.727  6.08  3.19  2.23  1.76  1.47  1.28  1.15  1.05  0.98  0.92
      37     1.343  6.25  4.77  3.65  2.95  2.47  2.16  1.93  1.77  1.64  1.54
     710   overall  1.78  1.08  0.77  0.64  0.55  0.50  0.46  0.43  0.40  0.39


 List of    128 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   23   21   16    1.72  166.76 -0.2166E+02  0.1100E+02 "alien"
   35   24    6    1.74  155.61 -0.2022E+02  0.8071E+01 "alien"
   33   32   16    1.37  116.89 -0.2676E+02  0.1566E+02 "alien"
   38    1   11    1.75   78.76 -0.1023E+02  0.6944E+01 "alien"
   37    2   13    1.70   71.99 -0.9353E+01  0.1055E+02 "alien"
   30   12   23    1.37   69.16 -0.1584E+02  0.7289E+01 "alien"
   37    7   18    1.48   69.07 -0.2211E+02  0.1315E+02 "alien"
   30   26   20    1.36   68.39 -0.1566E+02  0.1419E+02 "alien"
   40   31    5    1.49   62.84 -0.2012E+02  0.1436E+02 "alien"
   54    4    6    1.39   61.18 -0.1401E+02  0.8839E+01 "alien"
   31    4   16    1.72   59.49 -0.7728E+01  0.1035E+02 "alien"
   26   24   19    1.48   51.90 -0.1661E+02  0.9811E+01 "alien"
   51    5   11    1.38   50.81 -0.1163E+02  0.8475E+01 "alien"
   33    8   20    1.47   49.20 -0.1575E+02  0.1152E+02 "alien"
   45    7   16    1.38   49.16 -0.1126E+02  0.8976E+01 "alien"
   25   21   15    1.73   48.99 -0.6365E+01  0.1021E+02 "alien"
   37   31   12    1.43   47.82 -0.2875E+02  0.1529E+02 "alien"
   54    4    7    1.37   44.34 -0.1015E+02  0.9136E+01 "alien"
   51    2    3    1.49   42.71 -0.1367E+02  0.8661E+01 "alien"
   31   22   10    1.79   41.82 -0.1510E+02  0.1254E+02 "alien"
   38    2   12    1.71   41.54 -0.5396E+01  0.1209E+02 "alien"
   36   25    5    1.71   41.16 -0.5347E+01  0.7485E+01 "alien"
   46    4   12    1.48   41.02 -0.1313E+02  0.8389E+01 "alien"
   38    7   18    1.46   40.32 -0.1291E+02  0.7528E+01 "alien"
   54    5    8    1.36   37.96  0.2200E+04  0.7234E+02 "alien"
   51    2    5    1.48   34.62 -0.1108E+02  0.9096E+01 "alien"
   33   23    8    1.77   34.40 -0.1242E+02  0.7292E+01 "alien"
   26   20   13    1.84   34.03 -0.1229E+02  0.6561E+01 "alien"
   40   31    7    1.46   31.56 -0.1010E+02  0.9425E+01 "alien"
   31   30   18    1.37   30.30 -0.6937E+01  0.8717E+01 "alien"
   17   15   18    1.81   29.89 -0.1079E+02  0.5037E+01 "alien"
   41   33    6    1.42   29.67 -0.1783E+02  0.1532E+02 "alien"
   42    8   18    1.38   29.15 -0.6675E+01  0.8096E+01 "alien"
   36   24    1    1.77   29.11 -0.1051E+02  0.7714E+01 "alien"
   35   23    3    1.82   28.91 -0.1044E+02  0.7890E+01 "alien"
   52    5   10    1.37   28.28 -0.6475E+01  0.8977E+01 "alien"
   33   24    9    1.72   28.27 -0.3672E+01  0.7415E+01 "alien"
   34   31   15    1.40   27.45 -0.6286E+01  0.8401E+01 "alien"
   40   30    2    1.53   27.36 -0.2014E+02  0.9829E+01 "alien"
   36   24    3    1.76   27.14 -0.3526E+01  0.9283E+01 "alien"
   29    9   21    1.48   26.73 -0.8559E+01  0.6567E+01 "alien"
   41    1    4    1.84   26.45 -0.9548E+01  0.7997E+01 "alien"
   36   10   21    1.37   26.40 -0.6045E+01  0.7514E+01 "alien"
   30    8   20    1.52   25.89 -0.1906E+02  0.7111E+01 "alien"
   36   25    4    1.72   25.40 -0.3300E+01  0.7910E+01 "alien"
   30   28   16    1.48   25.26 -0.8086E+01  0.1089E+02 "alien"
   41   32    4    1.46   24.91 -0.7974E+01  0.9577E+01 "alien"
   37   31   11    1.45   23.73 -0.7597E+01  0.9501E+01 "alien"
   43    0    7    1.70   23.65 -0.1739E+02  0.1354E+02 "alien"
   35   23    2    1.83   23.58 -0.8512E+01  0.7560E+01 "alien"
   27   18   23    1.37   23.25 -0.5323E+01  0.1017E+02 "alien"
   44    0    4    1.72   23.09 -0.1618E+02  0.1363E+02 "alien"
   26   22   15    1.69   22.65 -0.2942E+01  0.7473E+01 "alien"
   34   29   13    1.49   22.47 -0.7193E+01  0.8121E+01 "alien"
   40    1    5    1.86   22.32 -0.8056E+01  0.1144E+02 "alien"
   36   30   12    1.46   22.08 -0.7068E+01  0.8299E+01 "alien"
   22   11   23    1.47   22.01 -0.7046E+01  0.7356E+01 "alien"
   28   21   12    1.81   21.06 -0.7601E+01  0.6963E+01 "alien"
   38   30    9    1.49   20.95 -0.6707E+01  0.8066E+01 "alien"
   44    1    0    1.75   20.63 -0.1177E+02  0.1259E+02 "alien"
   43    1    1    1.78   20.60 -0.7437E+01  0.8069E+01 "alien"
   30   18    2    2.18   20.08  0.1761E+02  0.8430E+01 "alien"
   41   31    2    1.49   19.56 -0.6262E+01  0.8369E+01 "alien"
   37    9   20    1.39   19.16 -0.4388E+01  0.7395E+01 "alien"
   54    3    4    1.41   19.15 -0.1151E+02  0.9319E+01 "alien"
   41   31    3    1.48   18.91 -0.6055E+01  0.8178E+01 "alien"
   37    8   19    1.43   18.88 -0.1135E+02  0.1079E+02 "alien"
   30   29   18    1.39   18.84 -0.4314E+01  0.1330E+02 "alien"
   41   31    0    1.49   18.36 -0.2150E+02  0.1146E+02 "alien"
   26   21   21    1.43   18.09 -0.1087E+02  0.9697E+01 "alien"
   49    5   12    1.40   17.98 -0.1081E+02  0.8679E+01 "alien"
   44    0    5    1.70   17.62 -0.7900E+01  0.7593E+01 "alien"
   37   33   13    1.37   17.23 -0.3945E+01  0.8495E+01 "alien"
   26   25   17    1.55   17.18 -0.1265E+02  0.7888E+01 "alien"
   38   31   10    1.45   17.12 -0.5482E+01  0.8375E+01 "alien"
   27    2    8    2.44   16.93  0.1937E+02  0.1145E+02 "alien"
   26    3    5    2.74   16.87 -0.2066E+02  0.7192E+01 "alien"
   53    3    6    1.41   16.64 -0.1000E+02  0.1562E+02 "alien"
   35    2    4    2.14   16.56  0.1452E+02  0.1627E+02 "alien"
   20   12    9    2.61   16.09 -0.1971E+02  0.1098E+02 "alien"
   27   20   12    1.86   16.05 -0.5795E+01  0.6971E+01 "alien"
   37   32   12    1.41   15.94 -0.9584E+01  0.8875E+01 "alien"
   54    3    3    1.41   15.05 -0.9048E+01  0.8608E+01 "alien"
   29   17    4    2.22   14.95  0.1310E+02  0.1411E+02 "alien"
   45    3   11    1.53   14.91 -0.1098E+02  0.8103E+01 "alien"
   33   23    7    1.80   14.90 -0.5378E+01  0.1379E+02 "alien"
   22   12    6    2.77   14.83 -0.1816E+02  0.1042E+02 "alien"
   31   28   15    1.49   14.70 -0.4707E+01  0.1324E+02 "alien"
   30    9   21    1.47   14.28 -0.4571E+01  0.9340E+01 "alien"
   39    0    8    1.82   14.09 -0.1350E+02  0.1169E+02 "alien"
   39   32   10    1.42   14.06 -0.8450E+01  0.9371E+01 "alien"
   26    9   21    1.53   13.90 -0.1024E+02  0.8221E+01 "alien"
   27    2    7    2.52   13.86 -0.1697E+02  0.1276E+02 "alien"
   28   27   17    1.49   13.83 -0.4428E+01  0.7510E+01 "alien"
   26   19   22    1.41   13.83 -0.8314E+01  0.9162E+01 "alien"
   21   12   22    1.53   13.62 -0.1002E+02  0.7977E+01 "alien"
   34   22    4    1.86   13.37 -0.4826E+01  0.9654E+01 "alien"
   50    1    2    1.53   13.25 -0.9755E+01  0.9608E+01 "alien"
   27    9   21    1.52   13.03 -0.9592E+01  0.7959E+01 "alien"
   34    3    0    2.25   12.98  0.2794E+02  0.1303E+02 "alien"
   28   16    1    2.38   12.89  0.1475E+02  0.7608E+01 "alien"
   50    5   12    1.38   12.85 -0.2943E+01  0.9730E+01 "alien"
   47    1    6    1.58   12.85 -0.1929E+02  0.1148E+02 "alien"
   42    2   11    1.62   12.82 -0.1522E+02  0.7504E+01 "alien"
   21   17   14    1.98   12.63  0.2813E+02  0.9465E+01 "alien"
   36   26    6    1.67   12.49 -0.1482E+02  0.1276E+02 "alien"
   28   17    5    2.24   12.29  0.1078E+02  0.9997E+01 "alien"
   43    1    3    1.77   11.86 -0.1541E+01  0.8469E+01 "alien"
   54    3    0    1.42   11.72 -0.1588E+02  0.1253E+02 "alien"
   35    2   13    1.76   11.71 -0.1521E+01  0.6968E+01 "alien"
   20   15   13    2.13   11.62  0.1019E+02  0.9321E+01 "alien"
   34   26   10    1.63   11.43 -0.1357E+02  0.7552E+01 "alien"
   36   32   13    1.40   11.41 -0.6856E+01  0.1476E+02 "alien"
   35    2    5    2.11   11.09 -0.1250E+02  0.1171E+02 "alien"
   46    2    9    1.55   11.04 -0.8127E+01  0.7879E+01 "alien"
   36    0   10    1.87   10.89 -0.8389E+01  0.7579E+01 "alien"
   34    1   12    1.85   10.80 -0.3900E+01  0.6695E+01 "alien"
   45    1    8    1.61   10.69 -0.1603E+02  0.9535E+01 "alien"
   30   18    3    2.16   10.64  0.9325E+01  0.8275E+01 "alien"
   30   28   19    1.37   10.62 -0.2431E+01  0.7740E+01 "alien"
   39   29    4    1.56   10.46 -0.1569E+02  0.8698E+01 "alien"
   30   18    1    2.19   10.43  0.9143E+01  0.9332E+01 "alien"
   50    2    4    1.52   10.41 -0.7667E+01  0.1526E+02 "alien"
   14   11   13    2.44   10.34  0.1183E+02  0.8478E+01 "alien"
   48    2    7    1.54   10.32 -0.7599E+01  0.8342E+01 "alien"
   27    3    5    2.65   10.31 -0.1263E+02  0.8045E+01 "alien"
   38    9   20    1.38   10.15 -0.2324E+01  0.1189E+02 "alien"
   40    2   12    1.65   10.13 -0.1202E+02  0.7734E+01 "alien"


 List of    128 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   54    5    8    1.36   37.96  0.2200E+04  0.7234E+02
   30   26   20    1.36   68.39 -0.1566E+02  0.1419E+02
   27   18   23    1.37   23.25 -0.5323E+01  0.1017E+02
   36   10   21    1.37   26.40 -0.6045E+01  0.7514E+01
   30   12   23    1.37   69.16 -0.1584E+02  0.7289E+01
   31   30   18    1.37   30.30 -0.6937E+01  0.8717E+01
   33   32   16    1.37  116.89 -0.2676E+02  0.1566E+02
   37   33   13    1.37   17.23 -0.3945E+01  0.8495E+01
   52    5   10    1.37   28.28 -0.6475E+01  0.8977E+01
   54    4    7    1.37   44.34 -0.1015E+02  0.9136E+01
   30   28   19    1.37   10.62 -0.2431E+01  0.7740E+01
   38    9   20    1.38   10.15 -0.2324E+01  0.1189E+02
   42    8   18    1.38   29.15 -0.6675E+01  0.8096E+01
   51    5   11    1.38   50.81 -0.1163E+02  0.8475E+01
   50    5   12    1.38   12.85 -0.2943E+01  0.9730E+01
   45    7   16    1.38   49.16 -0.1126E+02  0.8976E+01
   54    4    6    1.39   61.18 -0.1401E+02  0.8839E+01
   37    9   20    1.39   19.16 -0.4388E+01  0.7395E+01
   30   29   18    1.39   18.84 -0.4314E+01  0.1330E+02
   34   31   15    1.40   27.45 -0.6286E+01  0.8401E+01
   49    5   12    1.40   17.98 -0.1081E+02  0.8679E+01
   36   32   13    1.40   11.41 -0.6856E+01  0.1476E+02
   54    3    4    1.41   19.15 -0.1151E+02  0.9319E+01
   26   19   22    1.41   13.83 -0.8314E+01  0.9162E+01
   37   32   12    1.41   15.94 -0.9584E+01  0.8875E+01
   53    3    6    1.41   16.64 -0.1000E+02  0.1562E+02
   54    3    3    1.41   15.05 -0.9048E+01  0.8608E+01
   39   32   10    1.42   14.06 -0.8450E+01  0.9371E+01
   54    3    0    1.42   11.72 -0.1588E+02  0.1253E+02
   41   33    6    1.42   29.67 -0.1783E+02  0.1532E+02
   37   31   12    1.43   47.82 -0.2875E+02  0.1529E+02
   26   21   21    1.43   18.09 -0.1087E+02  0.9697E+01
   37    8   19    1.43   18.88 -0.1135E+02  0.1079E+02
   37   31   11    1.45   23.73 -0.7597E+01  0.9501E+01
   38   31   10    1.45   17.12 -0.5482E+01  0.8375E+01
   38    7   18    1.46   40.32 -0.1291E+02  0.7528E+01
   36   30   12    1.46   22.08 -0.7068E+01  0.8299E+01
   41   32    4    1.46   24.91 -0.7974E+01  0.9577E+01
   40   31    7    1.46   31.56 -0.1010E+02  0.9425E+01
   33    8   20    1.47   49.20 -0.1575E+02  0.1152E+02
   30    9   21    1.47   14.28 -0.4571E+01  0.9340E+01
   22   11   23    1.47   22.01 -0.7046E+01  0.7356E+01
   37    7   18    1.48   69.07 -0.2211E+02  0.1315E+02
   46    4   12    1.48   41.02 -0.1313E+02  0.8389E+01
   30   28   16    1.48   25.26 -0.8086E+01  0.1089E+02
   51    2    5    1.48   34.62 -0.1108E+02  0.9096E+01
   26   24   19    1.48   51.90 -0.1661E+02  0.9811E+01
   29    9   21    1.48   26.73 -0.8559E+01  0.6567E+01
   41   31    3    1.48   18.91 -0.6055E+01  0.8178E+01
   34   29   13    1.49   22.47 -0.7193E+01  0.8121E+01
   38   30    9    1.49   20.95 -0.6707E+01  0.8066E+01
   28   27   17    1.49   13.83 -0.4428E+01  0.7510E+01
   40   31    5    1.49   62.84 -0.2012E+02  0.1436E+02
   41   31    2    1.49   19.56 -0.6262E+01  0.8369E+01
   31   28   15    1.49   14.70 -0.4707E+01  0.1324E+02
   41   31    0    1.49   18.36 -0.2150E+02  0.1146E+02
   51    2    3    1.49   42.71 -0.1367E+02  0.8661E+01
   27    9   21    1.52   13.03 -0.9592E+01  0.7959E+01
   50    2    4    1.52   10.41 -0.7667E+01  0.1526E+02
   30    8   20    1.52   25.89 -0.1906E+02  0.7111E+01
   50    1    2    1.53   13.25 -0.9755E+01  0.9608E+01
   45    3   11    1.53   14.91 -0.1098E+02  0.8103E+01
   26    9   21    1.53   13.90 -0.1024E+02  0.8221E+01
   40   30    2    1.53   27.36 -0.2014E+02  0.9829E+01
   21   12   22    1.53   13.62 -0.1002E+02  0.7977E+01
   48    2    7    1.54   10.32 -0.7599E+01  0.8342E+01
   26   25   17    1.55   17.18 -0.1265E+02  0.7888E+01
   46    2    9    1.55   11.04 -0.8127E+01  0.7879E+01
   39   29    4    1.56   10.46 -0.1569E+02  0.8698E+01
   47    1    6    1.58   12.85 -0.1929E+02  0.1148E+02
   45    1    8    1.61   10.69 -0.1603E+02  0.9535E+01
   42    2   11    1.62   12.82 -0.1522E+02  0.7504E+01
   34   26   10    1.63   11.43 -0.1357E+02  0.7552E+01
   40    2   12    1.65   10.13 -0.1202E+02  0.7734E+01
   36   26    6    1.67   12.49 -0.1482E+02  0.1276E+02
   26   22   15    1.69   22.65 -0.2942E+01  0.7473E+01
   37    2   13    1.70   71.99 -0.9353E+01  0.1055E+02
   43    0    7    1.70   23.65 -0.1739E+02  0.1354E+02
   44    0    5    1.70   17.62 -0.7900E+01  0.7593E+01
   38    2   12    1.71   41.54 -0.5396E+01  0.1209E+02
   36   25    5    1.71   41.16 -0.5347E+01  0.7485E+01
   44    0    4    1.72   23.09 -0.1618E+02  0.1363E+02
   31    4   16    1.72   59.49 -0.7728E+01  0.1035E+02
   23   21   16    1.72  166.76 -0.2166E+02  0.1100E+02
   33   24    9    1.72   28.27 -0.3672E+01  0.7415E+01
   36   25    4    1.72   25.40 -0.3300E+01  0.7910E+01
   25   21   15    1.73   48.99 -0.6365E+01  0.1021E+02
   35   24    6    1.74  155.61 -0.2022E+02  0.8071E+01
   44    1    0    1.75   20.63 -0.1177E+02  0.1259E+02
   38    1   11    1.75   78.76 -0.1023E+02  0.6944E+01
   36   24    3    1.76   27.14 -0.3526E+01  0.9283E+01
   35    2   13    1.76   11.71 -0.1521E+01  0.6968E+01
   43    1    3    1.77   11.86 -0.1541E+01  0.8469E+01
   36   24    1    1.77   29.11 -0.1051E+02  0.7714E+01
   33   23    8    1.77   34.40 -0.1242E+02  0.7292E+01
   43    1    1    1.78   20.60 -0.7437E+01  0.8069E+01
   31   22   10    1.79   41.82 -0.1510E+02  0.1254E+02
   33   23    7    1.80   14.90 -0.5378E+01  0.1379E+02
   17   15   18    1.81   29.89 -0.1079E+02  0.5037E+01
   28   21   12    1.81   21.06 -0.7601E+01  0.6963E+01
   35   23    3    1.82   28.91 -0.1044E+02  0.7890E+01
   39    0    8    1.82   14.09 -0.1350E+02  0.1169E+02
   35   23    2    1.83   23.58 -0.8512E+01  0.7560E+01
   26   20   13    1.84   34.03 -0.1229E+02  0.6561E+01
   41    1    4    1.84   26.45 -0.9548E+01  0.7997E+01
   34    1   12    1.85   10.80 -0.3900E+01  0.6695E+01
   34   22    4    1.86   13.37 -0.4826E+01  0.9654E+01
   27   20   12    1.86   16.05 -0.5795E+01  0.6971E+01
   40    1    5    1.86   22.32 -0.8056E+01  0.1144E+02
   36    0   10    1.87   10.89 -0.8389E+01  0.7579E+01
   21   17   14    1.98   12.63  0.2813E+02  0.9465E+01
   35    2    5    2.11   11.09 -0.1250E+02  0.1171E+02
   20   15   13    2.13   11.62  0.1019E+02  0.9321E+01
   35    2    4    2.14   16.56  0.1452E+02  0.1627E+02
   30   18    3    2.16   10.64  0.9325E+01  0.8275E+01
   30   18    2    2.18   20.08  0.1761E+02  0.8430E+01
   30   18    1    2.19   10.43  0.9143E+01  0.9332E+01
   29   17    4    2.22   14.95  0.1310E+02  0.1411E+02
   28   17    5    2.24   12.29  0.1078E+02  0.9997E+01
   34    3    0    2.25   12.98  0.2794E+02  0.1303E+02
   28   16    1    2.38   12.89  0.1475E+02  0.7608E+01
   14   11   13    2.44   10.34  0.1183E+02  0.8478E+01
   27    2    8    2.44   16.93  0.1937E+02  0.1145E+02
   27    2    7    2.52   13.86 -0.1697E+02  0.1276E+02
   20   12    9    2.61   16.09 -0.1971E+02  0.1098E+02
   27    3    5    2.65   10.31 -0.1263E+02  0.8045E+01
   26    3    5    2.74   16.87 -0.2066E+02  0.7192E+01
   22   12    6    2.77   14.83 -0.1816E+02  0.1042E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        62
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1335
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1201
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      134


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
