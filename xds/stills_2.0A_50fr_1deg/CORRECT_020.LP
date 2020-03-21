
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.446    77.446    38.762  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      20      20
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      20
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000859  0.000696
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
 REFLECTING_RANGE_E.S.D.=   0.115
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

    1398 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      20 ...      20
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1398 REFLECTIONS ACCEPTED



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

 *  44        aP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  31        aP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  35        mP          0.0      77.4   38.8   77.4  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  34        mP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.5  109.5   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  32        oP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  23        oC          0.0     109.5  109.5   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.5  109.5   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.4   77.4   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.7   38.8   77.4  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.7   38.8   77.4  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.8  159.7   77.4  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.8  159.7   77.4  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.8  159.7   77.4  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.8  159.7   77.4  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    27        mC        500.0     159.7   38.8  109.5  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.8  109.5  116.2  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.8  159.7  159.7  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.5  116.2   38.8  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        998.0     173.2   77.4   38.8  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.0      77.4  173.2   38.8  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.0      77.4  173.2   38.8  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.0      77.4   77.4   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.2  116.2  116.2  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.6   86.6  116.2 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.6   86.6  109.5  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.6   86.6  116.2 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.5   86.6   86.6  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.6   86.6  109.5  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.8   86.6  248.0 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.6   86.6   77.4  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.6   86.6   77.4  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.8   77.4  177.5  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.6   86.6  177.5 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.6   86.6   77.4  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.6   86.6   86.6  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.5  109.5   38.8  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.8   77.4  177.5 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.6  177.5   77.4 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    5.0      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.788   -23.525     0.029
 UNIT_CELL_B-AXIS=    23.525    73.788     0.015
 UNIT_CELL_C-AXIS=    -0.016    -0.003    38.763
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       249 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.24
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.380    77.380    38.726  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  3.0E-02 3.0E-02 1.2E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012923  0.012923  0.025822  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.743   -23.445     0.035
 COORDINATES OF UNIT CELL B-AXIS    23.445    73.743     0.011
 COORDINATES OF UNIT CELL C-AXIS    -0.018     0.000    38.726
 CRYSTAL MOSAICITY (DEGREES)     0.115
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006408  0.001860
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.91   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1261.21   1242.20
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.14
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.38    77.38    38.73  90.000  90.000  90.000


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
  0.05  0.15        84        108.333       1.568
  0.15  0.25       171        198.758       0.855
  0.25  0.35       189        143.415       1.184
  0.35  0.45       151        105.245       1.614
  0.45  0.55       168        279.856       0.607
  0.55  0.65       195        102.059       1.664
  0.65  0.75       195        201.264       0.844
  0.75  0.85       158        167.964       1.011
  0.85  0.95        87        214.387       0.792
 
 -0.05  0.95      1398        169.842       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               223.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.909
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    19.1 XMAX=    19.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59560 NYBIN=    9
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        496


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               223.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.902
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     9.4 XMAX=  2455.2 NXBIN=    3
 YMIN=     8.0 YMAX=  2517.5 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        496


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               223.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.896
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    19.1 XMAX=    19.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1263
 DETECTOR_SURFACE_POSITION=    1648    1689
 DETECTOR_SURFACE_POSITION=     817    1689
 DETECTOR_SURFACE_POSITION=     817     837
 DETECTOR_SURFACE_POSITION=    1648     837
 DETECTOR_SURFACE_POSITION=    2172    1662
 DETECTOR_SURFACE_POSITION=    1622    2227
 DETECTOR_SURFACE_POSITION=     843    2227
 DETECTOR_SURFACE_POSITION=     292    1662
 DETECTOR_SURFACE_POSITION=     292     863
 DETECTOR_SURFACE_POSITION=     843     298
 DETECTOR_SURFACE_POSITION=    1622     298
 DETECTOR_SURFACE_POSITION=    2172     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        496


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
 9.130E-01  1.300E-03   29.03


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

   24.470   3.152    23.82   1.38      2.07      2.08      62      89      36
    3.152   2.238    11.25   1.28      4.86      4.06     182     214      16
    2.238   1.830     2.03   0.90     21.16     21.90     208     243       0
    1.830   1.586    -0.07   0.86    -99.90    -99.90     244     287       0
    1.586   1.419    -0.11   0.84    -99.90    -99.90     278     335       2
    1.419   1.296     0.12   0.64    -99.90    -99.90     124     173       2
   --------------------------------------------------------------------------
   24.470   1.296     3.72   0.94      6.53      6.62    1098    1341      56
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      20      20)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.682   0.9192E+01  0.1200E+02         0.77           1 
   47    0    0     1.646  -0.2014E+02  0.1191E+02        -1.69           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS   -11.4% OF MEAN INTENSITY


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

     3.84          60      39      1285        3.0%       2.1%      2.0%       42   29.64      2.9%    99.9*     0    0.000       0
     2.74          89      54      2098        2.6%       3.1%      2.7%       70   23.31      4.4%    99.6*     0    0.000       0
     2.24         154      88      2665        3.3%       5.8%      5.0%      132   12.73      8.1%    99.2*     0    0.000       0
     1.94         158      89      3101        2.9%      14.7%     15.0%      138    4.18     20.8%    96.6*     0    0.000       0
     1.74         168      98      3481        2.8%     -99.9%    -99.9%      140    0.00    -99.9%   -16.2      0    0.000       0
     1.59         204     117      3848        3.0%     -99.9%    -99.9%      174    0.00    -99.9%     9.3      0    0.000       0
     1.47         216     128      4143        3.1%     -99.9%    -99.9%      176    0.00    -99.9%    -6.7      0    0.000       0
     1.37         210     119      4429        2.7%     -99.9%    -99.9%      182    0.00    -99.9%    15.9      0    0.000       0
     1.30          82      60      4743        1.3%     -99.9%    -99.9%       44    0.00    -99.9%    -9.1      0    0.000       0
    total        1341     792     29793        2.7%       6.5%      6.6%     1098    4.89      9.2%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1398
 NUMBER OF REJECTED MISFITS                                56
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1341
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    792



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

   WILSON LINE (using all data) : A=   7.545 B=  25.973 CORRELATION=  0.91
      #      RES      SS        <I>       log(<I>)       BO
      17     7.257   0.005  2.0116E+03       7.607      -6.5
      18     4.585   0.012  2.6254E+03       7.873     -13.8
      21     3.511   0.020  1.0735E+03       6.979      14.0
      28     2.996   0.028  3.8948E+02       5.965      28.4
      36     2.649   0.036  3.3691E+02       5.820      24.2
      47     2.402   0.043  1.6823E+02       5.125      27.9
      38     2.208   0.051  1.0227E+02       4.628      28.5
      44     2.056   0.059  3.8288E+01       3.645      33.0
      44     1.932   0.067 -2.9597E-01    -999.999   -1000.0
      52     1.830   0.075 -4.8441E-01    -999.999   -1000.0
      48     1.740   0.083 -5.5835E-01    -999.999   -1000.0
      52     1.666   0.090 -1.7327E+00    -999.999   -1000.0
      52     1.597   0.098  9.5844E-01      -0.042      38.7
      65     1.536   0.106 -1.9626E+00    -999.999   -1000.0
      71     1.481   0.114 -1.4174E+00    -999.999   -1000.0
      55     1.434   0.122 -1.6866E+00    -999.999   -1000.0
      45     1.390   0.129  2.3678E+00       0.862      25.8
      41     1.352   0.137 -3.7795E-02    -999.999   -1000.0
      18     1.315   0.145  7.4958E+00       2.014      19.1


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.257        1.836        1.240        0.668
       3     4.585        1.074        0.605        0.254
       2     3.511        0.390        0.085        0.013
       4     2.996        1.072        0.466        0.153
       3     2.649        0.311        0.074        0.012
       5     2.402        0.992        0.551        0.237
       2     2.208        1.258        0.505        0.147
       7     2.056        0.387        0.150        0.046
       7     1.740       20.245       39.926      112.089
       3     1.597       83.639      245.006      795.847
       6     1.481       22.860       -3.721       61.824
       3     1.434       25.956       48.323       93.723
      82   overall        3.409        9.148       29.197


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.257        0.455        0.318        0.179
      15     4.585        1.026        1.050        0.898
      19     3.511        0.814        0.645        0.475
      24     2.996        0.779        0.590        0.383
      33     2.649        1.163        1.317        1.378
      42     2.402        0.889        0.732        0.513
      36     2.208        0.819        0.646        0.454
      37     2.056        1.387        1.793        2.217
      37     1.932       25.415      107.970      852.795
      49     1.597       43.210      -15.434      730.449
      45     1.390        7.854        8.788       37.046
      40     1.352       47.917       29.893     1111.643
      18     1.315        1.567        1.540        1.574
     710   overall        3.818       -0.174       53.080

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.257  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.32  0.30  0.29
       3     4.585  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       2     3.511  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       4     2.996  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.649  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       5     2.402  0.00  0.00  0.00  0.00  0.77  0.71  1.00  0.95  0.91  0.88
       2     2.208  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       7     2.056  1.73  1.24  1.03  1.21  1.10  1.02  0.96  0.91  0.87  1.05
       7     1.740  2.30  1.65  1.37  1.21  1.10  1.02  0.96  1.14  1.09  1.05
       3     1.597  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       6     1.481  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       3     1.434  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      82   overall  0.74  0.53  0.44  0.41  0.42  0.41  0.41  0.43  0.43  0.43


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.257  1.75  2.76  1.93  1.52  1.91  1.66  1.66  1.51  1.40  1.32
      15     4.585  1.40  1.10  1.03  0.81  1.36  1.18  1.06  1.33  1.24  1.16
      19     3.511  0.00  0.58  0.41  0.32  0.67  0.82  1.05  0.96  1.06  1.08
      24     2.996  1.31  1.61  1.13  1.01  1.16  1.11  1.08  1.21  1.12  1.12
      33     2.649  0.64  1.50  1.29  1.19  1.00  0.94  1.02  1.05  1.12  1.10
      42     2.402  1.25  0.92  1.01  1.01  1.09  1.11  1.04  1.04  1.12  1.05
      36     2.208  0.88  1.07  0.75  1.01  1.06  0.92  0.88  0.91  1.03  1.01
      37     2.056  1.70  1.04  0.83  0.82  0.96  1.08  0.97  0.93  0.96  0.90
      37     1.932  5.40  2.83  1.98  1.56  1.31  1.14  1.02  0.93  0.87  0.81
      49     1.597  4.07  2.25  1.57  1.30  1.14  1.00  0.89  0.82  0.76  0.71
      45     1.390  3.50  1.84  1.29  1.01  0.90  0.89  0.79  0.73  0.67  0.67
      40     1.352  4.99  2.62  1.83  1.44  1.21  1.05  0.94  0.86  0.80  0.75
      18     1.315  2.92  1.53  1.07  1.01  0.99  0.86  0.77  0.81  0.84  0.79
     710   overall  1.48  0.98  0.73  0.64  0.63  0.59  0.56  0.54  0.54  0.52


 List of     88 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   49    7   15    1.34  525.13 -0.1985E+02  0.1494E+02 "alien"
   50    6   13    1.37  485.55 -0.1835E+02  0.9421E+01 "alien"
   53    5   10    1.36  437.42 -0.1653E+02  0.9595E+01 "alien"
   52    5   11    1.37  374.79 -0.1417E+02  0.1509E+02 "alien"
   46    8   17    1.34  271.83 -0.1027E+02  0.9182E+01 "alien"
   32   12   23    1.35  268.54 -0.1015E+02  0.1011E+02 "alien"
   30   25   21    1.35  264.29 -0.9989E+01  0.1185E+02 "alien"
   28   19   23    1.36  183.45 -0.6933E+01  0.7570E+01 "alien"
   45    8   17    1.36  176.99 -0.6689E+01  0.1048E+02 "alien"
   32   30   19    1.33  165.95 -0.6272E+01  0.1363E+02 "alien"
   34   33   16    1.35  164.89 -0.6232E+01  0.9351E+01 "alien"
   53    6   11    1.34  147.17 -0.5562E+01  0.9723E+01 "alien"
   42    9   19    1.35  146.22 -0.5526E+01  0.8509E+01 "alien"
   37   34   14    1.35  135.08 -0.5105E+01  0.9702E+01 "alien"
   38   34   13    1.35   93.21 -0.3523E+01  0.9650E+01 "alien"
   28   17   24    1.33   91.87 -0.3472E+01  0.1245E+02 "alien"
   43    9   19    1.33   71.78 -0.2713E+01  0.1433E+02 "alien"
   39   34   12    1.36   71.26 -0.2693E+01  0.1520E+02 "alien"
   41    1    3    1.87   54.16 -0.2623E+02  0.1470E+02 "alien"
   31   22   10    1.80   43.49 -0.2107E+02  0.1379E+02 "alien"
   36   24    1    1.79   39.90 -0.1933E+02  0.8505E+01 "alien"
   30    2   14    1.88   33.79 -0.1000E+02  0.6826E+01 "alien"
   38    2    4    1.99   32.50 -0.9620E+01  0.9993E+01 "alien"
   40    2    1    1.93   30.53 -0.9035E+01  0.1034E+02 "alien"
   34    3   15    1.70   28.98 -0.1618E+02  0.1182E+02 "alien"
   35   33   15    1.37   26.51 -0.1002E+01  0.1408E+02 "alien"
   36    1    8    1.96   24.00 -0.7103E+01  0.8484E+01 "alien"
   33   33   17    1.34   23.13 -0.4736E+01  0.1324E+02 "alien"
   39    1    6    1.90   22.13 -0.6549E+01  0.8055E+01 "alien"
   32   23   10    1.75   20.91 -0.1168E+02  0.1520E+02 "alien"
   15   10   19    1.84   20.62 -0.9989E+01  0.6826E+01 "alien"
   35   24    6    1.75   20.41 -0.1140E+02  0.8528E+01 "alien"
   29   22   12    1.78   20.05 -0.1119E+02  0.7585E+01 "alien"
   29    2   14    1.92   19.34 -0.5723E+01  0.7320E+01 "alien"
   34   23    5    1.83   18.04 -0.8741E+01  0.8760E+01 "alien"
   39   29    4    1.57   17.51  0.1678E+02  0.1005E+02 "alien"
   30    1   13    1.95   17.15 -0.5077E+01  0.1126E+02 "alien"
   24   19   14    1.87   17.08 -0.8272E+01  0.1123E+02 "alien"
   28   22   13    1.76   16.96 -0.9471E+01  0.7441E+01 "alien"
   35   25    7    1.71   16.96 -0.9467E+01  0.8089E+01 "alien"
   15   12   18    1.90   16.88 -0.4995E+01  0.5141E+01 "alien"
   38   28    6    1.59   16.75  0.1605E+02  0.1266E+02 "alien"
   31   26   14    1.57   16.68  0.1599E+02  0.8460E+01 "alien"
   34    7   19    1.50   15.89 -0.2253E+02  0.1078E+02 "alien"
   34   22    1    1.91   15.70 -0.4647E+01  0.8509E+01 "alien"
   37    0    9    1.88   15.24 -0.1300E+02  0.9015E+01 "alien"
   31    3   15    1.79   15.04 -0.7283E+01  0.7162E+01 "alien"
   32   31   18    1.35   14.80 -0.5595E+00  0.1394E+02 "alien"
   33    0   11    1.95   14.59 -0.1187E+02  0.1315E+02 "alien"
   20   19   18    1.71   14.58 -0.8144E+01  0.6749E+01 "alien"
   37   25    1    1.73   14.47 -0.8080E+01  0.8647E+01 "alien"
   21   19   19    1.63   13.89 -0.2407E+02  0.1022E+02 "alien"
   30   28   16    1.49   13.85 -0.1963E+02  0.1259E+02 "alien"
   40    3   13    1.62   13.54  0.1297E+02  0.8489E+01 "alien"
   37   25    2    1.73   13.52 -0.7546E+01  0.9508E+01 "alien"
   33    1   12    1.90   13.37 -0.3958E+01  0.7775E+01 "alien"
   45    4   13    1.48   13.34 -0.1891E+02  0.8822E+01 "alien"
   37    8   19    1.44   13.22 -0.2229E+02  0.1200E+02 "alien"
   25   18   12    1.98   13.13 -0.3887E+01  0.1178E+02 "alien"
   37   27    7    1.62   13.09  0.1254E+02  0.8804E+01 "alien"
   35    0   10    1.92   13.08 -0.9573E+01  0.7905E+01 "alien"
   35   23    0    1.85   13.08 -0.1567E+02  0.1191E+02 "alien"
   40   31    5    1.50   12.73 -0.1804E+02  0.1683E+02 "alien"
   14    8   19    1.88   12.53 -0.6069E+01  0.5707E+01 "alien"
   26   19   12    1.93   12.46 -0.3686E+01  0.7713E+01 "alien"
   33   11   22    1.38   12.12  0.2869E+02  0.1122E+02 "alien"
   32   23    9    1.79   11.89 -0.5757E+01  0.8468E+01 "alien"
   42    1    2    1.83   11.71 -0.5671E+01  0.8891E+01 "alien"
   53    2    1    1.46   11.59 -0.1643E+02  0.1117E+02 "alien"
   43    1    8    1.69   11.24 -0.1948E+02  0.1399E+02 "alien"
   54    3    1    1.43   11.11 -0.1874E+02  0.9227E+01 "alien"
   40    2    0    1.93   11.07 -0.7077E+01  0.1299E+02 "alien"
   37    1    7    1.96   11.00 -0.3256E+01  0.8399E+01 "alien"
   46    5   14    1.43   10.97 -0.1850E+02  0.1369E+02 "alien"
   24   14   23    1.44   10.96 -0.1848E+02  0.8478E+01 "alien"
   45    2    9    1.60   10.93  0.1048E+02  0.9471E+01 "alien"
   50    1    2    1.54   10.90 -0.2138E+02  0.9235E+01 "alien"
   19   18   18    1.74   10.83 -0.6046E+01  0.6087E+01 "alien"
   39    0    8    1.84   10.53 -0.1065E+02  0.1362E+02 "alien"
   45    2   10    1.57   10.51  0.1007E+02  0.9625E+01 "alien"
   32   26   13    1.59   10.44  0.1000E+02  0.8540E+01 "alien"
   40   32    7    1.46   10.43 -0.1758E+02  0.9930E+01 "alien"
   36    2   13    1.74   10.39 -0.5799E+01  0.7297E+01 "alien"
   41    5   15    1.52   10.35 -0.2031E+02  0.1425E+02 "alien"
   39   28    0    1.61   10.30  0.2029E+02  0.1303E+02 "alien"
   40   30    5    1.52   10.28 -0.2018E+02  0.8887E+01 "alien"
   33   27   13    1.55   10.21 -0.2005E+02  0.1366E+02 "alien"
   41   31    3    1.50   10.08 -0.1429E+02  0.8998E+01 "alien"


 List of     88 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   28   17   24    1.33   91.87 -0.3472E+01  0.1245E+02
   43    9   19    1.33   71.78 -0.2713E+01  0.1433E+02
   32   30   19    1.33  165.95 -0.6272E+01  0.1363E+02
   49    7   15    1.34  525.13 -0.1985E+02  0.1494E+02
   46    8   17    1.34  271.83 -0.1027E+02  0.9182E+01
   33   33   17    1.34   23.13 -0.4736E+01  0.1324E+02
   53    6   11    1.34  147.17 -0.5562E+01  0.9723E+01
   37   34   14    1.35  135.08 -0.5105E+01  0.9702E+01
   42    9   19    1.35  146.22 -0.5526E+01  0.8509E+01
   30   25   21    1.35  264.29 -0.9989E+01  0.1185E+02
   32   12   23    1.35  268.54 -0.1015E+02  0.1011E+02
   32   31   18    1.35   14.80 -0.5595E+00  0.1394E+02
   38   34   13    1.35   93.21 -0.3523E+01  0.9650E+01
   34   33   16    1.35  164.89 -0.6232E+01  0.9351E+01
   28   19   23    1.36  183.45 -0.6933E+01  0.7570E+01
   39   34   12    1.36   71.26 -0.2693E+01  0.1520E+02
   45    8   17    1.36  176.99 -0.6689E+01  0.1048E+02
   53    5   10    1.36  437.42 -0.1653E+02  0.9595E+01
   35   33   15    1.37   26.51 -0.1002E+01  0.1408E+02
   52    5   11    1.37  374.79 -0.1417E+02  0.1509E+02
   50    6   13    1.37  485.55 -0.1835E+02  0.9421E+01
   33   11   22    1.38   12.12  0.2869E+02  0.1122E+02
   54    3    1    1.43   11.11 -0.1874E+02  0.9227E+01
   46    5   14    1.43   10.97 -0.1850E+02  0.1369E+02
   24   14   23    1.44   10.96 -0.1848E+02  0.8478E+01
   37    8   19    1.44   13.22 -0.2229E+02  0.1200E+02
   40   32    7    1.46   10.43 -0.1758E+02  0.9930E+01
   53    2    1    1.46   11.59 -0.1643E+02  0.1117E+02
   45    4   13    1.48   13.34 -0.1891E+02  0.8822E+01
   30   28   16    1.49   13.85 -0.1963E+02  0.1259E+02
   41   31    3    1.50   10.08 -0.1429E+02  0.8998E+01
   40   31    5    1.50   12.73 -0.1804E+02  0.1683E+02
   34    7   19    1.50   15.89 -0.2253E+02  0.1078E+02
   41    5   15    1.52   10.35 -0.2031E+02  0.1425E+02
   40   30    5    1.52   10.28 -0.2018E+02  0.8887E+01
   50    1    2    1.54   10.90 -0.2138E+02  0.9235E+01
   33   27   13    1.55   10.21 -0.2005E+02  0.1366E+02
   45    2   10    1.57   10.51  0.1007E+02  0.9625E+01
   39   29    4    1.57   17.51  0.1678E+02  0.1005E+02
   31   26   14    1.57   16.68  0.1599E+02  0.8460E+01
   32   26   13    1.59   10.44  0.1000E+02  0.8540E+01
   38   28    6    1.59   16.75  0.1605E+02  0.1266E+02
   45    2    9    1.60   10.93  0.1048E+02  0.9471E+01
   39   28    0    1.61   10.30  0.2029E+02  0.1303E+02
   37   27    7    1.62   13.09  0.1254E+02  0.8804E+01
   40    3   13    1.62   13.54  0.1297E+02  0.8489E+01
   21   19   19    1.63   13.89 -0.2407E+02  0.1022E+02
   43    1    8    1.69   11.24 -0.1948E+02  0.1399E+02
   34    3   15    1.70   28.98 -0.1618E+02  0.1182E+02
   20   19   18    1.71   14.58 -0.8144E+01  0.6749E+01
   35   25    7    1.71   16.96 -0.9467E+01  0.8089E+01
   37   25    2    1.73   13.52 -0.7546E+01  0.9508E+01
   37   25    1    1.73   14.47 -0.8080E+01  0.8647E+01
   19   18   18    1.74   10.83 -0.6046E+01  0.6087E+01
   36    2   13    1.74   10.39 -0.5799E+01  0.7297E+01
   32   23   10    1.75   20.91 -0.1168E+02  0.1520E+02
   35   24    6    1.75   20.41 -0.1140E+02  0.8528E+01
   28   22   13    1.76   16.96 -0.9471E+01  0.7441E+01
   29   22   12    1.78   20.05 -0.1119E+02  0.7585E+01
   32   23    9    1.79   11.89 -0.5757E+01  0.8468E+01
   36   24    1    1.79   39.90 -0.1933E+02  0.8505E+01
   31    3   15    1.79   15.04 -0.7283E+01  0.7162E+01
   31   22   10    1.80   43.49 -0.2107E+02  0.1379E+02
   34   23    5    1.83   18.04 -0.8741E+01  0.8760E+01
   42    1    2    1.83   11.71 -0.5671E+01  0.8891E+01
   39    0    8    1.84   10.53 -0.1065E+02  0.1362E+02
   15   10   19    1.84   20.62 -0.9989E+01  0.6826E+01
   35   23    0    1.85   13.08 -0.1567E+02  0.1191E+02
   24   19   14    1.87   17.08 -0.8272E+01  0.1123E+02
   41    1    3    1.87   54.16 -0.2623E+02  0.1470E+02
   14    8   19    1.88   12.53 -0.6069E+01  0.5707E+01
   37    0    9    1.88   15.24 -0.1300E+02  0.9015E+01
   30    2   14    1.88   33.79 -0.1000E+02  0.6826E+01
   39    1    6    1.90   22.13 -0.6549E+01  0.8055E+01
   33    1   12    1.90   13.37 -0.3958E+01  0.7775E+01
   15   12   18    1.90   16.88 -0.4995E+01  0.5141E+01
   34   22    1    1.91   15.70 -0.4647E+01  0.8509E+01
   29    2   14    1.92   19.34 -0.5723E+01  0.7320E+01
   35    0   10    1.92   13.08 -0.9573E+01  0.7905E+01
   26   19   12    1.93   12.46 -0.3686E+01  0.7713E+01
   40    2    1    1.93   30.53 -0.9035E+01  0.1034E+02
   40    2    0    1.93   11.07 -0.7077E+01  0.1299E+02
   30    1   13    1.95   17.15 -0.5077E+01  0.1126E+02
   33    0   11    1.95   14.59 -0.1187E+02  0.1315E+02
   37    1    7    1.96   11.00 -0.3256E+01  0.8399E+01
   36    1    8    1.96   24.00 -0.7103E+01  0.8484E+01
   25   18   12    1.98   13.13 -0.3887E+01  0.1178E+02
   38    2    4    1.99   32.50 -0.9620E+01  0.9993E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        33
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1398
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1291
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      107


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
