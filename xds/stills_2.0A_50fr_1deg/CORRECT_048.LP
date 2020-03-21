
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.658    77.658    38.874  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      48      48
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      48
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001030  0.000644
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
 REFLECTING_RANGE_E.S.D.=   0.206
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

    1074 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      48 ...      48
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1074 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               32
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

 *  44        aP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  31        aP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  35        mP          0.0      77.7   38.9   77.7  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  33        mP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.8  109.8   38.9  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.7   77.7   38.9  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     160.1   38.9   77.7  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     160.1   38.9   77.7  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.9  160.1   77.7  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.9  160.1   77.7  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.9  160.1   77.7  90.0  90.0  75.9    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.9  160.1   77.7  90.0  90.0  75.9    0  0  1  0  0  2  1  0 -1  0  0  0
    27        mC        500.0     160.1   38.9  109.8  90.0 133.3  75.9    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.9  109.8  116.5  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.9  160.1  160.1  86.6 104.1 104.1    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.8  116.5   38.9  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        997.7     173.6   77.7   38.9  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        997.7      77.7  173.6   38.9  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        997.7      77.7  173.6   38.9  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        997.7      77.7   77.7   38.9  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.5  116.5  116.5  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.8   86.8  116.5 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.8   86.8  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.8   86.8  116.5 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.8   86.8   86.8  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.8   86.8  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.9   86.8  248.6 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.8   86.8   77.7  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.8   86.8   77.7  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.9   77.7  177.9  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.8   86.8  177.9 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.8   86.8   77.7  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.8   86.8   86.8  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.9  109.8   38.9  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.9   77.7  177.9 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.8  177.9   77.7 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    4.9      12    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    74.016   -23.504     0.057
 UNIT_CELL_B-AXIS=    23.504    74.016     0.023
 UNIT_CELL_C-AXIS=    -0.031    -0.003    38.874
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM        91 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.29
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.801    77.801    38.930  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.0E-01 2.0E-01 7.4E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012853  0.012853  0.025687  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    74.171   -23.487     0.008
 COORDINATES OF UNIT CELL B-AXIS    23.487    74.171     0.006
 COORDINATES OF UNIT CELL C-AXIS    -0.005    -0.002    38.930
 CRYSTAL MOSAICITY (DEGREES)     0.206
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999976  0.006657  0.001691
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.85   1232.92
 DETECTOR ORIGIN (PIXELS) AT                     1261.39   1242.63
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.81
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.80    77.80    38.93  90.000  90.000  90.000


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
  0.05  0.15         9         88.841       0.970
  0.15  0.25       116         43.946       1.961
  0.25  0.35       161         93.029       0.927
  0.35  0.45       156         32.880       2.621
  0.45  0.55       162        167.220       0.515
  0.55  0.65       180         55.021       1.567
  0.65  0.75       175        117.389       0.734
  0.75  0.85       110         70.982       1.214
  0.85  0.95         5        244.357       0.353
 
 -0.05  0.95      1074         86.192       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        5
 DEGREES OF FREEDOM OF CHI^2 FIT               143.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.790
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    47.1 XMAX=    47.9 NXBIN=    1
 YMIN= 0.00165 YMAX= 0.59443 NYBIN=    5
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        298


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               143.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.789
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    87.0 XMAX=  2449.8 NXBIN=    2
 YMIN=     7.1 YMAX=  2517.5 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        298


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               143.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.788
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    47.1 XMAX=    47.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1268    1262
 DETECTOR_SURFACE_POSITION=    1670    1689
 DETECTOR_SURFACE_POSITION=     867    1689
 DETECTOR_SURFACE_POSITION=     867     836
 DETECTOR_SURFACE_POSITION=    1670     836
 DETECTOR_SURFACE_POSITION=    2177    1662
 DETECTOR_SURFACE_POSITION=    1645    2227
 DETECTOR_SURFACE_POSITION=     892    2227
 DETECTOR_SURFACE_POSITION=     360    1662
 DETECTOR_SURFACE_POSITION=     360     863
 DETECTOR_SURFACE_POSITION=     892     297
 DETECTOR_SURFACE_POSITION=    1645     297
 DETECTOR_SURFACE_POSITION=    2177     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        298


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
 8.072E-01  2.916E-03   20.61


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

   24.603   2.884    12.49   1.55      3.41      3.12      78     100      12
    2.884   2.047     1.83   0.88     28.39     29.62     176     218       2
    2.047   1.673     0.07   0.84    649.85    699.42     216     255       4
    1.673   1.450    -0.05   0.65    -99.90    -99.90     232     290       0
    1.450   1.297     0.01   0.96    880.08    864.12     136     188       4
   --------------------------------------------------------------------------
   24.603   1.297     1.57   0.88     11.04     11.39     838    1051      22
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      48      48)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.691   0.6279E+00  0.1310E+02         0.05           1 
   47    0    0     1.655  -0.9452E+01  0.1210E+02        -0.78           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS   -10.6% OF MEAN INTENSITY


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

     3.85          53      34      1297        2.6%       3.0%      2.9%       38   19.63      4.3%    99.8*     0    0.000       0
     2.74          68      40      2129        1.9%       7.3%      6.1%       56    9.66     10.3%    99.0*     0    0.000       0
     2.24         128      77      2699        2.9%      23.3%     26.1%      102    2.64     32.9%    89.4*     0    0.000       0
     1.94         112      66      3139        2.1%     155.7%    164.8%       92    0.45    220.2%    23.2      0    0.000       0
     1.74         150      86      3535        2.4%   21325.9%  22103.2%      128   -0.02  30159.4%     5.7      0    0.000       0
     1.59         152      89      3883        2.3%    1887.0%   2154.3%      126    0.03   2668.7%    23.2      0    0.000       0
     1.47         168     101      4205        2.4%     -99.9%    -99.9%      134    0.00    -99.9%    22.2      0    0.000       0
     1.38         150      90      4496        2.0%   11914.3%  12334.5%      120   -0.01  16849.4%     8.1      0    0.000       0
     1.30          70      49      4760        1.0%     335.1%    296.3%       42    0.06    473.9%   -25.9      0    0.000       0
    total        1051     632     30143        2.1%      11.0%     11.4%      838    2.03     15.6%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1074
 NUMBER OF REJECTED MISFITS                                22
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1051
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    632



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

   WILSON LINE (using all data) : A=   6.107 B=  31.523 CORRELATION=  0.81
      #      RES      SS        <I>       log(<I>)       BO
      12     7.531   0.004  1.7347E+03       7.459    -153.3
      18     4.752   0.011  1.3614E+03       7.216     -50.1
      15     3.618   0.019  4.1707E+02       6.033       1.9
      23     2.986   0.028  7.6415E+01       4.336      31.6
      32     2.634   0.036  5.7136E+01       4.045      28.6
      39     2.413   0.043  2.0443E+01       3.018      36.0
      32     2.215   0.051  9.8463E+00       2.287      37.5
      33     2.060   0.059  2.6537E+00       0.976      43.5
      34     1.932   0.067 -1.6278E+00    -999.999   -1000.0
      40     1.830   0.075  6.4068E-01      -0.445      43.9
      40     1.743   0.082  2.3752E+00       0.865      31.8
      45     1.668   0.090  7.7279E-03      -4.863      61.1
      46     1.599   0.098 -7.3341E-01    -999.999   -1000.0
      48     1.537   0.106 -1.0405E+00    -999.999   -1000.0
      48     1.484   0.114  3.7601E-01      -0.978      31.2
      47     1.435   0.121  2.7180E-01      -1.303      30.5
      36     1.393   0.129 -1.3587E+00    -999.999   -1000.0
      33     1.352   0.137  7.9310E-02      -2.534      31.6
      11     1.313   0.145  4.5127E+00       1.507      15.9


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       2     7.531        3.656        3.231        2.097
       3     4.752        1.218        0.727        0.323
       1     3.618        0.480        0.115        0.020
       3     2.986        3.077        1.996        0.953
       2     2.634        0.014        0.001        0.000
       3     2.413        1.623        1.101        0.559
       2     2.215        1.051        0.527        0.189
       5     2.060        6.714       12.495       17.480
       7     1.932       62.666       99.328      920.130
       8     1.743        6.024        4.708        4.884
       4     1.599       93.919      395.358     1924.318
       5     1.537       39.045      131.952      345.861
       3     1.484      204.745      828.228     4052.270
       1     1.393       28.077       51.536       67.568
      72   overall    23038.762   276503.511808352996.675


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.531        0.625        0.499        0.312
      15     4.752        1.261        1.596        1.625
      14     3.618        1.211        1.682        2.042
      20     2.986        0.498        0.318        0.176
      30     2.634        1.462        2.295        3.431
      36     2.413        0.970        0.800        0.599
      30     2.215        1.608        1.361        1.139
      28     2.060        7.040        4.756       22.964
      32     1.830      110.237      344.265     4226.633
      32     1.743        8.609        6.831       27.367
      34     1.668   938675.696**************************
      45     1.484      497.665     1268.678    88842.548
      44     1.435      822.922    -7401.411   427276.379
      32     1.352    13732.506  -260094.775112267022.556
      11     1.313        1.367        0.790        0.625
     560   overall    57887.878 -9474992.425*************

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       2     7.531  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.73
       3     4.752  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       1     3.618  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.986  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.634  2.01  1.45  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       3     2.413  0.00  0.00  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     2.215  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.060  1.61  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
       7     1.932  2.30  1.65  1.37  1.21  1.10  1.02  0.96  0.91  0.87  1.05
       8     1.743  1.51  1.08  0.90  0.79  0.72  0.67  0.63  0.60  0.57  0.55
       4     1.599  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     1.537  3.22  2.31  1.92  1.69  1.54  1.78  1.67  1.59  1.52  1.47
       3     1.484  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       1     1.393  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      72   overall  1.12  0.84  0.77  0.68  0.61  0.59  0.56  0.53  0.51  0.55


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.531  3.15  3.31  2.31  2.12  1.78  1.55  1.59  1.45  1.35  1.27
      15     4.752  2.80  1.47  1.29  1.01  1.02  1.03  1.32  1.33  1.46  1.37
      14     3.618  0.75  0.39  0.28  1.08  0.91  0.79  1.28  1.17  1.20  1.24
      20     2.986  0.53  1.10  1.16  1.36  1.14  1.33  1.39  1.36  1.35  1.27
      30     2.634  1.40  0.92  0.90  1.01  1.02  1.11  1.13  1.03  1.01  1.11
      36     2.413  2.04  1.23  1.18  1.01  0.92  0.99  0.94  0.96  0.94  0.92
      30     2.215  4.20  2.39  1.80  1.42  1.19  1.03  0.93  0.85  0.84  0.79
      28     2.060  4.50  2.36  1.65  1.41  1.27  1.11  1.06  0.97  0.90  0.85
      32     1.830  5.58  2.93  2.05  1.61  1.35  1.18  1.06  0.96  0.90  0.84
      32     1.743  4.27  2.24  1.57  1.33  1.19  1.11  1.06  0.96  0.90  0.84
      34     1.668  4.64  2.43  1.70  1.34  1.12  0.98  0.88  0.80  0.74  0.70
      45     1.484  5.14  2.70  1.89  1.55  1.30  1.13  1.02  0.93  0.86  0.81
      44     1.435  4.30  2.26  1.58  1.24  1.04  0.91  0.81  0.74  0.69  0.65
      32     1.352  4.93  2.59  1.81  1.42  1.19  1.04  0.93  0.85  0.79  0.74
      11     1.313  0.96  1.00  1.05  1.10  0.92  1.01  0.90  0.99  0.92  0.86
     560   overall  2.72  1.53  1.14  0.98  0.85  0.79  0.76  0.71  0.68  0.65


 List of    115 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   37    3   14    1.67 2522.67  0.1949E+02  0.8363E+01 "alien"
   38    3   14    1.65 2192.26  0.1694E+02  0.9318E+01 "alien"
   37   26    4    1.69 1698.04  0.1312E+02  0.9029E+01 "alien"
   21   19   19    1.64 1624.55  0.1255E+02  0.1097E+02 "alien"
   38   26    1    1.69 1539.01  0.1189E+02  0.1114E+02 "alien"
   42    2   11    1.64 1509.99  0.1167E+02  0.8693E+01 "alien"
   34   26   10    1.65 1281.99  0.9907E+01  0.1476E+02 "alien"
   36    3   14    1.70 1265.83  0.9782E+01  0.1386E+02 "alien"
   34    4   16    1.66 1236.07  0.9552E+01  0.8746E+01 "alien"
   45    1    7    1.65 1119.84  0.8654E+01  0.9207E+01 "alien"
   37   27    6    1.64  903.75  0.6984E+01  0.1389E+02 "alien"
   35   25    8    1.70  634.07  0.4900E+01  0.9969E+01 "alien"
   37   26    5    1.68  522.87  0.4041E+01  0.8676E+01 "alien"
   42    1    9    1.70  346.73  0.2679E+01  0.8951E+01 "alien"
   41   10   20    1.34  335.30  0.2659E+02  0.1061E+02 "alien"
   35   26    9    1.65  315.73  0.2440E+01  0.8481E+01 "alien"
   29   22   22    1.36  287.24  0.2278E+02  0.1437E+02 "alien"
   32   30   19    1.34  223.34  0.1771E+02  0.1453E+02 "alien"
   52    6   12    1.35  210.56  0.1670E+02  0.1464E+02 "alien"
   48    8   16    1.34  178.82  0.1418E+02  0.1603E+02 "alien"
   50    7   14    1.35  177.63  0.1409E+02  0.1208E+02 "alien"
   28   23   14    1.70  174.55  0.1349E+01  0.7984E+01 "alien"
   35   11   22    1.36  141.42  0.1122E+02  0.1389E+02 "alien"
   49    7   15    1.34  122.34  0.9703E+01  0.1709E+02 "alien"
   54    5    9    1.36  108.47  0.8603E+01  0.1047E+02 "alien"
   32   25   12    1.65   98.80  0.7635E+00  0.8697E+01 "alien"
   53    4    7    1.42   93.38  0.2538E+02  0.1205E+02 "alien"
   38   34   13    1.36   93.03  0.7378E+01  0.9802E+01 "alien"
   30   25   21    1.36   88.50  0.7019E+01  0.9683E+01 "alien"
   46    8   17    1.35   76.23  0.6046E+01  0.1039E+02 "alien"
   42    6   16    1.46   72.16  0.2713E+02  0.1377E+02 "alien"
   36   31   13    1.44   58.55  0.1591E+02  0.9940E+01 "alien"
   31   29   16    1.46   57.06  0.2146E+02  0.1480E+02 "alien"
   52    2    3    1.49   54.27  0.2041E+02  0.9945E+01 "alien"
   36   26    7    1.67   52.54  0.4060E+00  0.8536E+01 "alien"
   33   30   15    1.45   51.32  0.1395E+02  0.1027E+02 "alien"
   31   30   17    1.42   49.37  0.1342E+02  0.9989E+01 "alien"
   54    3    2    1.43   48.29  0.1312E+02  0.9920E+01 "alien"
   51    4   10    1.42   45.44  0.1235E+02  0.1819E+02 "alien"
   41   32    5    1.47   43.91  0.1651E+02  0.1045E+02 "alien"
   46    0    1    1.69   41.99  0.1050E+02  0.8987E+01 "alien"
   34   33   16    1.36   40.96  0.3248E+01  0.1673E+02 "alien"
   26   23   20    1.47   40.19  0.1511E+02  0.1289E+02 "alien"
   34   29   13    1.51   37.21  0.1399E+02  0.1625E+02 "alien"
   40   31    7    1.48   37.15  0.1397E+02  0.9947E+01 "alien"
   37    4   15    1.63   37.10 -0.2721E+02  0.1360E+02 "alien"
   45    0    5    1.69   35.73  0.4785E+01  0.9521E+01 "alien"
   52    3    6    1.46   34.98  0.9506E+01  0.9356E+01 "alien"
   39   31    8    1.49   34.21  0.1286E+02  0.1141E+02 "alien"
   33    8   20    1.48   33.90  0.1275E+02  0.1328E+02 "alien"
   47    4   12    1.47   33.28  0.1251E+02  0.1055E+02 "alien"
   24   20   15    1.80   32.65  0.2092E+02  0.1322E+02 "alien"
   48    1    3    1.61   30.78 -0.2258E+02  0.1577E+02 "alien"
   53    6   11    1.35   30.00  0.2379E+01  0.1134E+02 "alien"
   35    5   17    1.59   29.59 -0.2170E+02  0.1310E+02 "alien"
   45    6   15    1.43   28.97  0.7875E+01  0.1092E+02 "alien"
   38   29    7    1.56   28.78 -0.2995E+02  0.1625E+02 "alien"
   51    3    8    1.45   28.10  0.7638E+01  0.1636E+02 "alien"
   46    0    2    1.68   28.05  0.1816E+01  0.8675E+01 "alien"
   47    5   13    1.44   27.93  0.7592E+01  0.1023E+02 "alien"
   45    0    4    1.70   27.21  0.1632E+01  0.8387E+01 "alien"
   41   32    6    1.46   25.82  0.7018E+01  0.9939E+01 "alien"
   29    9   21    1.50   25.54  0.9603E+01  0.8507E+01 "alien"
   53    3    4    1.45   25.24  0.6860E+01  0.1044E+02 "alien"
   52    4    8    1.43   24.76  0.6728E+01  0.1146E+02 "alien"
   35   23    3    1.84   24.16  0.1548E+02  0.8889E+01 "alien"
   43    1    3    1.79   22.78  0.1460E+02  0.9855E+01 "alien"
   38    8   19    1.43   22.72  0.6176E+01  0.9326E+01 "alien"
   29   28   17    1.48   22.67  0.8525E+01  0.1006E+02 "alien"
   45    9   18    1.33   22.16  0.1758E+01  0.1484E+02 "alien"
   35   27   10    1.60   21.61 -0.1585E+02  0.8883E+01 "alien"
   38   32   11    1.43   21.40  0.5817E+01  0.1026E+02 "alien"
   25   20   21    1.47   20.87  0.7846E+01  0.1194E+02 "alien"
   40    1    6    1.86   20.48  0.1312E+02  0.1016E+02 "alien"
   43    1    2    1.80   19.85  0.1272E+02  0.8511E+01 "alien"
   31    3   15    1.80   18.78  0.1203E+02  0.1162E+02 "alien"
   42    5   15    1.50   18.59  0.6989E+01  0.9260E+01 "alien"
   39   29    5    1.57   18.45 -0.1353E+02  0.9114E+01 "alien"
   30   22   11    1.80   18.23  0.1168E+02  0.8332E+01 "alien"
   42    7   17    1.43   17.94  0.4875E+01  0.1413E+02 "alien"
   48    3   10    1.49   17.92  0.6737E+01  0.1534E+02 "alien"
   33   23    8    1.80   17.76  0.1138E+02  0.8308E+01 "alien"
   41   31    4    1.50   16.94  0.6368E+01  0.9682E+01 "alien"
   36   24    2    1.79   16.71  0.1071E+02  0.8774E+01 "alien"
   54    3    3    1.43   16.47  0.4476E+01  0.9619E+01 "alien"
   47    2    8    1.57   16.42 -0.1708E+02  0.8956E+01 "alien"
   38   28    6    1.60   16.33 -0.1197E+02  0.1237E+02 "alien"
   38   28    7    1.58   16.23 -0.1191E+02  0.1624E+02 "alien"
   37   29    9    1.55   16.18 -0.1683E+02  0.1515E+02 "alien"
   32   27   14    1.54   15.28 -0.1590E+02  0.1005E+02 "alien"
   44    3   12    1.55   15.15 -0.1576E+02  0.1460E+02 "alien"
   54    3    0    1.44   14.83  0.1128E+02  0.1541E+02 "alien"
   28    9   21    1.52   14.78 -0.1538E+02  0.1122E+02 "alien"
   48    2    8    1.54   14.77 -0.1537E+02  0.9882E+01 "alien"
   46    0    0    1.69   14.68  0.6279E+00  0.1310E+02 "alien"
   37    1   11    1.81   14.14  0.9058E+01  0.8966E+01 "alien"
   36   27    9    1.61   14.03 -0.1029E+02  0.9509E+01 "alien"
   30    1   13    1.96   14.03 -0.2283E+02  0.1374E+02 "alien"
   40   32    8    1.45   13.88  0.3771E+01  0.9551E+01 "alien"
   42   32    3    1.46   13.82  0.5197E+01  0.1126E+02 "alien"
   39   29    6    1.55   13.26 -0.1380E+02  0.1018E+02 "alien"
   28   28   18    1.45   12.98  0.8657E+01  0.9696E+01 "alien"
   30   28   19    1.39   12.69 -0.1725E+02  0.9962E+01 "alien"
   52    2    0    1.50   12.18  0.1064E+02  0.1336E+02 "alien"
   30   25   14    1.62   12.09 -0.8865E+01  0.8120E+01 "alien"
   33   32   18    1.33   12.08  0.9584E+00  0.1344E+02 "alien"
   42   32    0    1.47   11.90  0.1020E+02  0.1342E+02 "alien"
   33   26   12    1.61   11.78 -0.8643E+01  0.9871E+01 "alien"
   27   22   21    1.43   11.35  0.3086E+01  0.9301E+01 "alien"
   35   31   14    1.43   10.88  0.2957E+01  0.9381E+01 "alien"
   47    0    1    1.65   10.27  0.1627E+00  0.9360E+01 "alien"
   40   29    3    1.56   10.22 -0.1064E+02  0.1135E+02 "alien"
   36   27    8    1.63   10.18 -0.7464E+01  0.1551E+02 "alien"
   47    3   10    1.52   10.17 -0.1058E+02  0.9448E+01 "alien"
   39   33   11    1.40   10.14 -0.1378E+02  0.1005E+02 "alien"


 List of    115 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   33   32   18    1.33   12.08  0.9584E+00  0.1344E+02
   45    9   18    1.33   22.16  0.1758E+01  0.1484E+02
   48    8   16    1.34  178.82  0.1418E+02  0.1603E+02
   41   10   20    1.34  335.30  0.2659E+02  0.1061E+02
   32   30   19    1.34  223.34  0.1771E+02  0.1453E+02
   49    7   15    1.34  122.34  0.9703E+01  0.1709E+02
   46    8   17    1.35   76.23  0.6046E+01  0.1039E+02
   50    7   14    1.35  177.63  0.1409E+02  0.1208E+02
   53    6   11    1.35   30.00  0.2379E+01  0.1134E+02
   52    6   12    1.35  210.56  0.1670E+02  0.1464E+02
   30   25   21    1.36   88.50  0.7019E+01  0.9683E+01
   35   11   22    1.36  141.42  0.1122E+02  0.1389E+02
   38   34   13    1.36   93.03  0.7378E+01  0.9802E+01
   34   33   16    1.36   40.96  0.3248E+01  0.1673E+02
   54    5    9    1.36  108.47  0.8603E+01  0.1047E+02
   29   22   22    1.36  287.24  0.2278E+02  0.1437E+02
   30   28   19    1.39   12.69 -0.1725E+02  0.9962E+01
   39   33   11    1.40   10.14 -0.1378E+02  0.1005E+02
   53    4    7    1.42   93.38  0.2538E+02  0.1205E+02
   51    4   10    1.42   45.44  0.1235E+02  0.1819E+02
   31   30   17    1.42   49.37  0.1342E+02  0.9989E+01
   52    4    8    1.43   24.76  0.6728E+01  0.1146E+02
   27   22   21    1.43   11.35  0.3086E+01  0.9301E+01
   35   31   14    1.43   10.88  0.2957E+01  0.9381E+01
   42    7   17    1.43   17.94  0.4875E+01  0.1413E+02
   54    3    3    1.43   16.47  0.4476E+01  0.9619E+01
   45    6   15    1.43   28.97  0.7875E+01  0.1092E+02
   38   32   11    1.43   21.40  0.5817E+01  0.1026E+02
   38    8   19    1.43   22.72  0.6176E+01  0.9326E+01
   54    3    2    1.43   48.29  0.1312E+02  0.9920E+01
   36   31   13    1.44   58.55  0.1591E+02  0.9940E+01
   54    3    0    1.44   14.83  0.1128E+02  0.1541E+02
   47    5   13    1.44   27.93  0.7592E+01  0.1023E+02
   33   30   15    1.45   51.32  0.1395E+02  0.1027E+02
   53    3    4    1.45   25.24  0.6860E+01  0.1044E+02
   40   32    8    1.45   13.88  0.3771E+01  0.9551E+01
   51    3    8    1.45   28.10  0.7638E+01  0.1636E+02
   28   28   18    1.45   12.98  0.8657E+01  0.9696E+01
   52    3    6    1.46   34.98  0.9506E+01  0.9356E+01
   41   32    6    1.46   25.82  0.7018E+01  0.9939E+01
   31   29   16    1.46   57.06  0.2146E+02  0.1480E+02
   42   32    3    1.46   13.82  0.5197E+01  0.1126E+02
   42    6   16    1.46   72.16  0.2713E+02  0.1377E+02
   41   32    5    1.47   43.91  0.1651E+02  0.1045E+02
   26   23   20    1.47   40.19  0.1511E+02  0.1289E+02
   47    4   12    1.47   33.28  0.1251E+02  0.1055E+02
   42   32    0    1.47   11.90  0.1020E+02  0.1342E+02
   25   20   21    1.47   20.87  0.7846E+01  0.1194E+02
   29   28   17    1.48   22.67  0.8525E+01  0.1006E+02
   40   31    7    1.48   37.15  0.1397E+02  0.9947E+01
   33    8   20    1.48   33.90  0.1275E+02  0.1328E+02
   52    2    3    1.49   54.27  0.2041E+02  0.9945E+01
   39   31    8    1.49   34.21  0.1286E+02  0.1141E+02
   48    3   10    1.49   17.92  0.6737E+01  0.1534E+02
   52    2    0    1.50   12.18  0.1064E+02  0.1336E+02
   41   31    4    1.50   16.94  0.6368E+01  0.9682E+01
   42    5   15    1.50   18.59  0.6989E+01  0.9260E+01
   29    9   21    1.50   25.54  0.9603E+01  0.8507E+01
   34   29   13    1.51   37.21  0.1399E+02  0.1625E+02
   28    9   21    1.52   14.78 -0.1538E+02  0.1122E+02
   47    3   10    1.52   10.17 -0.1058E+02  0.9448E+01
   48    2    8    1.54   14.77 -0.1537E+02  0.9882E+01
   32   27   14    1.54   15.28 -0.1590E+02  0.1005E+02
   37   29    9    1.55   16.18 -0.1683E+02  0.1515E+02
   44    3   12    1.55   15.15 -0.1576E+02  0.1460E+02
   39   29    6    1.55   13.26 -0.1380E+02  0.1018E+02
   38   29    7    1.56   28.78 -0.2995E+02  0.1625E+02
   40   29    3    1.56   10.22 -0.1064E+02  0.1135E+02
   47    2    8    1.57   16.42 -0.1708E+02  0.8956E+01
   39   29    5    1.57   18.45 -0.1353E+02  0.9114E+01
   38   28    7    1.58   16.23 -0.1191E+02  0.1624E+02
   35    5   17    1.59   29.59 -0.2170E+02  0.1310E+02
   38   28    6    1.60   16.33 -0.1197E+02  0.1237E+02
   35   27   10    1.60   21.61 -0.1585E+02  0.8883E+01
   36   27    9    1.61   14.03 -0.1029E+02  0.9509E+01
   48    1    3    1.61   30.78 -0.2258E+02  0.1577E+02
   33   26   12    1.61   11.78 -0.8643E+01  0.9871E+01
   30   25   14    1.62   12.09 -0.8865E+01  0.8120E+01
   37    4   15    1.63   37.10 -0.2721E+02  0.1360E+02
   36   27    8    1.63   10.18 -0.7464E+01  0.1551E+02
   42    2   11    1.64 1509.99  0.1167E+02  0.8693E+01
   21   19   19    1.64 1624.55  0.1255E+02  0.1097E+02
   37   27    6    1.64  903.75  0.6984E+01  0.1389E+02
   38    3   14    1.65 2192.26  0.1694E+02  0.9318E+01
   34   26   10    1.65 1281.99  0.9907E+01  0.1476E+02
   35   26    9    1.65  315.73  0.2440E+01  0.8481E+01
   32   25   12    1.65   98.80  0.7635E+00  0.8697E+01
   45    1    7    1.65 1119.84  0.8654E+01  0.9207E+01
   47    0    1    1.65   10.27  0.1627E+00  0.9360E+01
   34    4   16    1.66 1236.07  0.9552E+01  0.8746E+01
   36   26    7    1.67   52.54  0.4060E+00  0.8536E+01
   37    3   14    1.67 2522.67  0.1949E+02  0.8363E+01
   37   26    5    1.68  522.87  0.4041E+01  0.8676E+01
   46    0    2    1.68   28.05  0.1816E+01  0.8675E+01
   45    0    5    1.69   35.73  0.4785E+01  0.9521E+01
   38   26    1    1.69 1539.01  0.1189E+02  0.1114E+02
   46    0    1    1.69   41.99  0.1050E+02  0.8987E+01
   46    0    0    1.69   14.68  0.6279E+00  0.1310E+02
   37   26    4    1.69 1698.04  0.1312E+02  0.9029E+01
   35   25    8    1.70  634.07  0.4900E+01  0.9969E+01
   28   23   14    1.70  174.55  0.1349E+01  0.7984E+01
   45    0    4    1.70   27.21  0.1632E+01  0.8387E+01
   42    1    9    1.70  346.73  0.2679E+01  0.8951E+01
   36    3   14    1.70 1265.83  0.9782E+01  0.1386E+02
   36   24    2    1.79   16.71  0.1071E+02  0.8774E+01
   43    1    3    1.79   22.78  0.1460E+02  0.9855E+01
   24   20   15    1.80   32.65  0.2092E+02  0.1322E+02
   33   23    8    1.80   17.76  0.1138E+02  0.8308E+01
   31    3   15    1.80   18.78  0.1203E+02  0.1162E+02
   30   22   11    1.80   18.23  0.1168E+02  0.8332E+01
   43    1    2    1.80   19.85  0.1272E+02  0.8511E+01
   37    1   11    1.81   14.14  0.9058E+01  0.8966E+01
   35   23    3    1.84   24.16  0.1548E+02  0.8889E+01
   40    1    6    1.86   20.48  0.1312E+02  0.1016E+02
   30    1   13    1.96   14.03 -0.2283E+02  0.1374E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        74
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1074
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)     929
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      145


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
