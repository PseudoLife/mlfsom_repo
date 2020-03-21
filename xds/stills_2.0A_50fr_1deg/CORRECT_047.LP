
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.644    77.644    38.870  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      47      47
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      47
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001030  0.000635
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
 REFLECTING_RANGE_E.S.D.=   0.204
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

    1087 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      47 ...      47
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1087 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               34
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
 *  25        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  33        mP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  20        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0   -1  1  0  0 -1 -1  0  0  0  0  1  0
 *  32        oP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  21        tP          0.0      77.6   77.6   38.9  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     160.1   38.9   77.6  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     160.1   38.9   77.6  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.9  160.1   77.6  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.9  160.1   77.6  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.9  160.1   77.6  90.0  90.0  75.9    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.9  160.1   77.6  90.0  90.0  75.9    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     160.1   38.9  109.8  90.0 133.3  75.9   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.9  109.8  116.5  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.9  160.1  160.1  86.6 104.1 104.1    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.8  116.5   38.9  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        997.5     173.6   77.6   38.9  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        997.5      77.6  173.6   38.9  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        997.5      77.6  173.6   38.9  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        997.5      77.6   77.6   38.9  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.5  116.5  116.5  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.8   86.8  116.5 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.8   86.8  116.5 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.8   86.8   86.8  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.8   86.8  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.9   86.8  248.6 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
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

  *   0.00        0    4.6      12    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    74.024   -23.427     0.052
 UNIT_CELL_B-AXIS=    23.427    74.024     0.023
 UNIT_CELL_C-AXIS=    -0.028    -0.003    38.871
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM        95 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.33
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.665    77.665    38.884  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.0E-01 2.0E-01 1.1E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012876  0.012876  0.025717  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    74.039   -23.455     0.004
 COORDINATES OF UNIT CELL B-AXIS    23.455    74.039     0.006
 COORDINATES OF UNIT CELL C-AXIS    -0.003    -0.002    38.884
 CRYSTAL MOSAICITY (DEGREES)     0.204
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006621  0.001687
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.88   1232.91
 DETECTOR ORIGIN (PIXELS) AT                     1261.43   1242.55
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.44
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.66    77.66    38.88  90.000  90.000  90.000


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
  0.05  0.15         7        470.980       0.189
  0.15  0.25       119         51.956       1.712
  0.25  0.35       161         83.098       1.071
  0.35  0.45       161         40.949       2.172
  0.45  0.55       161        168.557       0.528
  0.55  0.65       177         43.163       2.061
  0.65  0.75       177        124.194       0.716
  0.75  0.85       115         79.207       1.123
  0.85  0.95         9        153.135       0.581
 
 -0.05  0.95      1087         88.959       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        5
 DEGREES OF FREEDOM OF CHI^2 FIT               138.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.873
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    46.1 XMAX=    46.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59614 NYBIN=    5
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        282


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT               138.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.871
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    87.4 XMAX=  2450.4 NXBIN=    2
 YMIN=     6.2 YMAX=  2516.8 NYBIN=    2
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        282


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               131.2
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.874
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    46.1 XMAX=    46.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1269    1262
 DETECTOR_SURFACE_POSITION=    1670    1688
 DETECTOR_SURFACE_POSITION=     868    1688
 DETECTOR_SURFACE_POSITION=     868     835
 DETECTOR_SURFACE_POSITION=    1670     835
 DETECTOR_SURFACE_POSITION=    2177    1661
 DETECTOR_SURFACE_POSITION=    1645    2226
 DETECTOR_SURFACE_POSITION=     893    2226
 DETECTOR_SURFACE_POSITION=     361    1661
 DETECTOR_SURFACE_POSITION=     361     862
 DETECTOR_SURFACE_POSITION=     893     297
 DETECTOR_SURFACE_POSITION=    1645     297
 DETECTOR_SURFACE_POSITION=    2177     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        282


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
 8.399E-01  2.666E-03   21.13


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

   24.560   2.880    12.86   1.67      3.06      3.08      84     108       6
    2.880   2.044     1.84   0.80     23.79     27.70     182     222       2
    2.044   1.671    -0.06   0.65    -99.90    -99.90     218     258       2
    1.671   1.448    -0.05   0.79    -99.90    -99.90     232     290       2
    1.448   1.295     0.02   0.97   2519.95   2445.64     144     196       0
   --------------------------------------------------------------------------
   24.560   1.295     1.65   0.87     10.61     11.46     860    1074      12
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      47      47)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.688  -0.1353E+02  0.1426E+02        -0.95           1 
   47    0    0     1.652  -0.1396E+02  0.1290E+02        -1.08           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS   -14.7% OF MEAN INTENSITY


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

     3.84          59      38      1297        2.9%       2.6%      2.8%       42   19.92      3.7%    99.8*     0    0.000       0
     2.73          70      41      2129        1.9%       6.8%      5.9%       58   10.12      9.6%    98.9*     0    0.000       0
     2.24         131      77      2695        2.9%      24.4%     23.8%      108    2.70     34.6%    87.0*     0    0.000       0
     1.94         113      67      3136        2.1%     216.1%    316.3%       92    0.20    305.6%    49.0*     0    0.000       0
     1.74         150      86      3542        2.4%     -99.9%    -99.9%      128    0.00    -99.9%     4.2      0    0.000       0
     1.59         154      91      3877        2.3%     -99.9%    -99.9%      126    0.00    -99.9%     7.7      0    0.000       0
     1.47         167     100      4203        2.4%     -99.9%    -99.9%      134    0.00    -99.9%   -24.8      0    0.000       0
     1.37         157      94      4504        2.1%     673.8%    676.8%      126    0.18    952.9%     4.6      0    0.000       0
     1.30          71      49      4748        1.0%     -99.9%    -99.9%       44    0.00    -99.9%   -34.9      0    0.000       0
    total        1072     643     30131        2.1%      10.6%     11.4%      858    2.14     15.0%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1087
 NUMBER OF REJECTED MISFITS                                12
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1074
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    644



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

   WILSON LINE (using all data) : A=   6.287 B=  27.465 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      14     7.333   0.005  1.7608E+03       7.474    -127.6
      20     4.747   0.011  1.2684E+03       7.146     -38.7
      16     3.591   0.019  4.1423E+02       6.026       6.7
      23     2.981   0.028  8.0013E+01       4.382      33.9
      31     2.634   0.036  6.1135E+01       4.113      30.2
      40     2.411   0.043  1.9325E+01       2.961      38.7
      33     2.211   0.051  8.0612E+00       2.087      41.1
      33     2.055   0.059 -3.4750E-02    -999.999   -1000.0
      36     1.929   0.067  1.1593E+00       0.148      45.7
      39     1.826   0.075 -4.4999E-01    -999.999   -1000.0
      41     1.739   0.083 -2.4248E-01    -999.999   -1000.0
      46     1.665   0.090 -9.1255E-01    -999.999   -1000.0
      46     1.596   0.098  7.4575E-01      -0.293      33.5
      49     1.535   0.106 -1.2534E+00    -999.999   -1000.0
      47     1.482   0.114 -9.5446E-01    -999.999   -1000.0
      47     1.433   0.122  3.6032E+00       1.282      20.5
      38     1.392   0.129 -1.4564E-01    -999.999   -1000.0
      34     1.351   0.137 -2.5267E+00    -999.999   -1000.0
      11     1.311   0.145  5.7728E-01      -0.549      23.5


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       3     7.333        2.738        2.192        1.348
       3     4.747        1.481        0.978        0.480
       1     3.591        0.484        0.117        0.020
       3     2.981        3.231        2.131        1.034
       4     2.411        0.171        0.046        0.009
       2     2.211        1.093        0.513        0.186
       7     1.929       81.867      472.323     2530.958
       8     1.826       35.522      137.571      451.230
       4     1.596       56.763      137.430      428.022
       3     1.482       65.819      112.447      412.369
      73   overall       11.401       52.999      266.301


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.333        0.566        0.447        0.279
      17     4.747        1.296        1.748        1.905
      15     3.591        1.217        1.744        2.195
      20     2.981        0.492        0.317        0.175
      30     2.634        1.340        1.937        2.674
      36     2.411        1.155        1.102        0.946
      31     2.211        1.303        1.097        0.851
      28     2.055       38.277       89.622      794.841
      42     1.596       95.231      401.965     5588.726
      44     1.433        6.820        6.516       45.385
      37     1.392       84.081      223.969     3412.097
      11     1.311      219.204     -152.815    17925.547
     571   overall       13.780       29.407      784.884

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       3     7.333  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.49
       3     4.747  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       1     3.591  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.981  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       4     2.411  1.01  1.45  1.80  1.58  1.44  1.34  1.26  1.19  1.14  1.10
       2     2.211  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       7     1.929  1.73  1.24  1.03  0.91  0.82  0.76  0.72  0.68  0.65  0.63
       8     1.826  3.02  2.17  1.80  1.58  1.44  1.34  1.26  1.19  1.14  1.10
       4     1.596  1.01  0.72  0.60  0.53  0.48  0.89  0.84  0.80  0.76  0.73
       3     1.482  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      73   overall  0.83  0.63  0.56  0.49  0.45  0.44  0.41  0.39  0.38  0.38


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.333  2.87  3.01  2.46  2.21  1.85  1.61  1.63  1.49  1.38  1.29
      17     4.747  3.09  1.62  1.59  1.25  1.20  1.17  1.17  1.28  1.29  1.40
      15     3.591  0.70  0.37  0.51  0.61  0.85  0.74  1.32  1.21  1.24  1.27
      20     2.981  1.58  1.38  1.54  1.52  1.40  1.33  1.49  1.36  1.26  1.27
      30     2.634  0.70  0.92  1.03  1.01  0.93  1.11  1.13  1.09  1.07  1.05
      36     2.411  1.75  1.23  0.86  0.84  0.85  0.99  0.88  0.86  0.84  0.83
      31     2.211  2.71  1.60  1.24  1.08  0.98  0.86  0.83  0.88  0.82  0.77
      28     2.055  5.25  2.76  1.93  1.52  1.27  1.11  1.06  0.97  0.90  0.85
      42     1.596  6.00  3.15  2.20  1.73  1.45  1.27  1.14  1.04  0.96  0.94
      44     1.433  2.87  1.63  1.14  1.10  0.92  0.81  0.81  0.74  0.73  0.76
      37     1.392  4.83  2.53  1.77  1.39  1.17  1.02  0.91  0.83  0.77  0.73
      11     1.311  4.78  2.51  1.75  1.38  1.16  1.01  0.90  0.83  0.77  0.72
     571   overall  1.84  1.08  0.83  0.72  0.64  0.59  0.59  0.56  0.53  0.52


 List of     95 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   29   19    8    2.03  433.13 -0.1505E+02  0.9025E+01 "alien"
   37    2    4    2.05  374.97 -0.1303E+02  0.9089E+01 "alien"
   17   15   15    2.07  302.47 -0.1051E+02  0.1132E+02 "alien"
   36    2    5    2.08  235.16 -0.8172E+01  0.1344E+02 "alien"
   27    1   13    2.07  230.42 -0.8007E+01  0.7828E+01 "alien"
   26   18   11    2.02  151.14 -0.5252E+01  0.1403E+02 "alien"
   55    3    2    1.41  151.09 -0.2200E+02  0.1217E+02 "alien"
   36   32   14    1.39  148.20 -0.2158E+02  0.1681E+02 "alien"
   32   20    3    2.03  144.22 -0.5012E+01  0.9263E+01 "alien"
   26    1   13    2.11  135.21 -0.4699E+01  0.8084E+01 "alien"
   54    4    7    1.39  121.20 -0.1765E+02  0.1044E+02 "alien"
   32    1    9    2.12   97.46 -0.3387E+01  0.8773E+01 "alien"
   38    2    1    2.04   96.58 -0.3356E+01  0.9815E+01 "alien"
   21   20   17    1.74   96.13 -0.2331E+02  0.7949E+01 "alien"
   52    4    9    1.41   95.62 -0.1393E+02  0.1069E+02 "alien"
   28   23   21    1.40   89.10 -0.1298E+02  0.9409E+01 "alien"
   42    8   18    1.39   85.55 -0.1246E+02  0.1032E+02 "alien"
   36   10   21    1.38   75.48 -0.1099E+02  0.9588E+01 "alien"
   39   33   11    1.40   73.15 -0.1065E+02  0.1034E+02 "alien"
   45    7   16    1.40   69.93 -0.1018E+02  0.1566E+02 "alien"
   38   33   12    1.39   68.29 -0.9946E+01  0.9913E+01 "alien"
   38    2    2    2.03   58.93 -0.2048E+01  0.9275E+01 "alien"
   35   32   15    1.38   56.74 -0.8263E+01  0.1611E+02 "alien"
   50    5   12    1.40   56.04 -0.8161E+01  0.1136E+02 "alien"
   41    1    5    1.84   54.97 -0.2474E+02  0.1276E+02 "alien"
   50    6   13    1.37   49.92 -0.7270E+01  0.1006E+02 "alien"
   33   23    8    1.79   45.01 -0.2025E+02  0.8428E+01 "alien"
   30   19    6    2.07   40.41 -0.1404E+01  0.9204E+01 "alien"
   37    1   11    1.80   39.73 -0.1788E+02  0.8890E+01 "alien"
   36   27    8    1.63   39.00  0.2909E+02  0.1699E+02 "alien"
   55    3    1    1.41   38.08 -0.5545E+01  0.1109E+02 "alien"
   31    9   21    1.47   37.26 -0.3557E+02  0.1300E+02 "alien"
   32   24   11    1.70   37.04 -0.8982E+01  0.8879E+01 "alien"
   32    0   11    2.00   36.11 -0.2258E+02  0.9313E+01 "alien"
   51    5   11    1.39   35.40 -0.5156E+01  0.9962E+01 "alien"
   37   33   13    1.39   35.13 -0.5116E+01  0.1166E+02 "alien"
   40   11   21    1.32   34.86  0.2012E+02  0.1569E+02 "alien"
   20   19   18    1.71   33.79 -0.8194E+01  0.1166E+02 "alien"
   43    2   11    1.61   33.04  0.2464E+02  0.1546E+02 "alien"
   30    0   11    2.09   31.41 -0.1249E+02  0.8886E+01 "alien"
   39    4   15    1.57   30.33  0.2262E+02  0.1345E+02 "alien"
   26   22   15    1.71   29.17 -0.7073E+01  0.7988E+01 "alien"
   38    1   11    1.77   27.18 -0.6590E+01  0.8183E+01 "alien"
   30    8   20    1.54   26.86 -0.3367E+02  0.1215E+02 "alien"
   35   24    5    1.78   26.67 -0.6467E+01  0.9265E+01 "alien"
   42   10   20    1.32   25.67  0.1482E+02  0.1515E+02 "alien"
   43    8   18    1.37   25.27 -0.3680E+01  0.1056E+02 "alien"
   21   19   16    1.82   24.77 -0.1114E+02  0.1193E+02 "alien"
   33   22    6    1.87   23.94 -0.1077E+02  0.8441E+01 "alien"
   27   21   13    1.81   23.15 -0.1042E+02  0.7969E+01 "alien"
   34   23    6    1.82   22.80 -0.1026E+02  0.8456E+01 "alien"
   40    1   10    1.74   22.70 -0.5503E+01  0.8605E+01 "alien"
   37   27    6    1.64   20.95 -0.1912E+02  0.1450E+02 "alien"
   25   20   21    1.47   20.89 -0.1994E+02  0.1178E+02 "alien"
   39   28    2    1.61   20.64  0.1540E+02  0.9164E+01 "alien"
   44    1    0    1.76   20.58 -0.2181E+02  0.1461E+02 "alien"
   37   25    1    1.74   20.51 -0.4972E+01  0.8773E+01 "alien"
   49    8   16    1.32   19.42  0.1121E+02  0.1644E+02 "alien"
   25   22   20    1.49   19.40 -0.1851E+02  0.1347E+02 "alien"
   36   25    4    1.74   18.91 -0.4586E+01  0.1582E+02 "alien"
   50    3    8    1.48   17.82 -0.1700E+02  0.9380E+01 "alien"
   29   20   10    1.92   17.70  0.2052E+02  0.1387E+02 "alien"
   37   26    5    1.68   17.32 -0.1581E+02  0.8872E+01 "alien"
   35   23    1    1.85   17.12 -0.7703E+01  0.8996E+01 "alien"
   38   28    7    1.58   16.93  0.1263E+02  0.1583E+02 "alien"
   30    1   13    1.96   16.63  0.1928E+02  0.1452E+02 "alien"
   55    3    0    1.41   16.15 -0.7256E+01  0.1515E+02 "alien"
   41    1    4    1.86   15.23 -0.6854E+01  0.9106E+01 "alien"
   36    3   14    1.70   15.17 -0.1384E+02  0.1421E+02 "alien"
   39    2    1    1.99   14.10  0.1634E+02  0.9426E+01 "alien"
   34   22    0    1.92   13.96  0.4255E+02  0.1668E+02 "alien"
   35    5   17    1.58   13.93  0.1039E+02  0.1387E+02 "alien"
   23   21   19    1.58   13.83  0.1032E+02  0.1162E+02 "alien"
   33    8   20    1.48   13.72 -0.1309E+02  0.1347E+02 "alien"
   37   25    0    1.74   13.55 -0.8390E+01  0.1241E+02 "alien"
   32   26   13    1.59   13.26  0.9892E+01  0.9070E+01 "alien"
   40   31    7    1.48   13.20 -0.1260E+02  0.1057E+02 "alien"
   51    3    8    1.45   13.10  0.4720E+02  0.1803E+02 "alien"
   19   16   19    1.71   12.31 -0.2985E+01  0.1000E+02 "alien"
   31    0   11    2.04   12.27 -0.9967E+00  0.8356E+01 "alien"
   38    7   18    1.47   12.01 -0.1147E+02  0.1380E+02 "alien"
   48    2    8    1.53   11.64 -0.1459E+02  0.1021E+02 "alien"
   36   26    7    1.67   11.60 -0.1059E+02  0.8583E+01 "alien"
   28   21   12    1.83   11.52 -0.5185E+01  0.8132E+01 "alien"
   31   25   13    1.63   11.40 -0.1041E+02  0.8783E+01 "alien"
   51    6   13    1.35   11.17 -0.2822E+02  0.1820E+02 "alien"
   26   23   16    1.65   11.02 -0.1006E+02  0.1251E+02 "alien"
   48    3   10    1.49   10.65 -0.1017E+02  0.1477E+02 "alien"
   49    6   14    1.37   10.60 -0.2677E+02  0.1894E+02 "alien"
   47    1    5    1.62   10.48  0.7815E+01  0.9896E+01 "alien"
   27   20   12    1.88   10.47  0.1213E+02  0.8245E+01 "alien"
   33   11   22    1.39   10.32 -0.1504E+01  0.1386E+02 "alien"
   51    2    4    1.50   10.28 -0.9814E+01  0.9326E+01 "alien"
   43    0    6    1.74   10.04 -0.4921E+01  0.8482E+01 "alien"
   34    4   16    1.66   10.04 -0.9159E+01  0.8722E+01 "alien"


 List of     95 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   49    8   16    1.32   19.42  0.1121E+02  0.1644E+02
   40   11   21    1.32   34.86  0.2012E+02  0.1569E+02
   42   10   20    1.32   25.67  0.1482E+02  0.1515E+02
   51    6   13    1.35   11.17 -0.2822E+02  0.1820E+02
   49    6   14    1.37   10.60 -0.2677E+02  0.1894E+02
   50    6   13    1.37   49.92 -0.7270E+01  0.1006E+02
   43    8   18    1.37   25.27 -0.3680E+01  0.1056E+02
   36   10   21    1.38   75.48 -0.1099E+02  0.9588E+01
   35   32   15    1.38   56.74 -0.8263E+01  0.1611E+02
   33   11   22    1.39   10.32 -0.1504E+01  0.1386E+02
   37   33   13    1.39   35.13 -0.5116E+01  0.1166E+02
   54    4    7    1.39  121.20 -0.1765E+02  0.1044E+02
   42    8   18    1.39   85.55 -0.1246E+02  0.1032E+02
   51    5   11    1.39   35.40 -0.5156E+01  0.9962E+01
   38   33   12    1.39   68.29 -0.9946E+01  0.9913E+01
   36   32   14    1.39  148.20 -0.2158E+02  0.1681E+02
   50    5   12    1.40   56.04 -0.8161E+01  0.1136E+02
   45    7   16    1.40   69.93 -0.1018E+02  0.1566E+02
   39   33   11    1.40   73.15 -0.1065E+02  0.1034E+02
   28   23   21    1.40   89.10 -0.1298E+02  0.9409E+01
   55    3    2    1.41  151.09 -0.2200E+02  0.1217E+02
   52    4    9    1.41   95.62 -0.1393E+02  0.1069E+02
   55    3    1    1.41   38.08 -0.5545E+01  0.1109E+02
   55    3    0    1.41   16.15 -0.7256E+01  0.1515E+02
   51    3    8    1.45   13.10  0.4720E+02  0.1803E+02
   31    9   21    1.47   37.26 -0.3557E+02  0.1300E+02
   38    7   18    1.47   12.01 -0.1147E+02  0.1380E+02
   25   20   21    1.47   20.89 -0.1994E+02  0.1178E+02
   50    3    8    1.48   17.82 -0.1700E+02  0.9380E+01
   40   31    7    1.48   13.20 -0.1260E+02  0.1057E+02
   33    8   20    1.48   13.72 -0.1309E+02  0.1347E+02
   48    3   10    1.49   10.65 -0.1017E+02  0.1477E+02
   25   22   20    1.49   19.40 -0.1851E+02  0.1347E+02
   51    2    4    1.50   10.28 -0.9814E+01  0.9326E+01
   48    2    8    1.53   11.64 -0.1459E+02  0.1021E+02
   30    8   20    1.54   26.86 -0.3367E+02  0.1215E+02
   39    4   15    1.57   30.33  0.2262E+02  0.1345E+02
   38   28    7    1.58   16.93  0.1263E+02  0.1583E+02
   23   21   19    1.58   13.83  0.1032E+02  0.1162E+02
   35    5   17    1.58   13.93  0.1039E+02  0.1387E+02
   32   26   13    1.59   13.26  0.9892E+01  0.9070E+01
   43    2   11    1.61   33.04  0.2464E+02  0.1546E+02
   39   28    2    1.61   20.64  0.1540E+02  0.9164E+01
   47    1    5    1.62   10.48  0.7815E+01  0.9896E+01
   36   27    8    1.63   39.00  0.2909E+02  0.1699E+02
   31   25   13    1.63   11.40 -0.1041E+02  0.8783E+01
   37   27    6    1.64   20.95 -0.1912E+02  0.1450E+02
   26   23   16    1.65   11.02 -0.1006E+02  0.1251E+02
   34    4   16    1.66   10.04 -0.9159E+01  0.8722E+01
   36   26    7    1.67   11.60 -0.1059E+02  0.8583E+01
   37   26    5    1.68   17.32 -0.1581E+02  0.8872E+01
   36    3   14    1.70   15.17 -0.1384E+02  0.1421E+02
   32   24   11    1.70   37.04 -0.8982E+01  0.8879E+01
   26   22   15    1.71   29.17 -0.7073E+01  0.7988E+01
   19   16   19    1.71   12.31 -0.2985E+01  0.1000E+02
   20   19   18    1.71   33.79 -0.8194E+01  0.1166E+02
   40    1   10    1.74   22.70 -0.5503E+01  0.8605E+01
   37   25    1    1.74   20.51 -0.4972E+01  0.8773E+01
   37   25    0    1.74   13.55 -0.8390E+01  0.1241E+02
   21   20   17    1.74   96.13 -0.2331E+02  0.7949E+01
   43    0    6    1.74   10.04 -0.4921E+01  0.8482E+01
   36   25    4    1.74   18.91 -0.4586E+01  0.1582E+02
   44    1    0    1.76   20.58 -0.2181E+02  0.1461E+02
   38    1   11    1.77   27.18 -0.6590E+01  0.8183E+01
   35   24    5    1.78   26.67 -0.6467E+01  0.9265E+01
   33   23    8    1.79   45.01 -0.2025E+02  0.8428E+01
   37    1   11    1.80   39.73 -0.1788E+02  0.8890E+01
   27   21   13    1.81   23.15 -0.1042E+02  0.7969E+01
   34   23    6    1.82   22.80 -0.1026E+02  0.8456E+01
   21   19   16    1.82   24.77 -0.1114E+02  0.1193E+02
   28   21   12    1.83   11.52 -0.5185E+01  0.8132E+01
   41    1    5    1.84   54.97 -0.2474E+02  0.1276E+02
   35   23    1    1.85   17.12 -0.7703E+01  0.8996E+01
   41    1    4    1.86   15.23 -0.6854E+01  0.9106E+01
   33   22    6    1.87   23.94 -0.1077E+02  0.8441E+01
   27   20   12    1.88   10.47  0.1213E+02  0.8245E+01
   34   22    0    1.92   13.96  0.4255E+02  0.1668E+02
   29   20   10    1.92   17.70  0.2052E+02  0.1387E+02
   30    1   13    1.96   16.63  0.1928E+02  0.1452E+02
   39    2    1    1.99   14.10  0.1634E+02  0.9426E+01
   32    0   11    2.00   36.11 -0.2258E+02  0.9313E+01
   26   18   11    2.02  151.14 -0.5252E+01  0.1403E+02
   38    2    2    2.03   58.93 -0.2048E+01  0.9275E+01
   32   20    3    2.03  144.22 -0.5012E+01  0.9263E+01
   29   19    8    2.03  433.13 -0.1505E+02  0.9025E+01
   38    2    1    2.04   96.58 -0.3356E+01  0.9815E+01
   31    0   11    2.04   12.27 -0.9967E+00  0.8356E+01
   37    2    4    2.05  374.97 -0.1303E+02  0.9089E+01
   17   15   15    2.07  302.47 -0.1051E+02  0.1132E+02
   30   19    6    2.07   40.41 -0.1404E+01  0.9204E+01
   27    1   13    2.07  230.42 -0.8007E+01  0.7828E+01
   36    2    5    2.08  235.16 -0.8172E+01  0.1344E+02
   30    0   11    2.09   31.41 -0.1249E+02  0.8886E+01
   26    1   13    2.11  135.21 -0.4699E+01  0.8084E+01
   32    1    9    2.12   97.46 -0.3387E+01  0.8773E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        57
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1087
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)     980
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      107


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
