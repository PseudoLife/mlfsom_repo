
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.568    77.568    38.820  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      35      35
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      35
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000974  0.000678
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
 BEAM_DIVERGENCE_E.S.D.=   0.047
 REFLECTING_RANGE_E.S.D.=   0.159
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

    1250 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      35 ...      35
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1250 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               41
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
 *  31        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0  1  0  0  1  0  0  0
 *  34        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  35        mP          0.0      77.6   38.8   77.6  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  32        oP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    28        mC        250.0      38.8  159.9   77.6  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.8  159.9   77.6  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    37        mC        250.0     159.9   38.8   77.6  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.9   38.8   77.6  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.8  159.9   77.6  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.8  159.9   77.6  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    27        mC        500.0     159.9   38.8  109.7  90.0 133.3  76.0    0  2  1  0  0  0  1  0  1 -1  0  0
    19        oI        500.0      38.8  109.7  116.4  90.0  70.5  90.0    0  0  1  0  1 -1  0  0  1  1  1  0
    26        oF        625.0      38.8  159.9  159.9  86.6 104.0 104.0    0  0 -1  0  0  2  1  0  2  0  1  0
    18        tI        625.0     109.7  116.4   38.8  70.5  90.0  90.0    1 -1  0  0 -1 -1 -1  0  0  0 -1  0
    30        mC        998.1      77.6  173.4   38.8  90.0  90.0  63.4    0  1  0  0 -2  1  0  0  0  0  1  0
    41        mC        998.1     173.4   77.6   38.8  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    40        oC        998.1      77.6  173.4   38.8  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.1      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.4  116.4  116.4  96.4  96.4 141.0    1 -1 -1  0 -1  1 -1  0  1  1  1  0
     2        hR        999.0      86.7   86.7  116.4 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.7   86.7  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.7   86.7  116.4 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.7   86.7   86.7  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.7   86.7  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.7   86.7  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.8   86.7  248.3 102.1  99.0 116.6    0  0 -1  0  0  1  1  0  3 -1  1  0
    10        mC        999.0      86.7   86.7   77.6  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.7   86.7   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.8   77.6  177.7  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.7   86.7  177.7 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.7   86.7   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.7   86.7   86.7  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.7  109.7   38.8  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.8   77.6  177.7 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.7  177.7   77.6 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    4.3      17    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.887   -23.610     0.040
 UNIT_CELL_B-AXIS=    23.610    73.887     0.013
 UNIT_CELL_C-AXIS=    -0.021    -0.000    38.821
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       154 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.28
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.490    77.490    38.778  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.1E-01 1.1E-01 4.3E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012905  0.012905  0.025788  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.847   -23.479     0.058
 COORDINATES OF UNIT CELL B-AXIS    23.479    73.847     0.008
 COORDINATES OF UNIT CELL C-AXIS    -0.029     0.005    38.778
 CRYSTAL MOSAICITY (DEGREES)     0.159
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006340  0.002184
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.91
 DETECTOR ORIGIN (PIXELS) AT                     1260.69   1242.12
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.10
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.49    77.49    38.78  90.000  90.000  90.000


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
  0.05  0.15        38        101.801       1.082
  0.15  0.25       162        113.530       0.970
  0.25  0.35       172        104.827       1.050
  0.35  0.45       166         63.373       1.738
  0.45  0.55       165        193.924       0.568
  0.55  0.65       178         75.403       1.460
  0.65  0.75       194        122.403       0.900
  0.75  0.85       145        105.373       1.045
  0.85  0.95        30         79.629       1.383
 
 -0.05  0.95      1250        110.115       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        8
 DEGREES OF FREEDOM OF CHI^2 FIT               204.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.836
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    34.1 XMAX=    34.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59396 NYBIN=    8
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        434


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               204.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.831
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    13.4 XMAX=  2451.4 NXBIN=    3
 YMIN=     6.5 YMAX=  2520.2 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        434


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               204.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.825
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    34.1 XMAX=    34.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1263
 DETECTOR_SURFACE_POSITION=    1646    1690
 DETECTOR_SURFACE_POSITION=     818    1690
 DETECTOR_SURFACE_POSITION=     818     836
 DETECTOR_SURFACE_POSITION=    1646     836
 DETECTOR_SURFACE_POSITION=    2169    1664
 DETECTOR_SURFACE_POSITION=    1621    2230
 DETECTOR_SURFACE_POSITION=     844    2230
 DETECTOR_SURFACE_POSITION=     295    1664
 DETECTOR_SURFACE_POSITION=     295     863
 DETECTOR_SURFACE_POSITION=     844     297
 DETECTOR_SURFACE_POSITION=    1621     297
 DETECTOR_SURFACE_POSITION=    2169     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        434


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
 7.979E-01  1.124E-03   33.39


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

   24.504   3.156    22.99   1.20      2.02      1.90      58      80      26
    3.156   2.241     6.39   1.29      9.02      7.50     170     202       4
    2.241   1.832     0.56   0.85     82.12     88.43     182     209       0
    1.832   1.588    -0.01   0.86    -99.90    -99.90     228     262       2
    1.588   1.421    -0.03   0.86   2332.84   2488.07     240     297       2
    1.421   1.298     0.12   1.03    307.75    286.97     112     163       2
   --------------------------------------------------------------------------
   24.504   1.298     2.68   0.97      8.80      8.73     990    1213      36
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      35      35)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.685   0.2494E+02  0.1321E+02         1.89           1 
   47    0    0     1.649   0.7659E+01  0.1250E+02         0.61           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     6.8% OF MEAN INTENSITY


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

     3.85          54      34      1285        2.6%       1.9%      1.8%       40   31.46      2.7%    99.9*     0    0.000       0
     2.74          76      46      2098        2.2%       5.0%      3.8%       60   16.51      7.0%    99.3*     0    0.000       0
     2.24         152      88      2665        3.3%      11.1%      9.7%      128    6.92     15.7%    97.2*     0    0.000       0
     1.94         134      75      3101        2.4%      54.4%     56.3%      118    1.21     76.9%    58.8*     0    0.000       0
     1.74         153      86      3479        2.5%     -99.9%    -99.9%      134    0.00    -99.9%    19.5      0    0.000       0
     1.59         184     105      3851        2.7%     -99.9%    -99.9%      158    0.00    -99.9%    -6.2      0    0.000       0
     1.47         189     115      4143        2.8%    2388.5%   2631.1%      148   -0.05   3377.8%    15.0      0    0.000       0
     1.38         186     110      4427        2.5%    1026.3%   1048.4%      152    0.04   1451.4%     9.4      0    0.000       0
     1.30          83      57      4739        1.2%     219.5%    185.7%       52    0.28    310.4%    66.0*     0    0.000       0
    total        1211     716     29788        2.4%       8.8%      8.7%      990    3.54     12.4%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1250
 NUMBER OF REJECTED MISFITS                                36
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1213
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    718



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

   WILSON LINE (using all data) : A=   6.558 B=  23.788 CORRELATION=  0.85
      #      RES      SS        <I>       log(<I>)       BO
      14     7.399   0.005  1.8966E+03       7.548    -108.4
      16     4.711   0.011  1.8237E+03       7.509     -42.2
      19     3.598   0.019  6.5928E+02       6.491       1.7
      25     2.988   0.028  1.4916E+02       5.005      27.7
      32     2.637   0.036  1.4218E+02       4.957      22.3
      49     2.403   0.043  6.0392E+01       4.101      28.4
      34     2.210   0.051  2.6272E+01       3.269      32.1
      37     2.058   0.059  8.8512E+00       2.181      37.1
      38     1.935   0.067 -1.9134E+00    -999.999   -1000.0
      46     1.830   0.075 -8.2921E-02    -999.999   -1000.0
      44     1.740   0.083 -4.4442E-01    -999.999   -1000.0
      43     1.669   0.090  1.4842E+00       0.395      34.3
      50     1.599   0.098  4.8344E-01      -0.727      37.2
      56     1.540   0.105 -1.4286E+00    -999.999   -1000.0
      65     1.483   0.114 -6.2645E-01    -999.999   -1000.0
      49     1.435   0.121 -1.4623E+00    -999.999   -1000.0
      43     1.391   0.129  3.5008E+00       1.253      20.5
      39     1.353   0.137  2.4837E+00       0.910      20.7
      16     1.318   0.144  6.8346E+00       1.922      16.1


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.399        1.912        1.386        0.786
       3     4.711        1.202        0.718        0.319
       2     3.598        0.338        0.068        0.010
       3     2.988        2.690        1.598        0.693
       5     2.403        0.713        0.331        0.119
       1     2.210        2.036        1.007        0.355
       6     2.058        0.250        0.029        0.011
      10     1.830       31.394       -6.324      145.807
      11     1.669       16.669       12.194       51.886
       3     1.599       46.666       70.121      313.789
       3     1.540       35.399       84.336      198.747
       1     1.353        0.374        0.079        0.012
      74   overall        4.741        4.865       20.532


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.399        0.553        0.411        0.240
      13     4.711        1.242        1.502        1.438
      17     3.598        1.021        1.010        0.902
      22     2.988        0.690        0.578        0.415
      31     2.637        1.175        1.383        1.534
      44     2.403        0.958        0.847        0.648
      33     2.210        0.770        0.529        0.312
      31     2.058        1.647        1.800        2.015
      32     1.669       21.944       64.223      357.898
      47     1.599      195.720     -102.535    17990.706
      43     1.391        4.312        5.199        9.851
      38     1.353       34.376      398.665     4399.880
      16     1.318        2.184        2.908        4.009
     641   overall       18.154       20.066     1598.898

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.399  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.73
       3     4.711  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.598  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.73
       3     2.988  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       5     2.403  0.00  0.58  0.48  0.42  0.38  0.71  1.00  0.95  0.91  0.88
       1     2.210  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.058  2.01  1.45  1.20  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      10     1.830  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
      11     1.669  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.72  0.69  0.67
       3     1.599  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.540  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       1     1.353  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      74   overall  0.98  0.74  0.62  0.57  0.52  0.51  0.50  0.47  0.47  0.51


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.399  3.15  3.31  2.31  2.12  1.78  1.55  1.59  1.45  1.35  1.27
      13     4.711  2.42  1.27  1.48  1.17  1.17  1.53  1.38  1.54  1.43  1.34
      17     3.598  1.24  0.97  0.91  0.71  0.90  0.91  1.17  1.18  1.09  1.02
      22     2.988  2.39  2.26  1.58  1.24  1.39  1.41  1.35  1.32  1.30  1.22
      31     2.637  0.68  1.07  1.00  0.98  0.98  0.93  1.03  1.00  1.09  1.12
      44     2.403  0.48  0.75  1.05  1.10  1.10  1.16  1.08  1.03  1.11  1.08
      33     2.210  1.27  0.84  0.70  0.74  0.77  0.81  0.84  0.88  0.92  1.05
      31     2.058  2.37  1.60  1.12  0.98  0.82  0.86  0.77  0.70  0.71  0.66
      32     1.669  4.93  2.59  1.93  1.71  1.51  1.32  1.18  1.08  1.00  0.94
      47     1.599  4.47  2.35  1.72  1.36  1.19  1.04  0.93  0.89  0.82  0.77
      43     1.391  3.91  2.05  1.44  1.20  1.00  0.88  0.83  0.84  0.78  0.74
      38     1.353  6.64  3.63  2.54  2.08  1.81  1.57  1.41  1.29  1.20  1.12
      16     1.318  3.94  2.41  1.69  1.33  1.11  1.11  1.12  1.02  0.95  0.89
     641   overall  1.79  1.12  0.87  0.75  0.69  0.66  0.63  0.61  0.59  0.57


 List of     85 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   34   22    4    1.88  189.68 -0.1573E+02  0.9383E+01 "alien"
   34   23    5    1.83  185.37 -0.1537E+02  0.9404E+01 "alien"
   43    1    1    1.80  127.43 -0.1057E+02  0.8147E+01 "alien"
   20   18   16    1.85  124.24 -0.1030E+02  0.9061E+01 "alien"
   41    1    4    1.85  117.81 -0.9769E+01  0.8647E+01 "alien"
   28   21   12    1.83   95.04 -0.7881E+01  0.7465E+01 "alien"
   36   24    1    1.79   91.90 -0.7620E+01  0.8214E+01 "alien"
   34   23    6    1.81   91.34 -0.7574E+01  0.8035E+01 "alien"
   30   22   11    1.79   84.75 -0.7028E+01  0.7790E+01 "alien"
   35   23    1    1.85   76.84 -0.6372E+01  0.8205E+01 "alien"
   33   23    8    1.79   72.67 -0.6026E+01  0.7632E+01 "alien"
   27   20   12    1.88   67.15 -0.5568E+01  0.7323E+01 "alien"
   34    1   12    1.86   56.08 -0.4651E+01  0.7445E+01 "alien"
   54    5    9    1.36   44.65  0.1109E+03  0.1146E+02 "alien"
   35    2   13    1.78   42.35 -0.1882E+02  0.8398E+01 "alien"
   45    2    9    1.60   40.87  0.1976E+02  0.1418E+02 "alien"
   36   27    8    1.62   40.13  0.1940E+02  0.9606E+01 "alien"
   38   28    6    1.59   38.47  0.1860E+02  0.1184E+02 "alien"
   40    1    6    1.86   36.21 -0.3003E+01  0.8719E+01 "alien"
   40   31    7    1.48   33.61 -0.2106E+02  0.9431E+01 "alien"
   36   25    4    1.74   33.07 -0.1470E+02  0.8988E+01 "alien"
   21   19   16    1.81   30.68 -0.2544E+01  0.9685E+01 "alien"
   52    2    4    1.47   29.52 -0.1849E+02  0.1642E+02 "alien"
   25   20   21    1.47   29.23 -0.1831E+02  0.1021E+02 "alien"
   51    2    5    1.49   27.78 -0.1740E+02  0.9384E+01 "alien"
   37   31   11    1.46   25.82 -0.1618E+02  0.1636E+02 "alien"
   18   18   17    1.83   25.73 -0.1450E+02  0.9973E+01 "alien"
   43    1    0    1.80   25.68 -0.1441E+02  0.1178E+02 "alien"
   42   32    1    1.47   25.26 -0.1583E+02  0.8970E+01 "alien"
   23   21   19    1.58   24.99  0.1208E+02  0.7529E+01 "alien"
   36   25    5    1.72   24.72 -0.1099E+02  0.8061E+01 "alien"
   34   29   13    1.50   24.28 -0.1521E+02  0.9894E+01 "alien"
   37   25    2    1.73   24.07 -0.1070E+02  0.8945E+01 "alien"
   38   31   10    1.46   22.47 -0.1408E+02  0.8835E+01 "alien"
   31   26   14    1.58   21.28  0.1029E+02  0.8372E+01 "alien"
   33   25   11    1.65   21.02  0.3120E+02  0.1340E+02 "alien"
   39   28    2    1.61   20.57  0.9943E+01  0.8292E+01 "alien"
   33   26   12    1.60   20.48  0.9901E+01  0.9375E+01 "alien"
   39   29    4    1.57   20.45  0.9886E+01  0.1101E+02 "alien"
   47    4   12    1.46   19.98 -0.1252E+02  0.8974E+01 "alien"
   39    1   11    1.73   19.89 -0.8839E+01  0.1351E+02 "alien"
   30   28   16    1.49   19.66 -0.1232E+02  0.1432E+02 "alien"
   37   30   10    1.50   18.81 -0.1178E+02  0.1526E+02 "alien"
   37   29   10    1.52   18.70 -0.2672E+02  0.1503E+02 "alien"
   28   24   15    1.63   18.53  0.8959E+01  0.7664E+01 "alien"
   35    5   17    1.58   18.18  0.8790E+01  0.1305E+02 "alien"
   40    0    8    1.80   17.70 -0.5092E+01  0.7890E+01 "alien"
   53    3    6    1.42   17.59 -0.2572E+02  0.1138E+02 "alien"
   48    3    9    1.51   17.42 -0.1091E+02  0.1487E+02 "alien"
   35   23    0    1.85   17.23 -0.4781E+01  0.1205E+02 "alien"
   25   18   12    1.98   16.64 -0.3184E+02  0.1145E+02 "alien"
   34    3   15    1.71   16.55 -0.7355E+01  0.1318E+02 "alien"
   40    1   10    1.73   16.15 -0.7177E+01  0.7831E+01 "alien"
   35   24    6    1.76   16.05 -0.7134E+01  0.8532E+01 "alien"
   46    1    6    1.63   15.86  0.7666E+01  0.1224E+02 "alien"
   32   27   14    1.54   14.78 -0.2112E+02  0.9648E+01 "alien"
   46    4   12    1.49   14.50 -0.9081E+01  0.9573E+01 "alien"
   43    2   11    1.60   14.43  0.6974E+01  0.1000E+02 "alien"
   37    2   13    1.71   14.16 -0.6295E+01  0.8864E+01 "alien"
   45    2   10    1.57   14.06  0.6799E+01  0.9799E+01 "alien"
   29   23   13    1.71   13.71 -0.6093E+01  0.7396E+01 "alien"
   42    0    7    1.75   13.68 -0.1567E+02  0.8011E+01 "alien"
   28   23   14    1.69   13.62  0.2021E+02  0.7630E+01 "alien"
   39    0    8    1.84   13.35 -0.2786E+01  0.1504E+02 "alien"
   40   32    7    1.46   13.25 -0.1938E+02  0.1597E+02 "alien"
   50    3    8    1.47   13.05 -0.8175E+01  0.8868E+01 "alien"
   23   13   23    1.46   12.91 -0.8085E+01  0.9355E+01 "alien"
   30   25   14    1.61   12.71  0.6145E+01  0.7751E+01 "alien"
   42    4   14    1.53   12.66 -0.1808E+02  0.8408E+01 "alien"
   32   30   19    1.34   12.27  0.3048E+02  0.1435E+02 "alien"
   27   24   20    1.44   12.25 -0.1791E+02  0.1206E+02 "alien"
   26   25   17    1.56   12.20 -0.1742E+02  0.1274E+02 "alien"
   37   27    7    1.62   12.18  0.5887E+01  0.8455E+01 "alien"
   37   25    0    1.74   11.60 -0.1153E+02  0.1130E+02 "alien"
   33    3   15    1.73   11.53 -0.5126E+01  0.7425E+01 "alien"
   41    3   13    1.59   11.38  0.5501E+01  0.8343E+01 "alien"
   41   31    0    1.51   11.25 -0.1540E+02  0.1324E+02 "alien"
   21   17   20    1.61   11.19  0.5408E+01  0.9005E+01 "alien"
   35   31   14    1.42   10.76 -0.1573E+02  0.9096E+01 "alien"
   49    5   12    1.41   10.41 -0.1522E+02  0.9527E+01 "alien"
   49    4   11    1.44   10.40 -0.1520E+02  0.9326E+01 "alien"
   26   10   22    1.49   10.31 -0.6459E+01  0.9803E+01 "alien"
   39   31    8    1.48   10.20 -0.6391E+01  0.9474E+01 "alien"
   40   30    3    1.54   10.17 -0.1452E+02  0.9867E+01 "alien"
   38    1   11    1.76   10.03 -0.4456E+01  0.7641E+01 "alien"


 List of     85 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   32   30   19    1.34   12.27  0.3048E+02  0.1435E+02
   54    5    9    1.36   44.65  0.1109E+03  0.1146E+02
   49    5   12    1.41   10.41 -0.1522E+02  0.9527E+01
   35   31   14    1.42   10.76 -0.1573E+02  0.9096E+01
   53    3    6    1.42   17.59 -0.2572E+02  0.1138E+02
   27   24   20    1.44   12.25 -0.1791E+02  0.1206E+02
   49    4   11    1.44   10.40 -0.1520E+02  0.9326E+01
   40   32    7    1.46   13.25 -0.1938E+02  0.1597E+02
   37   31   11    1.46   25.82 -0.1618E+02  0.1636E+02
   23   13   23    1.46   12.91 -0.8085E+01  0.9355E+01
   38   31   10    1.46   22.47 -0.1408E+02  0.8835E+01
   47    4   12    1.46   19.98 -0.1252E+02  0.8974E+01
   42   32    1    1.47   25.26 -0.1583E+02  0.8970E+01
   25   20   21    1.47   29.23 -0.1831E+02  0.1021E+02
   52    2    4    1.47   29.52 -0.1849E+02  0.1642E+02
   50    3    8    1.47   13.05 -0.8175E+01  0.8868E+01
   40   31    7    1.48   33.61 -0.2106E+02  0.9431E+01
   39   31    8    1.48   10.20 -0.6391E+01  0.9474E+01
   26   10   22    1.49   10.31 -0.6459E+01  0.9803E+01
   46    4   12    1.49   14.50 -0.9081E+01  0.9573E+01
   30   28   16    1.49   19.66 -0.1232E+02  0.1432E+02
   51    2    5    1.49   27.78 -0.1740E+02  0.9384E+01
   34   29   13    1.50   24.28 -0.1521E+02  0.9894E+01
   37   30   10    1.50   18.81 -0.1178E+02  0.1526E+02
   41   31    0    1.51   11.25 -0.1540E+02  0.1324E+02
   48    3    9    1.51   17.42 -0.1091E+02  0.1487E+02
   37   29   10    1.52   18.70 -0.2672E+02  0.1503E+02
   42    4   14    1.53   12.66 -0.1808E+02  0.8408E+01
   40   30    3    1.54   10.17 -0.1452E+02  0.9867E+01
   32   27   14    1.54   14.78 -0.2112E+02  0.9648E+01
   26   25   17    1.56   12.20 -0.1742E+02  0.1274E+02
   45    2   10    1.57   14.06  0.6799E+01  0.9799E+01
   39   29    4    1.57   20.45  0.9886E+01  0.1101E+02
   31   26   14    1.58   21.28  0.1029E+02  0.8372E+01
   23   21   19    1.58   24.99  0.1208E+02  0.7529E+01
   35    5   17    1.58   18.18  0.8790E+01  0.1305E+02
   38   28    6    1.59   38.47  0.1860E+02  0.1184E+02
   41    3   13    1.59   11.38  0.5501E+01  0.8343E+01
   45    2    9    1.60   40.87  0.1976E+02  0.1418E+02
   33   26   12    1.60   20.48  0.9901E+01  0.9375E+01
   43    2   11    1.60   14.43  0.6974E+01  0.1000E+02
   21   17   20    1.61   11.19  0.5408E+01  0.9005E+01
   39   28    2    1.61   20.57  0.9943E+01  0.8292E+01
   30   25   14    1.61   12.71  0.6145E+01  0.7751E+01
   37   27    7    1.62   12.18  0.5887E+01  0.8455E+01
   36   27    8    1.62   40.13  0.1940E+02  0.9606E+01
   46    1    6    1.63   15.86  0.7666E+01  0.1224E+02
   28   24   15    1.63   18.53  0.8959E+01  0.7664E+01
   33   25   11    1.65   21.02  0.3120E+02  0.1340E+02
   28   23   14    1.69   13.62  0.2021E+02  0.7630E+01
   34    3   15    1.71   16.55 -0.7355E+01  0.1318E+02
   37    2   13    1.71   14.16 -0.6295E+01  0.8864E+01
   29   23   13    1.71   13.71 -0.6093E+01  0.7396E+01
   36   25    5    1.72   24.72 -0.1099E+02  0.8061E+01
   37   25    2    1.73   24.07 -0.1070E+02  0.8945E+01
   39    1   11    1.73   19.89 -0.8839E+01  0.1351E+02
   40    1   10    1.73   16.15 -0.7177E+01  0.7831E+01
   33    3   15    1.73   11.53 -0.5126E+01  0.7425E+01
   37   25    0    1.74   11.60 -0.1153E+02  0.1130E+02
   36   25    4    1.74   33.07 -0.1470E+02  0.8988E+01
   42    0    7    1.75   13.68 -0.1567E+02  0.8011E+01
   35   24    6    1.76   16.05 -0.7134E+01  0.8532E+01
   38    1   11    1.76   10.03 -0.4456E+01  0.7641E+01
   35    2   13    1.78   42.35 -0.1882E+02  0.8398E+01
   36   24    1    1.79   91.90 -0.7620E+01  0.8214E+01
   33   23    8    1.79   72.67 -0.6026E+01  0.7632E+01
   30   22   11    1.79   84.75 -0.7028E+01  0.7790E+01
   40    0    8    1.80   17.70 -0.5092E+01  0.7890E+01
   43    1    1    1.80  127.43 -0.1057E+02  0.8147E+01
   43    1    0    1.80   25.68 -0.1441E+02  0.1178E+02
   34   23    6    1.81   91.34 -0.7574E+01  0.8035E+01
   21   19   16    1.81   30.68 -0.2544E+01  0.9685E+01
   18   18   17    1.83   25.73 -0.1450E+02  0.9973E+01
   28   21   12    1.83   95.04 -0.7881E+01  0.7465E+01
   34   23    5    1.83  185.37 -0.1537E+02  0.9404E+01
   39    0    8    1.84   13.35 -0.2786E+01  0.1504E+02
   35   23    1    1.85   76.84 -0.6372E+01  0.8205E+01
   35   23    0    1.85   17.23 -0.4781E+01  0.1205E+02
   20   18   16    1.85  124.24 -0.1030E+02  0.9061E+01
   41    1    4    1.85  117.81 -0.9769E+01  0.8647E+01
   40    1    6    1.86   36.21 -0.3003E+01  0.8719E+01
   34    1   12    1.86   56.08 -0.4651E+01  0.7445E+01
   27   20   12    1.88   67.15 -0.5568E+01  0.7323E+01
   34   22    4    1.88  189.68 -0.1573E+02  0.9383E+01
   25   18   12    1.98   16.64 -0.3184E+02  0.1145E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        42
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1250
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1143
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      107


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
