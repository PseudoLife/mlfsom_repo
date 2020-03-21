
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.334    77.334    38.703  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=       8       8
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       8
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  1.000000 -0.000698  0.000704
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
 BEAM_DIVERGENCE_E.S.D.=   0.041
 REFLECTING_RANGE_E.S.D.=   0.088
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

    1485 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       8 ...       8
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1485 REFLECTIONS ACCEPTED



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
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    27        mC        500.0     159.4   38.7  109.4  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.7  109.4  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.4  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        998.1     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.1      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.1      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.1      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.4   86.5   86.5  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.7   86.5  247.6 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.5   86.5   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.2  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.5   86.5  177.2 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
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

  *   0.00        0    5.5      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.699   -23.430     0.022
 UNIT_CELL_B-AXIS=    23.430    73.699     0.009
 UNIT_CELL_C-AXIS=    -0.012    -0.001    38.703
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       330 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.21
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.314    77.314    38.690  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  8.8E-02 8.8E-02 4.4E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012934  0.012934  0.025846  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.681   -23.422     0.029
 COORDINATES OF UNIT CELL B-AXIS    23.422    73.681     0.010
 COORDINATES OF UNIT CELL C-AXIS    -0.015    -0.000    38.690
 CRYSTAL MOSAICITY (DEGREES)     0.088
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006341  0.002005
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1260.95   1242.10
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.23
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.31    77.31    38.69  90.000  90.000  90.000


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
  0.05  0.15       116        234.310       1.145
  0.15  0.25       179        247.329       1.085
  0.25  0.35       192        237.151       1.131
  0.35  0.45       146        158.824       1.689
  0.45  0.55       180        333.883       0.804
  0.55  0.65       196        260.205       1.031
  0.65  0.75       182        198.195       1.354
  0.75  0.85       170        255.161       1.052
  0.85  0.95       124        546.320       0.491
 
 -0.05  0.95      1485        268.327       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               234.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.032
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     7.1 XMAX=     7.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59666 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        500


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               234.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.029
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    30.6 XMAX=  2438.0 NXBIN=    3
 YMIN=     7.1 YMAX=  2518.9 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        500


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               234.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.021
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     7.1 XMAX=     7.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1234    1263
 DETECTOR_SURFACE_POSITION=    1643    1690
 DETECTOR_SURFACE_POSITION=     825    1690
 DETECTOR_SURFACE_POSITION=     825     836
 DETECTOR_SURFACE_POSITION=    1643     836
 DETECTOR_SURFACE_POSITION=    2160    1663
 DETECTOR_SURFACE_POSITION=    1618    2228
 DETECTOR_SURFACE_POSITION=     851    2228
 DETECTOR_SURFACE_POSITION=     309    1663
 DETECTOR_SURFACE_POSITION=     309     863
 DETECTOR_SURFACE_POSITION=     851     298
 DETECTOR_SURFACE_POSITION=    1618     298
 DETECTOR_SURFACE_POSITION=    2160     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        500


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
 1.388E+00  2.144E-03   18.33


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

   24.449   3.397    16.12   1.53      2.58      3.17      66      91      18
    3.397   2.414    12.71   1.49      4.74      3.75     130     170      12
    2.414   1.974     8.02   1.04      6.20      5.44     188     221       0
    1.974   1.711    -0.00   0.57   3516.93   4601.32     190     240       0
    1.711   1.531    -0.03   0.55    -99.90    -99.90     218     278       0
    1.531   1.398    -0.05   0.59    -99.90    -99.90     246     309       0
    1.398   1.295     0.01   0.57    -99.90    -99.90      98     143       2
   --------------------------------------------------------------------------
   24.449   1.295     3.70   0.81      5.58      5.98    1136    1452      32
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       8       8)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681   0.1116E+02  0.1380E+02         0.81           1 
   47    0    0     1.645  -0.9599E+01  0.1344E+02        -0.71           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -3.4% OF MEAN INTENSITY


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

     3.84          69      44      1285        3.4%       2.4%      3.2%       50   20.15      3.4%    99.9*     0    0.000       0
     2.73         105      67      2098        3.2%       4.4%      3.5%       76   17.60      6.2%    99.3*     0    0.000       0
     2.24         167     100      2665        3.8%       4.9%      4.1%      134   14.42      6.9%    98.8*     0    0.000       0
     1.94         166      94      3101        3.0%       7.4%      7.2%      144    7.09     10.4%    98.6*     0    0.000       0
     1.73         185     112      3479        3.2%    1075.7%   1303.4%      146    0.02   1521.3%    -8.3      0    0.000       0
     1.58         214     127      3850        3.3%     -99.9%    -99.9%      174    0.00    -99.9%   -16.2      0    0.000       0
     1.47         241     149      4144        3.6%    1940.6%   2532.7%      184   -0.06   2744.4%   -15.9      0    0.000       0
     1.37         220     131      4427        3.0%     -99.9%    -99.9%      178    0.00    -99.9%     9.8      0    0.000       0
     1.29          85      60      4743        1.3%     -99.9%    -99.9%       50    0.00    -99.9%     2.3      0    0.000       0
    total        1452     884     29792        3.0%       5.6%      6.0%     1136    4.70      7.9%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1485
 NUMBER OF REJECTED MISFITS                                32
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1452
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    884



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

   WILSON LINE (using all data) : A=   7.995 B=  31.843 CORRELATION=  0.82
      #      RES      SS        <I>       log(<I>)       BO
      19     7.377   0.005  2.9026E+03       7.973       2.4
      21     4.621   0.012  2.9643E+03       7.994       0.0
      30     3.495   0.020  1.3568E+03       7.213      19.1
      31     2.995   0.028  6.9156E+02       6.539      26.1
      41     2.645   0.036  7.0223E+02       6.554      20.2
      54     2.393   0.044  4.3274E+02       6.070      22.1
      43     2.204   0.051  2.9661E+02       5.692      22.4
      45     2.055   0.059  1.2856E+02       4.856      26.5
      47     1.929   0.067  1.7518E+00       0.561      55.4
      61     1.828   0.075  3.0765E-01      -1.179      61.3
      52     1.739   0.083 -1.9635E+00    -999.999   -1000.0
      55     1.664   0.090 -2.8560E-01    -999.999   -1000.0
      61     1.595   0.098 -3.2983E-01    -999.999   -1000.0
      76     1.535   0.106 -7.2222E-01    -999.999   -1000.0
      79     1.480   0.114 -3.7525E-02    -999.999   -1000.0
      58     1.431   0.122 -2.1810E+00    -999.999   -1000.0
      53     1.390   0.129  4.5990E-01      -0.777      33.9
      40     1.351   0.137 -3.2637E-01    -999.999   -1000.0
      18     1.314   0.145  5.5880E+00       1.721      21.7


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.377        1.068        0.529        0.206
       3     4.621        1.183        0.695        0.304
       3     3.495        0.360        0.092        0.017
       4     2.995        0.969        0.406        0.129
       3     2.645        0.247        0.058        0.010
       5     2.393        0.768        0.349        0.119
       3     2.204        1.231        0.486        0.140
       8     2.055        0.267        0.087        0.022
       7     1.929       16.003       35.414       61.170
      10     1.664       18.060       20.795       48.378
       3     1.595       31.282       65.889      103.294
       3     1.535       11.367       20.511       29.325
      88   overall       18.917      -52.255      571.323


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      14     7.377        1.421        2.431        3.403
      18     4.621        1.082        1.198        1.132
      27     3.495        0.847        0.640        0.444
      27     2.995        0.870        0.773        0.625
      38     2.645        1.151        1.251        1.227
      49     2.393        0.924        0.763        0.536
      40     2.204        0.774        0.560        0.351
      37     2.055        1.454        1.948        2.405
      40     1.929       11.304        6.188       48.579
      50     1.828      336.274     1634.378    60575.127
      73     1.480       30.004       -7.322      517.082
      53     1.390      237.821      -42.459    33528.469
      18     1.314        3.396        7.107       14.227
     796   overall       37.928      100.650     6040.499

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.377  0.00  0.00  0.00  0.00  0.00  0.36  0.33  0.32  0.61  0.59
       3     4.621  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.495  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.995  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.645  1.34  0.96  0.80  0.70  0.64  0.59  0.56  1.06  1.01  0.98
       5     2.393  0.00  0.00  0.00  0.00  0.77  1.07  1.00  0.95  0.91  0.88
       3     2.204  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.055  2.01  1.81  1.50  1.32  1.20  1.11  1.05  0.99  0.95  1.10
       7     1.929  1.15  0.83  0.69  0.60  0.55  0.51  0.72  0.68  0.65  0.63
      10     1.664  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       3     1.595  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       3     1.535  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      88   overall  0.87  0.66  0.55  0.48  0.48  0.51  0.49  0.49  0.48  0.50


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      14     7.377  2.25  2.76  1.93  2.17  2.00  1.74  1.56  1.43  1.32  1.24
      18     4.621  1.75  1.84  1.29  1.01  1.27  1.23  1.10  1.01  1.12  1.14
      27     3.495  0.78  0.82  0.57  0.67  0.66  0.74  0.88  1.01  0.94  0.88
      27     2.995  1.56  1.43  1.00  1.12  1.13  0.99  1.18  1.14  1.12  1.05
      38     2.645  0.83  1.16  1.22  1.20  1.00  0.93  1.05  1.05  1.06  1.04
      49     2.393  1.29  1.01  1.02  1.05  1.04  1.13  1.14  1.04  1.10  1.03
      40     2.204  1.58  0.97  0.96  0.83  0.95  0.89  0.89  1.00  1.01  1.03
      37     2.055  1.70  1.19  1.15  1.23  1.10  0.96  0.91  1.03  1.00  1.03
      40     1.929  4.73  2.62  1.83  1.44  1.21  1.05  0.94  0.86  0.93  0.87
      50     1.828  4.83  2.54  1.85  1.46  1.22  1.06  0.95  0.87  0.81  0.76
      73     1.480  5.76  3.02  2.11  1.66  1.39  1.21  1.09  1.00  0.92  0.87
      53     1.390  4.56  2.39  1.75  1.37  1.15  1.00  0.90  0.82  0.76  0.72
      18     1.314  4.67  2.45  1.71  1.52  1.41  1.23  1.10  1.01  1.03  1.05
     796   overall  1.91  1.17  0.90  0.78  0.71  0.65  0.62  0.60  0.59  0.57


 List of    147 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   49    3    8    1.50  667.00 -0.2503E+02  0.1812E+02 "alien"
   47    4   12    1.46  402.02 -0.1509E+02  0.9863E+01 "alien"
   39   31    8    1.48  318.43 -0.1195E+02  0.1026E+02 "alien"
   39    6   17    1.48  283.80 -0.1065E+02  0.8819E+01 "alien"
   26   23   20    1.46  273.96 -0.1028E+02  0.1325E+02 "alien"
   53    2    1    1.46  260.84 -0.9788E+01  0.1266E+02 "alien"
   52    3    5    1.46  255.59 -0.9591E+01  0.1107E+02 "alien"
   43    5   15    1.47  251.06 -0.9421E+01  0.1069E+02 "alien"
   33   29   14    1.48  224.55 -0.8426E+01  0.1027E+02 "alien"
   23   13   23    1.46  218.13 -0.8185E+01  0.9191E+01 "alien"
   26   10   22    1.49  197.50 -0.7411E+01  0.9987E+01 "alien"
   37    7   18    1.48  185.94 -0.6978E+01  0.1432E+02 "alien"
   51    2    3    1.50  180.54 -0.6775E+01  0.1033E+02 "alien"
   41   31    4    1.49  176.50 -0.6623E+01  0.1093E+02 "alien"
   41   31    2    1.50  169.96 -0.6378E+01  0.1053E+02 "alien"
   42   32    2    1.46  165.62 -0.6215E+01  0.1143E+02 "alien"
   42    5   15    1.49  149.45 -0.5608E+01  0.9473E+01 "alien"
   38    6   17    1.51  141.42 -0.5307E+01  0.8799E+01 "alien"
   32   29   15    1.47  138.99 -0.5216E+01  0.1020E+02 "alien"
   52    2    3    1.48  131.91 -0.4950E+01  0.1047E+02 "alien"
   28   27   17    1.50  131.45 -0.4933E+01  0.9143E+01 "alien"
   33    8   20    1.47  128.02 -0.4804E+01  0.1221E+02 "alien"
   38   31   10    1.46  127.27 -0.4776E+01  0.1052E+02 "alien"
   41   32    4    1.47   99.41 -0.3731E+01  0.1133E+02 "alien"
   32    8   20    1.49   85.86 -0.3222E+01  0.8535E+01 "alien"
   52    2    1    1.48   81.98 -0.3076E+01  0.9994E+01 "alien"
   35   24    4    1.79   72.08  0.2218E+02  0.1715E+02 "alien"
   41   32    5    1.46   71.86 -0.2697E+01  0.1075E+02 "alien"
   37   34   14    1.34   66.46 -0.2169E+02  0.1541E+02 "alien"
   39    4   15    1.57   65.85 -0.2172E+02  0.1272E+02 "alien"
   49    3    9    1.48   64.64 -0.2426E+01  0.1032E+02 "alien"
   37    5   16    1.57   60.41 -0.1992E+02  0.1426E+02 "alien"
   21   19   16    1.81   54.60  0.1680E+02  0.1065E+02 "alien"
   35   25    8    1.68   53.13 -0.1517E+02  0.1491E+02 "alien"
   31    2   14    1.85   50.71  0.1560E+02  0.8091E+01 "alien"
   52    4    9    1.40   50.58  0.2326E+02  0.1118E+02 "alien"
   39   33   10    1.41   49.96  0.2298E+02  0.1227E+02 "alien"
   42    1    3    1.82   49.65  0.1527E+02  0.9900E+01 "alien"
   37    3   14    1.66   47.39 -0.1353E+02  0.8411E+01 "alien"
   29   20   23    1.34   46.62 -0.1521E+02  0.1183E+02 "alien"
   31   28   15    1.50   45.16 -0.1695E+01  0.1546E+02 "alien"
   39   31    9    1.46   44.59 -0.1673E+01  0.1207E+02 "alien"
   36   32   14    1.39   43.37  0.1995E+02  0.1952E+02 "alien"
   35   33   15    1.36   41.44 -0.1352E+02  0.1572E+02 "alien"
   35   26    9    1.64   41.38 -0.1182E+02  0.9360E+01 "alien"
   43    9   19    1.33   41.34 -0.1349E+02  0.1541E+02 "alien"
   50    2    6    1.50   41.18 -0.1545E+01  0.9965E+01 "alien"
   45    2    9    1.59   40.89 -0.1349E+02  0.9835E+01 "alien"
   34   33   16    1.35   39.62 -0.1293E+02  0.1060E+02 "alien"
   31    3   15    1.79   39.20  0.1206E+02  0.7913E+01 "alien"
   21   19   19    1.63   38.18 -0.1091E+02  0.1081E+02 "alien"
   45    8   17    1.36   36.89 -0.1204E+02  0.1089E+02 "alien"
   32   25   12    1.64   36.83 -0.1052E+02  0.9294E+01 "alien"
   37   10   21    1.36   36.63 -0.1195E+02  0.1337E+02 "alien"
   34   25    9    1.69   36.24 -0.1035E+02  0.9257E+01 "alien"
   34   11   22    1.36   35.78 -0.1168E+02  0.8294E+01 "alien"
   40   31    5    1.50   35.52 -0.1333E+01  0.1211E+02 "alien"
   33   11   22    1.38   34.34  0.1579E+02  0.9235E+01 "alien"
   42    1    1    1.84   32.89  0.1012E+02  0.1035E+02 "alien"
   39   28    3    1.60   32.70 -0.1078E+02  0.1007E+02 "alien"
   52    5   11    1.36   31.50 -0.1028E+02  0.1764E+02 "alien"
   50    7   14    1.34   31.14 -0.1016E+02  0.1099E+02 "alien"
   41    8   18    1.40   29.72  0.1367E+02  0.1741E+02 "alien"
   27   27   18    1.47   29.60 -0.1073E+02  0.9828E+01 "alien"
   33   27   13    1.55   28.80 -0.2080E+02  0.1471E+02 "alien"
   35    5   17    1.58   28.66 -0.9452E+01  0.1315E+02 "alien"
   41   31    0    1.50   28.33 -0.9131E+01  0.1509E+02 "alien"
   33   23    8    1.79   28.00  0.8615E+01  0.9366E+01 "alien"
   36   28    9    1.58   27.05 -0.8922E+01  0.1770E+02 "alien"
   38   28    5    1.60   26.59 -0.8770E+01  0.9869E+01 "alien"
   26   20   13    1.85   26.07  0.8022E+01  0.8413E+01 "alien"
   22   14   22    1.51   25.21 -0.1821E+02  0.8801E+01 "alien"
   33   28   13    1.53   25.02 -0.1807E+02  0.1631E+02 "alien"
   39   33   11    1.39   24.83  0.1142E+02  0.1168E+02 "alien"
   28   13   24    1.36   24.74 -0.8075E+01  0.1282E+02 "alien"
   30    7   19    1.58   24.62 -0.8121E+01  0.1298E+02 "alien"
   52    2    0    1.49   24.53 -0.5626E+01  0.1413E+02 "alien"
   53    2    0    1.46   24.35 -0.5502E+01  0.1627E+02 "alien"
   30   26   15    1.55   24.22 -0.1749E+02  0.1343E+02 "alien"
   49    7   15    1.34   24.16 -0.7884E+01  0.1217E+02 "alien"
   47    1    6    1.59   23.74 -0.7831E+01  0.1466E+02 "alien"
   27   23   15    1.66   23.21 -0.6629E+01  0.8585E+01 "alien"
   21    9   21    1.62   23.11 -0.7623E+01  0.7997E+01 "alien"
   43    3   12    1.57   22.72 -0.7494E+01  0.9169E+01 "alien"
   32   24   11    1.69   22.24 -0.6352E+01  0.9523E+01 "alien"
   34    5   17    1.60   21.75 -0.7175E+01  0.8023E+01 "alien"
   43    1    8    1.68   21.28 -0.6078E+01  0.1051E+02 "alien"
   37    9   20    1.40   21.13  0.9716E+01  0.8817E+01 "alien"
   25   20   14    1.82   21.11  0.6495E+01  0.8778E+01 "alien"
   41    1    4    1.85   19.76  0.6080E+01  0.9617E+01 "alien"
   21   15   21    1.57   19.74 -0.6512E+01  0.9451E+01 "alien"
   37   26    4    1.68   19.35 -0.5528E+01  0.9788E+01 "alien"
   21   17   20    1.60   19.16 -0.6319E+01  0.9126E+01 "alien"
   34   26   10    1.64   19.12 -0.5460E+01  0.9453E+01 "alien"
   34   27   11    1.59   19.05 -0.6284E+01  0.9947E+01 "alien"
   34   28   12    1.54   18.91 -0.1366E+02  0.1033E+02 "alien"
   24   10   22    1.51   18.75 -0.1354E+02  0.8526E+01 "alien"
   28   12   23    1.40   18.73  0.8613E+01  0.1222E+02 "alien"
   45    1    7    1.64   18.28 -0.5221E+01  0.9368E+01 "alien"
   53    5    9    1.38   18.17  0.8355E+01  0.1159E+02 "alien"
   30   25   14    1.61   17.87 -0.5893E+01  0.8917E+01 "alien"
   54    4    6    1.39   17.13  0.7880E+01  0.1056E+02 "alien"
   18   10   21    1.65   17.00 -0.4854E+01  0.6778E+01 "alien"
   46    6   15    1.40   16.76  0.7706E+01  0.1028E+02 "alien"
   34    2   13    1.80   16.72  0.5144E+01  0.1486E+02 "alien"
   53    4    8    1.39   16.67  0.7668E+01  0.1069E+02 "alien"
   43    8   18    1.37   16.65 -0.5435E+01  0.1032E+02 "alien"
   42    9   19    1.35   16.47 -0.5377E+01  0.9399E+01 "alien"
   50    6   13    1.36   16.19 -0.5285E+01  0.1058E+02 "alien"
   29   27   16    1.52   16.15 -0.1167E+02  0.1451E+02 "alien"
   22   19   15    1.85   15.91  0.4894E+01  0.1352E+02 "alien"
   34    4   16    1.65   15.84 -0.4525E+01  0.8440E+01 "alien"
   27   24   16    1.60   15.42 -0.5085E+01  0.1492E+02 "alien"
   32   32   17    1.37   15.27 -0.1439E+02  0.1122E+02 "alien"
   34   31   15    1.41   15.07  0.6929E+01  0.1071E+02 "alien"
   26    9   21    1.54   15.05 -0.1087E+02  0.9470E+01 "alien"
   41   31    1    1.50   14.85 -0.5574E+00  0.1060E+02 "alien"
   31   22    9    1.84   14.44  0.4442E+01  0.1025E+02 "alien"
   35   23    2    1.84   14.33  0.4410E+01  0.9898E+01 "alien"
   52    5   10    1.38   14.18  0.6521E+01  0.1109E+02 "alien"
   38   29    8    1.53   14.10 -0.1019E+02  0.1042E+02 "alien"
   35    6   18    1.53   13.86 -0.1001E+02  0.1427E+02 "alien"
   22   18   20    1.58   13.86 -0.4570E+01  0.6898E+01 "alien"
   31   29   19    1.36   13.76 -0.4490E+01  0.1024E+02 "alien"
   50    5   12    1.39   13.66  0.6281E+01  0.1127E+02 "alien"
   41    1    5    1.83   13.62  0.4191E+01  0.9847E+01 "alien"
   53    4    7    1.41   13.46  0.6192E+01  0.1069E+02 "alien"
   27   20   12    1.87   13.14  0.4042E+01  0.8450E+01 "alien"
   32   30   19    1.33   12.95 -0.4228E+01  0.1488E+02 "alien"
   31   27   15    1.52   12.80 -0.9246E+01  0.1561E+02 "alien"
   48    1    4    1.59   12.72 -0.4197E+01  0.9852E+01 "alien"
   19    9   22    1.59   12.68 -0.4181E+01  0.6571E+01 "alien"
   44    0    6    1.70   12.65 -0.8668E+01  0.9453E+01 "alien"
   27    9   21    1.52   12.19 -0.8804E+01  0.7717E+01 "alien"
   19    9   21    1.65   12.06 -0.3443E+01  0.6772E+01 "alien"
   39   28    2    1.60   12.01 -0.3963E+01  0.9967E+01 "alien"
   30   28   19    1.38   11.35  0.5218E+01  0.9792E+01 "alien"
   35    1   12    1.82   11.31  0.3480E+01  0.8755E+01 "alien"
   18   18   17    1.82   11.25  0.7559E+01  0.9466E+01 "alien"
   39    2   12    1.69   11.15 -0.3184E+01  0.8593E+01 "alien"
   36    0   10    1.88   11.08  0.7367E+01  0.9085E+01 "alien"
   44    4   13    1.51   10.75 -0.7767E+01  0.9565E+01 "alien"
   30   26   20    1.37   10.72  0.4929E+01  0.1579E+02 "alien"
   25   23   19    1.52   10.63 -0.7678E+01  0.9389E+01 "alien"
   46    0    3    1.67   10.56 -0.6300E+01  0.9377E+01 "alien"
   45    2   10    1.57   10.36 -0.3416E+01  0.1056E+02 "alien"
   15   10   19    1.84   10.33  0.3177E+01  0.7105E+01 "alien"


 List of    147 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   43    9   19    1.33   41.34 -0.1349E+02  0.1541E+02
   32   30   19    1.33   12.95 -0.4228E+01  0.1488E+02
   29   20   23    1.34   46.62 -0.1521E+02  0.1183E+02
   49    7   15    1.34   24.16 -0.7884E+01  0.1217E+02
   50    7   14    1.34   31.14 -0.1016E+02  0.1099E+02
   37   34   14    1.34   66.46 -0.2169E+02  0.1541E+02
   42    9   19    1.35   16.47 -0.5377E+01  0.9399E+01
   34   33   16    1.35   39.62 -0.1293E+02  0.1060E+02
   28   13   24    1.36   24.74 -0.8075E+01  0.1282E+02
   31   29   19    1.36   13.76 -0.4490E+01  0.1024E+02
   45    8   17    1.36   36.89 -0.1204E+02  0.1089E+02
   37   10   21    1.36   36.63 -0.1195E+02  0.1337E+02
   35   33   15    1.36   41.44 -0.1352E+02  0.1572E+02
   52    5   11    1.36   31.50 -0.1028E+02  0.1764E+02
   50    6   13    1.36   16.19 -0.5285E+01  0.1058E+02
   34   11   22    1.36   35.78 -0.1168E+02  0.8294E+01
   43    8   18    1.37   16.65 -0.5435E+01  0.1032E+02
   32   32   17    1.37   15.27 -0.1439E+02  0.1122E+02
   30   26   20    1.37   10.72  0.4929E+01  0.1579E+02
   53    5    9    1.38   18.17  0.8355E+01  0.1159E+02
   33   11   22    1.38   34.34  0.1579E+02  0.9235E+01
   52    5   10    1.38   14.18  0.6521E+01  0.1109E+02
   30   28   19    1.38   11.35  0.5218E+01  0.9792E+01
   36   32   14    1.39   43.37  0.1995E+02  0.1952E+02
   50    5   12    1.39   13.66  0.6281E+01  0.1127E+02
   39   33   11    1.39   24.83  0.1142E+02  0.1168E+02
   53    4    8    1.39   16.67  0.7668E+01  0.1069E+02
   54    4    6    1.39   17.13  0.7880E+01  0.1056E+02
   46    6   15    1.40   16.76  0.7706E+01  0.1028E+02
   37    9   20    1.40   21.13  0.9716E+01  0.8817E+01
   52    4    9    1.40   50.58  0.2326E+02  0.1118E+02
   28   12   23    1.40   18.73  0.8613E+01  0.1222E+02
   41    8   18    1.40   29.72  0.1367E+02  0.1741E+02
   53    4    7    1.41   13.46  0.6192E+01  0.1069E+02
   34   31   15    1.41   15.07  0.6929E+01  0.1071E+02
   39   33   10    1.41   49.96  0.2298E+02  0.1227E+02
   53    2    1    1.46  260.84 -0.9788E+01  0.1266E+02
   53    2    0    1.46   24.35 -0.5502E+01  0.1627E+02
   52    3    5    1.46  255.59 -0.9591E+01  0.1107E+02
   23   13   23    1.46  218.13 -0.8185E+01  0.9191E+01
   39   31    9    1.46   44.59 -0.1673E+01  0.1207E+02
   41   32    5    1.46   71.86 -0.2697E+01  0.1075E+02
   38   31   10    1.46  127.27 -0.4776E+01  0.1052E+02
   42   32    2    1.46  165.62 -0.6215E+01  0.1143E+02
   26   23   20    1.46  273.96 -0.1028E+02  0.1325E+02
   47    4   12    1.46  402.02 -0.1509E+02  0.9863E+01
   43    5   15    1.47  251.06 -0.9421E+01  0.1069E+02
   41   32    4    1.47   99.41 -0.3731E+01  0.1133E+02
   32   29   15    1.47  138.99 -0.5216E+01  0.1020E+02
   27   27   18    1.47   29.60 -0.1073E+02  0.9828E+01
   33    8   20    1.47  128.02 -0.4804E+01  0.1221E+02
   52    2    3    1.48  131.91 -0.4950E+01  0.1047E+02
   39   31    8    1.48  318.43 -0.1195E+02  0.1026E+02
   49    3    9    1.48   64.64 -0.2426E+01  0.1032E+02
   33   29   14    1.48  224.55 -0.8426E+01  0.1027E+02
   52    2    1    1.48   81.98 -0.3076E+01  0.9994E+01
   37    7   18    1.48  185.94 -0.6978E+01  0.1432E+02
   39    6   17    1.48  283.80 -0.1065E+02  0.8819E+01
   26   10   22    1.49  197.50 -0.7411E+01  0.9987E+01
   52    2    0    1.49   24.53 -0.5626E+01  0.1413E+02
   41   31    4    1.49  176.50 -0.6623E+01  0.1093E+02
   42    5   15    1.49  149.45 -0.5608E+01  0.9473E+01
   32    8   20    1.49   85.86 -0.3222E+01  0.8535E+01
   28   27   17    1.50  131.45 -0.4933E+01  0.9143E+01
   49    3    8    1.50  667.00 -0.2503E+02  0.1812E+02
   40   31    5    1.50   35.52 -0.1333E+01  0.1211E+02
   41   31    2    1.50  169.96 -0.6378E+01  0.1053E+02
   50    2    6    1.50   41.18 -0.1545E+01  0.9965E+01
   41   31    1    1.50   14.85 -0.5574E+00  0.1060E+02
   31   28   15    1.50   45.16 -0.1695E+01  0.1546E+02
   41   31    0    1.50   28.33 -0.9131E+01  0.1509E+02
   51    2    3    1.50  180.54 -0.6775E+01  0.1033E+02
   38    6   17    1.51  141.42 -0.5307E+01  0.8799E+01
   44    4   13    1.51   10.75 -0.7767E+01  0.9565E+01
   22   14   22    1.51   25.21 -0.1821E+02  0.8801E+01
   24   10   22    1.51   18.75 -0.1354E+02  0.8526E+01
   25   23   19    1.52   10.63 -0.7678E+01  0.9389E+01
   29   27   16    1.52   16.15 -0.1167E+02  0.1451E+02
   31   27   15    1.52   12.80 -0.9246E+01  0.1561E+02
   27    9   21    1.52   12.19 -0.8804E+01  0.7717E+01
   35    6   18    1.53   13.86 -0.1001E+02  0.1427E+02
   33   28   13    1.53   25.02 -0.1807E+02  0.1631E+02
   38   29    8    1.53   14.10 -0.1019E+02  0.1042E+02
   26    9   21    1.54   15.05 -0.1087E+02  0.9470E+01
   34   28   12    1.54   18.91 -0.1366E+02  0.1033E+02
   33   27   13    1.55   28.80 -0.2080E+02  0.1471E+02
   30   26   15    1.55   24.22 -0.1749E+02  0.1343E+02
   39    4   15    1.57   65.85 -0.2172E+02  0.1272E+02
   43    3   12    1.57   22.72 -0.7494E+01  0.9169E+01
   45    2   10    1.57   10.36 -0.3416E+01  0.1056E+02
   21   15   21    1.57   19.74 -0.6512E+01  0.9451E+01
   37    5   16    1.57   60.41 -0.1992E+02  0.1426E+02
   22   18   20    1.58   13.86 -0.4570E+01  0.6898E+01
   35    5   17    1.58   28.66 -0.9452E+01  0.1315E+02
   36   28    9    1.58   27.05 -0.8922E+01  0.1770E+02
   30    7   19    1.58   24.62 -0.8121E+01  0.1298E+02
   19    9   22    1.59   12.68 -0.4181E+01  0.6571E+01
   48    1    4    1.59   12.72 -0.4197E+01  0.9852E+01
   34   27   11    1.59   19.05 -0.6284E+01  0.9947E+01
   47    1    6    1.59   23.74 -0.7831E+01  0.1466E+02
   45    2    9    1.59   40.89 -0.1349E+02  0.9835E+01
   39   28    3    1.60   32.70 -0.1078E+02  0.1007E+02
   34    5   17    1.60   21.75 -0.7175E+01  0.8023E+01
   38   28    5    1.60   26.59 -0.8770E+01  0.9869E+01
   21   17   20    1.60   19.16 -0.6319E+01  0.9126E+01
   27   24   16    1.60   15.42 -0.5085E+01  0.1492E+02
   39   28    2    1.60   12.01 -0.3963E+01  0.9967E+01
   30   25   14    1.61   17.87 -0.5893E+01  0.8917E+01
   21    9   21    1.62   23.11 -0.7623E+01  0.7997E+01
   21   19   19    1.63   38.18 -0.1091E+02  0.1081E+02
   34   26   10    1.64   19.12 -0.5460E+01  0.9453E+01
   35   26    9    1.64   41.38 -0.1182E+02  0.9360E+01
   32   25   12    1.64   36.83 -0.1052E+02  0.9294E+01
   45    1    7    1.64   18.28 -0.5221E+01  0.9368E+01
   19    9   21    1.65   12.06 -0.3443E+01  0.6772E+01
   34    4   16    1.65   15.84 -0.4525E+01  0.8440E+01
   18   10   21    1.65   17.00 -0.4854E+01  0.6778E+01
   37    3   14    1.66   47.39 -0.1353E+02  0.8411E+01
   27   23   15    1.66   23.21 -0.6629E+01  0.8585E+01
   46    0    3    1.67   10.56 -0.6300E+01  0.9377E+01
   37   26    4    1.68   19.35 -0.5528E+01  0.9788E+01
   43    1    8    1.68   21.28 -0.6078E+01  0.1051E+02
   35   25    8    1.68   53.13 -0.1517E+02  0.1491E+02
   34   25    9    1.69   36.24 -0.1035E+02  0.9257E+01
   39    2   12    1.69   11.15 -0.3184E+01  0.8593E+01
   32   24   11    1.69   22.24 -0.6352E+01  0.9523E+01
   44    0    6    1.70   12.65 -0.8668E+01  0.9453E+01
   33   23    8    1.79   28.00  0.8615E+01  0.9366E+01
   31    3   15    1.79   39.20  0.1206E+02  0.7913E+01
   35   24    4    1.79   72.08  0.2218E+02  0.1715E+02
   34    2   13    1.80   16.72  0.5144E+01  0.1486E+02
   21   19   16    1.81   54.60  0.1680E+02  0.1065E+02
   25   20   14    1.82   21.11  0.6495E+01  0.8778E+01
   18   18   17    1.82   11.25  0.7559E+01  0.9466E+01
   35    1   12    1.82   11.31  0.3480E+01  0.8755E+01
   42    1    3    1.82   49.65  0.1527E+02  0.9900E+01
   41    1    5    1.83   13.62  0.4191E+01  0.9847E+01
   35   23    2    1.84   14.33  0.4410E+01  0.9898E+01
   42    1    1    1.84   32.89  0.1012E+02  0.1035E+02
   31   22    9    1.84   14.44  0.4442E+01  0.1025E+02
   15   10   19    1.84   10.33  0.3177E+01  0.7105E+01
   26   20   13    1.85   26.07  0.8022E+01  0.8413E+01
   31    2   14    1.85   50.71  0.1560E+02  0.8091E+01
   41    1    4    1.85   19.76  0.6080E+01  0.9617E+01
   22   19   15    1.85   15.91  0.4894E+01  0.1352E+02
   27   20   12    1.87   13.14  0.4042E+01  0.8450E+01
   36    0   10    1.88   11.08  0.7367E+01  0.9085E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        89
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1485
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1309
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      176


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
