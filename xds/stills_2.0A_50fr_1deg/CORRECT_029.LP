
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.524    77.524    38.794  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      29      29
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      29
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000939  0.000689
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
 BEAM_DIVERGENCE_E.S.D.=   0.046
 REFLECTING_RANGE_E.S.D.=   0.137
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

    1323 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      29 ...      29
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1323 REFLECTIONS ACCEPTED



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

 *  44        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.5   38.8   77.5  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  25        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  33        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  32        oP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  23        oC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0   -1  1  0  0 -1 -1  0  0  0  0  1  0
 *  21        tP          0.0      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.8   38.8   77.5  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.8   38.8   77.5  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  159.8   77.5  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.8   77.5  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  159.8   77.5  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  159.8   77.5  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     159.8   38.8  109.6  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.6  116.3  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  159.8  159.8  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.6  116.3   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        998.4     173.4   77.5   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.4      77.5  173.4   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.4      77.5  173.4   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.4      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.3  116.3  116.3  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.7   86.7  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.6   86.7   86.7  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.7   86.7  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.7   86.7  109.6  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.7  248.2 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.7   86.7   77.5  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.5  177.6  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.7   86.7  177.6 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.7   86.7   86.7  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.6  109.6   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.5  177.6 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.7  177.6   77.5 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    4.1      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.878   -23.496     0.067
 UNIT_CELL_B-AXIS=    23.496    73.878     0.027
 UNIT_CELL_C-AXIS=    -0.036    -0.003    38.793
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       185 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.32
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.461    77.461    38.760  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.6E-01 1.6E-01 7.2E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012910  0.012910  0.025800  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.821   -23.466     0.069
 COORDINATES OF UNIT CELL B-AXIS    23.466    73.821     0.022
 COORDINATES OF UNIT CELL C-AXIS    -0.036    -0.000    38.760
 CRYSTAL MOSAICITY (DEGREES)     0.137
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006370  0.002156
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1260.74   1242.15
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.16
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
  0.05  0.15        59         75.869       1.641
  0.15  0.25       162        136.959       0.909
  0.25  0.35       186        114.130       1.091
  0.35  0.45       159         76.640       1.624
  0.45  0.55       166        220.127       0.565
  0.55  0.65       185         77.190       1.613
  0.65  0.75       200        151.440       0.822
  0.75  0.85       146        121.759       1.022
  0.85  0.95        60         95.415       1.305
 
 -0.05  0.95      1323        124.481       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               216.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.971
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    28.1 XMAX=    28.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59444 NYBIN=    9
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        456


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               216.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.970
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    11.0 XMAX=  2453.9 NXBIN=    3
 YMIN=     8.0 YMAX=  2519.4 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        456


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               216.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.969
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    28.1 XMAX=    28.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1647    1690
 DETECTOR_SURFACE_POSITION=     818    1690
 DETECTOR_SURFACE_POSITION=     818     837
 DETECTOR_SURFACE_POSITION=    1647     837
 DETECTOR_SURFACE_POSITION=    2171    1664
 DETECTOR_SURFACE_POSITION=    1621    2229
 DETECTOR_SURFACE_POSITION=     844    2229
 DETECTOR_SURFACE_POSITION=     294    1664
 DETECTOR_SURFACE_POSITION=     294     864
 DETECTOR_SURFACE_POSITION=     844     298
 DETECTOR_SURFACE_POSITION=    1621     298
 DETECTOR_SURFACE_POSITION=    2171     864
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        456


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
 8.623E-01  1.201E-03   31.07


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

   24.495   3.155    22.71   1.25      2.25      2.03      68      93      20
    3.155   2.240     7.80   1.41      7.69      6.12     182     211       6
    2.240   1.832     1.04   0.87     42.11     46.63     198     225       2
    1.832   1.587     0.08   0.89    517.04    561.93     240     274       0
    1.587   1.420    -0.03   0.85    -99.90    -99.90     260     316       2
    1.420   1.297    -0.08   0.80    -99.90    -99.90     126     173       0
   --------------------------------------------------------------------------
   24.495   1.297     3.09   0.98      8.20      8.12    1074    1292      30
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      29      29)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.684   0.2834E+01  0.1216E+02         0.23           1 
   47    0    0     1.648   0.8008E+01  0.1256E+02         0.64           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     6.2% OF MEAN INTENSITY


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

     3.85          57      35      1285        2.7%       2.1%      1.9%       44   30.70      3.0%   100.0*     0    0.000       0
     2.74          88      53      2098        2.5%       4.2%      3.4%       70   19.04      5.9%    99.5*     0    0.000       0
     2.24         159      91      2665        3.4%       9.2%      7.3%      136    9.10     13.0%    98.5*     0    0.000       0
     1.94         146      82      3100        2.6%      29.3%     31.8%      128    2.05     41.5%    87.0*     0    0.000       0
     1.74         159      88      3481        2.5%    1228.6%   1433.5%      142    0.02   1737.5%     1.5      0    0.000       0
     1.59         194     110      3850        2.9%     365.8%    385.7%      168    0.13    517.4%    -4.2      0    0.000       0
     1.47         204     121      4144        2.9%     902.3%   1052.4%      166    0.05   1276.0%     8.4      0    0.000       0
     1.38         199     117      4427        2.6%     -99.9%    -99.9%      164    0.00    -99.9%   -10.0      0    0.000       0
     1.30          84      57      4737        1.2%     590.1%    724.9%       54    0.12    834.6%    25.8      0    0.000       0
    total        1290     754     29787        2.5%       8.2%      8.1%     1072    4.10     11.6%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1323
 NUMBER OF REJECTED MISFITS                                30
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1292
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    755



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

   WILSON LINE (using all data) : A=   7.361 B=  37.524 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      15     7.258   0.005  1.8481E+03       7.522     -16.9
      16     4.702   0.011  2.0505E+03       7.626     -11.7
      25     3.526   0.020  7.8461E+02       6.665      17.3
      26     3.001   0.028  2.0687E+02       5.332      36.6
      35     2.636   0.036  1.6562E+02       5.110      31.3
      49     2.402   0.043  9.0947E+01       4.510      32.9
      35     2.211   0.051  4.3814E+01       3.780      35.0
      42     2.057   0.059  1.6510E+01       2.804      38.6
      39     1.935   0.067  3.5336E-01      -1.040      62.9
      48     1.829   0.075 -7.8877E-02    -999.999   -1000.0
      43     1.740   0.083  6.5044E-01      -0.430      47.2
      48     1.668   0.090  2.3403E+00       0.850      36.2
      51     1.598   0.098  7.9223E-02      -2.535      50.6
      58     1.538   0.106  1.3001E+00       0.262      33.6
      70     1.482   0.114  1.3829E-01      -1.978      41.0
      52     1.434   0.122 -2.3538E+00    -999.999   -1000.0
      46     1.392   0.129 -1.1170E+00    -999.999   -1000.0
      39     1.353   0.137  2.4494E+00       0.896      23.7
      17     1.316   0.144 -1.1446E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.258        2.070        1.555        0.914
       3     4.702        1.254        0.766        0.347
       2     3.526        0.352        0.072        0.011
       3     3.001        2.176        1.202        0.495
       2     2.636        0.486        0.166        0.040
       5     2.402        0.855        0.475        0.204
       1     2.211        1.041        0.368        0.093
       6     2.057        0.623        0.252        0.081
      11     1.829       25.715       -8.212      208.299
       7     1.740      109.324       53.901     2338.610
      11     1.668        7.125        4.731        9.279
       6     1.482     2169.418    42632.793  1021511.132
      78   overall      288.459      948.488   120547.574


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.258        0.542        0.410        0.247
      13     4.702        1.234        1.493        1.428
      23     3.526        0.893        0.762        0.603
      23     3.001        0.752        0.611        0.409
      33     2.636        0.988        0.867        0.658
      44     2.402        0.969        0.852        0.633
      34     2.211        0.830        0.622        0.419
      36     2.057        1.306        1.559        1.872
      33     1.935      182.708     1378.189    24101.596
      36     1.740       77.050      247.435     3941.777
      37     1.668        7.130       10.853       26.174
      48     1.598     6886.778  -216641.695 32258349.989
      55     1.538       21.236       41.440      287.627
      38     1.353       16.279       23.311      151.613
     676   overall      751.016   -13544.187  2632206.059

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.258  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.37
       3     4.702  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.526  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.73
       3     3.001  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.636  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.402  0.00  0.00  0.00  0.85  0.77  0.71  1.00  0.95  0.91  0.88
       1     2.211  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.057  0.00  0.48  0.40  0.35  0.32  0.59  0.84  0.80  0.76  0.73
      11     1.829  2.20  1.58  1.31  1.15  1.05  0.97  0.91  0.87  0.83  0.80
       7     1.740  1.73  1.24  1.03  0.91  0.82  0.76  0.72  0.68  0.65  0.63
      11     1.668  1.46  1.31  1.09  0.96  0.87  0.97  1.22  1.16  1.11  1.07
       6     1.482  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      78   overall  0.88  0.70  0.58  0.57  0.52  0.53  0.58  0.55  0.55  0.56


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.258  1.91  3.01  2.10  1.65  1.85  1.61  1.63  1.49  1.38  1.29
      13     4.702  1.62  2.12  1.48  1.17  1.17  1.53  1.38  1.54  1.43  1.34
      23     3.526  0.91  0.72  0.67  0.53  0.66  0.87  1.12  1.11  1.10  1.03
      23     3.001  1.83  1.68  1.51  1.32  1.55  1.45  1.38  1.34  1.25  1.17
      33     2.636  0.96  1.17  0.94  1.10  1.00  0.94  0.90  0.99  0.97  1.01
      44     2.402  0.72  0.88  0.88  1.17  1.10  1.06  1.08  1.16  1.07  1.08
      34     2.211  1.55  0.97  0.79  0.62  0.60  0.59  0.70  0.75  0.84  0.93
      36     2.057  1.75  1.23  1.07  1.10  0.99  0.92  0.99  0.96  0.94  0.88
      33     1.935  4.78  2.51  1.75  1.38  1.23  1.14  1.02  0.94  0.87  0.81
      36     1.740  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      37     1.668  4.54  2.53  1.77  1.39  1.24  1.20  1.13  1.08  1.00  0.94
      48     1.598  5.04  2.64  1.85  1.45  1.22  1.06  0.95  0.87  0.81  0.76
      55     1.538  4.97  2.61  1.82  1.43  1.20  1.05  0.94  0.86  0.83  0.78
      38     1.353  4.98  2.76  1.93  1.52  1.27  1.11  0.99  0.96  0.89  0.83
     676   overall  2.22  1.36  1.01  0.86  0.78  0.73  0.71  0.69  0.66  0.63


 List of    113 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   31   22    9    1.84  442.57 -0.3491E+02  0.1461E+02 "alien"
   33   22    6    1.87  273.75 -0.2159E+02  0.8259E+01 "alien"
   35   23    3    1.83  242.82 -0.1915E+02  0.8509E+01 "alien"
   30   28   16    1.49  221.62  0.3065E+02  0.1473E+02 "alien"
   33    5   17    1.63  214.32  0.1698E+02  0.7386E+01 "alien"
   33   23    8    1.79  200.40 -0.1581E+02  0.8132E+01 "alien"
   46    2    9    1.57  190.93  0.1513E+02  0.8654E+01 "alien"
   36   27    9    1.60  149.82  0.1187E+02  0.9588E+01 "alien"
   27   25   19    1.46  143.87  0.1990E+02  0.1343E+02 "alien"
   37   27    7    1.62  139.65  0.1106E+02  0.9044E+01 "alien"
   45    1    8    1.62  138.34  0.1096E+02  0.9793E+01 "alien"
   43    2   11    1.60  137.48  0.1089E+02  0.9806E+01 "alien"
   42    1    3    1.83  121.43 -0.9578E+01  0.8815E+01 "alien"
   26   23   20    1.46  120.71  0.1669E+02  0.1321E+02 "alien"
   21   15   21    1.57  118.10  0.9357E+01  0.9846E+01 "alien"
   31   25   13    1.63  108.19  0.8571E+01  0.8103E+01 "alien"
   28   24   15    1.63  107.89  0.8547E+01  0.8124E+01 "alien"
   40    1    6    1.85  106.70 -0.8416E+01  0.8533E+01 "alien"
   25   24   17    1.60  105.56  0.8363E+01  0.1208E+02 "alien"
   50    2    6    1.51  103.88  0.1437E+02  0.9338E+01 "alien"
   21   19   16    1.81  103.73 -0.8182E+01  0.9584E+01 "alien"
   38   30    9    1.50  102.45  0.1417E+02  0.9115E+01 "alien"
   45    2   10    1.57  100.65  0.7974E+01  0.9840E+01 "alien"
   39   31    8    1.48   95.99  0.1327E+02  0.9636E+01 "alien"
   34    5   17    1.60   90.77  0.7191E+01  0.7447E+01 "alien"
   29    9   21    1.50   90.66  0.1254E+02  0.6903E+01 "alien"
   20   18   16    1.85   90.24 -0.7118E+01  0.9353E+01 "alien"
   38   28    6    1.59   90.11  0.7139E+01  0.1253E+02 "alien"
   23   13   23    1.46   89.73  0.1241E+02  0.9305E+01 "alien"
   36   27    8    1.62   85.69  0.6788E+01  0.9274E+01 "alien"
   25   10   22    1.50   83.76  0.1158E+02  0.9198E+01 "alien"
   39   28    4    1.59   83.05  0.6580E+01  0.1091E+02 "alien"
   44    5   14    1.48   76.25  0.1055E+02  0.1628E+02 "alien"
   39    6   17    1.49   75.70  0.1047E+02  0.8285E+01 "alien"
   39   29    4    1.57   73.64  0.5834E+01  0.1065E+02 "alien"
   27   21   13    1.80   70.48 -0.5560E+01  0.7269E+01 "alien"
   35   27   10    1.60   64.18  0.5084E+01  0.8600E+01 "alien"
   38    7   18    1.47   61.98  0.8571E+01  0.8761E+01 "alien"
   42    5   15    1.49   61.57  0.8514E+01  0.8590E+01 "alien"
   38    1    6    1.94   60.77  0.2147E+02  0.1030E+02 "alien"
   36   30   12    1.47   60.25  0.8332E+01  0.9089E+01 "alien"
   41   31    4    1.49   58.39  0.8075E+01  0.9375E+01 "alien"
   35   30   13    1.46   56.67  0.7837E+01  0.9012E+01 "alien"
   17   15   18    1.82   52.57 -0.4147E+01  0.5863E+01 "alien"
   48    1    3    1.60   52.11  0.4128E+01  0.1044E+02 "alien"
   24   23   17    1.63   51.42  0.4074E+01  0.8495E+01 "alien"
   48    1    5    1.58   48.18  0.3817E+01  0.8669E+01 "alien"
   42   32    1    1.47   45.52  0.6295E+01  0.9665E+01 "alien"
   29   21   11    1.84   43.98 -0.3469E+01  0.7802E+01 "alien"
   39   28    2    1.61   40.41  0.3202E+01  0.8834E+01 "alien"
   34    3   15    1.71   40.05  0.2605E+02  0.1322E+02 "alien"
   52    2    1    1.49   38.71  0.5353E+01  0.9198E+01 "alien"
   26   19   12    1.93   38.49  0.1360E+02  0.7914E+01 "alien"
   23   21   19    1.58   38.06  0.3015E+01  0.7028E+01 "alien"
   42    1    4    1.81   34.19 -0.2697E+01  0.8784E+01 "alien"
   50    8   15    1.32   32.77 -0.3751E+02  0.1804E+02 "alien"
   32    8   20    1.49   30.15  0.4169E+01  0.7916E+01 "alien"
   35   23    0    1.85   28.43 -0.1946E+02  0.1199E+02 "alien"
   51    2    3    1.51   28.34  0.3919E+01  0.9673E+01 "alien"
   29    1   13    1.99   28.10  0.9929E+01  0.7937E+01 "alien"
   46    1    7    1.61   26.99  0.2138E+01  0.8651E+01 "alien"
   38    6   17    1.51   25.79  0.3567E+01  0.9057E+01 "alien"
   27   27   18    1.48   24.79  0.2145E+02  0.1389E+02 "alien"
   34   22    3    1.89   24.63  0.8705E+01  0.8550E+01 "alien"
   41   31    1    1.51   24.52  0.3390E+01  0.9452E+01 "alien"
   37    1   11    1.80   23.52 -0.1855E+01  0.7926E+01 "alien"
   40    1    5    1.88   23.17 -0.1828E+01  0.1156E+02 "alien"
   39    2    2    1.97   22.60  0.7987E+01  0.8954E+01 "alien"
   39   28    1    1.61   21.40  0.1695E+01  0.8733E+01 "alien"
   36   29   11    1.51   20.53  0.2668E+02  0.1392E+02 "alien"
   43    5   15    1.47   19.81  0.2740E+01  0.9699E+01 "alien"
   37   31   11    1.46   19.55  0.2704E+01  0.1538E+02 "alien"
   41   31    0    1.51   19.47  0.1076E+02  0.1369E+02 "alien"
   32   28   20    1.33   18.40 -0.2106E+02  0.1442E+02 "alien"
   41   31    3    1.50   17.18  0.2376E+01  0.9013E+01 "alien"
   29   26   20    1.39   17.15 -0.1915E+02  0.1570E+02 "alien"
   28   20   11    1.90   16.89  0.5967E+01  0.8020E+01 "alien"
   39   31    9    1.46   16.75  0.2317E+01  0.1053E+02 "alien"
   37   33   13    1.38   16.74 -0.1870E+02  0.9703E+01 "alien"
   45   10   19    1.30   16.65 -0.1905E+02  0.1299E+02 "alien"
   45    7   16    1.39   16.55 -0.1849E+02  0.1069E+02 "alien"
   17   17   17    1.86   16.19 -0.3949E+01  0.6630E+01 "alien"
   25   19   13    1.90   16.14  0.5704E+01  0.7953E+01 "alien"
   32   30   19    1.34   15.73  0.3852E+02  0.1496E+02 "alien"
   44    1    1    1.76   15.68  0.1020E+02  0.9585E+01 "alien"
   46    6   15    1.40   15.60 -0.1743E+02  0.9285E+01 "alien"
   36    0   10    1.88   15.42 -0.3553E+01  0.8114E+01 "alien"
   33   24    9    1.74   14.98  0.9742E+01  0.8446E+01 "alien"
   39    2    1    1.98   14.54  0.5140E+01  0.8895E+01 "alien"
   21   18   15    1.90   14.11  0.4985E+01  0.7951E+01 "alien"
   34   24    7    1.76   13.97  0.9088E+01  0.9360E+01 "alien"
   33    3   15    1.73   13.69  0.8905E+01  0.7430E+01 "alien"
   29   22   12    1.78   13.31  0.8655E+01  0.7771E+01 "alien"
   31    8   20    1.51   13.29  0.1728E+02  0.7268E+01 "alien"
   40   29    2    1.56   13.17  0.1712E+02  0.9666E+01 "alien"
   23   23   18    1.60   13.15  0.2587E+01  0.1085E+02 "alien"
   26   26   18    1.51   12.88  0.4343E+01  0.8136E+01 "alien"
   42    0    7    1.75   12.81  0.2022E+02  0.8498E+01 "alien"
   35    6   18    1.53   12.63  0.1642E+02  0.1435E+02 "alien"
   42   32    0    1.47   12.53  0.4127E+01  0.1311E+02 "alien"
   36   24    2    1.78   11.88  0.7726E+01  0.8608E+01 "alien"
   15   12   18    1.90   11.58  0.4091E+01  0.5407E+01 "alien"
   36   25    5    1.72   11.37  0.7393E+01  0.8331E+01 "alien"
   26   13   24    1.38   11.36 -0.1269E+02  0.9707E+01 "alien"
   40    0    8    1.80   11.22 -0.1930E+01  0.8065E+01 "alien"
   22   14   22    1.52   10.57  0.1374E+02  0.9337E+01 "alien"
   48    7   15    1.36   10.57  0.2589E+02  0.1031E+02 "alien"
   39   29    5    1.56   10.56  0.1373E+02  0.8789E+01 "alien"
   27   20   22    1.40   10.55 -0.1178E+02  0.7498E+01 "alien"
   17   11   20    1.73   10.53  0.6851E+01  0.7058E+01 "alien"
   38   33   12    1.39   10.46 -0.1169E+02  0.9650E+01 "alien"
   33    1   12    1.90   10.32  0.3647E+01  0.8363E+01 "alien"
   29   20   10    1.91   10.04  0.3546E+01  0.1394E+02 "alien"


 List of    113 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   45   10   19    1.30   16.65 -0.1905E+02  0.1299E+02
   50    8   15    1.32   32.77 -0.3751E+02  0.1804E+02
   32   28   20    1.33   18.40 -0.2106E+02  0.1442E+02
   32   30   19    1.34   15.73  0.3852E+02  0.1496E+02
   48    7   15    1.36   10.57  0.2589E+02  0.1031E+02
   26   13   24    1.38   11.36 -0.1269E+02  0.9707E+01
   37   33   13    1.38   16.74 -0.1870E+02  0.9703E+01
   29   26   20    1.39   17.15 -0.1915E+02  0.1570E+02
   38   33   12    1.39   10.46 -0.1169E+02  0.9650E+01
   45    7   16    1.39   16.55 -0.1849E+02  0.1069E+02
   27   20   22    1.40   10.55 -0.1178E+02  0.7498E+01
   46    6   15    1.40   15.60 -0.1743E+02  0.9285E+01
   37   31   11    1.46   19.55  0.2704E+01  0.1538E+02
   23   13   23    1.46   89.73  0.1241E+02  0.9305E+01
   39   31    9    1.46   16.75  0.2317E+01  0.1053E+02
   26   23   20    1.46  120.71  0.1669E+02  0.1321E+02
   35   30   13    1.46   56.67  0.7837E+01  0.9012E+01
   27   25   19    1.46  143.87  0.1990E+02  0.1343E+02
   42   32    1    1.47   45.52  0.6295E+01  0.9665E+01
   42   32    0    1.47   12.53  0.4127E+01  0.1311E+02
   38    7   18    1.47   61.98  0.8571E+01  0.8761E+01
   43    5   15    1.47   19.81  0.2740E+01  0.9699E+01
   36   30   12    1.47   60.25  0.8332E+01  0.9089E+01
   27   27   18    1.48   24.79  0.2145E+02  0.1389E+02
   44    5   14    1.48   76.25  0.1055E+02  0.1628E+02
   39   31    8    1.48   95.99  0.1327E+02  0.9636E+01
   52    2    1    1.49   38.71  0.5353E+01  0.9198E+01
   39    6   17    1.49   75.70  0.1047E+02  0.8285E+01
   30   28   16    1.49  221.62  0.3065E+02  0.1473E+02
   41   31    4    1.49   58.39  0.8075E+01  0.9375E+01
   42    5   15    1.49   61.57  0.8514E+01  0.8590E+01
   32    8   20    1.49   30.15  0.4169E+01  0.7916E+01
   29    9   21    1.50   90.66  0.1254E+02  0.6903E+01
   41   31    3    1.50   17.18  0.2376E+01  0.9013E+01
   38   30    9    1.50  102.45  0.1417E+02  0.9115E+01
   25   10   22    1.50   83.76  0.1158E+02  0.9198E+01
   50    2    6    1.51  103.88  0.1437E+02  0.9338E+01
   26   26   18    1.51   12.88  0.4343E+01  0.8136E+01
   41   31    1    1.51   24.52  0.3390E+01  0.9452E+01
   41   31    0    1.51   19.47  0.1076E+02  0.1369E+02
   51    2    3    1.51   28.34  0.3919E+01  0.9673E+01
   38    6   17    1.51   25.79  0.3567E+01  0.9057E+01
   31    8   20    1.51   13.29  0.1728E+02  0.7268E+01
   36   29   11    1.51   20.53  0.2668E+02  0.1392E+02
   22   14   22    1.52   10.57  0.1374E+02  0.9337E+01
   35    6   18    1.53   12.63  0.1642E+02  0.1435E+02
   39   29    5    1.56   10.56  0.1373E+02  0.8789E+01
   40   29    2    1.56   13.17  0.1712E+02  0.9666E+01
   46    2    9    1.57  190.93  0.1513E+02  0.8654E+01
   45    2   10    1.57  100.65  0.7974E+01  0.9840E+01
   21   15   21    1.57  118.10  0.9357E+01  0.9846E+01
   39   29    4    1.57   73.64  0.5834E+01  0.1065E+02
   23   21   19    1.58   38.06  0.3015E+01  0.7028E+01
   48    1    5    1.58   48.18  0.3817E+01  0.8669E+01
   38   28    6    1.59   90.11  0.7139E+01  0.1253E+02
   39   28    4    1.59   83.05  0.6580E+01  0.1091E+02
   25   24   17    1.60  105.56  0.8363E+01  0.1208E+02
   35   27   10    1.60   64.18  0.5084E+01  0.8600E+01
   23   23   18    1.60   13.15  0.2587E+01  0.1085E+02
   36   27    9    1.60  149.82  0.1187E+02  0.9588E+01
   48    1    3    1.60   52.11  0.4128E+01  0.1044E+02
   43    2   11    1.60  137.48  0.1089E+02  0.9806E+01
   34    5   17    1.60   90.77  0.7191E+01  0.7447E+01
   39   28    2    1.61   40.41  0.3202E+01  0.8834E+01
   46    1    7    1.61   26.99  0.2138E+01  0.8651E+01
   39   28    1    1.61   21.40  0.1695E+01  0.8733E+01
   37   27    7    1.62  139.65  0.1106E+02  0.9044E+01
   45    1    8    1.62  138.34  0.1096E+02  0.9793E+01
   36   27    8    1.62   85.69  0.6788E+01  0.9274E+01
   33    5   17    1.63  214.32  0.1698E+02  0.7386E+01
   31   25   13    1.63  108.19  0.8571E+01  0.8103E+01
   24   23   17    1.63   51.42  0.4074E+01  0.8495E+01
   28   24   15    1.63  107.89  0.8547E+01  0.8124E+01
   34    3   15    1.71   40.05  0.2605E+02  0.1322E+02
   36   25    5    1.72   11.37  0.7393E+01  0.8331E+01
   17   11   20    1.73   10.53  0.6851E+01  0.7058E+01
   33    3   15    1.73   13.69  0.8905E+01  0.7430E+01
   33   24    9    1.74   14.98  0.9742E+01  0.8446E+01
   42    0    7    1.75   12.81  0.2022E+02  0.8498E+01
   44    1    1    1.76   15.68  0.1020E+02  0.9585E+01
   34   24    7    1.76   13.97  0.9088E+01  0.9360E+01
   29   22   12    1.78   13.31  0.8655E+01  0.7771E+01
   36   24    2    1.78   11.88  0.7726E+01  0.8608E+01
   33   23    8    1.79  200.40 -0.1581E+02  0.8132E+01
   40    0    8    1.80   11.22 -0.1930E+01  0.8065E+01
   37    1   11    1.80   23.52 -0.1855E+01  0.7926E+01
   27   21   13    1.80   70.48 -0.5560E+01  0.7269E+01
   42    1    4    1.81   34.19 -0.2697E+01  0.8784E+01
   21   19   16    1.81  103.73 -0.8182E+01  0.9584E+01
   17   15   18    1.82   52.57 -0.4147E+01  0.5863E+01
   42    1    3    1.83  121.43 -0.9578E+01  0.8815E+01
   35   23    3    1.83  242.82 -0.1915E+02  0.8509E+01
   31   22    9    1.84  442.57 -0.3491E+02  0.1461E+02
   29   21   11    1.84   43.98 -0.3469E+01  0.7802E+01
   35   23    0    1.85   28.43 -0.1946E+02  0.1199E+02
   20   18   16    1.85   90.24 -0.7118E+01  0.9353E+01
   40    1    6    1.85  106.70 -0.8416E+01  0.8533E+01
   17   17   17    1.86   16.19 -0.3949E+01  0.6630E+01
   33   22    6    1.87  273.75 -0.2159E+02  0.8259E+01
   40    1    5    1.88   23.17 -0.1828E+01  0.1156E+02
   36    0   10    1.88   15.42 -0.3553E+01  0.8114E+01
   34   22    3    1.89   24.63  0.8705E+01  0.8550E+01
   28   20   11    1.90   16.89  0.5967E+01  0.8020E+01
   33    1   12    1.90   10.32  0.3647E+01  0.8363E+01
   21   18   15    1.90   14.11  0.4985E+01  0.7951E+01
   15   12   18    1.90   11.58  0.4091E+01  0.5407E+01
   25   19   13    1.90   16.14  0.5704E+01  0.7953E+01
   29   20   10    1.91   10.04  0.3546E+01  0.1394E+02
   26   19   12    1.93   38.49  0.1360E+02  0.7914E+01
   38    1    6    1.94   60.77  0.2147E+02  0.1030E+02
   39    2    2    1.97   22.60  0.7987E+01  0.8954E+01
   39    2    1    1.98   14.54  0.5140E+01  0.8895E+01
   29    1   13    1.99   28.10  0.9929E+01  0.7937E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        71
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1323
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1169
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      154


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
