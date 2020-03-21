
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.305    77.305    38.687  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=       4       4
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       4
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  1.000000 -0.000591  0.000708
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
 REFLECTING_RANGE_E.S.D.=   0.072
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

    1548 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       4 ...       4
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1548 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               49
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

 *  44        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  25        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  23        oC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  20        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    27        mC        500.0     159.4   38.7  109.3  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.7  109.3  116.0  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.3  116.0   38.7  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        998.2     172.9   77.3   38.7  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.2      77.3  172.9   38.7  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.2      77.3  172.9   38.7  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.2      77.3   77.3   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.3   86.4   86.4  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.4   86.4  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.4  247.5 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.4   86.4   77.3  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.7   77.3  177.1  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.4   86.4  177.1 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.4   86.4   86.4  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.1  109.3   38.7  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.7   77.3  177.1 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.4  177.1   77.3 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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
 UNIT_CELL_A-AXIS=    73.676   -23.404     0.029
 UNIT_CELL_B-AXIS=    23.404    73.676     0.011
 UNIT_CELL_C-AXIS=    -0.015    -0.001    38.688
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       355 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.22
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.242    77.242    38.657  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.3E-01 1.3E-01 6.4E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012946  0.012946  0.025868  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.613   -23.397     0.016
 COORDINATES OF UNIT CELL B-AXIS    23.397    73.613     0.004
 COORDINATES OF UNIT CELL C-AXIS    -0.008     0.000    38.657
 CRYSTAL MOSAICITY (DEGREES)     0.072
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006279  0.002017
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.85   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1260.92   1242.01
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.10
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.24    77.24    38.66  90.000  90.000  90.000


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
 -0.05  0.05         2          6.157      54.516
  0.05  0.15       136        456.514       0.735
  0.15  0.25       189        303.249       1.107
  0.25  0.35       181        306.639       1.095
  0.35  0.45       155        168.190       1.996
  0.45  0.55       186        370.145       0.907
  0.55  0.65       197        337.320       0.995
  0.65  0.75       178        208.711       1.608
  0.75  0.85       169        307.376       1.092
  0.85  0.95       155        607.638       0.552
 
 -0.05  0.95      1548        335.635       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               241.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.897
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     3.0 XMAX=     4.0 NXBIN=    1
 YMIN= 0.00168 YMAX= 0.59774 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        536


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               241.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.895
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    11.9 XMAX=  2454.1 NXBIN=    3
 YMIN=     6.4 YMAX=  2520.0 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        536


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               241.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.891
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     3.0 XMAX=     4.0 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1263
 DETECTOR_SURFACE_POSITION=    1648    1690
 DETECTOR_SURFACE_POSITION=     818    1690
 DETECTOR_SURFACE_POSITION=     818     836
 DETECTOR_SURFACE_POSITION=    1648     836
 DETECTOR_SURFACE_POSITION=    2172    1663
 DETECTOR_SURFACE_POSITION=    1622    2229
 DETECTOR_SURFACE_POSITION=     844    2229
 DETECTOR_SURFACE_POSITION=     294    1663
 DETECTOR_SURFACE_POSITION=     294     863
 DETECTOR_SURFACE_POSITION=     844     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2172     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        536


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
 1.007E+00  3.486E-03   16.87


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

   24.426   3.394    15.67   0.99      2.48      3.40      60      91      24
    3.394   2.411    13.25   1.17      4.45      3.69     124     167      20
    2.411   1.972     9.55   1.09      5.31      4.50     178     222      10
    1.972   1.709     0.01   0.73    -99.90    -99.90     198     250       0
    1.709   1.530     0.04   0.57    916.30   1183.83     216     293       0
    1.530   1.397    -0.06   0.72    -99.90    -99.90     250     319       0
    1.397   1.293     0.13   0.47     42.06     54.93     100     149       2
   --------------------------------------------------------------------------
   24.426   1.293     3.87   0.79      5.01      5.31    1126    1491      56
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       4       4)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.679  -0.6333E+01  0.1080E+02        -0.59           1 
   47    0    0     1.643   0.9790E+01  0.1134E+02         0.86           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     2.7% OF MEAN INTENSITY


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

     3.84          68      46      1285        3.6%       2.3%      3.4%       44   19.03      3.3%   100.0*     0    0.000       0
     2.73         106      68      2098        3.2%       4.3%      3.6%       76   17.41      6.1%    99.1*     0    0.000       0
     2.23         162     100      2665        3.8%       4.4%      3.8%      124   15.65      6.3%    99.1*     0    0.000       0
     1.94         170     101      3101        3.3%       5.8%      5.4%      138    8.80      8.2%    99.2*     0    0.000       0
     1.73         192     116      3481        3.3%   18105.3%  20813.2%      152    0.03  25604.8%     0.8      0    0.000       0
     1.58         225     138      3848        3.6%     -99.9%    -99.9%      174    0.00    -99.9%    18.1      0    0.000       0
     1.47         253     159      4143        3.8%     -99.9%    -99.9%      188    0.00    -99.9%    -4.8      0    0.000       0
     1.37         226     137      4429        3.1%     -99.9%    -99.9%      178    0.00    -99.9%    10.5      0    0.000       0
     1.29          89      63      4740        1.3%      24.9%     28.9%       52    0.40     35.2%    99.6*     0    0.000       0
    total        1491     928     29790        3.1%       5.0%      5.3%     1126    4.87      7.1%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1548
 NUMBER OF REJECTED MISFITS                                56
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1491
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    928



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

   WILSON LINE (using all data) : A=   8.537 B=  35.980 CORRELATION=  0.81
      #      RES      SS        <I>       log(<I>)       BO
      18     7.404   0.005  4.0194E+03       8.299      26.1
      23     4.641   0.012  3.3472E+03       8.116      18.1
      32     3.507   0.020  2.1076E+03       7.653      21.7
      31     2.990   0.028  8.4995E+02       6.745      32.0
      40     2.648   0.036  8.8797E+02       6.789      24.5
      55     2.394   0.044  6.0795E+02       6.410      24.4
      43     2.203   0.052  4.2014E+02       6.041      24.2
      49     2.054   0.059  2.1584E+02       5.375      26.7
      51     1.930   0.067  1.2696E+00       0.239      61.8
      63     1.825   0.075 -1.9038E-01    -999.999   -1000.0
      56     1.736   0.083  2.6666E-01      -1.322      59.4
      61     1.662   0.091 -2.7360E-01    -999.999   -1000.0
      64     1.594   0.098  1.5662E+00       0.449      41.1
      83     1.534   0.106 -1.0138E+00    -999.999   -1000.0
      81     1.478   0.114  3.6529E-01      -1.007      41.7
      61     1.430   0.122 -9.4859E-01    -999.999   -1000.0
      56     1.389   0.130 -6.9314E-01    -999.999   -1000.0
      41     1.349   0.137  1.7878E+01       2.884      20.6
      19     1.313   0.145 -7.4286E-01    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       6     7.404        1.532        1.090        0.618
       3     4.641        0.971        0.513        0.202
       4     3.507        4.622        6.583        6.864
       4     2.990        0.904        0.366        0.112
       3     2.648        0.230        0.052        0.009
       5     2.394        0.654        0.274        0.086
       3     2.203        1.119        0.415        0.111
       8     2.054        0.229        0.073        0.018
      11     1.662       20.729       27.204       82.852
       4     1.594       12.108       19.266       22.492
       3     1.534       12.755       27.205       41.566
       6     1.478      325.286     3568.841    30782.063
       1     1.389       31.490       61.215       84.998
      95   overall       36.452      128.360     2650.882


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.404        0.827        0.977        0.950
      20     4.641        0.984        0.984        0.868
      28     3.507        0.580        0.382        0.225
      27     2.990        0.821        0.633        0.427
      37     2.648        1.116        1.160        1.080
      50     2.394        0.952        0.781        0.537
      40     2.203        0.788        0.574        0.359
      41     2.054        1.300        1.436        1.431
      44     1.930       14.640       30.869      115.395
      52     1.825       20.486        0.919      144.780
      45     1.736      543.057     5391.356   200407.110
      60     1.594       10.050       14.747       71.225
      39     1.349       18.489      231.343     2187.816
     832   overall       48.475      246.746    12589.755

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       6     7.404  0.00  0.00  0.00  0.35  0.32  0.59  0.56  0.53  0.51  0.49
       3     4.641  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       4     3.507  1.01  0.72  0.60  0.53  0.48  0.45  0.42  0.40  0.38  1.10
       4     2.990  0.00  0.00  0.00  0.00  0.48  0.45  0.42  0.40  0.38  0.37
       3     2.648  1.34  0.96  0.80  0.70  0.64  0.59  1.12  1.06  1.01  0.98
       5     2.394  0.00  0.00  0.00  0.85  0.77  1.07  1.00  0.95  0.91  0.88
       3     2.203  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.054  2.01  1.81  1.50  1.32  1.20  1.11  1.05  1.19  1.14  1.10
      11     1.662  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.72  0.69  0.67
       4     1.594  0.00  0.00  0.00  0.00  0.00  0.00  0.42  0.40  0.38  0.37
       3     1.534  2.68  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       6     1.478  0.67  0.48  0.40  0.35  0.32  0.30  0.28  0.27  0.25  0.24
       1     1.389  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      95   overall  0.68  0.55  0.45  0.47  0.44  0.45  0.46  0.45  0.45  0.46


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.404  4.38  2.30  2.57  2.27  1.91  1.66  1.49  1.36  1.26  1.19
      20     4.641  1.05  1.38  0.96  0.91  1.02  1.11  0.99  1.09  1.01  1.03
      28     3.507  1.13  0.79  0.83  0.87  1.09  1.11  1.06  1.17  1.14  1.13
      27     2.990  1.56  1.43  1.14  1.12  0.94  0.99  1.18  1.08  1.06  1.00
      37     2.648  0.85  1.34  1.25  1.07  0.96  0.84  0.97  1.08  1.05  1.03
      50     2.394  1.47  1.43  1.08  1.15  1.02  1.11  1.11  1.05  1.04  0.98
      40     2.203  1.05  1.38  1.06  0.83  0.89  0.94  0.89  1.00  1.05  1.03
      41     2.054  1.54  1.08  1.32  1.11  0.99  0.86  0.97  0.97  0.95  0.93
      44     1.930  4.30  2.26  1.67  1.31  1.10  0.96  0.90  0.83  0.77  0.72
      52     1.825  5.66  2.97  2.08  1.63  1.37  1.19  1.07  0.98  0.91  0.85
      45     1.736  5.14  2.70  1.89  1.48  1.24  1.08  0.97  0.93  0.86  0.81
      60     1.594  4.73  2.57  1.80  1.42  1.23  1.07  0.96  0.91  0.84  0.79
      39     1.349  5.12  2.97  2.18  1.94  1.95  1.76  1.68  1.68  1.60  1.50
     832   overall  1.87  1.18  0.91  0.78  0.71  0.66  0.64  0.63  0.60  0.57


 List of    139 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   32   23   10    1.75  102.53  0.2734E+02  0.1311E+02 "alien"
   33   26   11    1.63   86.89 -0.2377E+02  0.1407E+02 "alien"
   34   24    7    1.76   78.17  0.2084E+02  0.8246E+01 "alien"
   24   20   15    1.78   65.72 -0.1251E+02  0.1117E+02 "alien"
   37   26    2    1.70   62.94  0.1678E+02  0.1554E+02 "alien"
   34    1   12    1.86   62.65 -0.1193E+02  0.6873E+01 "alien"
   39    1   11    1.72   55.40  0.1477E+02  0.1153E+02 "alien"
   42    1    4    1.81   49.40 -0.9406E+01  0.7637E+01 "alien"
   33   23    7    1.81   48.00 -0.9139E+01  0.1194E+02 "alien"
   25   20   14    1.82   46.36 -0.8826E+01  0.6902E+01 "alien"
   39    1   10    1.76   45.54  0.1214E+02  0.1075E+02 "alien"
   34   26   10    1.64   41.69 -0.1141E+02  0.7812E+01 "alien"
   22   19   15    1.85   40.15 -0.7645E+01  0.1036E+02 "alien"
   31   24   12    1.68   40.13 -0.1098E+02  0.7375E+01 "alien"
   45    1    7    1.64   40.08 -0.1096E+02  0.7597E+01 "alien"
   34   23    6    1.81   40.07 -0.7629E+01  0.8100E+01 "alien"
   41    1    4    1.85   39.52 -0.7524E+01  0.7843E+01 "alien"
   41   31    4    1.48   39.24  0.1434E+02  0.9192E+01 "alien"
   30    3   15    1.82   38.67 -0.7363E+01  0.8591E+01 "alien"
   54    5    9    1.35   37.83  0.6763E+03  0.3336E+02 "alien"
   18   10   21    1.65   36.41 -0.9962E+01  0.5496E+01 "alien"
   17   15   18    1.82   35.66 -0.6789E+01  0.4931E+01 "alien"
   31    2   14    1.85   35.33 -0.6727E+01  0.6364E+01 "alien"
   33    3   15    1.73   33.65  0.8972E+01  0.6891E+01 "alien"
   49    4   10    1.46   33.21  0.1213E+02  0.1311E+02 "alien"
   36   25    5    1.72   32.46  0.8655E+01  0.8066E+01 "alien"
   31   22    9    1.84   32.31 -0.6151E+01  0.8323E+01 "alien"
   33   25   11    1.65   32.23 -0.8817E+01  0.1130E+02 "alien"
   38    6   17    1.51   30.40  0.1110E+02  0.7434E+01 "alien"
   31   28   15    1.50   29.07  0.1062E+02  0.1225E+02 "alien"
   32   28   14    1.52   29.01 -0.2941E+02  0.1374E+02 "alien"
   32   31   17    1.38   28.84 -0.1999E+02  0.1560E+02 "alien"
   36   26    7    1.66   28.35 -0.7755E+01  0.8136E+01 "alien"
   44    1    8    1.65   27.86 -0.7623E+01  0.7513E+01 "alien"
   40    2   12    1.65   27.85 -0.7620E+01  0.7122E+01 "alien"
   33   11   22    1.38   27.84 -0.1930E+02  0.7115E+01 "alien"
   51    3    6    1.47   27.65  0.1010E+02  0.1486E+02 "alien"
   51    8   15    1.29   27.64 -0.2053E+02  0.1265E+02 "alien"
   34   23    5    1.83   27.56 -0.5248E+01  0.8225E+01 "alien"
   37    2   13    1.71   26.76  0.7137E+01  0.8069E+01 "alien"
   30   22   11    1.79   26.43 -0.5032E+01  0.7438E+01 "alien"
   42    1    2    1.83   26.39 -0.5024E+01  0.7997E+01 "alien"
   39    6   17    1.48   25.93  0.9472E+01  0.7553E+01 "alien"
   53    2    1    1.46   25.59  0.9346E+01  0.1071E+02 "alien"
   39   33   10    1.41   25.33 -0.1755E+02  0.9505E+01 "alien"
   52    3    5    1.46   24.52  0.8958E+01  0.8861E+01 "alien"
   16    9   20    1.76   24.14  0.6438E+01  0.5063E+01 "alien"
   27   23   15    1.66   24.04 -0.6576E+01  0.7060E+01 "alien"
   28   21   12    1.82   23.96 -0.4561E+01  0.7186E+01 "alien"
   29   26   20    1.38   23.42 -0.1623E+02  0.1299E+02 "alien"
   32   28   20    1.32   23.12 -0.1718E+02  0.1160E+02 "alien"
   53    7   12    1.32   23.06 -0.1713E+02  0.1601E+02 "alien"
   24   19   14    1.86   22.98 -0.4375E+01  0.9404E+01 "alien"
   43    2   10    1.63   22.96 -0.6283E+01  0.8707E+01 "alien"
   34    3   15    1.70   22.67  0.6046E+01  0.1046E+02 "alien"
   37   29   10    1.51   22.63 -0.2294E+02  0.1304E+02 "alien"
   26   23   16    1.64   22.61 -0.6185E+01  0.1172E+02 "alien"
   52    2    3    1.47   22.53  0.8230E+01  0.8600E+01 "alien"
   36   26    6    1.68   22.13 -0.6053E+01  0.1178E+02 "alien"
   42    1    1    1.84   22.11 -0.4209E+01  0.8029E+01 "alien"
   40   11   21    1.31   21.74 -0.1615E+02  0.1117E+02 "alien"
   38    9   20    1.38   21.19 -0.1469E+02  0.1054E+02 "alien"
   33   28   13    1.53   20.49 -0.2077E+02  0.1271E+02 "alien"
   32   11   22    1.39   20.43 -0.1416E+02  0.9281E+01 "alien"
   33    4   16    1.67   20.36 -0.5571E+01  0.6231E+01 "alien"
   31    3   15    1.79   20.08 -0.3824E+01  0.6286E+01 "alien"
   26   13   24    1.38   19.24 -0.1333E+02  0.7996E+01 "alien"
   50    3    8    1.47   19.04  0.6954E+01  0.8686E+01 "alien"
   36   25    4    1.73   18.80  0.5013E+01  0.8270E+01 "alien"
   35   26    8    1.66   18.70 -0.5116E+01  0.1290E+02 "alien"
   42   32    0    1.46   18.16  0.2386E+02  0.1353E+02 "alien"
   42   32    1    1.46   18.15  0.6628E+01  0.9486E+01 "alien"
   24   16   22    1.47   17.64  0.6443E+01  0.9058E+01 "alien"
   36    9   20    1.42   17.43 -0.1654E+02  0.9639E+01 "alien"
   26   24   19    1.49   17.37  0.6346E+01  0.1011E+02 "alien"
   21    7   24    1.46   17.30  0.6320E+01  0.5502E+01 "alien"
   50    4   10    1.43   17.09 -0.1621E+02  0.8530E+01 "alien"
   49    5   12    1.41   16.81 -0.1594E+02  0.8470E+01 "alien"
   49    3    9    1.48   16.53  0.6037E+01  0.8267E+01 "alien"
   36    0   10    1.88   16.26 -0.9631E+01  0.7286E+01 "alien"
   27   24   20    1.43   16.18 -0.1535E+02  0.1035E+02 "alien"
   54    4    6    1.39   16.12 -0.1117E+02  0.8892E+01 "alien"
   27   20   22    1.40   16.09 -0.1115E+02  0.7168E+01 "alien"
   43    5   15    1.47   15.26  0.5573E+01  0.8879E+01 "alien"
   20   17   19    1.67   15.20 -0.4159E+01  0.8005E+01 "alien"
   32   29   15    1.47   14.99  0.5475E+01  0.8268E+01 "alien"
   35   24    5    1.77   14.94  0.3985E+01  0.8078E+01 "alien"
   35   23    0    1.84   14.89 -0.7963E+01  0.1138E+02 "alien"
   29   27   16    1.52   14.81 -0.1502E+02  0.1153E+02 "alien"
   29   20   10    1.91   14.75  0.1873E+02  0.1165E+02 "alien"
   29   24   14    1.65   14.64 -0.4005E+01  0.1067E+02 "alien"
   52    8   14    1.30   14.41 -0.1070E+02  0.1519E+02 "alien"
   43    1    0    1.80   14.38 -0.7419E+01  0.1078E+02 "alien"
   36   32   14    1.39   14.31 -0.9920E+01  0.1476E+02 "alien"
   36   30   12    1.47   14.24  0.5203E+01  0.8693E+01 "alien"
   25   21   15    1.74   13.75  0.3666E+01  0.9194E+01 "alien"
   41    4   14    1.55   13.70 -0.1389E+02  0.7345E+01 "alien"
   35   30   13    1.46   13.70  0.5006E+01  0.8495E+01 "alien"
   41   32    4    1.47   13.63  0.4980E+01  0.9325E+01 "alien"
   40   31    5    1.50   13.54  0.4946E+01  0.9759E+01 "alien"
   22    8   24    1.45   13.47 -0.1278E+02  0.6534E+01 "alien"
   17    8   21    1.68   13.43 -0.3674E+01  0.5159E+01 "alien"
   25   13   24    1.39   13.36 -0.9259E+01  0.5306E+01 "alien"
   47    1    6    1.59   13.31  0.2084E+02  0.1213E+02 "alien"
   44    0    6    1.69   13.20 -0.9003E+01  0.7801E+01 "alien"
   41   11   21    1.29   13.13 -0.9750E+01  0.1169E+02 "alien"
   40    0    8    1.79   13.09 -0.6159E+01  0.7448E+01 "alien"
   28   19   10    1.97   12.90  0.1638E+02  0.1454E+02 "alien"
   51    3    8    1.44   12.89 -0.1222E+02  0.1024E+02 "alien"
   19   14   20    1.66   12.86 -0.3517E+01  0.4935E+01 "alien"
   49    6   13    1.38   12.72 -0.8815E+01  0.1224E+02 "alien"
   53    2    0    1.46   12.48  0.1081E+02  0.1393E+02 "alien"
   31   30   17    1.41   12.39 -0.8588E+01  0.9399E+01 "alien"
   39   30    6    1.53   12.32 -0.1249E+02  0.9039E+01 "alien"
   26   22   15    1.70   12.25  0.3267E+01  0.7437E+01 "alien"
   54    3    4    1.41   12.18 -0.1155E+02  0.9097E+01 "alien"
   33   24    9    1.73   12.17  0.3245E+01  0.7625E+01 "alien"
   35   25    7    1.71   12.14  0.3237E+01  0.8155E+01 "alien"
   47    0    2    1.64   12.06 -0.7603E+01  0.7794E+01 "alien"
   47    3   10    1.51   12.03 -0.1220E+02  0.8131E+01 "alien"
   53    3    3    1.45   11.98 -0.1137E+02  0.9390E+01 "alien"
   24   13   24    1.40   11.91 -0.8253E+01  0.7144E+01 "alien"
   42    3   12    1.59   11.84  0.1854E+02  0.1216E+02 "alien"
   50    5   12    1.39   11.58 -0.8027E+01  0.9073E+01 "alien"
   26   26   18    1.50   11.47  0.9283E+01  0.7735E+01 "alien"
   39    5   16    1.52   11.43 -0.1159E+02  0.7285E+01 "alien"
   53    5    9    1.37   11.41 -0.7907E+01  0.8787E+01 "alien"
   45    0    5    1.68   11.40 -0.6884E+01  0.7835E+01 "alien"
   25   19   13    1.90   11.34  0.1439E+02  0.7416E+01 "alien"
   24    9   21    1.57   11.34  0.1775E+02  0.7316E+01 "alien"
   38   30    9    1.50   11.28  0.4121E+01  0.9082E+01 "alien"
   40   32    8    1.44   10.95 -0.1039E+02  0.8867E+01 "alien"
   51    5   11    1.39   10.70 -0.7417E+01  0.8682E+01 "alien"
   29   28   17    1.47   10.68  0.3900E+01  0.7842E+01 "alien"
   52    3    6    1.45   10.61 -0.1007E+02  0.8558E+01 "alien"
   36   29   11    1.51   10.60 -0.1075E+02  0.8331E+01 "alien"
   34   27   11    1.59   10.56  0.1653E+02  0.8528E+01 "alien"
   37   33   13    1.38   10.55 -0.7311E+01  0.9026E+01 "alien"
   42    3   13    1.56   10.13 -0.1027E+02  0.1335E+02 "alien"


 List of    139 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   51    8   15    1.29   27.64 -0.2053E+02  0.1265E+02
   41   11   21    1.29   13.13 -0.9750E+01  0.1169E+02
   52    8   14    1.30   14.41 -0.1070E+02  0.1519E+02
   40   11   21    1.31   21.74 -0.1615E+02  0.1117E+02
   53    7   12    1.32   23.06 -0.1713E+02  0.1601E+02
   32   28   20    1.32   23.12 -0.1718E+02  0.1160E+02
   54    5    9    1.35   37.83  0.6763E+03  0.3336E+02
   53    5    9    1.37   11.41 -0.7907E+01  0.8787E+01
   26   13   24    1.38   19.24 -0.1333E+02  0.7996E+01
   33   11   22    1.38   27.84 -0.1930E+02  0.7115E+01
   32   31   17    1.38   28.84 -0.1999E+02  0.1560E+02
   37   33   13    1.38   10.55 -0.7311E+01  0.9026E+01
   38    9   20    1.38   21.19 -0.1469E+02  0.1054E+02
   29   26   20    1.38   23.42 -0.1623E+02  0.1299E+02
   49    6   13    1.38   12.72 -0.8815E+01  0.1224E+02
   51    5   11    1.39   10.70 -0.7417E+01  0.8682E+01
   36   32   14    1.39   14.31 -0.9920E+01  0.1476E+02
   50    5   12    1.39   11.58 -0.8027E+01  0.9073E+01
   25   13   24    1.39   13.36 -0.9259E+01  0.5306E+01
   32   11   22    1.39   20.43 -0.1416E+02  0.9281E+01
   54    4    6    1.39   16.12 -0.1117E+02  0.8892E+01
   27   20   22    1.40   16.09 -0.1115E+02  0.7168E+01
   24   13   24    1.40   11.91 -0.8253E+01  0.7144E+01
   31   30   17    1.41   12.39 -0.8588E+01  0.9399E+01
   39   33   10    1.41   25.33 -0.1755E+02  0.9505E+01
   49    5   12    1.41   16.81 -0.1594E+02  0.8470E+01
   54    3    4    1.41   12.18 -0.1155E+02  0.9097E+01
   36    9   20    1.42   17.43 -0.1654E+02  0.9639E+01
   50    4   10    1.43   17.09 -0.1621E+02  0.8530E+01
   27   24   20    1.43   16.18 -0.1535E+02  0.1035E+02
   40   32    8    1.44   10.95 -0.1039E+02  0.8867E+01
   51    3    8    1.44   12.89 -0.1222E+02  0.1024E+02
   52    3    6    1.45   10.61 -0.1007E+02  0.8558E+01
   53    3    3    1.45   11.98 -0.1137E+02  0.9390E+01
   22    8   24    1.45   13.47 -0.1278E+02  0.6534E+01
   53    2    1    1.46   25.59  0.9346E+01  0.1071E+02
   49    4   10    1.46   33.21  0.1213E+02  0.1311E+02
   53    2    0    1.46   12.48  0.1081E+02  0.1393E+02
   52    3    5    1.46   24.52  0.8958E+01  0.8861E+01
   35   30   13    1.46   13.70  0.5006E+01  0.8495E+01
   42   32    1    1.46   18.15  0.6628E+01  0.9486E+01
   21    7   24    1.46   17.30  0.6320E+01  0.5502E+01
   42   32    0    1.46   18.16  0.2386E+02  0.1353E+02
   29   28   17    1.47   10.68  0.3900E+01  0.7842E+01
   43    5   15    1.47   15.26  0.5573E+01  0.8879E+01
   36   30   12    1.47   14.24  0.5203E+01  0.8693E+01
   41   32    4    1.47   13.63  0.4980E+01  0.9325E+01
   50    3    8    1.47   19.04  0.6954E+01  0.8686E+01
   24   16   22    1.47   17.64  0.6443E+01  0.9058E+01
   32   29   15    1.47   14.99  0.5475E+01  0.8268E+01
   51    3    6    1.47   27.65  0.1010E+02  0.1486E+02
   52    2    3    1.47   22.53  0.8230E+01  0.8600E+01
   49    3    9    1.48   16.53  0.6037E+01  0.8267E+01
   39    6   17    1.48   25.93  0.9472E+01  0.7553E+01
   41   31    4    1.48   39.24  0.1434E+02  0.9192E+01
   26   24   19    1.49   17.37  0.6346E+01  0.1011E+02
   38   30    9    1.50   11.28  0.4121E+01  0.9082E+01
   40   31    5    1.50   13.54  0.4946E+01  0.9759E+01
   26   26   18    1.50   11.47  0.9283E+01  0.7735E+01
   31   28   15    1.50   29.07  0.1062E+02  0.1225E+02
   38    6   17    1.51   30.40  0.1110E+02  0.7434E+01
   36   29   11    1.51   10.60 -0.1075E+02  0.8331E+01
   47    3   10    1.51   12.03 -0.1220E+02  0.8131E+01
   37   29   10    1.51   22.63 -0.2294E+02  0.1304E+02
   29   27   16    1.52   14.81 -0.1502E+02  0.1153E+02
   32   28   14    1.52   29.01 -0.2941E+02  0.1374E+02
   39    5   16    1.52   11.43 -0.1159E+02  0.7285E+01
   39   30    6    1.53   12.32 -0.1249E+02  0.9039E+01
   33   28   13    1.53   20.49 -0.2077E+02  0.1271E+02
   41    4   14    1.55   13.70 -0.1389E+02  0.7345E+01
   42    3   13    1.56   10.13 -0.1027E+02  0.1335E+02
   24    9   21    1.57   11.34  0.1775E+02  0.7316E+01
   34   27   11    1.59   10.56  0.1653E+02  0.8528E+01
   47    1    6    1.59   13.31  0.2084E+02  0.1213E+02
   42    3   12    1.59   11.84  0.1854E+02  0.1216E+02
   43    2   10    1.63   22.96 -0.6283E+01  0.8707E+01
   33   26   11    1.63   86.89 -0.2377E+02  0.1407E+02
   34   26   10    1.64   41.69 -0.1141E+02  0.7812E+01
   26   23   16    1.64   22.61 -0.6185E+01  0.1172E+02
   47    0    2    1.64   12.06 -0.7603E+01  0.7794E+01
   45    1    7    1.64   40.08 -0.1096E+02  0.7597E+01
   29   24   14    1.65   14.64 -0.4005E+01  0.1067E+02
   33   25   11    1.65   32.23 -0.8817E+01  0.1130E+02
   44    1    8    1.65   27.86 -0.7623E+01  0.7513E+01
   18   10   21    1.65   36.41 -0.9962E+01  0.5496E+01
   40    2   12    1.65   27.85 -0.7620E+01  0.7122E+01
   36   26    7    1.66   28.35 -0.7755E+01  0.8136E+01
   35   26    8    1.66   18.70 -0.5116E+01  0.1290E+02
   27   23   15    1.66   24.04 -0.6576E+01  0.7060E+01
   19   14   20    1.66   12.86 -0.3517E+01  0.4935E+01
   20   17   19    1.67   15.20 -0.4159E+01  0.8005E+01
   33    4   16    1.67   20.36 -0.5571E+01  0.6231E+01
   45    0    5    1.68   11.40 -0.6884E+01  0.7835E+01
   36   26    6    1.68   22.13 -0.6053E+01  0.1178E+02
   17    8   21    1.68   13.43 -0.3674E+01  0.5159E+01
   31   24   12    1.68   40.13 -0.1098E+02  0.7375E+01
   44    0    6    1.69   13.20 -0.9003E+01  0.7801E+01
   34    3   15    1.70   22.67  0.6046E+01  0.1046E+02
   37   26    2    1.70   62.94  0.1678E+02  0.1554E+02
   26   22   15    1.70   12.25  0.3267E+01  0.7437E+01
   37    2   13    1.71   26.76  0.7137E+01  0.8069E+01
   35   25    7    1.71   12.14  0.3237E+01  0.8155E+01
   36   25    5    1.72   32.46  0.8655E+01  0.8066E+01
   39    1   11    1.72   55.40  0.1477E+02  0.1153E+02
   33    3   15    1.73   33.65  0.8972E+01  0.6891E+01
   33   24    9    1.73   12.17  0.3245E+01  0.7625E+01
   36   25    4    1.73   18.80  0.5013E+01  0.8270E+01
   25   21   15    1.74   13.75  0.3666E+01  0.9194E+01
   32   23   10    1.75  102.53  0.2734E+02  0.1311E+02
   16    9   20    1.76   24.14  0.6438E+01  0.5063E+01
   34   24    7    1.76   78.17  0.2084E+02  0.8246E+01
   39    1   10    1.76   45.54  0.1214E+02  0.1075E+02
   35   24    5    1.77   14.94  0.3985E+01  0.8078E+01
   24   20   15    1.78   65.72 -0.1251E+02  0.1117E+02
   31    3   15    1.79   20.08 -0.3824E+01  0.6286E+01
   30   22   11    1.79   26.43 -0.5032E+01  0.7438E+01
   40    0    8    1.79   13.09 -0.6159E+01  0.7448E+01
   43    1    0    1.80   14.38 -0.7419E+01  0.1078E+02
   42    1    4    1.81   49.40 -0.9406E+01  0.7637E+01
   34   23    6    1.81   40.07 -0.7629E+01  0.8100E+01
   33   23    7    1.81   48.00 -0.9139E+01  0.1194E+02
   17   15   18    1.82   35.66 -0.6789E+01  0.4931E+01
   25   20   14    1.82   46.36 -0.8826E+01  0.6902E+01
   30    3   15    1.82   38.67 -0.7363E+01  0.8591E+01
   28   21   12    1.82   23.96 -0.4561E+01  0.7186E+01
   34   23    5    1.83   27.56 -0.5248E+01  0.8225E+01
   42    1    2    1.83   26.39 -0.5024E+01  0.7997E+01
   42    1    1    1.84   22.11 -0.4209E+01  0.8029E+01
   31   22    9    1.84   32.31 -0.6151E+01  0.8323E+01
   35   23    0    1.84   14.89 -0.7963E+01  0.1138E+02
   31    2   14    1.85   35.33 -0.6727E+01  0.6364E+01
   41    1    4    1.85   39.52 -0.7524E+01  0.7843E+01
   22   19   15    1.85   40.15 -0.7645E+01  0.1036E+02
   34    1   12    1.86   62.65 -0.1193E+02  0.6873E+01
   24   19   14    1.86   22.98 -0.4375E+01  0.9404E+01
   36    0   10    1.88   16.26 -0.9631E+01  0.7286E+01
   25   19   13    1.90   11.34  0.1439E+02  0.7416E+01
   29   20   10    1.91   14.75  0.1873E+02  0.1165E+02
   28   19   10    1.97   12.90  0.1638E+02  0.1454E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        63
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1548
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1393
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      155


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.2 sec
