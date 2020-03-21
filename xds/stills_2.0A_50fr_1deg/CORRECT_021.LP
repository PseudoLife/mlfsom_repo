
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.453    77.453    38.764  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      21      21
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      21
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000860  0.000692
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
 BEAM_DIVERGENCE_E.S.D.=   0.045
 REFLECTING_RANGE_E.S.D.=   0.114
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

    1407 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      21 ...      21
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1407 REFLECTIONS ACCEPTED



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

 *  44        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.5   38.8   77.5  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  33        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  32        oP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.5  109.5   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.5  109.5   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.5  109.5   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.7   38.8   77.5  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.7   38.8   77.5  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  159.7   77.5  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.7   77.5  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  159.7   77.5  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  159.7   77.5  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     159.7   38.8  109.5  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.5  116.2  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  159.7  159.7  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.5  116.2   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        998.1     173.2   77.5   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.1      77.5  173.2   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.1      77.5  173.2   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.1      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.2  116.2  116.2  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.6   86.6  116.2 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.6   86.6  116.2 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.5   86.6   86.6  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.6   86.6  109.5  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.6  248.0 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.6   86.6   77.5  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.6   86.6   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.5  177.5  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.6   86.6  177.5 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.6   86.6   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.6   86.6   86.6  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.5  109.5   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.5  177.5 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.6  177.5   77.5 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    5.2      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.810   -23.473     0.044
 UNIT_CELL_B-AXIS=    23.473    73.810     0.019
 UNIT_CELL_C-AXIS=    -0.024    -0.002    38.764
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       241 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.26
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.426    77.426    38.747  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.2E-02 2.2E-02 1.3E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012915  0.012915  0.025809  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.786   -23.462     0.036
 COORDINATES OF UNIT CELL B-AXIS    23.462    73.786     0.014
 COORDINATES OF UNIT CELL C-AXIS    -0.019    -0.001    38.747
 CRYSTAL MOSAICITY (DEGREES)     0.114
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006553  0.001794
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.90   1232.90
 DETECTOR ORIGIN (PIXELS) AT                     1261.29   1242.44
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.26
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.43    77.43    38.75  90.000  90.000  90.000


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
  0.05  0.15        83        109.208       1.501
  0.15  0.25       172        200.523       0.817
  0.25  0.35       195        127.608       1.285
  0.35  0.45       147        109.750       1.494
  0.45  0.55       171        226.856       0.723
  0.55  0.65       192        146.402       1.120
  0.65  0.75       194        174.801       0.938
  0.75  0.85       166        130.543       1.256
  0.85  0.95        87        270.989       0.605
 
 -0.05  0.95      1407        163.916       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               238.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.891
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    20.1 XMAX=    20.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59493 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        500


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               238.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.887
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     9.2 XMAX=  2455.4 NXBIN=    3
 YMIN=     7.5 YMAX=  2520.6 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        500


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               238.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.881
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    20.1 XMAX=    20.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1648    1691
 DETECTOR_SURFACE_POSITION=     817    1691
 DETECTOR_SURFACE_POSITION=     817     837
 DETECTOR_SURFACE_POSITION=    1648     837
 DETECTOR_SURFACE_POSITION=    2173    1664
 DETECTOR_SURFACE_POSITION=    1622    2230
 DETECTOR_SURFACE_POSITION=     843    2230
 DETECTOR_SURFACE_POSITION=     292    1664
 DETECTOR_SURFACE_POSITION=     292     864
 DETECTOR_SURFACE_POSITION=     843     298
 DETECTOR_SURFACE_POSITION=    1622     298
 DETECTOR_SURFACE_POSITION=    2173     864
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
 9.335E-01  3.081E-03   18.65


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

   24.484   3.154    15.92   1.45      3.15      3.16      84     108      18
    3.154   2.239     8.89   1.15      5.37      4.98     190     223       8
    2.239   1.831     1.78   0.88     23.82     24.20     210     245       0
    1.831   1.587     0.00   0.62    -99.90    -99.90     238     281       4
    1.587   1.420    -0.09   0.87    -99.90    -99.90     284     339       2
    1.420   1.296     0.09   0.49   1522.36   2165.44     130     178       0
   --------------------------------------------------------------------------
   24.484   1.296     3.00   0.86      6.99      7.27    1136    1374      32
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      21      21)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.683   0.9003E+01  0.1254E+02         0.72           1 
   47    0    0     1.647  -0.1005E+02  0.1232E+02        -0.82           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -5.8% OF MEAN INTENSITY


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

     3.85          66      41      1285        3.2%       2.7%      3.1%       50   20.51      3.8%    99.8*     0    0.000       0
     2.74         102      60      2098        2.9%       4.6%      3.6%       84   16.93      6.5%    98.9*     0    0.000       0
     2.24         163      93      2665        3.5%       6.2%      5.7%      140   10.48      8.8%    98.8*     0    0.000       0
     1.94         159      90      3101        2.9%      17.2%     16.7%      138    3.57     24.3%    95.9*     0    0.000       0
     1.74         168      97      3479        2.8%     -99.9%    -99.9%      142    0.00    -99.9%    -7.3      0    0.000       0
     1.59         199     115      3851        3.0%     -99.9%    -99.9%      168    0.00    -99.9%    16.6      0    0.000       0
     1.47         218     129      4143        3.1%     -99.9%    -99.9%      178    0.00    -99.9%   -10.5      0    0.000       0
     1.37         213     120      4427        2.7%     -99.9%    -99.9%      186    0.00    -99.9%    22.0      0    0.000       0
     1.30          86      61      4743        1.3%     174.4%    221.2%       50    0.10    246.6%    73.7*     0    0.000       0
    total        1374     806     29792        2.7%       7.0%      7.3%     1136    3.90      9.9%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1407
 NUMBER OF REJECTED MISFITS                                32
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1374
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    806



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

   WILSON LINE (using all data) : A=   7.149 B=  29.031 CORRELATION=  0.82
      #      RES      SS        <I>       log(<I>)       BO
      17     7.261   0.005  1.9928E+03       7.597     -47.2
      20     4.586   0.012  2.2193E+03       7.705     -23.4
      26     3.519   0.020  9.8463E+02       6.892       6.4
      29     3.003   0.028  3.4952E+02       5.857      23.3
      37     2.649   0.036  3.0889E+02       5.733      19.9
      51     2.401   0.043  1.5982E+02       5.074      23.9
      38     2.210   0.051  9.3076E+01       4.533      25.5
      45     2.058   0.059  3.4626E+01       3.545      30.5
      44     1.933   0.067  2.7708E-01      -1.283      63.0
      52     1.831   0.075 -2.3862E-02    -999.999   -1000.0
      46     1.742   0.082  9.6707E-01      -0.033      43.6
      51     1.667   0.090 -1.0200E+00    -999.999   -1000.0
      52     1.598   0.098  6.6072E-01      -0.414      38.6
      65     1.537   0.106 -1.7310E-01    -999.999   -1000.0
      72     1.482   0.114 -1.6602E+00    -999.999   -1000.0
      55     1.434   0.122 -2.3273E+00    -999.999   -1000.0
      46     1.392   0.129  1.0921E+01       2.391      18.4
      42     1.352   0.137  2.8456E-01      -1.257      30.7
      18     1.316   0.144  2.3561E+00       0.857      21.8


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.261        1.840        1.249        0.677
       3     4.586        1.470        0.971        0.477
       2     3.519        0.412        0.092        0.015
       4     3.003        1.173        0.536        0.185
       3     2.649        0.375        0.099        0.019
       5     2.401        0.974        0.551        0.240
       2     2.210        1.312        0.549        0.169
       8     2.058        0.277        0.093        0.027
       7     1.933      417.844     5788.924    70179.821
      10     1.831       33.799       -2.527      207.153
       3     1.598      199.496     1036.510     4896.676
       3     1.434       19.844       28.801       63.926
       1     1.392      695.392     6352.368    41448.919
      83   overall       53.390      597.253     6606.306


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.261        0.459        0.323        0.184
      17     4.586        1.154        1.422        1.439
      24     3.519        0.826        0.657        0.494
      25     3.003        0.816        0.653        0.443
      34     2.649        1.185        1.416        1.569
      46     2.401        0.953        0.834        0.619
      36     2.210        0.841        0.688        0.497
      37     2.058        1.345        1.632        1.868
      39     1.742       22.470       31.919      209.472
      49     1.598       83.369       94.318     3195.562
      45     1.392        0.315       -0.050        0.061
      41     1.352     1816.609    96120.081  7135780.263
      18     1.316       19.145       18.077      184.179
     723   overall      130.412     5557.842   409114.877

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.261  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.32  0.30  0.29
       3     4.586  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     3.519  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       4     3.003  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.649  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       5     2.401  0.00  0.00  0.00  0.42  0.77  1.07  1.00  0.95  0.91  0.88
       2     2.210  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.058  2.01  1.81  1.50  1.32  1.44  1.34  1.26  1.19  1.14  1.10
       7     1.933  2.30  1.65  1.37  1.21  1.10  1.02  0.96  0.91  0.87  0.84
      10     1.831  2.82  2.02  1.68  1.48  1.34  1.25  1.17  1.11  1.07  1.03
       3     1.598  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       3     1.434  1.34  0.96  0.80  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       1     1.392  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      83   overall  0.92  0.70  0.58  0.56  0.56  0.56  0.52  0.52  0.50  0.48


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.261  2.63  2.76  1.93  1.52  1.91  1.66  1.66  1.51  1.40  1.32
      17     4.586  1.24  1.62  1.13  1.07  1.05  1.43  1.29  1.18  1.39  1.30
      24     3.519  0.88  0.69  0.64  0.51  0.74  0.74  0.99  0.98  0.91  0.99
      25     3.003  1.26  1.77  1.23  1.21  1.12  1.15  1.19  1.23  1.15  1.14
      34     2.649  0.62  1.14  1.13  1.16  0.97  0.98  1.05  1.12  1.14  1.16
      46     2.401  1.14  1.08  1.09  1.06  1.10  1.11  1.08  1.14  1.10  1.07
      36     2.210  0.88  0.92  0.96  1.01  0.85  0.86  1.05  0.96  1.03  1.01
      37     2.058  0.85  0.75  0.63  0.90  1.03  0.90  0.91  0.88  0.87  0.94
      39     1.742  4.58  2.40  1.68  1.32  1.11  0.97  0.87  0.79  0.73  0.73
      49     1.598  5.15  2.70  1.89  1.49  1.24  1.09  0.97  0.89  0.83  0.77
      45     1.392  5.37  3.06  2.31  2.09  1.81  1.67  1.55  1.53  1.54  1.51
      41     1.352  6.66  3.50  2.45  2.00  1.74  1.51  1.41  1.28  1.19  1.12
      18     1.316  5.25  2.76  1.93  1.52  1.27  1.11  0.99  1.01  0.94  0.88
     723   overall  1.77  1.14  0.88  0.78  0.72  0.67  0.66  0.64  0.63  0.62


 List of    108 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   36   24    1    1.79  775.33 -0.1850E+02  0.8683E+01 "alien"
   30   21   10    1.86  563.65 -0.1345E+02  0.1351E+02 "alien"
   42    1    3    1.82  531.00 -0.1267E+02  0.8865E+01 "alien"
   41    1    5    1.83  487.30 -0.1163E+02  0.8934E+01 "alien"
   15   10   19    1.84  396.95 -0.9472E+01  0.7209E+01 "alien"
   34   23    5    1.83  352.42 -0.8409E+01  0.8904E+01 "alien"
   54    5    9    1.36  285.46  0.8123E+02  0.1070E+02 "alien"
   40    1    5    1.88  277.67 -0.6626E+01  0.1184E+02 "alien"
   14    8   19    1.88  245.11 -0.5849E+01  0.5848E+01 "alien"
   43    1    2    1.79  227.11 -0.5419E+01  0.8563E+01 "alien"
   35   23    1    1.85  220.87 -0.5270E+01  0.8630E+01 "alien"
   32   23    9    1.79  215.67 -0.5146E+01  0.8617E+01 "alien"
   30   22   11    1.79  194.51 -0.4641E+01  0.7869E+01 "alien"
   35    1   12    1.82  179.62 -0.4286E+01  0.7883E+01 "alien"
   28   21   12    1.82  155.64 -0.3714E+01  0.7556E+01 "alien"
   41    1    3    1.87  144.92 -0.3458E+01  0.1629E+02 "alien"
   41    1    4    1.85  127.76 -0.3049E+01  0.8818E+01 "alien"
   33   27   13    1.55  127.27 -0.2203E+02  0.1450E+02 "alien"
   33   22    6    1.87  125.40 -0.2992E+01  0.8433E+01 "alien"
   34   22    4    1.88  121.69 -0.2904E+01  0.9967E+01 "alien"
   30   26   15    1.56  104.95 -0.1817E+02  0.1273E+02 "alien"
   37    1   11    1.80   97.77 -0.2333E+01  0.7807E+01 "alien"
   54    5    8    1.37   96.83  0.2755E+02  0.1905E+02 "alien"
   17   15   18    1.82   95.45 -0.2278E+01  0.5427E+01 "alien"
   31   27   15    1.52   92.79 -0.1606E+02  0.1571E+02 "alien"
   34   28   12    1.54   90.19 -0.1561E+02  0.9859E+01 "alien"
   29   20   23    1.34   89.53  0.2548E+02  0.1225E+02 "alien"
   26    9   21    1.54   84.63 -0.1465E+02  0.9335E+01 "alien"
   42    1    4    1.81   81.16 -0.1937E+01  0.8681E+01 "alien"
   31    3   15    1.79   77.81 -0.1857E+01  0.7459E+01 "alien"
   49    2    7    1.52   77.32 -0.1338E+02  0.9678E+01 "alien"
   29   20   10    1.91   75.03  0.2079E+02  0.1459E+02 "alien"
   53    6   11    1.34   70.46  0.2005E+02  0.1038E+02 "alien"
   47    2    8    1.56   68.00 -0.1177E+02  0.8809E+01 "alien"
   41   10   20    1.33   64.20  0.1827E+02  0.9197E+01 "alien"
   41    9   19    1.37   63.57  0.1809E+02  0.1470E+02 "alien"
   50    2    5    1.52   61.10 -0.1058E+02  0.8999E+01 "alien"
   23   17   21    1.52   60.57 -0.1049E+02  0.1052E+02 "alien"
   44    3   12    1.54   58.81 -0.1018E+02  0.8903E+01 "alien"
   38    1    6    1.94   57.49  0.1593E+02  0.9689E+01 "alien"
   26   19   12    1.93   53.47  0.1482E+02  0.8045E+01 "alien"
   21   19   16    1.81   51.79 -0.1236E+01  0.9434E+01 "alien"
   42    4   14    1.53   39.17 -0.6780E+01  0.8418E+01 "alien"
   35    6   18    1.53   37.64 -0.6515E+01  0.1368E+02 "alien"
   35   23    0    1.85   37.11 -0.1758E+02  0.1194E+02 "alien"
   40   29    1    1.57   35.27 -0.6104E+01  0.8897E+01 "alien"
   42    3   13    1.56   35.07 -0.6071E+01  0.1552E+02 "alien"
   32   12   23    1.35   34.61  0.9849E+01  0.8741E+01 "alien"
   38    0    9    1.84   33.98 -0.1186E+02  0.8010E+01 "alien"
   38   29    8    1.54   33.77 -0.5845E+01  0.9052E+01 "alien"
   39   29    4    1.57   33.18  0.2192E+02  0.1021E+02 "alien"
   24   10   22    1.52   33.01 -0.5713E+01  0.7917E+01 "alien"
   46    2    9    1.57   31.73 -0.5493E+01  0.8675E+01 "alien"
   14   10   18    1.94   30.75  0.8520E+01  0.6729E+01 "alien"
   51    2    2    1.51   30.13 -0.5216E+01  0.1108E+02 "alien"
   39    2    3    1.96   28.78  0.7975E+01  0.1033E+02 "alien"
   25   18   12    1.98   28.31  0.7844E+01  0.1251E+02 "alien"
   36   24    0    1.79   27.71 -0.5371E+01  0.1201E+02 "alien"
   50    1    1    1.55   27.35 -0.4734E+01  0.9063E+01 "alien"
   26   25   17    1.56   26.82 -0.4643E+01  0.8919E+01 "alien"
   32   27   14    1.54   26.40 -0.4570E+01  0.1332E+02 "alien"
   40   30    2    1.54   26.37 -0.4565E+01  0.1164E+02 "alien"
   35    5   17    1.58   25.97  0.1716E+02  0.1243E+02 "alien"
   39    2    2    1.97   25.29  0.7008E+01  0.8964E+01 "alien"
   45    8   17    1.36   24.65  0.7013E+01  0.1084E+02 "alien"
   38   28    6    1.59   23.28  0.1538E+02  0.1283E+02 "alien"
   37   29    9    1.54   22.79 -0.3945E+01  0.9068E+01 "alien"
   52    5   11    1.37   22.16  0.6305E+01  0.1563E+02 "alien"
   39    2    0    1.98   20.08  0.2336E+02  0.1300E+02 "alien"
   34   22    2    1.90   19.65  0.5444E+01  0.8785E+01 "alien"
   37   27    7    1.62   18.84  0.1245E+02  0.8997E+01 "alien"
   47    1    5    1.61   17.96  0.1187E+02  0.9119E+01 "alien"
   28   19   23    1.36   17.16  0.4882E+01  0.7876E+01 "alien"
   30    1   13    1.95   17.07  0.4731E+01  0.1208E+02 "alien"
   48    2    7    1.55   16.69 -0.2888E+01  0.9007E+01 "alien"
   38    4   15    1.59   16.30  0.1077E+02  0.8069E+01 "alien"
   39   30    7    1.51   16.19 -0.2802E+01  0.9143E+01 "alien"
   35   23    3    1.83   15.99 -0.3817E+00  0.8676E+01 "alien"
   45    2    9    1.60   15.85  0.1047E+02  0.9820E+01 "alien"
   34    3   15    1.70   15.75  0.1523E+02  0.1245E+02 "alien"
   24   24   18    1.57   15.63 -0.8022E+01  0.7914E+01 "alien"
   55    3    0    1.41   15.52  0.4988E+03  0.3491E+02 "alien"
   32   26   13    1.59   15.10  0.9978E+01  0.8741E+01 "alien"
   40    2    1    1.93   14.97  0.4149E+01  0.1002E+02 "alien"
   39   30    6    1.53   14.64 -0.2534E+01  0.1119E+02 "alien"
   45    3   11    1.54   14.59 -0.2525E+01  0.8876E+01 "alien"
   37   26    5    1.67   14.14 -0.1442E+02  0.8817E+01 "alien"
   53    7   13    1.30   14.01  0.3300E+02  0.1423E+02 "alien"
   39   28    0    1.61   13.33  0.2215E+02  0.1333E+02 "alien"
   27   19   11    1.95   13.28  0.3679E+01  0.8651E+01 "alien"
   38    1   11    1.76   12.60  0.1218E+02  0.8027E+01 "alien"
   34   25    9    1.69   12.51 -0.1276E+02  0.8173E+01 "alien"
   40    2    0    1.93   12.40  0.8105E+01  0.1346E+02 "alien"
   43    1    8    1.69   12.18 -0.1243E+02  0.1465E+02 "alien"
   38   31   10    1.46   11.69 -0.1940E+02  0.9263E+01 "alien"
   36    3   14    1.69   11.29 -0.1152E+02  0.1180E+02 "alien"
   25   25   18    1.54   11.27 -0.4270E+01  0.7743E+01 "alien"
   26   22   15    1.71   11.15  0.1078E+02  0.7749E+01 "alien"
   28   24   15    1.63   11.12  0.7349E+01  0.8031E+01 "alien"
   34   24    8    1.74   11.02  0.1066E+02  0.8410E+01 "alien"
   40   29    0    1.57   10.65  0.1479E+02  0.1316E+02 "alien"
   30   28   16    1.49   10.61 -0.1762E+02  0.1248E+02 "alien"
   29    1   13    1.99   10.50  0.2909E+01  0.7833E+01 "alien"
   51    2    3    1.51   10.40 -0.1726E+02  0.9270E+01 "alien"
   44    0    6    1.70   10.39 -0.2192E+02  0.8578E+01 "alien"
   39   31    8    1.48   10.39 -0.1725E+02  0.9040E+01 "alien"
   16    8   20    1.77   10.10  0.9770E+01  0.6242E+01 "alien"
   35   25    8    1.69   10.07 -0.1027E+02  0.1397E+02 "alien"


 List of    108 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   14.01  0.3300E+02  0.1423E+02
   41   10   20    1.33   64.20  0.1827E+02  0.9197E+01
   29   20   23    1.34   89.53  0.2548E+02  0.1225E+02
   53    6   11    1.34   70.46  0.2005E+02  0.1038E+02
   32   12   23    1.35   34.61  0.9849E+01  0.8741E+01
   54    5    9    1.36  285.46  0.8123E+02  0.1070E+02
   28   19   23    1.36   17.16  0.4882E+01  0.7876E+01
   45    8   17    1.36   24.65  0.7013E+01  0.1084E+02
   52    5   11    1.37   22.16  0.6305E+01  0.1563E+02
   41    9   19    1.37   63.57  0.1809E+02  0.1470E+02
   54    5    8    1.37   96.83  0.2755E+02  0.1905E+02
   55    3    0    1.41   15.52  0.4988E+03  0.3491E+02
   38   31   10    1.46   11.69 -0.1940E+02  0.9263E+01
   39   31    8    1.48   10.39 -0.1725E+02  0.9040E+01
   30   28   16    1.49   10.61 -0.1762E+02  0.1248E+02
   51    2    3    1.51   10.40 -0.1726E+02  0.9270E+01
   51    2    2    1.51   30.13 -0.5216E+01  0.1108E+02
   39   30    7    1.51   16.19 -0.2802E+01  0.9143E+01
   24   10   22    1.52   33.01 -0.5713E+01  0.7917E+01
   50    2    5    1.52   61.10 -0.1058E+02  0.8999E+01
   49    2    7    1.52   77.32 -0.1338E+02  0.9678E+01
   31   27   15    1.52   92.79 -0.1606E+02  0.1571E+02
   23   17   21    1.52   60.57 -0.1049E+02  0.1052E+02
   39   30    6    1.53   14.64 -0.2534E+01  0.1119E+02
   42    4   14    1.53   39.17 -0.6780E+01  0.8418E+01
   35    6   18    1.53   37.64 -0.6515E+01  0.1368E+02
   25   25   18    1.54   11.27 -0.4270E+01  0.7743E+01
   38   29    8    1.54   33.77 -0.5845E+01  0.9052E+01
   32   27   14    1.54   26.40 -0.4570E+01  0.1332E+02
   37   29    9    1.54   22.79 -0.3945E+01  0.9068E+01
   44    3   12    1.54   58.81 -0.1018E+02  0.8903E+01
   26    9   21    1.54   84.63 -0.1465E+02  0.9335E+01
   45    3   11    1.54   14.59 -0.2525E+01  0.8876E+01
   40   30    2    1.54   26.37 -0.4565E+01  0.1164E+02
   34   28   12    1.54   90.19 -0.1561E+02  0.9859E+01
   50    1    1    1.55   27.35 -0.4734E+01  0.9063E+01
   48    2    7    1.55   16.69 -0.2888E+01  0.9007E+01
   33   27   13    1.55  127.27 -0.2203E+02  0.1450E+02
   30   26   15    1.56  104.95 -0.1817E+02  0.1273E+02
   47    2    8    1.56   68.00 -0.1177E+02  0.8809E+01
   26   25   17    1.56   26.82 -0.4643E+01  0.8919E+01
   42    3   13    1.56   35.07 -0.6071E+01  0.1552E+02
   24   24   18    1.57   15.63 -0.8022E+01  0.7914E+01
   40   29    1    1.57   35.27 -0.6104E+01  0.8897E+01
   46    2    9    1.57   31.73 -0.5493E+01  0.8675E+01
   40   29    0    1.57   10.65  0.1479E+02  0.1316E+02
   39   29    4    1.57   33.18  0.2192E+02  0.1021E+02
   35    5   17    1.58   25.97  0.1716E+02  0.1243E+02
   32   26   13    1.59   15.10  0.9978E+01  0.8741E+01
   38   28    6    1.59   23.28  0.1538E+02  0.1283E+02
   38    4   15    1.59   16.30  0.1077E+02  0.8069E+01
   45    2    9    1.60   15.85  0.1047E+02  0.9820E+01
   47    1    5    1.61   17.96  0.1187E+02  0.9119E+01
   39   28    0    1.61   13.33  0.2215E+02  0.1333E+02
   37   27    7    1.62   18.84  0.1245E+02  0.8997E+01
   28   24   15    1.63   11.12  0.7349E+01  0.8031E+01
   37   26    5    1.67   14.14 -0.1442E+02  0.8817E+01
   43    1    8    1.69   12.18 -0.1243E+02  0.1465E+02
   35   25    8    1.69   10.07 -0.1027E+02  0.1397E+02
   34   25    9    1.69   12.51 -0.1276E+02  0.8173E+01
   36    3   14    1.69   11.29 -0.1152E+02  0.1180E+02
   44    0    6    1.70   10.39 -0.2192E+02  0.8578E+01
   34    3   15    1.70   15.75  0.1523E+02  0.1245E+02
   26   22   15    1.71   11.15  0.1078E+02  0.7749E+01
   34   24    8    1.74   11.02  0.1066E+02  0.8410E+01
   38    1   11    1.76   12.60  0.1218E+02  0.8027E+01
   16    8   20    1.77   10.10  0.9770E+01  0.6242E+01
   32   23    9    1.79  215.67 -0.5146E+01  0.8617E+01
   36   24    1    1.79  775.33 -0.1850E+02  0.8683E+01
   36   24    0    1.79   27.71 -0.5371E+01  0.1201E+02
   31    3   15    1.79   77.81 -0.1857E+01  0.7459E+01
   30   22   11    1.79  194.51 -0.4641E+01  0.7869E+01
   43    1    2    1.79  227.11 -0.5419E+01  0.8563E+01
   37    1   11    1.80   97.77 -0.2333E+01  0.7807E+01
   42    1    4    1.81   81.16 -0.1937E+01  0.8681E+01
   21   19   16    1.81   51.79 -0.1236E+01  0.9434E+01
   17   15   18    1.82   95.45 -0.2278E+01  0.5427E+01
   35    1   12    1.82  179.62 -0.4286E+01  0.7883E+01
   42    1    3    1.82  531.00 -0.1267E+02  0.8865E+01
   28   21   12    1.82  155.64 -0.3714E+01  0.7556E+01
   35   23    3    1.83   15.99 -0.3817E+00  0.8676E+01
   34   23    5    1.83  352.42 -0.8409E+01  0.8904E+01
   41    1    5    1.83  487.30 -0.1163E+02  0.8934E+01
   38    0    9    1.84   33.98 -0.1186E+02  0.8010E+01
   15   10   19    1.84  396.95 -0.9472E+01  0.7209E+01
   35   23    1    1.85  220.87 -0.5270E+01  0.8630E+01
   35   23    0    1.85   37.11 -0.1758E+02  0.1194E+02
   41    1    4    1.85  127.76 -0.3049E+01  0.8818E+01
   30   21   10    1.86  563.65 -0.1345E+02  0.1351E+02
   41    1    3    1.87  144.92 -0.3458E+01  0.1629E+02
   33   22    6    1.87  125.40 -0.2992E+01  0.8433E+01
   34   22    4    1.88  121.69 -0.2904E+01  0.9967E+01
   14    8   19    1.88  245.11 -0.5849E+01  0.5848E+01
   40    1    5    1.88  277.67 -0.6626E+01  0.1184E+02
   34   22    2    1.90   19.65  0.5444E+01  0.8785E+01
   29   20   10    1.91   75.03  0.2079E+02  0.1459E+02
   26   19   12    1.93   53.47  0.1482E+02  0.8045E+01
   40    2    1    1.93   14.97  0.4149E+01  0.1002E+02
   40    2    0    1.93   12.40  0.8105E+01  0.1346E+02
   14   10   18    1.94   30.75  0.8520E+01  0.6729E+01
   38    1    6    1.94   57.49  0.1593E+02  0.9689E+01
   30    1   13    1.95   17.07  0.4731E+01  0.1208E+02
   27   19   11    1.95   13.28  0.3679E+01  0.8651E+01
   39    2    3    1.96   28.78  0.7975E+01  0.1033E+02
   39    2    2    1.97   25.29  0.7008E+01  0.8964E+01
   39    2    0    1.98   20.08  0.2336E+02  0.1300E+02
   25   18   12    1.98   28.31  0.7844E+01  0.1251E+02
   29    1   13    1.99   10.50  0.2909E+01  0.7833E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        69
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1407
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1264
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      143


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
