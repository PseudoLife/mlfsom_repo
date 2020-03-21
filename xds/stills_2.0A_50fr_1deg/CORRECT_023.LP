
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.480    77.480    38.773  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      23      23
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      23
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000886  0.000690
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
 REFLECTING_RANGE_E.S.D.=   0.119
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

    1385 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      23 ...      23
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1385 REFLECTIONS ACCEPTED



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
 *  25        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.7   38.8   77.5  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.7   38.8   77.5  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    28        mC        250.0      38.8  159.7   77.5  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  159.7   77.5  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    36        oC        250.0      38.8  159.7   77.5  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.7   77.5  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    27        mC        500.0     159.7   38.8  109.6  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.6  116.2  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  159.7  159.7  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.6  116.2   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    30        mC        998.3      77.5  173.3   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.3      77.5  173.3   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    41        mC        998.3     173.3   77.5   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    22        hP        998.3      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.2  116.2  116.2  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.6   86.6  116.2 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.6   86.6  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.6   86.6  116.2 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.6   86.6   86.6  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.6   86.6  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.6   86.6  109.6  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.6  248.1 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.6   86.6   77.5  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.6   86.6   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.5  177.5  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.6   86.6  177.5 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.6   86.6   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.6   86.6   86.6  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.5  109.6   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
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

  *   0.00        0    5.5      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.844   -23.456     0.057
 UNIT_CELL_B-AXIS=    23.456    73.844     0.020
 UNIT_CELL_C-AXIS=    -0.030    -0.001    38.773
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       236 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.30
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.326    77.326    38.700  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.3E-01 1.3E-01 6.2E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012932  0.012932  0.025840  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.692   -23.428     0.052
 COORDINATES OF UNIT CELL B-AXIS    23.428    73.692     0.018
 COORDINATES OF UNIT CELL C-AXIS    -0.028    -0.001    38.700
 CRYSTAL MOSAICITY (DEGREES)     0.119
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006473  0.001722
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.91   1232.91
 DETECTOR ORIGIN (PIXELS) AT                     1261.41   1242.31
 CRYSTAL TO DETECTOR DISTANCE (mm)       249.88
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.33    77.33    38.70  90.000  90.000  90.000


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
  0.05  0.15        81        137.712       1.110
  0.15  0.25       172        174.218       0.877
  0.25  0.35       188        127.786       1.196
  0.35  0.45       146         94.430       1.618
  0.45  0.55       180        195.071       0.783
  0.55  0.65       184        145.400       1.051
  0.65  0.75       203        161.065       0.949
  0.75  0.85       153        166.999       0.915
  0.85  0.95        78        161.275       0.947
 
 -0.05  0.95      1385        152.797       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               250.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.854
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    22.1 XMAX=    22.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59643 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        516


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               250.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.848
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     9.3 XMAX=  2455.0 NXBIN=    3
 YMIN=     7.9 YMAX=  2519.9 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        516


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               250.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.836
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    22.1 XMAX=    22.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1648    1691
 DETECTOR_SURFACE_POSITION=     817    1691
 DETECTOR_SURFACE_POSITION=     817     837
 DETECTOR_SURFACE_POSITION=    1648     837
 DETECTOR_SURFACE_POSITION=    2172    1664
 DETECTOR_SURFACE_POSITION=    1622    2229
 DETECTOR_SURFACE_POSITION=     843    2229
 DETECTOR_SURFACE_POSITION=     292    1664
 DETECTOR_SURFACE_POSITION=     292     864
 DETECTOR_SURFACE_POSITION=     843     298
 DETECTOR_SURFACE_POSITION=    1622     298
 DETECTOR_SURFACE_POSITION=    2172     864
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        516


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
 1.049E+00  4.373E-03   14.76


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

   24.453   3.150    12.88   1.36      3.59      3.95      88     110      14
    3.150   2.236     7.31   1.12      6.41      6.10     200     228       0
    2.236   1.829     1.40   0.89     30.18     31.32     208     241       0
    1.829   1.585    -0.06   0.58    -99.90    -99.90     246     284       0
    1.585   1.418     0.02   0.60    682.20    910.26     282     330       0
    1.418   1.295    -0.00   0.64   2489.95   3085.34     134     175       2
   --------------------------------------------------------------------------
   24.453   1.295     2.49   0.80      7.57      8.51    1158    1368      16
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      23      23)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681   0.9264E+01  0.1340E+02         0.69           1 
   47    0    0     1.645   0.6022E+01  0.1327E+02         0.45           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     3.7% OF MEAN INTENSITY


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

     3.84          68      41      1285        3.2%       3.1%      3.9%       54   16.83      4.3%    99.8*     0    0.000       0
     2.73         106      61      2098        2.9%       5.5%      4.5%       90   13.51      7.8%    99.2*     0    0.000       0
     2.24         164      92      2665        3.5%       7.1%      6.9%      144    8.92     10.1%    98.7*     0    0.000       0
     1.94         156      88      3101        2.8%      20.9%     21.4%      136    2.86     29.5%    93.5*     0    0.000       0
     1.74         167      95      3481        2.7%     -99.9%    -99.9%      144    0.00    -99.9%    -7.8      0    0.000       0
     1.58         202     115      3848        3.0%     -99.9%    -99.9%      174    0.00    -99.9%    16.7      0    0.000       0
     1.47         212     123      4143        3.0%     497.8%    640.8%      178    0.04    703.9%    -3.0      0    0.000       0
     1.37         210     117      4429        2.6%    1909.2%   2455.0%      186   -0.03   2700.0%     5.0      0    0.000       0
     1.29          83      57      4743        1.2%    8028.7%  11323.7%       52    0.01  11354.2%    24.5      0    0.000       0
    total        1368     789     29793        2.6%       7.6%      8.5%     1158    3.25     10.7%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1385
 NUMBER OF REJECTED MISFITS                                16
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1368
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    789



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

   WILSON LINE (using all data) : A=   7.533 B=  33.799 CORRELATION=  0.83
      #      RES      SS        <I>       log(<I>)       BO
      19     7.290   0.005  1.7965E+03       7.494       4.2
      18     4.694   0.011  2.1748E+03       7.685      -6.7
      26     3.515   0.020  9.2697E+02       6.832      17.3
      31     2.993   0.028  3.0321E+02       5.714      32.6
      37     2.640   0.036  2.7677E+02       5.623      26.6
      49     2.398   0.043  1.4712E+02       4.991      29.2
      38     2.207   0.051  7.9158E+01       4.371      30.8
      44     2.054   0.059  2.9220E+01       3.375      35.1
      43     1.930   0.067 -1.6941E+00    -999.999   -1000.0
      50     1.828   0.075 -1.7167E+00    -999.999   -1000.0
      46     1.739   0.083  1.6970E-01      -1.774      56.3
      51     1.665   0.090 -7.3319E-01    -999.999   -1000.0
      52     1.596   0.098 -1.2102E+00    -999.999   -1000.0
      62     1.535   0.106  3.2569E-02      -3.424      51.7
      69     1.479   0.114 -3.2748E-01    -999.999   -1000.0
      53     1.432   0.122  9.1314E-01      -0.091      31.3
      45     1.390   0.129 -6.3433E-01    -999.999   -1000.0
      38     1.350   0.137 -2.3865E+00    -999.999   -1000.0
      18     1.313   0.145  5.0606E+00       1.621      20.4


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.290        2.182        1.625        0.970
       3     4.694        1.442        0.946        0.461
       2     3.515        0.401        0.088        0.014
       4     2.993        1.284        0.613        0.220
       2     2.640        0.289        0.076        0.014
       5     2.398        0.804        0.386        0.143
       2     2.207        1.282        0.572        0.192
       7     2.054        0.599        0.270        0.093
      10     1.828       63.583      -28.898      727.647
       3     1.596      100.895      175.565     1054.500
       1     1.350       21.833       35.339       40.858
      81   overall      843.414   -23003.609  2030818.138


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      14     7.290        0.482        0.372        0.233
      15     4.694        1.197        1.508        1.523
      24     3.515        0.838        0.680        0.522
      27     2.993        0.876        0.775        0.587
      35     2.640        1.151        1.304        1.370
      44     2.398        0.924        0.775        0.552
      36     2.207        0.783        0.579        0.379
      37     2.054        1.271        1.466        1.580
      39     1.739     1263.620   -14082.419  1515229.437
      59     1.535    30580.651  -989698.386503314474.433
      50     1.432       46.649       60.386      746.918
      18     1.313        4.014       12.444       33.758
     708   overall     2621.708   -83245.701 42026392.785

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.290  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.30  0.29
       3     4.694  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.515  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       4     2.993  0.00  0.00  0.00  0.00  0.00  0.00  0.42  0.40  0.38  0.37
       2     2.640  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.398  0.00  0.00  0.00  0.00  0.77  0.71  1.00  0.95  0.91  0.88
       2     2.207  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       7     2.054  1.15  0.83  1.03  0.91  1.10  1.02  0.96  1.14  1.09  1.05
      10     1.828  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
       3     1.596  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       1     1.350  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      81   overall  0.65  0.46  0.41  0.37  0.40  0.37  0.39  0.39  0.39  0.40


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      14     7.290  1.50  2.76  1.93  1.73  2.00  1.74  1.70  1.56  1.44  1.36
      15     4.694  1.40  1.47  1.03  1.01  1.02  1.63  1.46  1.45  1.46  1.37
      24     3.515  0.88  0.69  0.64  0.51  0.64  0.74  1.08  0.98  1.05  0.99
      27     2.993  1.17  1.84  1.43  1.35  1.22  1.31  1.25  1.21  1.19  1.11
      35     2.640  0.90  0.95  1.21  1.13  0.94  0.95  1.08  1.09  1.16  1.08
      44     2.398  1.19  0.88  0.88  1.10  1.10  1.11  1.04  1.07  1.03  1.04
      36     2.207  1.17  0.61  0.86  0.84  0.85  0.86  0.94  0.91  0.98  0.97
      37     2.054  2.56  1.49  1.04  1.07  0.96  0.96  0.97  0.88  0.96  0.94
      39     1.739  5.12  2.69  1.98  1.56  1.30  1.14  1.02  0.93  0.86  0.81
      59     1.535  4.81  2.52  1.77  1.39  1.16  1.01  0.91  0.83  0.77  0.72
      50     1.432  4.41  2.32  1.62  1.33  1.12  0.98  0.87  0.80  0.74  0.70
      18     1.313  4.09  2.15  1.71  1.52  1.41  1.23  1.10  1.01  0.94  0.97
     708   overall  1.54  0.97  0.76  0.68  0.62  0.60  0.59  0.56  0.55  0.53


 List of    123 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   51    2    2    1.51  565.78  0.1843E+02  0.1263E+02 "alien"
   36   29   11    1.51  407.76  0.1328E+02  0.1014E+02 "alien"
   33    7   19    1.52  354.16  0.1153E+02  0.7826E+01 "alien"
   40   30    5    1.52  340.60  0.1109E+02  0.9713E+01 "alien"
   48    2    8    1.53  326.78  0.1064E+02  0.1004E+02 "alien"
   20    8   23    1.52  323.09  0.1052E+02  0.8084E+01 "alien"
   50    1    1    1.54  299.87  0.9766E+01  0.9778E+01 "alien"
   39   29    6    1.54  298.82  0.9732E+01  0.1045E+02 "alien"
   26   25   17    1.56  259.12  0.8439E+01  0.9525E+01 "alien"
   40   29    2    1.56  247.93  0.8075E+01  0.1005E+02 "alien"
   22   10   22    1.54  246.13  0.8016E+01  0.8607E+01 "alien"
   32   27   14    1.54  244.24  0.7955E+01  0.1090E+02 "alien"
   34   28   12    1.54  238.88  0.7780E+01  0.1073E+02 "alien"
   31    7   19    1.56  180.03  0.5863E+01  0.1195E+02 "alien"
   40   30    3    1.54  169.09  0.5507E+01  0.9847E+01 "alien"
   47    3   10    1.51  154.55  0.5034E+01  0.9472E+01 "alien"
   36   28   10    1.55  145.50  0.4739E+01  0.9580E+01 "alien"
   28    9   21    1.51  135.84  0.4424E+01  0.7065E+01 "alien"
   47    2    8    1.56  113.33  0.3691E+01  0.9528E+01 "alien"
   24   10   22    1.51  110.03  0.3583E+01  0.8557E+01 "alien"
   39   29    5    1.56  104.58  0.3406E+01  0.9520E+01 "alien"
   34   24    7    1.76  100.11  0.1699E+02  0.9789E+01 "alien"
   46    3   11    1.51   96.49  0.3143E+01  0.9575E+01 "alien"
   30   26   15    1.55   92.58  0.3015E+01  0.1308E+02 "alien"
   38    5   16    1.55   90.90  0.2960E+01  0.8745E+01 "alien"
   16    9   20    1.76   89.51  0.1519E+02  0.7158E+01 "alien"
   37   25    2    1.72   80.90  0.1373E+02  0.1030E+02 "alien"
   21   10   22    1.55   75.00  0.2443E+01  0.8501E+01 "alien"
   22   14   22    1.51   69.07  0.2250E+01  0.9349E+01 "alien"
   34   24    8    1.73   68.91  0.1169E+02  0.8987E+01 "alien"
   36   25    5    1.72   62.18  0.1055E+02  0.8980E+01 "alien"
   49    2    7    1.52   59.77  0.1947E+01  0.1484E+02 "alien"
   35   25    7    1.71   57.28  0.9721E+01  0.9076E+01 "alien"
   34    3   15    1.70   54.01  0.9166E+01  0.1313E+02 "alien"
   46    2    9    1.56   51.39  0.1674E+01  0.9137E+01 "alien"
   19   18   18    1.74   51.38  0.8719E+01  0.6727E+01 "alien"
   41   32    5    1.46   48.36 -0.1584E+02  0.9937E+01 "alien"
   40   31    7    1.47   47.98 -0.1571E+02  0.1063E+02 "alien"
   33   25   11    1.65   47.90 -0.3512E+02  0.1473E+02 "alien"
   43    1    3    1.78   41.77  0.7087E+01  0.9380E+01 "alien"
   23   13   23    1.46   40.14 -0.1314E+02  0.9542E+01 "alien"
   28   22   13    1.75   40.09  0.6803E+01  0.8162E+01 "alien"
   33    8   20    1.47   38.75 -0.1269E+02  0.1253E+02 "alien"
   34    7   19    1.50   35.98 -0.1178E+02  0.1221E+02 "alien"
   38   31   10    1.46   34.41 -0.1127E+02  0.9764E+01 "alien"
   26   22   15    1.70   33.25  0.5642E+01  0.8175E+01 "alien"
   31    9   21    1.46   31.85 -0.1043E+02  0.1237E+02 "alien"
   36   32   14    1.39   31.45 -0.1995E+02  0.1943E+02 "alien"
   39   30    7    1.51   30.19  0.9832E+00  0.9547E+01 "alien"
   42   32    1    1.46   30.00 -0.9824E+01  0.1001E+02 "alien"
   25   25   18    1.53   28.14  0.7738E+01  0.8433E+01 "alien"
   37   31   11    1.46   27.95 -0.9154E+01  0.1174E+02 "alien"
   43    1    8    1.69   27.75 -0.2035E+02  0.1632E+02 "alien"
   30   23   12    1.73   27.19  0.4615E+01  0.8488E+01 "alien"
   53    4    8    1.39   27.09 -0.1718E+02  0.1003E+02 "alien"
   51    2    3    1.50   26.28 -0.8607E+01  0.9974E+01 "alien"
   31   30   17    1.41   26.00 -0.1649E+02  0.9914E+01 "alien"
   50    2    6    1.50   25.90 -0.8481E+01  0.9583E+01 "alien"
   25   21   15    1.74   25.68  0.4358E+01  0.8830E+01 "alien"
   52    2    4    1.47   24.60 -0.8056E+01  0.1265E+02 "alien"
   23   15   22    1.49   21.33 -0.6984E+01  0.9288E+01 "alien"
   39   28    4    1.59   20.79 -0.2516E+02  0.1171E+02 "alien"
   37   25    0    1.73   20.37  0.1486E+02  0.1328E+02 "alien"
   29    9   21    1.49   19.88 -0.6509E+01  0.7038E+01 "alien"
   30   28   16    1.49   19.47 -0.6375E+01  0.1305E+02 "alien"
   26   19   22    1.42   18.92  0.1727E+02  0.8465E+01 "alien"
   37   26    4    1.68   18.82 -0.1380E+02  0.9043E+01 "alien"
   31   29   16    1.46   18.54  0.1693E+02  0.9844E+01 "alien"
   30   12   23    1.38   18.07 -0.1146E+02  0.8697E+01 "alien"
   33   29   14    1.48   17.30 -0.5667E+01  0.1053E+02 "alien"
   49    4   11    1.44   17.15  0.1566E+02  0.1031E+02 "alien"
   50    3    8    1.47   17.14 -0.5612E+01  0.9806E+01 "alien"
   36   30   12    1.47   16.87 -0.5524E+01  0.9875E+01 "alien"
   52    2    2    1.48   16.51 -0.5407E+01  0.9794E+01 "alien"
   49    5   12    1.41   15.64  0.1428E+02  0.1038E+02 "alien"
   52    3    7    1.43   15.24  0.1392E+02  0.1116E+02 "alien"
   38   33   12    1.39   15.16 -0.9618E+01  0.9975E+01 "alien"
   37   27    6    1.63   15.13 -0.1110E+02  0.1326E+02 "alien"
   53    2    1    1.46   15.13 -0.4954E+01  0.1186E+02 "alien"
   32    8   20    1.49   15.05 -0.4930E+01  0.8075E+01 "alien"
   53    4    7    1.41   14.91 -0.9456E+01  0.1071E+02 "alien"
   17    9   21    1.68   14.91 -0.1093E+02  0.6518E+01 "alien"
   18   10   21    1.65   14.65 -0.1074E+02  0.7530E+01 "alien"
   37   30   11    1.47   14.63 -0.4790E+01  0.1147E+02 "alien"
   30    9   21    1.48   14.26 -0.4669E+01  0.7737E+01 "alien"
   46    5   14    1.43   13.97  0.1276E+02  0.1219E+02 "alien"
   38    1   11    1.76   13.79  0.2340E+01  0.8187E+01 "alien"
   37   29    9    1.54   13.75  0.4479E+00  0.9642E+01 "alien"
   39    2    3    1.96   13.66 -0.2314E+02  0.1130E+02 "alien"
   36   31   13    1.43   13.61  0.1243E+02  0.1022E+02 "alien"
   19    9   22    1.59   13.51 -0.1634E+02  0.8068E+01 "alien"
   45    0    3    1.70   13.23  0.5606E+01  0.9801E+01 "alien"
   48    2    7    1.55   13.19  0.4295E+00  0.9591E+01 "alien"
   38   26    0    1.68   13.18 -0.2404E+02  0.1684E+02 "alien"
   42    0    7    1.75   13.04  0.5454E+01  0.8693E+01 "alien"
   55    3    0    1.40   12.93 -0.2006E+02  0.1301E+02 "alien"
   42   32    0    1.46   12.68 -0.9989E+01  0.1415E+02 "alien"
   47    4   12    1.46   12.61 -0.4129E+01  0.9874E+01 "alien"
   48    1    4    1.59   12.33 -0.1492E+02  0.9582E+01 "alien"
   21   19   19    1.63   12.31 -0.9028E+01  0.1132E+02 "alien"
   27   24   16    1.60   12.31 -0.1489E+02  0.1559E+02 "alien"
   28   21   22    1.38   12.18 -0.7723E+01  0.8051E+01 "alien"
   50    4   10    1.43   12.07  0.1102E+02  0.1009E+02 "alien"
   41   31    1    1.50   12.00 -0.3929E+01  0.9927E+01 "alien"
   52    2    0    1.49   11.97 -0.8986E+01  0.1384E+02 "alien"
   21   17   20    1.60   11.81 -0.1430E+02  0.9692E+01 "alien"
   40   33    9    1.41   11.80 -0.7483E+01  0.1021E+02 "alien"
   27   22   14    1.73   11.65  0.1978E+01  0.7778E+01 "alien"
   35   25    8    1.69   11.49 -0.8424E+01  0.1073E+02 "alien"
   37    9   20    1.40   11.32 -0.7180E+01  0.8477E+01 "alien"
   46    0    1    1.68   11.26 -0.1805E+02  0.9289E+01 "alien"
   39   31    9    1.46   11.23 -0.3678E+01  0.1090E+02 "alien"
   39   34   11    1.38   11.02 -0.6989E+01  0.1226E+02 "alien"
   53    7   13    1.30   10.98  0.5555E+02  0.1613E+02 "alien"
   26   21   21    1.44   10.94  0.9993E+01  0.8639E+01 "alien"
   37   33   13    1.38   10.90 -0.6914E+01  0.1011E+02 "alien"
   26   26   18    1.50   10.82 -0.7526E+01  0.8583E+01 "alien"
   46    0    3    1.67   10.75 -0.1667E+02  0.9029E+01 "alien"
   15    8   20    1.78   10.73  0.1821E+01  0.6254E+01 "alien"
   20   18   16    1.85   10.52 -0.1805E+02  0.1008E+02 "alien"
   39   32   10    1.43   10.31  0.9419E+01  0.1050E+02 "alien"
   45    5   14    1.45   10.27  0.9377E+01  0.9514E+01 "alien"
   24   23   17    1.63   10.22 -0.1236E+02  0.8565E+01 "alien"


 List of    123 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   10.98  0.5555E+02  0.1613E+02
   39   34   11    1.38   11.02 -0.6989E+01  0.1226E+02
   28   21   22    1.38   12.18 -0.7723E+01  0.8051E+01
   30   12   23    1.38   18.07 -0.1146E+02  0.8697E+01
   37   33   13    1.38   10.90 -0.6914E+01  0.1011E+02
   38   33   12    1.39   15.16 -0.9618E+01  0.9975E+01
   36   32   14    1.39   31.45 -0.1995E+02  0.1943E+02
   53    4    8    1.39   27.09 -0.1718E+02  0.1003E+02
   37    9   20    1.40   11.32 -0.7180E+01  0.8477E+01
   55    3    0    1.40   12.93 -0.2006E+02  0.1301E+02
   53    4    7    1.41   14.91 -0.9456E+01  0.1071E+02
   31   30   17    1.41   26.00 -0.1649E+02  0.9914E+01
   40   33    9    1.41   11.80 -0.7483E+01  0.1021E+02
   49    5   12    1.41   15.64  0.1428E+02  0.1038E+02
   26   19   22    1.42   18.92  0.1727E+02  0.8465E+01
   39   32   10    1.43   10.31  0.9419E+01  0.1050E+02
   36   31   13    1.43   13.61  0.1243E+02  0.1022E+02
   46    5   14    1.43   13.97  0.1276E+02  0.1219E+02
   50    4   10    1.43   12.07  0.1102E+02  0.1009E+02
   52    3    7    1.43   15.24  0.1392E+02  0.1116E+02
   49    4   11    1.44   17.15  0.1566E+02  0.1031E+02
   26   21   21    1.44   10.94  0.9993E+01  0.8639E+01
   45    5   14    1.45   10.27  0.9377E+01  0.9514E+01
   31   29   16    1.46   18.54  0.1693E+02  0.9844E+01
   53    2    1    1.46   15.13 -0.4954E+01  0.1186E+02
   37   31   11    1.46   27.95 -0.9154E+01  0.1174E+02
   23   13   23    1.46   40.14 -0.1314E+02  0.9542E+01
   39   31    9    1.46   11.23 -0.3678E+01  0.1090E+02
   41   32    5    1.46   48.36 -0.1584E+02  0.9937E+01
   38   31   10    1.46   34.41 -0.1127E+02  0.9764E+01
   31    9   21    1.46   31.85 -0.1043E+02  0.1237E+02
   47    4   12    1.46   12.61 -0.4129E+01  0.9874E+01
   42   32    1    1.46   30.00 -0.9824E+01  0.1001E+02
   42   32    0    1.46   12.68 -0.9989E+01  0.1415E+02
   52    2    4    1.47   24.60 -0.8056E+01  0.1265E+02
   36   30   12    1.47   16.87 -0.5524E+01  0.9875E+01
   50    3    8    1.47   17.14 -0.5612E+01  0.9806E+01
   40   31    7    1.47   47.98 -0.1571E+02  0.1063E+02
   37   30   11    1.47   14.63 -0.4790E+01  0.1147E+02
   33    8   20    1.47   38.75 -0.1269E+02  0.1253E+02
   30    9   21    1.48   14.26 -0.4669E+01  0.7737E+01
   52    2    2    1.48   16.51 -0.5407E+01  0.9794E+01
   33   29   14    1.48   17.30 -0.5667E+01  0.1053E+02
   52    2    0    1.49   11.97 -0.8986E+01  0.1384E+02
   30   28   16    1.49   19.47 -0.6375E+01  0.1305E+02
   23   15   22    1.49   21.33 -0.6984E+01  0.9288E+01
   32    8   20    1.49   15.05 -0.4930E+01  0.8075E+01
   29    9   21    1.49   19.88 -0.6509E+01  0.7038E+01
   50    2    6    1.50   25.90 -0.8481E+01  0.9583E+01
   34    7   19    1.50   35.98 -0.1178E+02  0.1221E+02
   41   31    1    1.50   12.00 -0.3929E+01  0.9927E+01
   26   26   18    1.50   10.82 -0.7526E+01  0.8583E+01
   51    2    3    1.50   26.28 -0.8607E+01  0.9974E+01
   28    9   21    1.51  135.84  0.4424E+01  0.7065E+01
   51    2    2    1.51  565.78  0.1843E+02  0.1263E+02
   36   29   11    1.51  407.76  0.1328E+02  0.1014E+02
   47    3   10    1.51  154.55  0.5034E+01  0.9472E+01
   39   30    7    1.51   30.19  0.9832E+00  0.9547E+01
   22   14   22    1.51   69.07  0.2250E+01  0.9349E+01
   24   10   22    1.51  110.03  0.3583E+01  0.8557E+01
   46    3   11    1.51   96.49  0.3143E+01  0.9575E+01
   49    2    7    1.52   59.77  0.1947E+01  0.1484E+02
   40   30    5    1.52  340.60  0.1109E+02  0.9713E+01
   33    7   19    1.52  354.16  0.1153E+02  0.7826E+01
   20    8   23    1.52  323.09  0.1052E+02  0.8084E+01
   48    2    8    1.53  326.78  0.1064E+02  0.1004E+02
   25   25   18    1.53   28.14  0.7738E+01  0.8433E+01
   40   30    3    1.54  169.09  0.5507E+01  0.9847E+01
   32   27   14    1.54  244.24  0.7955E+01  0.1090E+02
   37   29    9    1.54   13.75  0.4479E+00  0.9642E+01
   22   10   22    1.54  246.13  0.8016E+01  0.8607E+01
   34   28   12    1.54  238.88  0.7780E+01  0.1073E+02
   39   29    6    1.54  298.82  0.9732E+01  0.1045E+02
   50    1    1    1.54  299.87  0.9766E+01  0.9778E+01
   48    2    7    1.55   13.19  0.4295E+00  0.9591E+01
   38    5   16    1.55   90.90  0.2960E+01  0.8745E+01
   36   28   10    1.55  145.50  0.4739E+01  0.9580E+01
   30   26   15    1.55   92.58  0.3015E+01  0.1308E+02
   21   10   22    1.55   75.00  0.2443E+01  0.8501E+01
   47    2    8    1.56  113.33  0.3691E+01  0.9528E+01
   39   29    5    1.56  104.58  0.3406E+01  0.9520E+01
   40   29    2    1.56  247.93  0.8075E+01  0.1005E+02
   26   25   17    1.56  259.12  0.8439E+01  0.9525E+01
   31    7   19    1.56  180.03  0.5863E+01  0.1195E+02
   46    2    9    1.56   51.39  0.1674E+01  0.9137E+01
   19    9   22    1.59   13.51 -0.1634E+02  0.8068E+01
   39   28    4    1.59   20.79 -0.2516E+02  0.1171E+02
   48    1    4    1.59   12.33 -0.1492E+02  0.9582E+01
   27   24   16    1.60   12.31 -0.1489E+02  0.1559E+02
   21   17   20    1.60   11.81 -0.1430E+02  0.9692E+01
   24   23   17    1.63   10.22 -0.1236E+02  0.8565E+01
   21   19   19    1.63   12.31 -0.9028E+01  0.1132E+02
   37   27    6    1.63   15.13 -0.1110E+02  0.1326E+02
   33   25   11    1.65   47.90 -0.3512E+02  0.1473E+02
   18   10   21    1.65   14.65 -0.1074E+02  0.7530E+01
   46    0    3    1.67   10.75 -0.1667E+02  0.9029E+01
   17    9   21    1.68   14.91 -0.1093E+02  0.6518E+01
   38   26    0    1.68   13.18 -0.2404E+02  0.1684E+02
   46    0    1    1.68   11.26 -0.1805E+02  0.9289E+01
   37   26    4    1.68   18.82 -0.1380E+02  0.9043E+01
   43    1    8    1.69   27.75 -0.2035E+02  0.1632E+02
   35   25    8    1.69   11.49 -0.8424E+01  0.1073E+02
   34    3   15    1.70   54.01  0.9166E+01  0.1313E+02
   45    0    3    1.70   13.23  0.5606E+01  0.9801E+01
   26   22   15    1.70   33.25  0.5642E+01  0.8175E+01
   35   25    7    1.71   57.28  0.9721E+01  0.9076E+01
   36   25    5    1.72   62.18  0.1055E+02  0.8980E+01
   37   25    2    1.72   80.90  0.1373E+02  0.1030E+02
   30   23   12    1.73   27.19  0.4615E+01  0.8488E+01
   27   22   14    1.73   11.65  0.1978E+01  0.7778E+01
   37   25    0    1.73   20.37  0.1486E+02  0.1328E+02
   34   24    8    1.73   68.91  0.1169E+02  0.8987E+01
   19   18   18    1.74   51.38  0.8719E+01  0.6727E+01
   25   21   15    1.74   25.68  0.4358E+01  0.8830E+01
   42    0    7    1.75   13.04  0.5454E+01  0.8693E+01
   28   22   13    1.75   40.09  0.6803E+01  0.8162E+01
   16    9   20    1.76   89.51  0.1519E+02  0.7158E+01
   38    1   11    1.76   13.79  0.2340E+01  0.8187E+01
   34   24    7    1.76  100.11  0.1699E+02  0.9789E+01
   15    8   20    1.78   10.73  0.1821E+01  0.6254E+01
   43    1    3    1.78   41.77  0.7087E+01  0.9380E+01
   20   18   16    1.85   10.52 -0.1805E+02  0.1008E+02
   39    2    3    1.96   13.66 -0.2314E+02  0.1130E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        63
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1385
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1262
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      123


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
