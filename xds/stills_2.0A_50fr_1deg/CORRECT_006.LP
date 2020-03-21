
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.317    77.317    38.692  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=       6       6
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       6
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000642  0.000712
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
 REFLECTING_RANGE_E.S.D.=   0.085
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

    1508 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       6 ...       6
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1508 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               47
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
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  25        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  23        oC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  20        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    27        mC        500.0     159.4   38.7  109.3  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.7  109.3  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.3  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        998.3     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.3      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.3      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.3      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.5   86.5  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.3   86.5   86.5  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.5   86.5  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.5   86.5  109.3  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.7   86.5  247.5 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.5   86.5   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.2  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.5   86.5  177.2 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.2  109.3   38.7  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
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

  *   0.00        0    6.1      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.686   -23.417     0.036
 UNIT_CELL_B-AXIS=    23.417    73.686     0.012
 UNIT_CELL_C-AXIS=    -0.019    -0.000    38.692
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       338 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.17
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.263    77.263    38.666  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  8.5E-02 8.5E-02 5.0E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012943  0.012943  0.025863  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.634   -23.401     0.036
 COORDINATES OF UNIT CELL B-AXIS    23.401    73.634     0.011
 COORDINATES OF UNIT CELL C-AXIS    -0.019     0.000    38.666
 CRYSTAL MOSAICITY (DEGREES)     0.085
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006406  0.002207
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.86   1232.90
 DETECTOR ORIGIN (PIXELS) AT                     1260.65   1242.21
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.10
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.26    77.26    38.67  90.000  90.000  90.000


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
  0.05  0.15       123        309.353       0.933
  0.15  0.25       182        242.689       1.190
  0.25  0.35       186        270.147       1.069
  0.35  0.45       157        163.938       1.761
  0.45  0.55       176        347.314       0.831
  0.55  0.65       199        311.798       0.926
  0.65  0.75       181        189.999       1.520
  0.75  0.85       169        282.531       1.022
  0.85  0.95       135        532.943       0.542
 
 -0.05  0.95      1508        288.770       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               236.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.031
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     5.1 XMAX=     5.9 NXBIN=    1
 YMIN= 0.00168 YMAX= 0.59744 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               236.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.028
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    12.2 XMAX=  2453.8 NXBIN=    3
 YMIN=     6.7 YMAX=  2519.5 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               236.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.022
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     5.1 XMAX=     5.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1263
 DETECTOR_SURFACE_POSITION=    1648    1690
 DETECTOR_SURFACE_POSITION=     818    1690
 DETECTOR_SURFACE_POSITION=     818     836
 DETECTOR_SURFACE_POSITION=    1648     836
 DETECTOR_SURFACE_POSITION=    2171    1663
 DETECTOR_SURFACE_POSITION=    1622    2229
 DETECTOR_SURFACE_POSITION=     844    2229
 DETECTOR_SURFACE_POSITION=     295    1663
 DETECTOR_SURFACE_POSITION=     295     863
 DETECTOR_SURFACE_POSITION=     844     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2171     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


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
 1.412E+00  2.537E-03   16.71


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

   24.433   3.395    15.08   1.33      2.70      3.46      64      90      20
    3.395   2.412    12.30   1.33      4.95      3.90     126     167      16
    2.412   1.973     8.26   1.14      5.76      5.19     182     220       6
    1.973   1.710    -0.05   0.51    -99.90    -99.90     192     246       0
    1.710   1.530    -0.00   0.61    -99.90    -99.90     214     282       0
    1.530   1.397     0.04   0.63    801.50   1000.69     244     312       0
    1.397   1.294     0.05   0.68    647.74    824.04     102     148       0
   --------------------------------------------------------------------------
   24.433   1.294     3.57   0.81      5.58      5.96    1124    1465      42
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       6       6)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.680   0.3831E+01  0.1330E+02         0.29           1 
   47    0    0     1.644  -0.8361E+01  0.1328E+02        -0.63           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -2.8% OF MEAN INTENSITY


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

     3.84          68      44      1285        3.4%       2.6%      3.4%       48   18.74      3.6%    99.9*     0    0.000       0
     2.73         107      68      2098        3.2%       4.3%      3.7%       78   16.49      6.1%    98.6*     0    0.000       0
     2.23         160      97      2665        3.6%       5.3%      4.1%      126   14.33      7.5%    98.6*     0    0.000       0
     1.94         167      97      3101        3.1%       6.3%      6.6%      140    7.34      8.9%    98.9*     0    0.000       0
     1.73         190     116      3479        3.3%     -99.9%    -99.9%      148    0.00    -99.9%   -22.9      0    0.000       0
     1.58         219     134      3850        3.5%     -99.9%    -99.9%      170    0.00    -99.9%     4.7      0    0.000       0
     1.47         243     152      4143        3.7%     316.8%    391.0%      182    0.11    448.0%   -21.9      0    0.000       0
     1.37         221     132      4428        3.0%     -99.9%    -99.9%      178    0.00    -99.9%     4.9      0    0.000       0
     1.29          90      63      4743        1.3%     329.0%    405.4%       54    0.11    465.3%   -25.5      0    0.000       0
    total        1465     903     29792        3.0%       5.6%      6.0%     1124    4.49      7.9%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1508
 NUMBER OF REJECTED MISFITS                                42
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1465
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    903



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

   WILSON LINE (using all data) : A=   8.543 B=  32.550 CORRELATION=  0.95
      #      RES      SS        <I>       log(<I>)       BO
      18     7.516   0.004  3.0954E+03       8.038      57.1
      22     4.644   0.012  3.1456E+03       8.054      21.1
      30     3.492   0.020  1.4472E+03       7.277      30.9
      32     2.993   0.028  7.3138E+02       6.595      34.9
      39     2.643   0.036  8.2001E+02       6.709      25.6
      54     2.392   0.044  4.8941E+02       6.193      26.9
      41     2.202   0.052  3.6319E+02       5.895      25.7
      47     2.055   0.059  1.7040E+02       5.138      28.8
      50     1.930   0.067 -4.2199E-01    -999.999   -1000.0
      63     1.825   0.075 -5.6974E-01    -999.999   -1000.0
      54     1.738   0.083 -3.6799E-01    -999.999   -1000.0
      60     1.663   0.090  1.0521E+00       0.051      47.0
      62     1.595   0.098 -8.8213E-01    -999.999   -1000.0
      78     1.533   0.106  1.4930E+00       0.401      38.3
      80     1.479   0.114  1.2651E+00       0.235      36.3
      58     1.431   0.122 -6.8667E-01    -999.999   -1000.0
      54     1.390   0.129 -3.4927E-01    -999.999   -1000.0
      41     1.350   0.137  2.2094E+00       0.793      28.3
      20     1.312   0.145  1.6108E+00       0.477      27.8


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.516        0.960        0.450        0.167
       3     4.644        1.071        0.599        0.249
       3     3.492        0.367        0.094        0.017
       4     2.993        1.060        0.462        0.152
       3     2.643        0.232        0.053        0.009
       5     2.392        0.723        0.320        0.107
       3     2.202        1.107        0.406        0.107
       8     2.055        0.261        0.087        0.022
      11     1.825       28.196       83.844      297.278
      11     1.738       48.986       61.981      529.434
       3     1.533       11.834       19.058       23.851
       6     1.479       10.739        7.308       15.182
       3     1.431       29.902       50.618      142.246
       2     1.350        9.585       11.904       11.098
      92   overall        3.383       -0.016        7.128


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      13     7.516        1.370        2.215        2.930
      19     4.644        1.045        1.100        1.003
      27     3.492        0.847        0.638        0.439
      28     2.993        0.879        0.792        0.645
      36     2.643        1.108        1.137        1.048
      49     2.392        0.946        0.802        0.583
      38     2.202        0.768        0.544        0.331
      39     2.055        1.307        1.486        1.529
      49     1.663       31.759       79.909      577.325
      75     1.533       16.629       11.098      158.350
      74     1.479       14.356        4.968       82.103
      39     1.350        7.806       11.589       46.185
      20     1.312       38.204      253.233     2127.743
     811   overall        6.394       13.416      111.988

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.516  0.00  0.00  0.00  0.00  0.00  0.36  0.33  0.64  0.61  0.59
       3     4.644  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.492  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.993  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.643  1.34  0.96  0.80  0.70  0.64  0.59  1.12  1.06  1.01  0.98
       5     2.392  0.00  0.00  0.00  0.00  0.77  1.07  1.00  0.95  0.91  0.88
       3     2.202  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.055  2.01  1.45  1.50  1.32  1.20  1.11  1.05  1.19  1.14  1.10
      11     1.825  2.20  1.58  1.31  1.15  1.05  0.97  0.91  0.87  0.83  0.80
      11     1.738  2.20  1.58  1.31  1.15  1.05  0.97  0.91  0.87  0.83  0.80
       3     1.533  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     1.479  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       3     1.431  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     1.350  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      92   overall  0.96  0.69  0.60  0.53  0.52  0.54  0.53  0.54  0.51  0.51


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      13     7.516  1.62  2.55  1.78  2.10  1.95  1.70  1.53  1.40  1.30  1.22
      19     4.644  1.66  1.74  1.22  1.12  1.20  1.17  1.05  1.05  1.06  1.08
      27     3.492  0.78  0.61  0.57  0.67  0.66  0.90  0.88  0.94  0.94  0.88
      28     2.993  1.50  1.58  1.24  1.19  1.18  1.11  1.14  1.17  1.08  1.07
      36     2.643  0.88  1.23  1.18  1.10  0.99  0.92  0.99  1.06  1.03  1.01
      49     2.392  1.29  1.01  1.10  1.05  1.04  1.13  1.09  1.04  1.10  1.03
      38     2.202  1.66  1.31  1.02  0.80  0.87  0.93  0.84  1.05  1.02  1.04
      39     2.055  1.62  1.13  1.09  1.17  1.04  0.97  0.92  0.93  0.95  0.93
      49     1.663  5.15  2.70  1.89  1.49  1.24  1.09  1.01  0.96  0.89  0.84
      75     1.533  5.18  2.80  1.95  1.54  1.29  1.12  1.01  0.92  0.88  0.82
      74     1.479  4.26  2.31  1.62  1.27  1.10  0.96  0.91  0.86  0.84  0.79
      39     1.350  4.85  2.55  1.78  1.40  1.17  1.02  0.92  0.88  0.86  0.81
      20     1.312  5.78  3.03  2.12  1.67  1.40  1.22  1.09  1.09  1.01  0.95
     811   overall  1.97  1.22  0.92  0.79  0.71  0.66  0.62  0.61  0.60  0.57


 List of     93 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   46    6   15    1.40   67.89 -0.2371E+02  0.9981E+01 "alien"
   30    1   13    1.95   50.27 -0.2122E+02  0.1210E+02 "alien"
   49    6   13    1.39   49.63 -0.1734E+02  0.1581E+02 "alien"
   39    1   10    1.76   43.16 -0.1588E+02  0.1329E+02 "alien"
   40   33    9    1.41   39.62 -0.1384E+02  0.1093E+02 "alien"
   53    4    8    1.39   39.56 -0.1382E+02  0.1047E+02 "alien"
   36    2   13    1.74   35.61 -0.1310E+02  0.8035E+01 "alien"
   31    4   16    1.73   34.67 -0.1276E+02  0.1087E+02 "alien"
   42    8   18    1.38   34.17 -0.1193E+02  0.9226E+01 "alien"
   27   20   22    1.40   32.69 -0.1142E+02  0.8431E+01 "alien"
   53    7   13    1.30   31.51  0.5075E+02  0.1547E+02 "alien"
   39    1   11    1.73   31.37 -0.1154E+02  0.1350E+02 "alien"
   37   32   12    1.42   31.13 -0.2138E+02  0.1048E+02 "alien"
   30   26   20    1.37   29.36 -0.1025E+02  0.1501E+02 "alien"
   33   22    5    1.89   28.12 -0.1186E+02  0.1077E+02 "alien"
   21   18   15    1.89   27.59 -0.1164E+02  0.7706E+01 "alien"
   26   13   24    1.38   27.27 -0.9526E+01  0.9695E+01 "alien"
   32   23   10    1.75   25.18 -0.9265E+01  0.1577E+02 "alien"
   43    1    3    1.78   25.00 -0.9199E+01  0.9291E+01 "alien"
   30   21   10    1.85   24.08 -0.1372E+02  0.1308E+02 "alien"
   27   24   20    1.43   23.18 -0.1592E+02  0.1324E+02 "alien"
   45    1    8    1.62   22.73 -0.2005E+02  0.1024E+02 "alien"
   37    1    7    1.95   22.24 -0.9387E+01  0.9429E+01 "alien"
   25   18   12    1.98   21.77 -0.9186E+01  0.1227E+02 "alien"
   35   25    7    1.71   21.75 -0.8003E+01  0.9582E+01 "alien"
   31    6   18    1.61   21.60 -0.1905E+02  0.1152E+02 "alien"
   42    7   17    1.42   21.53 -0.1479E+02  0.9101E+01 "alien"
   25   13   24    1.39   21.37 -0.7463E+01  0.6329E+01 "alien"
   34   22    3    1.89   21.01 -0.8866E+01  0.9767E+01 "alien"
   50    5   12    1.39   20.85 -0.7283E+01  0.1090E+02 "alien"
   35   26    8    1.66   20.81  0.2189E+02  0.1727E+02 "alien"
   37    9   20    1.40   20.80 -0.7264E+01  0.8747E+01 "alien"
   39   28    2    1.60   20.37 -0.1797E+02  0.1010E+02 "alien"
   27   21   13    1.80   20.34 -0.1159E+02  0.8159E+01 "alien"
   31    5   17    1.67   19.66  0.2069E+02  0.1363E+02 "alien"
   27   18   23    1.37   19.66 -0.6866E+01  0.1051E+02 "alien"
   27   24   16    1.60   19.57 -0.1727E+02  0.1348E+02 "alien"
   49    4   11    1.43   19.29 -0.1325E+02  0.9997E+01 "alien"
   19   16   19    1.70   19.02 -0.6998E+01  0.6110E+01 "alien"
   33   32   16    1.38   18.95 -0.6618E+01  0.1675E+02 "alien"
   38    1    7    1.91   18.16 -0.7664E+01  0.9405E+01 "alien"
   38    4   15    1.59   18.16 -0.1602E+02  0.8474E+01 "alien"
   30   29   18    1.40   17.64 -0.6161E+01  0.1463E+02 "alien"
   46    1    6    1.62   17.31 -0.1527E+02  0.1356E+02 "alien"
   54    4    7    1.38   17.17 -0.5997E+01  0.1081E+02 "alien"
   24   20   15    1.78   16.94 -0.9652E+01  0.1408E+02 "alien"
   34   22    4    1.87   16.79 -0.9565E+01  0.1069E+02 "alien"
   38    1   11    1.76   16.75 -0.6165E+01  0.8661E+01 "alien"
   43    0    7    1.71   16.63 -0.1956E+02  0.1400E+02 "alien"
   33   26   11    1.63   16.52  0.1738E+02  0.1819E+02 "alien"
   22   13   24    1.42   16.49 -0.1132E+02  0.7913E+01 "alien"
   30    7   19    1.58   16.16 -0.1426E+02  0.1205E+02 "alien"
   27   20   12    1.87   16.08 -0.9163E+01  0.8335E+01 "alien"
   36    1    8    1.96   15.94 -0.6725E+01  0.9384E+01 "alien"
   39   29    5    1.56   15.46  0.2308E+02  0.1080E+02 "alien"
   53    3    3    1.45   15.42 -0.1059E+02  0.1141E+02 "alien"
   48    5   13    1.41   15.23 -0.1046E+02  0.1060E+02 "alien"
   20   19   18    1.70   14.84 -0.5461E+01  0.7235E+01 "alien"
   35    1   12    1.82   14.79 -0.8429E+01  0.8987E+01 "alien"
   21   17   20    1.60   14.51 -0.1280E+02  0.8833E+01 "alien"
   33    0   11    1.95   14.46 -0.1662E+02  0.1378E+02 "alien"
   44    1    8    1.65   14.45  0.1520E+02  0.9381E+01 "alien"
   36   24    2    1.78   14.39 -0.5296E+01  0.9637E+01 "alien"
   20   17   15    1.94   14.14 -0.5965E+01  0.7512E+01 "alien"
   36   33   14    1.37   13.48 -0.4709E+01  0.1062E+02 "alien"
   48    7   15    1.35   13.42  0.2965E+02  0.1596E+02 "alien"
   54    4    6    1.39   12.96 -0.4528E+01  0.1053E+02 "alien"
   34   26   10    1.64   12.96  0.1363E+02  0.9592E+01 "alien"
   42    6   16    1.45   12.94 -0.8888E+01  0.9285E+01 "alien"
   26   21   14    1.77   12.50 -0.4601E+01  0.8009E+01 "alien"
   43    6   16    1.43   12.07 -0.8291E+01  0.9924E+01 "alien"
   44    0    4    1.73   11.94 -0.1004E+02  0.1090E+02 "alien"
   33   29   14    1.48   11.86  0.1501E+02  0.1019E+02 "alien"
   43    4   13    1.53   11.79  0.1760E+02  0.1596E+02 "alien"
   40    2   12    1.66   11.70  0.1231E+02  0.9089E+01 "alien"
   37   25    1    1.73   11.46 -0.4219E+01  0.9631E+01 "alien"
   33   27   13    1.55   11.39  0.1700E+02  0.1462E+02 "alien"
   30   26   15    1.55   11.38  0.1700E+02  0.1307E+02 "alien"
   39   28    0    1.61   11.34 -0.2199E+02  0.1429E+02 "alien"
   42    1    2    1.83   11.24 -0.6404E+01  0.9565E+01 "alien"
   53    4    7    1.41   11.19 -0.3909E+01  0.1079E+02 "alien"
   36    9   20    1.42   10.98 -0.7539E+01  0.1170E+02 "alien"
   21   19   19    1.63   10.94  0.1151E+02  0.1072E+02 "alien"
   33   30   15    1.44   10.91 -0.7490E+01  0.9917E+01 "alien"
   38    1    6    1.94   10.84 -0.4575E+01  0.9649E+01 "alien"
   34   23    5    1.83   10.52 -0.5993E+01  0.9889E+01 "alien"
   36   31   13    1.43   10.45 -0.7178E+01  0.1045E+02 "alien"
   32   24   11    1.69   10.29  0.1083E+02  0.9866E+01 "alien"
   41    8   18    1.40   10.26 -0.3584E+01  0.1542E+02 "alien"
   33   24    9    1.73   10.23 -0.3764E+01  0.9180E+01 "alien"
   44    1    0    1.76   10.12 -0.7568E+01  0.1302E+02 "alien"
   34   31   15    1.41   10.10 -0.3528E+01  0.1061E+02 "alien"
   44    0    5    1.71   10.07 -0.7510E+01  0.9022E+01 "alien"


 List of     93 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   31.51  0.5075E+02  0.1547E+02
   48    7   15    1.35   13.42  0.2965E+02  0.1596E+02
   30   26   20    1.37   29.36 -0.1025E+02  0.1501E+02
   36   33   14    1.37   13.48 -0.4709E+01  0.1062E+02
   27   18   23    1.37   19.66 -0.6866E+01  0.1051E+02
   26   13   24    1.38   27.27 -0.9526E+01  0.9695E+01
   33   32   16    1.38   18.95 -0.6618E+01  0.1675E+02
   54    4    7    1.38   17.17 -0.5997E+01  0.1081E+02
   42    8   18    1.38   34.17 -0.1193E+02  0.9226E+01
   49    6   13    1.39   49.63 -0.1734E+02  0.1581E+02
   50    5   12    1.39   20.85 -0.7283E+01  0.1090E+02
   25   13   24    1.39   21.37 -0.7463E+01  0.6329E+01
   53    4    8    1.39   39.56 -0.1382E+02  0.1047E+02
   54    4    6    1.39   12.96 -0.4528E+01  0.1053E+02
   27   20   22    1.40   32.69 -0.1142E+02  0.8431E+01
   46    6   15    1.40   67.89 -0.2371E+02  0.9981E+01
   37    9   20    1.40   20.80 -0.7264E+01  0.8747E+01
   41    8   18    1.40   10.26 -0.3584E+01  0.1542E+02
   30   29   18    1.40   17.64 -0.6161E+01  0.1463E+02
   53    4    7    1.41   11.19 -0.3909E+01  0.1079E+02
   34   31   15    1.41   10.10 -0.3528E+01  0.1061E+02
   40   33    9    1.41   39.62 -0.1384E+02  0.1093E+02
   48    5   13    1.41   15.23 -0.1046E+02  0.1060E+02
   36    9   20    1.42   10.98 -0.7539E+01  0.1170E+02
   37   32   12    1.42   31.13 -0.2138E+02  0.1048E+02
   42    7   17    1.42   21.53 -0.1479E+02  0.9101E+01
   22   13   24    1.42   16.49 -0.1132E+02  0.7913E+01
   36   31   13    1.43   10.45 -0.7178E+01  0.1045E+02
   43    6   16    1.43   12.07 -0.8291E+01  0.9924E+01
   27   24   20    1.43   23.18 -0.1592E+02  0.1324E+02
   49    4   11    1.43   19.29 -0.1325E+02  0.9997E+01
   33   30   15    1.44   10.91 -0.7490E+01  0.9917E+01
   53    3    3    1.45   15.42 -0.1059E+02  0.1141E+02
   42    6   16    1.45   12.94 -0.8888E+01  0.9285E+01
   33   29   14    1.48   11.86  0.1501E+02  0.1019E+02
   43    4   13    1.53   11.79  0.1760E+02  0.1596E+02
   33   27   13    1.55   11.39  0.1700E+02  0.1462E+02
   30   26   15    1.55   11.38  0.1700E+02  0.1307E+02
   39   29    5    1.56   15.46  0.2308E+02  0.1080E+02
   30    7   19    1.58   16.16 -0.1426E+02  0.1205E+02
   38    4   15    1.59   18.16 -0.1602E+02  0.8474E+01
   27   24   16    1.60   19.57 -0.1727E+02  0.1348E+02
   21   17   20    1.60   14.51 -0.1280E+02  0.8833E+01
   39   28    2    1.60   20.37 -0.1797E+02  0.1010E+02
   39   28    0    1.61   11.34 -0.2199E+02  0.1429E+02
   31    6   18    1.61   21.60 -0.1905E+02  0.1152E+02
   45    1    8    1.62   22.73 -0.2005E+02  0.1024E+02
   46    1    6    1.62   17.31 -0.1527E+02  0.1356E+02
   33   26   11    1.63   16.52  0.1738E+02  0.1819E+02
   21   19   19    1.63   10.94  0.1151E+02  0.1072E+02
   34   26   10    1.64   12.96  0.1363E+02  0.9592E+01
   44    1    8    1.65   14.45  0.1520E+02  0.9381E+01
   40    2   12    1.66   11.70  0.1231E+02  0.9089E+01
   35   26    8    1.66   20.81  0.2189E+02  0.1727E+02
   31    5   17    1.67   19.66  0.2069E+02  0.1363E+02
   32   24   11    1.69   10.29  0.1083E+02  0.9866E+01
   19   16   19    1.70   19.02 -0.6998E+01  0.6110E+01
   20   19   18    1.70   14.84 -0.5461E+01  0.7235E+01
   35   25    7    1.71   21.75 -0.8003E+01  0.9582E+01
   43    0    7    1.71   16.63 -0.1956E+02  0.1400E+02
   44    0    5    1.71   10.07 -0.7510E+01  0.9022E+01
   39    1   11    1.73   31.37 -0.1154E+02  0.1350E+02
   44    0    4    1.73   11.94 -0.1004E+02  0.1090E+02
   31    4   16    1.73   34.67 -0.1276E+02  0.1087E+02
   37   25    1    1.73   11.46 -0.4219E+01  0.9631E+01
   33   24    9    1.73   10.23 -0.3764E+01  0.9180E+01
   36    2   13    1.74   35.61 -0.1310E+02  0.8035E+01
   32   23   10    1.75   25.18 -0.9265E+01  0.1577E+02
   44    1    0    1.76   10.12 -0.7568E+01  0.1302E+02
   38    1   11    1.76   16.75 -0.6165E+01  0.8661E+01
   39    1   10    1.76   43.16 -0.1588E+02  0.1329E+02
   26   21   14    1.77   12.50 -0.4601E+01  0.8009E+01
   36   24    2    1.78   14.39 -0.5296E+01  0.9637E+01
   43    1    3    1.78   25.00 -0.9199E+01  0.9291E+01
   24   20   15    1.78   16.94 -0.9652E+01  0.1408E+02
   27   21   13    1.80   20.34 -0.1159E+02  0.8159E+01
   35    1   12    1.82   14.79 -0.8429E+01  0.8987E+01
   34   23    5    1.83   10.52 -0.5993E+01  0.9889E+01
   42    1    2    1.83   11.24 -0.6404E+01  0.9565E+01
   30   21   10    1.85   24.08 -0.1372E+02  0.1308E+02
   27   20   12    1.87   16.08 -0.9163E+01  0.8335E+01
   34   22    4    1.87   16.79 -0.9565E+01  0.1069E+02
   34   22    3    1.89   21.01 -0.8866E+01  0.9767E+01
   33   22    5    1.89   28.12 -0.1186E+02  0.1077E+02
   21   18   15    1.89   27.59 -0.1164E+02  0.7706E+01
   38    1    7    1.91   18.16 -0.7664E+01  0.9405E+01
   38    1    6    1.94   10.84 -0.4575E+01  0.9649E+01
   20   17   15    1.94   14.14 -0.5965E+01  0.7512E+01
   30    1   13    1.95   50.27 -0.2122E+02  0.1210E+02
   33    0   11    1.95   14.46 -0.1662E+02  0.1378E+02
   37    1    7    1.95   22.24 -0.9387E+01  0.9429E+01
   36    1    8    1.96   15.94 -0.6725E+01  0.9384E+01
   25   18   12    1.98   21.77 -0.9186E+01  0.1227E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        34
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1508
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1412
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       96


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
