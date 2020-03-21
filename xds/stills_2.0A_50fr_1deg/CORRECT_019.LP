
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.435    77.435    38.755  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      19      19
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      19
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000845  0.000685
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
 REFLECTING_RANGE_E.S.D.=   0.112
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

    1405 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      19 ...      19
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1405 REFLECTIONS ACCEPTED



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
 *  31        aP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0  1  0  0  1  0  0  0
 *  34        mP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  35        mP          0.0      77.4   38.8   77.4  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  32        oP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.5  109.5   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.5  109.5   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.5  109.5   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.4   77.4   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    28        mC        250.0      38.8  159.6   77.4  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.8  159.6   77.4  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    37        mC        250.0     159.6   38.8   77.4  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.6   38.8   77.4  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.8  159.6   77.4  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.8  159.6   77.4  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    27        mC        500.0     159.6   38.8  109.5  90.0 133.3  76.0    0  2  1  0  0  0  1  0  1 -1  0  0
    19        oI        500.0      38.8  109.5  116.2  90.0  70.5  90.0    0  0  1  0  1 -1  0  0  1  1  1  0
    26        oF        625.0      38.8  159.6  159.6  86.6 104.0 104.0    0  0 -1  0  0  2  1  0  2  0  1  0
    18        tI        625.0     109.5  116.2   38.8  70.5  90.0  90.0    1 -1  0  0 -1 -1 -1  0  0  0 -1  0
    30        mC        998.1      77.4  173.2   38.8  90.0  90.0  63.4    0  1  0  0 -2  1  0  0  0  0  1  0
    41        mC        998.1     173.2   77.4   38.8  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    40        oC        998.1      77.4  173.2   38.8  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.1      77.4   77.4   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.2  116.2  116.2  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.6   86.6  116.2 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.6   86.6  109.5  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.6   86.6  116.2 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.5   86.6   86.6  78.4  50.8  50.8    1  1  0  0  1  0 -1  0  0  1 -1  0
     7        tI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    1  0 -1  0  0  1 -1  0  1  1  0  0
     8        oI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    0 -1  1  0 -1  0  1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.6  247.9 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.6   86.6   77.4  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.6   86.6   77.4  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.8   77.4  177.4  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.6   86.6  177.4 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.6   86.6   77.4  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.6   86.6   86.6  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.4  109.5   38.8  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.8   77.4  177.4 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.6  177.4   77.4 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    5.4      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.801   -23.444     0.035
 UNIT_CELL_B-AXIS=    23.444    73.801     0.016
 UNIT_CELL_C-AXIS=    -0.019    -0.002    38.755
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       258 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.27
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.354    77.354    38.715  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.2E-01 1.2E-01 5.0E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012928  0.012928  0.025830  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.723   -23.421     0.033
 COORDINATES OF UNIT CELL B-AXIS    23.421    73.723     0.011
 COORDINATES OF UNIT CELL C-AXIS    -0.018    -0.000    38.715
 CRYSTAL MOSAICITY (DEGREES)     0.112
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006611  0.001813
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.91   1232.91
 DETECTOR ORIGIN (PIXELS) AT                     1261.28   1242.53
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.08
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.35    77.35    38.71  90.000  90.000  90.000


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
  0.05  0.15        86         99.460       1.755
  0.15  0.25       167        208.692       0.836
  0.25  0.35       196        132.794       1.314
  0.35  0.45       149        103.970       1.678
  0.45  0.55       169        246.601       0.708
  0.55  0.65       195        158.091       1.104
  0.65  0.75       187        182.094       0.958
  0.75  0.85       172        150.435       1.160
  0.85  0.95        84        331.397       0.527
 
 -0.05  0.95      1405        174.514       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               240.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.063
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    18.1 XMAX=    18.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59599 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        508


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               240.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.061
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     8.7 XMAX=  2456.0 NXBIN=    3
 YMIN=     6.9 YMAX=  2516.6 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        508


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               240.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.056
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    18.1 XMAX=    18.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1262
 DETECTOR_SURFACE_POSITION=    1648    1688
 DETECTOR_SURFACE_POSITION=     817    1688
 DETECTOR_SURFACE_POSITION=     817     835
 DETECTOR_SURFACE_POSITION=    1648     835
 DETECTOR_SURFACE_POSITION=    2173    1661
 DETECTOR_SURFACE_POSITION=    1622    2226
 DETECTOR_SURFACE_POSITION=     843    2226
 DETECTOR_SURFACE_POSITION=     292    1661
 DETECTOR_SURFACE_POSITION=     292     862
 DETECTOR_SURFACE_POSITION=     843     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2173     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        508


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
 9.913E-01  1.377E-03   27.07


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

   24.461   3.151    21.96   1.57      2.31      2.23      66      94      30
    3.151   2.237    11.00   1.43      5.14      4.06     186     221      10
    2.237   1.829     2.14   0.97     20.36     20.53     208     244       0
    1.829   1.585    -0.04   0.68    -99.90    -99.90     242     288       0
    1.585   1.419     0.06   0.62    486.25    625.69     282     341       0
    1.419   1.295     0.03   0.85    719.31    774.23     126     176       0
   --------------------------------------------------------------------------
   24.461   1.295     3.69   0.92      6.42      6.58    1110    1364      40
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      19      19)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.682  -0.1342E+02  0.1205E+02        -1.11           1 
   47    0    0     1.646   0.6254E+01  0.1180E+02         0.53           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     3.4% OF MEAN INTENSITY


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

     3.84          61      40      1285        3.1%       2.3%      2.2%       42   27.64      3.2%    99.9*     0    0.000       0
     2.73          95      57      2098        2.7%       2.7%      2.9%       76   21.25      3.8%    99.7*     0    0.000       0
     2.24         159      92      2665        3.5%       6.8%      4.7%      134   12.89      9.7%    98.8*     0    0.000       0
     1.94         158      90      3101        2.9%      16.0%     14.0%      136    4.27     22.6%    96.0*     0    0.000       0
     1.74         170      99      3481        2.8%     511.2%    748.7%      142    0.06    723.0%    14.4      0    0.000       0
     1.59         204     118      3847        3.1%     -99.9%    -99.9%      172    0.00    -99.9%    -5.8      0    0.000       0
     1.47         222     132      4142        3.2%     768.4%    903.8%      180    0.03   1086.7%    12.4      0    0.000       0
     1.37         210     120      4431        2.7%     708.3%    937.7%      180    0.08   1001.6%    20.6      0    0.000       0
     1.30          85      61      4743        1.3%     221.0%    222.4%       48    0.03    312.6%    63.8*     0    0.000       0
    total        1364     809     29793        2.7%       6.4%      6.6%     1110    4.81      9.1%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1405
 NUMBER OF REJECTED MISFITS                                40
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1364
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    809



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

   WILSON LINE (using all data) : A=   7.834 B=  33.195 CORRELATION=  0.93
      #      RES      SS        <I>       log(<I>)       BO
      17     7.254   0.005  2.0311E+03       7.616      22.9
      19     4.559   0.012  2.5541E+03       7.845      -0.5
      22     3.541   0.020  1.1096E+03       7.012      20.6
      30     2.991   0.028  3.6809E+02       5.908      34.5
      37     2.648   0.036  3.5780E+02       5.880      27.4
      51     2.399   0.043  1.8621E+02       5.227      30.0
      37     2.206   0.051  1.0838E+02       4.686      30.7
      45     2.057   0.059  4.3846E+01       3.781      34.3
      44     1.931   0.067 -8.6711E-01    -999.999   -1000.0
      53     1.829   0.075  1.2148E+00       0.195      51.1
      49     1.739   0.083 -1.4112E-01    -999.999   -1000.0
      52     1.665   0.090 -3.4901E-01    -999.999   -1000.0
      52     1.597   0.098 -3.7024E-01    -999.999   -1000.0
      68     1.536   0.106  1.9304E-01      -1.645      44.7
      72     1.480   0.114 -2.2999E-01    -999.999   -1000.0
      55     1.433   0.122  1.7109E+00       0.537      30.0
      46     1.391   0.129  5.4739E-01      -0.603      32.6
      42     1.351   0.137 -6.2619E-01    -999.999   -1000.0
      18     1.315   0.145  1.3038E+00       0.265      26.2


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.254        1.834        1.234        0.662
       3     4.559        1.195        0.711        0.314
       2     3.541        0.370        0.078        0.012
       4     2.991        1.262        0.595        0.213
       3     2.648        0.304        0.072        0.012
       5     2.399        0.987        0.554        0.239
       2     2.206        1.361        0.557        0.164
       8     2.057        0.419        0.157        0.050
      10     1.829        9.888       -0.358       16.224
      11     1.665       35.054      153.815      653.882
       3     1.433        6.513        6.988        5.752
      83   overall       26.132     -250.159     1953.597


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.254        0.456        0.319        0.180
      16     4.559        1.061        1.145        1.032
      20     3.541        0.821        0.636        0.457
      26     2.991        0.840        0.695        0.486
      34     2.648        1.159        1.320        1.404
      46     2.399        0.934        0.814        0.604
      35     2.206        0.779        0.580        0.377
      37     2.057        1.431        1.975        2.531
      43     1.829       22.210       12.559      164.454
      42     1.739       36.579      -15.531      765.436
      49     1.597       36.229        2.389      689.319
      65     1.536     1147.622     1892.329   616044.250
      52     1.433       15.122        5.285      128.577
      46     1.391      138.372      547.466    13550.943
      41     1.351      131.244     2514.885    56557.371
      18     1.315       53.557       85.219     1784.334
     726   overall      115.550      207.662    56077.573

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.254  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.32  0.30  0.29
       3     4.559  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.541  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       4     2.991  0.00  0.00  0.00  0.00  0.00  0.00  0.42  0.40  0.38  0.37
       3     2.648  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       5     2.399  0.00  0.00  0.00  0.42  0.77  1.07  1.00  0.95  0.91  0.88
       2     2.206  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.057  1.01  0.72  0.90  0.79  0.96  0.89  1.05  0.99  0.95  0.92
      10     1.829  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.73
      11     1.665  2.20  1.58  1.31  1.15  1.05  0.97  0.91  0.87  0.83  0.80
       3     1.433  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      83   overall  0.68  0.49  0.43  0.41  0.42  0.41  0.42  0.42  0.40  0.42


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.254  1.75  2.76  1.93  1.52  1.91  1.66  1.66  1.51  1.40  1.32
      16     4.559  1.31  1.38  1.21  0.95  1.43  1.25  1.12  1.36  1.26  1.19
      20     3.541  0.53  0.83  0.58  0.45  0.64  0.78  0.99  0.91  0.93  1.11
      26     2.991  2.02  1.91  1.34  1.17  1.27  1.28  1.22  1.19  1.17  1.10
      34     2.648  0.62  1.14  1.25  1.16  0.97  0.98  1.11  1.12  1.14  1.07
      46     2.399  0.91  0.96  0.92  1.12  1.16  1.06  1.08  1.11  1.14  1.07
      35     2.206  0.90  0.95  0.77  0.78  0.80  1.01  0.91  0.88  1.06  1.13
      37     2.057  1.42  1.49  1.04  1.07  1.10  1.02  1.07  1.03  1.09  1.07
      43     1.829  4.89  2.57  1.79  1.41  1.18  1.03  0.97  0.89  0.82  0.77
      42     1.739  6.51  3.42  2.39  1.88  1.57  1.37  1.23  1.12  1.04  0.98
      49     1.597  5.36  2.81  1.97  1.55  1.30  1.13  1.01  0.93  0.86  0.81
      65     1.536  5.01  2.63  1.84  1.45  1.21  1.06  0.98  0.89  0.83  0.78
      52     1.433  4.45  2.33  1.63  1.34  1.12  1.02  0.92  0.94  0.91  0.85
      46     1.391  5.71  3.00  2.10  1.65  1.44  1.25  1.12  1.07  0.99  0.93
      41     1.351  4.36  2.29  1.60  1.26  1.05  0.92  0.82  0.75  0.70  0.66
      18     1.315  4.67  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
     726   overall  2.87  1.73  1.26  1.05  0.96  0.88  0.83  0.80  0.78  0.75


 List of    146 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   36   25    5    1.72  142.70 -0.2014E+02  0.8636E+01 "alien"
   41    4   14    1.55  117.98  0.2277E+02  0.8905E+01 "alien"
   36    2   13    1.74  109.66 -0.1547E+02  0.7569E+01 "alien"
   31    4   16    1.73  108.70 -0.1534E+02  0.1243E+02 "alien"
   35   25    7    1.71   98.10 -0.1384E+02  0.8687E+01 "alien"
   32   23   10    1.75   96.39 -0.1360E+02  0.1560E+02 "alien"
   40    1   10    1.73   88.17 -0.1244E+02  0.8118E+01 "alien"
   30   26   15    1.56   86.45  0.1669E+02  0.1241E+02 "alien"
   48    2    8    1.53   85.29  0.1646E+02  0.9435E+01 "alien"
   50    1    2    1.54   79.94  0.1543E+02  0.9589E+01 "alien"
   15    8   20    1.78   78.62 -0.1110E+02  0.5813E+01 "alien"
   29    8   20    1.55   75.01  0.1448E+02  0.1152E+02 "alien"
   37   29   10    1.51   74.87  0.1445E+02  0.1652E+02 "alien"
   23   13   23    1.46   72.36 -0.1664E+02  0.8823E+01 "alien"
   51    2    5    1.49   69.08 -0.1589E+02  0.9213E+01 "alien"
   47    3   10    1.51   68.16  0.1316E+02  0.9240E+01 "alien"
   25   20   21    1.47   66.92 -0.1539E+02  0.1085E+02 "alien"
   44    3   12    1.54   65.81  0.1270E+02  0.9176E+01 "alien"
   41   31    2    1.50   64.83 -0.1491E+02  0.9592E+01 "alien"
   50    2    5    1.52   64.17  0.1239E+02  0.9237E+01 "alien"
   39   31    8    1.48   56.44 -0.1298E+02  0.9260E+01 "alien"
   27   24   16    1.60   56.27 -0.2083E+02  0.1431E+02 "alien"
   49    3    9    1.48   55.36 -0.1273E+02  0.9096E+01 "alien"
   33   25   11    1.65   55.07 -0.1922E+02  0.1375E+02 "alien"
   52    2    2    1.48   50.92 -0.1171E+02  0.9415E+01 "alien"
   33   25   10    1.68   50.90 -0.1776E+02  0.1309E+02 "alien"
   38    5   16    1.55   50.33  0.9715E+01  0.8082E+01 "alien"
   21   10   23    1.50   50.21 -0.1155E+02  0.8121E+01 "alien"
   39    5   16    1.53   50.20  0.9690E+01  0.8317E+01 "alien"
   36   26    6    1.68   49.41 -0.1724E+02  0.1489E+02 "alien"
   39    8   19    1.41   49.08  0.2687E+02  0.1341E+02 "alien"
   32   27   14    1.54   48.96  0.9451E+01  0.1374E+02 "alien"
   37    2   13    1.71   44.84 -0.6328E+01  0.8527E+01 "alien"
   52    5   10    1.38   44.77  0.2451E+02  0.1015E+02 "alien"
   44    4   13    1.51   44.39  0.8570E+01  0.8867E+01 "alien"
   42   32    1    1.46   44.25 -0.1018E+02  0.9950E+01 "alien"
   31   25   13    1.63   43.81 -0.1622E+02  0.8169E+01 "alien"
   32   24   11    1.69   43.57 -0.1521E+02  0.8713E+01 "alien"
   47    1    5    1.61   42.88 -0.1588E+02  0.8982E+01 "alien"
   39    2   12    1.69   42.44 -0.1481E+02  0.7959E+01 "alien"
   42    4   14    1.53   41.70  0.8050E+01  0.8686E+01 "alien"
   49    2    7    1.52   40.52  0.7822E+01  0.9855E+01 "alien"
   16    8   20    1.77   38.67 -0.5457E+01  0.5586E+01 "alien"
   48    2    7    1.55   37.23  0.7188E+01  0.8968E+01 "alien"
   22   10   22    1.54   36.47  0.7040E+01  0.7916E+01 "alien"
   37   26    5    1.67   34.21 -0.1194E+02  0.8645E+01 "alien"
   35    6   18    1.53   33.83  0.6531E+01  0.1332E+02 "alien"
   39   31    9    1.46   33.74 -0.7759E+01  0.1044E+02 "alien"
   25   24   17    1.59   33.23 -0.1230E+02  0.8612E+01 "alien"
   44    5   14    1.48   32.79 -0.7541E+01  0.1481E+02 "alien"
   20   13   21    1.60   32.36 -0.1198E+02  0.8249E+01 "alien"
   52    2    1    1.49   31.90 -0.7336E+01  0.9260E+01 "alien"
   34    6   18    1.55   31.76  0.6132E+01  0.7664E+01 "alien"
   27   22   14    1.73   31.72 -0.4477E+01  0.7353E+01 "alien"
   33   12   23    1.34   31.37 -0.1964E+02  0.1385E+02 "alien"
   30   23   12    1.73   29.25 -0.4128E+01  0.8298E+01 "alien"
   46    1    7    1.61   27.72 -0.1026E+02  0.8616E+01 "alien"
   50    3    8    1.47   27.68 -0.6366E+01  0.9334E+01 "alien"
   48    3    9    1.51   27.19 -0.6253E+01  0.1043E+02 "alien"
   53    7   13    1.30   26.66  0.3475E+02  0.1450E+02 "alien"
   49    6   14    1.36   26.18 -0.1640E+02  0.1633E+02 "alien"
   35   33   15    1.36   26.01 -0.1629E+02  0.1370E+02 "alien"
   51    3    7    1.46   25.88 -0.5952E+01  0.9198E+01 "alien"
   38   34   13    1.35   25.87 -0.1620E+02  0.1000E+02 "alien"
   36   32   14    1.39   25.66  0.1405E+02  0.1813E+02 "alien"
   36   27    9    1.60   25.61 -0.9482E+01  0.9559E+01 "alien"
   39   30    6    1.53   25.57  0.4936E+01  0.1145E+02 "alien"
   32    8   20    1.49   24.22 -0.5571E+01  0.7561E+01 "alien"
   48    7   15    1.36   24.20 -0.1515E+02  0.1292E+02 "alien"
   40   31    5    1.50   23.97 -0.5513E+01  0.1674E+02 "alien"
   47    2    9    1.54   23.90  0.4614E+01  0.1587E+02 "alien"
   38   29    8    1.53   23.89  0.4612E+01  0.9125E+01 "alien"
   39   33   10    1.41   23.67  0.1296E+02  0.1183E+02 "alien"
   38    9   20    1.38   23.04  0.1261E+02  0.9461E+01 "alien"
   42    5   15    1.49   23.03 -0.5297E+01  0.8680E+01 "alien"
   53    2    0    1.46   22.64 -0.2706E+02  0.1505E+02 "alien"
   31   30   18    1.38   22.51  0.1232E+02  0.9928E+01 "alien"
   21   10   22    1.56   22.35  0.4315E+01  0.7944E+01 "alien"
   33   22    5    1.89   22.14 -0.1919E+02  0.1511E+02 "alien"
   44    1    0    1.76   22.00 -0.1528E+02  0.1283E+02 "alien"
   38   31   10    1.46   21.50 -0.4944E+01  0.9348E+01 "alien"
   19   18   18    1.74   21.26 -0.3000E+01  0.6282E+01 "alien"
   35   25    8    1.69   21.19 -0.7397E+01  0.1408E+02 "alien"
   33   29   14    1.49   21.04 -0.4839E+01  0.1018E+02 "alien"
   34   27   11    1.59   20.26 -0.7501E+01  0.9459E+01 "alien"
   44    0    5    1.71   19.78 -0.1142E+02  0.8413E+01 "alien"
   21   17   14    1.99   19.68 -0.1706E+02  0.8117E+01 "alien"
   37   30   11    1.47   19.65 -0.4518E+01  0.1096E+02 "alien"
   25   10   22    1.50   19.63 -0.4515E+01  0.8605E+01 "alien"
   51    2    2    1.51   19.52  0.3768E+01  0.1085E+02 "alien"
   38   29    7    1.55   19.48  0.3761E+01  0.9319E+01 "alien"
   38   28    5    1.60   19.02 -0.7042E+01  0.9141E+01 "alien"
   45    8   17    1.36   18.94 -0.1186E+02  0.1065E+02 "alien"
   29   23   13    1.71   18.66 -0.2633E+01  0.7774E+01 "alien"
   17   11   20    1.73   18.62 -0.2628E+01  0.6458E+01 "alien"
   28   21   22    1.38   18.18  0.9951E+01  0.8070E+01 "alien"
   39   29    4    1.57   17.85 -0.6608E+01  0.9554E+01 "alien"
   41    9   19    1.37   17.65 -0.1105E+02  0.1369E+02 "alien"
   38    6   17    1.51   17.56 -0.4039E+01  0.8365E+01 "alien"
   50    7   14    1.34   17.50 -0.1096E+02  0.1021E+02 "alien"
   41   32    4    1.47   17.19 -0.3953E+01  0.1068E+02 "alien"
   24   10   22    1.51   17.08  0.3296E+01  0.8184E+01 "alien"
   33    5   17    1.62   16.82 -0.6228E+01  0.7041E+01 "alien"
   43    3   12    1.57   16.69 -0.6181E+01  0.8503E+01 "alien"
   27   23   15    1.67   16.20 -0.5654E+01  0.7679E+01 "alien"
   19    9   21    1.65   15.95 -0.5568E+01  0.6454E+01 "alien"
   30   20    8    1.96   15.68 -0.1359E+02  0.8636E+01 "alien"
   32   21    6    1.93   15.58 -0.1351E+02  0.8607E+01 "alien"
   22   16   21    1.55   15.50  0.2992E+01  0.6219E+01 "alien"
   33   33   17    1.34   15.10 -0.2699E+02  0.1370E+02 "alien"
   27   20   22    1.40   14.67  0.8028E+01  0.7602E+01 "alien"
   52    5   11    1.36   14.64 -0.9166E+01  0.1647E+02 "alien"
   26   10   22    1.49   14.57 -0.3351E+01  0.9295E+01 "alien"
   37    3   14    1.66   14.26 -0.4977E+01  0.7608E+01 "alien"
   28   17   24    1.33   14.24 -0.8919E+01  0.1258E+02 "alien"
   43    1    9    1.66   14.06 -0.4908E+01  0.8575E+01 "alien"
   30    1   13    1.95   13.85 -0.1201E+02  0.1200E+02 "alien"
   36    1    8    1.96   13.52 -0.1172E+02  0.8527E+01 "alien"
   40   29    0    1.57   13.45 -0.1262E+02  0.1288E+02 "alien"
   31   29   19    1.36   13.28 -0.8315E+01  0.9077E+01 "alien"
   31   24   12    1.68   13.21 -0.4612E+01  0.7851E+01 "alien"
   42    0    7    1.75   13.08 -0.4562E+01  0.8480E+01 "alien"
   33    8   20    1.48   13.07 -0.3007E+01  0.1172E+02 "alien"
   34   29   13    1.50   13.00 -0.2991E+01  0.9257E+01 "alien"
   52    2    0    1.49   12.84 -0.7178E+01  0.1309E+02 "alien"
   43    2   11    1.60   12.83 -0.4750E+01  0.9677E+01 "alien"
   31   22   10    1.80   12.81  0.1556E+02  0.1430E+02 "alien"
   41   31    1    1.50   12.81 -0.2945E+01  0.9249E+01 "alien"
   31   26   21    1.33   12.73  0.1660E+02  0.1426E+02 "alien"
   29   26   20    1.39   12.71  0.6955E+01  0.1568E+02 "alien"
   28    2   14    1.95   12.51 -0.1085E+02  0.1154E+02 "alien"
   39   32    9    1.44   12.41  0.2123E+02  0.1070E+02 "alien"
   34   24    8    1.74   12.26 -0.1730E+01  0.8417E+01 "alien"
   37   26    3    1.70   11.90 -0.4154E+01  0.9178E+01 "alien"
   26   20   13    1.85   11.81  0.1435E+02  0.7720E+01 "alien"
   23   23   18    1.60   11.64 -0.9665E+01  0.1040E+02 "alien"
   23   17   21    1.52   11.55  0.2229E+01  0.1054E+02 "alien"
   49    7   15    1.34   11.46 -0.7174E+01  0.1501E+02 "alien"
   48    1    4    1.59   11.33 -0.4194E+01  0.8823E+01 "alien"
   37    8   19    1.44   11.02  0.1885E+02  0.1245E+02 "alien"
   50    2    6    1.50   10.88 -0.2502E+01  0.9070E+01 "alien"
   42   32    3    1.46   10.86  0.1857E+02  0.1247E+02 "alien"
   22   11   23    1.48   10.60 -0.2438E+01  0.9089E+01 "alien"
   53    3    6    1.42   10.58  0.1810E+02  0.1153E+02 "alien"
   45    0    5    1.68   10.33 -0.7420E+01  0.8669E+01 "alien"
   26   13   24    1.38   10.24  0.5608E+01  0.9471E+01 "alien"


 List of    146 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   26.66  0.3475E+02  0.1450E+02
   31   26   21    1.33   12.73  0.1660E+02  0.1426E+02
   28   17   24    1.33   14.24 -0.8919E+01  0.1258E+02
   49    7   15    1.34   11.46 -0.7174E+01  0.1501E+02
   33   12   23    1.34   31.37 -0.1964E+02  0.1385E+02
   33   33   17    1.34   15.10 -0.2699E+02  0.1370E+02
   50    7   14    1.34   17.50 -0.1096E+02  0.1021E+02
   38   34   13    1.35   25.87 -0.1620E+02  0.1000E+02
   48    7   15    1.36   24.20 -0.1515E+02  0.1292E+02
   31   29   19    1.36   13.28 -0.8315E+01  0.9077E+01
   45    8   17    1.36   18.94 -0.1186E+02  0.1065E+02
   49    6   14    1.36   26.18 -0.1640E+02  0.1633E+02
   35   33   15    1.36   26.01 -0.1629E+02  0.1370E+02
   52    5   11    1.36   14.64 -0.9166E+01  0.1647E+02
   41    9   19    1.37   17.65 -0.1105E+02  0.1369E+02
   28   21   22    1.38   18.18  0.9951E+01  0.8070E+01
   31   30   18    1.38   22.51  0.1232E+02  0.9928E+01
   26   13   24    1.38   10.24  0.5608E+01  0.9471E+01
   52    5   10    1.38   44.77  0.2451E+02  0.1015E+02
   38    9   20    1.38   23.04  0.1261E+02  0.9461E+01
   29   26   20    1.39   12.71  0.6955E+01  0.1568E+02
   36   32   14    1.39   25.66  0.1405E+02  0.1813E+02
   27   20   22    1.40   14.67  0.8028E+01  0.7602E+01
   39    8   19    1.41   49.08  0.2687E+02  0.1341E+02
   39   33   10    1.41   23.67  0.1296E+02  0.1183E+02
   53    3    6    1.42   10.58  0.1810E+02  0.1153E+02
   37    8   19    1.44   11.02  0.1885E+02  0.1245E+02
   39   32    9    1.44   12.41  0.2123E+02  0.1070E+02
   42   32    3    1.46   10.86  0.1857E+02  0.1247E+02
   53    2    0    1.46   22.64 -0.2706E+02  0.1505E+02
   23   13   23    1.46   72.36 -0.1664E+02  0.8823E+01
   51    3    7    1.46   25.88 -0.5952E+01  0.9198E+01
   39   31    9    1.46   33.74 -0.7759E+01  0.1044E+02
   38   31   10    1.46   21.50 -0.4944E+01  0.9348E+01
   42   32    1    1.46   44.25 -0.1018E+02  0.9950E+01
   25   20   21    1.47   66.92 -0.1539E+02  0.1085E+02
   41   32    4    1.47   17.19 -0.3953E+01  0.1068E+02
   50    3    8    1.47   27.68 -0.6366E+01  0.9334E+01
   37   30   11    1.47   19.65 -0.4518E+01  0.1096E+02
   33    8   20    1.48   13.07 -0.3007E+01  0.1172E+02
   44    5   14    1.48   32.79 -0.7541E+01  0.1481E+02
   39   31    8    1.48   56.44 -0.1298E+02  0.9260E+01
   49    3    9    1.48   55.36 -0.1273E+02  0.9096E+01
   52    2    2    1.48   50.92 -0.1171E+02  0.9415E+01
   22   11   23    1.48   10.60 -0.2438E+01  0.9089E+01
   33   29   14    1.49   21.04 -0.4839E+01  0.1018E+02
   52    2    1    1.49   31.90 -0.7336E+01  0.9260E+01
   26   10   22    1.49   14.57 -0.3351E+01  0.9295E+01
   52    2    0    1.49   12.84 -0.7178E+01  0.1309E+02
   51    2    5    1.49   69.08 -0.1589E+02  0.9213E+01
   42    5   15    1.49   23.03 -0.5297E+01  0.8680E+01
   32    8   20    1.49   24.22 -0.5571E+01  0.7561E+01
   34   29   13    1.50   13.00 -0.2991E+01  0.9257E+01
   40   31    5    1.50   23.97 -0.5513E+01  0.1674E+02
   41   31    2    1.50   64.83 -0.1491E+02  0.9592E+01
   25   10   22    1.50   19.63 -0.4515E+01  0.8605E+01
   21   10   23    1.50   50.21 -0.1155E+02  0.8121E+01
   50    2    6    1.50   10.88 -0.2502E+01  0.9070E+01
   41   31    1    1.50   12.81 -0.2945E+01  0.9249E+01
   48    3    9    1.51   27.19 -0.6253E+01  0.1043E+02
   38    6   17    1.51   17.56 -0.4039E+01  0.8365E+01
   44    4   13    1.51   44.39  0.8570E+01  0.8867E+01
   51    2    2    1.51   19.52  0.3768E+01  0.1085E+02
   47    3   10    1.51   68.16  0.1316E+02  0.9240E+01
   37   29   10    1.51   74.87  0.1445E+02  0.1652E+02
   24   10   22    1.51   17.08  0.3296E+01  0.8184E+01
   50    2    5    1.52   64.17  0.1239E+02  0.9237E+01
   49    2    7    1.52   40.52  0.7822E+01  0.9855E+01
   23   17   21    1.52   11.55  0.2229E+01  0.1054E+02
   39    5   16    1.53   50.20  0.9690E+01  0.8317E+01
   39   30    6    1.53   25.57  0.4936E+01  0.1145E+02
   48    2    8    1.53   85.29  0.1646E+02  0.9435E+01
   42    4   14    1.53   41.70  0.8050E+01  0.8686E+01
   35    6   18    1.53   33.83  0.6531E+01  0.1332E+02
   38   29    8    1.53   23.89  0.4612E+01  0.9125E+01
   47    2    9    1.54   23.90  0.4614E+01  0.1587E+02
   32   27   14    1.54   48.96  0.9451E+01  0.1374E+02
   44    3   12    1.54   65.81  0.1270E+02  0.9176E+01
   50    1    2    1.54   79.94  0.1543E+02  0.9589E+01
   22   10   22    1.54   36.47  0.7040E+01  0.7916E+01
   48    2    7    1.55   37.23  0.7188E+01  0.8968E+01
   29    8   20    1.55   75.01  0.1448E+02  0.1152E+02
   22   16   21    1.55   15.50  0.2992E+01  0.6219E+01
   38    5   16    1.55   50.33  0.9715E+01  0.8082E+01
   34    6   18    1.55   31.76  0.6132E+01  0.7664E+01
   38   29    7    1.55   19.48  0.3761E+01  0.9319E+01
   41    4   14    1.55  117.98  0.2277E+02  0.8905E+01
   30   26   15    1.56   86.45  0.1669E+02  0.1241E+02
   21   10   22    1.56   22.35  0.4315E+01  0.7944E+01
   40   29    0    1.57   13.45 -0.1262E+02  0.1288E+02
   43    3   12    1.57   16.69 -0.6181E+01  0.8503E+01
   39   29    4    1.57   17.85 -0.6608E+01  0.9554E+01
   48    1    4    1.59   11.33 -0.4194E+01  0.8823E+01
   34   27   11    1.59   20.26 -0.7501E+01  0.9459E+01
   25   24   17    1.59   33.23 -0.1230E+02  0.8612E+01
   23   23   18    1.60   11.64 -0.9665E+01  0.1040E+02
   36   27    9    1.60   25.61 -0.9482E+01  0.9559E+01
   43    2   11    1.60   12.83 -0.4750E+01  0.9677E+01
   20   13   21    1.60   32.36 -0.1198E+02  0.8249E+01
   38   28    5    1.60   19.02 -0.7042E+01  0.9141E+01
   27   24   16    1.60   56.27 -0.2083E+02  0.1431E+02
   46    1    7    1.61   27.72 -0.1026E+02  0.8616E+01
   47    1    5    1.61   42.88 -0.1588E+02  0.8982E+01
   33    5   17    1.62   16.82 -0.6228E+01  0.7041E+01
   31   25   13    1.63   43.81 -0.1622E+02  0.8169E+01
   19    9   21    1.65   15.95 -0.5568E+01  0.6454E+01
   33   25   11    1.65   55.07 -0.1922E+02  0.1375E+02
   43    1    9    1.66   14.06 -0.4908E+01  0.8575E+01
   37    3   14    1.66   14.26 -0.4977E+01  0.7608E+01
   27   23   15    1.67   16.20 -0.5654E+01  0.7679E+01
   37   26    5    1.67   34.21 -0.1194E+02  0.8645E+01
   45    0    5    1.68   10.33 -0.7420E+01  0.8669E+01
   36   26    6    1.68   49.41 -0.1724E+02  0.1489E+02
   33   25   10    1.68   50.90 -0.1776E+02  0.1309E+02
   31   24   12    1.68   13.21 -0.4612E+01  0.7851E+01
   35   25    8    1.69   21.19 -0.7397E+01  0.1408E+02
   39    2   12    1.69   42.44 -0.1481E+02  0.7959E+01
   32   24   11    1.69   43.57 -0.1521E+02  0.8713E+01
   37   26    3    1.70   11.90 -0.4154E+01  0.9178E+01
   37    2   13    1.71   44.84 -0.6328E+01  0.8527E+01
   35   25    7    1.71   98.10 -0.1384E+02  0.8687E+01
   29   23   13    1.71   18.66 -0.2633E+01  0.7774E+01
   44    0    5    1.71   19.78 -0.1142E+02  0.8413E+01
   36   25    5    1.72  142.70 -0.2014E+02  0.8636E+01
   17   11   20    1.73   18.62 -0.2628E+01  0.6458E+01
   30   23   12    1.73   29.25 -0.4128E+01  0.8298E+01
   40    1   10    1.73   88.17 -0.1244E+02  0.8118E+01
   31    4   16    1.73  108.70 -0.1534E+02  0.1243E+02
   27   22   14    1.73   31.72 -0.4477E+01  0.7353E+01
   34   24    8    1.74   12.26 -0.1730E+01  0.8417E+01
   19   18   18    1.74   21.26 -0.3000E+01  0.6282E+01
   36    2   13    1.74  109.66 -0.1547E+02  0.7569E+01
   42    0    7    1.75   13.08 -0.4562E+01  0.8480E+01
   32   23   10    1.75   96.39 -0.1360E+02  0.1560E+02
   44    1    0    1.76   22.00 -0.1528E+02  0.1283E+02
   16    8   20    1.77   38.67 -0.5457E+01  0.5586E+01
   15    8   20    1.78   78.62 -0.1110E+02  0.5813E+01
   31   22   10    1.80   12.81  0.1556E+02  0.1430E+02
   26   20   13    1.85   11.81  0.1435E+02  0.7720E+01
   33   22    5    1.89   22.14 -0.1919E+02  0.1511E+02
   32   21    6    1.93   15.58 -0.1351E+02  0.8607E+01
   30    1   13    1.95   13.85 -0.1201E+02  0.1200E+02
   28    2   14    1.95   12.51 -0.1085E+02  0.1154E+02
   30   20    8    1.96   15.68 -0.1359E+02  0.8636E+01
   36    1    8    1.96   13.52 -0.1172E+02  0.8527E+01
   21   17   14    1.99   19.68 -0.1706E+02  0.8117E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        85
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1405
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1228
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      177


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
