
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.469    77.469    38.766  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      22      22
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      22
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000873  0.000701
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
 REFLECTING_RANGE_E.S.D.=   0.117
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

    1394 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      22 ...      22
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1394 REFLECTIONS ACCEPTED



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
 *  25        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  31        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  20        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  35        mP          0.0      77.5   38.8   77.5  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  33        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  21        tP          0.0      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
 *  32        oP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    37        mC        250.0     159.7   38.8   77.5  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.7   38.8   77.5  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  159.7   77.5  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.7   77.5  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  159.7   77.5  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.8  159.7   77.5  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    19        oI        500.0      38.8  109.6  116.2  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    27        mC        500.0     159.7   38.8  109.6  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    18        tI        625.0     109.6  116.2   38.8  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    26        oF        625.0      38.8  159.7  159.7  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    22        hP        998.4      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    30        mC        998.4      77.5  173.2   38.8  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.4      77.5  173.2   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    41        mC        998.4     173.2   77.5   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.2  116.2  116.2  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.6   86.6  116.2 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.6   86.6  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.6   86.6  116.2 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.6   86.6   86.6  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.6   86.6  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.6   86.6  109.6  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.6  248.0 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.6   86.6   77.5  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
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
 UNIT_CELL_A-AXIS=    73.789   -23.594    -0.016
 UNIT_CELL_B-AXIS=    23.594    73.789    -0.003
 UNIT_CELL_C-AXIS=     0.008    -0.001    38.765
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       239 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.26
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.452    77.452    38.756  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.3E-01 1.3E-01 5.6E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012911  0.012911  0.025802  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.810   -23.471     0.032
 COORDINATES OF UNIT CELL B-AXIS    23.471    73.810     0.009
 COORDINATES OF UNIT CELL C-AXIS    -0.016     0.001    38.756
 CRYSTAL MOSAICITY (DEGREES)     0.117
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006511  0.001789
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.91   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1261.31   1242.36
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.33
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.45    77.45    38.76  90.000  90.000  90.000


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
  0.05  0.15        83         89.143       1.784
  0.15  0.25       169        194.512       0.818
  0.25  0.35       194        140.154       1.135
  0.35  0.45       146         89.707       1.773
  0.45  0.55       171        221.413       0.718
  0.55  0.65       192        141.648       1.123
  0.65  0.75       194        176.133       0.903
  0.75  0.85       165        122.517       1.298
  0.85  0.95        80        271.104       0.587
 
 -0.05  0.95      1394        159.031       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               237.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.897
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    21.1 XMAX=    21.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59458 NYBIN=    9
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        496


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               237.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.895
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     9.1 XMAX=  2453.8 NXBIN=    3
 YMIN=     7.6 YMAX=  2517.4 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        496


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               237.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.884
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    21.1 XMAX=    21.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1231    1263
 DETECTOR_SURFACE_POSITION=    1647    1689
 DETECTOR_SURFACE_POSITION=     816    1689
 DETECTOR_SURFACE_POSITION=     816     836
 DETECTOR_SURFACE_POSITION=    1647     836
 DETECTOR_SURFACE_POSITION=    2171    1662
 DETECTOR_SURFACE_POSITION=    1621    2227
 DETECTOR_SURFACE_POSITION=     842    2227
 DETECTOR_SURFACE_POSITION=     292    1662
 DETECTOR_SURFACE_POSITION=     292     863
 DETECTOR_SURFACE_POSITION=     842     298
 DETECTOR_SURFACE_POSITION=    1621     298
 DETECTOR_SURFACE_POSITION=    2171     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        496


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
 1.188E+00  3.017E-03   16.71


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

   24.492   3.155    14.21   1.34      3.30      3.53      84     110      14
    3.155   2.240     7.72   1.06      6.28      5.79     198     227       2
    2.240   1.831     1.46   0.78     27.97     30.31     206     241       2
    1.831   1.587    -0.06   0.76    -99.90    -99.90     242     282       0
    1.587   1.420    -0.01   0.61   1452.53   1860.76     286     339       0
    1.420   1.297     0.00   0.72    -99.90    -99.90     132     174       2
   --------------------------------------------------------------------------
   24.492   1.297     2.66   0.82      7.74      8.46    1148    1373      20
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      22      22)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.684  -0.3934E+01  0.1401E+02        -0.28           1 
   47    0    0     1.648   0.1590E+02  0.1439E+02         1.11           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     9.3% OF MEAN INTENSITY


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

     3.85          69      43      1285        3.3%       3.0%      3.5%       52   18.36      4.2%    99.8*     0    0.000       0
     2.74         104      61      2098        2.9%       4.7%      4.2%       86   14.52      6.6%    98.9*     0    0.000       0
     2.24         164      92      2665        3.5%       7.5%      6.6%      144    9.36     10.6%    98.3*     0    0.000       0
     1.94         155      88      3101        2.8%      21.4%     20.6%      134    2.97     30.2%    93.3*     0    0.000       0
     1.74         167      96      3480        2.8%     -99.9%    -99.9%      142    0.00    -99.9%     6.4      0    0.000       0
     1.59         201     115      3850        3.0%     -99.9%    -99.9%      172    0.00    -99.9%     0.2      0    0.000       0
     1.47         218     128      4143        3.1%    2413.9%   3092.0%      180   -0.06   3413.8%    14.7      0    0.000       0
     1.38         210     117      4428        2.6%    1084.6%   1405.2%      186    0.04   1533.9%    -2.3      0    0.000       0
     1.30          85      59      4742        1.2%     -99.9%    -99.9%       52    0.00    -99.9%     1.9      0    0.000       0
    total        1373     799     29792        2.7%       7.7%      8.5%     1148    3.48     11.0%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1394
 NUMBER OF REJECTED MISFITS                                20
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1373
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    799



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

   WILSON LINE (using all data) : A=   7.271 B=  31.732 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      19     7.301   0.005  1.8171E+03       7.505     -24.9
      20     4.579   0.012  2.1950E+03       7.694     -17.7
      26     3.520   0.020  9.5454E+02       6.861      10.2
      31     2.999   0.028  3.1946E+02       5.767      27.1
      37     2.644   0.036  2.9041E+02       5.671      22.4
      49     2.402   0.043  1.5834E+02       5.065      25.5
      38     2.210   0.051  8.5600E+01       4.450      27.6
      43     2.057   0.059  3.2173E+01       3.471      32.2
      44     1.934   0.067  6.9723E-01      -0.361      57.1
      51     1.831   0.075  8.1841E-01      -0.200      50.1
      46     1.742   0.082 -1.5985E+00    -999.999   -1000.0
      51     1.668   0.090  8.1019E-01      -0.210      41.6
      53     1.598   0.098 -5.3436E-01    -999.999   -1000.0
      65     1.537   0.106 -2.8738E+00    -999.999   -1000.0
      70     1.482   0.114  1.4769E+00       0.390      30.2
      55     1.435   0.121  1.0174E+00       0.017      29.9
      43     1.392   0.129  3.5947E-01      -1.023      32.1
      40     1.353   0.137 -1.3498E+00    -999.999   -1000.0
      18     1.316   0.144  2.4616E+00       0.901      22.1


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.301        2.192        1.626        0.964
       3     4.579        1.483        0.987        0.487
       2     3.520        0.394        0.086        0.013
       4     2.999        1.342        0.658        0.243
       2     2.644        0.221        0.051        0.008
       5     2.402        0.802        0.385        0.143
       2     2.210        1.180        0.480        0.144
       7     2.057        0.314        0.113        0.032
       7     1.934       68.266      329.810     1532.261
      10     1.831       29.446       16.602      121.791
       3     1.598       71.213       27.072      576.996
      81   overall       15.065       37.743      182.191


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      14     7.301        0.475        0.364        0.226
      17     4.579        1.143        1.399        1.410
      24     3.520        0.840        0.673        0.505
      27     2.999        0.855        0.736        0.540
      35     2.644        1.195        1.412        1.555
      44     2.402        0.905        0.753        0.532
      36     2.210        0.821        0.628        0.419
      36     2.057        1.460        2.082        2.707
      37     1.934       83.275     -107.343     2974.200
      41     1.831       60.359      139.373     2097.220
      40     1.668       70.579       53.110     2377.825
      64     1.482       23.961       37.871      268.904
      52     1.435       25.116       20.941      285.627
      43     1.392      359.925      345.164    67137.300
      18     1.316       15.245       56.148      273.545
     718   overall       37.887       32.698     4478.104

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.301  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.30  0.29
       3     4.579  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.520  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       4     2.999  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       2     2.644  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.402  0.00  0.00  0.00  0.00  0.77  0.71  1.00  0.95  0.91  0.88
       2     2.210  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       7     2.057  1.15  1.24  1.03  1.21  1.10  1.02  1.20  1.14  1.09  1.05
       7     1.934  1.73  1.24  1.03  0.91  0.82  0.76  0.72  0.68  0.65  0.63
      10     1.831  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.59
       3     1.598  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      81   overall  0.65  0.50  0.41  0.39  0.40  0.40  0.41  0.39  0.39  0.40


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      14     7.301  1.50  2.76  1.93  1.73  1.82  1.74  1.70  1.56  1.44  1.36
      17     4.579  1.24  1.62  1.13  1.07  1.05  1.43  1.29  1.28  1.39  1.30
      24     3.520  0.88  0.92  0.64  0.51  0.53  0.83  1.08  0.98  0.98  0.99
      27     2.999  1.17  1.84  1.43  1.35  1.22  1.31  1.18  1.21  1.19  1.17
      35     2.644  0.90  1.26  1.21  1.13  0.94  0.95  1.02  1.09  1.11  1.13
      44     2.402  0.96  0.75  0.96  0.97  1.10  1.01  1.04  1.11  1.11  1.04
      36     2.210  1.17  0.92  0.86  0.84  0.85  0.86  1.05  0.96  0.98  1.10
      36     2.057  0.58  1.07  0.96  1.10  1.20  1.11  1.10  1.06  1.08  1.01
      37     1.934  5.11  2.68  1.88  1.56  1.31  1.14  1.02  0.93  0.87  0.81
      41     1.831  5.38  2.83  1.98  1.55  1.30  1.14  1.02  0.93  0.86  0.81
      40     1.668  5.52  2.90  2.03  1.59  1.40  1.22  1.14  1.04  0.97  0.91
      64     1.482  4.27  2.24  1.57  1.28  1.07  0.94  0.84  0.77  0.71  0.67
      52     1.435  4.85  2.55  1.78  1.40  1.22  1.07  0.96  0.87  0.84  0.79
      43     1.392  5.13  2.69  1.88  1.48  1.24  1.13  1.02  0.93  0.86  0.81
      18     1.316  4.67  2.45  1.71  1.35  1.13  0.99  0.99  0.91  0.84  0.88
     718   overall  2.36  1.47  1.10  0.93  0.85  0.80  0.77  0.73  0.71  0.69


 List of     77 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   35   10   21    1.39   55.56  0.1997E+02  0.1123E+02 "alien"
   52    5   10    1.38   54.35  0.1954E+02  0.1187E+02 "alien"
   39   34   11    1.38   53.70  0.1930E+02  0.1335E+02 "alien"
   38   33   12    1.39   51.20  0.1841E+02  0.1111E+02 "alien"
   27   24   16    1.61   47.30 -0.2528E+02  0.1561E+02 "alien"
   27   20   22    1.40   35.81  0.1287E+02  0.8684E+01 "alien"
   54    4    7    1.38   32.10  0.1154E+02  0.1101E+02 "alien"
   41    3   13    1.59   31.02 -0.1657E+02  0.9591E+01 "alien"
   42    1    2    1.84   29.54  0.2418E+02  0.1081E+02 "alien"
   41    1    3    1.87   28.32  0.2318E+02  0.1863E+02 "alien"
   47    1    5    1.61   27.54 -0.1471E+02  0.9971E+01 "alien"
   35    5   17    1.58   27.07 -0.1446E+02  0.1367E+02 "alien"
   36   33   14    1.38   23.60  0.8485E+01  0.1145E+02 "alien"
   43    1    9    1.66   23.21  0.1881E+02  0.1011E+02 "alien"
   27   19   11    1.95   22.58  0.1575E+02  0.9959E+01 "alien"
   29   20   10    1.91   22.21  0.1549E+02  0.1569E+02 "alien"
   32   25   12    1.64   21.19  0.1717E+02  0.9730E+01 "alien"
   35   26    9    1.64   20.90  0.1694E+02  0.9977E+01 "alien"
   36   26    6    1.68   20.03  0.1623E+02  0.1725E+02 "alien"
   36   10   21    1.38   19.77  0.7107E+01  0.9273E+01 "alien"
   45    2    9    1.60   19.30 -0.1031E+02  0.1102E+02 "alien"
   40   33    9    1.41   18.81  0.6762E+01  0.1081E+02 "alien"
   37    1    7    1.96   18.65  0.1300E+02  0.9723E+01 "alien"
   45    8   17    1.36   18.57 -0.2506E+02  0.1244E+02 "alien"
   52    5   11    1.37   18.47 -0.2494E+02  0.1775E+02 "alien"
   31   30   18    1.38   18.41  0.6616E+01  0.1114E+02 "alien"
   28   24   15    1.63   18.40 -0.9831E+01  0.8775E+01 "alien"
   33    5   17    1.63   18.33 -0.9796E+01  0.8055E+01 "alien"
   53    7   13    1.30   18.32  0.4509E+02  0.1665E+02 "alien"
   28   17   24    1.33   17.83 -0.2407E+02  0.1480E+02 "alien"
   34   25    9    1.69   17.82  0.1444E+02  0.9577E+01 "alien"
   40    2   12    1.66   17.57  0.1424E+02  0.9565E+01 "alien"
   49    5   12    1.41   16.99  0.1728E+02  0.1132E+02 "alien"
   19    9   22    1.59   16.79 -0.8972E+01  0.8097E+01 "alien"
   39    2    3    1.96   16.68  0.1163E+02  0.1212E+02 "alien"
   42    1    9    1.69   16.67  0.1350E+02  0.1000E+02 "alien"
   33    8   20    1.48   16.25  0.2399E+02  0.1417E+02 "alien"
   37   30   11    1.48   16.15  0.2385E+02  0.1303E+02 "alien"
   43    2   11    1.60   15.91 -0.8502E+01  0.1090E+02 "alien"
   29    2   14    1.92   15.54  0.1083E+02  0.8873E+01 "alien"
   51    6   13    1.35   15.54 -0.2097E+02  0.1872E+02 "alien"
   29   12   23    1.39   15.11  0.5430E+01  0.1293E+02 "alien"
   31   22    9    1.84   14.32  0.1172E+02  0.1154E+02 "alien"
   25   10   22    1.50   14.09  0.2080E+02  0.1061E+02 "alien"
   34   29   13    1.50   14.00  0.2068E+02  0.1064E+02 "alien"
   34   22    1    1.91   13.93  0.9711E+01  0.1018E+02 "alien"
   31   22   10    1.80   13.80  0.1130E+02  0.1559E+02 "alien"
   31   30   17    1.41   13.29  0.4776E+01  0.1087E+02 "alien"
   33    0   11    1.95   13.08  0.2253E+02  0.1597E+02 "alien"
   41   33    6    1.43   12.99  0.1322E+02  0.1375E+02 "alien"
   29   27   19    1.41   12.96  0.4660E+01  0.9743E+01 "alien"
   23   23   18    1.60   12.88 -0.1677E+02  0.1179E+02 "alien"
   50    7   14    1.34   12.82 -0.1730E+02  0.1156E+02 "alien"
   43    9   19    1.33   12.81 -0.1729E+02  0.1622E+02 "alien"
   35    1   12    1.83   12.76  0.1044E+02  0.9409E+01 "alien"
   39    7   18    1.45   12.67  0.1289E+02  0.9628E+01 "alien"
   48    1    3    1.60   12.24 -0.6542E+01  0.1102E+02 "alien"
   43    1    2    1.79   12.23  0.1001E+02  0.9901E+01 "alien"
   34   22    2    1.90   12.18  0.8493E+01  0.9845E+01 "alien"
   37   30   10    1.50   12.16  0.1796E+02  0.1169E+02 "alien"
   37   27    6    1.64   12.06  0.9769E+01  0.1424E+02 "alien"
   46    1    7    1.61   12.03 -0.6431E+01  0.9776E+01 "alien"
   36   27    8    1.62   11.70 -0.6252E+01  0.9938E+01 "alien"
   43    6   16    1.44   11.41  0.1161E+02  0.1028E+02 "alien"
   52    3    5    1.46   11.23  0.1659E+02  0.1894E+02 "alien"
   37    1   11    1.80   11.22  0.9186E+01  0.9310E+01 "alien"
   25   18   12    1.98   11.08  0.7728E+01  0.1381E+02 "alien"
   36   24    3    1.77   10.99 -0.1757E+02  0.1094E+02 "alien"
   40   31    6    1.49   10.72  0.1583E+02  0.1068E+02 "alien"
   31   29   16    1.46   10.68  0.1087E+02  0.9979E+01 "alien"
   21   17   14    1.99   10.60  0.7388E+01  0.9224E+01 "alien"
   44    0    6    1.70   10.52  0.1778E+02  0.1020E+02 "alien"
   21   19   19    1.64   10.45  0.8470E+01  0.1177E+02 "alien"
   39   31    9    1.46   10.40  0.1536E+02  0.1206E+02 "alien"
   45    1    7    1.64   10.37  0.8398E+01  0.9812E+01 "alien"
   40    1    6    1.85   10.22  0.8365E+01  0.9846E+01 "alien"
   22   13   24    1.43   10.08  0.1026E+02  0.9450E+01 "alien"


 List of     77 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   18.32  0.4509E+02  0.1665E+02
   28   17   24    1.33   17.83 -0.2407E+02  0.1480E+02
   43    9   19    1.33   12.81 -0.1729E+02  0.1622E+02
   50    7   14    1.34   12.82 -0.1730E+02  0.1156E+02
   51    6   13    1.35   15.54 -0.2097E+02  0.1872E+02
   45    8   17    1.36   18.57 -0.2506E+02  0.1244E+02
   52    5   11    1.37   18.47 -0.2494E+02  0.1775E+02
   36   33   14    1.38   23.60  0.8485E+01  0.1145E+02
   39   34   11    1.38   53.70  0.1930E+02  0.1335E+02
   36   10   21    1.38   19.77  0.7107E+01  0.9273E+01
   31   30   18    1.38   18.41  0.6616E+01  0.1114E+02
   52    5   10    1.38   54.35  0.1954E+02  0.1187E+02
   54    4    7    1.38   32.10  0.1154E+02  0.1101E+02
   38   33   12    1.39   51.20  0.1841E+02  0.1111E+02
   29   12   23    1.39   15.11  0.5430E+01  0.1293E+02
   35   10   21    1.39   55.56  0.1997E+02  0.1123E+02
   27   20   22    1.40   35.81  0.1287E+02  0.8684E+01
   31   30   17    1.41   13.29  0.4776E+01  0.1087E+02
   40   33    9    1.41   18.81  0.6762E+01  0.1081E+02
   29   27   19    1.41   12.96  0.4660E+01  0.9743E+01
   49    5   12    1.41   16.99  0.1728E+02  0.1132E+02
   22   13   24    1.43   10.08  0.1026E+02  0.9450E+01
   41   33    6    1.43   12.99  0.1322E+02  0.1375E+02
   43    6   16    1.44   11.41  0.1161E+02  0.1028E+02
   39    7   18    1.45   12.67  0.1289E+02  0.9628E+01
   31   29   16    1.46   10.68  0.1087E+02  0.9979E+01
   52    3    5    1.46   11.23  0.1659E+02  0.1894E+02
   39   31    9    1.46   10.40  0.1536E+02  0.1206E+02
   37   30   11    1.48   16.15  0.2385E+02  0.1303E+02
   33    8   20    1.48   16.25  0.2399E+02  0.1417E+02
   40   31    6    1.49   10.72  0.1583E+02  0.1068E+02
   34   29   13    1.50   14.00  0.2068E+02  0.1064E+02
   37   30   10    1.50   12.16  0.1796E+02  0.1169E+02
   25   10   22    1.50   14.09  0.2080E+02  0.1061E+02
   35    5   17    1.58   27.07 -0.1446E+02  0.1367E+02
   19    9   22    1.59   16.79 -0.8972E+01  0.8097E+01
   41    3   13    1.59   31.02 -0.1657E+02  0.9591E+01
   45    2    9    1.60   19.30 -0.1031E+02  0.1102E+02
   23   23   18    1.60   12.88 -0.1677E+02  0.1179E+02
   48    1    3    1.60   12.24 -0.6542E+01  0.1102E+02
   43    2   11    1.60   15.91 -0.8502E+01  0.1090E+02
   27   24   16    1.61   47.30 -0.2528E+02  0.1561E+02
   46    1    7    1.61   12.03 -0.6431E+01  0.9776E+01
   47    1    5    1.61   27.54 -0.1471E+02  0.9971E+01
   36   27    8    1.62   11.70 -0.6252E+01  0.9938E+01
   33    5   17    1.63   18.33 -0.9796E+01  0.8055E+01
   28   24   15    1.63   18.40 -0.9831E+01  0.8775E+01
   21   19   19    1.64   10.45  0.8470E+01  0.1177E+02
   37   27    6    1.64   12.06  0.9769E+01  0.1424E+02
   35   26    9    1.64   20.90  0.1694E+02  0.9977E+01
   32   25   12    1.64   21.19  0.1717E+02  0.9730E+01
   45    1    7    1.64   10.37  0.8398E+01  0.9812E+01
   40    2   12    1.66   17.57  0.1424E+02  0.9565E+01
   43    1    9    1.66   23.21  0.1881E+02  0.1011E+02
   36   26    6    1.68   20.03  0.1623E+02  0.1725E+02
   34   25    9    1.69   17.82  0.1444E+02  0.9577E+01
   42    1    9    1.69   16.67  0.1350E+02  0.1000E+02
   44    0    6    1.70   10.52  0.1778E+02  0.1020E+02
   36   24    3    1.77   10.99 -0.1757E+02  0.1094E+02
   43    1    2    1.79   12.23  0.1001E+02  0.9901E+01
   37    1   11    1.80   11.22  0.9186E+01  0.9310E+01
   31   22   10    1.80   13.80  0.1130E+02  0.1559E+02
   35    1   12    1.83   12.76  0.1044E+02  0.9409E+01
   42    1    2    1.84   29.54  0.2418E+02  0.1081E+02
   31   22    9    1.84   14.32  0.1172E+02  0.1154E+02
   40    1    6    1.85   10.22  0.8365E+01  0.9846E+01
   41    1    3    1.87   28.32  0.2318E+02  0.1863E+02
   34   22    2    1.90   12.18  0.8493E+01  0.9845E+01
   34   22    1    1.91   13.93  0.9711E+01  0.1018E+02
   29   20   10    1.91   22.21  0.1549E+02  0.1569E+02
   29    2   14    1.92   15.54  0.1083E+02  0.8873E+01
   27   19   11    1.95   22.58  0.1575E+02  0.9959E+01
   33    0   11    1.95   13.08  0.2253E+02  0.1597E+02
   37    1    7    1.96   18.65  0.1300E+02  0.9723E+01
   39    2    3    1.96   16.68  0.1163E+02  0.1212E+02
   25   18   12    1.98   11.08  0.7728E+01  0.1381E+02
   21   17   14    1.99   10.60  0.7388E+01  0.9224E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        19
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1394
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1341
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       53


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
