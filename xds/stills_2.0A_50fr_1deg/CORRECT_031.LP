
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.543    77.543    38.806  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      31      31
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      31
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000959  0.000684
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
 REFLECTING_RANGE_E.S.D.=   0.147
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

    1284 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      31 ...      31
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1284 REFLECTIONS ACCEPTED



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
 *  33        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  32        oP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.9   38.8   77.5  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.9   38.8   77.5  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    28        mC        250.0      38.8  159.9   77.5  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  159.9   77.5  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    36        oC        250.0      38.8  159.9   77.5  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.9   77.5  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    27        mC        500.0     159.9   38.8  109.7  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.7  116.3  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  159.9  159.9  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.7  116.3   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        998.2     173.4   77.5   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.2      77.5  173.4   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.2      77.5  173.4   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.2      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.3  116.3  116.3  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.7   86.7  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.7   86.7   86.7  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.7   86.7  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.7   86.7  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.7  248.3 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.7   86.7   77.5  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.5  177.7  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.7   86.7  177.7 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.7   86.7   86.7  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.7  109.7   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.5  177.7 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.7  177.7   77.5 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    4.5      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.895   -23.505     0.058
 UNIT_CELL_B-AXIS=    23.505    73.895     0.019
 UNIT_CELL_C-AXIS=    -0.030    -0.000    38.807
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       173 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.27
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.438    77.438    38.755  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.3E-01 1.3E-01 6.2E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012914  0.012914  0.025803  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.800   -23.455     0.068
 COORDINATES OF UNIT CELL B-AXIS    23.455    73.800     0.023
 COORDINATES OF UNIT CELL C-AXIS    -0.036    -0.001    38.755
 CRYSTAL MOSAICITY (DEGREES)     0.147
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006541  0.002003
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.90   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1260.99   1242.40
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.04
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.44    77.44    38.75  90.000  90.000  90.000


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
  0.05  0.15        51         98.445       1.221
  0.15  0.25       166        127.897       0.940
  0.25  0.35       178        124.101       0.969
  0.35  0.45       160         51.979       2.312
  0.45  0.55       162        220.907       0.544
  0.55  0.65       182         83.675       1.436
  0.65  0.75       201        130.915       0.918
  0.75  0.85       141        116.428       1.032
  0.85  0.95        43         91.387       1.315
 
 -0.05  0.95      1284        120.198       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        8
 DEGREES OF FREEDOM OF CHI^2 FIT               211.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.966
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    30.1 XMAX=    30.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59473 NYBIN=    8
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        442


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               211.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.966
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    11.3 XMAX=  2452.8 NXBIN=    3
 YMIN=     8.7 YMAX=  2519.3 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        442


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               211.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.963
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    30.1 XMAX=    30.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1647    1690
 DETECTOR_SURFACE_POSITION=     817    1690
 DETECTOR_SURFACE_POSITION=     817     838
 DETECTOR_SURFACE_POSITION=    1647     838
 DETECTOR_SURFACE_POSITION=    2170    1664
 DETECTOR_SURFACE_POSITION=    1621    2229
 DETECTOR_SURFACE_POSITION=     843    2229
 DETECTOR_SURFACE_POSITION=     294    1664
 DETECTOR_SURFACE_POSITION=     294     864
 DETECTOR_SURFACE_POSITION=     843     299
 DETECTOR_SURFACE_POSITION=    1621     299
 DETECTOR_SURFACE_POSITION=    2170     864
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        442


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
 8.892E-01  1.489E-03   27.48


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

   24.488   3.154    19.89   1.68      2.33      2.26      68      89      18
    3.154   2.240     7.17   1.24      7.80      6.35     188     213       2
    2.240   1.831     0.80   0.80     57.51     62.69     198     220       0
    1.831   1.587     0.00   0.80   1821.18   1897.77     236     264       2
    1.587   1.420    -0.08   0.80    -99.90    -99.90     254     303       4
    1.420   1.297     0.03   0.75   3449.22   3950.68     120     162       6
   --------------------------------------------------------------------------
   24.488   1.297     2.76   0.93      8.43      8.53    1064    1251      32
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      31      31)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.683  -0.1584E+02  0.1290E+02        -1.23           1 
   47    0    0     1.648  -0.5887E+01  0.1266E+02        -0.47           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -4.6% OF MEAN INTENSITY


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

     3.85          61      37      1285        2.9%       2.2%      2.2%       48   26.63      3.2%    99.9*     0    0.000       0
     2.74          84      48      2098        2.3%       5.2%      3.7%       72   17.18      7.4%    99.3*     0    0.000       0
     2.24         157      89      2665        3.3%       9.0%      8.4%      136    8.04     12.7%    98.1*     0    0.000       0
     1.94         142      78      3101        2.5%      39.2%     45.7%      128    1.51     55.4%    86.3*     0    0.000       0
     1.74         156      85      3479        2.4%     481.9%    467.0%      142    0.15    681.6%     0.1      0    0.000       0
     1.59         186     104      3850        2.7%    1247.9%   1342.7%      164   -0.03   1764.9%     1.1      0    0.000       0
     1.47         192     113      4143        2.7%     -99.9%    -99.9%      158    0.00    -99.9%     4.6      0    0.000       0
     1.38         192     111      4428        2.5%    2582.0%   2851.3%      162    0.01   3651.6%     0.3      0    0.000       0
     1.30          81      54      4744        1.1%     -99.9%    -99.9%       54    0.00    -99.9%    11.4      0    0.000       0
    total        1251     719     29793        2.4%       8.4%      8.5%     1064    3.66     11.9%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1284
 NUMBER OF REJECTED MISFITS                                32
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1251
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    719



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

   WILSON LINE (using all data) : A=   7.254 B=  35.555 CORRELATION=  0.85
      #      RES      SS        <I>       log(<I>)       BO
      15     7.256   0.005  1.8386E+03       7.517     -27.7
      18     4.721   0.011  1.8052E+03       7.498     -10.9
      20     3.571   0.020  7.5212E+02       6.623      16.1
      26     2.990   0.028  1.9842E+02       5.290      35.1
      33     2.638   0.036  1.8097E+02       5.198      28.6
      49     2.401   0.043  8.0755E+01       4.391      33.0
      35     2.209   0.051  3.3671E+01       3.517      36.5
      39     2.055   0.059  1.3030E+01       2.567      39.6
      38     1.934   0.067 -1.3338E+00    -999.999   -1000.0
      46     1.829   0.075  2.5037E+00       0.918      42.4
      42     1.739   0.083  3.4296E-01      -1.070      50.4
      42     1.667   0.090  1.4100E-01      -1.959      51.2
      50     1.598   0.098 -2.9873E-01    -999.999   -1000.0
      54     1.538   0.106 -2.1559E+00    -999.999   -1000.0
      66     1.482   0.114 -1.9065E+00    -999.999   -1000.0
      50     1.435   0.121 -8.7845E-02    -999.999   -1000.0
      43     1.392   0.129 -1.1934E-01    -999.999   -1000.0
      37     1.352   0.137  3.4511E+00       1.239      22.0
      16     1.315   0.145 -3.8354E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.256        2.072        1.562        0.921
       3     4.721        1.532        1.033        0.518
       2     3.571        0.335        0.068        0.010
       3     2.990        2.162        1.185        0.481
       2     2.638        0.382        0.116        0.025
       5     2.401        0.677        0.337        0.132
       1     2.209        0.773        0.235        0.051
       6     2.055        0.557       -0.003        0.061
       6     1.934       23.431       49.032      128.137
       9     1.829        4.834        5.281        4.415
       7     1.739      256.943     -726.275    13591.716
       1     1.392       33.861       68.255       98.275
       1     1.352        9.690       10.450        8.049
      75   overall      144.285     -773.244    20878.386


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.256        0.550        0.420        0.255
      15     4.721        1.293        1.700        1.777
      18     3.571        0.960        0.897        0.762
      23     2.990        0.775        0.653        0.456
      31     2.638        1.189        1.460        1.690
      44     2.401        0.949        0.809        0.584
      34     2.209        0.921        0.755        0.559
      33     2.055        2.000        2.634        3.356
      37     1.829        5.758        3.355       10.600
      31     1.667     1727.982    16398.099  1096416.877
      47     1.598       42.336      -75.323     1146.233
      47     1.435       43.225     -112.479     1125.818
      36     1.352        3.363        2.855        4.469
     644   overall       94.303      809.499    53555.195

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.256  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.73
       3     4.721  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     3.571  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.73
       3     2.990  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.638  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.401  0.00  0.00  0.00  0.42  0.77  1.07  1.00  0.95  0.91  0.88
       1     2.209  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.055  1.34  1.45  1.20  1.06  1.28  1.19  1.12  1.06  1.01  0.98
       6     1.934  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.27  1.22
       9     1.829  0.89  0.64  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       7     1.739  0.58  0.41  0.34  0.30  0.27  0.25  0.24  0.23  0.22  0.21
       1     1.392  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       1     1.352  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      75   overall  0.64  0.50  0.45  0.42  0.44  0.43  0.40  0.38  0.41  0.43


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.256  2.87  3.01  2.10  1.93  1.85  1.61  1.63  1.49  1.38  1.29
      15     4.721  2.80  1.84  1.54  1.21  1.02  1.63  1.46  1.33  1.46  1.37
      18     3.571  1.17  0.61  0.64  0.67  0.71  0.86  1.10  1.11  1.12  1.14
      23     2.990  1.37  1.92  1.51  1.58  1.55  1.35  1.30  1.34  1.25  1.17
      31     2.638  0.34  1.07  1.24  1.17  0.98  0.93  1.03  1.05  1.14  1.12
      44     2.401  0.48  0.75  0.96  1.10  1.16  1.06  1.04  1.11  1.07  1.04
      34     2.209  1.55  0.97  0.91  0.80  0.82  0.98  0.93  0.91  0.89  0.93
      33     2.055  2.87  1.84  1.29  1.10  0.92  1.01  0.90  0.88  0.87  0.91
      37     1.829  5.40  2.83  1.98  1.56  1.31  1.20  1.07  0.98  0.91  0.86
      31     1.667  5.42  2.85  1.99  1.57  1.31  1.14  1.03  0.94  0.87  0.82
      47     1.598  5.59  2.93  2.05  1.61  1.35  1.18  1.06  0.97  0.90  0.84
      47     1.435  6.26  3.29  2.38  1.87  1.57  1.37  1.23  1.12  1.04  0.98
      36     1.352  3.79  1.99  1.50  1.26  1.13  0.99  0.94  0.86  0.80  0.83
     644   overall  2.12  1.29  1.00  0.86  0.76  0.72  0.68  0.66  0.63  0.62


 List of     89 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   53    3    6    1.42  328.11 -0.2882E+02  0.1141E+02 "alien"
   48    5   13    1.41  234.78 -0.2062E+02  0.1506E+02 "alien"
   39    8   19    1.41  215.35 -0.2570E+02  0.1288E+02 "alien"
   28   26   19    1.44  188.83 -0.1659E+02  0.1322E+02 "alien"
   49    4   11    1.44  173.94 -0.1528E+02  0.9540E+01 "alien"
   39   33   11    1.39  171.13 -0.2042E+02  0.1019E+02 "alien"
   31   29   16    1.46  165.14 -0.1451E+02  0.9019E+01 "alien"
   39   34   11    1.38  148.22 -0.1769E+02  0.1754E+02 "alien"
   46    7   16    1.37  137.55 -0.1642E+02  0.9245E+01 "alien"
   30   28   19    1.39  130.16 -0.1553E+02  0.8885E+01 "alien"
   42    1    9    1.69  117.52  0.1657E+02  0.8763E+01 "alien"
   45    5   14    1.45  117.35 -0.1031E+02  0.8954E+01 "alien"
   37   27    6    1.64  114.97  0.1621E+02  0.1319E+02 "alien"
   51    4    9    1.43  112.78 -0.9908E+01  0.9768E+01 "alien"
   26   13   24    1.38  103.89 -0.1240E+02  0.1031E+02 "alien"
   40   32    7    1.46  103.39 -0.9082E+01  0.1158E+02 "alien"
   35   10   21    1.39   99.26 -0.1185E+02  0.1296E+02 "alien"
   42   32    3    1.46   93.20 -0.8187E+01  0.1172E+02 "alien"
   47    1    6    1.60   90.29 -0.2697E+02  0.1234E+02 "alien"
   44    1    8    1.65   86.79  0.1224E+02  0.8653E+01 "alien"
   40    2   12    1.66   81.89  0.1155E+02  0.8337E+01 "alien"
   34   25    9    1.69   75.59  0.1066E+02  0.8446E+01 "alien"
   39    2   12    1.69   73.13  0.1031E+02  0.8212E+01 "alien"
   39    7   18    1.45   70.58 -0.6200E+01  0.8337E+01 "alien"
   38   28    6    1.59   70.58 -0.2108E+02  0.1242E+02 "alien"
   52    4    8    1.42   70.30 -0.6175E+01  0.1044E+02 "alien"
   39   28    4    1.59   68.73 -0.2053E+02  0.1099E+02 "alien"
   37   32   12    1.42   66.08 -0.5805E+01  0.1091E+02 "alien"
   43    1    9    1.66   60.08  0.8472E+01  0.8760E+01 "alien"
   53    3    5    1.43   58.79 -0.5164E+01  0.9560E+01 "alien"
   38    3   14    1.64   56.25  0.7932E+01  0.8207E+01 "alien"
   34   26   10    1.64   55.85  0.7875E+01  0.9306E+01 "alien"
   36   10   21    1.38   55.28 -0.6597E+01  0.8428E+01 "alien"
   42    6   16    1.46   52.77 -0.4635E+01  0.9452E+01 "alien"
   52    4    9    1.40   52.62 -0.6279E+01  0.9821E+01 "alien"
   32   24   11    1.70   51.57  0.7271E+01  0.8814E+01 "alien"
   35    2   13    1.77   46.01  0.1578E+02  0.8505E+01 "alien"
   48    6   14    1.39   45.91 -0.5480E+01  0.1360E+02 "alien"
   53    5    9    1.38   45.85 -0.5472E+01  0.1825E+02 "alien"
   43    7   17    1.40   45.38 -0.5415E+01  0.9000E+01 "alien"
   32   25   12    1.64   39.81  0.5613E+01  0.8526E+01 "alien"
   54    4    6    1.40   39.65 -0.4732E+01  0.1009E+02 "alien"
   28   24   15    1.63   39.40 -0.1177E+02  0.7893E+01 "alien"
   37    3   14    1.67   38.36  0.5409E+01  0.7626E+01 "alien"
   50    5   12    1.39   37.08 -0.4425E+01  0.1013E+02 "alien"
   39   28    2    1.61   34.54 -0.1032E+02  0.8711E+01 "alien"
   51    4   10    1.41   33.07 -0.3947E+01  0.1215E+02 "alien"
   39   28    1    1.61   32.15 -0.9603E+01  0.8987E+01 "alien"
   37   33   13    1.38   30.66 -0.3659E+01  0.9701E+01 "alien"
   36   25    5    1.72   29.86  0.1024E+02  0.8294E+01 "alien"
   17   11   20    1.73   28.49  0.9773E+01  0.7489E+01 "alien"
   40    3   13    1.62   27.68 -0.8269E+01  0.8820E+01 "alien"
   37   25    2    1.73   26.96  0.9246E+01  0.9853E+01 "alien"
   38   28    5    1.61   26.46 -0.7903E+01  0.1029E+02 "alien"
   27   18   23    1.38   26.43 -0.3154E+01  0.8089E+01 "alien"
   29   23   13    1.71   26.37  0.9042E+01  0.7756E+01 "alien"
   51    5   11    1.39   26.09 -0.3113E+01  0.9539E+01 "alien"
   46    1    6    1.63   24.57 -0.7341E+01  0.1241E+02 "alien"
   34   24    8    1.74   23.89  0.8192E+01  0.8796E+01 "alien"
   21   17   20    1.61   23.55 -0.7036E+01  0.9131E+01 "alien"
   28   22   13    1.76   21.99  0.7542E+01  0.7648E+01 "alien"
   54    3    0    1.43   21.99 -0.9490E+01  0.1433E+02 "alien"
   36   31   13    1.43   21.46 -0.1886E+01  0.9825E+01 "alien"
   38   33   12    1.39   21.08 -0.2516E+01  0.9884E+01 "alien"
   44    0    6    1.70   19.81  0.1147E+02  0.8896E+01 "alien"
   28   28   18    1.45   17.98 -0.5603E+01  0.9956E+01 "alien"
   33   24    9    1.74   17.34  0.5945E+01  0.8589E+01 "alien"
   40   29    0    1.57   17.31 -0.1742E+02  0.1273E+02 "alien"
   23   23   18    1.60   16.70 -0.1604E+02  0.1077E+02 "alien"
   26   22   15    1.71   16.06  0.5507E+01  0.7919E+01 "alien"
   38    2    4    1.99   15.17 -0.2023E+02  0.1082E+02 "alien"
   33    6   18    1.57   14.71 -0.4393E+01  0.7186E+01 "alien"
   25   21   15    1.75   14.15  0.4853E+01  0.8710E+01 "alien"
   26   21   21    1.44   13.91 -0.1222E+01  0.8102E+01 "alien"
   39    4   15    1.57   13.41 -0.4006E+01  0.1277E+02 "alien"
   44    1    0    1.76   13.23  0.1132E+02  0.1304E+02 "alien"
   22   18   20    1.58   13.22 -0.3949E+01  0.7014E+01 "alien"
   22   13   24    1.43   13.20 -0.1160E+01  0.9313E+01 "alien"
   27   22   14    1.73   13.05  0.4476E+01  0.7490E+01 "alien"
   53    4    8    1.40   12.60 -0.1504E+01  0.9786E+01 "alien"
   45    0    4    1.69   12.59  0.4241E+01  0.8599E+01 "alien"
   22   21   17    1.70   11.78  0.1662E+01  0.9508E+01 "alien"
   48    3    9    1.51   11.68 -0.2226E+02  0.1149E+02 "alien"
   37   25    0    1.73   11.62  0.8924E+01  0.1245E+02 "alien"
   43    0    7    1.71   11.59  0.8888E+01  0.1016E+02 "alien"
   34   22    1    1.91   11.04 -0.1473E+02  0.9244E+01 "alien"
   36    1    8    1.97   10.78 -0.1437E+02  0.8827E+01 "alien"
   39    2    2    1.97   10.30 -0.1374E+02  0.9023E+01 "alien"
   30    2   14    1.89   10.06 -0.1342E+02  0.7187E+01 "alien"


 List of     89 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   46    7   16    1.37  137.55 -0.1642E+02  0.9245E+01
   27   18   23    1.38   26.43 -0.3154E+01  0.8089E+01
   39   34   11    1.38  148.22 -0.1769E+02  0.1754E+02
   53    5    9    1.38   45.85 -0.5472E+01  0.1825E+02
   36   10   21    1.38   55.28 -0.6597E+01  0.8428E+01
   26   13   24    1.38  103.89 -0.1240E+02  0.1031E+02
   37   33   13    1.38   30.66 -0.3659E+01  0.9701E+01
   30   28   19    1.39  130.16 -0.1553E+02  0.8885E+01
   48    6   14    1.39   45.91 -0.5480E+01  0.1360E+02
   51    5   11    1.39   26.09 -0.3113E+01  0.9539E+01
   38   33   12    1.39   21.08 -0.2516E+01  0.9884E+01
   50    5   12    1.39   37.08 -0.4425E+01  0.1013E+02
   39   33   11    1.39  171.13 -0.2042E+02  0.1019E+02
   35   10   21    1.39   99.26 -0.1185E+02  0.1296E+02
   53    4    8    1.40   12.60 -0.1504E+01  0.9786E+01
   54    4    6    1.40   39.65 -0.4732E+01  0.1009E+02
   43    7   17    1.40   45.38 -0.5415E+01  0.9000E+01
   52    4    9    1.40   52.62 -0.6279E+01  0.9821E+01
   39    8   19    1.41  215.35 -0.2570E+02  0.1288E+02
   51    4   10    1.41   33.07 -0.3947E+01  0.1215E+02
   48    5   13    1.41  234.78 -0.2062E+02  0.1506E+02
   52    4    8    1.42   70.30 -0.6175E+01  0.1044E+02
   37   32   12    1.42   66.08 -0.5805E+01  0.1091E+02
   53    3    6    1.42  328.11 -0.2882E+02  0.1141E+02
   22   13   24    1.43   13.20 -0.1160E+01  0.9313E+01
   51    4    9    1.43  112.78 -0.9908E+01  0.9768E+01
   36   31   13    1.43   21.46 -0.1886E+01  0.9825E+01
   54    3    0    1.43   21.99 -0.9490E+01  0.1433E+02
   53    3    5    1.43   58.79 -0.5164E+01  0.9560E+01
   28   26   19    1.44  188.83 -0.1659E+02  0.1322E+02
   49    4   11    1.44  173.94 -0.1528E+02  0.9540E+01
   26   21   21    1.44   13.91 -0.1222E+01  0.8102E+01
   39    7   18    1.45   70.58 -0.6200E+01  0.8337E+01
   28   28   18    1.45   17.98 -0.5603E+01  0.9956E+01
   45    5   14    1.45  117.35 -0.1031E+02  0.8954E+01
   31   29   16    1.46  165.14 -0.1451E+02  0.9019E+01
   42   32    3    1.46   93.20 -0.8187E+01  0.1172E+02
   42    6   16    1.46   52.77 -0.4635E+01  0.9452E+01
   40   32    7    1.46  103.39 -0.9082E+01  0.1158E+02
   48    3    9    1.51   11.68 -0.2226E+02  0.1149E+02
   40   29    0    1.57   17.31 -0.1742E+02  0.1273E+02
   39    4   15    1.57   13.41 -0.4006E+01  0.1277E+02
   33    6   18    1.57   14.71 -0.4393E+01  0.7186E+01
   22   18   20    1.58   13.22 -0.3949E+01  0.7014E+01
   38   28    6    1.59   70.58 -0.2108E+02  0.1242E+02
   39   28    4    1.59   68.73 -0.2053E+02  0.1099E+02
   47    1    6    1.60   90.29 -0.2697E+02  0.1234E+02
   23   23   18    1.60   16.70 -0.1604E+02  0.1077E+02
   38   28    5    1.61   26.46 -0.7903E+01  0.1029E+02
   21   17   20    1.61   23.55 -0.7036E+01  0.9131E+01
   39   28    2    1.61   34.54 -0.1032E+02  0.8711E+01
   39   28    1    1.61   32.15 -0.9603E+01  0.8987E+01
   40    3   13    1.62   27.68 -0.8269E+01  0.8820E+01
   46    1    6    1.63   24.57 -0.7341E+01  0.1241E+02
   28   24   15    1.63   39.40 -0.1177E+02  0.7893E+01
   37   27    6    1.64  114.97  0.1621E+02  0.1319E+02
   38    3   14    1.64   56.25  0.7932E+01  0.8207E+01
   34   26   10    1.64   55.85  0.7875E+01  0.9306E+01
   32   25   12    1.64   39.81  0.5613E+01  0.8526E+01
   44    1    8    1.65   86.79  0.1224E+02  0.8653E+01
   40    2   12    1.66   81.89  0.1155E+02  0.8337E+01
   43    1    9    1.66   60.08  0.8472E+01  0.8760E+01
   37    3   14    1.67   38.36  0.5409E+01  0.7626E+01
   34   25    9    1.69   75.59  0.1066E+02  0.8446E+01
   39    2   12    1.69   73.13  0.1031E+02  0.8212E+01
   45    0    4    1.69   12.59  0.4241E+01  0.8599E+01
   42    1    9    1.69  117.52  0.1657E+02  0.8763E+01
   32   24   11    1.70   51.57  0.7271E+01  0.8814E+01
   44    0    6    1.70   19.81  0.1147E+02  0.8896E+01
   22   21   17    1.70   11.78  0.1662E+01  0.9508E+01
   26   22   15    1.71   16.06  0.5507E+01  0.7919E+01
   29   23   13    1.71   26.37  0.9042E+01  0.7756E+01
   43    0    7    1.71   11.59  0.8888E+01  0.1016E+02
   36   25    5    1.72   29.86  0.1024E+02  0.8294E+01
   37   25    2    1.73   26.96  0.9246E+01  0.9853E+01
   17   11   20    1.73   28.49  0.9773E+01  0.7489E+01
   27   22   14    1.73   13.05  0.4476E+01  0.7490E+01
   37   25    0    1.73   11.62  0.8924E+01  0.1245E+02
   33   24    9    1.74   17.34  0.5945E+01  0.8589E+01
   34   24    8    1.74   23.89  0.8192E+01  0.8796E+01
   25   21   15    1.75   14.15  0.4853E+01  0.8710E+01
   28   22   13    1.76   21.99  0.7542E+01  0.7648E+01
   44    1    0    1.76   13.23  0.1132E+02  0.1304E+02
   35    2   13    1.77   46.01  0.1578E+02  0.8505E+01
   30    2   14    1.89   10.06 -0.1342E+02  0.7187E+01
   34   22    1    1.91   11.04 -0.1473E+02  0.9244E+01
   36    1    8    1.97   10.78 -0.1437E+02  0.8827E+01
   39    2    2    1.97   10.30 -0.1374E+02  0.9023E+01
   38    2    4    1.99   15.17 -0.2023E+02  0.1082E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        64
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1284
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1139
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      145


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
