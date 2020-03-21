
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.630    77.630    38.867  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      46      46
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      46
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001038  0.000642
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
 BEAM_DIVERGENCE_E.S.D.=   0.050
 REFLECTING_RANGE_E.S.D.=   0.203
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

    1094 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      46 ...      46
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1094 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               36
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

 *  44        aP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  31        aP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0  1  0  0  1  0  0  0
 *  34        mP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  35        mP          0.0      77.6   38.9   77.6  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  32        oP          0.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.8  109.8   38.9  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.6   77.6   38.9  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    28        mC        250.0      38.9  160.1   77.6  90.0  90.0  75.9    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.9  160.1   77.6  90.0  90.0  75.9    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    36        oC        250.0      38.9  160.1   77.6  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.9  160.1   77.6  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    37        mC        250.0     160.1   38.9   77.6  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     160.1   38.9   77.6  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    27        mC        500.0     160.1   38.9  109.8  90.0 133.3  75.9    0  2  1  0  0  0  1  0  1 -1  0  0
    19        oI        500.0      38.9  109.8  116.5  90.0  70.5  90.0    0  0  1  0  1 -1  0  0  1  1  1  0
    26        oF        625.0      38.9  160.1  160.1  86.6 104.1 104.1    0  0 -1  0  0  2  1  0  2  0  1  0
    18        tI        625.0     109.8  116.5   38.9  70.5  90.0  90.0    1 -1  0  0 -1 -1 -1  0  0  0 -1  0
    30        mC        997.3      77.6  173.6   38.9  90.0  90.0  63.4    0  1  0  0 -2  1  0  0  0  0  1  0
    41        mC        997.3     173.6   77.6   38.9  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    40        oC        997.3      77.6  173.6   38.9  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        997.3      77.6   77.6   38.9  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.5  116.5  116.5  96.4  96.4 141.0    1 -1 -1  0 -1  1 -1  0  1  1  1  0
     2        hR        999.0      86.8   86.8  116.5 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.8   86.8  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.8   86.8  116.5 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.8   86.8   86.8  78.4  50.8  50.8    1  1  0  0  1  0 -1  0  0  1 -1  0
     7        tI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    1  0 -1  0  0  1 -1  0  1  1  0  0
     8        oI        999.0      86.8   86.8  109.8  50.8  50.8  78.4    0 -1  1  0 -1  0  1  0 -1 -1  0  0
     9        hR        999.0      38.9   86.8  248.5 102.1  99.0 116.6    0  0 -1  0  0  1  1  0  3 -1  1  0
    10        mC        999.0      86.8   86.8   77.6  90.0  90.0 126.8    0  1 -1  0  0 -1 -1  0 -1  0  0  0
    11        tP        999.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.9   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.9   77.6  177.9  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.8   86.8  177.9 107.0 119.2  53.2    0 -1  1  0  0 -1 -1  0  2  1 -1  0
    14        mC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.8   86.8   86.8  78.4 101.6  53.2    0 -1 -1  0  0 -1  1  0 -1  0  1  0
    24        hR        999.0     177.9  109.8   38.9  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.9   77.6  177.9 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.8  177.9   77.6 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    4.0      12    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.985   -23.509     0.022
 UNIT_CELL_B-AXIS=    23.509    73.985     0.022
 UNIT_CELL_C-AXIS=    -0.014    -0.007    38.867
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       101 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.27
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.13
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.627    77.627    38.850  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  5.3E-02 5.3E-02 3.6E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012882  0.012882  0.025740  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.984   -23.502     0.062
 COORDINATES OF UNIT CELL B-AXIS    23.502    73.984     0.017
 COORDINATES OF UNIT CELL C-AXIS    -0.032     0.001    38.850
 CRYSTAL MOSAICITY (DEGREES)     0.203
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006354  0.002224
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.90   1232.91
 DETECTOR ORIGIN (PIXELS) AT                     1260.66   1242.16
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.28
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.63    77.63    38.85  90.000  90.000  90.000


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
  0.05  0.15         7          8.348      11.230
  0.15  0.25       133         86.661       1.082
  0.25  0.35       156        103.505       0.906
  0.35  0.45       160         36.697       2.555
  0.45  0.55       172        138.711       0.676
  0.55  0.65       170         69.020       1.358
  0.65  0.75       172        139.711       0.671
  0.75  0.85       115         77.092       1.216
  0.85  0.95         9         52.349       1.791
 
 -0.05  0.95      1094         93.749       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        5
 DEGREES OF FREEDOM OF CHI^2 FIT               126.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.714
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    45.1 XMAX=    45.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59182 NYBIN=    5
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        272


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT               126.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.713
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    78.2 XMAX=  2450.5 NXBIN=    2
 YMIN=     6.5 YMAX=  2518.0 NYBIN=    2
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        272


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               120.5
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.723
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    45.1 XMAX=    45.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1264    1262
 DETECTOR_SURFACE_POSITION=    1667    1689
 DETECTOR_SURFACE_POSITION=     861    1689
 DETECTOR_SURFACE_POSITION=     861     836
 DETECTOR_SURFACE_POSITION=    1667     836
 DETECTOR_SURFACE_POSITION=    2176    1662
 DETECTOR_SURFACE_POSITION=    1642    2228
 DETECTOR_SURFACE_POSITION=     887    2228
 DETECTOR_SURFACE_POSITION=     353    1662
 DETECTOR_SURFACE_POSITION=     353     862
 DETECTOR_SURFACE_POSITION=     887     297
 DETECTOR_SURFACE_POSITION=    1642     297
 DETECTOR_SURFACE_POSITION=    2176     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        272


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
 7.818E-01  8.996E-04   37.71


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

   24.548   2.890    18.90   1.09      2.21      1.96      72     102      18
    2.890   2.051     2.13   1.02     30.64     29.69     160     213       2
    2.051   1.677    -0.04   0.83    -99.90    -99.90     208     263       0
    1.677   1.453     0.00   1.10   3670.30   3491.85     224     291       2
    1.453   1.300     0.12   1.01    423.59    423.37     140     200       2
   --------------------------------------------------------------------------
   24.548   1.300     2.24   1.00     10.54     10.30     804    1069      24
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      46      46)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.688  -0.1796E+01  0.1389E+02        -0.13           1 
   47    0    0     1.652  -0.1130E+01  0.1262E+02        -0.09           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -1.2% OF MEAN INTENSITY


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

     3.86          51      35      1285        2.7%       1.9%      1.7%       32   30.53      2.6%    99.9*     0    0.000       0
     2.74          69      44      2098        2.1%       5.2%      4.7%       50   13.35      7.4%    99.4*     0    0.000       0
     2.25         130      80      2665        3.0%      22.3%     23.2%      100    2.99     31.6%    89.9*     0    0.000       0
     1.95         120      73      3101        2.4%     225.0%    199.0%       94    0.45    318.2%     1.6      0    0.000       0
     1.74         137      83      3481        2.4%     -99.9%    -99.9%      108    0.00    -99.9%    -8.7      0    0.000       0
     1.59         169     103      3848        2.7%     -99.9%    -99.9%      132    0.00    -99.9%   -19.1      0    0.000       0
     1.47         160      98      4143        2.4%     -99.9%    -99.9%      124    0.00    -99.9%    17.4      0    0.000       0
     1.38         161     101      4429        2.3%     228.1%    229.6%      120    0.20    322.6%    17.3      0    0.000       0
     1.30          72      50      4743        1.1%     -99.9%    -99.9%       44    0.00    -99.9%   -35.5      0    0.000       0
    total        1069     667     29793        2.2%      10.5%     10.3%      804    2.90     14.9%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1094
 NUMBER OF REJECTED MISFITS                                24
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1069
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    667



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

   WILSON LINE (using all data) : A=   5.920 B=  26.596 CORRELATION=  0.79
      #      RES      SS        <I>       log(<I>)       BO
      14     7.396   0.005  1.8437E+03       7.520    -175.0
      17     4.756   0.011  1.4119E+03       7.253     -60.3
      18     3.595   0.019  4.8216E+02       6.178      -6.7
      24     2.984   0.028  9.7839E+01       4.583      23.8
      30     2.643   0.036  6.8756E+01       4.231      23.6
      43     2.409   0.043  2.3557E+01       3.159      32.0
      34     2.214   0.051  9.9109E+00       2.294      35.6
      35     2.062   0.059  3.4084E+00       1.226      39.9
      36     1.939   0.066  1.5655E-01      -1.854      58.5
      45     1.831   0.075 -1.5861E+00    -999.999   -1000.0
      43     1.743   0.082  4.7030E-01      -0.754      40.6
      42     1.672   0.089 -1.4908E+00    -999.999   -1000.0
      48     1.604   0.097  7.2785E-01      -0.318      32.1
      47     1.543   0.105 -6.5709E-01    -999.999   -1000.0
      56     1.486   0.113  1.0722E+00       0.070      25.8
      45     1.437   0.121  9.2176E-01      -0.081      24.8
      41     1.396   0.128  5.9547E+00       1.784      16.1
      34     1.356   0.136  1.6785E+00       0.518      19.9
      15     1.318   0.144 -4.3613E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.396        2.003        1.555        0.940
       3     4.756        1.263        0.764        0.344
       1     3.595        0.450        0.104        0.017
       3     2.984        2.557        1.511        0.660
       1     2.643        0.000        0.000        0.000
       3     2.409        0.234        0.046        0.007
       1     2.214        0.759        0.229        0.049
       6     2.062        1.340       -0.229        0.333
       6     1.939     2348.881    56167.897  1401471.051
       3     1.543       25.868       39.014       95.184
       1     1.396       34.483       70.145      101.920
       1     1.356       30.071       57.124       77.510
      71   overall      216.537     4697.944   118619.849


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.396        0.571        0.429        0.252
      14     4.756        1.299        1.656        1.675
      17     3.595        1.062        1.228        1.305
      21     2.984        0.631        0.502        0.349
      29     2.643        1.477        2.382        3.634
      40     2.409        1.124        1.063        0.883
      33     2.214        1.473        1.306        1.150
      29     2.062        5.105        6.681       15.286
      36     1.743      164.741      570.570    16873.286
      44     1.604       90.662      194.838     8110.055
      51     1.486       84.195      266.880     4722.219
      42     1.437       60.054      163.789     1838.459
      40     1.396        2.715        2.215        4.174
      33     1.356       23.183       -7.495      299.466
     596   overall      124.100     -302.065    76252.952

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.396  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.38  0.73
       3     4.756  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       1     3.595  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.984  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     2.643  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       3     2.409  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       1     2.214  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.062  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       6     1.939  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       3     1.543  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       1     1.396  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     1.356  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      71   overall  0.57  0.41  0.34  0.30  0.27  0.25  0.24  0.25  0.24  0.29


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.396  4.20  3.31  2.31  2.12  1.78  1.55  1.59  1.45  1.35  1.27
      14     4.756  3.00  1.58  1.65  1.30  1.09  1.11  1.42  1.43  1.44  1.36
      17     3.595  1.24  0.65  0.68  0.71  0.75  0.78  1.17  1.18  1.19  1.12
      21     2.984  1.00  1.05  1.47  1.30  1.09  1.37  1.42  1.30  1.28  1.21
      29     2.643  1.45  1.14  1.20  1.05  1.14  1.07  1.10  1.06  1.05  1.09
      40     2.409  1.58  1.24  1.06  1.14  1.02  1.00  0.94  0.95  1.01  0.95
      33     2.214  2.55  1.50  1.29  1.10  1.00  0.94  1.08  0.99  0.92  0.91
      29     2.062  3.26  1.71  1.33  1.15  1.05  0.92  0.82  0.75  0.70  0.76
      36     1.743  4.38  2.30  1.61  1.35  1.13  1.05  0.94  0.86  0.80  0.75
      44     1.604  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      51     1.486  4.74  2.49  1.74  1.37  1.15  1.00  0.90  0.89  0.83  0.78
      42     1.437  4.50  2.36  1.65  1.30  1.15  1.00  0.90  0.82  0.76  0.72
      40     1.396  3.15  2.07  1.54  1.44  1.21  1.11  0.99  0.91  0.84  0.87
      33     1.356  4.46  2.34  1.64  1.29  1.08  0.94  0.84  0.83  0.77  0.72
     596   overall  2.52  1.44  1.11  0.95  0.83  0.77  0.74  0.70  0.67  0.65


 List of     77 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   38    1    6    1.95  121.07  0.1895E+02  0.1523E+02 "alien"
   25   18   12    1.99   93.51  0.1464E+02  0.1264E+02 "alien"
   34   22    2    1.91   89.32  0.1398E+02  0.9193E+01 "alien"
   39    2    3    1.96   87.19  0.1365E+02  0.1685E+02 "alien"
   33    1   12    1.90   73.64  0.1153E+02  0.1256E+02 "alien"
   21   20   17    1.74   50.33  0.2367E+02  0.1049E+02 "alien"
   38   28    5    1.61   47.53  0.3459E+02  0.1743E+02 "alien"
   21   17   14    2.00   42.51  0.6656E+01  0.1256E+02 "alien"
   25   22   20    1.49   41.25  0.4423E+02  0.1365E+02 "alien"
   36   25    6    1.71   40.28  0.1894E+02  0.1544E+02 "alien"
   38    5   16    1.56   37.28 -0.2450E+02  0.1182E+02 "alien"
   30   20    8    1.97   36.88  0.5773E+01  0.8669E+01 "alien"
   32   27   14    1.54   36.51 -0.2399E+02  0.1319E+02 "alien"
   48    2    7    1.55   31.51 -0.2071E+02  0.9918E+01 "alien"
   33    7   19    1.53   28.84 -0.1895E+02  0.8263E+01 "alien"
   39   29    5    1.56   28.01 -0.1840E+02  0.9228E+01 "alien"
   29   22   12    1.78   27.70  0.1303E+02  0.8211E+01 "alien"
   46    5   14    1.44   26.79  0.2469E+02  0.1585E+02 "alien"
   33    8   20    1.48   26.35  0.2825E+02  0.1340E+02 "alien"
   47    5   13    1.44   25.94  0.2391E+02  0.1075E+02 "alien"
   34    0   10    1.97   25.51  0.2661E+02  0.1446E+02 "alien"
   30    1   13    1.96   24.78  0.3879E+01  0.1311E+02 "alien"
   33   29   14    1.49   24.29  0.2604E+02  0.1692E+02 "alien"
   20   19   18    1.71   23.67  0.1113E+02  0.1070E+02 "alien"
   25   19   13    1.90   23.29  0.3647E+01  0.8580E+01 "alien"
   34   29   13    1.50   23.25  0.2493E+02  0.1566E+02 "alien"
   40    1    5    1.88   22.44 -0.3560E+02  0.1251E+02 "alien"
   39    6   17    1.49   22.37  0.2399E+02  0.9066E+01 "alien"
   45    3   11    1.55   21.76 -0.1430E+02  0.9808E+01 "alien"
   49    5   12    1.42   21.39  0.1972E+02  0.1080E+02 "alien"
   46    1    6    1.63   20.48  0.1490E+02  0.1380E+02 "alien"
   40   30    5    1.52   20.43 -0.1343E+02  0.9494E+01 "alien"
   22   16   21    1.55   19.71 -0.1295E+02  0.1047E+02 "alien"
   34   22    0    1.92   18.83  0.1119E+02  0.1516E+02 "alien"
   42   32    3    1.46   18.60  0.1715E+02  0.1188E+02 "alien"
   34   22    3    1.90   18.49  0.2895E+01  0.8526E+01 "alien"
   32   26   13    1.59   17.86  0.1300E+02  0.9231E+01 "alien"
   29   20   10    1.92   17.25  0.2700E+01  0.1402E+02 "alien"
   43    1    3    1.79   17.23  0.8104E+01  0.9737E+01 "alien"
   41   32    6    1.45   16.35  0.1507E+02  0.1019E+02 "alien"
   25   21   15    1.75   16.15  0.7597E+01  0.1096E+02 "alien"
   36   28   10    1.56   15.57 -0.1023E+02  0.9011E+01 "alien"
   36   25    5    1.73   15.56  0.7316E+01  0.8346E+01 "alien"
   39    2    0    1.99   15.52  0.7153E+01  0.1233E+02 "alien"
   36   31   13    1.43   15.24  0.1405E+02  0.1003E+02 "alien"
   41   31    2    1.51   14.88  0.1596E+02  0.9693E+01 "alien"
   44    2   10    1.60   14.78  0.1075E+02  0.8782E+01 "alien"
   30   25   14    1.62   14.47  0.1053E+02  0.8287E+01 "alien"
   28   19   23    1.36   14.14  0.2374E+02  0.1384E+02 "alien"
   53    2    1    1.46   13.98  0.1499E+02  0.1740E+02 "alien"
   34    3   15    1.71   13.92  0.6549E+01  0.1331E+02 "alien"
   24   19   14    1.87   13.65 -0.2164E+02  0.1180E+02 "alien"
   38   29    7    1.56   13.50 -0.8869E+01  0.1500E+02 "alien"
   45    4   13    1.49   13.25  0.1420E+02  0.1097E+02 "alien"
   33   26   12    1.60   12.98  0.9447E+01  0.1006E+02 "alien"
   29   23   13    1.72   12.58  0.5916E+01  0.8190E+01 "alien"
   50    7   14    1.34   12.56  0.2109E+02  0.1610E+02 "alien"
   20   17   15    1.95   12.46  0.1950E+01  0.7433E+01 "alien"
   36   27    9    1.60   12.46  0.9066E+01  0.9931E+01 "alien"
   36    2   13    1.75   12.20  0.5739E+01  0.7944E+01 "alien"
   39    5   16    1.53   12.11 -0.7954E+01  0.8904E+01 "alien"
   52    2    2    1.49   11.98  0.1284E+02  0.9169E+01 "alien"
   37   26    3    1.70   11.94 -0.1780E+02  0.1558E+02 "alien"
   37    4   15    1.62   11.90  0.8659E+01  0.1363E+02 "alien"
   53    3    6    1.43   11.59  0.1068E+02  0.1682E+02 "alien"
   45    2   10    1.58   11.37  0.8275E+01  0.1076E+02 "alien"
   32   23    9    1.79   10.79 -0.1711E+02  0.1424E+02 "alien"
   37   25    1    1.74   10.59  0.4979E+01  0.8975E+01 "alien"
   35   25    8    1.69   10.49 -0.1564E+02  0.1392E+02 "alien"
   51    2    5    1.49   10.40  0.1115E+02  0.1024E+02 "alien"
   42   32    2    1.47   10.39  0.1114E+02  0.1005E+02 "alien"
   36   25    4    1.74   10.17  0.4782E+01  0.1055E+02 "alien"
   54    3    3    1.43   10.17  0.9371E+01  0.9578E+01 "alien"
   38   31   10    1.47   10.14  0.1087E+02  0.9748E+01 "alien"
   36   24    3    1.78   10.10  0.4749E+01  0.1010E+02 "alien"
   34   26   10    1.64   10.03 -0.1495E+02  0.1383E+02 "alien"
   24   23   17    1.63   10.02  0.7291E+01  0.1245E+02 "alien"


 List of     77 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   50    7   14    1.34   12.56  0.2109E+02  0.1610E+02
   28   19   23    1.36   14.14  0.2374E+02  0.1384E+02
   49    5   12    1.42   21.39  0.1972E+02  0.1080E+02
   53    3    6    1.43   11.59  0.1068E+02  0.1682E+02
   54    3    3    1.43   10.17  0.9371E+01  0.9578E+01
   36   31   13    1.43   15.24  0.1405E+02  0.1003E+02
   46    5   14    1.44   26.79  0.2469E+02  0.1585E+02
   47    5   13    1.44   25.94  0.2391E+02  0.1075E+02
   41   32    6    1.45   16.35  0.1507E+02  0.1019E+02
   42   32    3    1.46   18.60  0.1715E+02  0.1188E+02
   53    2    1    1.46   13.98  0.1499E+02  0.1740E+02
   38   31   10    1.47   10.14  0.1087E+02  0.9748E+01
   42   32    2    1.47   10.39  0.1114E+02  0.1005E+02
   33    8   20    1.48   26.35  0.2825E+02  0.1340E+02
   52    2    2    1.49   11.98  0.1284E+02  0.9169E+01
   45    4   13    1.49   13.25  0.1420E+02  0.1097E+02
   33   29   14    1.49   24.29  0.2604E+02  0.1692E+02
   39    6   17    1.49   22.37  0.2399E+02  0.9066E+01
   25   22   20    1.49   41.25  0.4423E+02  0.1365E+02
   51    2    5    1.49   10.40  0.1115E+02  0.1024E+02
   34   29   13    1.50   23.25  0.2493E+02  0.1566E+02
   41   31    2    1.51   14.88  0.1596E+02  0.9693E+01
   40   30    5    1.52   20.43 -0.1343E+02  0.9494E+01
   33    7   19    1.53   28.84 -0.1895E+02  0.8263E+01
   39    5   16    1.53   12.11 -0.7954E+01  0.8904E+01
   32   27   14    1.54   36.51 -0.2399E+02  0.1319E+02
   45    3   11    1.55   21.76 -0.1430E+02  0.9808E+01
   48    2    7    1.55   31.51 -0.2071E+02  0.9918E+01
   22   16   21    1.55   19.71 -0.1295E+02  0.1047E+02
   38    5   16    1.56   37.28 -0.2450E+02  0.1182E+02
   38   29    7    1.56   13.50 -0.8869E+01  0.1500E+02
   36   28   10    1.56   15.57 -0.1023E+02  0.9011E+01
   39   29    5    1.56   28.01 -0.1840E+02  0.9228E+01
   45    2   10    1.58   11.37  0.8275E+01  0.1076E+02
   32   26   13    1.59   17.86  0.1300E+02  0.9231E+01
   36   27    9    1.60   12.46  0.9066E+01  0.9931E+01
   33   26   12    1.60   12.98  0.9447E+01  0.1006E+02
   44    2   10    1.60   14.78  0.1075E+02  0.8782E+01
   38   28    5    1.61   47.53  0.3459E+02  0.1743E+02
   30   25   14    1.62   14.47  0.1053E+02  0.8287E+01
   37    4   15    1.62   11.90  0.8659E+01  0.1363E+02
   46    1    6    1.63   20.48  0.1490E+02  0.1380E+02
   24   23   17    1.63   10.02  0.7291E+01  0.1245E+02
   34   26   10    1.64   10.03 -0.1495E+02  0.1383E+02
   35   25    8    1.69   10.49 -0.1564E+02  0.1392E+02
   37   26    3    1.70   11.94 -0.1780E+02  0.1558E+02
   36   25    6    1.71   40.28  0.1894E+02  0.1544E+02
   34    3   15    1.71   13.92  0.6549E+01  0.1331E+02
   20   19   18    1.71   23.67  0.1113E+02  0.1070E+02
   29   23   13    1.72   12.58  0.5916E+01  0.8190E+01
   36   25    5    1.73   15.56  0.7316E+01  0.8346E+01
   37   25    1    1.74   10.59  0.4979E+01  0.8975E+01
   21   20   17    1.74   50.33  0.2367E+02  0.1049E+02
   36   25    4    1.74   10.17  0.4782E+01  0.1055E+02
   36    2   13    1.75   12.20  0.5739E+01  0.7944E+01
   25   21   15    1.75   16.15  0.7597E+01  0.1096E+02
   36   24    3    1.78   10.10  0.4749E+01  0.1010E+02
   29   22   12    1.78   27.70  0.1303E+02  0.8211E+01
   43    1    3    1.79   17.23  0.8104E+01  0.9737E+01
   32   23    9    1.79   10.79 -0.1711E+02  0.1424E+02
   24   19   14    1.87   13.65 -0.2164E+02  0.1180E+02
   40    1    5    1.88   22.44 -0.3560E+02  0.1251E+02
   34   22    3    1.90   18.49  0.2895E+01  0.8526E+01
   33    1   12    1.90   73.64  0.1153E+02  0.1256E+02
   25   19   13    1.90   23.29  0.3647E+01  0.8580E+01
   34   22    2    1.91   89.32  0.1398E+02  0.9193E+01
   29   20   10    1.92   17.25  0.2700E+01  0.1402E+02
   34   22    0    1.92   18.83  0.1119E+02  0.1516E+02
   38    1    6    1.95  121.07  0.1895E+02  0.1523E+02
   20   17   15    1.95   12.46  0.1950E+01  0.7433E+01
   30    1   13    1.96   24.78  0.3879E+01  0.1311E+02
   39    2    3    1.96   87.19  0.1365E+02  0.1685E+02
   30   20    8    1.97   36.88  0.5773E+01  0.8669E+01
   34    0   10    1.97   25.51  0.2661E+02  0.1446E+02
   39    2    0    1.99   15.52  0.7153E+01  0.1233E+02
   25   18   12    1.99   93.51  0.1464E+02  0.1264E+02
   21   17   14    2.00   42.51  0.6656E+01  0.1256E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        32
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1094
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1026
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       68


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
