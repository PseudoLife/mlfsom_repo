
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)   6-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.251    77.251    38.666  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/fake/myfake_???.img SMV
 DATA_RANGE=       1       1
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       1
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  1.000000 -0.000469  0.000730
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
 BEAM_DIVERGENCE_E.S.D.=   0.043
 REFLECTING_RANGE_E.S.D.=   0.090
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

    1495 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       1 ...       1
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1495 REFLECTIONS ACCEPTED



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
 UNIT_CELL_CONSTANTS=    77.00    77.00    38.00  90.000  90.000  90.000
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

 *  44        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.2  109.2   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  23        oC          0.0     109.2  109.2   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  20        mC          0.0     109.2  109.2   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    37        mC        250.0     159.3   38.7   77.3  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.3   38.7   77.3  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.7  159.3   77.3  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.7  159.3   77.3  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.7  159.3   77.3  90.0  90.0  75.9    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.7  159.3   77.3  90.0  90.0  75.9    0  0  1  0  0  2  1  0 -1  0  0  0
    19        oI        500.0      38.7  109.2  115.9  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    27        mC        500.0     159.3   38.7  109.2  90.0 133.3  75.9    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    18        tI        625.0     109.2  115.9   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    26        oF        625.0      38.7  159.3  159.3  86.6 104.1 104.1    0  0  1  0  0  2 -1  0 -2  0 -1  0
    41        mC        997.9     172.7   77.3   38.7  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    22        hP        997.9      77.3   77.3   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    30        mC        997.9      77.3  172.7   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        997.9      77.3  172.7   38.7  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
     1        cF        999.0     115.9  115.9  115.9  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.4   86.4  115.9 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.4   86.4  109.2  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.4   86.4  115.9 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.2   86.4   86.4  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.4   86.4  109.2  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.4   86.4  109.2  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.4  247.3 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.4   86.4   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.7   77.3  177.0  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.4   86.4  177.0 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.4   86.4   86.4  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.0  109.2   38.7  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.7   77.3  177.0 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.4  177.0   77.3 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    6.4      19    0.02  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.00    77.00    38.00  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.617   -23.413     0.030
 UNIT_CELL_B-AXIS=    23.413    73.617     0.010
 UNIT_CELL_C-AXIS=    -0.016    -0.000    38.667
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       568 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.37
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.224    77.224    38.654  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.1E-02 2.1E-02 7.1E-03 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012949  0.012949  0.025871  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.597   -23.388     0.018
 COORDINATES OF UNIT CELL B-AXIS    23.388    73.597     0.006
 COORDINATES OF UNIT CELL C-AXIS    -0.009     0.000    38.654
 CRYSTAL MOSAICITY (DEGREES)     0.090
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006383  0.002020
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.85   1232.92
 DETECTOR ORIGIN (PIXELS) AT                     1260.92   1242.21
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.14
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.22    77.22    38.65  90.000  90.000  90.000


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
  0.05  0.15       118        644.359       0.975
  0.15  0.25       194        532.596       1.180
  0.25  0.35       170        614.318       1.023
  0.35  0.45       157        428.300       1.467
  0.45  0.55       179        829.818       0.757
  0.55  0.65       193        553.040       1.136
  0.65  0.75       186        535.798       1.173
  0.75  0.85       167        613.970       1.023
  0.85  0.95       131        998.717       0.629
 
 -0.05  0.95      1495        628.316       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       20
 DEGREES OF FREEDOM OF CHI^2 FIT               505.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.893
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     0.1 XMAX=     0.9 NXBIN=    1
 YMIN= 0.00168 YMAX= 0.59795 NYBIN=   20
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS       1062


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       25
 DEGREES OF FREEDOM OF CHI^2 FIT               504.9
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.887
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    10.9 XMAX=  2439.9 NXBIN=    5
 YMIN=     5.9 YMAX=  2520.9 NYBIN=    5
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS       1062


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               505.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.885
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     0.1 XMAX=     0.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1225    1263
 DETECTOR_SURFACE_POSITION=    1638    1691
 DETECTOR_SURFACE_POSITION=     813    1691
 DETECTOR_SURFACE_POSITION=     813     836
 DETECTOR_SURFACE_POSITION=    1638     836
 DETECTOR_SURFACE_POSITION=    2159    1664
 DETECTOR_SURFACE_POSITION=    1612    2230
 DETECTOR_SURFACE_POSITION=     839    2230
 DETECTOR_SURFACE_POSITION=     292    1664
 DETECTOR_SURFACE_POSITION=     292     863
 DETECTOR_SURFACE_POSITION=     839     297
 DETECTOR_SURFACE_POSITION=    1612     297
 DETECTOR_SURFACE_POSITION=    2159     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS       1062


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
 1.615E+00  1.918E-03   17.96


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

   24.420   3.393    16.65   1.59      2.73      3.19      62      88      22
    3.393   2.411    14.26   1.47      4.03      3.44     128     169      12
    2.411   1.972    11.61   1.05      4.06      3.98     176     215      12
    1.972   1.709     7.01   0.89      5.93      6.63     192     241       2
    1.709   1.529     4.56   0.71      8.38      9.85     214     281       0
    1.529   1.396     2.39   0.58     17.64     23.33     240     306       0
    1.396   1.293     1.38   0.44     25.26     38.06      98     142       4
   --------------------------------------------------------------------------
   24.420   1.293     7.12   0.88      4.16      4.42    1110    1442      52
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       1       1)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.679   0.6124E+02  0.1790E+02         3.42           1 
   47    0    0     1.643   0.2531E+02  0.1660E+02         1.52           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     3.9% OF MEAN INTENSITY


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

     3.84          66      43      1285        3.3%       2.6%      3.2%       46   20.57      3.7%    99.9*     0    0.000       0
     2.73         106      68      2098        3.2%       3.7%      3.3%       76   18.59      5.2%    99.1*     0    0.000       0
     2.23         161      97      2664        3.6%       3.8%      3.6%      128   17.31      5.3%    99.2*     0    0.000       0
     1.94         163      96      3102        3.1%       4.7%      4.3%      134   13.57      6.6%    99.4*     0    0.000       0
     1.73         186     111      3481        3.2%       6.8%      7.1%      150    8.71      9.6%    98.9*     0    0.000       0
     1.58         216     132      3848        3.4%       7.0%      8.7%      168    6.35     10.0%    98.8*     0    0.000       0
     1.47         244     150      4142        3.6%      12.4%     15.9%      188    4.08     17.5%    95.5*     0    0.000       0
     1.37         215     130      4433        2.9%      22.1%     28.6%      170    2.54     31.3%    81.7*     0    0.000       0
     1.29          85      60      4737        1.3%      25.9%     50.0%       50    1.43     36.6%    86.9*     0    0.000       0
    total        1442     887     29790        3.0%       4.2%      4.4%     1110    8.98      5.9%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1495
 NUMBER OF REJECTED MISFITS                                52
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1442
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    887



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

   WILSON LINE (using all data) : A=   8.680 B=  20.370 CORRELATION=  0.99
      #      RES      SS        <I>       log(<I>)       BO
      18     7.513   0.004  4.8260E+03       8.482      22.4
      21     4.667   0.011  5.6968E+03       8.648       1.4
      30     3.491   0.021  2.6665E+03       7.889      19.3
      33     2.995   0.028  1.4770E+03       7.298      24.8
      39     2.640   0.036  1.5370E+03       7.338      18.7
      52     2.392   0.044  1.0479E+03       6.955      19.8
      41     2.204   0.051  8.4070E+02       6.734      18.9
      47     2.054   0.059  4.6153E+02       6.135      21.5
      49     1.929   0.067  3.3659E+02       5.819      21.3
      60     1.826   0.075  2.1960E+02       5.392      21.9
      54     1.736   0.083  1.4264E+02       4.960      22.4
      58     1.662   0.090  1.3626E+02       4.915      20.8
      61     1.593   0.098  9.8171E+01       4.587      20.8
      78     1.533   0.106  7.3307E+01       4.295      20.6
      77     1.478   0.114  5.8826E+01       4.075      20.1
      58     1.430   0.122  3.5678E+01       3.575      20.9
      53     1.389   0.130  3.7809E+01       3.633      19.5
      38     1.349   0.137  2.7007E+01       3.296      19.6
      20     1.312   0.145  1.8079E+01       2.895      19.9


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.513        0.966        0.427        0.145
       3     4.667        0.998        0.535        0.214
       3     3.491        0.400        0.109        0.021
       4     2.995        0.984        0.421        0.138
       3     2.640        0.282        0.075        0.015
       5     2.392        0.898        0.452        0.171
       3     2.204        1.179        0.456        0.128
       8     2.054        0.224        0.060        0.012
       7     1.929        1.851        1.562        0.964
      11     1.826        0.642        0.344        0.147
      11     1.736        1.070        0.748        0.397
      10     1.662        0.584        0.394        0.209
       3     1.593        0.209        0.052        0.010
       3     1.533        0.024        0.002        0.000
       6     1.478        0.835        0.345        0.110
       3     1.430        0.767        0.380        0.139
       2     1.349        0.033        0.003        0.000
      91   overall        0.762        0.443        0.212


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      13     7.513        1.192        1.695        1.996
      18     4.667        1.006        1.014        0.887
      27     3.491        0.832        0.603        0.392
      29     2.995        0.865        0.739        0.568
      36     2.640        1.106        1.182        1.177
      47     2.392        0.895        0.710        0.477
      38     2.204        0.779        0.573        0.364
      39     2.054        1.182        1.251        1.220
      42     1.929        0.803        0.718        0.623
      49     1.826        1.090        1.058        0.885
      43     1.736        0.873        0.701        0.500
      48     1.662        1.279        1.505        1.604
      58     1.593        1.018        1.048        1.081
      75     1.533        1.057        1.078        1.070
      71     1.478        0.921        0.845        0.717
      55     1.430        0.950        0.818        0.644
      52     1.389        1.195        1.385        1.604
      36     1.349        1.255        1.289        1.154
      20     1.312        2.174        2.950        4.241
     796   overall        1.045        1.046        1.008

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.513  0.00  0.00  0.00  0.00  0.00  0.36  0.33  0.32  0.61  0.59
       3     4.667  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.491  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.995  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.640  1.34  0.96  0.80  0.70  0.64  0.59  1.12  1.06  1.01  0.98
       5     2.392  0.00  0.00  0.00  0.42  0.77  1.07  1.00  0.95  0.91  0.88
       3     2.204  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.054  1.01  1.08  0.90  1.06  0.96  0.89  0.84  0.99  1.14  1.28
       7     1.929  0.00  0.00  0.69  0.91  1.10  1.02  0.96  0.91  1.09  1.05
      11     1.826  1.10  1.05  0.87  0.77  0.87  1.13  1.07  1.01  0.97  0.93
      11     1.736  0.73  1.05  1.09  0.96  0.87  0.97  0.91  1.01  1.11  1.07
      10     1.662  1.21  1.16  1.20  1.06  1.15  1.07  1.00  0.95  1.22  1.32
       3     1.593  0.00  0.96  0.80  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.533  1.34  0.96  1.60  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       6     1.478  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.53  0.51  0.73
       3     1.430  0.00  0.00  0.00  0.70  1.28  1.19  1.12  1.06  1.01  0.98
       2     1.349  2.01  1.45  1.20  1.06  1.92  1.78  1.67  1.59  1.52  1.47
      91   overall  0.66  0.67  0.69  0.74  0.80  0.86  0.83  0.86  0.92  0.95


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      13     7.513  1.62  2.12  1.78  2.33  1.95  1.70  1.53  1.40  1.30  1.22
      18     4.667  1.17  1.53  1.07  1.18  1.13  1.11  0.99  1.01  1.12  1.05
      27     3.491  1.17  0.61  0.57  0.79  0.66  0.90  0.88  0.94  0.94  0.88
      29     2.995  1.45  1.52  1.06  1.05  1.14  1.15  1.16  1.13  1.05  1.04
      36     2.640  0.88  0.92  1.29  1.18  0.99  0.92  0.99  1.06  1.03  1.01
      47     2.392  1.34  1.06  1.07  1.10  1.03  1.04  1.06  1.00  1.04  1.04
      38     2.204  0.83  1.31  1.02  0.80  1.00  0.93  0.94  1.00  1.06  1.08
      39     2.054  1.08  1.13  1.09  0.93  0.91  0.80  0.87  0.93  0.99  0.97
      42     1.929  0.50  1.18  1.10  1.16  1.15  1.00  1.04  0.95  1.00  0.98
      49     1.826  1.07  0.90  1.18  1.18  1.30  1.13  1.01  0.93  0.96  1.03
      43     1.736  1.22  0.77  0.81  0.92  1.00  1.13  1.11  1.01  0.98  1.03
      48     1.662  1.09  1.26  1.45  1.26  1.11  1.06  1.03  1.10  1.09  1.02
      58     1.593  0.72  0.95  0.86  0.99  1.01  0.92  0.92  0.94  0.93  0.93
      75     1.533  1.12  0.88  0.82  0.81  0.88  0.92  0.95  0.97  0.97  1.01
      71     1.478  0.59  0.70  0.98  1.03  1.11  1.06  1.09  1.05  1.02  0.98
      55     1.430  1.53  1.20  0.98  1.16  1.02  0.97  0.90  0.86  0.89  0.83
      52     1.389  2.02  1.27  0.96  0.93  0.98  1.07  1.07  1.01  1.07  1.00
      36     1.349  2.04  1.38  1.29  1.26  1.20  1.05  1.05  1.01  0.94  0.92
      20     1.312  3.15  2.21  1.54  1.21  1.02  0.89  0.89  1.00  0.93  0.95
     796   overall  1.20  1.10  1.05  1.06  1.05  1.01  1.01  1.00  1.00  0.99

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT         0
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1495
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1443
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       52


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
