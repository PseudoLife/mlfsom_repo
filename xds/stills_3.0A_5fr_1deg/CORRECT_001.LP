
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.301    77.301    38.674  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_3.0A_5fr_1deg/stills_3.0A_5fr_1deg_???_001.img SMV
 DATA_RANGE=       1       1
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       1
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000952  0.000646
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
 BEAM_DIVERGENCE_E.S.D.=   0.036
 REFLECTING_RANGE_E.S.D.=   0.078
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

    1528 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       1 ...       1
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1528 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               46
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
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    27        mC        500.0     159.4   38.7  109.3  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.7  109.3  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.3  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    30        mC        998.8      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.8      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    41        mC        998.8     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    22        hP        998.8      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.5    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.4   86.4  109.3  50.8  50.8  78.5   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.5    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.3   86.4   86.4  78.5  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.4   86.4  109.3  50.8  50.8  78.5   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.4   86.4  109.3  50.8  50.8  78.5    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.7   86.4  247.5 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.4   86.4   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.1  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.4   86.4  177.1 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.4   86.4   86.4  78.5 101.5  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.1  109.3   38.7  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.7   77.3  177.1 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.4  177.1   77.3 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    4.9      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.657   -23.457     0.045
 UNIT_CELL_B-AXIS=    23.457    73.657     0.014
 UNIT_CELL_C-AXIS=    -0.024     0.000    38.675
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       148 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.10
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.217    77.217    38.644  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  3.2E-02 3.2E-02 1.4E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012951  0.012951  0.025877  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.593   -23.377     0.027
 COORDINATES OF UNIT CELL B-AXIS    23.377    73.593     0.010
 COORDINATES OF UNIT CELL C-AXIS    -0.014    -0.001    38.644
 CRYSTAL MOSAICITY (DEGREES)     0.078
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006437  0.001783
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.88   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1261.29   1242.25
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.08
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.22    77.22    38.64  90.000  90.000  90.000


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
  0.05  0.15       128        171.451       1.151
  0.15  0.25       196        161.144       1.225
  0.25  0.35       178        183.296       1.077
  0.35  0.45       146         70.888       2.784
  0.45  0.55       184        342.268       0.577
  0.55  0.65       190        119.794       1.647
  0.65  0.75       190        233.024       0.847
  0.75  0.85       169        120.695       1.635
  0.85  0.95       147        371.739       0.531
 
 -0.05  0.95      1528        197.357       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               144.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.582
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     0.1 XMAX=     0.9 NXBIN=    1
 YMIN= 0.00168 YMAX= 0.59814 NYBIN=    6
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        312


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               144.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.581
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    10.4 XMAX=  2439.8 NXBIN=    2
 YMIN=     6.6 YMAX=  2520.6 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        312


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               144.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.577
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     0.1 XMAX=     0.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1225    1264
 DETECTOR_SURFACE_POSITION=    1638    1691
 DETECTOR_SURFACE_POSITION=     812    1691
 DETECTOR_SURFACE_POSITION=     812     837
 DETECTOR_SURFACE_POSITION=    1638     837
 DETECTOR_SURFACE_POSITION=    2159    1664
 DETECTOR_SURFACE_POSITION=    1612    2230
 DETECTOR_SURFACE_POSITION=     838    2230
 DETECTOR_SURFACE_POSITION=     291    1664
 DETECTOR_SURFACE_POSITION=     291     863
 DETECTOR_SURFACE_POSITION=     838     297
 DETECTOR_SURFACE_POSITION=    1612     297
 DETECTOR_SURFACE_POSITION=    2159     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        312


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
 8.514E-01  4.426E-03   16.29


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

   24.418   3.393    15.41   1.00      2.73      3.51      68      89      24
    3.393   2.411     3.27   0.95      6.66      5.53     150     184       0
    2.411   1.971     0.03   0.80    714.67    807.40     194     228       2
    1.971   1.709    -0.05   0.99    -99.90    -99.90     206     250       0
    1.709   1.529     0.02   0.82   2205.38   2420.57     232     289       0
    1.529   1.396     0.02   0.79    -99.90    -99.90     258     315       0
    1.396   1.293     0.04   1.05   1633.51   1541.41     102     146       0
   --------------------------------------------------------------------------
   24.418   1.293     1.32   0.88      5.96      6.61    1210    1501      26
 


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

   46    0    0     1.679   0.5924E+01  0.1001E+02         0.59           1 
   47    0    0     1.643  -0.2793E+01  0.9928E+01        -0.28           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -1.3% OF MEAN INTENSITY


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

     3.84          67      42      1285        3.3%       2.4%      3.5%       50   19.22      3.4%    99.8*     0    0.000       0
     2.73         114      70      2098        3.3%       4.6%      3.8%       88   10.41      6.5%    99.3*     0    0.000       0
     2.23         173      99      2665        3.7%    2774.5%   2974.1%      148    0.04   3923.8%    13.5      0    0.000       0
     1.94         173      99      3101        3.2%     760.6%    877.4%      148    0.07   1075.6%     5.7      0    0.000       0
     1.73         192     113      3479        3.2%     -99.9%    -99.9%      158    0.00    -99.9%    -3.5      0    0.000       0
     1.58         223     132      3850        3.4%     679.9%    726.9%      182    0.09    961.5%    -6.2      0    0.000       0
     1.47         250     150      4143        3.6%     -99.9%    -99.9%      200    0.00    -99.9%     1.0      0    0.000       0
     1.37         219     128      4428        2.9%    1450.4%   1545.1%      182    0.09   2051.2%     7.6      0    0.000       0
     1.29          88      62      4741        1.3%     750.6%    893.9%       52    0.09   1061.5%   -10.4      0    0.000       0
    total        1499     895     29790        3.0%       6.0%      6.6%     1208    1.74      8.4%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1528
 NUMBER OF REJECTED MISFITS                                26
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1501
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    896



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

   WILSON LINE (using all data) : A=   5.148 B=  25.231 CORRELATION=  0.63
      #      RES      SS        <I>       log(<I>)       BO
      18     7.512   0.004  3.2295E+03       8.080    -330.9
      20     4.639   0.012  3.4889E+03       8.157    -129.5
      29     3.493   0.020  1.7174E+03       7.449     -56.1
      36     2.985   0.028  3.8366E+02       5.950     -14.3
      41     2.635   0.036  1.7379E-02      -4.053     127.8
      52     2.389   0.044  2.6733E-01      -1.319      73.8
      43     2.201   0.052  1.5664E+00       0.449      45.5
      48     2.054   0.059 -5.0615E-01    -999.999   -1000.0
      49     1.928   0.067 -3.2553E-01    -999.999   -1000.0
      62     1.825   0.075  5.7046E-01      -0.561      38.0
      54     1.736   0.083 -1.2573E+00    -999.999   -1000.0
      58     1.662   0.091  1.2693E+00       0.238      27.1
      61     1.593   0.098  2.0833E-01      -1.569      34.1
      78     1.532   0.107 -1.1040E+00    -999.999   -1000.0
      78     1.478   0.115  1.1233E+00       0.116      22.0
      56     1.429   0.122  8.1489E-01      -0.205      21.9
      53     1.389   0.130  1.4497E+00       0.371      18.4
      39     1.350   0.137  1.6246E-01      -1.817      25.4
      19     1.313   0.145  9.7398E-01      -0.026      17.8


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.512        0.948        0.444        0.164
       3     4.639        0.987        0.528        0.211
       3     3.493        0.374        0.098        0.018
       4     2.985        1.393        1.139        0.666
       3     2.635    74112.660  9817329.473*************
       3     2.201        7.934        9.473        8.143
       7     1.928       31.775       61.470      230.142
      11     1.736       31.116       21.218      141.755
      11     1.662        9.028       12.583       15.297
       3     1.593      775.515     4799.331    70298.112
       6     1.478       27.935      -13.401      147.286
       3     1.429       41.345      135.748      355.566
      91   overall     2521.942   323310.292 33147671.220


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      13     7.512        1.430        2.393        3.266
      17     4.639        1.057        1.116        1.009
      26     3.493        0.813        0.580        0.374
      32     2.985        1.943        2.979        3.685
      47     2.389      370.582       18.790    60568.491
      40     2.201        7.241       -0.395       35.348
      51     1.825       78.918      -99.533     4772.852
      47     1.662       15.005        0.521      125.843
      58     1.593      502.129     -422.368   103650.168
      72     1.478       23.310       45.350      325.087
      53     1.429       46.632       74.156     1206.485
      53     1.389       13.042       10.335       64.803
      37     1.350      823.107    11096.691   405571.704
      19     1.313      128.143     1109.851    17051.692
     803   overall     3686.491   -56579.369160776057.204

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.512  0.00  0.00  0.00  0.00  0.00  0.36  0.33  0.64  0.61  0.59
       3     4.639  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       3     3.493  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.985  3.02  2.17  1.80  1.58  1.44  1.34  1.26  1.19  1.14  1.10
       3     2.635  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.201  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       7     1.928  2.30  1.65  1.37  1.21  1.10  1.02  0.96  0.91  0.87  0.84
      11     1.736  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.87  0.83  0.80
      11     1.662  1.10  0.79  0.65  0.58  0.52  0.49  0.46  0.43  0.42  0.40
       3     1.593  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       6     1.478  2.01  1.45  1.20  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.429  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      91   overall  1.02  0.73  0.61  0.56  0.51  0.49  0.46  0.47  0.47  0.45


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      13     7.512  1.62  2.12  2.08  2.33  1.95  1.70  1.53  1.40  1.30  1.22
      17     4.639  1.24  1.62  1.13  1.25  1.20  1.04  1.05  1.07  1.09  1.12
      26     3.493  0.81  0.42  0.45  0.58  0.68  0.85  0.92  0.91  0.91  0.97
      32     2.985  6.57  3.45  2.41  1.90  1.67  1.45  1.30  1.25  1.21  1.14
      47     2.389  4.92  2.58  1.81  1.42  1.19  1.04  0.93  0.85  0.79  0.74
      40     2.201  4.47  2.34  1.64  1.36  1.33  1.16  1.09  1.00  0.93  0.91
      51     1.825  5.15  2.70  1.89  1.49  1.25  1.09  1.01  0.96  0.89  0.84
      47     1.662  4.47  2.35  1.64  1.29  1.08  0.94  0.85  0.85  0.82  0.84
      58     1.593  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      72     1.478  4.82  2.61  1.82  1.47  1.24  1.11  0.99  0.93  0.87  0.83
      53     1.429  5.16  2.81  1.97  1.55  1.29  1.13  1.01  0.96  0.89  0.84
      53     1.389  4.36  2.50  1.75  1.43  1.25  1.13  1.05  0.96  0.92  0.87
      37     1.350  4.83  2.53  1.77  1.39  1.17  1.02  0.91  0.83  0.77  0.73
      19     1.313  6.08  3.19  2.23  1.76  1.47  1.28  1.15  1.05  0.98  0.92
     803   overall  3.25  1.77  1.25  1.02  0.88  0.78  0.72  0.67  0.64  0.61


 List of    146 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   27    2    7    2.53 1064.48  0.1850E+02  0.9586E+01 "alien"
   29    3    3    2.59  521.71  0.9067E+01  0.6318E+01 "alien"
   26   14    1    2.61  517.90  0.9000E+01  0.6363E+01 "alien"
   20    1   10    2.73  485.14  0.8431E+01  0.6363E+01 "alien"
   21    1   10    2.66  393.18  0.6833E+01  0.4776E+01 "alien"
   24   13    3    2.76  375.37  0.6523E+01  0.6065E+01 "alien"
   26   14    2    2.59  349.26  0.6070E+01  0.7578E+01 "alien"
   15   11   12    2.54  231.48  0.4023E+01  0.4158E+01 "alien"
   30    3    3    2.51  209.65  0.3643E+01  0.6376E+01 "alien"
   16   11   11    2.63  177.12  0.3078E+01  0.6226E+01 "alien"
   24   13    4    2.71  167.10  0.2904E+01  0.6177E+01 "alien"
   22    1   10    2.60  154.13  0.2679E+01  0.4766E+01 "alien"
   25   14    4    2.60  144.68  0.2514E+01  0.6266E+01 "alien"
   12   10   13    2.55  126.44  0.2197E+01  0.3488E+01 "alien"
   28   19   23    1.35  104.24  0.1694E+02  0.8312E+01 "alien"
   35   33   15    1.36   99.26  0.1613E+02  0.1126E+02 "alien"
   52    6   12    1.34   89.67  0.1457E+02  0.7933E+01 "alien"
   39   28    4    1.59   74.41  0.1550E+02  0.1231E+02 "alien"
   28    3    5    2.58   73.86  0.1284E+01  0.1009E+02 "alien"
   27   15    3    2.45   63.70  0.1703E+02  0.7738E+01 "alien"
   45    8   17    1.36   61.68  0.1002E+02  0.7423E+01 "alien"
   27   24   16    1.60   57.48  0.1197E+02  0.1106E+02 "alien"
   37    5   16    1.57   56.93  0.1186E+02  0.1010E+02 "alien"
   29   20   23    1.33   53.22  0.8645E+01  0.8326E+01 "alien"
   53    7   13    1.30   51.84  0.5049E+02  0.1123E+02 "alien"
   40    2    1    1.93   49.03 -0.1596E+02  0.8237E+01 "alien"
   27    2    8    2.46   46.84  0.1252E+02  0.8424E+01 "alien"
   38   34   13    1.35   46.26  0.7515E+01  0.7789E+01 "alien"
   10    9   13    2.64   44.77  0.7780E+00  0.3411E+01 "alien"
   32    3    2    2.38   41.24  0.1102E+02  0.6748E+01 "alien"
   39   28    1    1.61   39.45  0.8220E+01  0.7095E+01 "alien"
   46    1    6    1.62   38.73  0.8069E+01  0.9925E+01 "alien"
   49    7   15    1.33   38.60  0.6271E+01  0.8048E+01 "alien"
   39    1    6    1.89   38.27 -0.1246E+02  0.6706E+01 "alien"
   32    3    3    2.36   38.14  0.1020E+02  0.6799E+01 "alien"
   25   24   17    1.59   37.44  0.7800E+01  0.6097E+01 "alien"
   28   16    3    2.35   37.25  0.9958E+01  0.6901E+01 "alien"
   26   14    0    2.61   37.22  0.1297E+02  0.9072E+01 "alien"
   34   27   11    1.59   36.20  0.7543E+01  0.6773E+01 "alien"
   33   22    5    1.89   35.89 -0.1168E+02  0.7169E+01 "alien"
   21   15   21    1.57   35.77  0.7452E+01  0.7671E+01 "alien"
   39    2    2    1.97   33.62 -0.1094E+02  0.7032E+01 "alien"
   31    3    4    2.40   32.60  0.8714E+01  0.1069E+02 "alien"
   45    2   10    1.57   31.93  0.6651E+01  0.7392E+01 "alien"
   48    1    5    1.57   31.72  0.6609E+01  0.7216E+01 "alien"
   38   28    6    1.59   31.37  0.6536E+01  0.9963E+01 "alien"
   21   17   14    1.99   31.06 -0.1011E+02  0.5508E+01 "alien"
   25   18   12    1.98   30.50 -0.9930E+01  0.8613E+01 "alien"
   35    2    6    2.08   30.11 -0.1524E+02  0.1166E+02 "alien"
   41    3   13    1.59   29.95  0.6239E+01  0.6424E+01 "alien"
   28   16    1    2.39   29.79  0.7963E+01  0.6721E+01 "alien"
   28   24   15    1.62   29.37  0.6118E+01  0.6160E+01 "alien"
   19   13   11    2.43   28.34  0.7577E+01  0.5817E+01 "alien"
   17    0   11    2.78   28.18  0.4149E+01  0.6002E+01 "alien"
    8    8   13    2.73   27.86  0.3988E+01  0.3572E+01 "alien"
   28    2    7    2.46   27.82  0.7437E+01  0.6240E+01 "alien"
   39   28    3    1.60   27.53  0.5736E+01  0.7040E+01 "alien"
   31    3    1    2.47   26.24  0.7016E+01  0.6534E+01 "alien"
   24    1   10    2.47   25.81  0.6900E+01  0.5118E+01 "alien"
   38    8   19    1.42   25.04  0.2040E+02  0.6757E+01 "alien"
   43    1    2    1.79   25.01  0.1427E+02  0.7170E+01 "alien"
   33   23    7    1.81   24.76  0.1412E+02  0.1075E+02 "alien"
   40    1    4    1.89   23.94 -0.7793E+01  0.1260E+02 "alien"
   35   24    4    1.79   23.88  0.1362E+02  0.8473E+01 "alien"
   20   13   10    2.48   23.70  0.6336E+01  0.5426E+01 "alien"
   21   18   15    1.89   23.17 -0.7543E+01  0.5244E+01 "alien"
   25   15   23    1.42   22.39  0.1824E+02  0.8128E+01 "alien"
   34    5   17    1.60   22.36  0.4658E+01  0.5787E+01 "alien"
   49    4   10    1.46   22.33  0.2509E+02  0.1205E+02 "alien"
   39   32   10    1.42   22.33  0.1819E+02  0.7812E+01 "alien"
   45    1    8    1.62   21.97  0.4578E+01  0.7357E+01 "alien"
   26   16    7    2.30   21.40  0.5722E+01  0.6543E+01 "alien"
   35    1   12    1.82   21.10  0.1204E+02  0.6370E+01 "alien"
   37    2    2    2.07   21.10 -0.1068E+02  0.7659E+01 "alien"
   31   19    3    2.10   21.06 -0.1066E+02  0.7248E+01 "alien"
   39    2    3    1.95   20.74 -0.6751E+01  0.8356E+01 "alien"
   51    6   12    1.36   20.20  0.3282E+01  0.8273E+01 "alien"
   41    9   19    1.36   20.16  0.3275E+01  0.7264E+01 "alien"
   38    4   15    1.59   19.98  0.4163E+01  0.6043E+01 "alien"
   36    1    8    1.96   19.89 -0.6475E+01  0.6533E+01 "alien"
   48    1    3    1.60   19.41  0.4043E+01  0.6993E+01 "alien"
   42    1    4    1.81   19.23  0.1097E+02  0.7063E+01 "alien"
   43    2   11    1.60   19.22  0.4003E+01  0.7787E+01 "alien"
   30   21    9    1.89   18.47 -0.6012E+01  0.7624E+01 "alien"
   34   22    1    1.90   18.41 -0.5992E+01  0.6847E+01 "alien"
   23   23   18    1.59   18.11  0.1353E+02  0.8066E+01 "alien"
   29   22   22    1.35   17.78  0.2889E+01  0.6267E+01 "alien"
   14   11   13    2.45   17.33  0.4632E+01  0.6448E+01 "alien"
   46    1    7    1.61   17.20  0.3583E+01  0.6944E+01 "alien"
   38   29    8    1.53   17.00 -0.1877E+02  0.7047E+01 "alien"
   19   12   10    2.57   16.98  0.2951E+00  0.5454E+01 "alien"
   31    3    3    2.43   16.74  0.4476E+01  0.6382E+01 "alien"
   32   20    2    2.03   16.61 -0.8407E+01  0.7181E+01 "alien"
   30   19    6    2.06   16.52 -0.8360E+01  0.7091E+01 "alien"
   20   18   16    1.85   16.43  0.9373E+01  0.7106E+01 "alien"
   33   32   18    1.32   16.22  0.1580E+02  0.1102E+02 "alien"
   43    5   15    1.47   15.77  0.1771E+02  0.7990E+01 "alien"
   27    1   13    2.06   15.49 -0.7842E+01  0.5487E+01 "alien"
   31   25   13    1.62   15.41  0.3210E+01  0.6205E+01 "alien"
   37   32   12    1.42   15.33  0.1249E+02  0.7495E+01 "alien"
   34    2   13    1.80   15.29  0.8725E+01  0.1009E+02 "alien"
   42    9   19    1.35   15.13  0.2458E+01  0.7064E+01 "alien"
   32   30   19    1.33   14.68  0.2386E+01  0.1042E+02 "alien"
   25   19   13    1.89   14.44 -0.4702E+01  0.5882E+01 "alien"
   48    9   17    1.30   14.27  0.1389E+02  0.1089E+02 "alien"
   35   23    1    1.84   14.26  0.8135E+01  0.6772E+01 "alien"
   21   20   18    1.67   14.13  0.1793E+02  0.8573E+01 "alien"
   34    7   19    1.50   14.13  0.1587E+02  0.8925E+01 "alien"
   33   29   14    1.48   14.01  0.1574E+02  0.7254E+01 "alien"
   21    9   21    1.62   13.97  0.2911E+01  0.5713E+01 "alien"
   36   34   15    1.33   13.89  0.2256E+01  0.1045E+02 "alien"
   38   28    5    1.60   13.75  0.2864E+01  0.7226E+01 "alien"
   31   22    9    1.84   13.54  0.7725E+01  0.6852E+01 "alien"
   41    1    4    1.85   13.39  0.7636E+01  0.7021E+01 "alien"
   31   27   15    1.52   13.32 -0.1471E+02  0.1119E+02 "alien"
   52    3    6    1.44   13.08  0.1066E+02  0.7681E+01 "alien"
   31    3    0    2.48   13.04  0.8595E+01  0.9358E+01 "alien"
   39   28    2    1.60   12.95  0.2697E+01  0.7084E+01 "alien"
   34   22    0    1.91   12.66 -0.9897E+01  0.9658E+01 "alien"
   25   20   21    1.46   12.65  0.1422E+02  0.8549E+01 "alien"
   16    2   14    2.39   12.43  0.3323E+01  0.3582E+01 "alien"
   39   28    0    1.61   12.41  0.6107E+01  0.9991E+01 "alien"
   36   32   13    1.41   12.21  0.9951E+01  0.1139E+02 "alien"
   13    5   17    2.10   12.15 -0.6148E+01  0.3754E+01 "alien"
   32   26   13    1.58   12.11  0.2523E+01  0.6641E+01 "alien"
   26    3    5    2.76   12.02  0.2090E+00  0.5905E+01 "alien"
   32   28   14    1.52   11.85 -0.1308E+02  0.1187E+02 "alien"
   25    1   10    2.41   11.71  0.3132E+01  0.5295E+01 "alien"
   20    1   13    2.35   11.54  0.3086E+01  0.4237E+01 "alien"
   36   26    6    1.68   11.47  0.1456E+02  0.1044E+02 "alien"
   30    3    4    2.48   11.40  0.3047E+01  0.6459E+01 "alien"
   39   11   21    1.32   11.36  0.1107E+02  0.1124E+02 "alien"
   26   20   13    1.85   11.30  0.6448E+01  0.5892E+01 "alien"
   43    3   12    1.57   11.26  0.2345E+01  0.6585E+01 "alien"
   26   21   21    1.44   11.09  0.9036E+01  0.6676E+01 "alien"
   41    8   18    1.40   11.09  0.1607E+02  0.1213E+02 "alien"
   54    3    3    1.42   10.92  0.8899E+01  0.7680E+01 "alien"
   37    0    9    1.88   10.88 -0.7551E+01  0.6390E+01 "alien"
   46    2    9    1.56   10.70 -0.1181E+02  0.6891E+01 "alien"
   27   16    6    2.30   10.40  0.2779E+01  0.6578E+01 "alien"
   30   26   20    1.37   10.36  0.1502E+02  0.1136E+02 "alien"
   40   30    4    1.52   10.22 -0.1129E+02  0.7275E+01 "alien"
   35   34   16    1.32   10.18  0.9917E+01  0.1060E+02 "alien"
   37    2    3    2.06   10.10 -0.5112E+01  0.7291E+01 "alien"
   29    2   14    1.91   10.06 -0.3274E+01  0.5407E+01 "alien"
   22   11   23    1.48   10.05  0.1129E+02  0.5908E+01 "alien"


 List of    146 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   48    9   17    1.30   14.27  0.1389E+02  0.1089E+02
   53    7   13    1.30   51.84  0.5049E+02  0.1123E+02
   33   32   18    1.32   16.22  0.1580E+02  0.1102E+02
   35   34   16    1.32   10.18  0.9917E+01  0.1060E+02
   39   11   21    1.32   11.36  0.1107E+02  0.1124E+02
   32   30   19    1.33   14.68  0.2386E+01  0.1042E+02
   29   20   23    1.33   53.22  0.8645E+01  0.8326E+01
   36   34   15    1.33   13.89  0.2256E+01  0.1045E+02
   49    7   15    1.33   38.60  0.6271E+01  0.8048E+01
   52    6   12    1.34   89.67  0.1457E+02  0.7933E+01
   42    9   19    1.35   15.13  0.2458E+01  0.7064E+01
   38   34   13    1.35   46.26  0.7515E+01  0.7789E+01
   29   22   22    1.35   17.78  0.2889E+01  0.6267E+01
   28   19   23    1.35  104.24  0.1694E+02  0.8312E+01
   45    8   17    1.36   61.68  0.1002E+02  0.7423E+01
   35   33   15    1.36   99.26  0.1613E+02  0.1126E+02
   51    6   12    1.36   20.20  0.3282E+01  0.8273E+01
   41    9   19    1.36   20.16  0.3275E+01  0.7264E+01
   30   26   20    1.37   10.36  0.1502E+02  0.1136E+02
   41    8   18    1.40   11.09  0.1607E+02  0.1213E+02
   36   32   13    1.41   12.21  0.9951E+01  0.1139E+02
   37   32   12    1.42   15.33  0.1249E+02  0.7495E+01
   25   15   23    1.42   22.39  0.1824E+02  0.8128E+01
   54    3    3    1.42   10.92  0.8899E+01  0.7680E+01
   38    8   19    1.42   25.04  0.2040E+02  0.6757E+01
   39   32   10    1.42   22.33  0.1819E+02  0.7812E+01
   26   21   21    1.44   11.09  0.9036E+01  0.6676E+01
   52    3    6    1.44   13.08  0.1066E+02  0.7681E+01
   49    4   10    1.46   22.33  0.2509E+02  0.1205E+02
   25   20   21    1.46   12.65  0.1422E+02  0.8549E+01
   43    5   15    1.47   15.77  0.1771E+02  0.7990E+01
   22   11   23    1.48   10.05  0.1129E+02  0.5908E+01
   33   29   14    1.48   14.01  0.1574E+02  0.7254E+01
   34    7   19    1.50   14.13  0.1587E+02  0.8925E+01
   32   28   14    1.52   11.85 -0.1308E+02  0.1187E+02
   31   27   15    1.52   13.32 -0.1471E+02  0.1119E+02
   40   30    4    1.52   10.22 -0.1129E+02  0.7275E+01
   38   29    8    1.53   17.00 -0.1877E+02  0.7047E+01
   46    2    9    1.56   10.70 -0.1181E+02  0.6891E+01
   43    3   12    1.57   11.26  0.2345E+01  0.6585E+01
   45    2   10    1.57   31.93  0.6651E+01  0.7392E+01
   21   15   21    1.57   35.77  0.7452E+01  0.7671E+01
   37    5   16    1.57   56.93  0.1186E+02  0.1010E+02
   48    1    5    1.57   31.72  0.6609E+01  0.7216E+01
   32   26   13    1.58   12.11  0.2523E+01  0.6641E+01
   38   28    6    1.59   31.37  0.6536E+01  0.9963E+01
   39   28    4    1.59   74.41  0.1550E+02  0.1231E+02
   34   27   11    1.59   36.20  0.7543E+01  0.6773E+01
   41    3   13    1.59   29.95  0.6239E+01  0.6424E+01
   38    4   15    1.59   19.98  0.4163E+01  0.6043E+01
   25   24   17    1.59   37.44  0.7800E+01  0.6097E+01
   23   23   18    1.59   18.11  0.1353E+02  0.8066E+01
   48    1    3    1.60   19.41  0.4043E+01  0.6993E+01
   39   28    3    1.60   27.53  0.5736E+01  0.7040E+01
   43    2   11    1.60   19.22  0.4003E+01  0.7787E+01
   34    5   17    1.60   22.36  0.4658E+01  0.5787E+01
   38   28    5    1.60   13.75  0.2864E+01  0.7226E+01
   27   24   16    1.60   57.48  0.1197E+02  0.1106E+02
   39   28    2    1.60   12.95  0.2697E+01  0.7084E+01
   46    1    7    1.61   17.20  0.3583E+01  0.6944E+01
   39   28    1    1.61   39.45  0.8220E+01  0.7095E+01
   39   28    0    1.61   12.41  0.6107E+01  0.9991E+01
   21    9   21    1.62   13.97  0.2911E+01  0.5713E+01
   45    1    8    1.62   21.97  0.4578E+01  0.7357E+01
   31   25   13    1.62   15.41  0.3210E+01  0.6205E+01
   46    1    6    1.62   38.73  0.8069E+01  0.9925E+01
   28   24   15    1.62   29.37  0.6118E+01  0.6160E+01
   21   20   18    1.67   14.13  0.1793E+02  0.8573E+01
   36   26    6    1.68   11.47  0.1456E+02  0.1044E+02
   43    1    2    1.79   25.01  0.1427E+02  0.7170E+01
   35   24    4    1.79   23.88  0.1362E+02  0.8473E+01
   34    2   13    1.80   15.29  0.8725E+01  0.1009E+02
   42    1    4    1.81   19.23  0.1097E+02  0.7063E+01
   33   23    7    1.81   24.76  0.1412E+02  0.1075E+02
   35    1   12    1.82   21.10  0.1204E+02  0.6370E+01
   31   22    9    1.84   13.54  0.7725E+01  0.6852E+01
   35   23    1    1.84   14.26  0.8135E+01  0.6772E+01
   26   20   13    1.85   11.30  0.6448E+01  0.5892E+01
   20   18   16    1.85   16.43  0.9373E+01  0.7106E+01
   41    1    4    1.85   13.39  0.7636E+01  0.7021E+01
   37    0    9    1.88   10.88 -0.7551E+01  0.6390E+01
   33   22    5    1.89   35.89 -0.1168E+02  0.7169E+01
   39    1    6    1.89   38.27 -0.1246E+02  0.6706E+01
   40    1    4    1.89   23.94 -0.7793E+01  0.1260E+02
   30   21    9    1.89   18.47 -0.6012E+01  0.7624E+01
   21   18   15    1.89   23.17 -0.7543E+01  0.5244E+01
   25   19   13    1.89   14.44 -0.4702E+01  0.5882E+01
   34   22    1    1.90   18.41 -0.5992E+01  0.6847E+01
   34   22    0    1.91   12.66 -0.9897E+01  0.9658E+01
   29    2   14    1.91   10.06 -0.3274E+01  0.5407E+01
   40    2    1    1.93   49.03 -0.1596E+02  0.8237E+01
   39    2    3    1.95   20.74 -0.6751E+01  0.8356E+01
   36    1    8    1.96   19.89 -0.6475E+01  0.6533E+01
   39    2    2    1.97   33.62 -0.1094E+02  0.7032E+01
   25   18   12    1.98   30.50 -0.9930E+01  0.8613E+01
   21   17   14    1.99   31.06 -0.1011E+02  0.5508E+01
   32   20    2    2.03   16.61 -0.8407E+01  0.7181E+01
   37    2    3    2.06   10.10 -0.5112E+01  0.7291E+01
   27    1   13    2.06   15.49 -0.7842E+01  0.5487E+01
   30   19    6    2.06   16.52 -0.8360E+01  0.7091E+01
   37    2    2    2.07   21.10 -0.1068E+02  0.7659E+01
   35    2    6    2.08   30.11 -0.1524E+02  0.1166E+02
   31   19    3    2.10   21.06 -0.1066E+02  0.7248E+01
   13    5   17    2.10   12.15 -0.6148E+01  0.3754E+01
   27   16    6    2.30   10.40  0.2779E+01  0.6578E+01
   26   16    7    2.30   21.40  0.5722E+01  0.6543E+01
   28   16    3    2.35   37.25  0.9958E+01  0.6901E+01
   20    1   13    2.35   11.54  0.3086E+01  0.4237E+01
   32    3    3    2.36   38.14  0.1020E+02  0.6799E+01
   32    3    2    2.38   41.24  0.1102E+02  0.6748E+01
   28   16    1    2.39   29.79  0.7963E+01  0.6721E+01
   16    2   14    2.39   12.43  0.3323E+01  0.3582E+01
   31    3    4    2.40   32.60  0.8714E+01  0.1069E+02
   25    1   10    2.41   11.71  0.3132E+01  0.5295E+01
   19   13   11    2.43   28.34  0.7577E+01  0.5817E+01
   31    3    3    2.43   16.74  0.4476E+01  0.6382E+01
   14   11   13    2.45   17.33  0.4632E+01  0.6448E+01
   27   15    3    2.45   63.70  0.1703E+02  0.7738E+01
   27    2    8    2.46   46.84  0.1252E+02  0.8424E+01
   28    2    7    2.46   27.82  0.7437E+01  0.6240E+01
   24    1   10    2.47   25.81  0.6900E+01  0.5118E+01
   31    3    1    2.47   26.24  0.7016E+01  0.6534E+01
   30    3    4    2.48   11.40  0.3047E+01  0.6459E+01
   31    3    0    2.48   13.04  0.8595E+01  0.9358E+01
   20   13   10    2.48   23.70  0.6336E+01  0.5426E+01
   30    3    3    2.51  209.65  0.3643E+01  0.6376E+01
   27    2    7    2.53 1064.48  0.1850E+02  0.9586E+01
   15   11   12    2.54  231.48  0.4023E+01  0.4158E+01
   12   10   13    2.55  126.44  0.2197E+01  0.3488E+01
   19   12   10    2.57   16.98  0.2951E+00  0.5454E+01
   28    3    5    2.58   73.86  0.1284E+01  0.1009E+02
   26   14    2    2.59  349.26  0.6070E+01  0.7578E+01
   29    3    3    2.59  521.71  0.9067E+01  0.6318E+01
   25   14    4    2.60  144.68  0.2514E+01  0.6266E+01
   22    1   10    2.60  154.13  0.2679E+01  0.4766E+01
   26   14    1    2.61  517.90  0.9000E+01  0.6363E+01
   26   14    0    2.61   37.22  0.1297E+02  0.9072E+01
   16   11   11    2.63  177.12  0.3078E+01  0.6226E+01
   10    9   13    2.64   44.77  0.7780E+00  0.3411E+01
   21    1   10    2.66  393.18  0.6833E+01  0.4776E+01
   24   13    4    2.71  167.10  0.2904E+01  0.6177E+01
    8    8   13    2.73   27.86  0.3988E+01  0.3572E+01
   20    1   10    2.73  485.14  0.8431E+01  0.6363E+01
   26    3    5    2.76   12.02  0.2090E+00  0.5905E+01
   24   13    3    2.76  375.37  0.6523E+01  0.6065E+01
   17    0   11    2.78   28.18  0.4149E+01  0.6002E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        72
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1528
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1381
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      147


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
