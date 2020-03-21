
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.509    77.509    38.788  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      27      27
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      27
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000911  0.000679
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
 REFLECTING_RANGE_E.S.D.=   0.131
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

    1337 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      27 ...      27
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1337 REFLECTIONS ACCEPTED



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
 *  25        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0   -1  1  0  0 -1 -1  0  0  0  0  1  0
 *  21        tP          0.0      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.8   38.8   77.5  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.8   38.8   77.5  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  159.8   77.5  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.8   77.5  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  159.8   77.5  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  159.8   77.5  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     159.8   38.8  109.6  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.6  116.3  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  159.8  159.8  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.6  116.3   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        998.3     173.3   77.5   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.3      77.5  173.3   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.3      77.5  173.3   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.3      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.3  116.3  116.3  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.7   86.7  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.6   86.7   86.7  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.7   86.7  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.7   86.7  109.6  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.7  248.2 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.7   86.7   77.5  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.5  177.6  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.7   86.7  177.6 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.7   86.7   86.7  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.6  109.6   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.5  177.6 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.7  177.6   77.5 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    4.2      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.861   -23.503     0.041
 UNIT_CELL_B-AXIS=    23.502    73.861     0.019
 UNIT_CELL_C-AXIS=    -0.022    -0.003    38.788
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       194 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.26
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.471    77.471    38.774  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.5E-01 1.5E-01 6.9E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012908  0.012908  0.025790  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.830   -23.468     0.049
 COORDINATES OF UNIT CELL B-AXIS    23.468    73.830     0.018
 COORDINATES OF UNIT CELL C-AXIS    -0.026    -0.001    38.774
 CRYSTAL MOSAICITY (DEGREES)     0.131
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006365  0.001978
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.92   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1261.04   1242.15
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.28
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.47    77.47    38.77  90.000  90.000  90.000


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
  0.05  0.15        59         94.985       1.393
  0.15  0.25       169        141.199       0.937
  0.25  0.35       186        119.357       1.108
  0.35  0.45       154         86.630       1.527
  0.45  0.55       170        174.896       0.756
  0.55  0.65       188        123.512       1.071
  0.65  0.75       197        167.034       0.792
  0.75  0.85       148        128.493       1.030
  0.85  0.95        66        105.931       1.249
 
 -0.05  0.95      1337        132.292       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               208.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.902
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    26.1 XMAX=    26.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59418 NYBIN=    9
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        450


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               208.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.900
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    10.4 XMAX=  2453.8 NXBIN=    3
 YMIN=     8.0 YMAX=  2520.1 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        450


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               208.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.892
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    26.1 XMAX=    26.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1647    1691
 DETECTOR_SURFACE_POSITION=     817    1691
 DETECTOR_SURFACE_POSITION=     817     837
 DETECTOR_SURFACE_POSITION=    1647     837
 DETECTOR_SURFACE_POSITION=    2171    1664
 DETECTOR_SURFACE_POSITION=    1621    2230
 DETECTOR_SURFACE_POSITION=     843    2230
 DETECTOR_SURFACE_POSITION=     293    1664
 DETECTOR_SURFACE_POSITION=     293     864
 DETECTOR_SURFACE_POSITION=     843     299
 DETECTOR_SURFACE_POSITION=    1621     299
 DETECTOR_SURFACE_POSITION=    2171     864
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        450


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
 7.840E-01  1.175E-03   32.94


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

   24.498   3.156    24.46   1.38      2.13      1.89      64      87      28
    3.156   2.241     9.13   1.12      5.44      5.06     186     214       6
    2.241   1.832     1.22   1.10     41.86     39.28     198     228       2
    1.832   1.588    -0.04   0.91    -99.90    -99.90     240     273       2
    1.588   1.421     0.01   0.84    -99.90    -99.90     270     322       0
    1.421   1.297     0.01   0.75    792.50    892.84     128     174       0
   --------------------------------------------------------------------------
   24.498   1.297     3.36   0.97      7.49      7.38    1086    1298      38
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      27      27)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.684   0.7790E+01  0.1140E+02         0.68           1 
   47    0    0     1.648  -0.9269E+01  0.1117E+02        -0.83           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -7.0% OF MEAN INTENSITY


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

     3.85          53      33      1285        2.6%       2.0%      1.8%       40   33.31      2.8%    99.9*     0    0.000       0
     2.74          91      53      2098        2.5%       4.0%      3.0%       76   20.84      5.7%    99.5*     0    0.000       0
     2.24         157      90      2665        3.4%       6.2%      6.5%      134   10.48      8.7%    99.0*     0    0.000       0
     1.94         146      83      3101        2.7%      27.6%     26.0%      126    2.48     39.0%    89.9*     0    0.000       0
     1.74         155      86      3481        2.5%    1266.5%   1227.8%      138   -0.00   1791.1%     2.0      0    0.000       0
     1.59         196     111      3847        2.9%    2014.2%   2162.7%      170   -0.06   2848.6%   -18.3      0    0.000       0
     1.47         207     121      4142        2.9%   10273.2%  10517.2%      172    0.06  14528.5%    -9.3      0    0.000       0
     1.38         203     118      4433        2.7%     -99.9%    -99.9%      170    0.00    -99.9%    15.3      0    0.000       0
     1.30          84      57      4741        1.2%     128.0%    144.6%       54    0.20    181.0%    64.3*     0    0.000       0
    total        1292     752     29793        2.5%       7.4%      7.3%     1080    4.46     10.5%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1337
 NUMBER OF REJECTED MISFITS                                38
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1298
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    755



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

   WILSON LINE (using all data) : A=   7.234 B=  28.697 CORRELATION=  0.93
      #      RES      SS        <I>       log(<I>)       BO
      14     7.398   0.005  1.9426E+03       7.572     -37.0
      15     4.683   0.011  2.2648E+03       7.725     -21.6
      23     3.518   0.020  7.6388E+02       6.638      14.7
      28     2.993   0.028  2.3209E+02       5.447      32.0
      33     2.643   0.036  2.0405E+02       5.318      26.8
      49     2.402   0.043  1.0904E+02       4.692      29.3
      37     2.211   0.051  5.2541E+01       3.962      32.0
      41     2.060   0.059  2.0183E+01       3.005      35.9
      40     1.935   0.067 -2.7235E-01    -999.999   -1000.0
      48     1.831   0.075 -1.0657E+00    -999.999   -1000.0
      42     1.739   0.083  7.9249E-01      -0.233      45.2
      49     1.669   0.090 -6.4047E-01    -999.999   -1000.0
      50     1.599   0.098 -8.0474E-01    -999.999   -1000.0
      60     1.538   0.106  2.2424E+00       0.808      30.4
      69     1.482   0.114 -3.4868E-01    -999.999   -1000.0
      53     1.435   0.121  6.6930E-01      -0.402      31.4
      46     1.393   0.129 -1.5266E+00    -999.999   -1000.0
      39     1.353   0.137  1.3895E+00       0.329      25.3
      17     1.316   0.144  2.1106E+00       0.747      22.5


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.398        1.914        1.385        0.784
       3     4.683        1.127        0.653        0.281
       2     3.518        0.453        0.106        0.018
       3     2.993        2.106        1.165        0.482
       2     2.643        0.406        0.127        0.028
       5     2.402        0.745        0.380        0.152
       1     2.211        1.051        0.373        0.095
       6     2.060        0.615        0.285        0.110
      10     1.831       31.578       42.202      106.715
       7     1.739       29.556       61.793      136.271
       3     1.538        7.660        9.274       11.095
       3     1.435       48.337      194.991      572.649
       1     1.353       11.494       13.498       11.323
      77   overall        5.376       13.952       35.370


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.398        0.538        0.394        0.227
      12     4.683        1.173        1.332        1.197
      21     3.518        0.965        0.954        0.873
      25     2.993        0.741        0.571        0.366
      31     2.643        1.264        1.709        2.154
      44     2.402        0.975        0.857        0.634
      36     2.211        0.888        0.731        0.533
      35     2.060        1.298        1.579        1.773
      34     1.935       52.359      -37.217     1430.153
      57     1.538        5.801        5.500       17.046
      50     1.435      108.896       32.773     5079.101
      38     1.353       99.384      843.840    17118.220
      17     1.316       15.938       23.077       96.306
     676   overall       17.288       44.502     1425.816

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.398  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.38  0.73
       3     4.683  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       2     3.518  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.993  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.643  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.402  0.00  0.00  0.00  0.00  0.38  1.07  1.00  0.95  0.91  0.88
       1     2.211  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.060  0.67  0.48  0.40  0.35  0.96  0.89  0.84  0.80  1.01  0.98
      10     1.831  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
       7     1.739  0.58  0.41  0.34  0.30  0.27  0.25  0.24  0.23  0.22  0.21
       3     1.538  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.435  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     1.353  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      77   overall  0.63  0.45  0.37  0.33  0.37  0.39  0.37  0.37  0.40  0.40


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.398  3.15  3.31  2.31  1.82  1.78  1.55  1.59  1.45  1.35  1.27
      12     4.683  1.75  1.84  1.29  1.01  1.69  1.48  1.32  1.51  1.40  1.32
      21     3.518  1.00  0.79  0.92  0.72  0.73  0.84  0.95  1.21  1.12  1.05
      25     2.993  1.26  1.99  1.39  1.21  1.42  1.42  1.27  1.23  1.21  1.20
      31     2.643  1.02  1.07  1.00  1.17  1.07  1.00  1.09  1.05  1.14  1.07
      44     2.402  0.96  0.88  0.96  0.90  1.27  1.16  1.08  1.07  1.07  1.08
      36     2.211  1.46  0.77  0.96  0.76  0.78  0.86  0.88  0.96  1.03  1.01
      35     2.060  2.10  1.26  1.10  1.04  0.94  0.95  0.96  1.04  1.01  1.04
      34     1.935  6.18  3.25  2.27  1.78  1.49  1.30  1.17  1.07  0.99  0.93
      57     1.538  4.42  2.32  1.69  1.38  1.16  1.01  0.94  0.89  0.83  0.80
      50     1.435  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.82
      38     1.353  4.70  2.47  1.73  1.44  1.20  1.05  0.94  0.86  0.80  0.79
      17     1.316  4.95  2.60  1.82  1.43  1.20  1.17  1.05  0.96  0.89  0.84
     676   overall  1.91  1.16  0.90  0.75  0.72  0.67  0.63  0.62  0.60  0.59


 List of     91 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   30    1   13    1.95   92.23 -0.2512E+02  0.1159E+02 "alien"
   28    2   14    1.95   70.99 -0.1933E+02  0.1170E+02 "alien"
   54    5    9    1.36   61.00  0.8476E+02  0.9733E+01 "alien"
   39    2    3    1.96   60.87 -0.1658E+02  0.1006E+02 "alien"
   33    1   12    1.90   56.84 -0.1548E+02  0.7743E+01 "alien"
   38   30    9    1.50   39.39 -0.1373E+02  0.8497E+01 "alien"
   52    3    5    1.46   37.60 -0.1311E+02  0.1590E+02 "alien"
   37   30   11    1.48   37.52 -0.1308E+02  0.1014E+02 "alien"
   26   19   12    1.93   36.70 -0.9995E+01  0.7272E+01 "alien"
   41   31    3    1.50   35.56 -0.1240E+02  0.8504E+01 "alien"
   34   29   13    1.50   35.55 -0.1240E+02  0.8915E+01 "alien"
   39   31    9    1.46   35.44 -0.1236E+02  0.9765E+01 "alien"
   29    9   21    1.50   34.46 -0.1201E+02  0.6206E+01 "alien"
   30   28   16    1.49   33.20 -0.1158E+02  0.1230E+02 "alien"
   37   32   12    1.42   32.10  0.2148E+02  0.1021E+02 "alien"
   37   31   11    1.46   31.13 -0.1085E+02  0.1089E+02 "alien"
   15   12   18    1.90   27.38 -0.7457E+01  0.4860E+01 "alien"
   54    3    1    1.43   27.22  0.1822E+02  0.9338E+01 "alien"
   24   19   14    1.87   26.33 -0.2806E+02  0.1141E+02 "alien"
   52    3    7    1.44   26.28  0.1759E+02  0.9742E+01 "alien"
   29    2   14    1.92   26.19 -0.7132E+01  0.7289E+01 "alien"
   41   32    5    1.46   25.99 -0.9061E+01  0.8913E+01 "alien"
   33    8   20    1.48   25.75 -0.8977E+01  0.1114E+02 "alien"
   40   31    6    1.49   25.06 -0.8740E+01  0.8573E+01 "alien"
   30   20    8    1.96   23.89 -0.6505E+01  0.7827E+01 "alien"
   38   31   10    1.46   23.76 -0.8284E+01  0.8674E+01 "alien"
   41   31    4    1.49   23.39 -0.8157E+01  0.9001E+01 "alien"
   29    1   13    1.99   22.94 -0.6248E+01  0.7301E+01 "alien"
   51    6   12    1.37   22.92  0.3185E+02  0.1684E+02 "alien"
   36   30   12    1.47   22.37 -0.7798E+01  0.8560E+01 "alien"
   53    3    5    1.43   22.25  0.1489E+02  0.8988E+01 "alien"
   41   32    6    1.45   21.64  0.1448E+02  0.9010E+01 "alien"
   35    2   13    1.78   21.29  0.1687E+02  0.7775E+01 "alien"
   50    2    6    1.51   20.76 -0.7240E+01  0.8522E+01 "alien"
   21   18   15    1.90   19.37 -0.5275E+01  0.7139E+01 "alien"
   43    1    8    1.69   19.32 -0.1238E+02  0.1412E+02 "alien"
   34    0   10    1.96   19.08 -0.2011E+02  0.9108E+01 "alien"
   39   29    4    1.57   18.94 -0.1524E+02  0.9579E+01 "alien"
   53    2    1    1.46   18.41 -0.6418E+01  0.1072E+02 "alien"
   52    4    8    1.42   18.07  0.1210E+02  0.9730E+01 "alien"
   33   23    8    1.79   18.05 -0.1923E+02  0.7558E+01 "alien"
   36    9   20    1.42   17.98  0.1203E+02  0.1321E+02 "alien"
   40   31    7    1.48   17.48 -0.6094E+01  0.9190E+01 "alien"
   52    2    1    1.49   17.24 -0.6013E+01  0.8530E+01 "alien"
   31   22   10    1.80   17.24 -0.1838E+02  0.1247E+02 "alien"
   25   24   17    1.60   16.71 -0.1345E+02  0.1050E+02 "alien"
   48    5   13    1.41   16.30  0.1091E+02  0.9475E+01 "alien"
   47    4   12    1.46   16.00 -0.5580E+01  0.8581E+01 "alien"
   38    8   19    1.43   15.98  0.1069E+02  0.7762E+01 "alien"
   28   26   19    1.44   15.64  0.1047E+02  0.1179E+02 "alien"
   34   22    0    1.91   15.35 -0.1215E+02  0.1194E+02 "alien"
   18   10   21    1.66   14.47 -0.9268E+01  0.6760E+01 "alien"
   27   13   24    1.37   14.46  0.2009E+02  0.1005E+02 "alien"
   52    2    0    1.49   14.39 -0.1359E+02  0.1242E+02 "alien"
   34   22    1    1.91   14.33 -0.3902E+01  0.8497E+01 "alien"
   52    4    9    1.40   14.28 -0.2179E+02  0.9031E+01 "alien"
   32    8   20    1.50   14.04 -0.4894E+01  0.7302E+01 "alien"
   34    4   16    1.65   13.55 -0.8681E+01  0.7067E+01 "alien"
   43    1    2    1.79   13.44 -0.1433E+02  0.7949E+01 "alien"
   37   10   21    1.36   13.11  0.1822E+02  0.1212E+02 "alien"
   40   32    8    1.44   13.10  0.8768E+01  0.8719E+01 "alien"
   27   24   16    1.61   13.03 -0.1049E+02  0.1314E+02 "alien"
   45    2   10    1.57   12.87 -0.1035E+02  0.8798E+01 "alien"
   29   24   14    1.65   12.84 -0.8221E+01  0.7310E+01 "alien"
   40    2    1    1.93   12.35 -0.3365E+01  0.9909E+01 "alien"
   29   21   11    1.84   12.34 -0.1315E+02  0.7278E+01 "alien"
   36    3   14    1.70   12.32 -0.7893E+01  0.1198E+02 "alien"
   42    1    9    1.70   12.32 -0.7888E+01  0.7535E+01 "alien"
   41   31    1    1.51   12.08 -0.4212E+01  0.8384E+01 "alien"
   49    3    9    1.48   11.98 -0.4176E+01  0.8484E+01 "alien"
   45    6   15    1.42   11.97  0.8009E+01  0.9139E+01 "alien"
   36   24    2    1.78   11.71  0.9281E+01  0.8081E+01 "alien"
   39   34   12    1.36   11.64  0.1617E+02  0.1109E+02 "alien"
   45    5   14    1.46   11.55  0.7727E+01  0.8329E+01 "alien"
   26   21   14    1.78   11.45  0.9071E+01  0.6743E+01 "alien"
   47    1    6    1.60   11.44 -0.9207E+01  0.1147E+02 "alien"
   28   23   21    1.40   11.30 -0.1724E+02  0.7459E+01 "alien"
   35   30   13    1.46   11.19 -0.3902E+01  0.8529E+01 "alien"
   37   26    5    1.67   11.15 -0.7143E+01  0.7942E+01 "alien"
   36   27    8    1.62   11.13 -0.8959E+01  0.8559E+01 "alien"
   32   25   12    1.64   11.09 -0.7103E+01  0.7696E+01 "alien"
   30   25   14    1.61   10.88 -0.8756E+01  0.7422E+01 "alien"
   41   11   21    1.30   10.87  0.2294E+02  0.1328E+02 "alien"
   28   23   14    1.69   10.64 -0.6815E+01  0.6979E+01 "alien"
   34   22    2    1.90   10.60 -0.2887E+01  0.8251E+01 "alien"
   34    7   19    1.51   10.40 -0.3627E+01  0.1080E+02 "alien"
   37   26    3    1.70   10.36 -0.6637E+01  0.8870E+01 "alien"
   50    8   15    1.32   10.28  0.2169E+02  0.1639E+02 "alien"
   40   32    7    1.46   10.26  0.6868E+01  0.1037E+02 "alien"
   29   29   18    1.42   10.13  0.1379E+02  0.1461E+02 "alien"
   19   16   19    1.71   10.01  0.7933E+01  0.5564E+01 "alien"


 List of     91 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   41   11   21    1.30   10.87  0.2294E+02  0.1328E+02
   50    8   15    1.32   10.28  0.2169E+02  0.1639E+02
   54    5    9    1.36   61.00  0.8476E+02  0.9733E+01
   39   34   12    1.36   11.64  0.1617E+02  0.1109E+02
   37   10   21    1.36   13.11  0.1822E+02  0.1212E+02
   51    6   12    1.37   22.92  0.3185E+02  0.1684E+02
   27   13   24    1.37   14.46  0.2009E+02  0.1005E+02
   28   23   21    1.40   11.30 -0.1724E+02  0.7459E+01
   52    4    9    1.40   14.28 -0.2179E+02  0.9031E+01
   48    5   13    1.41   16.30  0.1091E+02  0.9475E+01
   52    4    8    1.42   18.07  0.1210E+02  0.9730E+01
   29   29   18    1.42   10.13  0.1379E+02  0.1461E+02
   36    9   20    1.42   17.98  0.1203E+02  0.1321E+02
   37   32   12    1.42   32.10  0.2148E+02  0.1021E+02
   45    6   15    1.42   11.97  0.8009E+01  0.9139E+01
   38    8   19    1.43   15.98  0.1069E+02  0.7762E+01
   54    3    1    1.43   27.22  0.1822E+02  0.9338E+01
   53    3    5    1.43   22.25  0.1489E+02  0.8988E+01
   52    3    7    1.44   26.28  0.1759E+02  0.9742E+01
   28   26   19    1.44   15.64  0.1047E+02  0.1179E+02
   40   32    8    1.44   13.10  0.8768E+01  0.8719E+01
   41   32    6    1.45   21.64  0.1448E+02  0.9010E+01
   45    5   14    1.46   11.55  0.7727E+01  0.8329E+01
   40   32    7    1.46   10.26  0.6868E+01  0.1037E+02
   53    2    1    1.46   18.41 -0.6418E+01  0.1072E+02
   37   31   11    1.46   31.13 -0.1085E+02  0.1089E+02
   52    3    5    1.46   37.60 -0.1311E+02  0.1590E+02
   39   31    9    1.46   35.44 -0.1236E+02  0.9765E+01
   41   32    5    1.46   25.99 -0.9061E+01  0.8913E+01
   38   31   10    1.46   23.76 -0.8284E+01  0.8674E+01
   47    4   12    1.46   16.00 -0.5580E+01  0.8581E+01
   35   30   13    1.46   11.19 -0.3902E+01  0.8529E+01
   36   30   12    1.47   22.37 -0.7798E+01  0.8560E+01
   40   31    7    1.48   17.48 -0.6094E+01  0.9190E+01
   37   30   11    1.48   37.52 -0.1308E+02  0.1014E+02
   33    8   20    1.48   25.75 -0.8977E+01  0.1114E+02
   49    3    9    1.48   11.98 -0.4176E+01  0.8484E+01
   52    2    1    1.49   17.24 -0.6013E+01  0.8530E+01
   52    2    0    1.49   14.39 -0.1359E+02  0.1242E+02
   30   28   16    1.49   33.20 -0.1158E+02  0.1230E+02
   41   31    4    1.49   23.39 -0.8157E+01  0.9001E+01
   40   31    6    1.49   25.06 -0.8740E+01  0.8573E+01
   32    8   20    1.50   14.04 -0.4894E+01  0.7302E+01
   29    9   21    1.50   34.46 -0.1201E+02  0.6206E+01
   41   31    3    1.50   35.56 -0.1240E+02  0.8504E+01
   34   29   13    1.50   35.55 -0.1240E+02  0.8915E+01
   38   30    9    1.50   39.39 -0.1373E+02  0.8497E+01
   50    2    6    1.51   20.76 -0.7240E+01  0.8522E+01
   34    7   19    1.51   10.40 -0.3627E+01  0.1080E+02
   41   31    1    1.51   12.08 -0.4212E+01  0.8384E+01
   45    2   10    1.57   12.87 -0.1035E+02  0.8798E+01
   39   29    4    1.57   18.94 -0.1524E+02  0.9579E+01
   25   24   17    1.60   16.71 -0.1345E+02  0.1050E+02
   47    1    6    1.60   11.44 -0.9207E+01  0.1147E+02
   27   24   16    1.61   13.03 -0.1049E+02  0.1314E+02
   30   25   14    1.61   10.88 -0.8756E+01  0.7422E+01
   36   27    8    1.62   11.13 -0.8959E+01  0.8559E+01
   32   25   12    1.64   11.09 -0.7103E+01  0.7696E+01
   29   24   14    1.65   12.84 -0.8221E+01  0.7310E+01
   34    4   16    1.65   13.55 -0.8681E+01  0.7067E+01
   18   10   21    1.66   14.47 -0.9268E+01  0.6760E+01
   37   26    5    1.67   11.15 -0.7143E+01  0.7942E+01
   43    1    8    1.69   19.32 -0.1238E+02  0.1412E+02
   28   23   14    1.69   10.64 -0.6815E+01  0.6979E+01
   42    1    9    1.70   12.32 -0.7888E+01  0.7535E+01
   36    3   14    1.70   12.32 -0.7893E+01  0.1198E+02
   37   26    3    1.70   10.36 -0.6637E+01  0.8870E+01
   19   16   19    1.71   10.01  0.7933E+01  0.5564E+01
   35    2   13    1.78   21.29  0.1687E+02  0.7775E+01
   26   21   14    1.78   11.45  0.9071E+01  0.6743E+01
   36   24    2    1.78   11.71  0.9281E+01  0.8081E+01
   33   23    8    1.79   18.05 -0.1923E+02  0.7558E+01
   43    1    2    1.79   13.44 -0.1433E+02  0.7949E+01
   31   22   10    1.80   17.24 -0.1838E+02  0.1247E+02
   29   21   11    1.84   12.34 -0.1315E+02  0.7278E+01
   24   19   14    1.87   26.33 -0.2806E+02  0.1141E+02
   33    1   12    1.90   56.84 -0.1548E+02  0.7743E+01
   21   18   15    1.90   19.37 -0.5275E+01  0.7139E+01
   15   12   18    1.90   27.38 -0.7457E+01  0.4860E+01
   34   22    2    1.90   10.60 -0.2887E+01  0.8251E+01
   34   22    1    1.91   14.33 -0.3902E+01  0.8497E+01
   34   22    0    1.91   15.35 -0.1215E+02  0.1194E+02
   29    2   14    1.92   26.19 -0.7132E+01  0.7289E+01
   26   19   12    1.93   36.70 -0.9995E+01  0.7272E+01
   40    2    1    1.93   12.35 -0.3365E+01  0.9909E+01
   30    1   13    1.95   92.23 -0.2512E+02  0.1159E+02
   28    2   14    1.95   70.99 -0.1933E+02  0.1170E+02
   39    2    3    1.96   60.87 -0.1658E+02  0.1006E+02
   30   20    8    1.96   23.89 -0.6505E+01  0.7827E+01
   34    0   10    1.96   19.08 -0.2011E+02  0.9108E+01
   29    1   13    1.99   22.94 -0.6248E+01  0.7301E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        34
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1337
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1238
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       99


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
