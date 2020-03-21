
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.550    77.550    38.810  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      32      32
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      32
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000954  0.000671
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
 BEAM_DIVERGENCE_E.S.D.=   0.047
 REFLECTING_RANGE_E.S.D.=   0.149
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

    1283 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      32 ...      32
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1283 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               41
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

 *  44        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  31        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0  1  0  0  1  0  0  0
 *  34        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  35        mP          0.0      77.6   38.8   77.6  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  32        oP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.9   38.8   77.6  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.9   38.8   77.6  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    28        mC        250.0      38.8  159.9   77.6  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.8  159.9   77.6  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    36        oC        250.0      38.8  159.9   77.6  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.8  159.9   77.6  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    27        mC        500.0     159.9   38.8  109.7  90.0 133.3  76.0    0  2  1  0  0  0  1  0  1 -1  0  0
    19        oI        500.0      38.8  109.7  116.3  90.0  70.5  90.0    0  0  1  0  1 -1  0  0  1  1  1  0
    26        oF        625.0      38.8  159.9  159.9  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.7  116.3   38.8  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        998.3     173.4   77.6   38.8  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.3      77.6  173.4   38.8  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.3      77.6  173.4   38.8  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.3      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.3  116.3  116.3  96.4  96.4 141.0    1 -1 -1  0 -1  1 -1  0  1  1  1  0
     2        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.7   86.7  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.7   86.7   86.7  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.7   86.7  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.7   86.7  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.8   86.7  248.3 102.1  99.0 116.6    0  0 -1  0  0  1  1  0  3 -1  1  0
    10        mC        999.0      86.7   86.7   77.6  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.7   86.7   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.8   77.6  177.7  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.7   86.7  177.7 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.7   86.7   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.7   86.7   86.7  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.7  109.7   38.8  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.8   77.6  177.7 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.7  177.7   77.6 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    4.2      17    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.905   -23.498     0.065
 UNIT_CELL_B-AXIS=    23.498    73.905     0.023
 UNIT_CELL_C-AXIS=    -0.034    -0.001    38.809
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       169 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.29
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.458    77.458    38.764  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.9E-01 1.9E-01 8.4E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012910  0.012910  0.025797  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.823   -23.449     0.033
 COORDINATES OF UNIT CELL B-AXIS    23.449    73.823     0.009
 COORDINATES OF UNIT CELL C-AXIS    -0.017     0.001    38.764
 CRYSTAL MOSAICITY (DEGREES)     0.149
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006473  0.001904
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.89   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1261.12   1242.29
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.08
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.46    77.46    38.76  90.000  90.000  90.000


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
  0.05  0.15        44         97.657       1.205
  0.15  0.25       171        117.496       1.002
  0.25  0.35       177        107.054       1.099
  0.35  0.45       163         73.857       1.593
  0.45  0.55       161        211.136       0.557
  0.55  0.65       185         87.727       1.341
  0.65  0.75       190        130.811       0.900
  0.75  0.85       151        113.965       1.033
  0.85  0.95        41         81.077       1.451
 
 -0.05  0.95      1283        117.682       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        7
 DEGREES OF FREEDOM OF CHI^2 FIT               184.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.886
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    31.1 XMAX=    31.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59443 NYBIN=    7
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        396


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               184.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.878
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    12.2 XMAX=  2453.3 NXBIN=    3
 YMIN=     8.8 YMAX=  2519.1 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        396


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               184.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.872
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    31.1 XMAX=    31.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1264
 DETECTOR_SURFACE_POSITION=    1647    1690
 DETECTOR_SURFACE_POSITION=     818    1690
 DETECTOR_SURFACE_POSITION=     818     838
 DETECTOR_SURFACE_POSITION=    1647     838
 DETECTOR_SURFACE_POSITION=    2171    1664
 DETECTOR_SURFACE_POSITION=    1621    2229
 DETECTOR_SURFACE_POSITION=     844    2229
 DETECTOR_SURFACE_POSITION=     294    1664
 DETECTOR_SURFACE_POSITION=     294     864
 DETECTOR_SURFACE_POSITION=     844     299
 DETECTOR_SURFACE_POSITION=    1621     299
 DETECTOR_SURFACE_POSITION=    2171     864
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        396


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
 7.812E-01  1.274E-03   31.69


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

   24.494   3.155    22.29   1.47      2.19      1.99      66      91      18
    3.155   2.240     7.34   1.16      7.17      6.30     170     202      10
    2.240   1.832     0.78   0.75     55.17     63.32     184     213       4
    1.832   1.587    -0.01   0.89   2127.85   2206.09     230     264       6
    1.587   1.420    -0.13   0.81    -99.90    -99.90     238     300       4
    1.420   1.297     0.00   1.13    -99.90    -99.90     114     166       4
   --------------------------------------------------------------------------
   24.494   1.297     2.94   0.96      8.02      7.95    1002    1236      46
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      32      32)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.684  -0.8338E+01  0.1243E+02        -0.67           1 
   47    0    0     1.648   0.7279E+01  0.1197E+02         0.61           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     5.9% OF MEAN INTENSITY


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

     3.85          58      36      1285        2.8%       2.0%      1.9%       44   30.01      2.8%    99.9*     0    0.000       0
     2.74          84      51      2098        2.4%       4.8%      3.4%       66   18.68      6.8%    99.3*     0    0.000       0
     2.24         151      88      2665        3.3%       8.6%      8.4%      126    8.04     12.1%    98.4*     0    0.000       0
     1.94         138      78      3101        2.5%      37.4%     42.9%      120    1.55     52.9%    85.0*     0    0.000       0
     1.74         153      86      3479        2.5%    1644.0%   1796.3%      134    0.01   2325.0%   -20.0      0    0.000       0
     1.59         186     106      3850        2.8%    1712.9%   1769.9%      160   -0.07   2422.4%    -2.7      0    0.000       0
     1.47         191     118      4143        2.8%     -99.9%    -99.9%      146    0.00    -99.9%    16.6      0    0.000       0
     1.38         192     114      4428        2.6%     -99.9%    -99.9%      156    0.00    -99.9%    -2.5      0    0.000       0
     1.30          83      58      4744        1.2%     -99.9%    -99.9%       50    0.00    -99.9%   -11.2      0    0.000       0
    total        1236     735     29793        2.5%       8.0%      8.0%     1002    3.84     11.3%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1283
 NUMBER OF REJECTED MISFITS                                46
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1236
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    735



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

   WILSON LINE (using all data) : A=   7.236 B=  32.760 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      15     7.258   0.005  1.8564E+03       7.526     -30.7
      17     4.726   0.011  1.8421E+03       7.519     -12.7
      23     3.536   0.020  7.0976E+02       6.565      16.8
      26     2.993   0.028  2.0186E+02       5.308      34.5
      31     2.645   0.036  1.7506E+02       5.165      29.0
      50     2.400   0.043  7.2371E+01       4.282      34.0
      35     2.209   0.051  3.4003E+01       3.526      36.2
      38     2.057   0.059  1.1935E+01       2.479      40.2
      39     1.935   0.067 -9.8650E-02    -999.999   -1000.0
      46     1.828   0.075 -3.6537E-01    -999.999   -1000.0
      44     1.739   0.083  3.1760E-01      -1.147      50.7
      44     1.668   0.090 -4.1993E-02    -999.999   -1000.0
      49     1.598   0.098 -2.0986E+00    -999.999   -1000.0
      58     1.538   0.106 -2.8701E+00    -999.999   -1000.0
      67     1.482   0.114 -2.8145E-01    -999.999   -1000.0
      50     1.435   0.121 -1.5657E+00    -999.999   -1000.0
      46     1.392   0.129 -9.4409E-01    -999.999   -1000.0
      40     1.352   0.137  2.7722E-02      -3.586      39.6
      17     1.316   0.144  4.5804E+00       1.522      19.8


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.258        2.070        1.544        0.900
       3     4.726        1.387        0.889        0.423
       2     3.536        0.383        0.083        0.013
       3     2.993        2.022        1.046        0.399
       2     2.645        0.364        0.107        0.023
       5     2.400        0.665        0.297        0.103
       1     2.209        1.491        0.631        0.191
       6     2.057        0.751        0.046        0.086
       6     1.935       25.600       -8.866       89.002
      11     1.828       14.787        3.155       46.550
      11     1.668       18.789       28.378       66.319
       3     1.598      114.088      553.422     2570.288
       3     1.538       13.276       10.363       19.975
       5     1.482       43.724       98.998      397.085
       3     1.435      105.048      368.832     1704.346
       1     1.392       38.168       81.686      124.870
      77   overall      282.179   -13022.352   454915.378


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.258        0.546        0.416        0.252
      14     4.726        1.289        1.646        1.665
      21     3.536        0.951        0.882        0.750
      23     2.993        0.743        0.573        0.363
      29     2.645        1.168        1.356        1.462
      45     2.400        0.962        0.834        0.613
      34     2.209        0.842        0.635        0.424
      32     2.057        1.634        2.057        2.418
      37     1.739      287.828      859.811    31244.653
      39     1.352    78607.869  4474524.093*************
      17     1.316        3.789        6.296        9.889
     658   overall     4675.745   265256.206218052961.573

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.258  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.73
       3     4.726  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.536  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.73
       3     2.993  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.645  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.400  0.00  0.00  0.00  0.42  0.77  1.07  1.00  0.95  0.91  0.88
       1     2.209  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.057  1.34  0.96  0.80  0.70  0.96  0.89  0.84  0.80  1.01  0.98
       6     1.935  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      11     1.828  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.72  0.83  0.80
      11     1.668  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.72  0.69  0.67
       3     1.598  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.538  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       5     1.482  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.59
       3     1.435  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       1     1.392  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      77   overall  1.36  0.98  0.81  0.74  0.72  0.69  0.65  0.62  0.65  0.69


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.258  2.87  3.01  2.10  1.93  1.85  1.61  1.63  1.49  1.38  1.29
      14     4.726  3.00  1.97  1.65  1.30  1.27  1.58  1.42  1.30  1.44  1.36
      21     3.536  1.00  0.79  0.73  0.72  0.73  0.95  1.23  1.21  1.12  1.05
      23     2.993  1.83  1.92  1.51  1.32  1.44  1.54  1.38  1.26  1.17  1.10
      29     2.645  0.72  0.95  1.46  1.15  1.14  0.99  1.10  1.00  1.05  1.09
      45     2.400  0.93  0.74  0.69  1.21  1.30  1.13  1.06  1.05  1.01  1.02
      34     2.209  0.93  0.81  0.79  0.98  0.90  0.85  0.99  0.91  0.89  0.93
      32     2.057  1.31  1.21  0.96  0.95  0.87  0.97  0.87  0.91  0.95  0.94
      37     1.739  5.68  2.98  2.09  1.64  1.37  1.20  1.07  0.98  0.91  0.86
      39     1.352  5.12  2.69  1.88  1.48  1.24  1.08  0.97  0.88  0.82  0.77
      17     1.316  5.56  2.92  2.04  1.61  1.35  1.30  1.17  1.18  1.09  1.02
     658   overall  1.18  0.78  0.63  0.58  0.54  0.52  0.51  0.48  0.47  0.45


 List of    126 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   28   17   24    1.33 1143.67  0.3170E+02  0.1234E+02 "alien"
   41    9   19    1.37  683.05  0.1894E+02  0.1483E+02 "alien"
   36   34   15    1.34  681.37  0.1889E+02  0.1358E+02 "alien"
   29   22   22    1.36  543.10  0.1506E+02  0.8337E+01 "alien"
   32   12   23    1.35  489.80  0.1358E+02  0.1148E+02 "alien"
   53    5   10    1.36  441.59  0.1224E+02  0.9670E+01 "alien"
   26   23   16    1.64  424.24 -0.1782E+02  0.1224E+02 "alien"
   49    7   15    1.34  411.28  0.1140E+02  0.1531E+02 "alien"
   37    3   14    1.67  366.29 -0.1538E+02  0.7443E+01 "alien"
   37   26    3    1.70  314.29 -0.1320E+02  0.9829E+01 "alien"
   21   19   19    1.64  262.06 -0.1100E+02  0.9449E+01 "alien"
   30   25   21    1.35  236.78  0.6564E+01  0.1128E+02 "alien"
   43    9   19    1.33  232.53  0.6446E+01  0.1504E+02 "alien"
   45    8   17    1.36  222.51  0.6168E+01  0.1446E+02 "alien"
   27   13   24    1.37  214.26  0.5940E+01  0.1060E+02 "alien"
   32   30   19    1.34  209.14  0.5798E+01  0.1308E+02 "alien"
   39    2    3    1.96  157.97 -0.1558E+02  0.1056E+02 "alien"
   50    7   14    1.34  156.54  0.4339E+01  0.1054E+02 "alien"
   39    2    2    1.97  155.95 -0.1538E+02  0.8414E+01 "alien"
   44    1    8    1.65  149.94 -0.6297E+01  0.8189E+01 "alien"
   30   20    8    1.96  147.91 -0.1459E+02  0.7968E+01 "alien"
   42    2   11    1.63  131.96 -0.5542E+01  0.8012E+01 "alien"
   52    5   11    1.37  120.34  0.3336E+01  0.1600E+02 "alien"
   29    1   13    1.99  117.42 -0.1158E+02  0.7769E+01 "alien"
   36    3   14    1.70   94.61 -0.3973E+01  0.1215E+02 "alien"
   33   12   23    1.34   89.15  0.2472E+01  0.1387E+02 "alien"
   30   28   16    1.49   86.99 -0.2449E+02  0.1421E+02 "alien"
   29   20   10    1.91   83.64 -0.8251E+01  0.1279E+02 "alien"
   27   19   11    1.95   79.75 -0.7867E+01  0.8133E+01 "alien"
   31   22    9    1.84   78.86 -0.2881E+02  0.1425E+02 "alien"
   28    2   14    1.95   69.58 -0.6865E+01  0.1254E+02 "alien"
   42    1    9    1.69   69.46 -0.2917E+01  0.7945E+01 "alien"
   53    2    1    1.46   69.08 -0.1944E+02  0.1126E+02 "alien"
   43    1    9    1.66   65.79 -0.2763E+01  0.8202E+01 "alien"
   31   24   12    1.69   63.03 -0.2647E+01  0.7566E+01 "alien"
   51    6   13    1.35   62.33  0.1728E+01  0.1593E+02 "alien"
   31   21    8    1.90   61.91 -0.6107E+01  0.7870E+01 "alien"
   34   22    1    1.91   58.50 -0.5771E+01  0.8494E+01 "alien"
   33   25   11    1.65   56.66 -0.2379E+01  0.1262E+02 "alien"
   40   31    6    1.49   51.70 -0.1455E+02  0.8812E+01 "alien"
   34   25    9    1.69   48.93 -0.2055E+01  0.7865E+01 "alien"
   38    6   17    1.51   48.78 -0.1373E+02  0.8845E+01 "alien"
   39   31    8    1.48   44.57 -0.1254E+02  0.9016E+01 "alien"
   34    3   15    1.71   43.93  0.1395E+02  0.1291E+02 "alien"
   32   23   10    1.75   42.83  0.1360E+02  0.1400E+02 "alien"
   35   30   13    1.46   42.26 -0.1189E+02  0.8717E+01 "alien"
   31   29   19    1.36   41.52  0.1151E+01  0.8612E+01 "alien"
   41   32    4    1.47   41.12 -0.1157E+02  0.1705E+02 "alien"
   27   25   19    1.47   40.92 -0.1152E+02  0.1291E+02 "alien"
   52    2    3    1.48   39.95 -0.1124E+02  0.9112E+01 "alien"
   19   14   20    1.67   39.74 -0.1669E+01  0.5649E+01 "alien"
   30   23   12    1.73   38.56  0.1225E+02  0.7717E+01 "alien"
   26   23   20    1.46   37.96 -0.1068E+02  0.1171E+02 "alien"
   44    1    1    1.76   37.81  0.1201E+02  0.9720E+01 "alien"
   38    1    7    1.91   35.33 -0.3486E+01  0.8133E+01 "alien"
   37   25    2    1.73   35.19  0.1117E+02  0.9170E+01 "alien"
   34   24    8    1.74   32.95  0.1047E+02  0.8109E+01 "alien"
   28   27   17    1.50   29.89 -0.8413E+01  0.8701E+01 "alien"
   35    2   13    1.78   29.12  0.9249E+01  0.8108E+01 "alien"
   51    3    7    1.46   29.06 -0.8180E+01  0.8762E+01 "alien"
   24   19   14    1.87   28.94 -0.1057E+02  0.1076E+02 "alien"
   20   19   18    1.71   28.60  0.9082E+01  0.9121E+01 "alien"
   48    6   14    1.39   27.58 -0.2604E+02  0.1274E+02 "alien"
   47    0    1    1.65   26.61 -0.8216E+01  0.8537E+01 "alien"
   35   24    6    1.76   26.47  0.8408E+01  0.8349E+01 "alien"
   39    1   11    1.73   25.92  0.8232E+01  0.1361E+02 "alien"
   43    1    2    1.79   25.45 -0.9298E+01  0.8011E+01 "alien"
   51    2    5    1.49   25.06 -0.7053E+01  0.9175E+01 "alien"
   34   33   16    1.36   24.69  0.6846E+00  0.1010E+02 "alien"
   37    1   11    1.80   24.25 -0.8859E+01  0.7820E+01 "alien"
   42   32    2    1.46   24.21 -0.6813E+01  0.9265E+01 "alien"
   52    2    2    1.48   24.08 -0.6777E+01  0.8930E+01 "alien"
   39    2    0    1.98   24.07 -0.1394E+02  0.1168E+02 "alien"
   30    9   21    1.48   23.66 -0.6658E+01  0.7314E+01 "alien"
   46    0    4    1.66   23.65 -0.5627E+01  0.9227E+01 "alien"
   37   27    6    1.64   23.59 -0.9908E+00  0.1164E+02 "alien"
   33   29   14    1.49   23.57 -0.6633E+01  0.1399E+02 "alien"
   39   31    9    1.46   23.56 -0.6632E+01  0.1016E+02 "alien"
   34    1   12    1.86   23.33 -0.8524E+01  0.7307E+01 "alien"
   21   18   15    1.90   23.22 -0.2291E+01  0.9656E+01 "alien"
   35   23    1    1.85   23.12 -0.8449E+01  0.8130E+01 "alien"
   45    4   13    1.49   22.77 -0.6408E+01  0.9140E+01 "alien"
   31   30   18    1.38   22.45 -0.2119E+02  0.1334E+02 "alien"
   21   20   17    1.73   22.13  0.7028E+01  0.7270E+01 "alien"
   29   21   11    1.84   21.77 -0.7954E+01  0.7575E+01 "alien"
   35   32   15    1.38   21.40 -0.2020E+02  0.1514E+02 "alien"
   46    0    0    1.68   21.33 -0.8338E+01  0.1243E+02 "alien"
   39    1    6    1.90   21.22 -0.2094E+01  0.8030E+01 "alien"
   25   20   21    1.47   20.94 -0.5893E+01  0.1017E+02 "alien"
   25   10   22    1.50   20.88 -0.5875E+01  0.8838E+01 "alien"
   41   31    1    1.51   20.15 -0.5671E+01  0.9071E+01 "alien"
   44    0    6    1.70   19.79 -0.3407E+01  0.8482E+01 "alien"
   47    0    2    1.64   19.79 -0.3405E+01  0.8901E+01 "alien"
   25   19   13    1.90   19.16 -0.1890E+01  0.7446E+01 "alien"
   29   24   14    1.65   18.62 -0.7818E+00  0.7570E+01 "alien"
   32   21    6    1.93   18.24 -0.1800E+01  0.8103E+01 "alien"
   34   23    5    1.83   17.98 -0.6570E+01  0.9013E+01 "alien"
   35   33   15    1.37   17.69  0.4905E+00  0.1428E+02 "alien"
   36    2   13    1.74   17.67  0.5612E+01  0.7040E+01 "alien"
   34   24    7    1.76   17.09  0.5427E+01  0.9643E+01 "alien"
   38   26    0    1.68   16.86 -0.2304E+01  0.1459E+02 "alien"
   45    1    7    1.64   16.30 -0.6845E+00  0.8332E+01 "alien"
   34    0   10    1.96   16.29 -0.5008E+01  0.9618E+01 "alien"
   49    4   11    1.44   15.97 -0.2500E+02  0.9267E+01 "alien"
   35   23    2    1.84   15.93 -0.5820E+01  0.8100E+01 "alien"
   38    7   18    1.47   15.54 -0.4373E+01  0.8962E+01 "alien"
   52    3    6    1.45   15.39 -0.2409E+02  0.9050E+01 "alien"
   37    8   19    1.44   15.25 -0.2387E+02  0.1276E+02 "alien"
   53    3    6    1.42   15.18 -0.2376E+02  0.1116E+02 "alien"
   30   22   11    1.79   14.15 -0.5170E+01  0.7533E+01 "alien"
   41   31    0    1.51   13.98 -0.1035E+02  0.1294E+02 "alien"
   38   31   10    1.46   13.46 -0.3789E+01  0.8759E+01 "alien"
   42   32    0    1.47   12.89 -0.8856E+01  0.1246E+02 "alien"
   43    1    0    1.80   12.79 -0.1132E+02  0.1115E+02 "alien"
   33   24    9    1.74   12.61  0.4006E+01  0.8275E+01 "alien"
   52    5   10    1.38   12.33 -0.1164E+02  0.1028E+02 "alien"
   55    3    1    1.41   11.82 -0.1116E+02  0.8456E+01 "alien"
   51    4   10    1.41   11.75 -0.1110E+02  0.1557E+02 "alien"
   33   27   13    1.55   11.49 -0.3297E+02  0.1259E+02 "alien"
   24   13   24    1.40   11.18 -0.1056E+02  0.9823E+01 "alien"
   21   17   20    1.61   11.16 -0.2341E+02  0.8341E+01 "alien"
   36    0   10    1.88   11.10 -0.8780E+01  0.7882E+01 "alien"
   37    4   15    1.62   10.96 -0.2301E+02  0.1191E+02 "alien"
   44    0    5    1.72   10.93  0.7431E+01  0.8139E+01 "alien"
   31    2   14    1.85   10.41 -0.3805E+01  0.6936E+01 "alien"
   28   26   19    1.44   10.10 -0.1582E+02  0.1235E+02 "alien"


 List of    126 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   28   17   24    1.33 1143.67  0.3170E+02  0.1234E+02
   43    9   19    1.33  232.53  0.6446E+01  0.1504E+02
   32   30   19    1.34  209.14  0.5798E+01  0.1308E+02
   36   34   15    1.34  681.37  0.1889E+02  0.1358E+02
   49    7   15    1.34  411.28  0.1140E+02  0.1531E+02
   33   12   23    1.34   89.15  0.2472E+01  0.1387E+02
   50    7   14    1.34  156.54  0.4339E+01  0.1054E+02
   51    6   13    1.35   62.33  0.1728E+01  0.1593E+02
   30   25   21    1.35  236.78  0.6564E+01  0.1128E+02
   32   12   23    1.35  489.80  0.1358E+02  0.1148E+02
   34   33   16    1.36   24.69  0.6846E+00  0.1010E+02
   29   22   22    1.36  543.10  0.1506E+02  0.8337E+01
   31   29   19    1.36   41.52  0.1151E+01  0.8612E+01
   45    8   17    1.36  222.51  0.6168E+01  0.1446E+02
   53    5   10    1.36  441.59  0.1224E+02  0.9670E+01
   35   33   15    1.37   17.69  0.4905E+00  0.1428E+02
   52    5   11    1.37  120.34  0.3336E+01  0.1600E+02
   41    9   19    1.37  683.05  0.1894E+02  0.1483E+02
   27   13   24    1.37  214.26  0.5940E+01  0.1060E+02
   31   30   18    1.38   22.45 -0.2119E+02  0.1334E+02
   35   32   15    1.38   21.40 -0.2020E+02  0.1514E+02
   52    5   10    1.38   12.33 -0.1164E+02  0.1028E+02
   48    6   14    1.39   27.58 -0.2604E+02  0.1274E+02
   24   13   24    1.40   11.18 -0.1056E+02  0.9823E+01
   55    3    1    1.41   11.82 -0.1116E+02  0.8456E+01
   51    4   10    1.41   11.75 -0.1110E+02  0.1557E+02
   53    3    6    1.42   15.18 -0.2376E+02  0.1116E+02
   28   26   19    1.44   10.10 -0.1582E+02  0.1235E+02
   49    4   11    1.44   15.97 -0.2500E+02  0.9267E+01
   37    8   19    1.44   15.25 -0.2387E+02  0.1276E+02
   52    3    6    1.45   15.39 -0.2409E+02  0.9050E+01
   53    2    1    1.46   69.08 -0.1944E+02  0.1126E+02
   51    3    7    1.46   29.06 -0.8180E+01  0.8762E+01
   39   31    9    1.46   23.56 -0.6632E+01  0.1016E+02
   38   31   10    1.46   13.46 -0.3789E+01  0.8759E+01
   42   32    2    1.46   24.21 -0.6813E+01  0.9265E+01
   26   23   20    1.46   37.96 -0.1068E+02  0.1171E+02
   35   30   13    1.46   42.26 -0.1189E+02  0.8717E+01
   27   25   19    1.47   40.92 -0.1152E+02  0.1291E+02
   42   32    0    1.47   12.89 -0.8856E+01  0.1246E+02
   38    7   18    1.47   15.54 -0.4373E+01  0.8962E+01
   25   20   21    1.47   20.94 -0.5893E+01  0.1017E+02
   41   32    4    1.47   41.12 -0.1157E+02  0.1705E+02
   52    2    3    1.48   39.95 -0.1124E+02  0.9112E+01
   30    9   21    1.48   23.66 -0.6658E+01  0.7314E+01
   39   31    8    1.48   44.57 -0.1254E+02  0.9016E+01
   52    2    2    1.48   24.08 -0.6777E+01  0.8930E+01
   45    4   13    1.49   22.77 -0.6408E+01  0.9140E+01
   33   29   14    1.49   23.57 -0.6633E+01  0.1399E+02
   30   28   16    1.49   86.99 -0.2449E+02  0.1421E+02
   51    2    5    1.49   25.06 -0.7053E+01  0.9175E+01
   40   31    6    1.49   51.70 -0.1455E+02  0.8812E+01
   28   27   17    1.50   29.89 -0.8413E+01  0.8701E+01
   25   10   22    1.50   20.88 -0.5875E+01  0.8838E+01
   41   31    1    1.51   20.15 -0.5671E+01  0.9071E+01
   41   31    0    1.51   13.98 -0.1035E+02  0.1294E+02
   38    6   17    1.51   48.78 -0.1373E+02  0.8845E+01
   33   27   13    1.55   11.49 -0.3297E+02  0.1259E+02
   21   17   20    1.61   11.16 -0.2341E+02  0.8341E+01
   37    4   15    1.62   10.96 -0.2301E+02  0.1191E+02
   42    2   11    1.63  131.96 -0.5542E+01  0.8012E+01
   21   19   19    1.64  262.06 -0.1100E+02  0.9449E+01
   37   27    6    1.64   23.59 -0.9908E+00  0.1164E+02
   26   23   16    1.64  424.24 -0.1782E+02  0.1224E+02
   47    0    2    1.64   19.79 -0.3405E+01  0.8901E+01
   45    1    7    1.64   16.30 -0.6845E+00  0.8332E+01
   47    0    1    1.65   26.61 -0.8216E+01  0.8537E+01
   29   24   14    1.65   18.62 -0.7818E+00  0.7570E+01
   33   25   11    1.65   56.66 -0.2379E+01  0.1262E+02
   44    1    8    1.65  149.94 -0.6297E+01  0.8189E+01
   46    0    4    1.66   23.65 -0.5627E+01  0.9227E+01
   43    1    9    1.66   65.79 -0.2763E+01  0.8202E+01
   37    3   14    1.67  366.29 -0.1538E+02  0.7443E+01
   19   14   20    1.67   39.74 -0.1669E+01  0.5649E+01
   38   26    0    1.68   16.86 -0.2304E+01  0.1459E+02
   46    0    0    1.68   21.33 -0.8338E+01  0.1243E+02
   31   24   12    1.69   63.03 -0.2647E+01  0.7566E+01
   34   25    9    1.69   48.93 -0.2055E+01  0.7865E+01
   42    1    9    1.69   69.46 -0.2917E+01  0.7945E+01
   36    3   14    1.70   94.61 -0.3973E+01  0.1215E+02
   37   26    3    1.70  314.29 -0.1320E+02  0.9829E+01
   44    0    6    1.70   19.79 -0.3407E+01  0.8482E+01
   34    3   15    1.71   43.93  0.1395E+02  0.1291E+02
   20   19   18    1.71   28.60  0.9082E+01  0.9121E+01
   44    0    5    1.72   10.93  0.7431E+01  0.8139E+01
   37   25    2    1.73   35.19  0.1117E+02  0.9170E+01
   39    1   11    1.73   25.92  0.8232E+01  0.1361E+02
   30   23   12    1.73   38.56  0.1225E+02  0.7717E+01
   21   20   17    1.73   22.13  0.7028E+01  0.7270E+01
   33   24    9    1.74   12.61  0.4006E+01  0.8275E+01
   34   24    8    1.74   32.95  0.1047E+02  0.8109E+01
   36    2   13    1.74   17.67  0.5612E+01  0.7040E+01
   32   23   10    1.75   42.83  0.1360E+02  0.1400E+02
   35   24    6    1.76   26.47  0.8408E+01  0.8349E+01
   44    1    1    1.76   37.81  0.1201E+02  0.9720E+01
   34   24    7    1.76   17.09  0.5427E+01  0.9643E+01
   35    2   13    1.78   29.12  0.9249E+01  0.8108E+01
   30   22   11    1.79   14.15 -0.5170E+01  0.7533E+01
   43    1    2    1.79   25.45 -0.9298E+01  0.8011E+01
   37    1   11    1.80   24.25 -0.8859E+01  0.7820E+01
   43    1    0    1.80   12.79 -0.1132E+02  0.1115E+02
   34   23    5    1.83   17.98 -0.6570E+01  0.9013E+01
   35   23    2    1.84   15.93 -0.5820E+01  0.8100E+01
   31   22    9    1.84   78.86 -0.2881E+02  0.1425E+02
   29   21   11    1.84   21.77 -0.7954E+01  0.7575E+01
   35   23    1    1.85   23.12 -0.8449E+01  0.8130E+01
   31    2   14    1.85   10.41 -0.3805E+01  0.6936E+01
   34    1   12    1.86   23.33 -0.8524E+01  0.7307E+01
   24   19   14    1.87   28.94 -0.1057E+02  0.1076E+02
   36    0   10    1.88   11.10 -0.8780E+01  0.7882E+01
   39    1    6    1.90   21.22 -0.2094E+01  0.8030E+01
   21   18   15    1.90   23.22 -0.2291E+01  0.9656E+01
   25   19   13    1.90   19.16 -0.1890E+01  0.7446E+01
   31   21    8    1.90   61.91 -0.6107E+01  0.7870E+01
   34   22    1    1.91   58.50 -0.5771E+01  0.8494E+01
   38    1    7    1.91   35.33 -0.3486E+01  0.8133E+01
   29   20   10    1.91   83.64 -0.8251E+01  0.1279E+02
   32   21    6    1.93   18.24 -0.1800E+01  0.8103E+01
   27   19   11    1.95   79.75 -0.7867E+01  0.8133E+01
   28    2   14    1.95   69.58 -0.6865E+01  0.1254E+02
   39    2    3    1.96  157.97 -0.1558E+02  0.1056E+02
   30   20    8    1.96  147.91 -0.1459E+02  0.7968E+01
   34    0   10    1.96   16.29 -0.5008E+01  0.9618E+01
   39    2    2    1.97  155.95 -0.1538E+02  0.8414E+01
   39    2    0    1.98   24.07 -0.1394E+02  0.1168E+02
   29    1   13    1.99  117.42 -0.1158E+02  0.7769E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        91
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1283
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1094
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      189


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
