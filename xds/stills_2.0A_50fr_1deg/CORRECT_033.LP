
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.550    77.550    38.807  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      33      33
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      33
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000970  0.000669
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
 REFLECTING_RANGE_E.S.D.=   0.150
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

    1281 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      33 ...      33
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1281 REFLECTIONS ACCEPTED



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

 *  31        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0  1  0  0  1  0  0  0
 *  44        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  35        mP          0.0      77.6   38.8   77.6  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  32        oP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  20        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  25        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    28        mC        250.0      38.8  159.9   77.6  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.8  159.9   77.6  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    37        mC        250.0     159.9   38.8   77.6  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.9   38.8   77.6  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.8  159.9   77.6  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.8  159.9   77.6  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    27        mC        500.0     159.9   38.8  109.7  90.0 133.3  76.0    0  2  1  0  0  0  1  0  1 -1  0  0
    19        oI        500.0      38.8  109.7  116.3  90.0  70.5  90.0    0  0  1  0  1 -1  0  0  1  1  1  0
    26        oF        625.0      38.8  159.9  159.9  86.6 104.0 104.0    0  0 -1  0  0  2  1  0  2  0  1  0
    18        tI        625.0     109.7  116.3   38.8  70.5  90.0  90.0    1 -1  0  0 -1 -1 -1  0  0  0 -1  0
    30        mC        998.4      77.6  173.4   38.8  90.0  90.0  63.4    0  1  0  0 -2  1  0  0  0  0  1  0
    40        oC        998.4      77.6  173.4   38.8  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    41        mC        998.4     173.4   77.6   38.8  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    22        hP        998.4      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.3  116.3  116.3  96.4  96.4 141.0    1 -1 -1  0 -1  1 -1  0  1  1  1  0
     2        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.7   86.7  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.7   86.7   86.7  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.7   86.7  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.7   86.7  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.8   86.7  248.3 102.1  99.0 116.6    0  0 -1  0  0  1  1  0  3 -1  1  0
    10        mC        999.0      86.7   86.7   77.6  90.0  90.0 126.8    0  1 -1  0  0 -1 -1  0 -1  0  0  0
    11        tP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.7   86.7   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.8   77.6  177.7  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.7   86.7  177.7 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.7   86.7   77.6  90.0  90.0  53.2    0 -1  1  0  0 -1 -1  0  1  0  0  0
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

  *   0.00        0    4.3      17    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.884   -23.565     0.094
 UNIT_CELL_B-AXIS=    23.565    73.884     0.034
 UNIT_CELL_C-AXIS=    -0.050    -0.002    38.807
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       167 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.27
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.11
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.520    77.520    38.793  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  7.2E-02 7.2E-02 2.8E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012900  0.012900  0.025778  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.870   -23.504     0.057
 COORDINATES OF UNIT CELL B-AXIS    23.504    73.870     0.021
 COORDINATES OF UNIT CELL C-AXIS    -0.030    -0.002    38.793
 CRYSTAL MOSAICITY (DEGREES)     0.150
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999979  0.005988  0.002311
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.89   1232.86
 DETECTOR ORIGIN (PIXELS) AT                     1260.53   1241.58
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.28
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.52    77.52    38.79  90.000  90.000  90.000


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
  0.05  0.15        47         96.593       1.202
  0.15  0.25       165        121.327       0.957
  0.25  0.35       179        110.565       1.050
  0.35  0.45       162         62.110       1.870
  0.45  0.55       163        204.335       0.568
  0.55  0.65       185         86.563       1.342
  0.65  0.75       196        123.953       0.937
  0.75  0.85       139        114.199       1.017
  0.85  0.95        45        108.435       1.071
 
 -0.05  0.95      1281        116.144       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        8
 DEGREES OF FREEDOM OF CHI^2 FIT               193.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.035
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    32.1 XMAX=    32.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59350 NYBIN=    8
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        408


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               193.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.019
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    13.2 XMAX=  2453.0 NXBIN=    3
 YMIN=     6.4 YMAX=  2520.2 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        408


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               191.1
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.007
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    32.1 XMAX=    32.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1263
 DETECTOR_SURFACE_POSITION=    1647    1690
 DETECTOR_SURFACE_POSITION=     819    1690
 DETECTOR_SURFACE_POSITION=     819     836
 DETECTOR_SURFACE_POSITION=    1647     836
 DETECTOR_SURFACE_POSITION=    2171    1664
 DETECTOR_SURFACE_POSITION=    1622    2229
 DETECTOR_SURFACE_POSITION=     845    2229
 DETECTOR_SURFACE_POSITION=     295    1664
 DETECTOR_SURFACE_POSITION=     295     863
 DETECTOR_SURFACE_POSITION=     845     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2171     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        408


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
 7.858E-01  8.677E-04   38.30


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

   24.514   3.157    25.04   1.44      1.82      1.71      64      90      20
    3.157   2.242     7.40   1.27      7.74      6.44     178     209       2
    2.242   1.833     0.69   1.23     76.62     66.82     188     215       2
    1.833   1.589    -0.05   0.73    -99.90    -99.90     236     269       0
    1.589   1.422    -0.10   0.78    -99.90    -99.90     244     302       2
    1.422   1.298    -0.07   0.86    -99.90    -99.90     112     165       4
   --------------------------------------------------------------------------
   24.514   1.298     3.12   0.99      8.21      8.15    1022    1250      30
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      33      33)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.685  -0.1664E+02  0.1198E+02        -1.39           1 
   47    0    0     1.649   0.6914E+01  0.1190E+02         0.58           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     5.7% OF MEAN INTENSITY


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

     3.85          57      36      1285        2.8%       1.6%      1.6%       42   34.39      2.3%    99.9*     0    0.000       0
     2.74          88      53      2098        2.5%       4.6%      3.4%       70   19.34      6.5%    99.4*     0    0.000       0
     2.24         154      89      2665        3.3%       9.8%      8.6%      130    8.01     13.9%    98.2*     0    0.000       0
     1.94         140      79      3101        2.5%      54.1%     46.4%      122    1.34     76.6%    65.9*     0    0.000       0
     1.74         154      85      3479        2.4%     910.3%    956.3%      138   -0.01   1287.4%     2.0      0    0.000       0
     1.59         188     107      3851        2.8%     -99.9%    -99.9%      162    0.00    -99.9%     7.8      0    0.000       0
     1.47         193     116      4143        2.8%     -99.9%    -99.9%      154    0.00    -99.9%    21.5      0    0.000       0
     1.38         187     113      4427        2.6%     -99.9%    -99.9%      148    0.00    -99.9%    19.6      0    0.000       0
     1.30          83      57      4739        1.2%    4928.5%   4573.3%       52    0.11   6970.0%    36.8      0    0.000       0
    total        1244     735     29788        2.5%       8.2%      8.1%     1018    4.15     11.6%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1281
 NUMBER OF REJECTED MISFITS                                30
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1250
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    739



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

   WILSON LINE (using all data) : A=   6.937 B=  27.233 CORRELATION=  0.89
      #      RES      SS        <I>       log(<I>)       BO
      14     7.402   0.005  1.9178E+03       7.559     -68.1
      18     4.711   0.011  1.7912E+03       7.491     -24.6
      23     3.549   0.020  6.7776E+02       6.519      10.5
      28     2.995   0.028  1.8819E+02       5.237      30.5
      32     2.642   0.036  1.6329E+02       5.096      25.7
      50     2.402   0.043  7.0282E+01       4.253      31.0
      35     2.211   0.051  2.9049E+01       3.369      34.9
      39     2.059   0.059  1.1022E+01       2.400      38.5
      39     1.936   0.067 -8.9583E-01    -999.999   -1000.0
      45     1.831   0.075 -2.4555E-01    -999.999   -1000.0
      44     1.740   0.083  2.3731E+00       0.864      36.8
      44     1.670   0.090 -9.3295E-01    -999.999   -1000.0
      51     1.600   0.098 -2.3911E+00    -999.999   -1000.0
      56     1.540   0.105 -2.6554E+00    -999.999   -1000.0
      66     1.484   0.114  6.8537E-01      -0.378      32.2
      52     1.436   0.121 -1.4569E+00    -999.999   -1000.0
      43     1.392   0.129 -8.6814E-02    -999.999   -1000.0
      38     1.354   0.136 -5.7280E-01    -999.999   -1000.0
      17     1.318   0.144  5.0120E+00       1.612      18.5


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.402        1.890        1.347        0.750
       3     4.711        1.407        0.911        0.438
       2     3.549        0.396        0.086        0.013
       3     2.995        2.000        1.017        0.378
       2     2.642        0.396        0.122        0.027
       5     2.402        0.938        0.575        0.272
       1     2.211        1.423        0.588        0.174
       6     2.059        0.409        0.132        0.052
       6     1.936       49.930      170.768      797.969
       7     1.740        2.996        3.027        2.424
       3     1.600       28.015       66.526      124.632
      75   overall       15.090      -83.358      648.526


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.402        0.550        0.408        0.238
      15     4.711        1.248        1.580        1.605
      21     3.549        0.981        0.939        0.821
      25     2.995        0.778        0.643        0.449
      30     2.642        1.180        1.373        1.483
      45     2.402        0.897        0.725        0.498
      34     2.211        0.901        0.764        0.573
      33     2.059        2.179        2.949        3.798
      35     1.831       36.039      -52.372      714.972
      37     1.740        6.556        6.838       16.817
      61     1.484      113.214     -103.100     6212.877
      17     1.318        2.305        3.365        4.830
     659   overall       11.275       -8.677      576.565

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.402  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.73
       3     4.711  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.549  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.995  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.642  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.402  0.00  0.00  0.96  0.85  0.77  0.71  0.67  0.64  0.91  0.88
       1     2.211  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.059  3.36  2.41  2.00  1.76  1.60  1.48  1.40  1.33  1.27  1.22
       6     1.936  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       7     1.740  2.30  1.65  1.37  1.51  1.37  1.27  1.20  1.14  1.09  1.05
       3     1.600  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      75   overall  0.86  0.62  0.58  0.54  0.49  0.45  0.42  0.40  0.43  0.45


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.402  3.15  3.31  2.31  2.12  1.78  1.55  1.59  1.45  1.35  1.27
      15     4.711  2.80  1.47  1.54  1.21  1.19  1.48  1.32  1.21  1.35  1.37
      21     3.549  1.00  1.05  0.92  0.87  0.85  0.95  1.23  1.21  1.12  1.05
      25     2.995  1.26  1.77  1.39  1.33  1.32  1.33  1.27  1.23  1.21  1.14
      30     2.642  1.05  1.29  1.16  1.11  1.02  1.03  0.99  1.09  1.07  1.11
      45     2.402  0.70  0.61  0.60  0.88  1.13  1.13  1.06  1.09  1.09  1.05
      34     2.211  1.24  1.14  1.13  0.98  0.82  0.85  0.82  0.85  0.84  0.98
      33     2.059  2.87  1.67  1.17  1.01  1.08  1.07  0.96  0.94  0.87  0.81
      35     1.831  6.91  3.63  2.54  1.99  1.67  1.46  1.31  1.19  1.11  1.04
      37     1.740  5.11  2.68  1.88  1.48  1.24  1.08  1.02  0.93  0.87  0.81
      61     1.484  4.31  2.26  1.71  1.34  1.12  0.98  0.88  0.80  0.77  0.75
      17     1.318  4.33  2.27  1.82  1.43  1.20  1.04  0.93  0.85  0.79  0.74
     659   overall  1.66  1.04  0.81  0.70  0.64  0.62  0.58  0.56  0.54  0.53


 List of     93 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   53    5    9    1.38  330.40 -0.2868E+02  0.1653E+02 "alien"
   36   33   14    1.38  239.67 -0.2081E+02  0.1086E+02 "alien"
   35   32   15    1.38  212.95 -0.1849E+02  0.1609E+02 "alien"
   51    4   10    1.41  204.76 -0.1778E+02  0.1556E+02 "alien"
   45    7   16    1.39  202.03 -0.1754E+02  0.1362E+02 "alien"
   25   13   24    1.39  182.01 -0.1580E+02  0.7281E+01 "alien"
   40   33    9    1.41  127.80 -0.1109E+02  0.9079E+01 "alien"
   34   31   15    1.41  123.61 -0.1073E+02  0.9022E+01 "alien"
   26   13   24    1.38  123.61 -0.1073E+02  0.9828E+01 "alien"
   48    6   14    1.39  109.01 -0.9463E+01  0.1275E+02 "alien"
   24   19   14    1.87   95.21 -0.2338E+02  0.1122E+02 "alien"
   38   33   12    1.39   94.03 -0.8163E+01  0.9508E+01 "alien"
   38    9   20    1.39   86.85 -0.7540E+01  0.1178E+02 "alien"
   30   28   19    1.39   65.47 -0.5683E+01  0.8560E+01 "alien"
   35    1   12    1.83   60.44 -0.1484E+02  0.7980E+01 "alien"
   17   15   18    1.82   59.82 -0.1469E+02  0.6083E+01 "alien"
   27   18   23    1.38   50.16 -0.4354E+01  0.7949E+01 "alien"
   51    5   11    1.39   44.86 -0.3895E+01  0.9318E+01 "alien"
   31    3   15    1.79   44.70 -0.1098E+02  0.7686E+01 "alien"
   46    7   16    1.37   44.55 -0.3868E+01  0.9091E+01 "alien"
   34   22    4    1.88   43.53 -0.1069E+02  0.9236E+01 "alien"
   39   34   12    1.36   42.69 -0.2445E+02  0.1120E+02 "alien"
   31    2   14    1.85   39.77 -0.9766E+01  0.6943E+01 "alien"
   53    4    8    1.40   35.90 -0.3116E+01  0.9503E+01 "alien"
   41    9   19    1.37   30.05 -0.1721E+02  0.1461E+02 "alien"
   43    5   15    1.47   29.13  0.1997E+02  0.9962E+01 "alien"
   32   30   19    1.34   28.20 -0.1615E+02  0.1358E+02 "alien"
   29   22   22    1.36   27.94 -0.1600E+02  0.8317E+01 "alien"
   50    7   14    1.34   27.41 -0.1570E+02  0.1048E+02 "alien"
   20   18   16    1.86   26.67 -0.6550E+01  0.9183E+01 "alien"
   46    8   17    1.34   25.92 -0.1485E+02  0.9179E+01 "alien"
   50    3    8    1.47   25.77  0.1766E+02  0.9105E+01 "alien"
   39   34   11    1.38   25.30 -0.2197E+01  0.1719E+02 "alien"
   49    6   14    1.37   24.73 -0.1417E+02  0.1689E+02 "alien"
   52    2    2    1.49   24.63  0.1688E+02  0.9066E+01 "alien"
   40   31    7    1.48   23.94  0.1641E+02  0.9636E+01 "alien"
   52    2    1    1.49   22.78  0.1561E+02  0.8850E+01 "alien"
   39    0    8    1.84   21.81 -0.2592E+02  0.1504E+02 "alien"
   42    1    4    1.81   21.64 -0.5314E+01  0.8309E+01 "alien"
   33    1   12    1.90   21.39 -0.1916E+02  0.8398E+01 "alien"
   52    3    5    1.46   21.22  0.1455E+02  0.1713E+02 "alien"
   51    6   13    1.35   20.37 -0.1167E+02  0.1618E+02 "alien"
   49    7   15    1.34   20.16 -0.1155E+02  0.1488E+02 "alien"
   28   21   12    1.83   19.98 -0.4907E+01  0.7359E+01 "alien"
   50    6   13    1.37   19.68 -0.1128E+02  0.9208E+01 "alien"
   44    5   14    1.48   19.62  0.1345E+02  0.1571E+02 "alien"
   37   34   14    1.35   19.38 -0.1110E+02  0.9403E+01 "alien"
   27   13   24    1.37   19.12 -0.1095E+02  0.1053E+02 "alien"
   48    5   13    1.41   18.19 -0.2649E+02  0.1467E+02 "alien"
   53    2    1    1.46   17.91  0.1227E+02  0.1134E+02 "alien"
   33   25   11    1.65   17.51 -0.1633E+02  0.1295E+02 "alien"
   46    5   14    1.43   16.84 -0.2453E+02  0.1450E+02 "alien"
   39    6   17    1.49   16.57  0.1136E+02  0.8009E+01 "alien"
   39    0    9    1.81   16.39 -0.1264E+02  0.8669E+01 "alien"
   32   29   15    1.47   15.62  0.1071E+02  0.1345E+02 "alien"
   29   24   14    1.65   15.57 -0.1453E+02  0.7581E+01 "alien"
   38    6   17    1.51   15.42  0.1057E+02  0.9210E+01 "alien"
   53    6   11    1.34   15.36 -0.8797E+01  0.1046E+02 "alien"
   21   17   14    1.99   15.24 -0.1365E+02  0.1106E+02 "alien"
   26   23   20    1.46   14.54  0.9967E+01  0.1240E+02 "alien"
   52    2    3    1.48   14.39  0.9862E+01  0.9274E+01 "alien"
   41   31    2    1.50   14.10  0.9663E+01  0.8967E+01 "alien"
   36   31   13    1.43   14.03 -0.2044E+02  0.9294E+01 "alien"
   27   25   19    1.47   13.96  0.9568E+01  0.1329E+02 "alien"
   43    8   18    1.37   13.63 -0.7806E+01  0.9641E+01 "alien"
   51    2    4    1.50   13.61  0.9329E+01  0.8849E+01 "alien"
   29   20   10    1.91   13.20 -0.1183E+02  0.1344E+02 "alien"
   31   21    8    1.90   13.17 -0.1180E+02  0.7987E+01 "alien"
   27   23   15    1.67   12.91 -0.1204E+02  0.7525E+01 "alien"
   33   33   17    1.34   12.82 -0.1782E+02  0.1328E+02 "alien"
   26   10   22    1.49   12.74  0.8729E+01  0.9840E+01 "alien"
   35   23    0    1.85   12.51 -0.7303E+01  0.1186E+02 "alien"
   30   20    8    1.97   12.48 -0.1118E+02  0.8189E+01 "alien"
   34    7   19    1.51   12.27  0.8407E+01  0.1162E+02 "alien"
   36   26    7    1.67   12.23 -0.1141E+02  0.7909E+01 "alien"
   43    1    8    1.69   12.10 -0.1129E+02  0.1555E+02 "alien"
   38    8   19    1.43   11.89 -0.1732E+02  0.7875E+01 "alien"
   28   23   14    1.69   11.87 -0.1108E+02  0.7322E+01 "alien"
   42   32    3    1.46   11.80 -0.1719E+02  0.1094E+02 "alien"
   41   32    5    1.46   11.55  0.7914E+01  0.9045E+01 "alien"
   51    2    3    1.51   11.48  0.7871E+01  0.9208E+01 "alien"
   25   19   13    1.90   11.46 -0.1027E+02  0.7646E+01 "alien"
   40    0    8    1.80   11.38 -0.6156E+01  0.7995E+01 "alien"
   35   29   12    1.51   11.30  0.7744E+01  0.9278E+01 "alien"
   45    6   15    1.42   11.29 -0.1644E+02  0.9503E+01 "alien"
   28   21   22    1.38   11.23 -0.9747E+00  0.8121E+01 "alien"
   25   20   21    1.47   11.22  0.7693E+01  0.1070E+02 "alien"
   28    2   14    1.96   11.20 -0.1003E+02  0.1248E+02 "alien"
   38   28    6    1.59   11.05 -0.2643E+02  0.1165E+02 "alien"
   53    2    0    1.46   11.04  0.1631E+02  0.1568E+02 "alien"
   43    1    0    1.80   10.99 -0.5796E+01  0.1133E+02 "alien"
   34   25    9    1.69   10.13 -0.9446E+01  0.7878E+01 "alien"
   39    2    3    1.96   10.11 -0.9059E+01  0.1108E+02 "alien"


 List of     93 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   32   30   19    1.34   28.20 -0.1615E+02  0.1358E+02
   49    7   15    1.34   20.16 -0.1155E+02  0.1488E+02
   46    8   17    1.34   25.92 -0.1485E+02  0.9179E+01
   33   33   17    1.34   12.82 -0.1782E+02  0.1328E+02
   50    7   14    1.34   27.41 -0.1570E+02  0.1048E+02
   53    6   11    1.34   15.36 -0.8797E+01  0.1046E+02
   51    6   13    1.35   20.37 -0.1167E+02  0.1618E+02
   37   34   14    1.35   19.38 -0.1110E+02  0.9403E+01
   29   22   22    1.36   27.94 -0.1600E+02  0.8317E+01
   39   34   12    1.36   42.69 -0.2445E+02  0.1120E+02
   49    6   14    1.37   24.73 -0.1417E+02  0.1689E+02
   50    6   13    1.37   19.68 -0.1128E+02  0.9208E+01
   43    8   18    1.37   13.63 -0.7806E+01  0.9641E+01
   41    9   19    1.37   30.05 -0.1721E+02  0.1461E+02
   27   13   24    1.37   19.12 -0.1095E+02  0.1053E+02
   46    7   16    1.37   44.55 -0.3868E+01  0.9091E+01
   36   33   14    1.38  239.67 -0.2081E+02  0.1086E+02
   27   18   23    1.38   50.16 -0.4354E+01  0.7949E+01
   39   34   11    1.38   25.30 -0.2197E+01  0.1719E+02
   28   21   22    1.38   11.23 -0.9747E+00  0.8121E+01
   53    5    9    1.38  330.40 -0.2868E+02  0.1653E+02
   35   32   15    1.38  212.95 -0.1849E+02  0.1609E+02
   26   13   24    1.38  123.61 -0.1073E+02  0.9828E+01
   30   28   19    1.39   65.47 -0.5683E+01  0.8560E+01
   48    6   14    1.39  109.01 -0.9463E+01  0.1275E+02
   38    9   20    1.39   86.85 -0.7540E+01  0.1178E+02
   51    5   11    1.39   44.86 -0.3895E+01  0.9318E+01
   38   33   12    1.39   94.03 -0.8163E+01  0.9508E+01
   45    7   16    1.39  202.03 -0.1754E+02  0.1362E+02
   25   13   24    1.39  182.01 -0.1580E+02  0.7281E+01
   53    4    8    1.40   35.90 -0.3116E+01  0.9503E+01
   51    4   10    1.41  204.76 -0.1778E+02  0.1556E+02
   34   31   15    1.41  123.61 -0.1073E+02  0.9022E+01
   40   33    9    1.41  127.80 -0.1109E+02  0.9079E+01
   48    5   13    1.41   18.19 -0.2649E+02  0.1467E+02
   45    6   15    1.42   11.29 -0.1644E+02  0.9503E+01
   38    8   19    1.43   11.89 -0.1732E+02  0.7875E+01
   36   31   13    1.43   14.03 -0.2044E+02  0.9294E+01
   46    5   14    1.43   16.84 -0.2453E+02  0.1450E+02
   42   32    3    1.46   11.80 -0.1719E+02  0.1094E+02
   53    2    1    1.46   17.91  0.1227E+02  0.1134E+02
   53    2    0    1.46   11.04  0.1631E+02  0.1568E+02
   52    3    5    1.46   21.22  0.1455E+02  0.1713E+02
   41   32    5    1.46   11.55  0.7914E+01  0.9045E+01
   26   23   20    1.46   14.54  0.9967E+01  0.1240E+02
   27   25   19    1.47   13.96  0.9568E+01  0.1329E+02
   25   20   21    1.47   11.22  0.7693E+01  0.1070E+02
   43    5   15    1.47   29.13  0.1997E+02  0.9962E+01
   50    3    8    1.47   25.77  0.1766E+02  0.9105E+01
   32   29   15    1.47   15.62  0.1071E+02  0.1345E+02
   40   31    7    1.48   23.94  0.1641E+02  0.9636E+01
   52    2    3    1.48   14.39  0.9862E+01  0.9274E+01
   44    5   14    1.48   19.62  0.1345E+02  0.1571E+02
   52    2    2    1.49   24.63  0.1688E+02  0.9066E+01
   52    2    1    1.49   22.78  0.1561E+02  0.8850E+01
   39    6   17    1.49   16.57  0.1136E+02  0.8009E+01
   26   10   22    1.49   12.74  0.8729E+01  0.9840E+01
   51    2    4    1.50   13.61  0.9329E+01  0.8849E+01
   41   31    2    1.50   14.10  0.9663E+01  0.8967E+01
   34    7   19    1.51   12.27  0.8407E+01  0.1162E+02
   35   29   12    1.51   11.30  0.7744E+01  0.9278E+01
   51    2    3    1.51   11.48  0.7871E+01  0.9208E+01
   38    6   17    1.51   15.42  0.1057E+02  0.9210E+01
   38   28    6    1.59   11.05 -0.2643E+02  0.1165E+02
   29   24   14    1.65   15.57 -0.1453E+02  0.7581E+01
   33   25   11    1.65   17.51 -0.1633E+02  0.1295E+02
   36   26    7    1.67   12.23 -0.1141E+02  0.7909E+01
   27   23   15    1.67   12.91 -0.1204E+02  0.7525E+01
   43    1    8    1.69   12.10 -0.1129E+02  0.1555E+02
   34   25    9    1.69   10.13 -0.9446E+01  0.7878E+01
   28   23   14    1.69   11.87 -0.1108E+02  0.7322E+01
   31    3   15    1.79   44.70 -0.1098E+02  0.7686E+01
   40    0    8    1.80   11.38 -0.6156E+01  0.7995E+01
   43    1    0    1.80   10.99 -0.5796E+01  0.1133E+02
   39    0    9    1.81   16.39 -0.1264E+02  0.8669E+01
   42    1    4    1.81   21.64 -0.5314E+01  0.8309E+01
   17   15   18    1.82   59.82 -0.1469E+02  0.6083E+01
   35    1   12    1.83   60.44 -0.1484E+02  0.7980E+01
   28   21   12    1.83   19.98 -0.4907E+01  0.7359E+01
   39    0    8    1.84   21.81 -0.2592E+02  0.1504E+02
   35   23    0    1.85   12.51 -0.7303E+01  0.1186E+02
   31    2   14    1.85   39.77 -0.9766E+01  0.6943E+01
   20   18   16    1.86   26.67 -0.6550E+01  0.9183E+01
   24   19   14    1.87   95.21 -0.2338E+02  0.1122E+02
   34   22    4    1.88   43.53 -0.1069E+02  0.9236E+01
   33    1   12    1.90   21.39 -0.1916E+02  0.8398E+01
   25   19   13    1.90   11.46 -0.1027E+02  0.7646E+01
   31   21    8    1.90   13.17 -0.1180E+02  0.7987E+01
   29   20   10    1.91   13.20 -0.1183E+02  0.1344E+02
   28    2   14    1.96   11.20 -0.1003E+02  0.1248E+02
   39    2    3    1.96   10.11 -0.9059E+01  0.1108E+02
   30   20    8    1.97   12.48 -0.1118E+02  0.8189E+01
   21   17   14    1.99   15.24 -0.1365E+02  0.1106E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        18
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1281
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1220
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       61


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
