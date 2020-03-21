
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.600    77.600    38.836  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      39      39
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      39
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000992  0.000662
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
 BEAM_DIVERGENCE_E.S.D.=   0.048
 REFLECTING_RANGE_E.S.D.=   0.182
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

    1162 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      39 ...      39
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1162 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               38
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
 *  31        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  35        mP          0.0      77.6   38.8   77.6  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  25        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  34        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  32        oP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  20        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     160.0   38.8   77.6  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     160.0   38.8   77.6  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.8  160.0   77.6  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.8  160.0   77.6  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.8  160.0   77.6  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.8  160.0   77.6  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    27        mC        500.0     160.0   38.8  109.7  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.8  109.7  116.4  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.8  160.0  160.0  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.7  116.4   38.8  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    30        mC        998.1      77.6  173.5   38.8  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.1      77.6  173.5   38.8  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    41        mC        998.1     173.5   77.6   38.8  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    22        hP        998.1      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.4  116.4  116.4  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.8   86.8  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.7   86.8   86.8  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.8   86.8  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.8   86.8  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.8   86.8  248.4 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.8   86.8   77.6  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.8   77.6  177.8  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.8   86.8  177.8 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.8   86.8   86.8  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.8  109.7   38.8  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.8   77.6  177.8 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.8  177.8   77.6 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    4.1      14    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.972   -23.448     0.058
 UNIT_CELL_B-AXIS=    23.448    73.972     0.023
 UNIT_CELL_C-AXIS=    -0.031    -0.002    38.836
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       127 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.33
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.569    77.569    38.827  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.1E-01 1.1E-01 5.9E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012892  0.012892  0.025755  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.920   -23.513     0.026
 COORDINATES OF UNIT CELL B-AXIS    23.513    73.920     0.009
 COORDINATES OF UNIT CELL C-AXIS    -0.014    -0.001    38.827
 CRYSTAL MOSAICITY (DEGREES)     0.182
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999979  0.006269  0.001733
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.92   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1261.40   1242.01
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.29
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.57    77.57    38.83  90.000  90.000  90.000


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
  0.05  0.15        22         26.590       3.891
  0.15  0.25       149         98.813       1.047
  0.25  0.35       159        109.552       0.944
  0.35  0.45       163         52.718       1.962
  0.45  0.55       166        134.860       0.767
  0.55  0.65       165        107.085       0.966
  0.65  0.75       194        107.887       0.959
  0.75  0.85       132        130.893       0.790
  0.85  0.95        12         51.839       1.996
 
 -0.05  0.95      1162        103.447       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        6
 DEGREES OF FREEDOM OF CHI^2 FIT               151.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.835
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    38.1 XMAX=    38.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59263 NYBIN=    6
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        328


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               151.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.829
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    13.4 XMAX=  2452.2 NXBIN=    3
 YMIN=     5.6 YMAX=  2518.0 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        328


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               148.7
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.826
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    38.1 XMAX=    38.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1262
 DETECTOR_SURFACE_POSITION=    1647    1689
 DETECTOR_SURFACE_POSITION=     819    1689
 DETECTOR_SURFACE_POSITION=     819     835
 DETECTOR_SURFACE_POSITION=    1647     835
 DETECTOR_SURFACE_POSITION=    2170    1662
 DETECTOR_SURFACE_POSITION=    1621    2227
 DETECTOR_SURFACE_POSITION=     845    2227
 DETECTOR_SURFACE_POSITION=     295    1662
 DETECTOR_SURFACE_POSITION=     295     862
 DETECTOR_SURFACE_POSITION=     845     296
 DETECTOR_SURFACE_POSITION=    1621     296
 DETECTOR_SURFACE_POSITION=    2170     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        328


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
 7.085E-01  7.027E-04   44.82


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

   24.530   2.888    22.70   1.34      1.95      1.65      76     105      22
    2.888   2.050     3.50   1.14     18.43     16.35     182     224       4
    2.050   1.675     0.04   0.91   1082.77   1078.05     218     264       6
    1.675   1.452    -0.05   0.94    -99.90    -99.90     248     317       0
    1.452   1.299     0.06   0.80    625.14    697.67     162     215       4
   --------------------------------------------------------------------------
   24.530   1.299     2.82   0.98      8.98      8.73     886    1125      36
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      39      39)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.686  -0.1458E+02  0.1203E+02        -1.21           1 
   47    0    0     1.650  -0.1436E+02  0.1109E+02        -1.29           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS   -13.6% OF MEAN INTENSITY


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

     3.85          54      35      1285        2.7%       1.8%      1.4%       38   36.86      2.5%   100.0*     0    0.000       0
     2.74          71      45      2098        2.1%       4.0%      4.1%       52   16.90      5.6%    99.7*     0    0.000       0
     2.24         136      80      2664        3.0%      15.7%     12.8%      112    5.02     22.3%    95.5*     0    0.000       0
     1.95         126      73      3102        2.4%      81.5%     76.1%      106    0.91    115.2%    49.6*     0    0.000       0
     1.74         133      80      3481        2.3%     -99.9%    -99.9%      106    0.00    -99.9%   -11.5      0    0.000       0
     1.59         177     103      3848        2.7%    1541.7%   1472.2%      148    0.06   2180.3%     2.5      0    0.000       0
     1.47         176     109      4142        2.6%     -99.9%    -99.9%      134    0.00    -99.9%    24.5      0    0.000       0
     1.38         173     103      4433        2.3%   22898.2%  26530.9%      140    0.02  32383.0%    16.5      0    0.000       0
     1.30          79      54      4737        1.1%     272.6%    272.8%       50    0.11    385.5%    -7.5      0    0.000       0
    total        1125     682     29790        2.3%       9.0%      8.7%      886    3.70     12.7%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1162
 NUMBER OF REJECTED MISFITS                                36
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1125
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    682



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

   WILSON LINE (using all data) : A=   6.412 B=  27.768 CORRELATION=  0.82
      #      RES      SS        <I>       log(<I>)       BO
      14     7.501   0.004  1.8185E+03       7.506    -123.1
      17     4.753   0.011  1.6101E+03       7.384     -43.9
      17     3.574   0.020  5.9755E+02       6.393       0.5
      26     2.988   0.028  1.3739E+02       4.923      26.6
      29     2.641   0.036  9.9470E+01       4.600      25.3
      45     2.407   0.043  3.7492E+01       3.624      32.3
      33     2.211   0.051  2.0557E+01       3.023      33.1
      35     2.061   0.059  5.2232E+00       1.653      40.4
      36     1.939   0.066 -1.4117E+00    -999.999   -1000.0
      43     1.831   0.075  2.8043E-01      -1.271      51.5
      41     1.742   0.082  2.5628E+00       0.941      33.2
      43     1.671   0.090 -6.6031E-01    -999.999   -1000.0
      48     1.602   0.097  7.4100E-01      -0.300      34.4
      53     1.541   0.105 -1.2430E+00    -999.999   -1000.0
      62     1.485   0.113 -1.1509E-01    -999.999   -1000.0
      46     1.437   0.121  1.8580E+00       0.619      23.9
      41     1.395   0.128 -7.0997E-01    -999.999   -1000.0
      37     1.354   0.136 -2.3693E-02    -999.999   -1000.0
      16     1.317   0.144  4.0325E+00       1.394      17.4


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.501        2.065        1.580        0.944
       3     4.753        1.400        0.908        0.437
       2     3.574        0.336        0.068        0.010
       3     2.988        2.170        1.195        0.490
       1     2.641        0.015        0.001        0.000
       4     2.407        0.415        0.118        0.026
       1     2.211        0.737        0.219        0.047
       6     2.061        4.105        3.744        2.654
       7     1.742        3.367        3.244        2.820
       3     1.602       60.291      277.503      920.571
       3     1.541        9.375       12.839       17.112
       5     1.485       56.815       99.351      351.681
       3     1.437       34.808       52.666      115.971
       1     1.354       26.716       47.834       61.176
      74   overall      121.353     -529.822    22165.805


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.501        0.592        0.473        0.291
      14     4.753        1.276        1.616        1.623
      15     3.574        1.034        1.084        1.030
      23     2.988        0.740        0.600        0.404
      28     2.641        1.442        2.350        3.561
      41     2.407        1.147        1.168        1.030
      32     2.211        0.916        0.662        0.406
      29     2.061        1.893        1.159        1.515
      33     1.831      422.392    -4591.152   194791.752
      34     1.742        6.202        5.296       18.151
      32     1.671       46.041      133.441     1579.411
      45     1.602       88.572       79.919     2865.409
      43     1.437       12.322        9.385       83.601
      16     1.317        6.869        1.453       12.992
     608   overall       31.255     -241.906    10792.336

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.501  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.73
       3     4.753  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.574  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.73
       3     2.988  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     2.641  0.00  0.00  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       4     2.407  0.00  0.00  0.00  0.00  0.00  0.89  0.84  0.80  0.76  0.73
       1     2.211  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.061  0.67  0.48  0.40  0.35  0.32  0.30  0.28  0.27  0.25  0.24
       7     1.742  0.58  0.41  0.69  0.60  0.55  0.51  0.48  0.45  0.43  0.42
       3     1.602  1.34  0.96  0.80  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.541  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       5     1.485  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
       3     1.437  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       1     1.354  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      74   overall  0.60  0.43  0.42  0.40  0.36  0.39  0.36  0.34  0.35  0.40


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.501  4.20  3.86  2.70  2.12  2.03  1.77  1.59  1.45  1.35  1.27
      14     4.753  3.00  1.58  1.65  1.30  1.27  1.58  1.42  1.30  1.44  1.36
      15     3.574  1.40  0.74  0.77  0.61  0.85  0.74  1.19  1.21  1.12  1.05
      23     2.988  1.37  1.44  1.51  1.45  1.66  1.45  1.30  1.26  1.17  1.24
      28     2.641  1.50  1.18  1.10  1.08  0.91  1.03  1.06  0.97  1.08  1.07
      41     2.407  1.54  1.35  1.13  1.04  1.18  1.03  1.02  1.02  1.07  1.08
      32     2.211  1.64  1.21  0.96  0.95  0.79  0.83  0.87  0.91  0.95  0.99
      29     2.061  3.26  1.90  1.33  1.15  1.05  0.99  0.96  0.94  0.87  0.87
      33     1.831  4.78  2.51  1.75  1.38  1.16  1.07  0.96  0.88  0.82  0.77
      34     1.742  4.33  2.60  1.93  1.52  1.35  1.17  1.05  0.96  0.89  0.84
      32     1.671  5.91  3.10  2.17  1.80  1.51  1.32  1.18  1.08  1.00  0.94
      45     1.602  4.90  2.57  1.80  1.42  1.19  1.03  0.93  0.85  0.79  0.74
      43     1.437  4.64  2.44  1.79  1.41  1.18  1.08  1.02  0.93  0.86  0.81
      16     1.317  3.28  2.07  1.45  1.14  0.95  0.83  0.74  0.68  0.63  0.59
     608   overall  2.23  1.33  1.02  0.85  0.77  0.71  0.68  0.64  0.62  0.61


 List of    113 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   34   11   22    1.37  764.88 -0.1812E+02  0.8056E+01 "alien"
   49    7   15    1.34  564.46 -0.1337E+02  0.1476E+02 "alien"
   35   11   22    1.35  442.93 -0.1049E+02  0.1134E+02 "alien"
   36   34   15    1.34  417.73 -0.9897E+01  0.1342E+02 "alien"
   37   10   21    1.37  401.86 -0.9522E+01  0.1211E+02 "alien"
   37   34   14    1.35  296.29 -0.7020E+01  0.9135E+01 "alien"
   28   19   23    1.36  268.05 -0.6351E+01  0.8095E+01 "alien"
   48    3    9    1.51  258.39 -0.2974E+02  0.1591E+02 "alien"
   45    8   17    1.36  247.29 -0.5859E+01  0.1479E+02 "alien"
   35   33   15    1.37  245.09 -0.5807E+01  0.1442E+02 "alien"
   39   31    9    1.46  239.26 -0.2754E+02  0.9695E+01 "alien"
   30   25   21    1.35  231.75 -0.5491E+01  0.1176E+02 "alien"
   52    2    1    1.49  211.91 -0.2439E+02  0.8454E+01 "alien"
   28   17   24    1.34  211.35 -0.5008E+01  0.1218E+02 "alien"
   52    5   11    1.37  205.18 -0.4861E+01  0.1517E+02 "alien"
   46    4   12    1.49  187.04 -0.2153E+02  0.9569E+01 "alien"
   53    5   10    1.36  179.51 -0.4253E+01  0.9290E+01 "alien"
   34   29   13    1.50  173.14 -0.1993E+02  0.9850E+01 "alien"
   31   29   19    1.36  144.61 -0.3426E+01  0.8508E+01 "alien"
   42    9   19    1.35  144.26 -0.3418E+01  0.8765E+01 "alien"
   27   25   19    1.47  135.37 -0.1558E+02  0.1298E+02 "alien"
   30   28   16    1.49  135.04 -0.1554E+02  0.1375E+02 "alien"
   49    6   14    1.37  129.60 -0.3071E+01  0.1542E+02 "alien"
   42   32    1    1.47  117.02 -0.1347E+02  0.8962E+01 "alien"
   41   31    2    1.50  107.28 -0.1235E+02  0.8540E+01 "alien"
   41   31    3    1.50   93.45 -0.1075E+02  0.8233E+01 "alien"
   40   31    6    1.49   68.97 -0.7938E+01  0.8485E+01 "alien"
   38   31   10    1.46   66.78 -0.7685E+01  0.8330E+01 "alien"
   32   31   18    1.35   66.57 -0.1577E+01  0.8920E+01 "alien"
   35   29   12    1.51   60.85 -0.7003E+01  0.8976E+01 "alien"
   52    2    4    1.47   58.11 -0.6688E+01  0.1618E+02 "alien"
   20   18   16    1.86   57.88  0.1623E+02  0.9499E+01 "alien"
   29    9   21    1.50   55.05 -0.6336E+01  0.7106E+01 "alien"
   38    7   18    1.47   50.84 -0.5852E+01  0.1176E+02 "alien"
   35   23    3    1.83   38.49  0.1079E+02  0.8135E+01 "alien"
   47    4   12    1.47   37.47 -0.4312E+01  0.8694E+01 "alien"
   30   21   10    1.86   33.04  0.9264E+01  0.1273E+02 "alien"
   30   22   11    1.80   32.12  0.9007E+01  0.7630E+01 "alien"
   41    1    5    1.84   30.84  0.8648E+01  0.1136E+02 "alien"
   31   24   12    1.69   30.56 -0.2018E+02  0.7308E+01 "alien"
   34    7   19    1.51   30.40 -0.3498E+01  0.1099E+02 "alien"
   38    6   17    1.51   30.32 -0.3490E+01  0.1177E+02 "alien"
   43    1    1    1.80   29.32  0.8222E+01  0.7808E+01 "alien"
   39   28    3    1.60   28.08  0.2080E+02  0.8585E+01 "alien"
   34   23    5    1.84   28.06  0.7869E+01  0.9363E+01 "alien"
   34   27   11    1.59   26.65  0.1975E+02  0.1355E+02 "alien"
   30    9   21    1.48   26.25 -0.3021E+01  0.9993E+01 "alien"
   43    7   17    1.40   25.59 -0.1817E+02  0.9016E+01 "alien"
   40    1    5    1.88   24.50  0.6869E+01  0.1170E+02 "alien"
   32   12   23    1.35   24.13 -0.5717E+00  0.1116E+02 "alien"
   48    1    3    1.60   23.89  0.1770E+02  0.1484E+02 "alien"
   54    4    6    1.40   23.73 -0.1685E+02  0.9872E+01 "alien"
   51    5   11    1.39   22.49 -0.1597E+02  0.8900E+01 "alien"
   41   31    0    1.51   22.20 -0.1279E+02  0.1290E+02 "alien"
   42    1    4    1.81   22.10  0.6199E+01  0.8252E+01 "alien"
   45    7   16    1.39   21.39 -0.1518E+02  0.1473E+02 "alien"
   36   32   14    1.39   20.71 -0.1470E+02  0.1604E+02 "alien"
   34    1   12    1.86   20.46  0.5738E+01  0.7275E+01 "alien"
   47    1    6    1.60   20.27  0.1502E+02  0.1140E+02 "alien"
   26   10   22    1.49   19.76 -0.2274E+01  0.9754E+01 "alien"
   41    0    7    1.79   19.66  0.2234E+02  0.1360E+02 "alien"
   41    4   14    1.56   19.34 -0.2404E+02  0.9043E+01 "alien"
   33    8   20    1.48   19.23 -0.2213E+01  0.1125E+02 "alien"
   37   30   11    1.48   18.07 -0.2080E+01  0.1024E+02 "alien"
   26   23   20    1.47   17.54 -0.2019E+01  0.1184E+02 "alien"
   53    2    0    1.46   17.49 -0.6876E+01  0.1504E+02 "alien"
   32   27   14    1.54   17.17 -0.2134E+02  0.1283E+02 "alien"
   39    4   15    1.57   16.74  0.1240E+02  0.1120E+02 "alien"
   34   31   15    1.41   16.74 -0.1188E+02  0.8656E+01 "alien"
   17   17   17    1.86   16.05  0.1378E+02  0.8771E+01 "alien"
   21   19   19    1.64   15.55 -0.1027E+02  0.9892E+01 "alien"
   51    4   10    1.41   15.50 -0.1100E+02  0.1486E+02 "alien"
   37   26    4    1.69   15.20 -0.1004E+02  0.7946E+01 "alien"
   32   24   11    1.70   15.17 -0.1001E+02  0.8269E+01 "alien"
   35   10   21    1.40   15.16 -0.1076E+02  0.1176E+02 "alien"
   26   25   17    1.57   15.02 -0.1867E+02  0.1199E+02 "alien"
   27   18   23    1.38   14.97 -0.1062E+02  0.1113E+02 "alien"
   37   26    3    1.70   14.87 -0.9817E+01  0.9854E+01 "alien"
   34    4   16    1.66   14.68 -0.9696E+01  0.7372E+01 "alien"
   38    9   20    1.39   14.63 -0.1039E+02  0.1158E+02 "alien"
   33    5   17    1.63   14.42  0.1069E+02  0.7483E+01 "alien"
   44    3   12    1.55   14.41 -0.1791E+02  0.8533E+01 "alien"
   34   25    9    1.69   14.26 -0.9417E+01  0.7546E+01 "alien"
   36   24    0    1.79   14.24  0.1070E+02  0.1112E+02 "alien"
   36   27    9    1.60   13.98  0.1036E+02  0.9223E+01 "alien"
   39   29    5    1.56   13.73 -0.1707E+02  0.8119E+01 "alien"
   45    1    8    1.62   13.33  0.9877E+01  0.9426E+01 "alien"
   33   27   13    1.55   12.75 -0.1585E+02  0.1353E+02 "alien"
   36   27    8    1.62   12.45  0.9224E+01  0.9648E+01 "alien"
   27    9   21    1.53   12.31 -0.1530E+02  0.9739E+01 "alien"
   28   27   17    1.50   12.12 -0.1395E+01  0.8797E+01 "alien"
   46    3   11    1.52   12.09 -0.1503E+02  0.8338E+01 "alien"
   55    3    0    1.41   11.82 -0.1904E+02  0.1319E+02 "alien"
   40    1    6    1.86   11.79  0.3306E+01  0.8586E+01 "alien"
   31    2   14    1.86   11.78  0.3303E+01  0.6791E+01 "alien"
   45    2   10    1.57   11.77  0.8721E+01  0.9395E+01 "alien"
   47    0    2    1.64   11.74 -0.1750E+02  0.8869E+01 "alien"
   46    5   14    1.43   11.66  0.2167E+02  0.1409E+02 "alien"
   46    1    7    1.61   11.63  0.8614E+01  0.8056E+01 "alien"
   46    6   15    1.40   11.40 -0.8091E+01  0.8816E+01 "alien"
   52    3    7    1.44   11.36  0.2111E+02  0.1026E+02 "alien"
   41   31    4    1.49   11.10 -0.1277E+01  0.9001E+01 "alien"
   39    8   19    1.41   11.00 -0.7813E+01  0.1240E+02 "alien"
   39   28    0    1.62   10.88  0.1720E+02  0.1189E+02 "alien"
   46    0    4    1.66   10.68 -0.1484E+02  0.8902E+01 "alien"
   27   20   22    1.40   10.53 -0.7479E+01  0.7791E+01 "alien"
   37   31   12    1.44   10.50  0.1951E+02  0.1545E+02 "alien"
   37   27    7    1.62   10.49  0.7775E+01  0.8383E+01 "alien"
   17   15   18    1.82   10.39  0.2914E+01  0.7864E+01 "alien"
   29    2   14    1.92   10.39 -0.1466E+02  0.1046E+02 "alien"
   38    0    9    1.85   10.34  0.5979E+01  0.7781E+01 "alien"
   42    4   14    1.53   10.32 -0.1283E+02  0.8037E+01 "alien"
   44    0    6    1.70   10.10 -0.1354E+02  0.8428E+01 "alien"


 List of    113 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   28   17   24    1.34  211.35 -0.5008E+01  0.1218E+02
   36   34   15    1.34  417.73 -0.9897E+01  0.1342E+02
   49    7   15    1.34  564.46 -0.1337E+02  0.1476E+02
   37   34   14    1.35  296.29 -0.7020E+01  0.9135E+01
   42    9   19    1.35  144.26 -0.3418E+01  0.8765E+01
   30   25   21    1.35  231.75 -0.5491E+01  0.1176E+02
   32   12   23    1.35   24.13 -0.5717E+00  0.1116E+02
   32   31   18    1.35   66.57 -0.1577E+01  0.8920E+01
   35   11   22    1.35  442.93 -0.1049E+02  0.1134E+02
   28   19   23    1.36  268.05 -0.6351E+01  0.8095E+01
   31   29   19    1.36  144.61 -0.3426E+01  0.8508E+01
   45    8   17    1.36  247.29 -0.5859E+01  0.1479E+02
   53    5   10    1.36  179.51 -0.4253E+01  0.9290E+01
   37   10   21    1.37  401.86 -0.9522E+01  0.1211E+02
   49    6   14    1.37  129.60 -0.3071E+01  0.1542E+02
   35   33   15    1.37  245.09 -0.5807E+01  0.1442E+02
   52    5   11    1.37  205.18 -0.4861E+01  0.1517E+02
   34   11   22    1.37  764.88 -0.1812E+02  0.8056E+01
   27   18   23    1.38   14.97 -0.1062E+02  0.1113E+02
   38    9   20    1.39   14.63 -0.1039E+02  0.1158E+02
   51    5   11    1.39   22.49 -0.1597E+02  0.8900E+01
   36   32   14    1.39   20.71 -0.1470E+02  0.1604E+02
   45    7   16    1.39   21.39 -0.1518E+02  0.1473E+02
   35   10   21    1.40   15.16 -0.1076E+02  0.1176E+02
   54    4    6    1.40   23.73 -0.1685E+02  0.9872E+01
   27   20   22    1.40   10.53 -0.7479E+01  0.7791E+01
   43    7   17    1.40   25.59 -0.1817E+02  0.9016E+01
   46    6   15    1.40   11.40 -0.8091E+01  0.8816E+01
   55    3    0    1.41   11.82 -0.1904E+02  0.1319E+02
   39    8   19    1.41   11.00 -0.7813E+01  0.1240E+02
   51    4   10    1.41   15.50 -0.1100E+02  0.1486E+02
   34   31   15    1.41   16.74 -0.1188E+02  0.8656E+01
   46    5   14    1.43   11.66  0.2167E+02  0.1409E+02
   52    3    7    1.44   11.36  0.2111E+02  0.1026E+02
   37   31   12    1.44   10.50  0.1951E+02  0.1545E+02
   53    2    0    1.46   17.49 -0.6876E+01  0.1504E+02
   39   31    9    1.46  239.26 -0.2754E+02  0.9695E+01
   38   31   10    1.46   66.78 -0.7685E+01  0.8330E+01
   26   23   20    1.47   17.54 -0.2019E+01  0.1184E+02
   47    4   12    1.47   37.47 -0.4312E+01  0.8694E+01
   27   25   19    1.47  135.37 -0.1558E+02  0.1298E+02
   42   32    1    1.47  117.02 -0.1347E+02  0.8962E+01
   38    7   18    1.47   50.84 -0.5852E+01  0.1176E+02
   52    2    4    1.47   58.11 -0.6688E+01  0.1618E+02
   37   30   11    1.48   18.07 -0.2080E+01  0.1024E+02
   33    8   20    1.48   19.23 -0.2213E+01  0.1125E+02
   30    9   21    1.48   26.25 -0.3021E+01  0.9993E+01
   52    2    1    1.49  211.91 -0.2439E+02  0.8454E+01
   26   10   22    1.49   19.76 -0.2274E+01  0.9754E+01
   46    4   12    1.49  187.04 -0.2153E+02  0.9569E+01
   41   31    4    1.49   11.10 -0.1277E+01  0.9001E+01
   30   28   16    1.49  135.04 -0.1554E+02  0.1375E+02
   40   31    6    1.49   68.97 -0.7938E+01  0.8485E+01
   29    9   21    1.50   55.05 -0.6336E+01  0.7106E+01
   41   31    3    1.50   93.45 -0.1075E+02  0.8233E+01
   34   29   13    1.50  173.14 -0.1993E+02  0.9850E+01
   28   27   17    1.50   12.12 -0.1395E+01  0.8797E+01
   41   31    2    1.50  107.28 -0.1235E+02  0.8540E+01
   34    7   19    1.51   30.40 -0.3498E+01  0.1099E+02
   41   31    0    1.51   22.20 -0.1279E+02  0.1290E+02
   35   29   12    1.51   60.85 -0.7003E+01  0.8976E+01
   48    3    9    1.51  258.39 -0.2974E+02  0.1591E+02
   38    6   17    1.51   30.32 -0.3490E+01  0.1177E+02
   46    3   11    1.52   12.09 -0.1503E+02  0.8338E+01
   27    9   21    1.53   12.31 -0.1530E+02  0.9739E+01
   42    4   14    1.53   10.32 -0.1283E+02  0.8037E+01
   32   27   14    1.54   17.17 -0.2134E+02  0.1283E+02
   44    3   12    1.55   14.41 -0.1791E+02  0.8533E+01
   33   27   13    1.55   12.75 -0.1585E+02  0.1353E+02
   41    4   14    1.56   19.34 -0.2404E+02  0.9043E+01
   39   29    5    1.56   13.73 -0.1707E+02  0.8119E+01
   26   25   17    1.57   15.02 -0.1867E+02  0.1199E+02
   39    4   15    1.57   16.74  0.1240E+02  0.1120E+02
   45    2   10    1.57   11.77  0.8721E+01  0.9395E+01
   34   27   11    1.59   26.65  0.1975E+02  0.1355E+02
   47    1    6    1.60   20.27  0.1502E+02  0.1140E+02
   36   27    9    1.60   13.98  0.1036E+02  0.9223E+01
   48    1    3    1.60   23.89  0.1770E+02  0.1484E+02
   39   28    3    1.60   28.08  0.2080E+02  0.8585E+01
   46    1    7    1.61   11.63  0.8614E+01  0.8056E+01
   39   28    0    1.62   10.88  0.1720E+02  0.1189E+02
   37   27    7    1.62   10.49  0.7775E+01  0.8383E+01
   45    1    8    1.62   13.33  0.9877E+01  0.9426E+01
   36   27    8    1.62   12.45  0.9224E+01  0.9648E+01
   33    5   17    1.63   14.42  0.1069E+02  0.7483E+01
   21   19   19    1.64   15.55 -0.1027E+02  0.9892E+01
   47    0    2    1.64   11.74 -0.1750E+02  0.8869E+01
   34    4   16    1.66   14.68 -0.9696E+01  0.7372E+01
   46    0    4    1.66   10.68 -0.1484E+02  0.8902E+01
   31   24   12    1.69   30.56 -0.2018E+02  0.7308E+01
   37   26    4    1.69   15.20 -0.1004E+02  0.7946E+01
   34   25    9    1.69   14.26 -0.9417E+01  0.7546E+01
   32   24   11    1.70   15.17 -0.1001E+02  0.8269E+01
   37   26    3    1.70   14.87 -0.9817E+01  0.9854E+01
   44    0    6    1.70   10.10 -0.1354E+02  0.8428E+01
   41    0    7    1.79   19.66  0.2234E+02  0.1360E+02
   36   24    0    1.79   14.24  0.1070E+02  0.1112E+02
   30   22   11    1.80   32.12  0.9007E+01  0.7630E+01
   43    1    1    1.80   29.32  0.8222E+01  0.7808E+01
   42    1    4    1.81   22.10  0.6199E+01  0.8252E+01
   17   15   18    1.82   10.39  0.2914E+01  0.7864E+01
   35   23    3    1.83   38.49  0.1079E+02  0.8135E+01
   34   23    5    1.84   28.06  0.7869E+01  0.9363E+01
   41    1    5    1.84   30.84  0.8648E+01  0.1136E+02
   38    0    9    1.85   10.34  0.5979E+01  0.7781E+01
   31    2   14    1.86   11.78  0.3303E+01  0.6791E+01
   20   18   16    1.86   57.88  0.1623E+02  0.9499E+01
   40    1    6    1.86   11.79  0.3306E+01  0.8586E+01
   30   21   10    1.86   33.04  0.9264E+01  0.1273E+02
   34    1   12    1.86   20.46  0.5738E+01  0.7275E+01
   17   17   17    1.86   16.05  0.1378E+02  0.8771E+01
   40    1    5    1.88   24.50  0.6869E+01  0.1170E+02
   29    2   14    1.92   10.39 -0.1466E+02  0.1046E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        59
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1162
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1037
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      125


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
