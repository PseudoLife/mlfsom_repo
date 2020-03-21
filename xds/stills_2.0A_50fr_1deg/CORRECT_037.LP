
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.582    77.582    38.824  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      37      37
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      37
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000986  0.000665
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
 REFLECTING_RANGE_E.S.D.=   0.173
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

    1192 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      37 ...      37
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1192 REFLECTIONS ACCEPTED



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
 *  31        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  35        mP          0.0      77.6   38.8   77.6  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  33        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.9   38.8   77.6  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.9   38.8   77.6  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.8  159.9   77.6  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.8  159.9   77.6  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.8  159.9   77.6  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.8  159.9   77.6  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    27        mC        500.0     159.9   38.8  109.7  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.8  109.7  116.4  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.8  159.9  159.9  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.7  116.4   38.8  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    30        mC        998.3      77.6  173.5   38.8  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.3      77.6  173.5   38.8  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    41        mC        998.3     173.5   77.6   38.8  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    22        hP        998.3      77.6   77.6   38.8  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
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

  *   0.00        0    3.7      14    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.945   -23.475     0.058
 UNIT_CELL_B-AXIS=    23.475    73.945     0.024
 UNIT_CELL_C-AXIS=    -0.031    -0.002    38.824
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       139 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.33
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.523    77.523    38.803  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  9.8E-02 9.8E-02 3.5E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012899  0.012899  0.025771  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.886   -23.466     0.054
 COORDINATES OF UNIT CELL B-AXIS    23.466    73.886     0.021
 COORDINATES OF UNIT CELL C-AXIS    -0.029    -0.002    38.803
 CRYSTAL MOSAICITY (DEGREES)     0.173
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999975  0.006735  0.002082
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.88   1232.92
 DETECTOR ORIGIN (PIXELS) AT                     1260.85   1242.72
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.17
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.52    77.52    38.80  90.000  90.000  90.000


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
  0.05  0.15        25         54.112       1.989
  0.15  0.25       155        102.223       1.053
  0.25  0.35       168        108.252       0.994
  0.35  0.45       165         57.528       1.871
  0.45  0.55       164        145.196       0.741
  0.55  0.65       168        111.617       0.964
  0.65  0.75       200        125.980       0.854
  0.75  0.85       131        107.843       0.998
  0.85  0.95        16         97.186       1.108
 
 -0.05  0.95      1192        107.649       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        7
 DEGREES OF FREEDOM OF CHI^2 FIT               182.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.788
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    36.1 XMAX=    36.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59335 NYBIN=    7
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        388


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               182.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.783
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    13.2 XMAX=  2452.5 NXBIN=    3
 YMIN=     5.8 YMAX=  2518.9 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        388


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               182.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.778
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    36.1 XMAX=    36.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1262
 DETECTOR_SURFACE_POSITION=    1647    1689
 DETECTOR_SURFACE_POSITION=     819    1689
 DETECTOR_SURFACE_POSITION=     819     836
 DETECTOR_SURFACE_POSITION=    1647     836
 DETECTOR_SURFACE_POSITION=    2170    1662
 DETECTOR_SURFACE_POSITION=    1621    2228
 DETECTOR_SURFACE_POSITION=     844    2228
 DETECTOR_SURFACE_POSITION=     295    1662
 DETECTOR_SURFACE_POSITION=     295     862
 DETECTOR_SURFACE_POSITION=     844     296
 DETECTOR_SURFACE_POSITION=    1621     296
 DETECTOR_SURFACE_POSITION=    2170     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        388


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
 7.218E-01  9.252E-04   38.70


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

   24.515   2.887    20.89   1.16      1.81      1.81      76     104      24
    2.887   2.048     4.39   1.11     11.63     11.26     198     237       2
    2.048   1.674     0.00   0.97   2445.16   2450.80     224     271       2
    1.674   1.451    -0.02   0.95   8744.40   8814.67     254     320       4
    1.451   1.298     0.01   0.79    -99.90    -99.90     162     219       8
   --------------------------------------------------------------------------
   24.515   1.298     2.79   0.98      8.18      8.37     914    1151      40
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      37      37)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.685  -0.1675E+01  0.1199E+02        -0.14           1 
   47    0    0     1.649   0.3141E+02  0.1211E+02         2.60           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    28.4% OF MEAN INTENSITY


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

     3.85          54      35      1285        2.7%       1.6%      1.6%       38   33.22      2.3%    99.9*     0    0.000       0
     2.74          72      44      2098        2.1%       3.4%      3.7%       56   17.03      4.8%    99.7*     0    0.000       0
     2.24         146      85      2664        3.2%      11.9%     10.6%      122    6.53     16.9%    96.4*     0    0.000       0
     1.94         126      73      3102        2.4%      63.1%     68.1%      106    0.89     89.2%    62.8*     0    0.000       0
     1.74         139      83      3481        2.4%     -99.9%    -99.9%      112    0.00    -99.9%     7.3      0    0.000       0
     1.59         177     103      3848        2.7%     680.1%    628.4%      148   -0.01    961.9%     3.4      0    0.000       0
     1.47         176     108      4142        2.6%     -99.9%    -99.9%      136    0.00    -99.9%     0.6      0    0.000       0
     1.38         178     107      4433        2.4%     -99.9%    -99.9%      142    0.00    -99.9%    14.5      0    0.000       0
     1.30          81      55      4733        1.2%     -99.9%    -99.9%       52    0.00    -99.9%    -7.3      0    0.000       0
    total        1149     693     29786        2.3%       8.1%      8.3%      912    3.64     11.5%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1192
 NUMBER OF REJECTED MISFITS                                40
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1151
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    694



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

   WILSON LINE (using all data) : A=   6.929 B=  32.441 CORRELATION=  0.94
      #      RES      SS        <I>       log(<I>)       BO
      14     7.403   0.005  1.9115E+03       7.556     -68.7
      18     4.749   0.011  1.5976E+03       7.376     -20.2
      17     3.565   0.020  6.2688E+02       6.441      12.4
      24     2.977   0.028  1.3678E+02       4.918      35.6
      32     2.638   0.036  1.2772E+02       4.850      28.9
      46     2.406   0.043  5.2023E+01       3.952      34.5
      34     2.211   0.051  1.9927E+01       2.992      38.5
      35     2.059   0.059  6.4575E+00       1.865      42.9
      38     1.936   0.067 -1.4454E+00    -999.999   -1000.0
      45     1.830   0.075 -1.1022E+00    -999.999   -1000.0
      40     1.741   0.082  6.5939E-01      -0.416      44.5
      42     1.670   0.090  9.2043E-01      -0.083      39.1
      50     1.600   0.098 -4.0322E-01    -999.999   -1000.0
      52     1.540   0.105 -1.3904E+00    -999.999   -1000.0
      62     1.484   0.114  4.2014E-01      -0.867      34.3
      50     1.436   0.121  3.0712E-01      -1.181      33.5
      40     1.393   0.129  1.5267E-01      -1.879      34.2
      37     1.353   0.137  3.4339E-01      -1.069      29.3
      17     1.317   0.144  1.1285E+00       0.121      23.6


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.403        1.912        1.395        0.796
       3     4.749        1.559        1.064        0.539
       2     3.565        0.346        0.071        0.010
       3     2.977        2.685        1.603        0.704
       1     2.638        0.000        0.000        0.000
       4     2.406        0.280        0.074        0.015
       1     2.211        0.381        0.082        0.012
       6     2.059        1.714        0.359        0.493
       6     1.936       53.788       33.022      452.257
       7     1.741      112.674      614.523     3496.242
      11     1.670       21.120        5.834      101.832
       1     1.393    46442.029  3467028.395184873889.233
      74   overall      667.114    46714.126  2500102.581


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.403        0.558        0.416        0.244
      15     4.749        1.326        1.770        1.874
      15     3.565        1.040        1.092        1.036
      21     2.977        0.753        0.718        0.571
      31     2.638        1.210        1.508        1.789
      42     2.406        1.020        0.979        0.837
      33     2.211        1.069        0.905        0.665
      29     2.059        2.456        2.107        3.005
      35     1.830       44.831       50.682      891.908
      31     1.670       59.906       95.758     1510.880
      57     1.484      460.558      609.434   153440.127
      47     1.436      470.669     1257.161   185042.098
      36     1.353      651.574     8076.345   249308.505
      17     1.317       86.366      134.679     2502.583
     619   overall      303.225      509.924   297193.477

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.403  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.73
       3     4.749  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     3.565  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.73
       3     2.977  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     2.638  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       4     2.406  1.01  0.72  0.60  0.53  0.96  0.89  0.84  0.80  1.14  1.10
       1     2.211  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.059  2.01  1.45  1.20  1.06  1.28  1.19  1.12  1.06  1.01  0.98
       6     1.936  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       7     1.741  2.30  1.65  1.37  1.21  1.10  1.02  0.96  0.91  0.87  0.84
      11     1.670  1.46  1.05  0.87  0.77  0.70  0.65  0.61  0.58  0.55  0.53
       1     1.393  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      74   overall  0.98  0.70  0.58  0.51  0.52  0.48  0.45  0.43  0.45  0.48


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.403  3.15  3.31  2.31  2.12  1.78  1.55  1.59  1.45  1.35  1.27
      15     4.749  3.50  1.84  1.80  1.42  1.19  1.18  1.46  1.33  1.46  1.37
      15     3.565  1.40  0.74  0.77  0.81  0.85  0.89  1.19  1.09  1.12  1.05
      21     2.977  2.00  2.10  1.47  1.59  1.45  1.48  1.42  1.38  1.36  1.28
      31     2.638  0.68  0.89  1.00  1.08  1.07  1.00  1.03  1.00  1.03  1.12
      42     2.406  0.00  0.13  0.37  0.94  1.03  1.16  1.09  0.99  1.04  1.13
      33     2.211  2.23  1.34  0.94  0.83  0.85  0.87  0.96  0.99  1.02  1.01
      29     2.059  2.90  1.71  1.20  1.05  0.96  0.92  0.82  0.75  0.76  0.82
      35     1.830  5.10  2.68  1.87  1.47  1.23  1.08  0.96  0.88  0.82  0.81
      31     1.670  5.08  2.67  1.87  1.47  1.23  1.07  0.96  0.88  0.82  0.77
      57     1.484  5.16  2.71  1.90  1.49  1.25  1.09  0.98  0.89  0.83  0.78
      47     1.436  4.47  2.35  1.72  1.36  1.14  0.99  0.89  0.81  0.75  0.71
      36     1.353  5.84  3.06  2.14  1.69  1.41  1.23  1.10  1.01  0.94  0.92
      17     1.317  5.56  2.92  2.04  1.61  1.35  1.17  1.05  0.96  0.89  0.84
     619   overall  2.38  1.36  1.02  0.89  0.79  0.74  0.71  0.66  0.64  0.63


 List of    116 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   48    6   14    1.39  200.81  0.3066E+02  0.1301E+02 "alien"
   39    8   19    1.41  114.37  0.1746E+02  0.1226E+02 "alien"
   28   21   22    1.38  111.12  0.1696E+02  0.8459E+01 "alien"
   41   33    7    1.42  104.60  0.3212E+02  0.1466E+02 "alien"
   43    9   19    1.34  102.61  0.3524E+02  0.1491E+02 "alien"
   52    2    4    1.47  102.47  0.4305E+02  0.1703E+02 "alien"
   38    9   20    1.39   90.43  0.1381E+02  0.1203E+02 "alien"
   28   17   24    1.33   84.52  0.2902E+02  0.1250E+02 "alien"
   51    4   10    1.41   73.98  0.1129E+02  0.1528E+02 "alien"
   27   18   23    1.38   72.52  0.1107E+02  0.1099E+02 "alien"
   37    4   15    1.62   62.68 -0.2527E+02  0.1209E+02 "alien"
   29   20   23    1.34   61.38  0.2108E+02  0.1275E+02 "alien"
   38   33   12    1.39   59.47  0.9080E+01  0.9186E+01 "alien"
   27   24   16    1.61   53.51 -0.2158E+02  0.1249E+02 "alien"
   29   27   19    1.41   48.38  0.7386E+01  0.9185E+01 "alien"
   34   33   16    1.36   48.32  0.1659E+02  0.1040E+02 "alien"
   46    1    6    1.63   47.99 -0.1935E+02  0.1191E+02 "alien"
   49    7   15    1.34   45.85  0.1574E+02  0.1516E+02 "alien"
   27   24   20    1.44   43.45  0.1335E+02  0.1235E+02 "alien"
   49    3    9    1.48   42.52  0.1787E+02  0.8697E+01 "alien"
   33   30   15    1.44   42.09  0.1293E+02  0.8634E+01 "alien"
   36   10   21    1.38   41.91  0.6398E+01  0.8026E+01 "alien"
   43    5   15    1.47   41.59  0.1748E+02  0.9566E+01 "alien"
   52    5   10    1.39   40.14  0.6128E+01  0.1005E+02 "alien"
   38   34   13    1.35   38.87  0.1335E+02  0.9242E+01 "alien"
   37   30   11    1.48   38.49  0.1617E+02  0.1059E+02 "alien"
   45    1    8    1.62   38.29 -0.1544E+02  0.9426E+01 "alien"
   30   25   21    1.35   37.95  0.1303E+02  0.1226E+02 "alien"
   30   28   16    1.49   36.21  0.1521E+02  0.1371E+02 "alien"
   38    8   19    1.43   36.05  0.1107E+02  0.8221E+01 "alien"
   35   10   21    1.40   35.97  0.5492E+01  0.1168E+02 "alien"
   27   13   24    1.37   34.36  0.1180E+02  0.1070E+02 "alien"
   31    9   21    1.46   34.21  0.1437E+02  0.1175E+02 "alien"
   54    3    1    1.43   33.21  0.1020E+02  0.9527E+01 "alien"
   47    1    5    1.61   32.50 -0.1311E+02  0.8651E+01 "alien"
   55    3    0    1.41   31.71  0.5699E+02  0.1294E+02 "alien"
   29   28   17    1.47   31.27  0.1314E+02  0.8344E+01 "alien"
   45    2   10    1.57   28.98 -0.1168E+02  0.9287E+01 "alien"
   49    4   11    1.44   28.36  0.8710E+01  0.9241E+01 "alien"
   36   31   13    1.43   28.17  0.8652E+01  0.9478E+01 "alien"
   34   27   11    1.59   28.10 -0.1133E+02  0.1351E+02 "alien"
   31   25   13    1.63   27.67 -0.1116E+02  0.7660E+01 "alien"
   33   29   14    1.49   27.59  0.1159E+02  0.1421E+02 "alien"
   41   32    5    1.46   27.47  0.1154E+02  0.8867E+01 "alien"
   42    8   18    1.39   27.25  0.4161E+01  0.9405E+01 "alien"
   39   31    8    1.48   26.72  0.1123E+02  0.9442E+01 "alien"
   45    5   14    1.46   26.12  0.8023E+01  0.8491E+01 "alien"
   29    9   21    1.50   25.97  0.1091E+02  0.7218E+01 "alien"
   40    3   13    1.62   25.95 -0.1046E+02  0.8454E+01 "alien"
   36    3   14    1.70   25.64  0.2360E+02  0.1263E+02 "alien"
   32   30   16    1.43   25.62  0.7869E+01  0.8396E+01 "alien"
   30    9   21    1.48   25.39  0.1067E+02  0.9910E+01 "alien"
   39    6   17    1.49   24.91  0.1047E+02  0.8028E+01 "alien"
   25   20   21    1.47   24.28  0.1020E+02  0.1076E+02 "alien"
   36   24    2    1.78   23.26  0.1534E+02  0.8112E+01 "alien"
   41   10   20    1.33   23.22  0.7974E+01  0.9039E+01 "alien"
   52    3    6    1.45   23.07  0.7085E+01  0.8756E+01 "alien"
   53    3    4    1.44   23.03  0.7073E+01  0.9594E+01 "alien"
   48    9   17    1.30   22.68  0.2559E+02  0.1653E+02 "alien"
   34    7   19    1.51   21.97  0.9231E+01  0.1113E+02 "alien"
   42    9   19    1.35   21.74  0.7465E+01  0.8771E+01 "alien"
   49    5   12    1.42   21.68  0.6657E+01  0.9204E+01 "alien"
   38    4   15    1.60   21.52 -0.8679E+01  0.7374E+01 "alien"
   30   28   19    1.39   21.00  0.3205E+01  0.8552E+01 "alien"
   33   32   18    1.33   20.43  0.2306E+02  0.1258E+02 "alien"
   40   11   21    1.31   20.19  0.2278E+02  0.1352E+02 "alien"
   39   34   12    1.36   19.29  0.6622E+01  0.1119E+02 "alien"
   41    3   13    1.59   19.24 -0.7757E+01  0.1159E+02 "alien"
   30   25   14    1.61   19.12 -0.7710E+01  0.7485E+01 "alien"
   27   20   22    1.40   18.46  0.2818E+01  0.7806E+01 "alien"
   51    2    5    1.49   18.44  0.7749E+01  0.9162E+01 "alien"
   39   31    9    1.46   18.40  0.7730E+01  0.9904E+01 "alien"
   32   25   12    1.64   18.28  0.1682E+02  0.8233E+01 "alien"
   53    2    1    1.46   18.09  0.7600E+01  0.1141E+02 "alien"
   33    8   20    1.48   17.98  0.7556E+01  0.1164E+02 "alien"
   50    4   10    1.44   17.88  0.5491E+01  0.8856E+01 "alien"
   40    2   12    1.66   17.88  0.1645E+02  0.7958E+01 "alien"
   38   28    5    1.61   17.85 -0.7196E+01  0.1029E+02 "alien"
   54    3    2    1.43   17.35  0.5328E+01  0.9020E+01 "alien"
   28    9   21    1.51   17.24 -0.2397E+02  0.9386E+01 "alien"
   39    2   12    1.69   17.12  0.1575E+02  0.7795E+01 "alien"
   50    2    5    1.52   16.92 -0.2352E+02  0.8771E+01 "alien"
   36   24    3    1.77   16.49  0.1088E+02  0.9332E+01 "alien"
   37    2   13    1.71   15.95  0.1052E+02  0.8642E+01 "alien"
   51    2    3    1.51   15.50  0.6511E+01  0.9249E+01 "alien"
   35    1   12    1.83   15.46 -0.1704E+02  0.7964E+01 "alien"
   26   19   22    1.42   15.31  0.4701E+01  0.7891E+01 "alien"
   45    9   18    1.33   15.24  0.1720E+02  0.1389E+02 "alien"
   21   17   20    1.61   15.20 -0.6130E+01  0.8641E+01 "alien"
   40    1    6    1.86   14.97 -0.1650E+02  0.8474E+01 "alien"
   35   25    7    1.71   14.90  0.9825E+01  0.8032E+01 "alien"
   50    2    6    1.51   14.87  0.6247E+01  0.8820E+01 "alien"
   44    1    0    1.76   14.42  0.2582E+02  0.1313E+02 "alien"
   53    3    5    1.44   14.19  0.4359E+01  0.8739E+01 "alien"
   34    5   17    1.60   14.13 -0.5697E+01  0.7067E+01 "alien"
   45    4   13    1.49   13.82  0.5806E+01  0.9092E+01 "alien"
   36   34   15    1.34   13.60  0.4669E+01  0.1299E+02 "alien"
   37   30   10    1.50   12.93  0.5432E+01  0.1478E+02 "alien"
   21   17   14    1.99   12.76 -0.1845E+02  0.1103E+02 "alien"
   34    4   16    1.66   12.58  0.1158E+02  0.7520E+01 "alien"
   36    2   13    1.74   12.34  0.8138E+01  0.7087E+01 "alien"
   28   27   17    1.50   12.33  0.5178E+01  0.8943E+01 "alien"
   23   22   17    1.67   12.25  0.1128E+02  0.9953E+01 "alien"
   35   24    5    1.78   12.15  0.8009E+01  0.8055E+01 "alien"
   39   29    4    1.57   11.97 -0.4828E+01  0.1476E+02 "alien"
   41    1    5    1.84   11.82 -0.1303E+02  0.1116E+02 "alien"
   36   25    5    1.72   11.81  0.7788E+01  0.7707E+01 "alien"
   35   23    2    1.84   11.52 -0.1269E+02  0.7839E+01 "alien"
   54    4    6    1.40   11.39  0.1738E+01  0.9747E+01 "alien"
   39   32   10    1.43   11.32  0.3478E+01  0.9532E+01 "alien"
   39    1   11    1.73   11.07  0.7297E+01  0.1341E+02 "alien"
   52    5   11    1.37   11.02  0.3783E+01  0.1527E+02 "alien"
   27    9   21    1.53   10.95 -0.1523E+02  0.9392E+01 "alien"
   29    2   14    1.92   10.80 -0.1562E+02  0.1009E+02 "alien"
   41    0    7    1.79   10.46 -0.2395E+02  0.1340E+02 "alien"
   30    8   20    1.53   10.46 -0.1454E+02  0.1022E+02 "alien"


 List of    116 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   48    9   17    1.30   22.68  0.2559E+02  0.1653E+02
   40   11   21    1.31   20.19  0.2278E+02  0.1352E+02
   33   32   18    1.33   20.43  0.2306E+02  0.1258E+02
   45    9   18    1.33   15.24  0.1720E+02  0.1389E+02
   41   10   20    1.33   23.22  0.7974E+01  0.9039E+01
   28   17   24    1.33   84.52  0.2902E+02  0.1250E+02
   43    9   19    1.34  102.61  0.3524E+02  0.1491E+02
   29   20   23    1.34   61.38  0.2108E+02  0.1275E+02
   36   34   15    1.34   13.60  0.4669E+01  0.1299E+02
   49    7   15    1.34   45.85  0.1574E+02  0.1516E+02
   42    9   19    1.35   21.74  0.7465E+01  0.8771E+01
   30   25   21    1.35   37.95  0.1303E+02  0.1226E+02
   38   34   13    1.35   38.87  0.1335E+02  0.9242E+01
   34   33   16    1.36   48.32  0.1659E+02  0.1040E+02
   39   34   12    1.36   19.29  0.6622E+01  0.1119E+02
   52    5   11    1.37   11.02  0.3783E+01  0.1527E+02
   27   13   24    1.37   34.36  0.1180E+02  0.1070E+02
   27   18   23    1.38   72.52  0.1107E+02  0.1099E+02
   28   21   22    1.38  111.12  0.1696E+02  0.8459E+01
   36   10   21    1.38   41.91  0.6398E+01  0.8026E+01
   52    5   10    1.39   40.14  0.6128E+01  0.1005E+02
   30   28   19    1.39   21.00  0.3205E+01  0.8552E+01
   48    6   14    1.39  200.81  0.3066E+02  0.1301E+02
   38    9   20    1.39   90.43  0.1381E+02  0.1203E+02
   42    8   18    1.39   27.25  0.4161E+01  0.9405E+01
   38   33   12    1.39   59.47  0.9080E+01  0.9186E+01
   35   10   21    1.40   35.97  0.5492E+01  0.1168E+02
   54    4    6    1.40   11.39  0.1738E+01  0.9747E+01
   27   20   22    1.40   18.46  0.2818E+01  0.7806E+01
   55    3    0    1.41   31.71  0.5699E+02  0.1294E+02
   39    8   19    1.41  114.37  0.1746E+02  0.1226E+02
   51    4   10    1.41   73.98  0.1129E+02  0.1528E+02
   29   27   19    1.41   48.38  0.7386E+01  0.9185E+01
   49    5   12    1.42   21.68  0.6657E+01  0.9204E+01
   26   19   22    1.42   15.31  0.4701E+01  0.7891E+01
   41   33    7    1.42  104.60  0.3212E+02  0.1466E+02
   38    8   19    1.43   36.05  0.1107E+02  0.8221E+01
   32   30   16    1.43   25.62  0.7869E+01  0.8396E+01
   39   32   10    1.43   11.32  0.3478E+01  0.9532E+01
   54    3    2    1.43   17.35  0.5328E+01  0.9020E+01
   36   31   13    1.43   28.17  0.8652E+01  0.9478E+01
   54    3    1    1.43   33.21  0.1020E+02  0.9527E+01
   53    3    5    1.44   14.19  0.4359E+01  0.8739E+01
   50    4   10    1.44   17.88  0.5491E+01  0.8856E+01
   27   24   20    1.44   43.45  0.1335E+02  0.1235E+02
   49    4   11    1.44   28.36  0.8710E+01  0.9241E+01
   33   30   15    1.44   42.09  0.1293E+02  0.8634E+01
   53    3    4    1.44   23.03  0.7073E+01  0.9594E+01
   52    3    6    1.45   23.07  0.7085E+01  0.8756E+01
   45    5   14    1.46   26.12  0.8023E+01  0.8491E+01
   53    2    1    1.46   18.09  0.7600E+01  0.1141E+02
   39   31    9    1.46   18.40  0.7730E+01  0.9904E+01
   41   32    5    1.46   27.47  0.1154E+02  0.8867E+01
   31    9   21    1.46   34.21  0.1437E+02  0.1175E+02
   25   20   21    1.47   24.28  0.1020E+02  0.1076E+02
   29   28   17    1.47   31.27  0.1314E+02  0.8344E+01
   43    5   15    1.47   41.59  0.1748E+02  0.9566E+01
   52    2    4    1.47  102.47  0.4305E+02  0.1703E+02
   37   30   11    1.48   38.49  0.1617E+02  0.1059E+02
   33    8   20    1.48   17.98  0.7556E+01  0.1164E+02
   30    9   21    1.48   25.39  0.1067E+02  0.9910E+01
   39   31    8    1.48   26.72  0.1123E+02  0.9442E+01
   49    3    9    1.48   42.52  0.1787E+02  0.8697E+01
   45    4   13    1.49   13.82  0.5806E+01  0.9092E+01
   33   29   14    1.49   27.59  0.1159E+02  0.1421E+02
   39    6   17    1.49   24.91  0.1047E+02  0.8028E+01
   30   28   16    1.49   36.21  0.1521E+02  0.1371E+02
   51    2    5    1.49   18.44  0.7749E+01  0.9162E+01
   29    9   21    1.50   25.97  0.1091E+02  0.7218E+01
   37   30   10    1.50   12.93  0.5432E+01  0.1478E+02
   28   27   17    1.50   12.33  0.5178E+01  0.8943E+01
   50    2    6    1.51   14.87  0.6247E+01  0.8820E+01
   34    7   19    1.51   21.97  0.9231E+01  0.1113E+02
   51    2    3    1.51   15.50  0.6511E+01  0.9249E+01
   28    9   21    1.51   17.24 -0.2397E+02  0.9386E+01
   50    2    5    1.52   16.92 -0.2352E+02  0.8771E+01
   27    9   21    1.53   10.95 -0.1523E+02  0.9392E+01
   30    8   20    1.53   10.46 -0.1454E+02  0.1022E+02
   45    2   10    1.57   28.98 -0.1168E+02  0.9287E+01
   39   29    4    1.57   11.97 -0.4828E+01  0.1476E+02
   34   27   11    1.59   28.10 -0.1133E+02  0.1351E+02
   41    3   13    1.59   19.24 -0.7757E+01  0.1159E+02
   38    4   15    1.60   21.52 -0.8679E+01  0.7374E+01
   34    5   17    1.60   14.13 -0.5697E+01  0.7067E+01
   38   28    5    1.61   17.85 -0.7196E+01  0.1029E+02
   27   24   16    1.61   53.51 -0.2158E+02  0.1249E+02
   21   17   20    1.61   15.20 -0.6130E+01  0.8641E+01
   47    1    5    1.61   32.50 -0.1311E+02  0.8651E+01
   30   25   14    1.61   19.12 -0.7710E+01  0.7485E+01
   40    3   13    1.62   25.95 -0.1046E+02  0.8454E+01
   37    4   15    1.62   62.68 -0.2527E+02  0.1209E+02
   45    1    8    1.62   38.29 -0.1544E+02  0.9426E+01
   46    1    6    1.63   47.99 -0.1935E+02  0.1191E+02
   31   25   13    1.63   27.67 -0.1116E+02  0.7660E+01
   32   25   12    1.64   18.28  0.1682E+02  0.8233E+01
   34    4   16    1.66   12.58  0.1158E+02  0.7520E+01
   40    2   12    1.66   17.88  0.1645E+02  0.7958E+01
   23   22   17    1.67   12.25  0.1128E+02  0.9953E+01
   39    2   12    1.69   17.12  0.1575E+02  0.7795E+01
   36    3   14    1.70   25.64  0.2360E+02  0.1263E+02
   37    2   13    1.71   15.95  0.1052E+02  0.8642E+01
   35   25    7    1.71   14.90  0.9825E+01  0.8032E+01
   36   25    5    1.72   11.81  0.7788E+01  0.7707E+01
   39    1   11    1.73   11.07  0.7297E+01  0.1341E+02
   36    2   13    1.74   12.34  0.8138E+01  0.7087E+01
   44    1    0    1.76   14.42  0.2582E+02  0.1313E+02
   36   24    3    1.77   16.49  0.1088E+02  0.9332E+01
   35   24    5    1.78   12.15  0.8009E+01  0.8055E+01
   36   24    2    1.78   23.26  0.1534E+02  0.8112E+01
   41    0    7    1.79   10.46 -0.2395E+02  0.1340E+02
   35    1   12    1.83   15.46 -0.1704E+02  0.7964E+01
   41    1    5    1.84   11.82 -0.1303E+02  0.1116E+02
   35   23    2    1.84   11.52 -0.1269E+02  0.7839E+01
   40    1    6    1.86   14.97 -0.1650E+02  0.8474E+01
   29    2   14    1.92   10.80 -0.1562E+02  0.1009E+02
   21   17   14    1.99   12.76 -0.1845E+02  0.1103E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        81
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1192
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1025
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      167


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
