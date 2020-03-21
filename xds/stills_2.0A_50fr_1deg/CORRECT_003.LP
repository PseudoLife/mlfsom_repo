
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.295    77.295    38.681  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=       3       3
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       3
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  1.000000 -0.000574  0.000714
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
 BEAM_DIVERGENCE_E.S.D.=   0.039
 REFLECTING_RANGE_E.S.D.=   0.072
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

    1546 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       3 ...       3
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1546 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               48
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
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    27        mC        500.0     159.4   38.7  109.3  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.7  109.3  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.3  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        998.2     172.8   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.2      77.3  172.8   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.2      77.3  172.8   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.2      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.4   86.4  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.3   86.4   86.4  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.4   86.4  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.7   86.4  247.5 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.4   86.4   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.1  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.4   86.4  177.1 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.4   86.4   86.4  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
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
 UNIT_CELL_A-AXIS=    73.667   -23.402     0.034
 UNIT_CELL_B-AXIS=    23.402    73.667     0.011
 UNIT_CELL_C-AXIS=    -0.018    -0.000    38.682
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       351 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.19
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.249    77.249    38.661  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.6E-01 1.6E-01 7.5E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012945  0.012945  0.025866  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.621   -23.394     0.035
 COORDINATES OF UNIT CELL B-AXIS    23.394    73.621     0.011
 COORDINATES OF UNIT CELL C-AXIS    -0.018    -0.000    38.661
 CRYSTAL MOSAICITY (DEGREES)     0.072
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006297  0.002236
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.84   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1260.59   1242.05
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.15
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000


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
  0.05  0.15       141        323.058       1.043
  0.15  0.25       187        319.109       1.056
  0.25  0.35       184        290.733       1.159
  0.35  0.45       158        211.412       1.593
  0.45  0.55       180        439.580       0.766
  0.55  0.65       193        287.547       1.172
  0.65  0.75       188        234.468       1.437
  0.75  0.85       167        330.957       1.018
  0.85  0.95       148        639.847       0.526
 
 -0.05  0.95      1546        336.863       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               240.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.971
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     2.1 XMAX=     2.9 NXBIN=    1
 YMIN= 0.00168 YMAX= 0.59763 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        522


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               240.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.971
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    11.6 XMAX=  2454.3 NXBIN=    3
 YMIN=     6.2 YMAX=  2520.1 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        522


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               240.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.969
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     2.1 XMAX=     2.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1263
 DETECTOR_SURFACE_POSITION=    1648    1690
 DETECTOR_SURFACE_POSITION=     818    1690
 DETECTOR_SURFACE_POSITION=     818     836
 DETECTOR_SURFACE_POSITION=    1648     836
 DETECTOR_SURFACE_POSITION=    2172    1663
 DETECTOR_SURFACE_POSITION=    1622    2229
 DETECTOR_SURFACE_POSITION=     844    2229
 DETECTOR_SURFACE_POSITION=     294    1663
 DETECTOR_SURFACE_POSITION=     294     863
 DETECTOR_SURFACE_POSITION=     844     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2172     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        522


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
 1.297E+00  3.028E-03   15.96


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

   24.428   3.394    14.50   1.16      2.48      3.60      66      95      20
    3.394   2.412    12.37   1.33      4.70      3.93     128     171      16
    2.412   1.972     8.79   1.01      5.52      4.82     178     219      12
    1.972   1.709    -0.03   0.59    -99.90    -99.90     200     251       0
    1.709   1.530    -0.05   0.71   7539.80   8928.69     220     293       0
    1.530   1.397    -0.05   0.62    -99.90    -99.90     248     318       0
    1.397   1.294     0.06   0.91    465.06    468.08     100     148       2
   --------------------------------------------------------------------------
   24.428   1.294     3.60   0.83      5.22      5.65    1140    1495      50
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       3       3)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.679   0.1069E+02  0.1263E+02         0.85           1 
   47    0    0     1.644  -0.1579E+02  0.1306E+02        -1.21           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -4.4% OF MEAN INTENSITY


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

     3.84          72      47      1285        3.7%       2.2%      3.6%       50   17.76      3.2%    99.9*     0    0.000       0
     2.73         106      68      2098        3.2%       4.7%      3.8%       76   16.39      6.6%    98.4*     0    0.000       0
     2.23         165     101      2665        3.8%       4.5%      4.1%      128   14.56      6.4%    98.7*     0    0.000       0
     1.94         168      99      3101        3.2%       6.1%      5.8%      138    7.99      8.7%    99.0*     0    0.000       0
     1.73         193     116      3481        3.3%     -99.9%    -99.9%      154    0.00    -99.9%   -17.0      0    0.000       0
     1.58         225     137      3848        3.6%    2256.1%   2758.6%      176   -0.05   3190.6%     2.7      0    0.000       0
     1.47         253     158      4143        3.8%     -99.9%    -99.9%      190    0.00    -99.9%    -4.8      0    0.000       0
     1.37         224     136      4429        3.1%     -99.9%    -99.9%      176    0.00    -99.9%     2.1      0    0.000       0
     1.29          89      63      4743        1.3%     293.9%    259.0%       52    0.07    415.6%   -25.6      0    0.000       0
    total        1495     925     29793        3.1%       5.2%      5.7%     1140    4.52      7.4%    99.7*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1546
 NUMBER OF REJECTED MISFITS                                50
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1495
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    925



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

   WILSON LINE (using all data) : A=   8.705 B=  34.553 CORRELATION=  0.91
      #      RES      SS        <I>       log(<I>)       BO
      18     7.515   0.004  3.1029E+03       8.040      75.1
      24     4.620   0.012  3.1029E+03       8.040      28.4
      32     3.507   0.020  2.2036E+03       7.698      24.8
      31     2.991   0.028  8.8654E+02       6.787      34.3
      41     2.642   0.036  9.2903E+02       6.834      26.1
      55     2.395   0.044  6.2797E+02       6.442      25.9
      40     2.207   0.051  4.3370E+02       6.072      25.6
      50     2.055   0.059  2.3409E+02       5.456      27.5
      51     1.930   0.067 -1.0922E+00    -999.999   -1000.0
      63     1.825   0.075  5.3187E-01      -0.631      62.2
      56     1.736   0.083 -1.0977E+00    -999.999   -1000.0
      60     1.662   0.090 -4.1740E-01    -999.999   -1000.0
      63     1.594   0.098 -8.5729E-01    -999.999   -1000.0
      83     1.534   0.106 -1.0552E+00    -999.999   -1000.0
      81     1.478   0.114  6.5283E-01      -0.426      39.9
      61     1.430   0.122 -1.6211E+00    -999.999   -1000.0
      55     1.390   0.129 -6.4140E-01    -999.999   -1000.0
      41     1.350   0.137  1.2257E+00       0.204      31.0
      19     1.313   0.145  7.9716E-01      -0.227      30.8


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.515        0.994        0.479        0.183
       3     4.620        1.197        0.704        0.309
       4     3.507        4.729        6.812        7.185
       4     2.991        0.910        0.369        0.113
       3     2.642        0.243        0.058        0.010
       5     2.395        0.712        0.313        0.103
       3     2.207        1.176        0.449        0.124
       8     2.055        0.233        0.074        0.018
      11     1.662       31.958      -55.596      312.921
       3     1.534       43.386      158.114      426.368
       6     1.478       49.352     -100.364      641.079
       3     1.430        9.494        5.159       12.394
       1     1.390      186.195      880.125     2971.606
       2     1.350       10.477       12.062       10.077
      94   overall       13.611       -0.677      160.719


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      13     7.515        1.281        1.955        2.455
      21     4.620        1.096        1.248        1.238
      28     3.507        0.573        0.371        0.214
      27     2.991        0.806        0.607        0.399
      38     2.642        1.135        1.202        1.136
      50     2.395        0.964        0.818        0.590
      37     2.207        0.823        0.647        0.438
      42     2.055        1.299        1.440        1.461
      52     1.825       81.735      -84.626     4351.287
      75     1.478       71.217       88.169     2355.874
      39     1.350       40.688       17.576      928.492
      19     1.313      191.900     2921.321    59602.030
     830   overall       18.167       70.665     1893.774

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.515  0.00  0.00  0.00  0.00  0.00  0.36  0.33  0.64  0.61  0.59
       3     4.620  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       4     3.507  1.01  0.72  0.60  0.53  0.48  0.45  0.42  0.40  0.76  0.73
       4     2.991  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.642  1.34  0.96  0.80  0.70  0.64  0.59  1.12  1.06  1.01  0.98
       5     2.395  0.00  0.00  0.00  0.85  0.77  1.07  1.00  0.95  0.91  0.88
       3     2.207  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.055  1.51  1.45  1.50  1.32  1.20  1.11  1.05  1.19  1.14  1.10
      11     1.662  2.93  2.10  1.75  1.54  1.40  1.30  1.22  1.16  1.11  1.07
       3     1.534  2.68  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       6     1.478  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       3     1.430  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       1     1.390  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       2     1.350  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      94   overall  0.94  0.74  0.64  0.61  0.55  0.57  0.55  0.56  0.55  0.55


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      13     7.515  1.62  2.12  1.78  2.33  1.95  1.70  1.53  1.40  1.30  1.22
      21     4.620  2.00  1.84  1.29  1.16  1.21  1.16  1.14  1.04  1.12  1.05
      28     3.507  1.13  0.79  0.83  0.87  1.09  1.11  1.06  1.23  1.14  1.13
      27     2.991  1.56  1.43  1.14  1.12  0.94  0.99  1.10  1.08  1.06  1.00
      38     2.642  0.83  1.45  1.12  1.12  1.00  0.87  0.99  1.10  1.06  1.04
      50     2.395  1.47  1.32  1.08  1.21  1.07  1.11  1.07  1.09  1.08  1.01
      37     2.207  1.42  1.34  1.15  0.90  0.96  1.02  0.97  1.03  1.05  1.07
      42     2.055  2.00  1.31  1.19  1.16  1.03  0.90  0.95  0.95  0.96  0.94
      52     1.825  4.45  2.33  1.63  1.34  1.12  1.02  0.92  0.84  0.78  0.73
      75     1.478  5.04  2.65  1.85  1.46  1.22  1.06  0.95  0.87  0.83  0.78
      39     1.350  4.85  2.55  1.88  1.56  1.30  1.14  1.02  0.93  0.86  0.81
      19     1.313  6.08  3.19  2.23  1.76  1.47  1.28  1.15  1.05  0.98  0.92
     830   overall  1.56  1.00  0.76  0.68  0.61  0.57  0.54  0.53  0.52  0.50


 List of     94 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   71.50  0.5700E+02  0.1497E+02 "alien"
   48    6   14    1.38   38.66 -0.2480E+02  0.1315E+02 "alien"
   35   25    8    1.68   37.06 -0.1547E+02  0.1316E+02 "alien"
   34   25    9    1.68   34.66 -0.1447E+02  0.8731E+01 "alien"
   40   33    9    1.41   32.37 -0.2076E+02  0.1033E+02 "alien"
   53    2    1    1.46   31.69  0.2069E+02  0.1247E+02 "alien"
   52    2    4    1.47   30.86  0.2015E+02  0.1861E+02 "alien"
   38    9   20    1.38   30.02 -0.1925E+02  0.1202E+02 "alien"
   30   21   10    1.85   27.75  0.1476E+02  0.1252E+02 "alien"
   39   30    6    1.53   27.35 -0.2886E+02  0.1051E+02 "alien"
   32   24   11    1.69   26.68 -0.1114E+02  0.8930E+01 "alien"
   21   19   19    1.63   25.04 -0.1045E+02  0.9308E+01 "alien"
   34    4   16    1.65   24.78 -0.1034E+02  0.7855E+01 "alien"
   38    3   14    1.63   24.65 -0.1029E+02  0.8128E+01 "alien"
   33   22    6    1.86   24.16  0.1285E+02  0.9489E+01 "alien"
   36    1   11    1.83   24.11  0.1282E+02  0.1223E+02 "alien"
   33   23    8    1.78   22.82  0.1214E+02  0.9253E+01 "alien"
   33   25   11    1.65   22.71 -0.9478E+01  0.1276E+02 "alien"
   42    5   15    1.49   22.11  0.1443E+02  0.8851E+01 "alien"
   54    5    8    1.37   21.60  0.2647E+02  0.1153E+02 "alien"
   37    5   16    1.57   21.59 -0.1851E+02  0.1162E+02 "alien"
   31    5   17    1.67   20.79 -0.8679E+01  0.1104E+02 "alien"
   36   28    9    1.58   20.35 -0.1744E+02  0.1594E+02 "alien"
   33   32   16    1.38   20.19 -0.1295E+02  0.1569E+02 "alien"
   42    1    1    1.84   19.21  0.1022E+02  0.9199E+01 "alien"
   49    6   14    1.36   18.88  0.2314E+02  0.1699E+02 "alien"
   39   31    9    1.46   18.66  0.1218E+02  0.1104E+02 "alien"
   31   30   17    1.41   18.19 -0.1166E+02  0.1029E+02 "alien"
   49    2    6    1.53   17.29 -0.1824E+02  0.1410E+02 "alien"
   41    9   19    1.36   17.23  0.2111E+02  0.1021E+02 "alien"
   54    4    7    1.38   17.22 -0.1105E+02  0.1003E+02 "alien"
   44    0    6    1.69   17.00 -0.2333E+02  0.8790E+01 "alien"
   36   32   13    1.41   16.99 -0.2754E+02  0.1497E+02 "alien"
   33    4   16    1.67   16.51 -0.6891E+01  0.7096E+01 "alien"
   39    3   13    1.65   16.42 -0.6854E+01  0.1487E+02 "alien"
   30   22   11    1.79   16.36  0.8703E+01  0.8530E+01 "alien"
   22   21   17    1.69   16.36 -0.6827E+01  0.9625E+01 "alien"
   42   10   20    1.31   16.25  0.1296E+02  0.1339E+02 "alien"
   37   30   11    1.47   16.17  0.1055E+02  0.1408E+02 "alien"
   26   10   22    1.48   16.15  0.1054E+02  0.9177E+01 "alien"
   36   32   14    1.39   16.10 -0.1032E+02  0.1694E+02 "alien"
   32   28   14    1.52   15.96 -0.1684E+02  0.1483E+02 "alien"
   38   30    9    1.50   15.45  0.1009E+02  0.9826E+01 "alien"
   36   33   14    1.37   15.38 -0.9863E+01  0.9923E+01 "alien"
   46    2    8    1.58   15.11 -0.1296E+02  0.1654E+02 "alien"
   30    7   19    1.58   15.01 -0.1287E+02  0.9985E+01 "alien"
   28   19   10    1.97   14.96 -0.1634E+02  0.1558E+02 "alien"
   24   19   14    1.86   14.71  0.7822E+01  0.1106E+02 "alien"
   37   26    4    1.68   14.38 -0.6003E+01  0.9235E+01 "alien"
   37    4   15    1.62   14.27 -0.1223E+02  0.8198E+01 "alien"
   36   29   11    1.51   14.03 -0.1481E+02  0.9370E+01 "alien"
   23   15   22    1.49   14.01  0.9146E+01  0.6790E+01 "alien"
   52    5   10    1.38   13.88 -0.8903E+01  0.9859E+01 "alien"
   54    4    6    1.39   13.69 -0.8781E+01  0.1011E+02 "alien"
   30    9   21    1.48   13.60  0.8879E+01  0.9896E+01 "alien"
   39   30    7    1.51   13.20 -0.1393E+02  0.9945E+01 "alien"
   46    6   15    1.40   13.19 -0.8460E+01  0.9385E+01 "alien"
   49    3    8    1.50   13.06  0.8526E+01  0.1499E+02 "alien"
   53    4    8    1.39   12.97 -0.8322E+01  0.9893E+01 "alien"
   22   10   22    1.54   12.93 -0.1364E+02  0.7211E+01 "alien"
   30   28   19    1.38   12.87 -0.8252E+01  0.9096E+01 "alien"
   30   26   20    1.37   12.82 -0.8222E+01  0.1359E+02 "alien"
   27    2   14    1.98   12.79 -0.1397E+02  0.1126E+02 "alien"
   34    3   15    1.70   12.62 -0.1386E+02  0.1188E+02 "alien"
   37   26    3    1.69   12.62 -0.5267E+01  0.9395E+01 "alien"
   30   28   16    1.48   12.42  0.8105E+01  0.9528E+01 "alien"
   36   30   12    1.47   12.33  0.8050E+01  0.9808E+01 "alien"
   29   28   17    1.47   12.21  0.7971E+01  0.8943E+01 "alien"
   43    2   11    1.60   12.19 -0.1045E+02  0.1239E+02 "alien"
   38   28    6    1.59   11.99 -0.1028E+02  0.1374E+02 "alien"
   42    3   12    1.59   11.90 -0.1020E+02  0.1282E+02 "alien"
   29    9   21    1.49   11.86  0.7741E+01  0.6912E+01 "alien"
   35   23    0    1.84   11.80  0.1422E+02  0.1334E+02 "alien"
   38   31   10    1.46   11.69  0.7633E+01  0.9876E+01 "alien"
   38   29    7    1.55   11.61 -0.1225E+02  0.9922E+01 "alien"
   43    5   15    1.47   11.60  0.7575E+01  0.1025E+02 "alien"
   24   16   22    1.47   11.39  0.7433E+01  0.1023E+02 "alien"
   42    2   11    1.63   11.36 -0.4743E+01  0.8579E+01 "alien"
   39   28    4    1.59   11.34 -0.9717E+01  0.1189E+02 "alien"
   17    8   21    1.68   11.20 -0.4675E+01  0.6113E+01 "alien"
   40    1    4    1.89   11.11 -0.1213E+02  0.1567E+02 "alien"
   43    1    9    1.66   11.09 -0.4628E+01  0.8485E+01 "alien"
   30   26   15    1.55   10.98 -0.1159E+02  0.1205E+02 "alien"
   35   25    7    1.71   10.98 -0.1205E+02  0.9002E+01 "alien"
   46    0    2    1.67   10.91 -0.9728E+01  0.8892E+01 "alien"
   44    1    1    1.75   10.78 -0.1183E+02  0.9060E+01 "alien"
   36    1    8    1.96   10.69 -0.1167E+02  0.8766E+01 "alien"
   51    2    4    1.50   10.59  0.6911E+01  0.9768E+01 "alien"
   42    4   14    1.53   10.48 -0.1106E+02  0.8658E+01 "alien"
   39    0    9    1.80   10.33  0.1132E+02  0.8651E+01 "alien"
   28   27   17    1.50   10.33  0.6741E+01  0.8586E+01 "alien"
   41    8   18    1.40   10.32 -0.6619E+01  0.1356E+02 "alien"
   35   23    2    1.84   10.30  0.5476E+01  0.9382E+01 "alien"
   51    2    3    1.50   10.09  0.6585E+01  0.9983E+01 "alien"


 List of     94 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   71.50  0.5700E+02  0.1497E+02
   42   10   20    1.31   16.25  0.1296E+02  0.1339E+02
   49    6   14    1.36   18.88  0.2314E+02  0.1699E+02
   41    9   19    1.36   17.23  0.2111E+02  0.1021E+02
   54    5    8    1.37   21.60  0.2647E+02  0.1153E+02
   30   26   20    1.37   12.82 -0.8222E+01  0.1359E+02
   36   33   14    1.37   15.38 -0.9863E+01  0.9923E+01
   33   32   16    1.38   20.19 -0.1295E+02  0.1569E+02
   52    5   10    1.38   13.88 -0.8903E+01  0.9859E+01
   54    4    7    1.38   17.22 -0.1105E+02  0.1003E+02
   30   28   19    1.38   12.87 -0.8252E+01  0.9096E+01
   48    6   14    1.38   38.66 -0.2480E+02  0.1315E+02
   38    9   20    1.38   30.02 -0.1925E+02  0.1202E+02
   36   32   14    1.39   16.10 -0.1032E+02  0.1694E+02
   53    4    8    1.39   12.97 -0.8322E+01  0.9893E+01
   54    4    6    1.39   13.69 -0.8781E+01  0.1011E+02
   46    6   15    1.40   13.19 -0.8460E+01  0.9385E+01
   41    8   18    1.40   10.32 -0.6619E+01  0.1356E+02
   31   30   17    1.41   18.19 -0.1166E+02  0.1029E+02
   40   33    9    1.41   32.37 -0.2076E+02  0.1033E+02
   36   32   13    1.41   16.99 -0.2754E+02  0.1497E+02
   53    2    1    1.46   31.69  0.2069E+02  0.1247E+02
   39   31    9    1.46   18.66  0.1218E+02  0.1104E+02
   38   31   10    1.46   11.69  0.7633E+01  0.9876E+01
   29   28   17    1.47   12.21  0.7971E+01  0.8943E+01
   43    5   15    1.47   11.60  0.7575E+01  0.1025E+02
   52    2    4    1.47   30.86  0.2015E+02  0.1861E+02
   36   30   12    1.47   12.33  0.8050E+01  0.9808E+01
   24   16   22    1.47   11.39  0.7433E+01  0.1023E+02
   37   30   11    1.47   16.17  0.1055E+02  0.1408E+02
   30    9   21    1.48   13.60  0.8879E+01  0.9896E+01
   26   10   22    1.48   16.15  0.1054E+02  0.9177E+01
   30   28   16    1.48   12.42  0.8105E+01  0.9528E+01
   42    5   15    1.49   22.11  0.1443E+02  0.8851E+01
   23   15   22    1.49   14.01  0.9146E+01  0.6790E+01
   29    9   21    1.49   11.86  0.7741E+01  0.6912E+01
   51    2    4    1.50   10.59  0.6911E+01  0.9768E+01
   38   30    9    1.50   15.45  0.1009E+02  0.9826E+01
   28   27   17    1.50   10.33  0.6741E+01  0.8586E+01
   49    3    8    1.50   13.06  0.8526E+01  0.1499E+02
   51    2    3    1.50   10.09  0.6585E+01  0.9983E+01
   36   29   11    1.51   14.03 -0.1481E+02  0.9370E+01
   39   30    7    1.51   13.20 -0.1393E+02  0.9945E+01
   32   28   14    1.52   15.96 -0.1684E+02  0.1483E+02
   39   30    6    1.53   27.35 -0.2886E+02  0.1051E+02
   42    4   14    1.53   10.48 -0.1106E+02  0.8658E+01
   49    2    6    1.53   17.29 -0.1824E+02  0.1410E+02
   22   10   22    1.54   12.93 -0.1364E+02  0.7211E+01
   38   29    7    1.55   11.61 -0.1225E+02  0.9922E+01
   30   26   15    1.55   10.98 -0.1159E+02  0.1205E+02
   37    5   16    1.57   21.59 -0.1851E+02  0.1162E+02
   36   28    9    1.58   20.35 -0.1744E+02  0.1594E+02
   30    7   19    1.58   15.01 -0.1287E+02  0.9985E+01
   46    2    8    1.58   15.11 -0.1296E+02  0.1654E+02
   38   28    6    1.59   11.99 -0.1028E+02  0.1374E+02
   39   28    4    1.59   11.34 -0.9717E+01  0.1189E+02
   42    3   12    1.59   11.90 -0.1020E+02  0.1282E+02
   43    2   11    1.60   12.19 -0.1045E+02  0.1239E+02
   37    4   15    1.62   14.27 -0.1223E+02  0.8198E+01
   42    2   11    1.63   11.36 -0.4743E+01  0.8579E+01
   21   19   19    1.63   25.04 -0.1045E+02  0.9308E+01
   38    3   14    1.63   24.65 -0.1029E+02  0.8128E+01
   39    3   13    1.65   16.42 -0.6854E+01  0.1487E+02
   33   25   11    1.65   22.71 -0.9478E+01  0.1276E+02
   34    4   16    1.65   24.78 -0.1034E+02  0.7855E+01
   43    1    9    1.66   11.09 -0.4628E+01  0.8485E+01
   31    5   17    1.67   20.79 -0.8679E+01  0.1104E+02
   46    0    2    1.67   10.91 -0.9728E+01  0.8892E+01
   33    4   16    1.67   16.51 -0.6891E+01  0.7096E+01
   17    8   21    1.68   11.20 -0.4675E+01  0.6113E+01
   37   26    4    1.68   14.38 -0.6003E+01  0.9235E+01
   35   25    8    1.68   37.06 -0.1547E+02  0.1316E+02
   34   25    9    1.68   34.66 -0.1447E+02  0.8731E+01
   32   24   11    1.69   26.68 -0.1114E+02  0.8930E+01
   37   26    3    1.69   12.62 -0.5267E+01  0.9395E+01
   44    0    6    1.69   17.00 -0.2333E+02  0.8790E+01
   22   21   17    1.69   16.36 -0.6827E+01  0.9625E+01
   34    3   15    1.70   12.62 -0.1386E+02  0.1188E+02
   35   25    7    1.71   10.98 -0.1205E+02  0.9002E+01
   44    1    1    1.75   10.78 -0.1183E+02  0.9060E+01
   33   23    8    1.78   22.82  0.1214E+02  0.9253E+01
   30   22   11    1.79   16.36  0.8703E+01  0.8530E+01
   39    0    9    1.80   10.33  0.1132E+02  0.8651E+01
   36    1   11    1.83   24.11  0.1282E+02  0.1223E+02
   35   23    2    1.84   10.30  0.5476E+01  0.9382E+01
   42    1    1    1.84   19.21  0.1022E+02  0.9199E+01
   35   23    0    1.84   11.80  0.1422E+02  0.1334E+02
   30   21   10    1.85   27.75  0.1476E+02  0.1252E+02
   24   19   14    1.86   14.71  0.7822E+01  0.1106E+02
   33   22    6    1.86   24.16  0.1285E+02  0.9489E+01
   40    1    4    1.89   11.11 -0.1213E+02  0.1567E+02
   36    1    8    1.96   10.69 -0.1167E+02  0.8766E+01
   28   19   10    1.97   14.96 -0.1634E+02  0.1558E+02
   27    2   14    1.98   12.79 -0.1397E+02  0.1126E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        25
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1546
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1460
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       86


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.2 sec
