
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.358    77.358    38.714  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      11      11
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      11
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000743  0.000702
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
 BEAM_DIVERGENCE_E.S.D.=   0.041
 REFLECTING_RANGE_E.S.D.=   0.091
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

    1484 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      11 ...      11
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1484 REFLECTIONS ACCEPTED



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

 *  44        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  31        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  35        mP          0.0      77.4   38.7   77.4  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  33        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  21        tP          0.0      77.4   77.4   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
 *  32        oP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    37        mC        250.0     159.5   38.7   77.4  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.5   38.7   77.4  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.7  159.5   77.4  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.7  159.5   77.4  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.7  159.5   77.4  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.7  159.5   77.4  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    19        oI        500.0      38.7  109.4  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    27        mC        500.0     159.5   38.7  109.4  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    18        tI        625.0     109.4  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    26        oF        625.0      38.7  159.5  159.5  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    41        mC        998.2     173.0   77.4   38.7  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    22        hP        998.2      77.4   77.4   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    30        mC        998.2      77.4  173.0   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.2      77.4  173.0   38.7  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.4   86.5   86.5  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.5  247.7 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.5   86.5   77.4  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.5   86.5   77.4  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.7   77.4  177.3  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.5   86.5  177.3 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.5   86.5   77.4  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.3  109.4   38.7  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.7   77.4  177.3 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.5  177.3   77.4 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    6.8      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.728   -23.417     0.037
 UNIT_CELL_B-AXIS=    23.417    73.728     0.014
 UNIT_CELL_C-AXIS=    -0.020    -0.001    38.714
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       317 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.24
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.306    77.306    38.689  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  7.4E-02 7.4E-02 4.1E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012936  0.012936  0.025847  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.679   -23.402     0.035
 COORDINATES OF UNIT CELL B-AXIS    23.402    73.679     0.012
 COORDINATES OF UNIT CELL C-AXIS    -0.018    -0.000    38.689
 CRYSTAL MOSAICITY (DEGREES)     0.091
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999976  0.006619  0.002034
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.88   1232.91
 DETECTOR ORIGIN (PIXELS) AT                     1260.92   1242.54
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.11
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.31    77.31    38.69  90.000  90.000  90.000


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
  0.05  0.15       119        201.326       1.126
  0.15  0.25       170        209.861       1.080
  0.25  0.35       198        219.987       1.030
  0.35  0.45       142        133.265       1.701
  0.45  0.55       173        315.960       0.717
  0.55  0.65       198        209.817       1.080
  0.65  0.75       186        204.326       1.109
  0.75  0.85       170        212.628       1.066
  0.85  0.95       128        342.450       0.662
 
 -0.05  0.95      1484        226.621       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               240.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.993
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    10.1 XMAX=    10.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59675 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               240.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.990
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     7.3 XMAX=  2437.5 NXBIN=    3
 YMIN=     7.6 YMAX=  2518.3 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               240.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.983
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    10.1 XMAX=    10.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1222    1263
 DETECTOR_SURFACE_POSITION=    1635    1689
 DETECTOR_SURFACE_POSITION=     810    1689
 DETECTOR_SURFACE_POSITION=     810     837
 DETECTOR_SURFACE_POSITION=    1635     837
 DETECTOR_SURFACE_POSITION=    2156    1663
 DETECTOR_SURFACE_POSITION=    1609    2228
 DETECTOR_SURFACE_POSITION=     835    2228
 DETECTOR_SURFACE_POSITION=     288    1663
 DETECTOR_SURFACE_POSITION=     288     863
 DETECTOR_SURFACE_POSITION=     835     298
 DETECTOR_SURFACE_POSITION=    1609     298
 DETECTOR_SURFACE_POSITION=    2156     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        510


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
 1.105E+00  2.427E-03   19.31


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

   24.446   3.396    16.84   1.52      2.87      3.02      66      89      20
    3.396   2.413    12.89   1.45      4.59      3.64     130     171      10
    2.413   1.974     7.60   0.92      6.01      5.77     188     220       0
    1.974   1.711     0.04   0.72   1187.47   1381.84     190     241       0
    1.711   1.531     0.04   0.61  18853.94  23555.73     224     281       0
    1.531   1.398    -0.00   0.60    -99.90    -99.90     242     306       0
    1.398   1.295     0.00   0.51    -99.90    -99.90      96     143       2
   --------------------------------------------------------------------------
   24.446   1.295     3.72   0.82      5.73      5.95    1136    1451      32
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      11      11)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681  -0.1014E+02  0.1240E+02        -0.82           1 
   47    0    0     1.645   0.1043E+02  0.1258E+02         0.83           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     4.5% OF MEAN INTENSITY


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

     3.84          66      42      1285        3.3%       2.9%      3.0%       48   21.28      4.1%    99.8*     0    0.000       0
     2.73         108      68      2098        3.2%       3.5%      3.4%       80   18.04      5.0%    99.4*     0    0.000       0
     2.24         165      99      2665        3.7%       5.2%      4.1%      132   14.38      7.4%    98.9*     0    0.000       0
     1.94         166      94      3101        3.0%       7.8%      8.1%      144    6.52     11.0%    98.9*     0    0.000       0
     1.73         185     112      3481        3.2%    1104.6%   1297.3%      146    0.06   1562.1%   -16.6      0    0.000       0
     1.58         218     129      3848        3.4%    2052.4%   2653.5%      178    0.06   2902.5%    -0.4      0    0.000       0
     1.47         238     147      4143        3.5%     -99.9%    -99.9%      182    0.00    -99.9%    -1.5      0    0.000       0
     1.37         219     130      4429        2.9%    6184.8%   7445.4%      178    0.00   8746.7%    -6.6      0    0.000       0
     1.29          86      62      4743        1.3%     -99.9%    -99.9%       48    0.00    -99.9%    43.0      0    0.000       0
    total        1451     883     29793        3.0%       5.7%      6.0%     1136    4.73      8.1%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1484
 NUMBER OF REJECTED MISFITS                                32
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1451
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    883



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

   WILSON LINE (using all data) : A=   7.838 B=  32.896 CORRELATION=  0.84
      #      RES      SS        <I>       log(<I>)       BO
      17     7.250   0.005  2.1592E+03       7.678      16.9
      21     4.621   0.012  2.7601E+03       7.923      -3.6
      30     3.494   0.020  1.2281E+03       7.113      17.7
      32     2.987   0.028  5.8520E+02       6.372      26.1
      41     2.646   0.036  5.9157E+02       6.383      20.4
      53     2.394   0.044  3.3740E+02       5.821      23.1
      43     2.204   0.051  2.2476E+02       5.415      23.5
      45     2.055   0.059  9.0894E+01       4.510      28.1
      47     1.929   0.067  5.4090E-01      -0.615      62.9
      61     1.827   0.075 -2.0734E-02    -999.999   -1000.0
      53     1.739   0.083  1.6958E+00       0.528      44.2
      56     1.664   0.090  7.3286E-01      -0.311      45.1
      60     1.595   0.098 -3.4453E-01    -999.999   -1000.0
      76     1.535   0.106  1.0881E+00       0.084      36.5
      78     1.480   0.114 -1.4698E-02    -999.999   -1000.0
      58     1.431   0.122  8.5026E-01      -0.162      32.8
      52     1.390   0.129 -1.0783E+00    -999.999   -1000.0
      42     1.350   0.137 -1.8692E+00    -999.999   -1000.0
      18     1.314   0.145  6.5228E+00       1.875      20.6


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.250        1.828        1.196        0.618
       3     4.621        1.211        0.721        0.320
       3     3.494        0.337        0.083        0.015
       4     2.987        1.035        0.443        0.144
       3     2.646        0.257        0.060        0.010
       5     2.394        0.789        0.370        0.131
       3     2.204        1.257        0.492        0.138
       8     2.055        0.280        0.101        0.028
       7     1.929       50.290      193.195      666.114
       9     1.739        2.193        1.684        1.799
       3     1.535       36.277      117.688      298.761
       6     1.480       35.002       94.035      202.626
       3     1.431       75.451      332.984     1125.356
      89   overall       10.637       25.965      115.099


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.250        0.443        0.304        0.168
      18     4.621        1.091        1.229        1.174
      27     3.494        0.855        0.663        0.478
      28     2.987        0.916        0.862        0.732
      38     2.646        1.125        1.199        1.151
      48     2.394        0.947        0.804        0.579
      40     2.204        0.769        0.572        0.375
      37     2.055        1.561        2.258        3.006
      40     1.929       77.016       27.179     2924.492
      50     1.827       22.328        0.299      226.394
      44     1.739        6.751        6.684       19.153
      45     1.664       62.663      180.534     1744.635
      73     1.535       33.253       21.249      857.745
      55     1.431       34.439       22.142      483.441
      18     1.314        3.950       10.099       24.839
     794   overall       13.653       16.013      360.500

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.250  0.00  0.00  0.00  0.00  0.00  0.00  0.33  0.32  0.30  0.29
       3     4.621  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.494  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.987  0.00  0.00  0.00  0.00  0.00  0.00  0.42  0.40  0.38  0.37
       3     2.646  1.34  0.96  0.80  0.70  0.64  0.59  0.56  1.06  1.01  0.98
       5     2.394  0.00  0.00  0.00  0.42  0.77  1.07  1.00  0.95  0.91  0.88
       3     2.204  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.055  2.01  1.81  1.50  1.32  1.20  1.11  1.05  0.99  1.14  1.10
       7     1.929  1.73  1.24  1.03  0.91  0.82  0.76  0.72  0.68  0.65  0.63
       9     1.739  2.24  1.61  1.33  1.17  1.07  0.99  0.93  0.88  0.85  0.81
       3     1.535  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       6     1.480  3.36  2.41  2.00  1.76  1.60  1.48  1.40  1.33  1.27  1.22
       3     1.431  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      89   overall  1.00  0.75  0.62  0.57  0.54  0.52  0.53  0.52  0.51  0.51


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.250  2.63  2.76  1.93  1.52  1.91  1.85  1.66  1.51  1.40  1.32
      18     4.621  1.17  1.84  1.29  1.01  1.13  1.23  1.10  1.01  1.22  1.14
      27     3.494  0.78  0.82  0.57  0.67  0.66  0.74  0.88  0.94  0.94  0.88
      28     2.987  1.50  1.77  1.24  1.19  1.09  1.19  1.21  1.17  1.08  1.07
      38     2.646  0.83  1.16  1.32  1.12  0.94  0.93  1.05  1.10  1.06  1.04
      48     2.394  1.31  1.15  1.04  1.01  1.11  1.15  1.08  1.06  1.09  1.05
      40     2.204  0.79  0.97  0.87  0.91  0.95  0.94  0.89  1.00  1.01  1.07
      37     2.055  1.70  1.19  1.15  1.15  1.10  0.96  0.91  1.03  1.00  1.03
      40     1.929  5.52  2.90  2.03  1.59  1.33  1.16  1.04  1.00  0.93  0.87
      50     1.827  5.67  2.98  2.08  1.64  1.37  1.20  1.07  0.98  0.91  0.85
      44     1.739  4.78  2.51  1.75  1.45  1.27  1.16  1.04  0.95  0.88  0.83
      45     1.664  5.14  2.82  1.97  1.62  1.41  1.23  1.10  1.01  0.94  0.88
      73     1.535  4.75  2.57  1.80  1.41  1.18  1.06  0.98  0.90  0.83  0.80
      55     1.431  5.54  2.91  2.03  1.60  1.34  1.17  1.05  0.96  0.89  0.83
      18     1.314  4.09  2.45  1.71  1.35  1.13  0.99  0.88  0.81  0.75  0.79
     794   overall  2.49  1.53  1.13  0.95  0.86  0.80  0.75  0.72  0.69  0.67


 List of    126 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   49    3    8    1.50 1593.03 -0.2341E+02  0.1765E+02 "alien"
   37   30   11    1.47 1353.80 -0.1990E+02  0.1413E+02 "alien"
   44    5   14    1.48 1027.49 -0.1510E+02  0.1331E+02 "alien"
   43    5   15    1.47  916.92 -0.1348E+02  0.9806E+01 "alien"
   41   32    4    1.47  867.72 -0.1275E+02  0.9988E+01 "alien"
   40   31    6    1.49  841.08 -0.1236E+02  0.9450E+01 "alien"
   39   31    8    1.48  836.87 -0.1230E+02  0.9217E+01 "alien"
   52    3    5    1.46  826.46 -0.1215E+02  0.1032E+02 "alien"
   28   27   17    1.50  807.16 -0.1186E+02  0.8087E+01 "alien"
   31   22   10    1.80  778.13 -0.1613E+02  0.1267E+02 "alien"
   45    4   13    1.48  766.00 -0.1126E+02  0.9068E+01 "alien"
   42   32    2    1.46  631.39 -0.9280E+01  0.9779E+01 "alien"
   35   24    4    1.79  611.03 -0.1267E+02  0.1502E+02 "alien"
   41   32    5    1.46  590.46 -0.8678E+01  0.9758E+01 "alien"
   51    3    7    1.46  584.31 -0.8588E+01  0.9296E+01 "alien"
   27   10   22    1.47  562.65 -0.8270E+01  0.1013E+02 "alien"
   31    9   21    1.46  561.76 -0.8257E+01  0.1079E+02 "alien"
   34   23    5    1.83  536.57 -0.1113E+02  0.8880E+01 "alien"
   51    2    4    1.50  521.42 -0.7664E+01  0.9220E+01 "alien"
   37   31   11    1.46  492.78 -0.7243E+01  0.9787E+01 "alien"
   35   23    3    1.83  483.07 -0.1002E+02  0.8851E+01 "alien"
   30   21   10    1.85  477.68 -0.9904E+01  0.1201E+02 "alien"
   23   13   23    1.46  465.24 -0.6838E+01  0.8462E+01 "alien"
   35   23    1    1.84  456.47 -0.9465E+01  0.8823E+01 "alien"
   36   30   12    1.47  429.58 -0.6314E+01  0.9400E+01 "alien"
   27   25   19    1.46  418.28 -0.6148E+01  0.9399E+01 "alien"
   15   10   19    1.84  400.34 -0.8301E+01  0.6313E+01 "alien"
   42   32    1    1.46  387.18 -0.5691E+01  0.9724E+01 "alien"
   49    3    9    1.48  386.40 -0.5679E+01  0.8985E+01 "alien"
   52    2    2    1.48  353.97 -0.5203E+01  0.9247E+01 "alien"
   40   31    5    1.50  334.53 -0.4917E+01  0.1132E+02 "alien"
   32    8   20    1.49  317.84 -0.4672E+01  0.7694E+01 "alien"
   51    2    3    1.50  309.85 -0.4554E+01  0.9373E+01 "alien"
   41   31    4    1.49  307.87 -0.4525E+01  0.9789E+01 "alien"
   50    2    6    1.50  291.12 -0.4279E+01  0.9167E+01 "alien"
   34   22    4    1.87  281.43 -0.5835E+01  0.9872E+01 "alien"
   42    5   15    1.49  265.75 -0.3906E+01  0.8462E+01 "alien"
   14    8   19    1.87  203.04 -0.4210E+01  0.5340E+01 "alien"
   42    1    3    1.82  171.85 -0.3563E+01  0.8620E+01 "alien"
   39   31    9    1.46  162.89 -0.2394E+01  0.1088E+02 "alien"
   26   10   22    1.49  149.78 -0.2201E+01  0.8739E+01 "alien"
   46    4   12    1.49  143.78 -0.2113E+01  0.8915E+01 "alien"
   20   18   16    1.85  140.44 -0.2912E+01  0.8785E+01 "alien"
   43    1    2    1.79  139.77 -0.2898E+01  0.8378E+01 "alien"
   41    1    4    1.85  138.26 -0.2867E+01  0.8943E+01 "alien"
   17   15   18    1.82  136.14 -0.2823E+01  0.5358E+01 "alien"
   29   21   11    1.84  130.32 -0.2702E+01  0.8213E+01 "alien"
   22   19   15    1.85  125.50 -0.2602E+01  0.1289E+02 "alien"
   41    1    5    1.83  123.11 -0.2553E+01  0.8573E+01 "alien"
   31    2   14    1.85  119.83 -0.2484E+01  0.6914E+01 "alien"
   36    1   11    1.83  115.77 -0.2400E+01  0.1390E+02 "alien"
   37    1   11    1.80  113.43 -0.2352E+01  0.7991E+01 "alien"
   40   31    7    1.47  113.05 -0.1662E+01  0.1032E+02 "alien"
   36   28    9    1.58   94.08 -0.3241E+02  0.1628E+02 "alien"
   31    3   15    1.79   89.72 -0.1860E+01  0.7329E+01 "alien"
   34   29   13    1.50   86.92 -0.1278E+01  0.8938E+01 "alien"
   39    6   17    1.48   86.41 -0.1270E+01  0.8112E+01 "alien"
   30   22   11    1.79   63.89 -0.1325E+01  0.8024E+01 "alien"
   39   28    4    1.59   50.09 -0.1726E+02  0.1555E+02 "alien"
   37   27    7    1.61   41.43 -0.1427E+02  0.9040E+01 "alien"
   38   28    6    1.59   39.79 -0.1371E+02  0.1279E+02 "alien"
   43    1    0    1.80   35.10 -0.1186E+02  0.1198E+02 "alien"
   36    0   10    1.88   35.06 -0.1180E+02  0.8008E+01 "alien"
   17   17   17    1.86   33.57 -0.9785E+01  0.7972E+01 "alien"
   41    3   13    1.59   32.23 -0.1111E+02  0.8189E+01 "alien"
   27   27   18    1.47   30.26 -0.4567E+01  0.9109E+01 "alien"
   45    2    9    1.59   30.24 -0.1042E+02  0.9082E+01 "alien"
   18   18   17    1.82   30.17 -0.6374E+01  0.8563E+01 "alien"
   40    0    8    1.79   28.68 -0.5274E+01  0.8357E+01 "alien"
   30   21    9    1.89   27.90  0.1509E+02  0.1337E+02 "alien"
   42    3   13    1.56   27.73  0.3017E+02  0.1654E+02 "alien"
   32   24   11    1.69   27.00  0.1979E+02  0.9068E+01 "alien"
   39    0    9    1.80   26.70 -0.4105E+01  0.8476E+01 "alien"
   45    2   10    1.57   23.88 -0.8229E+01  0.9591E+01 "alien"
   39    2    1    1.98   23.51  0.1272E+02  0.9032E+01 "alien"
   28   20   11    1.89   22.26  0.1204E+02  0.8403E+01 "alien"
   26   23   16    1.64   21.79  0.1597E+02  0.1432E+02 "alien"
   35   23    0    1.85   20.79 -0.1918E+01  0.1253E+02 "alien"
   34    5   17    1.60   20.61 -0.7102E+01  0.7409E+01 "alien"
   27   24   16    1.60   20.42 -0.7036E+01  0.1293E+02 "alien"
   48    1    4    1.59   20.35 -0.7010E+01  0.8680E+01 "alien"
   21   17   20    1.60   20.14 -0.6940E+01  0.8422E+01 "alien"
   35   26    8    1.66   20.09  0.1472E+02  0.1618E+02 "alien"
   21   20   18    1.67   19.98  0.1464E+02  0.1124E+02 "alien"
   46    5   14    1.43   19.73  0.1678E+02  0.1405E+02 "alien"
   42    6   16    1.46   19.38  0.1648E+02  0.9144E+01 "alien"
   24   23   17    1.63   19.07 -0.6571E+01  0.7668E+01 "alien"
   24   13   24    1.40   17.97 -0.1938E+02  0.7995E+01 "alien"
   37   27    6    1.63   17.54  0.1285E+02  0.1268E+02 "alien"
   36    1    8    1.96   16.89  0.9136E+01  0.8723E+01 "alien"
   21   17   14    1.99   16.81  0.9092E+01  0.7880E+01 "alien"
   29   20    9    1.95   16.58  0.8966E+01  0.1333E+02 "alien"
   22   21   17    1.70   16.13  0.1182E+02  0.1016E+02 "alien"
   21   15   21    1.57   15.91 -0.5483E+01  0.9032E+01 "alien"
   39    2   12    1.69   15.89  0.1164E+02  0.8124E+01 "alien"
   17    9   21    1.67   15.85  0.1161E+02  0.6024E+01 "alien"
   50    2    4    1.53   15.52  0.1688E+02  0.1590E+02 "alien"
   43    4   13    1.53   15.48  0.1684E+02  0.1606E+02 "alien"
   40   29    0    1.56   15.22 -0.1510E+02  0.1333E+02 "alien"
   48    1    3    1.60   15.08 -0.5195E+01  0.9117E+01 "alien"
   52    5   11    1.36   14.99 -0.2802E+02  0.1543E+02 "alien"
   48    6   14    1.38   14.49 -0.1562E+02  0.1298E+02 "alien"
   36    3   14    1.69   14.33  0.1050E+02  0.8800E+01 "alien"
   40    2   12    1.66   14.05  0.1030E+02  0.8399E+01 "alien"
   33    7   19    1.52   13.75  0.1496E+02  0.7512E+01 "alien"
   40    7   18    1.43   13.55  0.1152E+02  0.1588E+02 "alien"
   45    7   16    1.39   13.39 -0.1444E+02  0.9467E+01 "alien"
   35   27   10    1.59   13.32 -0.4591E+01  0.8704E+01 "alien"
   36   32   14    1.39   13.10 -0.1413E+02  0.1725E+02 "alien"
   15   12   18    1.90   12.67  0.6856E+01  0.6519E+01 "alien"
   51    2    2    1.51   12.28  0.1336E+02  0.1026E+02 "alien"
   16   15   17    1.91   12.16  0.6579E+01  0.1018E+02 "alien"
   36   26    7    1.66   11.93  0.8742E+01  0.9056E+01 "alien"
   37   31   12    1.43   11.67  0.9927E+01  0.1547E+02 "alien"
   40    2    0    1.93   11.55  0.1392E+02  0.1358E+02 "alien"
   28   28   18    1.45   11.18  0.2069E+02  0.1233E+02 "alien"
   27    9   21    1.52   10.83  0.1178E+02  0.8537E+01 "alien"
   19    6   23    1.54   10.75  0.1170E+02  0.6310E+01 "alien"
   40   33    9    1.41   10.44 -0.1126E+02  0.9857E+01 "alien"
   15    7   20    1.79   10.34 -0.2144E+00  0.5508E+01 "alien"
   35   33   15    1.36   10.34 -0.1932E+02  0.1370E+02 "alien"
   28   12   23    1.40   10.27 -0.1108E+02  0.1173E+02 "alien"
   33   28   13    1.53   10.21  0.1111E+02  0.1459E+02 "alien"
   35    6   18    1.53   10.21  0.1111E+02  0.1329E+02 "alien"
   34    4   16    1.65   10.15  0.7440E+01  0.8078E+01 "alien"
   53    7   13    1.30   10.14  0.6612E+02  0.1505E+02 "alien"


 List of    126 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   10.14  0.6612E+02  0.1505E+02
   35   33   15    1.36   10.34 -0.1932E+02  0.1370E+02
   52    5   11    1.36   14.99 -0.2802E+02  0.1543E+02
   48    6   14    1.38   14.49 -0.1562E+02  0.1298E+02
   36   32   14    1.39   13.10 -0.1413E+02  0.1725E+02
   45    7   16    1.39   13.39 -0.1444E+02  0.9467E+01
   24   13   24    1.40   17.97 -0.1938E+02  0.7995E+01
   28   12   23    1.40   10.27 -0.1108E+02  0.1173E+02
   40   33    9    1.41   10.44 -0.1126E+02  0.9857E+01
   40    7   18    1.43   13.55  0.1152E+02  0.1588E+02
   46    5   14    1.43   19.73  0.1678E+02  0.1405E+02
   37   31   12    1.43   11.67  0.9927E+01  0.1547E+02
   28   28   18    1.45   11.18  0.2069E+02  0.1233E+02
   42    6   16    1.46   19.38  0.1648E+02  0.9144E+01
   37   31   11    1.46  492.78 -0.7243E+01  0.9787E+01
   52    3    5    1.46  826.46 -0.1215E+02  0.1032E+02
   23   13   23    1.46  465.24 -0.6838E+01  0.8462E+01
   51    3    7    1.46  584.31 -0.8588E+01  0.9296E+01
   39   31    9    1.46  162.89 -0.2394E+01  0.1088E+02
   41   32    5    1.46  590.46 -0.8678E+01  0.9758E+01
   42   32    2    1.46  631.39 -0.9280E+01  0.9779E+01
   31    9   21    1.46  561.76 -0.8257E+01  0.1079E+02
   27   25   19    1.46  418.28 -0.6148E+01  0.9399E+01
   42   32    1    1.46  387.18 -0.5691E+01  0.9724E+01
   43    5   15    1.47  916.92 -0.1348E+02  0.9806E+01
   36   30   12    1.47  429.58 -0.6314E+01  0.9400E+01
   41   32    4    1.47  867.72 -0.1275E+02  0.9988E+01
   27   10   22    1.47  562.65 -0.8270E+01  0.1013E+02
   40   31    7    1.47  113.05 -0.1662E+01  0.1032E+02
   37   30   11    1.47 1353.80 -0.1990E+02  0.1413E+02
   27   27   18    1.47   30.26 -0.4567E+01  0.9109E+01
   44    5   14    1.48 1027.49 -0.1510E+02  0.1331E+02
   39   31    8    1.48  836.87 -0.1230E+02  0.9217E+01
   49    3    9    1.48  386.40 -0.5679E+01  0.8985E+01
   52    2    2    1.48  353.97 -0.5203E+01  0.9247E+01
   45    4   13    1.48  766.00 -0.1126E+02  0.9068E+01
   39    6   17    1.48   86.41 -0.1270E+01  0.8112E+01
   26   10   22    1.49  149.78 -0.2201E+01  0.8739E+01
   46    4   12    1.49  143.78 -0.2113E+01  0.8915E+01
   41   31    4    1.49  307.87 -0.4525E+01  0.9789E+01
   40   31    6    1.49  841.08 -0.1236E+02  0.9450E+01
   42    5   15    1.49  265.75 -0.3906E+01  0.8462E+01
   32    8   20    1.49  317.84 -0.4672E+01  0.7694E+01
   34   29   13    1.50   86.92 -0.1278E+01  0.8938E+01
   51    2    4    1.50  521.42 -0.7664E+01  0.9220E+01
   28   27   17    1.50  807.16 -0.1186E+02  0.8087E+01
   49    3    8    1.50 1593.03 -0.2341E+02  0.1765E+02
   40   31    5    1.50  334.53 -0.4917E+01  0.1132E+02
   50    2    6    1.50  291.12 -0.4279E+01  0.9167E+01
   51    2    3    1.50  309.85 -0.4554E+01  0.9373E+01
   51    2    2    1.51   12.28  0.1336E+02  0.1026E+02
   33    7   19    1.52   13.75  0.1496E+02  0.7512E+01
   27    9   21    1.52   10.83  0.1178E+02  0.8537E+01
   50    2    4    1.53   15.52  0.1688E+02  0.1590E+02
   35    6   18    1.53   10.21  0.1111E+02  0.1329E+02
   33   28   13    1.53   10.21  0.1111E+02  0.1459E+02
   43    4   13    1.53   15.48  0.1684E+02  0.1606E+02
   19    6   23    1.54   10.75  0.1170E+02  0.6310E+01
   42    3   13    1.56   27.73  0.3017E+02  0.1654E+02
   40   29    0    1.56   15.22 -0.1510E+02  0.1333E+02
   45    2   10    1.57   23.88 -0.8229E+01  0.9591E+01
   21   15   21    1.57   15.91 -0.5483E+01  0.9032E+01
   36   28    9    1.58   94.08 -0.3241E+02  0.1628E+02
   38   28    6    1.59   39.79 -0.1371E+02  0.1279E+02
   48    1    4    1.59   20.35 -0.7010E+01  0.8680E+01
   39   28    4    1.59   50.09 -0.1726E+02  0.1555E+02
   41    3   13    1.59   32.23 -0.1111E+02  0.8189E+01
   35   27   10    1.59   13.32 -0.4591E+01  0.8704E+01
   45    2    9    1.59   30.24 -0.1042E+02  0.9082E+01
   48    1    3    1.60   15.08 -0.5195E+01  0.9117E+01
   34    5   17    1.60   20.61 -0.7102E+01  0.7409E+01
   27   24   16    1.60   20.42 -0.7036E+01  0.1293E+02
   21   17   20    1.60   20.14 -0.6940E+01  0.8422E+01
   37   27    7    1.61   41.43 -0.1427E+02  0.9040E+01
   24   23   17    1.63   19.07 -0.6571E+01  0.7668E+01
   37   27    6    1.63   17.54  0.1285E+02  0.1268E+02
   26   23   16    1.64   21.79  0.1597E+02  0.1432E+02
   34    4   16    1.65   10.15  0.7440E+01  0.8078E+01
   40    2   12    1.66   14.05  0.1030E+02  0.8399E+01
   36   26    7    1.66   11.93  0.8742E+01  0.9056E+01
   35   26    8    1.66   20.09  0.1472E+02  0.1618E+02
   21   20   18    1.67   19.98  0.1464E+02  0.1124E+02
   17    9   21    1.67   15.85  0.1161E+02  0.6024E+01
   39    2   12    1.69   15.89  0.1164E+02  0.8124E+01
   36    3   14    1.69   14.33  0.1050E+02  0.8800E+01
   32   24   11    1.69   27.00  0.1979E+02  0.9068E+01
   22   21   17    1.70   16.13  0.1182E+02  0.1016E+02
   15    7   20    1.79   10.34 -0.2144E+00  0.5508E+01
   31    3   15    1.79   89.72 -0.1860E+01  0.7329E+01
   30   22   11    1.79   63.89 -0.1325E+01  0.8024E+01
   43    1    2    1.79  139.77 -0.2898E+01  0.8378E+01
   35   24    4    1.79  611.03 -0.1267E+02  0.1502E+02
   40    0    8    1.79   28.68 -0.5274E+01  0.8357E+01
   37    1   11    1.80  113.43 -0.2352E+01  0.7991E+01
   43    1    0    1.80   35.10 -0.1186E+02  0.1198E+02
   39    0    9    1.80   26.70 -0.4105E+01  0.8476E+01
   31   22   10    1.80  778.13 -0.1613E+02  0.1267E+02
   17   15   18    1.82  136.14 -0.2823E+01  0.5358E+01
   18   18   17    1.82   30.17 -0.6374E+01  0.8563E+01
   42    1    3    1.82  171.85 -0.3563E+01  0.8620E+01
   35   23    3    1.83  483.07 -0.1002E+02  0.8851E+01
   34   23    5    1.83  536.57 -0.1113E+02  0.8880E+01
   41    1    5    1.83  123.11 -0.2553E+01  0.8573E+01
   36    1   11    1.83  115.77 -0.2400E+01  0.1390E+02
   15   10   19    1.84  400.34 -0.8301E+01  0.6313E+01
   29   21   11    1.84  130.32 -0.2702E+01  0.8213E+01
   35   23    1    1.84  456.47 -0.9465E+01  0.8823E+01
   35   23    0    1.85   20.79 -0.1918E+01  0.1253E+02
   31    2   14    1.85  119.83 -0.2484E+01  0.6914E+01
   20   18   16    1.85  140.44 -0.2912E+01  0.8785E+01
   41    1    4    1.85  138.26 -0.2867E+01  0.8943E+01
   22   19   15    1.85  125.50 -0.2602E+01  0.1289E+02
   30   21   10    1.85  477.68 -0.9904E+01  0.1201E+02
   17   17   17    1.86   33.57 -0.9785E+01  0.7972E+01
   34   22    4    1.87  281.43 -0.5835E+01  0.9872E+01
   14    8   19    1.87  203.04 -0.4210E+01  0.5340E+01
   36    0   10    1.88   35.06 -0.1180E+02  0.8008E+01
   28   20   11    1.89   22.26  0.1204E+02  0.8403E+01
   30   21    9    1.89   27.90  0.1509E+02  0.1337E+02
   15   12   18    1.90   12.67  0.6856E+01  0.6519E+01
   16   15   17    1.91   12.16  0.6579E+01  0.1018E+02
   40    2    0    1.93   11.55  0.1392E+02  0.1358E+02
   29   20    9    1.95   16.58  0.8966E+01  0.1333E+02
   36    1    8    1.96   16.89  0.9136E+01  0.8723E+01
   39    2    1    1.98   23.51  0.1272E+02  0.9032E+01
   21   17   14    1.99   16.81  0.9092E+01  0.7880E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        83
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1484
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1314
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      170


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
