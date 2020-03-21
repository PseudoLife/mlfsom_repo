
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.402    77.402    38.748  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      16      16
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      16
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000819  0.000682
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
 REFLECTING_RANGE_E.S.D.=   0.098
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

    1452 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      16 ...      16
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1452 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               44
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
 *  25        mC          0.0     109.5  109.5   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  31        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  23        oC          0.0     109.5  109.5   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  35        mP          0.0      77.4   38.7   77.4  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  20        mC          0.0     109.5  109.5   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  33        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  21        tP          0.0      77.4   77.4   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
 *  32        oP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    37        mC        250.0     159.6   38.7   77.4  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.6   38.7   77.4  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.7  159.6   77.4  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.7  159.6   77.4  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.7  159.6   77.4  90.0  90.0  75.9    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.7  159.6   77.4  90.0  90.0  75.9    0  0  1  0  0  2  1  0 -1  0  0  0
    19        oI        500.0      38.7  109.5  116.1  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    27        mC        500.0     159.6   38.7  109.5  90.0 133.3  75.9    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    18        tI        625.0     109.5  116.1   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    26        oF        625.0      38.7  159.6  159.6  86.6 104.1 104.1    0  0  1  0  0  2 -1  0 -2  0 -1  0
    41        mC        997.5     173.1   77.4   38.7  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    22        hP        997.5      77.4   77.4   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    30        mC        997.5      77.4  173.1   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        997.5      77.4  173.1   38.7  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
     1        cF        999.0     116.1  116.1  116.1  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.6   86.6  116.1 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.6   86.6  116.1 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.5   86.6   86.6  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.6   86.6  109.5  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.6  247.8 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.6   86.6   77.4  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.6   86.6   77.4  90.0  90.0  53.2   -1  0  1  0 -1  0 -1  0  0 -1  0  0
    15        tI        999.0      38.7   77.4  177.4  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.6   86.6  177.4 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.6   86.6   77.4  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.6   86.6   86.6  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.4  109.5   38.7  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.7   77.4  177.4 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.6  177.4   77.4 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    5.5      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.766   -23.442     0.023
 UNIT_CELL_B-AXIS=    23.442    73.766     0.014
 UNIT_CELL_C-AXIS=    -0.013    -0.003    38.749
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       289 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.26
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.389    77.389    38.728  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  4.4E-02 4.4E-02 1.9E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012922  0.012922  0.025821  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.748   -23.458     0.024
 COORDINATES OF UNIT CELL B-AXIS    23.458    73.748     0.007
 COORDINATES OF UNIT CELL C-AXIS    -0.013     0.001    38.728
 CRYSTAL MOSAICITY (DEGREES)     0.098
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006351  0.001819
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.92   1232.90
 DETECTOR ORIGIN (PIXELS) AT                     1261.27   1242.14
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.27
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.39    77.39    38.73  90.000  90.000  90.000


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
  0.05  0.15       106        158.435       1.211
  0.15  0.25       168        212.551       0.903
  0.25  0.35       201        162.508       1.181
  0.35  0.45       141        117.737       1.630
  0.45  0.55       170        269.118       0.713
  0.55  0.65       202        173.151       1.108
  0.65  0.75       183        178.507       1.075
  0.75  0.85       173        174.501       1.100
  0.85  0.95       108        308.537       0.622
 
 -0.05  0.95      1452        191.923       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               237.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.009
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    15.1 XMAX=    15.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59550 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        514


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               237.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.006
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     8.6 XMAX=  2456.4 NXBIN=    3
 YMIN=     6.4 YMAX=  2517.5 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        514


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               237.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.001
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    15.1 XMAX=    15.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1262
 DETECTOR_SURFACE_POSITION=    1648    1688
 DETECTOR_SURFACE_POSITION=     817    1688
 DETECTOR_SURFACE_POSITION=     817     835
 DETECTOR_SURFACE_POSITION=    1648     835
 DETECTOR_SURFACE_POSITION=    2173    1662
 DETECTOR_SURFACE_POSITION=    1622    2227
 DETECTOR_SURFACE_POSITION=     843    2227
 DETECTOR_SURFACE_POSITION=     292    1662
 DETECTOR_SURFACE_POSITION=     292     862
 DETECTOR_SURFACE_POSITION=     843     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2173     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        514


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
 9.896E-01  1.912E-03   22.99


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

   24.473   3.400    19.96   1.34      2.43      2.56      62      83      26
    3.400   2.416    13.31   1.60      4.37      3.48     136     174      10
    2.416   1.976     6.11   0.98      7.81      7.70     182     209       6
    1.976   1.712    -0.00   0.74   1826.13   2090.46     190     231       2
    1.712   1.532    -0.04   0.66    -99.90    -99.90     220     268       0
    1.532   1.399     0.09   0.78    822.69    931.22     242     296       2
    1.399   1.296     0.08   0.43    373.75    568.80      96     144       0
   --------------------------------------------------------------------------
   24.473   1.296     3.76   0.88      5.92      6.17    1128    1405      46
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      16      16)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.682   0.1382E+02  0.1217E+02         1.14           1 
   47    0    0     1.647   0.8751E+01  0.1203E+02         0.73           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     4.3% OF MEAN INTENSITY


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

     3.84          60      39      1285        3.0%       2.2%      2.5%       42   25.17      3.1%    99.9*     0    0.000       0
     2.74         108      66      2098        3.1%       3.6%      3.0%       84   20.24      5.1%    99.5*     0    0.000       0
     2.24         163      96      2665        3.6%       5.3%      4.4%      134   13.34      7.4%    99.0*     0    0.000       0
     1.94         160      90      3101        2.9%      11.0%     11.2%      140    5.11     15.6%    97.7*     0    0.000       0
     1.74         176     103      3479        3.0%    1635.6%   2045.5%      146   -0.02   2313.0%     9.4      0    0.000       0
     1.59         209     122      3851        3.2%     -99.9%    -99.9%      174    0.00    -99.9%    -3.9      0    0.000       0
     1.47         229     137      4143        3.3%     -99.9%    -99.9%      184    0.00    -99.9%    18.6      0    0.000       0
     1.37         215     126      4427        2.8%     573.5%    675.9%      178    0.15    811.0%   -15.6      0    0.000       0
     1.30          85      62      4743        1.3%     552.2%    731.3%       46   -0.01    781.0%    45.5      0    0.000       0
    total        1405     841     29792        2.8%       5.9%      6.2%     1128    4.85      8.4%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1452
 NUMBER OF REJECTED MISFITS                                46
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1405
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    841



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

   WILSON LINE (using all data) : A=   7.423 B=  30.680 CORRELATION=  0.84
      #      RES      SS        <I>       log(<I>)       BO
      16     7.404   0.005  2.1759E+03       7.685     -28.7
      19     4.637   0.012  2.7586E+03       7.922     -21.5
      29     3.497   0.020  1.0655E+03       6.971      11.1
      32     2.991   0.028  4.7015E+02       6.153      22.7
      41     2.643   0.036  4.3996E+02       6.087      18.7
      50     2.401   0.043  2.1295E+02       5.361      23.8
      39     2.207   0.051  1.3746E+02       4.923      24.3
      45     2.057   0.059  5.3353E+01       3.977      29.2
      44     1.930   0.067  3.8583E-01      -0.952      62.4
      55     1.830   0.075  3.5708E-01      -1.030      56.6
      51     1.741   0.083 -1.0176E+00    -999.999   -1000.0
      53     1.666   0.090 -3.9172E-01    -999.999   -1000.0
      54     1.598   0.098 -4.9712E-01    -999.999   -1000.0
      72     1.537   0.106  7.6808E-01      -0.264      36.3
      75     1.481   0.114  1.5283E+00       0.424      30.7
      57     1.433   0.122  1.1267E+00       0.119      30.0
      48     1.390   0.129  2.9824E+00       1.093      24.5
      43     1.352   0.137 -1.6594E+00    -999.999   -1000.0
      18     1.315   0.145  2.2769E+00       0.823      22.8


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.404        1.647        1.038        0.520
       3     4.637        1.056        0.590        0.245
       3     3.497        0.325        0.079        0.014
       4     2.991        1.021        0.433        0.139
       3     2.643        0.246        0.053        0.008
       5     2.401        1.007        0.563        0.244
       2     2.207        1.209        0.466        0.129
       8     2.057        0.370        0.147        0.045
      10     1.830      179.462      481.701     7811.159
       9     1.741       32.829       38.095      293.354
      11     1.666       12.937        2.279       17.717
       3     1.598       20.616       29.830       54.992
       3     1.537       24.680       12.031       80.173
       6     1.481       13.744       19.045       21.382
       3     1.433       56.615      149.609      531.977
      87   overall       38.229        4.896     1384.558


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.404        0.431        0.280        0.146
      16     4.637        1.017        1.046        0.908
      26     3.497        0.865        0.695        0.530
      28     2.991        0.942        0.917        0.801
      38     2.643        1.138        1.211        1.157
      45     2.401        0.922        0.785        0.568
      37     2.207        0.816        0.634        0.433
      37     2.057        1.520        2.164        2.797
      37     1.930      139.907       73.654     9585.925
      45     1.830      239.068      -68.471    23520.184
      69     1.537       56.756       95.632     1648.029
      69     1.481       16.035        5.723      136.317
      54     1.433       26.939       64.051      614.015
      48     1.390        3.225        2.808        6.091
      18     1.315       15.194       39.153      231.268
     754   overall       30.612       14.833     2087.614

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.404  0.00  0.00  0.00  0.00  0.00  0.00  0.33  0.32  0.30  0.59
       3     4.637  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       3     3.497  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.991  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.643  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       5     2.401  0.00  0.00  0.00  0.00  0.38  1.07  1.00  0.95  0.91  0.88
       2     2.207  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.057  1.01  1.45  1.20  1.32  1.20  1.11  1.05  0.99  0.95  1.10
      10     1.830  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       9     1.741  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
      11     1.666  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.72  0.69  0.67
       3     1.598  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.537  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       6     1.481  0.67  0.48  0.40  0.35  0.32  0.30  0.28  0.27  0.25  0.24
       3     1.433  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      87   overall  1.25  0.96  0.80  0.73  0.68  0.70  0.67  0.64  0.63  0.66


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.404  1.91  2.51  1.75  1.38  1.85  1.81  1.63  1.49  1.38  1.29
      16     4.637  0.66  1.38  0.96  0.76  1.27  1.11  0.99  1.25  1.16  1.19
      26     3.497  0.81  0.64  0.59  0.58  0.59  0.68  0.84  0.98  0.91  0.85
      28     2.991  1.88  1.58  1.38  1.19  1.18  1.19  1.14  1.17  1.08  1.07
      38     2.643  0.55  1.16  1.42  1.12  0.94  0.93  0.99  1.05  1.11  1.08
      45     2.401  1.17  1.23  1.03  1.08  1.13  1.13  1.10  1.09  1.09  1.09
      37     2.207  0.57  1.04  1.04  0.98  0.89  0.90  0.86  0.98  1.00  1.07
      37     2.057  2.27  1.49  1.15  0.98  1.10  1.08  0.97  0.98  1.00  1.03
      37     1.930  3.98  2.09  1.46  1.15  0.96  0.90  0.81  0.74  0.68  0.64
      45     1.830  4.67  2.45  1.71  1.35  1.13  0.99  0.88  0.81  0.75  0.70
      69     1.537  5.48  2.96  2.07  1.63  1.36  1.19  1.07  0.97  0.93  0.87
      69     1.481  5.03  2.64  1.90  1.49  1.25  1.09  1.01  0.92  0.85  0.80
      54     1.433  4.67  2.45  1.79  1.40  1.22  1.07  0.96  0.87  0.81  0.76
      48     1.390  3.28  1.84  1.37  1.07  0.95  0.83  0.74  0.72  0.70  0.76
      18     1.315  4.67  2.45  1.71  1.35  1.13  1.11  0.99  0.91  0.84  0.79
     754   overall  2.47  1.51  1.15  0.94  0.86  0.80  0.74  0.73  0.70  0.68


 List of     90 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   33   26   12    1.60   46.79 -0.2326E+02  0.9554E+01 "alien"
   42    1    9    1.69   45.48 -0.1781E+02  0.8140E+01 "alien"
   36   24    1    1.79   42.23  0.1508E+02  0.8600E+01 "alien"
   15   10   19    1.84   41.60  0.1485E+02  0.6907E+01 "alien"
   41    1    3    1.87   37.48  0.1338E+02  0.1046E+02 "alien"
   35   25    8    1.69   35.70 -0.1398E+02  0.1304E+02 "alien"
   25   18   12    1.98   35.11  0.1355E+02  0.1228E+02 "alien"
   40    1    5    1.88   34.86  0.1245E+02  0.1163E+02 "alien"
   37   26    5    1.67   34.70 -0.1359E+02  0.8561E+01 "alien"
   29   20   10    1.91   32.63  0.1259E+02  0.1380E+02 "alien"
   31   26   14    1.57   31.64 -0.1573E+02  0.8242E+01 "alien"
   43    1    1    1.80   31.46  0.1123E+02  0.8368E+01 "alien"
   27   19   11    1.95   29.79  0.1149E+02  0.8550E+01 "alien"
   21    9   21    1.62   29.68 -0.1476E+02  0.8219E+01 "alien"
   35   26    9    1.64   28.54 -0.1118E+02  0.8302E+01 "alien"
   42   32    3    1.46   27.83  0.3136E+02  0.1721E+02 "alien"
   19    9   21    1.65   27.70 -0.1085E+02  0.6136E+01 "alien"
   26   19   12    1.93   26.69  0.1030E+02  0.7872E+01 "alien"
   45    1    8    1.62   26.27 -0.1306E+02  0.9118E+01 "alien"
   35    5   17    1.58   26.20 -0.1302E+02  0.1140E+02 "alien"
   47    3   10    1.51   26.16  0.2009E+02  0.9108E+01 "alien"
   33    4   16    1.68   25.73 -0.1008E+02  0.6693E+01 "alien"
   29   27   16    1.52   25.71  0.1975E+02  0.1434E+02 "alien"
   48    1    5    1.58   25.36 -0.1261E+02  0.8543E+01 "alien"
   27   20   12    1.87   25.14  0.8979E+01  0.7504E+01 "alien"
   29   21   11    1.84   24.77  0.8843E+01  0.7816E+01 "alien"
   47    1    6    1.60   22.27 -0.1107E+02  0.1199E+02 "alien"
   36   27    8    1.62   22.01 -0.1094E+02  0.8460E+01 "alien"
   46    1    6    1.63   21.18 -0.1053E+02  0.1179E+02 "alien"
   27   21   13    1.80   21.02  0.7505E+01  0.7519E+01 "alien"
   38    5   16    1.55   20.37  0.1565E+02  0.8123E+01 "alien"
   51    2    2    1.51   20.10  0.1543E+02  0.1015E+02 "alien"
   30    8   20    1.53   19.51  0.1499E+02  0.7708E+01 "alien"
   37   29   10    1.52   18.96  0.1456E+02  0.1489E+02 "alien"
   39    4   15    1.57   18.05 -0.8975E+01  0.1098E+02 "alien"
   28   23   14    1.69   17.87 -0.6999E+01  0.7412E+01 "alien"
   26   20   13    1.85   17.74  0.6335E+01  0.7453E+01 "alien"
   22   18   20    1.58   17.72 -0.8807E+01  0.6259E+01 "alien"
   31   21    8    1.90   17.53  0.6764E+01  0.8388E+01 "alien"
   53    7   13    1.30   17.20  0.3917E+02  0.1362E+02 "alien"
   31    2   14    1.85   17.01  0.6075E+01  0.7132E+01 "alien"
   38   34   13    1.35   16.95 -0.2812E+02  0.9650E+01 "alien"
   33   22    5    1.89   16.74  0.6460E+01  0.1436E+02 "alien"
   47    2    9    1.54   15.75  0.1210E+02  0.1502E+02 "alien"
   37    1   11    1.80   15.19  0.5425E+01  0.7696E+01 "alien"
   46    2    9    1.57   15.16  0.1164E+02  0.8587E+01 "alien"
   43    1    8    1.69   14.79 -0.5794E+01  0.1373E+02 "alien"
   41    5   15    1.52   14.73  0.1131E+02  0.1397E+02 "alien"
   36   24    0    1.79   14.73  0.1460E+02  0.1245E+02 "alien"
   40    0    8    1.80   14.72  0.1458E+02  0.8131E+01 "alien"
   21   19   16    1.81   14.53  0.5189E+01  0.9404E+01 "alien"
   31    8   20    1.51   14.18  0.1089E+02  0.6904E+01 "alien"
   42    1    2    1.83   14.06  0.5021E+01  0.8891E+01 "alien"
   34    7   19    1.50   14.02  0.2143E+02  0.1120E+02 "alien"
   34   27   11    1.59   13.53 -0.6724E+01  0.8988E+01 "alien"
   31   22    9    1.84   13.51  0.4825E+01  0.9707E+01 "alien"
   40    2    1    1.93   13.49  0.5206E+01  0.1006E+02 "alien"
   39    2    1    1.98   13.49  0.5206E+01  0.8653E+01 "alien"
   30   28   16    1.49   13.49  0.2062E+02  0.1226E+02 "alien"
   39   29    4    1.57   13.12 -0.6523E+01  0.9305E+01 "alien"
   36   24    3    1.77   12.98 -0.1320E+02  0.9444E+01 "alien"
   14   10   18    1.94   12.60  0.4863E+01  0.6070E+01 "alien"
   39    1   10    1.77   12.55 -0.1277E+02  0.1304E+02 "alien"
   40    2   12    1.66   12.52 -0.4905E+01  0.7933E+01 "alien"
   52    4    8    1.42   12.52  0.1411E+02  0.9792E+01 "alien"
   35   23    1    1.85   12.50  0.4464E+01  0.8914E+01 "alien"
   30   25   14    1.61   12.42 -0.6174E+01  0.8147E+01 "alien"
   21   10   22    1.56   12.33  0.9471E+01  0.7800E+01 "alien"
   36   26    7    1.66   11.99 -0.4697E+01  0.8764E+01 "alien"
   32   23   10    1.75   11.97 -0.1218E+02  0.1468E+02 "alien"
   48    1    4    1.59   11.75 -0.5839E+01  0.8602E+01 "alien"
   38   28    6    1.59   11.58 -0.5756E+01  0.1243E+02 "alien"
   38    1    7    1.91   11.55  0.4456E+01  0.8278E+01 "alien"
   45    2   10    1.57   11.54 -0.5736E+01  0.9064E+01 "alien"
   40    2    0    1.93   11.29  0.9540E+01  0.1319E+02 "alien"
   38    1    6    1.94   11.22  0.4330E+01  0.8966E+01 "alien"
   39   28    1    1.61   10.99 -0.5461E+01  0.8843E+01 "alien"
   36   32   13    1.41   10.92  0.1230E+02  0.1597E+02 "alien"
   36    2   13    1.74   10.88 -0.1107E+02  0.7298E+01 "alien"
   41   32    6    1.45   10.84  0.1222E+02  0.9739E+01 "alien"
   32   24   11    1.70   10.82 -0.4240E+01  0.8270E+01 "alien"
   38   29    8    1.54   10.76  0.8263E+01  0.9184E+01 "alien"
   34   23    5    1.83   10.71  0.3824E+01  0.8926E+01 "alien"
   27   22   14    1.73   10.69 -0.1087E+02  0.7191E+01 "alien"
   37   10   21    1.36   10.39 -0.1723E+02  0.1154E+02 "alien"
   39    2    3    1.96   10.35  0.3992E+01  0.1030E+02 "alien"
   39   34   12    1.36   10.25 -0.1701E+02  0.1555E+02 "alien"
   34   22    2    1.90   10.22  0.3943E+01  0.8636E+01 "alien"
   38   26    0    1.68   10.15 -0.8096E+01  0.1596E+02 "alien"
   53    5   10    1.36   10.13 -0.1680E+02  0.9676E+01 "alien"


 List of     90 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   17.20  0.3917E+02  0.1362E+02
   38   34   13    1.35   16.95 -0.2812E+02  0.9650E+01
   39   34   12    1.36   10.25 -0.1701E+02  0.1555E+02
   53    5   10    1.36   10.13 -0.1680E+02  0.9676E+01
   37   10   21    1.36   10.39 -0.1723E+02  0.1154E+02
   36   32   13    1.41   10.92  0.1230E+02  0.1597E+02
   52    4    8    1.42   12.52  0.1411E+02  0.9792E+01
   41   32    6    1.45   10.84  0.1222E+02  0.9739E+01
   42   32    3    1.46   27.83  0.3136E+02  0.1721E+02
   30   28   16    1.49   13.49  0.2062E+02  0.1226E+02
   34    7   19    1.50   14.02  0.2143E+02  0.1120E+02
   31    8   20    1.51   14.18  0.1089E+02  0.6904E+01
   51    2    2    1.51   20.10  0.1543E+02  0.1015E+02
   47    3   10    1.51   26.16  0.2009E+02  0.9108E+01
   37   29   10    1.52   18.96  0.1456E+02  0.1489E+02
   41    5   15    1.52   14.73  0.1131E+02  0.1397E+02
   29   27   16    1.52   25.71  0.1975E+02  0.1434E+02
   30    8   20    1.53   19.51  0.1499E+02  0.7708E+01
   38   29    8    1.54   10.76  0.8263E+01  0.9184E+01
   47    2    9    1.54   15.75  0.1210E+02  0.1502E+02
   38    5   16    1.55   20.37  0.1565E+02  0.8123E+01
   21   10   22    1.56   12.33  0.9471E+01  0.7800E+01
   46    2    9    1.57   15.16  0.1164E+02  0.8587E+01
   39    4   15    1.57   18.05 -0.8975E+01  0.1098E+02
   45    2   10    1.57   11.54 -0.5736E+01  0.9064E+01
   39   29    4    1.57   13.12 -0.6523E+01  0.9305E+01
   31   26   14    1.57   31.64 -0.1573E+02  0.8242E+01
   22   18   20    1.58   17.72 -0.8807E+01  0.6259E+01
   48    1    5    1.58   25.36 -0.1261E+02  0.8543E+01
   35    5   17    1.58   26.20 -0.1302E+02  0.1140E+02
   38   28    6    1.59   11.58 -0.5756E+01  0.1243E+02
   48    1    4    1.59   11.75 -0.5839E+01  0.8602E+01
   34   27   11    1.59   13.53 -0.6724E+01  0.8988E+01
   47    1    6    1.60   22.27 -0.1107E+02  0.1199E+02
   33   26   12    1.60   46.79 -0.2326E+02  0.9554E+01
   39   28    1    1.61   10.99 -0.5461E+01  0.8843E+01
   30   25   14    1.61   12.42 -0.6174E+01  0.8147E+01
   21    9   21    1.62   29.68 -0.1476E+02  0.8219E+01
   45    1    8    1.62   26.27 -0.1306E+02  0.9118E+01
   36   27    8    1.62   22.01 -0.1094E+02  0.8460E+01
   46    1    6    1.63   21.18 -0.1053E+02  0.1179E+02
   35   26    9    1.64   28.54 -0.1118E+02  0.8302E+01
   19    9   21    1.65   27.70 -0.1085E+02  0.6136E+01
   40    2   12    1.66   12.52 -0.4905E+01  0.7933E+01
   36   26    7    1.66   11.99 -0.4697E+01  0.8764E+01
   37   26    5    1.67   34.70 -0.1359E+02  0.8561E+01
   33    4   16    1.68   25.73 -0.1008E+02  0.6693E+01
   38   26    0    1.68   10.15 -0.8096E+01  0.1596E+02
   43    1    8    1.69   14.79 -0.5794E+01  0.1373E+02
   35   25    8    1.69   35.70 -0.1398E+02  0.1304E+02
   28   23   14    1.69   17.87 -0.6999E+01  0.7412E+01
   42    1    9    1.69   45.48 -0.1781E+02  0.8140E+01
   32   24   11    1.70   10.82 -0.4240E+01  0.8270E+01
   27   22   14    1.73   10.69 -0.1087E+02  0.7191E+01
   36    2   13    1.74   10.88 -0.1107E+02  0.7298E+01
   32   23   10    1.75   11.97 -0.1218E+02  0.1468E+02
   39    1   10    1.77   12.55 -0.1277E+02  0.1304E+02
   36   24    3    1.77   12.98 -0.1320E+02  0.9444E+01
   36   24    1    1.79   42.23  0.1508E+02  0.8600E+01
   36   24    0    1.79   14.73  0.1460E+02  0.1245E+02
   40    0    8    1.80   14.72  0.1458E+02  0.8131E+01
   43    1    1    1.80   31.46  0.1123E+02  0.8368E+01
   37    1   11    1.80   15.19  0.5425E+01  0.7696E+01
   27   21   13    1.80   21.02  0.7505E+01  0.7519E+01
   21   19   16    1.81   14.53  0.5189E+01  0.9404E+01
   34   23    5    1.83   10.71  0.3824E+01  0.8926E+01
   42    1    2    1.83   14.06  0.5021E+01  0.8891E+01
   31   22    9    1.84   13.51  0.4825E+01  0.9707E+01
   15   10   19    1.84   41.60  0.1485E+02  0.6907E+01
   29   21   11    1.84   24.77  0.8843E+01  0.7816E+01
   35   23    1    1.85   12.50  0.4464E+01  0.8914E+01
   26   20   13    1.85   17.74  0.6335E+01  0.7453E+01
   31    2   14    1.85   17.01  0.6075E+01  0.7132E+01
   41    1    3    1.87   37.48  0.1338E+02  0.1046E+02
   27   20   12    1.87   25.14  0.8979E+01  0.7504E+01
   40    1    5    1.88   34.86  0.1245E+02  0.1163E+02
   33   22    5    1.89   16.74  0.6460E+01  0.1436E+02
   31   21    8    1.90   17.53  0.6764E+01  0.8388E+01
   34   22    2    1.90   10.22  0.3943E+01  0.8636E+01
   38    1    7    1.91   11.55  0.4456E+01  0.8278E+01
   29   20   10    1.91   32.63  0.1259E+02  0.1380E+02
   26   19   12    1.93   26.69  0.1030E+02  0.7872E+01
   40    2    1    1.93   13.49  0.5206E+01  0.1006E+02
   40    2    0    1.93   11.29  0.9540E+01  0.1319E+02
   14   10   18    1.94   12.60  0.4863E+01  0.6070E+01
   38    1    6    1.94   11.22  0.4330E+01  0.8966E+01
   27   19   11    1.95   29.79  0.1149E+02  0.8550E+01
   39    2    3    1.96   10.35  0.3992E+01  0.1030E+02
   39    2    1    1.98   13.49  0.5206E+01  0.8653E+01
   25   18   12    1.98   35.11  0.1355E+02  0.1228E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        32
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1452
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1354
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       98


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
