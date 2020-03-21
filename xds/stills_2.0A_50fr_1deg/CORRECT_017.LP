
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.414    77.414    38.751  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      17      17
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      17
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000834  0.000674
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
 REFLECTING_RANGE_E.S.D.=   0.099
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

    1440 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      17 ...      17
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1440 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               45
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

 *  44        aP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  1  0  0  0  0 -1  0  0
 *  35        mP          0.0      77.4   38.8   77.4  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  25        mC          0.0     109.5  109.5   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  34        mP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  33        mP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  23        oC          0.0     109.5  109.5   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  32        oP          0.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  20        mC          0.0     109.5  109.5   38.8  90.0  90.0  90.0    1 -1  0  0  1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.4   77.4   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.6   38.8   77.4  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.6   38.8   77.4  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  159.6   77.4  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.6   77.4  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  159.6   77.4  90.0  90.0  75.9    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  159.6   77.4  90.0  90.0  75.9    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     159.6   38.8  109.5  90.0 133.3  75.9   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.5  116.1  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  159.6  159.6  86.6 104.1 104.1    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.5  116.1   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        997.7     173.1   77.4   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        997.7      77.4  173.1   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        997.7      77.4  173.1   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        997.7      77.4   77.4   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.1  116.1  116.1  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.6   86.6  116.1 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.6   86.6  116.1 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.5   86.6   86.6  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.6   86.6  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.6   86.6  109.5  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.6  247.9 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.6   86.6   77.4  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.6   86.6   77.4  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.4  177.4  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.6   86.6  177.4 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.6   86.6   77.4  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.6   86.6   86.6  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.4  109.5   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.4  177.4 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
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

  *   0.00        0    6.1      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.771   -23.468     0.027
 UNIT_CELL_B-AXIS=    23.468    73.771     0.020
 UNIT_CELL_C-AXIS=    -0.016    -0.005    38.751
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       273 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.23
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.361    77.361    38.717  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.1E-01 1.1E-01 5.4E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012926  0.012926  0.025829  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.728   -23.429     0.027
 COORDINATES OF UNIT CELL B-AXIS    23.429    73.728     0.011
 COORDINATES OF UNIT CELL C-AXIS    -0.014    -0.001    38.717
 CRYSTAL MOSAICITY (DEGREES)     0.099
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999976  0.006600  0.001927
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.91
 DETECTOR ORIGIN (PIXELS) AT                     1261.07   1242.50
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.15
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.36    77.36    38.72  90.000  90.000  90.000


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
  0.05  0.15       104        157.516       1.161
  0.15  0.25       162        187.724       0.974
  0.25  0.35       195        162.071       1.129
  0.35  0.45       149         97.383       1.878
  0.45  0.55       168        277.191       0.660
  0.55  0.65       210        180.753       1.012
  0.65  0.75       176        154.119       1.187
  0.75  0.85       166        158.586       1.153
  0.85  0.95       110        295.650       0.619
 
 -0.05  0.95      1440        182.920       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               238.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.875
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    16.1 XMAX=    16.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59590 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        522


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               238.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.870
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     9.1 XMAX=  2456.4 NXBIN=    3
 YMIN=     6.8 YMAX=  2517.4 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        522


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               238.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.862
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    16.1 XMAX=    16.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1262
 DETECTOR_SURFACE_POSITION=    1648    1689
 DETECTOR_SURFACE_POSITION=     817    1689
 DETECTOR_SURFACE_POSITION=     817     836
 DETECTOR_SURFACE_POSITION=    1648     836
 DETECTOR_SURFACE_POSITION=    2173    1662
 DETECTOR_SURFACE_POSITION=    1622    2227
 DETECTOR_SURFACE_POSITION=     843    2227
 DETECTOR_SURFACE_POSITION=     292    1662
 DETECTOR_SURFACE_POSITION=     292     862
 DETECTOR_SURFACE_POSITION=     843     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2173     862
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
 8.638E-01  1.929E-03   24.50


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

   24.464   3.151    20.74   1.66      2.43      2.42      66      99      30
    3.151   2.237    11.49   1.14      4.51      3.98     178     224      14
    2.237   1.829     2.47   0.80     16.21     17.13     210     249       2
    1.829   1.585     0.01   0.80  10956.65  12245.74     242     296       0
    1.585   1.419     0.01   0.80   3183.62   3639.32     276     346       0
    1.419   1.295    -0.09   0.55    -99.90    -99.90     120     177       2
   --------------------------------------------------------------------------
   24.464   1.295     3.76   0.88      5.99      6.27    1092    1391      48
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      17      17)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.682  -0.5683E+01  0.1116E+02        -0.51           1 
   47    0    0     1.646  -0.6621E+01  0.1138E+02        -0.58           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -3.5% OF MEAN INTENSITY


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

     3.84          62      41      1285        3.2%       2.4%      2.4%       42   26.09      3.4%    99.9*     0    0.000       0
     2.73          99      64      2098        3.1%       3.2%      2.9%       70   20.29      4.5%    99.7*     0    0.000       0
     2.24         162      96      2665        3.6%       5.0%      4.5%      132   13.51      7.1%    99.3*     0    0.000       0
     1.94         160      91      3101        2.9%      11.6%     11.9%      138    5.01     16.3%    97.7*     0    0.000       0
     1.74         175     104      3481        3.0%    5291.6%   5574.1%      142   -0.02   7483.4%     6.7      0    0.000       0
     1.59         210     124      3848        3.2%    2754.3%   3251.7%      172    0.03   3895.2%    16.4      0    0.000       0
     1.47         225     136      4143        3.3%    5297.3%   6154.6%      178    0.06   7491.6%     7.4      0    0.000       0
     1.37         214     128      4429        2.9%    1883.4%   2232.6%      172   -0.10   2663.6%    15.2      0    0.000       0
     1.30          82      60      4743        1.3%     -99.9%    -99.9%       44    0.00    -99.9%    36.6      0    0.000       0
    total        1389     844     29793        2.8%       6.0%      6.3%     1090    4.88      8.4%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1440
 NUMBER OF REJECTED MISFITS                                48
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1391
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    845



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

   WILSON LINE (using all data) : A=   7.805 B=  37.057 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      17     7.255   0.005  2.0292E+03       7.615      19.9
      20     4.648   0.012  2.5701E+03       7.852      -2.0
      27     3.514   0.020  1.0611E+03       6.967      20.7
      31     2.991   0.028  4.3789E+02       6.082      30.8
      40     2.649   0.036  3.8412E+02       5.951      26.0
      52     2.397   0.043  2.1357E+02       5.364      28.1
      39     2.206   0.051  1.2388E+02       4.819      29.1
      45     2.057   0.059  4.9564E+01       3.903      33.0
      44     1.932   0.067  7.7304E-01      -0.257      60.2
      57     1.829   0.075 -3.8054E-01    -999.999   -1000.0
      52     1.740   0.083  2.3970E-01      -1.428      55.9
      53     1.665   0.090  9.5014E-01      -0.051      43.6
      55     1.598   0.098  3.3160E-01      -1.104      45.5
      70     1.536   0.106  1.2033E+00       0.185      35.9
      76     1.480   0.114  9.0449E-01      -0.100      34.7
      57     1.433   0.122 -1.4653E+00    -999.999   -1000.0
      50     1.390   0.129 -1.9153E+00    -999.999   -1000.0
      41     1.351   0.137 -9.2864E-01    -999.999   -1000.0
      18     1.315   0.145  1.4726E+00       0.387      25.6


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.255        1.850        1.238        0.658
       3     4.648        1.143        0.665        0.287
       3     3.514        0.304        0.071        0.012
       4     2.991        1.064        0.458        0.148
       3     2.649        0.329        0.082        0.015
       5     2.397        0.924        0.485        0.196
       2     2.206        1.733        0.809        0.273
       8     2.057        0.435        0.183        0.059
      11     1.829       26.699        9.088      117.279
      11     1.665       22.488        8.277       87.311
       4     1.598      613.933     4720.003    57521.461
       3     1.536        8.697        7.694        9.771
      88   overall       75.826      -34.230     4438.039


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.255        0.446        0.308        0.171
      17     4.648        1.056        1.142        1.044
      24     3.514        0.870        0.688        0.511
      27     2.991        0.987        1.020        0.939
      37     2.649        1.176        1.360        1.442
      47     2.397        0.959        0.835        0.611
      37     2.206        0.795        0.611        0.414
      37     2.057        1.334        1.721        2.034
      37     1.932       46.294       42.405      934.197
      43     1.740      597.678    -3656.156   174152.007
      42     1.665       38.563      -38.711     1038.411
      67     1.536       28.168       50.514      383.665
      70     1.480       39.329       25.601      729.042
      18     1.315       49.953      275.607     2587.261
     756   overall       61.783     -219.769    11562.995

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.255  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.32  0.30  0.29
       3     4.648  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       3     3.514  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.991  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.649  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       5     2.397  0.00  0.00  0.00  0.42  0.38  0.71  1.00  0.95  0.91  0.88
       2     2.206  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.057  2.01  1.45  1.20  1.32  1.20  1.11  1.05  0.99  0.95  0.92
      11     1.829  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.72  0.69  0.67
      11     1.665  1.10  0.79  0.65  0.58  0.52  0.49  0.46  0.43  0.42  0.40
       4     1.598  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       3     1.536  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      88   overall  0.82  0.59  0.49  0.48  0.44  0.45  0.44  0.43  0.43  0.42


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.255  2.63  2.76  1.93  1.52  1.91  1.66  1.66  1.51  1.40  1.32
      17     4.648  1.24  1.62  1.13  0.89  1.35  1.17  1.05  1.18  1.19  1.12
      24     3.514  0.88  0.69  0.64  0.63  0.74  0.74  0.83  0.91  0.84  0.86
      27     2.991  1.56  1.84  1.29  1.35  1.22  1.23  1.32  1.21  1.19  1.11
      37     2.649  0.85  1.34  1.36  1.23  1.03  0.96  1.02  1.08  1.05  1.07
      47     2.397  1.57  1.06  1.07  1.23  1.08  1.13  1.10  1.12  1.11  1.08
      37     2.206  1.70  1.04  0.83  0.90  0.89  0.90  0.91  1.03  1.00  1.07
      37     2.057  0.85  0.75  0.94  0.98  0.96  1.08  1.02  1.03  1.14  1.07
      37     1.932  4.26  2.24  1.56  1.23  1.03  0.90  0.81  0.74  0.73  0.77
      43     1.740  3.67  1.92  1.35  1.06  0.89  0.77  0.69  0.63  0.59  0.55
      42     1.665  4.75  2.50  1.84  1.44  1.21  1.06  0.95  0.86  0.80  0.75
      67     1.536  5.18  2.72  1.96  1.54  1.29  1.19  1.07  0.98  0.93  0.90
      70     1.480  4.65  2.44  1.71  1.39  1.16  1.01  0.91  0.83  0.77  0.72
      18     1.315  5.84  3.06  2.14  1.69  1.55  1.35  1.21  1.11  1.03  0.97
     756   overall  2.13  1.28  0.98  0.85  0.76  0.71  0.68  0.66  0.64  0.62


 List of    106 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   32   23   10    1.75   52.64  0.1262E+02  0.1340E+02 "alien"
   19   18   18    1.74   52.32  0.1254E+02  0.6131E+01 "alien"
   39    1   11    1.73   51.53  0.1235E+02  0.1139E+02 "alien"
   42    1    2    1.83   50.20 -0.1910E+02  0.8190E+01 "alien"
   34    3   15    1.70   49.89  0.1196E+02  0.1100E+02 "alien"
   33   23    8    1.79   47.45 -0.1806E+02  0.7603E+01 "alien"
   23   21   16    1.73   46.10  0.1105E+02  0.1180E+02 "alien"
   41    3   13    1.59   39.42  0.1307E+02  0.8096E+01 "alien"
   28   24   15    1.63   37.15  0.1232E+02  0.7449E+01 "alien"
   31    4   16    1.73   36.13  0.8660E+01  0.1099E+02 "alien"
   34   24    8    1.74   34.62  0.8299E+01  0.7870E+01 "alien"
   41    1    4    1.85   34.52 -0.1314E+02  0.8015E+01 "alien"
   32   23    9    1.79   34.06 -0.1296E+02  0.7900E+01 "alien"
   39    4   15    1.57   33.83  0.1122E+02  0.1056E+02 "alien"
   30    3   15    1.82   32.67 -0.1243E+02  0.1049E+02 "alien"
   53    7   13    1.30   31.80  0.4683E+02  0.1276E+02 "alien"
   33    3   15    1.73   30.81  0.7385E+01  0.6747E+01 "alien"
   26   22   15    1.71   30.73  0.7366E+01  0.7055E+01 "alien"
   37   25    2    1.73   30.48  0.7306E+01  0.8912E+01 "alien"
   38   28    6    1.59   30.41  0.1008E+02  0.1221E+02 "alien"
   37   27    7    1.62   30.14  0.9993E+01  0.8294E+01 "alien"
   47    1    6    1.59   29.16  0.9671E+01  0.1212E+02 "alien"
   38    2   12    1.72   27.63  0.6622E+01  0.1240E+02 "alien"
   34   22    4    1.87   27.04 -0.1029E+02  0.9102E+01 "alien"
   24   19   14    1.87   25.80 -0.9817E+01  0.1049E+02 "alien"
   52    5   11    1.36   25.66 -0.2383E+02  0.1382E+02 "alien"
   19   16   19    1.71   25.52  0.6116E+01  0.5305E+01 "alien"
   34   24    7    1.76   25.28  0.6059E+01  0.8446E+01 "alien"
   31   25   13    1.63   25.22  0.8363E+01  0.7568E+01 "alien"
   35   23    2    1.84   25.00 -0.9515E+01  0.7967E+01 "alien"
   35    5   17    1.58   24.53  0.8134E+01  0.1053E+02 "alien"
   36   25    4    1.74   24.22  0.5805E+01  0.8078E+01 "alien"
   33   25   10    1.68   23.96  0.2276E+02  0.1225E+02 "alien"
   52    3    5    1.46   23.09  0.2089E+02  0.1053E+02 "alien"
   29   20   10    1.91   23.03  0.1780E+02  0.1238E+02 "alien"
   36   27    8    1.62   22.33  0.7404E+01  0.8238E+01 "alien"
   30   25   21    1.35   21.88 -0.2032E+02  0.1091E+02 "alien"
   26   21   14    1.78   21.41  0.5133E+01  0.6830E+01 "alien"
   38    1   11    1.76   21.32  0.5111E+01  0.7390E+01 "alien"
   25   20   14    1.82   20.62 -0.7849E+01  0.6995E+01 "alien"
   47    3   10    1.51   20.02  0.2409E+02  0.8621E+01 "alien"
   31   26   14    1.57   19.87  0.6589E+01  0.7965E+01 "alien"
   41    1    3    1.87   19.25 -0.7326E+01  0.9689E+01 "alien"
   27   22   14    1.73   19.05  0.4565E+01  0.6863E+01 "alien"
   45    1    8    1.62   18.64  0.6181E+01  0.8698E+01 "alien"
   39   28    4    1.59   18.52  0.6141E+01  0.1078E+02 "alien"
   35   24    5    1.77   18.43  0.4419E+01  0.8087E+01 "alien"
   22   22   18    1.63   18.32  0.2215E+02  0.1215E+02 "alien"
   34   22    1    1.91   17.91  0.1384E+02  0.8596E+01 "alien"
   16    9   20    1.76   17.71  0.4245E+01  0.5407E+01 "alien"
   47    2    9    1.54   17.66  0.2125E+02  0.1390E+02 "alien"
   33   29   14    1.49   17.29  0.1564E+02  0.9018E+01 "alien"
   53    6   11    1.34   16.98 -0.1577E+02  0.9130E+01 "alien"
   41   32    6    1.45   16.77 -0.2457E+02  0.9203E+01 "alien"
   17   11   20    1.73   16.72  0.4008E+01  0.5857E+01 "alien"
   51    6   12    1.37   16.43 -0.1526E+02  0.1417E+02 "alien"
   33   24    9    1.73   15.85  0.3798E+01  0.7825E+01 "alien"
   47    1    5    1.61   15.48  0.5133E+01  0.8088E+01 "alien"
   51    2    5    1.49   15.37  0.1390E+02  0.8703E+01 "alien"
   37    1    7    1.96   14.74  0.1139E+02  0.8092E+01 "alien"
   27   24   16    1.60   14.70  0.4875E+01  0.1182E+02 "alien"
   39    0    8    1.84   14.15 -0.1434E+02  0.1235E+02 "alien"
   41   32    5    1.46   13.91  0.1258E+02  0.9047E+01 "alien"
   40   29    0    1.57   13.89  0.1205E+02  0.1270E+02 "alien"
   25    9   21    1.56   13.44  0.1618E+02  0.9008E+01 "alien"
   40   32    7    1.46   13.42 -0.1966E+02  0.9364E+01 "alien"
   39   28    3    1.60   13.30  0.4409E+01  0.8511E+01 "alien"
   49    7   15    1.34   13.28 -0.1233E+02  0.1409E+02 "alien"
   39   29    4    1.57   13.11  0.4348E+01  0.8920E+01 "alien"
   17   15   18    1.82   13.03 -0.4958E+01  0.5090E+01 "alien"
   39   31    9    1.46   12.97  0.1173E+02  0.9918E+01 "alien"
   53    2    1    1.46   12.95  0.1172E+02  0.1059E+02 "alien"
   33    6   18    1.57   12.90  0.4279E+01  0.6625E+01 "alien"
   23   17   21    1.52   12.78  0.1538E+02  0.9373E+01 "alien"
   40   29    1    1.56   12.75  0.1534E+02  0.8619E+01 "alien"
   33   28   13    1.53   12.74  0.1533E+02  0.1468E+02 "alien"
   42    3   13    1.56   12.74  0.1533E+02  0.1375E+02 "alien"
   38    5   16    1.55   12.68  0.1526E+02  0.7692E+01 "alien"
   50    2    4    1.53   12.58  0.1514E+02  0.1576E+02 "alien"
   42    0    7    1.75   12.55  0.7166E+01  0.7607E+01 "alien"
   41   11   21    1.30   12.42  0.1830E+02  0.1220E+02 "alien"
   21   18   15    1.90   12.40  0.9588E+01  0.6841E+01 "alien"
   41   31    4    1.49   12.26  0.1109E+02  0.9064E+01 "alien"
   37   31   11    1.46   12.03  0.1088E+02  0.9746E+01 "alien"
   18   18   17    1.82   11.81 -0.1018E+02  0.7751E+01 "alien"
   31   29   19    1.36   11.68 -0.1085E+02  0.8450E+01 "alien"
   29   27   16    1.52   11.64  0.1400E+02  0.1355E+02 "alien"
   18    9   21    1.66   11.57  0.1099E+02  0.6312E+01 "alien"
   34   25    9    1.69   11.50  0.1093E+02  0.7907E+01 "alien"
   25   24   17    1.59   11.46  0.3801E+01  0.7844E+01 "alien"
   33   25   11    1.65   11.45  0.1088E+02  0.1209E+02 "alien"
   32   30   19    1.33   11.41 -0.1060E+02  0.1277E+02 "alien"
   26   19   12    1.93   11.37  0.8793E+01  0.7352E+01 "alien"
   28   27   17    1.50   11.31  0.1023E+02  0.7858E+01 "alien"
   27   25   19    1.46   11.13  0.1007E+02  0.1111E+02 "alien"
   51    4    9    1.43   10.85 -0.1590E+02  0.8722E+01 "alien"
   53    3    3    1.45   10.81 -0.1584E+02  0.1556E+02 "alien"
   52    3    7    1.43   10.79 -0.1582E+02  0.9243E+01 "alien"
   42   32    1    1.46   10.64  0.9627E+01  0.9551E+01 "alien"
   38    1    6    1.94   10.64  0.8227E+01  0.8567E+01 "alien"
   15    7   19    1.87   10.62 -0.4042E+01  0.5533E+01 "alien"
   36    2   13    1.74   10.56  0.2531E+01  0.6940E+01 "alien"
   43    1    2    1.79   10.53 -0.4009E+01  0.7956E+01 "alien"
   36    9   20    1.42   10.37 -0.1519E+02  0.1079E+02 "alien"
   42   32    0    1.47   10.26  0.1904E+02  0.1343E+02 "alien"
   27   18   23    1.38   10.17 -0.1949E+02  0.9126E+01 "alien"


 List of    106 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   41   11   21    1.30   12.42  0.1830E+02  0.1220E+02
   53    7   13    1.30   31.80  0.4683E+02  0.1276E+02
   32   30   19    1.33   11.41 -0.1060E+02  0.1277E+02
   49    7   15    1.34   13.28 -0.1233E+02  0.1409E+02
   53    6   11    1.34   16.98 -0.1577E+02  0.9130E+01
   30   25   21    1.35   21.88 -0.2032E+02  0.1091E+02
   31   29   19    1.36   11.68 -0.1085E+02  0.8450E+01
   52    5   11    1.36   25.66 -0.2383E+02  0.1382E+02
   51    6   12    1.37   16.43 -0.1526E+02  0.1417E+02
   27   18   23    1.38   10.17 -0.1949E+02  0.9126E+01
   36    9   20    1.42   10.37 -0.1519E+02  0.1079E+02
   51    4    9    1.43   10.85 -0.1590E+02  0.8722E+01
   52    3    7    1.43   10.79 -0.1582E+02  0.9243E+01
   53    3    3    1.45   10.81 -0.1584E+02  0.1556E+02
   41   32    6    1.45   16.77 -0.2457E+02  0.9203E+01
   40   32    7    1.46   13.42 -0.1966E+02  0.9364E+01
   53    2    1    1.46   12.95  0.1172E+02  0.1059E+02
   37   31   11    1.46   12.03  0.1088E+02  0.9746E+01
   52    3    5    1.46   23.09  0.2089E+02  0.1053E+02
   39   31    9    1.46   12.97  0.1173E+02  0.9918E+01
   41   32    5    1.46   13.91  0.1258E+02  0.9047E+01
   27   25   19    1.46   11.13  0.1007E+02  0.1111E+02
   42   32    1    1.46   10.64  0.9627E+01  0.9551E+01
   42   32    0    1.47   10.26  0.1904E+02  0.1343E+02
   33   29   14    1.49   17.29  0.1564E+02  0.9018E+01
   41   31    4    1.49   12.26  0.1109E+02  0.9064E+01
   51    2    5    1.49   15.37  0.1390E+02  0.8703E+01
   28   27   17    1.50   11.31  0.1023E+02  0.7858E+01
   47    3   10    1.51   20.02  0.2409E+02  0.8621E+01
   29   27   16    1.52   11.64  0.1400E+02  0.1355E+02
   23   17   21    1.52   12.78  0.1538E+02  0.9373E+01
   50    2    4    1.53   12.58  0.1514E+02  0.1576E+02
   33   28   13    1.53   12.74  0.1533E+02  0.1468E+02
   47    2    9    1.54   17.66  0.2125E+02  0.1390E+02
   38    5   16    1.55   12.68  0.1526E+02  0.7692E+01
   25    9   21    1.56   13.44  0.1618E+02  0.9008E+01
   42    3   13    1.56   12.74  0.1533E+02  0.1375E+02
   40   29    1    1.56   12.75  0.1534E+02  0.8619E+01
   40   29    0    1.57   13.89  0.1205E+02  0.1270E+02
   39    4   15    1.57   33.83  0.1122E+02  0.1056E+02
   39   29    4    1.57   13.11  0.4348E+01  0.8920E+01
   31   26   14    1.57   19.87  0.6589E+01  0.7965E+01
   33    6   18    1.57   12.90  0.4279E+01  0.6625E+01
   35    5   17    1.58   24.53  0.8134E+01  0.1053E+02
   38   28    6    1.59   30.41  0.1008E+02  0.1221E+02
   39   28    4    1.59   18.52  0.6141E+01  0.1078E+02
   41    3   13    1.59   39.42  0.1307E+02  0.8096E+01
   25   24   17    1.59   11.46  0.3801E+01  0.7844E+01
   47    1    6    1.59   29.16  0.9671E+01  0.1212E+02
   39   28    3    1.60   13.30  0.4409E+01  0.8511E+01
   27   24   16    1.60   14.70  0.4875E+01  0.1182E+02
   47    1    5    1.61   15.48  0.5133E+01  0.8088E+01
   37   27    7    1.62   30.14  0.9993E+01  0.8294E+01
   45    1    8    1.62   18.64  0.6181E+01  0.8698E+01
   36   27    8    1.62   22.33  0.7404E+01  0.8238E+01
   22   22   18    1.63   18.32  0.2215E+02  0.1215E+02
   31   25   13    1.63   25.22  0.8363E+01  0.7568E+01
   28   24   15    1.63   37.15  0.1232E+02  0.7449E+01
   33   25   11    1.65   11.45  0.1088E+02  0.1209E+02
   18    9   21    1.66   11.57  0.1099E+02  0.6312E+01
   33   25   10    1.68   23.96  0.2276E+02  0.1225E+02
   34   25    9    1.69   11.50  0.1093E+02  0.7907E+01
   34    3   15    1.70   49.89  0.1196E+02  0.1100E+02
   26   22   15    1.71   30.73  0.7366E+01  0.7055E+01
   19   16   19    1.71   25.52  0.6116E+01  0.5305E+01
   38    2   12    1.72   27.63  0.6622E+01  0.1240E+02
   37   25    2    1.73   30.48  0.7306E+01  0.8912E+01
   17   11   20    1.73   16.72  0.4008E+01  0.5857E+01
   39    1   11    1.73   51.53  0.1235E+02  0.1139E+02
   31    4   16    1.73   36.13  0.8660E+01  0.1099E+02
   33    3   15    1.73   30.81  0.7385E+01  0.6747E+01
   27   22   14    1.73   19.05  0.4565E+01  0.6863E+01
   23   21   16    1.73   46.10  0.1105E+02  0.1180E+02
   33   24    9    1.73   15.85  0.3798E+01  0.7825E+01
   34   24    8    1.74   34.62  0.8299E+01  0.7870E+01
   36   25    4    1.74   24.22  0.5805E+01  0.8078E+01
   19   18   18    1.74   52.32  0.1254E+02  0.6131E+01
   36    2   13    1.74   10.56  0.2531E+01  0.6940E+01
   42    0    7    1.75   12.55  0.7166E+01  0.7607E+01
   32   23   10    1.75   52.64  0.1262E+02  0.1340E+02
   16    9   20    1.76   17.71  0.4245E+01  0.5407E+01
   38    1   11    1.76   21.32  0.5111E+01  0.7390E+01
   34   24    7    1.76   25.28  0.6059E+01  0.8446E+01
   35   24    5    1.77   18.43  0.4419E+01  0.8087E+01
   26   21   14    1.78   21.41  0.5133E+01  0.6830E+01
   32   23    9    1.79   34.06 -0.1296E+02  0.7900E+01
   33   23    8    1.79   47.45 -0.1806E+02  0.7603E+01
   43    1    2    1.79   10.53 -0.4009E+01  0.7956E+01
   17   15   18    1.82   13.03 -0.4958E+01  0.5090E+01
   25   20   14    1.82   20.62 -0.7849E+01  0.6995E+01
   30    3   15    1.82   32.67 -0.1243E+02  0.1049E+02
   18   18   17    1.82   11.81 -0.1018E+02  0.7751E+01
   42    1    2    1.83   50.20 -0.1910E+02  0.8190E+01
   39    0    8    1.84   14.15 -0.1434E+02  0.1235E+02
   35   23    2    1.84   25.00 -0.9515E+01  0.7967E+01
   41    1    4    1.85   34.52 -0.1314E+02  0.8015E+01
   24   19   14    1.87   25.80 -0.9817E+01  0.1049E+02
   41    1    3    1.87   19.25 -0.7326E+01  0.9689E+01
   15    7   19    1.87   10.62 -0.4042E+01  0.5533E+01
   34   22    4    1.87   27.04 -0.1029E+02  0.9102E+01
   21   18   15    1.90   12.40  0.9588E+01  0.6841E+01
   34   22    1    1.91   17.91  0.1384E+02  0.8596E+01
   29   20   10    1.91   23.03  0.1780E+02  0.1238E+02
   26   19   12    1.93   11.37  0.8793E+01  0.7352E+01
   38    1    6    1.94   10.64  0.8227E+01  0.8567E+01
   37    1    7    1.96   14.74  0.1139E+02  0.8092E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        42
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1440
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1326
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      114


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.2 sec
