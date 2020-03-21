
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.389    77.389    38.730  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      14      14
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      14
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000790  0.000697
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
 BEAM_DIVERGENCE_E.S.D.=   0.042
 REFLECTING_RANGE_E.S.D.=   0.097
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

    1455 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      14 ...      14
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1455 REFLECTIONS ACCEPTED



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

 *  44        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.4   38.7   77.4  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  33        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  32        oP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.4   77.4   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.5   38.7   77.4  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.5   38.7   77.4  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    28        mC        250.0      38.7  159.5   77.4  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.7  159.5   77.4  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    36        oC        250.0      38.7  159.5   77.4  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.7  159.5   77.4  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    27        mC        500.0     159.5   38.7  109.4  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.7  109.4  116.1  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.7  159.5  159.5  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.4  116.1   38.7  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        998.1     173.0   77.4   38.7  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.1      77.4  173.0   38.7  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.1      77.4  173.0   38.7  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.1      77.4   77.4   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.1  116.1  116.1  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.5   86.5  116.1 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.5   86.5  116.1 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.4   86.5   86.5  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.5  247.8 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.5   86.5   77.4  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
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

  *   0.00        0    6.0      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.756   -23.429     0.027
 UNIT_CELL_B-AXIS=    23.429    73.756     0.015
 UNIT_CELL_C-AXIS=    -0.015    -0.003    38.730
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       297 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.24
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.304    77.304    38.688  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.9E-01 1.9E-01 8.7E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012936  0.012936  0.025848  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.667   -23.434     0.030
 COORDINATES OF UNIT CELL B-AXIS    23.434    73.667     0.010
 COORDINATES OF UNIT CELL C-AXIS    -0.016    -0.000    38.688
 CRYSTAL MOSAICITY (DEGREES)     0.097
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006269  0.002062
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.89   1232.87
 DETECTOR ORIGIN (PIXELS) AT                     1260.89   1241.98
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.05
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.30    77.30    38.69  90.000  90.000  90.000


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
  0.05  0.15       107        171.868       1.190
  0.15  0.25       170        213.280       0.959
  0.25  0.35       199        164.731       1.241
  0.35  0.45       145        151.940       1.346
  0.45  0.55       168        297.557       0.687
  0.55  0.65       198        187.114       1.093
  0.65  0.75       187        188.968       1.082
  0.75  0.85       168        178.004       1.149
  0.85  0.95       113        316.558       0.646
 
 -0.05  0.95      1455        204.475       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               233.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.057
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    13.1 XMAX=    13.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59678 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        508


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               233.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.048
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     7.8 XMAX=  2436.7 NXBIN=    3
 YMIN=     6.1 YMAX=  2517.7 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        508


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               233.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.035
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    13.1 XMAX=    13.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1222    1262
 DETECTOR_SURFACE_POSITION=    1635    1688
 DETECTOR_SURFACE_POSITION=     810    1688
 DETECTOR_SURFACE_POSITION=     810     835
 DETECTOR_SURFACE_POSITION=    1635     835
 DETECTOR_SURFACE_POSITION=    2156    1662
 DETECTOR_SURFACE_POSITION=    1609    2227
 DETECTOR_SURFACE_POSITION=     836    2227
 DETECTOR_SURFACE_POSITION=     289    1662
 DETECTOR_SURFACE_POSITION=     289     862
 DETECTOR_SURFACE_POSITION=     836     297
 DETECTOR_SURFACE_POSITION=    1609     297
 DETECTOR_SURFACE_POSITION=    2156     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        508


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
 1.039E+00  1.294E-03   27.28


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

   24.446   3.396    22.17   1.56      2.18      2.20      58      82      24
    3.396   2.413    14.87   1.54      4.20      3.12     124     166      14
    2.413   1.974     7.03   1.15      7.31      6.81     184     214       4
    1.974   1.711     0.04   0.62   2318.87   2870.92     182     230       4
    1.711   1.531     0.00   0.69   1154.43   1345.12     220     273       0
    1.531   1.398    -0.02   0.82    -99.90    -99.90     240     299       0
    1.398   1.294     0.12   0.75    375.05    410.97      94     144       0
   --------------------------------------------------------------------------
   24.446   1.294     4.13   0.93      5.84      5.82    1102    1408      46
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      14      14)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681  -0.4472E+01  0.1260E+02        -0.35           1 
   47    0    0     1.645   0.2206E+02  0.1229E+02         1.80           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    10.5% OF MEAN INTENSITY


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

     3.84          59      39      1285        3.0%       2.2%      2.2%       40   27.56      3.2%    99.9*     0    0.000       0
     2.73         107      68      2098        3.2%       3.1%      2.7%       78   22.31      4.4%    99.6*     0    0.000       0
     2.24         160      96      2665        3.6%       4.7%      4.0%      128   14.94      6.7%    99.3*     0    0.000       0
     1.94         161      91      3101        2.9%      10.8%     10.1%      140    5.67     15.3%    98.1*     0    0.000       0
     1.73         177     106      3481        3.0%    1539.4%   1914.4%      142    0.07   2177.1%    -3.9      0    0.000       0
     1.58         208     123      3848        3.2%   25190.1%  29090.1%      170   -0.02  35624.2%    -5.8      0    0.000       0
     1.47         233     142      4143        3.4%   50434.5%  52461.6%      182    0.01  71325.2%    -1.0      0    0.000       0
     1.37         218     129      4429        2.9%     -99.9%    -99.9%      178    0.00    -99.9%     1.5      0    0.000       0
     1.29          85      63      4743        1.3%     122.1%    152.6%       44    0.22    172.7%    48.0      0    0.000       0
    total        1408     857     29793        2.9%       5.8%      5.8%     1102    5.32      8.3%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1455
 NUMBER OF REJECTED MISFITS                                46
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1408
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    857



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

   WILSON LINE (using all data) : A=   7.834 B=  34.944 CORRELATION=  0.86
      #      RES      SS        <I>       log(<I>)       BO
      16     7.154   0.005  2.2323E+03       7.711      12.6
      19     4.628   0.012  2.7142E+03       7.906      -3.1
      30     3.494   0.020  1.1206E+03       7.022      19.8
      32     2.981   0.028  5.1359E+02       6.241      28.3
      39     2.650   0.036  4.8124E+02       6.176      23.3
      53     2.394   0.044  2.6450E+02       5.578      25.9
      39     2.205   0.051  1.5629E+02       5.052      27.0
      45     2.055   0.059  6.8744E+01       4.230      30.4
      46     1.930   0.067 -7.5622E-01    -999.999   -1000.0
      57     1.827   0.075  1.9847E+00       0.685      47.7
      51     1.739   0.083  1.0725E-01      -2.233      60.9
      54     1.664   0.090 -7.9760E-01    -999.999   -1000.0
      57     1.595   0.098  9.5619E-01      -0.045      40.1
      73     1.535   0.106  1.0605E-01      -2.244      47.5
      76     1.479   0.114  5.0074E-01      -0.692      37.3
      59     1.431   0.122 -9.8832E-01    -999.999   -1000.0
      50     1.389   0.130  1.9494E-01      -1.635      36.6
      43     1.351   0.137  7.7923E-01      -0.249      29.5
      18     1.314   0.145  7.2677E+00       1.983      20.2


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.154        1.659        1.036        0.511
       3     4.628        1.186        0.701        0.308
       3     3.494        0.338        0.083        0.015
       4     2.981        1.042        0.452        0.149
       3     2.650        0.295        0.070        0.012
       5     2.394        0.819        0.402        0.151
       2     2.205        1.523        0.656        0.202
       8     2.055        0.309        0.115        0.033
      10     1.827        6.062        6.561       12.399
       9     1.739     2706.485     2546.863  1223327.816
       3     1.595       46.938       10.845      239.067
       3     1.535      527.344     4206.063    24003.472
       6     1.479      156.239      845.048     3471.132
       3     1.431       16.613       29.861       39.417
      86   overall      317.720      464.186   129136.997


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.154        0.436        0.286        0.150
      16     4.628        1.122        1.254        1.166
      27     3.494        0.868        0.696        0.525
      28     2.981        0.960        0.943        0.833
      36     2.650        1.157        1.311        1.362
      48     2.394        0.953        0.795        0.552
      37     2.205        0.792        0.599        0.405
      37     2.055        1.379        1.810        2.199
      47     1.827        7.533       13.171       37.322
      42     1.739     1306.184     4562.481   658595.430
      44     1.664       39.638       66.227     1376.710
      54     1.595       42.811       55.251     1023.545
      50     1.389      940.971      555.647   406581.389
      41     1.351      144.528      981.686    14958.774
      18     1.314        1.566        1.044        0.936
     771   overall      395.922     3591.705   347651.385

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.154  0.00  0.00  0.00  0.00  0.00  0.00  0.33  0.32  0.30  0.59
       3     4.628  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.494  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.981  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.650  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       5     2.394  0.00  0.00  0.00  0.42  0.77  1.07  1.00  0.95  0.91  0.88
       2     2.205  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.055  1.01  1.08  1.50  1.32  1.20  1.11  1.05  1.19  1.14  1.10
      10     1.827  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
       9     1.739  1.79  1.29  1.07  0.94  0.85  0.79  0.74  0.71  0.68  0.65
       3     1.595  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       3     1.535  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     1.479  0.67  0.48  0.40  0.35  0.32  0.30  0.28  0.27  0.25  0.24
       3     1.431  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      86   overall  0.94  0.71  0.64  0.59  0.56  0.56  0.55  0.54  0.51  0.55


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.154  2.87  2.51  1.75  1.38  1.85  1.81  1.63  1.49  1.38  1.29
      16     4.628  1.31  2.07  1.45  1.14  1.27  1.25  1.12  1.13  1.26  1.19
      27     3.494  0.78  0.61  0.57  0.67  0.66  0.82  0.88  1.01  0.94  0.88
      28     2.981  1.50  1.77  1.38  1.30  1.27  1.27  1.21  1.17  1.08  1.07
      36     2.650  0.88  1.38  1.29  1.18  0.99  0.99  1.05  1.11  1.03  1.05
      48     2.394  1.53  1.15  1.04  1.20  1.11  1.11  1.08  1.06  1.12  1.05
      37     2.205  1.70  1.04  0.94  0.74  1.03  0.96  0.86  0.88  1.00  0.98
      37     2.055  0.57  1.04  1.36  1.15  1.10  0.96  1.02  0.98  0.96  0.98
      47     1.827  5.37  2.82  1.97  1.61  1.35  1.18  1.06  0.97  0.90  0.84
      42     1.739  6.00  3.15  2.20  1.73  1.45  1.27  1.14  1.04  0.96  0.90
      44     1.664  5.25  2.76  1.93  1.52  1.27  1.16  1.04  0.95  0.88  0.83
      54     1.595  5.64  2.96  2.07  1.69  1.41  1.23  1.10  1.01  0.94  0.88
      50     1.389  5.04  2.65  1.85  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      41     1.351  4.87  2.56  1.88  1.48  1.24  1.08  0.97  0.89  0.86  0.81
      18     1.314  1.75  0.92  1.29  1.01  0.85  0.99  0.88  0.81  0.84  0.88
     771   overall  2.37  1.44  1.11  0.93  0.84  0.78  0.73  0.69  0.67  0.65


 List of    134 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   40   30    3    1.54  177.87  0.1886E+02  0.9680E+01 "alien"
   31   27   15    1.52  161.08  0.1708E+02  0.1454E+02 "alien"
   51    2    2    1.51  157.60  0.1671E+02  0.1040E+02 "alien"
   35   28   11    1.55  154.78  0.1641E+02  0.9102E+01 "alien"
   38   29    7    1.55  150.03  0.1591E+02  0.9342E+01 "alien"
   23   21   16    1.73  118.67  0.1273E+02  0.1308E+02 "alien"
   23   17   21    1.52  108.92  0.1155E+02  0.9995E+01 "alien"
   31    4   16    1.73   97.60  0.1047E+02  0.1130E+02 "alien"
   40   29    1    1.56   97.31  0.1032E+02  0.9514E+01 "alien"
   40   30    2    1.54   91.75  0.9730E+01  0.1067E+02 "alien"
   33    7   19    1.52   91.08  0.9659E+01  0.7238E+01 "alien"
   42    3   13    1.56   90.62  0.9610E+01  0.1539E+02 "alien"
   50    5   11    1.41   90.35  0.1761E+02  0.1854E+02 "alien"
   32   11   22    1.39   85.24  0.1662E+02  0.1326E+02 "alien"
   39   29    5    1.56   82.83  0.8784E+01  0.9498E+01 "alien"
   50    1    2    1.54   79.01  0.8379E+01  0.9419E+01 "alien"
   35    6   18    1.53   77.74  0.8244E+01  0.1258E+02 "alien"
   53    4    7    1.41   73.54  0.1434E+02  0.9883E+01 "alien"
   35   24    6    1.75   73.49  0.7882E+01  0.8735E+01 "alien"
   48    6   14    1.38   72.00  0.1404E+02  0.1391E+02 "alien"
   47    2    9    1.53   71.75  0.7609E+01  0.1504E+02 "alien"
   38    1   11    1.76   71.69  0.7689E+01  0.7868E+01 "alien"
   26   13   24    1.38   70.46  0.1374E+02  0.9245E+01 "alien"
   29   26   20    1.39   69.69  0.1359E+02  0.1504E+02 "alien"
   22   14   22    1.51   67.54  0.7162E+01  0.8212E+01 "alien"
   38    2   12    1.72   63.28  0.6787E+01  0.1287E+02 "alien"
   36   24    3    1.77   62.43  0.6695E+01  0.9772E+01 "alien"
   53    2    1    1.46   60.81  0.3045E+02  0.1172E+02 "alien"
   54    5    9    1.35   54.72  0.4264E+02  0.9325E+01 "alien"
   53    4    8    1.39   51.50  0.1004E+02  0.9792E+01 "alien"
   36   25    5    1.72   51.32  0.5504E+01  0.8730E+01 "alien"
   42    8   18    1.38   48.96  0.9545E+01  0.8907E+01 "alien"
   27    9   21    1.52   46.66  0.4948E+01  0.8330E+01 "alien"
   44    1    1    1.75   46.64  0.5002E+01  0.9079E+01 "alien"
   35   24    5    1.77   44.96  0.4822E+01  0.8669E+01 "alien"
   17   11   20    1.73   44.66  0.4789E+01  0.6131E+01 "alien"
   29   23   13    1.71   43.78  0.4695E+01  0.7708E+01 "alien"
   33    3   15    1.73   43.49  0.4665E+01  0.7171E+01 "alien"
   35   25    7    1.71   41.22  0.4421E+01  0.8593E+01 "alien"
   17   12   23    1.53   41.12  0.4361E+01  0.6525E+01 "alien"
   26   21   14    1.77   39.30  0.4215E+01  0.7512E+01 "alien"
   40   33    9    1.41   38.26  0.7458E+01  0.9891E+01 "alien"
   33   12   23    1.34   37.49  0.2921E+02  0.1412E+02 "alien"
   46    2    9    1.56   35.90  0.3807E+01  0.8720E+01 "alien"
   37   30   10    1.50   34.38  0.1722E+02  0.9938E+01 "alien"
   27   10   22    1.47   33.18  0.1661E+02  0.1036E+02 "alien"
   32    8   20    1.49   33.03  0.1654E+02  0.7821E+01 "alien"
   47    3   10    1.51   32.98  0.3497E+01  0.8831E+01 "alien"
   52    3    5    1.46   32.51  0.1628E+02  0.1112E+02 "alien"
   36   10   21    1.38   32.14  0.6265E+01  0.8000E+01 "alien"
   49    2    7    1.52   32.06  0.3400E+01  0.8890E+01 "alien"
   33   27   13    1.55   31.53  0.3344E+01  0.1271E+02 "alien"
   54    5    8    1.37   30.44  0.2372E+02  0.1632E+02 "alien"
   52    5   11    1.36   30.05  0.2342E+02  0.1564E+02 "alien"
   52    2    3    1.48   27.65  0.1385E+02  0.9408E+01 "alien"
   39   33   11    1.39   27.59  0.5378E+01  0.1029E+02 "alien"
   46    7   16    1.37   27.07  0.5278E+01  0.9183E+01 "alien"
   43    7   17    1.40   26.60  0.5185E+01  0.9093E+01 "alien"
   20   17   19    1.67   26.54 -0.2117E+02  0.1071E+02 "alien"
   47    1    6    1.59   25.24  0.2413E+02  0.1295E+02 "alien"
   52    2    4    1.47   25.06  0.1255E+02  0.1208E+02 "alien"
   39   33   10    1.41   24.92  0.4858E+01  0.1109E+02 "alien"
   37   25    0    1.73   24.82  0.1669E+02  0.1266E+02 "alien"
   19    6   23    1.54   24.44  0.2592E+01  0.5990E+01 "alien"
   29   20   10    1.91   24.02 -0.1817E+02  0.1296E+02 "alien"
   30   28   19    1.38   23.90  0.4658E+01  0.8833E+01 "alien"
   43    2   10    1.63   23.90 -0.1906E+02  0.1336E+02 "alien"
   40    7   18    1.43   23.33 -0.2306E+02  0.1507E+02 "alien"
   39   34   12    1.36   23.20  0.1808E+02  0.1633E+02 "alien"
   39    2    1    1.98   23.12 -0.1748E+02  0.8801E+01 "alien"
   26   10   22    1.49   22.98  0.1151E+02  0.9120E+01 "alien"
   27   18   23    1.37   22.90  0.4464E+01  0.9889E+01 "alien"
   19   18   18    1.74   22.79  0.2445E+01  0.6548E+01 "alien"
   38   33   12    1.39   22.10  0.4308E+01  0.9561E+01 "alien"
   50    6   13    1.36   21.99  0.1714E+02  0.9731E+01 "alien"
   21   12   22    1.54   21.63  0.2294E+01  0.8642E+01 "alien"
   39    1   10    1.76   20.49  0.2198E+01  0.1277E+02 "alien"
   28   23   21    1.39   20.12  0.3922E+01  0.8302E+01 "alien"
   37   31   12    1.43   19.88 -0.1964E+02  0.1473E+02 "alien"
   47    1    5    1.61   19.46  0.1860E+02  0.9092E+01 "alien"
   28   22   13    1.75   19.39  0.2080E+01  0.7788E+01 "alien"
   21   10   23    1.50   19.24  0.9634E+01  0.7770E+01 "alien"
   37    5   16    1.57   19.02  0.1819E+02  0.1454E+02 "alien"
   45    0    3    1.70   19.02  0.7857E+01  0.8446E+01 "alien"
   37   10   21    1.36   18.23  0.1421E+02  0.1218E+02 "alien"
   41   32    5    1.46   18.16  0.9095E+01  0.9694E+01 "alien"
   47    2    8    1.56   18.02  0.1911E+01  0.8548E+01 "alien"
   42    0    6    1.77   17.55  0.6455E+01  0.1515E+02 "alien"
   21   20   18    1.67   17.46 -0.1393E+02  0.1146E+02 "alien"
   43    0    6    1.73   17.28  0.6225E+01  0.8308E+01 "alien"
   52    4    9    1.40   17.06  0.3326E+01  0.9563E+01 "alien"
   34   29   13    1.50   16.78  0.8401E+01  0.9108E+01 "alien"
   41   31    4    1.49   16.44  0.8230E+01  0.9792E+01 "alien"
   50    3    8    1.47   16.32  0.8174E+01  0.9387E+01 "alien"
   33    5   17    1.62   16.32  0.1561E+02  0.7422E+01 "alien"
   37   34   14    1.34   16.23  0.1265E+02  0.1013E+02 "alien"
   34    7   19    1.50   15.97  0.7997E+01  0.1070E+02 "alien"
   39   32   10    1.42   15.40 -0.1522E+02  0.1023E+02 "alien"
   46    5   14    1.43   15.09 -0.1491E+02  0.1359E+02 "alien"
   50    1    0    1.55   14.85  0.4414E+01  0.1278E+02 "alien"
   28    2   14    1.95   14.79 -0.1119E+02  0.1058E+02 "alien"
   25   25   18    1.53   14.58  0.4246E+01  0.7950E+01 "alien"
   29   27   16    1.52   14.48  0.1535E+01  0.1401E+02 "alien"
   36    9   20    1.42   14.18 -0.1401E+02  0.1190E+02 "alien"
   24   24   18    1.56   14.12  0.3983E+01  0.8099E+01 "alien"
   38    7   18    1.46   14.10  0.7061E+01  0.8166E+01 "alien"
   30    1   13    1.95   14.08 -0.1065E+02  0.1112E+02 "alien"
   28    9   21    1.51   13.62  0.1445E+01  0.6337E+01 "alien"
   49    4   11    1.44   13.39 -0.1323E+02  0.9233E+01 "alien"
   39    1    6    1.89   13.35 -0.1009E+02  0.8332E+01 "alien"
   38    8   19    1.42   13.27 -0.1311E+02  0.7950E+01 "alien"
   53    2    0    1.46   13.18  0.1642E+02  0.1518E+02 "alien"
   51    2    3    1.50   12.82  0.6417E+01  0.9109E+01 "alien"
   22   13   24    1.42   12.73 -0.1258E+02  0.7537E+01 "alien"
   40    1    5    1.87   12.72  0.2525E+02  0.1240E+02 "alien"
   45    2    9    1.59   12.51  0.1196E+02  0.9261E+01 "alien"
   52    3    6    1.45   12.23 -0.1209E+02  0.9187E+01 "alien"
   42   32    0    1.46   12.19  0.1419E+02  0.1373E+02 "alien"
   30    9   21    1.48   12.11  0.6063E+01  0.9507E+01 "alien"
   39   28    3    1.60   12.10  0.1157E+02  0.9389E+01 "alien"
   47    4   12    1.46   12.01  0.6014E+01  0.8997E+01 "alien"
   41   32    4    1.47   11.50  0.5758E+01  0.1042E+02 "alien"
   54    3    4    1.41   11.37 -0.1123E+02  0.9509E+01 "alien"
   41    5   15    1.51   11.32  0.1200E+01  0.1316E+02 "alien"
   34   27   11    1.59   10.69  0.1023E+02  0.9391E+01 "alien"
   54    3    1    1.43   10.47 -0.1035E+02  0.9309E+01 "alien"
   42    4   14    1.53   10.45  0.1108E+01  0.8360E+01 "alien"
   41   31    0    1.50   10.39  0.1075E+02  0.1369E+02 "alien"
   33   25   11    1.65   10.30 -0.8218E+01  0.1302E+02 "alien"
   29   22   22    1.35   10.29  0.8015E+01  0.8190E+01 "alien"
   19   14   20    1.67   10.24 -0.8165E+01  0.5442E+01 "alien"
   21   18   15    1.90   10.09 -0.7631E+01  0.7106E+01 "alien"
   24   14   23    1.44   10.04 -0.9924E+01  0.6582E+01 "alien"
   27   27   18    1.47   10.04  0.1016E+02  0.9362E+01 "alien"


 List of    134 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   33   12   23    1.34   37.49  0.2921E+02  0.1412E+02
   37   34   14    1.34   16.23  0.1265E+02  0.1013E+02
   54    5    9    1.35   54.72  0.4264E+02  0.9325E+01
   29   22   22    1.35   10.29  0.8015E+01  0.8190E+01
   39   34   12    1.36   23.20  0.1808E+02  0.1633E+02
   37   10   21    1.36   18.23  0.1421E+02  0.1218E+02
   52    5   11    1.36   30.05  0.2342E+02  0.1564E+02
   50    6   13    1.36   21.99  0.1714E+02  0.9731E+01
   54    5    8    1.37   30.44  0.2372E+02  0.1632E+02
   46    7   16    1.37   27.07  0.5278E+01  0.9183E+01
   27   18   23    1.37   22.90  0.4464E+01  0.9889E+01
   36   10   21    1.38   32.14  0.6265E+01  0.8000E+01
   26   13   24    1.38   70.46  0.1374E+02  0.9245E+01
   30   28   19    1.38   23.90  0.4658E+01  0.8833E+01
   48    6   14    1.38   72.00  0.1404E+02  0.1391E+02
   42    8   18    1.38   48.96  0.9545E+01  0.8907E+01
   29   26   20    1.39   69.69  0.1359E+02  0.1504E+02
   38   33   12    1.39   22.10  0.4308E+01  0.9561E+01
   39   33   11    1.39   27.59  0.5378E+01  0.1029E+02
   53    4    8    1.39   51.50  0.1004E+02  0.9792E+01
   32   11   22    1.39   85.24  0.1662E+02  0.1326E+02
   28   23   21    1.39   20.12  0.3922E+01  0.8302E+01
   43    7   17    1.40   26.60  0.5185E+01  0.9093E+01
   52    4    9    1.40   17.06  0.3326E+01  0.9563E+01
   53    4    7    1.41   73.54  0.1434E+02  0.9883E+01
   40   33    9    1.41   38.26  0.7458E+01  0.9891E+01
   39   33   10    1.41   24.92  0.4858E+01  0.1109E+02
   50    5   11    1.41   90.35  0.1761E+02  0.1854E+02
   54    3    4    1.41   11.37 -0.1123E+02  0.9509E+01
   36    9   20    1.42   14.18 -0.1401E+02  0.1190E+02
   22   13   24    1.42   12.73 -0.1258E+02  0.7537E+01
   38    8   19    1.42   13.27 -0.1311E+02  0.7950E+01
   39   32   10    1.42   15.40 -0.1522E+02  0.1023E+02
   40    7   18    1.43   23.33 -0.2306E+02  0.1507E+02
   54    3    1    1.43   10.47 -0.1035E+02  0.9309E+01
   46    5   14    1.43   15.09 -0.1491E+02  0.1359E+02
   37   31   12    1.43   19.88 -0.1964E+02  0.1473E+02
   49    4   11    1.44   13.39 -0.1323E+02  0.9233E+01
   24   14   23    1.44   10.04 -0.9924E+01  0.6582E+01
   52    3    6    1.45   12.23 -0.1209E+02  0.9187E+01
   53    2    1    1.46   60.81  0.3045E+02  0.1172E+02
   53    2    0    1.46   13.18  0.1642E+02  0.1518E+02
   52    3    5    1.46   32.51  0.1628E+02  0.1112E+02
   41   32    5    1.46   18.16  0.9095E+01  0.9694E+01
   47    4   12    1.46   12.01  0.6014E+01  0.8997E+01
   42   32    0    1.46   12.19  0.1419E+02  0.1373E+02
   38    7   18    1.46   14.10  0.7061E+01  0.8166E+01
   52    2    4    1.47   25.06  0.1255E+02  0.1208E+02
   41   32    4    1.47   11.50  0.5758E+01  0.1042E+02
   50    3    8    1.47   16.32  0.8174E+01  0.9387E+01
   27   10   22    1.47   33.18  0.1661E+02  0.1036E+02
   27   27   18    1.47   10.04  0.1016E+02  0.9362E+01
   52    2    3    1.48   27.65  0.1385E+02  0.9408E+01
   30    9   21    1.48   12.11  0.6063E+01  0.9507E+01
   26   10   22    1.49   22.98  0.1151E+02  0.9120E+01
   41   31    4    1.49   16.44  0.8230E+01  0.9792E+01
   32    8   20    1.49   33.03  0.1654E+02  0.7821E+01
   34   29   13    1.50   16.78  0.8401E+01  0.9108E+01
   37   30   10    1.50   34.38  0.1722E+02  0.9938E+01
   21   10   23    1.50   19.24  0.9634E+01  0.7770E+01
   34    7   19    1.50   15.97  0.7997E+01  0.1070E+02
   41   31    0    1.50   10.39  0.1075E+02  0.1369E+02
   51    2    3    1.50   12.82  0.6417E+01  0.9109E+01
   28    9   21    1.51   13.62  0.1445E+01  0.6337E+01
   51    2    2    1.51  157.60  0.1671E+02  0.1040E+02
   47    3   10    1.51   32.98  0.3497E+01  0.8831E+01
   22   14   22    1.51   67.54  0.7162E+01  0.8212E+01
   41    5   15    1.51   11.32  0.1200E+01  0.1316E+02
   49    2    7    1.52   32.06  0.3400E+01  0.8890E+01
   29   27   16    1.52   14.48  0.1535E+01  0.1401E+02
   31   27   15    1.52  161.08  0.1708E+02  0.1454E+02
   33    7   19    1.52   91.08  0.9659E+01  0.7238E+01
   23   17   21    1.52  108.92  0.1155E+02  0.9995E+01
   27    9   21    1.52   46.66  0.4948E+01  0.8330E+01
   42    4   14    1.53   10.45  0.1108E+01  0.8360E+01
   35    6   18    1.53   77.74  0.8244E+01  0.1258E+02
   17   12   23    1.53   41.12  0.4361E+01  0.6525E+01
   25   25   18    1.53   14.58  0.4246E+01  0.7950E+01
   47    2    9    1.53   71.75  0.7609E+01  0.1504E+02
   40   30    3    1.54  177.87  0.1886E+02  0.9680E+01
   21   12   22    1.54   21.63  0.2294E+01  0.8642E+01
   50    1    2    1.54   79.01  0.8379E+01  0.9419E+01
   40   30    2    1.54   91.75  0.9730E+01  0.1067E+02
   19    6   23    1.54   24.44  0.2592E+01  0.5990E+01
   50    1    0    1.55   14.85  0.4414E+01  0.1278E+02
   33   27   13    1.55   31.53  0.3344E+01  0.1271E+02
   35   28   11    1.55  154.78  0.1641E+02  0.9102E+01
   38   29    7    1.55  150.03  0.1591E+02  0.9342E+01
   47    2    8    1.56   18.02  0.1911E+01  0.8548E+01
   39   29    5    1.56   82.83  0.8784E+01  0.9498E+01
   42    3   13    1.56   90.62  0.9610E+01  0.1539E+02
   24   24   18    1.56   14.12  0.3983E+01  0.8099E+01
   40   29    1    1.56   97.31  0.1032E+02  0.9514E+01
   46    2    9    1.56   35.90  0.3807E+01  0.8720E+01
   37    5   16    1.57   19.02  0.1819E+02  0.1454E+02
   34   27   11    1.59   10.69  0.1023E+02  0.9391E+01
   47    1    6    1.59   25.24  0.2413E+02  0.1295E+02
   45    2    9    1.59   12.51  0.1196E+02  0.9261E+01
   39   28    3    1.60   12.10  0.1157E+02  0.9389E+01
   47    1    5    1.61   19.46  0.1860E+02  0.9092E+01
   33    5   17    1.62   16.32  0.1561E+02  0.7422E+01
   43    2   10    1.63   23.90 -0.1906E+02  0.1336E+02
   33   25   11    1.65   10.30 -0.8218E+01  0.1302E+02
   19   14   20    1.67   10.24 -0.8165E+01  0.5442E+01
   21   20   18    1.67   17.46 -0.1393E+02  0.1146E+02
   20   17   19    1.67   26.54 -0.2117E+02  0.1071E+02
   45    0    3    1.70   19.02  0.7857E+01  0.8446E+01
   35   25    7    1.71   41.22  0.4421E+01  0.8593E+01
   29   23   13    1.71   43.78  0.4695E+01  0.7708E+01
   38    2   12    1.72   63.28  0.6787E+01  0.1287E+02
   36   25    5    1.72   51.32  0.5504E+01  0.8730E+01
   17   11   20    1.73   44.66  0.4789E+01  0.6131E+01
   31    4   16    1.73   97.60  0.1047E+02  0.1130E+02
   33    3   15    1.73   43.49  0.4665E+01  0.7171E+01
   37   25    0    1.73   24.82  0.1669E+02  0.1266E+02
   43    0    6    1.73   17.28  0.6225E+01  0.8308E+01
   23   21   16    1.73  118.67  0.1273E+02  0.1308E+02
   19   18   18    1.74   22.79  0.2445E+01  0.6548E+01
   35   24    6    1.75   73.49  0.7882E+01  0.8735E+01
   28   22   13    1.75   19.39  0.2080E+01  0.7788E+01
   44    1    1    1.75   46.64  0.5002E+01  0.9079E+01
   38    1   11    1.76   71.69  0.7689E+01  0.7868E+01
   39    1   10    1.76   20.49  0.2198E+01  0.1277E+02
   36   24    3    1.77   62.43  0.6695E+01  0.9772E+01
   42    0    6    1.77   17.55  0.6455E+01  0.1515E+02
   35   24    5    1.77   44.96  0.4822E+01  0.8669E+01
   26   21   14    1.77   39.30  0.4215E+01  0.7512E+01
   40    1    5    1.87   12.72  0.2525E+02  0.1240E+02
   39    1    6    1.89   13.35 -0.1009E+02  0.8332E+01
   21   18   15    1.90   10.09 -0.7631E+01  0.7106E+01
   29   20   10    1.91   24.02 -0.1817E+02  0.1296E+02
   30    1   13    1.95   14.08 -0.1065E+02  0.1112E+02
   28    2   14    1.95   14.79 -0.1119E+02  0.1058E+02
   39    2    1    1.98   23.12 -0.1748E+02  0.8801E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        78
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1455
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1287
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      168


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
