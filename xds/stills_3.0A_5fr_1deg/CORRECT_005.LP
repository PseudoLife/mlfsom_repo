
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.349    77.349    38.720  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_3.0A_5fr_1deg/stills_3.0A_5fr_1deg_???_001.img SMV
 DATA_RANGE=       5       5
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       5
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001033  0.000610
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
 REFLECTING_RANGE_E.S.D.=   0.202
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

    1051 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       5 ...       5
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1051 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               33
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
 *  31        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  35        mP          0.0      77.4   38.7   77.4  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  33        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  21        tP          0.0      77.4   77.4   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
 *  32        oP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    37        mC        250.0     159.5   38.7   77.4  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.5   38.7   77.4  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.7  159.5   77.4  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.7  159.5   77.4  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.7  159.5   77.4  90.0  90.0  75.9    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    29        mC        250.0      38.7  159.5   77.4  90.0  90.0  75.9    0  0  1  0  2  0  1  0  0  1  0  0
    19        oI        500.0      38.7  109.4  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    27        mC        500.0     159.5   38.7  109.4  90.0 133.3  75.9    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    18        tI        625.0     109.4  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    26        oF        625.0      38.7  159.5  159.5  86.6 104.1 104.1    0  0  1  0  0  2 -1  0 -2  0 -1  0
    41        mC        997.7     173.0   77.4   38.7  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    22        hP        997.7      77.4   77.4   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    30        mC        997.7      77.4  173.0   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        997.7      77.4  173.0   38.7  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     3        cP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.4   86.5   86.5  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.5  247.6 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.5   86.5   77.4  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.5   86.5   77.4  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.7   77.4  177.2  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.5   86.5  177.2 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.5   86.5   77.4  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.2  109.4   38.7  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.7   77.4  177.2 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.5  177.2   77.4 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    3.1       9    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.719   -23.423     0.032
 UNIT_CELL_B-AXIS=    23.423    73.719     0.011
 UNIT_CELL_C-AXIS=    -0.017    -0.000    38.720
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM        89 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.19
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.316    77.316    38.696  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.5E-01 2.5E-01 1.1E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012934  0.012934  0.025843  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.688   -23.404     0.007
 COORDINATES OF UNIT CELL B-AXIS    23.404    73.688    -0.000
 COORDINATES OF UNIT CELL C-AXIS    -0.003     0.001    38.696
 CRYSTAL MOSAICITY (DEGREES)     0.202
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006208  0.002302
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.83   1232.87
 DETECTOR ORIGIN (PIXELS) AT                     1260.48   1241.90
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.20
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.32    77.32    38.70  90.000  90.000  90.000


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
  0.05  0.15         8          0.173     583.122
  0.15  0.25       115         73.121       1.377
  0.25  0.35       152        116.249       0.866
  0.35  0.45       149         24.904       4.044
  0.45  0.55       166        170.723       0.590
  0.55  0.65       163         78.587       1.282
  0.65  0.75       170        111.584       0.903
  0.75  0.85       123        126.175       0.798
  0.85  0.95         5         85.667       1.176
 
 -0.05  0.95      1051        100.721       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        3
 DEGREES OF FREEDOM OF CHI^2 FIT                68.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.565
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     4.1 XMAX=     4.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59658 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        158


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT                68.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.564
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    63.3 XMAX=  2438.6 NXBIN=    2
 YMIN=     7.3 YMAX=  2518.9 NYBIN=    2
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        158


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT                63.9
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.573
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     4.1 XMAX=     4.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1251    1263
 DETECTOR_SURFACE_POSITION=    1654    1690
 DETECTOR_SURFACE_POSITION=     848    1690
 DETECTOR_SURFACE_POSITION=     848     837
 DETECTOR_SURFACE_POSITION=    1654     837
 DETECTOR_SURFACE_POSITION=    2164    1663
 DETECTOR_SURFACE_POSITION=    1629    2228
 DETECTOR_SURFACE_POSITION=     873    2228
 DETECTOR_SURFACE_POSITION=     338    1663
 DETECTOR_SURFACE_POSITION=     338     863
 DETECTOR_SURFACE_POSITION=     873     298
 DETECTOR_SURFACE_POSITION=    1629     298
 DETECTOR_SURFACE_POSITION=    2164     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        158


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
 7.316E-01  5.965E-04   47.87


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

   24.449   2.879    22.61   1.31      1.67      1.51      50      97      22
    2.879   2.043     0.00   0.86   3397.53   3634.73     148     210       0
    2.043   1.670    -0.02   0.99   2392.37   2439.60     192     255       0
    1.670   1.447    -0.03   1.12    -99.90    -99.90     196     283       2
    1.447   1.295     0.04   0.79    -99.90    -99.90     106     179       2
   --------------------------------------------------------------------------
   24.449   1.295     2.14   0.99      9.13      9.05     692    1024      26
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       5       5)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681  -0.8766E+01  0.1186E+02        -0.74           1 
   47    0    0     1.645  -0.1911E+02  0.1216E+02        -1.57           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS   -18.8% OF MEAN INTENSITY


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

     3.84          50      37      1285        2.9%       1.5%      1.3%       26   36.29      2.2%   100.0*     0    0.000       0
     2.73          66      49      2098        2.3%       6.2%      8.0%       34   11.32      8.8%    99.7*     0    0.000       0
     2.24         126      80      2665        3.0%     489.8%    511.8%       92    0.09    692.6%    -1.5      0    0.000       0
     1.94         117      75      3101        2.4%     -99.9%    -99.9%       84    0.00    -99.9%    11.7      0    0.000       0
     1.74         137      85      3481        2.4%    2086.3%   2401.1%      104   -0.08   2950.4%    13.4      0    0.000       0
     1.58         158     101      3847        2.6%     659.7%    655.1%      114    0.10    933.0%    -3.9      0    0.000       0
     1.47         157     104      4142        2.5%     -99.9%    -99.9%      106    0.00    -99.9%    -9.7      0    0.000       0
     1.37         151      98      4431        2.2%     -99.9%    -99.9%      106    0.00    -99.9%   -14.8      0    0.000       0
     1.29          62      49      4739        1.0%     530.3%    572.9%       26    0.21    750.0%    25.9      0    0.000       0
    total        1024     678     29789        2.3%       9.1%      9.0%      692    2.79     12.9%   100.0*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1051
 NUMBER OF REJECTED MISFITS                                26
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1024
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    678



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

   WILSON LINE (using all data) : A=   5.633 B=  22.373 CORRELATION=  0.74
      #      RES      SS        <I>       log(<I>)       BO
      14     7.639   0.004  2.0437E+03       7.623    -232.2
      19     4.745   0.011  1.4527E+03       7.281     -74.2
      22     3.525   0.020  5.0323E+02       6.221     -14.6
      24     2.972   0.028  6.5939E+01       4.189      25.5
      31     2.639   0.036  3.0028E+00       1.100      63.2
      44     2.400   0.043 -7.9918E-01    -999.999   -1000.0
      34     2.204   0.051  8.6396E-01      -0.146      56.2
      37     2.053   0.059 -6.4650E-01    -999.999   -1000.0
      34     1.930   0.067 -1.4988E+00    -999.999   -1000.0
      49     1.828   0.075 -3.7101E-01    -999.999   -1000.0
      39     1.738   0.083 -7.4724E-01    -999.999   -1000.0
      44     1.665   0.090  9.5310E-01      -0.048      31.5
      47     1.597   0.098  2.3129E+00       0.839      24.4
      50     1.536   0.106 -8.2633E-01    -999.999   -1000.0
      59     1.480   0.114 -1.4956E+00    -999.999   -1000.0
      46     1.432   0.122 -1.1818E+00    -999.999   -1000.0
      37     1.389   0.130  5.3113E-01      -0.633      24.2
      33     1.351   0.137  6.6612E-01      -0.406      22.0
      15     1.315   0.145  1.0524E+01       2.354      11.3


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.639        2.214        1.845        1.196
       3     4.745        1.395        0.871        0.407
       2     3.525        0.597        0.163        0.032
       3     2.972        3.164        3.376        2.574
       2     2.639        2.370        1.605        0.810
       5     2.400       36.629      -12.406      217.214
       5     2.053       32.544      119.629      362.923
      10     1.828       24.072        7.404       70.971
      10     1.665       20.619       37.284       79.120
       3     1.597       18.055       37.140       58.796
       5     1.480       21.665       22.974       59.791
      72   overall       10.283      -27.296      153.376


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.639        0.490        0.341        0.186
      16     4.745        1.442        2.077        2.370
      20     3.525        1.057        1.216        1.284
      21     2.972        2.380        4.830        8.015
      29     2.639        3.929        3.091        7.697
      33     2.204       74.656       54.558     2610.212
      34     1.665       48.242        0.341      728.640
      44     1.597       10.746       16.897       67.171
      37     1.389      177.806     -224.323    12536.552
      31     1.351      178.638      349.487    11763.017
      15     1.315        9.789       49.878      206.671
     606   overall       28.141       10.049     1560.939

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.639  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.38  0.37
       3     4.745  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     3.525  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.972  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       2     2.639  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       5     2.400  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
       5     2.053  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.59
      10     1.828  2.42  1.74  1.44  1.27  1.15  1.07  1.00  0.95  0.91  0.88
      10     1.665  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.95  0.91  0.88
       3     1.597  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       5     1.480  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.59
      72   overall  1.17  0.84  0.70  0.62  0.56  0.52  0.49  0.51  0.49  0.47


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.639  2.10  3.31  2.31  2.12  1.78  1.77  1.59  1.45  1.35  1.27
      16     4.745  2.63  2.07  1.69  1.33  1.59  1.66  1.49  1.36  1.26  1.19
      20     3.525  0.53  1.10  0.77  0.76  1.02  1.00  1.29  1.27  1.18  1.19
      21     2.972  7.51  3.94  2.76  2.17  1.82  1.58  1.42  1.30  1.28  1.21
      29     2.639  3.26  2.09  1.46  1.15  1.05  0.92  0.89  0.81  0.76  0.71
      33     2.204  3.82  2.17  1.52  1.19  1.00  0.87  0.78  0.77  0.71  0.67
      34     1.665  4.33  2.27  1.59  1.25  1.05  0.91  0.82  0.85  0.79  0.74
      44     1.597  5.49  3.01  2.10  1.65  1.44  1.26  1.13  1.03  1.00  0.93
      37     1.389  4.83  2.68  1.88  1.48  1.24  1.08  0.97  0.88  0.82  0.77
      31     1.351  5.42  2.85  1.99  1.57  1.31  1.14  1.03  0.94  0.87  0.82
      15     1.315  5.60  2.94  2.06  1.82  1.69  1.48  1.32  1.33  1.24  1.16
     606   overall  2.10  1.23  0.87  0.70  0.62  0.56  0.51  0.49  0.46  0.43


 List of    113 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   21   19   16    1.81   59.42 -0.2205E+02  0.1021E+02 "alien"
   35   23    3    1.83   56.34 -0.2090E+02  0.8345E+01 "alien"
   30   22   11    1.79   50.19 -0.1862E+02  0.7767E+01 "alien"
   35   33   15    1.36   40.39  0.2690E+02  0.1328E+02 "alien"
   32   20    1    2.05   39.94 -0.2582E+02  0.1509E+02 "alien"
   35    2    5    2.12   38.41 -0.2483E+02  0.1148E+02 "alien"
   21   16   13    2.09   38.17 -0.2468E+02  0.1008E+02 "alien"
   37    1   11    1.80   36.62 -0.1359E+02  0.7982E+01 "alien"
   35   10   21    1.39   35.89  0.1906E+02  0.1227E+02 "alien"
   24   15    9    2.31   34.34 -0.2744E+02  0.1188E+02 "alien"
   27   20   12    1.87   33.76 -0.1253E+02  0.7549E+01 "alien"
   34   33   16    1.35   33.57  0.2236E+02  0.1021E+02 "alien"
   39    1   10    1.76   31.90 -0.2383E+02  0.1373E+02 "alien"
   45    7   16    1.39   30.43  0.1616E+02  0.1346E+02 "alien"
   40    1    5    1.87   29.61 -0.1098E+02  0.1125E+02 "alien"
   37   10   21    1.36   29.44  0.1961E+02  0.1346E+02 "alien"
   28    1   10    2.25   28.84  0.2491E+02  0.1310E+02 "alien"
   29   19    8    2.03   28.21 -0.1824E+02  0.7978E+01 "alien"
   41    5   15    1.52   28.17 -0.2328E+02  0.1376E+02 "alien"
   52    5   10    1.38   26.38  0.1401E+02  0.9980E+01 "alien"
   50    6   13    1.36   24.34  0.1621E+02  0.9929E+01 "alien"
   30   21   10    1.85   23.15 -0.8589E+01  0.1319E+02 "alien"
   28   21   22    1.38   22.57  0.1199E+02  0.1266E+02 "alien"
   24   22   19    1.55   22.22 -0.1836E+02  0.1042E+02 "alien"
   31   29   19    1.36   22.16  0.1476E+02  0.1256E+02 "alien"
   31    7   19    1.56   21.54 -0.1780E+02  0.1033E+02 "alien"
   29   17    2    2.28   20.14  0.1740E+02  0.8348E+01 "alien"
   52    6   12    1.34   20.07  0.1337E+02  0.9571E+01 "alien"
   29   21   11    1.84   19.89 -0.7380E+01  0.7896E+01 "alien"
   38    9   20    1.38   19.55  0.1038E+02  0.1205E+02 "alien"
   39   34   11    1.38   19.11  0.1015E+02  0.1481E+02 "alien"
   38   33   12    1.39   18.97  0.1007E+02  0.8776E+01 "alien"
   31   22    9    1.84   18.68 -0.6931E+01  0.1336E+02 "alien"
   24   19   14    1.86   18.43 -0.6838E+01  0.1176E+02 "alien"
   27   20   22    1.40   18.42  0.9781E+01  0.1171E+02 "alien"
   52    4    9    1.40   17.93  0.9522E+01  0.1068E+02 "alien"
   23   14    8    2.47   17.88 -0.1429E+02  0.7838E+01 "alien"
   53    4    7    1.41   17.28  0.9178E+01  0.1005E+02 "alien"
   25   23   19    1.52   17.04 -0.1408E+02  0.1090E+02 "alien"
   29   22   22    1.35   16.88  0.1124E+02  0.1333E+02 "alien"
   34    1    8    2.06   16.65 -0.1076E+02  0.8958E+01 "alien"
   33   29   14    1.48   16.55 -0.2476E+02  0.1271E+02 "alien"
   22   16   12    2.13   16.54  0.1429E+02  0.1170E+02 "alien"
   53    7   13    1.30   16.51  0.1738E+03  0.1759E+02 "alien"
   53    3    4    1.44   16.42 -0.1941E+02  0.9192E+01 "alien"
   50    1    2    1.54   15.96 -0.1319E+02  0.1087E+02 "alien"
   51    5   11    1.39   15.71  0.8342E+01  0.9609E+01 "alien"
   51    6   12    1.36   15.66  0.1043E+02  0.1575E+02 "alien"
   36    9   20    1.42   15.64 -0.1849E+02  0.1187E+02 "alien"
   20   15   13    2.14   15.56  0.1344E+02  0.1011E+02 "alien"
   30    3   15    1.82   15.45 -0.5731E+01  0.1125E+02 "alien"
   32   20    3    2.02   15.42 -0.9971E+01  0.8152E+01 "alien"
   38    2   12    1.72   15.36 -0.1148E+02  0.1252E+02 "alien"
   37    2   13    1.71   15.35 -0.1147E+02  0.1171E+02 "alien"
   48    6   14    1.38   15.00  0.7968E+01  0.1234E+02 "alien"
   30    1    9    2.21   14.83  0.1281E+02  0.1193E+02 "alien"
   43    1    2    1.79   14.79 -0.5487E+01  0.8001E+01 "alien"
   47    1    6    1.59   14.48  0.3349E+02  0.1323E+02 "alien"
   50    1    0    1.55   14.45 -0.3249E+02  0.1251E+02 "alien"
   28   23   14    1.69   14.37  0.1369E+02  0.7584E+01 "alien"
   36   31   13    1.43   14.37 -0.1698E+02  0.1286E+02 "alien"
   33   11   22    1.38   14.25  0.7568E+01  0.1194E+02 "alien"
   40    2   12    1.66   14.19  0.1353E+02  0.8738E+01 "alien"
   22   21   17    1.70   14.19  0.1353E+02  0.1045E+02 "alien"
   33   25   11    1.65   14.13  0.1347E+02  0.1381E+02 "alien"
   30   18    2    2.20   14.03  0.1212E+02  0.9147E+01 "alien"
   32   24   11    1.69   14.01  0.1335E+02  0.1164E+02 "alien"
   20   18   16    1.85   13.93 -0.5168E+01  0.9235E+01 "alien"
   18   18   17    1.82   13.86 -0.1343E+02  0.9063E+01 "alien"
   28    1   13    2.02   13.79 -0.8917E+01  0.7647E+01 "alien"
   28   26   19    1.44   13.73 -0.1623E+02  0.1186E+02 "alien"
   25   19   13    1.90   13.73 -0.2058E+02  0.1032E+02 "alien"
   43    5   15    1.47   13.56 -0.2027E+02  0.1330E+02 "alien"
   29   20   23    1.34   13.38  0.8913E+01  0.1291E+02 "alien"
   37    1    7    1.95   13.36 -0.2002E+02  0.7966E+01 "alien"
   29    8   20    1.55   13.34 -0.1103E+02  0.1072E+02 "alien"
   27   16    6    2.30   13.31 -0.1064E+02  0.9128E+01 "alien"
   20   13   10    2.48   13.29 -0.1062E+02  0.9962E+01 "alien"
   22    1   13    2.27   13.04  0.1126E+02  0.9030E+01 "alien"
   30    2    7    2.33   12.81 -0.1024E+02  0.7687E+01 "alien"
   30   18    4    2.15   12.79  0.1105E+02  0.1014E+02 "alien"
   32   27   14    1.54   12.47 -0.1030E+02  0.1309E+02 "alien"
   47    2    8    1.56   12.47 -0.1030E+02  0.8182E+01 "alien"
   51    2    5    1.49   12.44 -0.1861E+02  0.1013E+02 "alien"
   35   11   22    1.35   12.44  0.8285E+01  0.1268E+02 "alien"
   45    3   11    1.54   12.11 -0.1001E+02  0.8417E+01 "alien"
   54    3    0    1.43   12.08 -0.3295E+02  0.1350E+02 "alien"
   46    5   14    1.43   11.92 -0.1408E+02  0.1497E+02 "alien"
   37    3   14    1.66   11.73  0.1118E+02  0.7660E+01 "alien"
   29    9   21    1.49   11.51 -0.1722E+02  0.1071E+02 "alien"
   26   23   20    1.46   11.49 -0.1719E+02  0.1200E+02 "alien"
   40   10   20    1.35   11.37  0.7575E+01  0.1444E+02 "alien"
   28   16    1    2.39   11.35 -0.9070E+01  0.8070E+01 "alien"
   43    0    7    1.71   11.23 -0.1830E+02  0.1421E+02 "alien"
   44    1    8    1.65   11.09  0.1057E+02  0.7956E+01 "alien"
   18   12   11    2.51   11.08 -0.8858E+01  0.8843E+01 "alien"
   25   15    7    2.39   11.03 -0.8814E+01  0.7265E+01 "alien"
   25   21   15    1.74   11.03 -0.8240E+01  0.1061E+02 "alien"
   40   32    7    1.46   10.96 -0.1296E+02  0.1432E+02 "alien"
   48    7   15    1.36   10.88  0.7247E+01  0.1286E+02 "alien"
   31    8   20    1.51   10.87 -0.8986E+01  0.8174E+01 "alien"
   17    2   14    2.36   10.86 -0.8680E+01  0.7750E+01 "alien"
   53    4    8    1.39   10.63  0.5648E+01  0.9604E+01 "alien"
   47    5   13    1.43   10.63 -0.1256E+02  0.8937E+01 "alien"
   35   25    8    1.68   10.62  0.1012E+02  0.1407E+02 "alien"
   30    1   13    1.95   10.57 -0.1584E+02  0.1256E+02 "alien"
   40    1   10    1.73   10.34 -0.7728E+01  0.7425E+01 "alien"
   33    4   16    1.68   10.30  0.9816E+01  0.7455E+01 "alien"
   42    1    2    1.83   10.24 -0.3798E+01  0.8783E+01 "alien"
   35   24    5    1.77   10.16 -0.7593E+01  0.7804E+01 "alien"
   42    1    9    1.69   10.09  0.9616E+01  0.7916E+01 "alien"
   31    3    0    2.48   10.07 -0.1630E+02  0.1383E+02 "alien"
   22   15   11    2.24   10.03  0.8666E+01  0.7346E+01 "alien"


 List of    113 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   16.51  0.1738E+03  0.1759E+02
   29   20   23    1.34   13.38  0.8913E+01  0.1291E+02
   52    6   12    1.34   20.07  0.1337E+02  0.9571E+01
   40   10   20    1.35   11.37  0.7575E+01  0.1444E+02
   35   11   22    1.35   12.44  0.8285E+01  0.1268E+02
   34   33   16    1.35   33.57  0.2236E+02  0.1021E+02
   29   22   22    1.35   16.88  0.1124E+02  0.1333E+02
   48    7   15    1.36   10.88  0.7247E+01  0.1286E+02
   31   29   19    1.36   22.16  0.1476E+02  0.1256E+02
   37   10   21    1.36   29.44  0.1961E+02  0.1346E+02
   35   33   15    1.36   40.39  0.2690E+02  0.1328E+02
   51    6   12    1.36   15.66  0.1043E+02  0.1575E+02
   50    6   13    1.36   24.34  0.1621E+02  0.9929E+01
   39   34   11    1.38   19.11  0.1015E+02  0.1481E+02
   28   21   22    1.38   22.57  0.1199E+02  0.1266E+02
   33   11   22    1.38   14.25  0.7568E+01  0.1194E+02
   52    5   10    1.38   26.38  0.1401E+02  0.9980E+01
   48    6   14    1.38   15.00  0.7968E+01  0.1234E+02
   38    9   20    1.38   19.55  0.1038E+02  0.1205E+02
   51    5   11    1.39   15.71  0.8342E+01  0.9609E+01
   38   33   12    1.39   18.97  0.1007E+02  0.8776E+01
   45    7   16    1.39   30.43  0.1616E+02  0.1346E+02
   35   10   21    1.39   35.89  0.1906E+02  0.1227E+02
   53    4    8    1.39   10.63  0.5648E+01  0.9604E+01
   27   20   22    1.40   18.42  0.9781E+01  0.1171E+02
   52    4    9    1.40   17.93  0.9522E+01  0.1068E+02
   53    4    7    1.41   17.28  0.9178E+01  0.1005E+02
   36    9   20    1.42   15.64 -0.1849E+02  0.1187E+02
   36   31   13    1.43   14.37 -0.1698E+02  0.1286E+02
   54    3    0    1.43   12.08 -0.3295E+02  0.1350E+02
   46    5   14    1.43   11.92 -0.1408E+02  0.1497E+02
   47    5   13    1.43   10.63 -0.1256E+02  0.8937E+01
   28   26   19    1.44   13.73 -0.1623E+02  0.1186E+02
   53    3    4    1.44   16.42 -0.1941E+02  0.9192E+01
   40   32    7    1.46   10.96 -0.1296E+02  0.1432E+02
   26   23   20    1.46   11.49 -0.1719E+02  0.1200E+02
   43    5   15    1.47   13.56 -0.2027E+02  0.1330E+02
   33   29   14    1.48   16.55 -0.2476E+02  0.1271E+02
   51    2    5    1.49   12.44 -0.1861E+02  0.1013E+02
   29    9   21    1.49   11.51 -0.1722E+02  0.1071E+02
   31    8   20    1.51   10.87 -0.8986E+01  0.8174E+01
   41    5   15    1.52   28.17 -0.2328E+02  0.1376E+02
   25   23   19    1.52   17.04 -0.1408E+02  0.1090E+02
   32   27   14    1.54   12.47 -0.1030E+02  0.1309E+02
   45    3   11    1.54   12.11 -0.1001E+02  0.8417E+01
   50    1    2    1.54   15.96 -0.1319E+02  0.1087E+02
   29    8   20    1.55   13.34 -0.1103E+02  0.1072E+02
   24   22   19    1.55   22.22 -0.1836E+02  0.1042E+02
   50    1    0    1.55   14.45 -0.3249E+02  0.1251E+02
   47    2    8    1.56   12.47 -0.1030E+02  0.8182E+01
   31    7   19    1.56   21.54 -0.1780E+02  0.1033E+02
   47    1    6    1.59   14.48  0.3349E+02  0.1323E+02
   33   25   11    1.65   14.13  0.1347E+02  0.1381E+02
   44    1    8    1.65   11.09  0.1057E+02  0.7956E+01
   40    2   12    1.66   14.19  0.1353E+02  0.8738E+01
   37    3   14    1.66   11.73  0.1118E+02  0.7660E+01
   33    4   16    1.68   10.30  0.9816E+01  0.7455E+01
   35   25    8    1.68   10.62  0.1012E+02  0.1407E+02
   28   23   14    1.69   14.37  0.1369E+02  0.7584E+01
   42    1    9    1.69   10.09  0.9616E+01  0.7916E+01
   32   24   11    1.69   14.01  0.1335E+02  0.1164E+02
   22   21   17    1.70   14.19  0.1353E+02  0.1045E+02
   37    2   13    1.71   15.35 -0.1147E+02  0.1171E+02
   43    0    7    1.71   11.23 -0.1830E+02  0.1421E+02
   38    2   12    1.72   15.36 -0.1148E+02  0.1252E+02
   40    1   10    1.73   10.34 -0.7728E+01  0.7425E+01
   25   21   15    1.74   11.03 -0.8240E+01  0.1061E+02
   39    1   10    1.76   31.90 -0.2383E+02  0.1373E+02
   35   24    5    1.77   10.16 -0.7593E+01  0.7804E+01
   30   22   11    1.79   50.19 -0.1862E+02  0.7767E+01
   43    1    2    1.79   14.79 -0.5487E+01  0.8001E+01
   37    1   11    1.80   36.62 -0.1359E+02  0.7982E+01
   21   19   16    1.81   59.42 -0.2205E+02  0.1021E+02
   30    3   15    1.82   15.45 -0.5731E+01  0.1125E+02
   18   18   17    1.82   13.86 -0.1343E+02  0.9063E+01
   35   23    3    1.83   56.34 -0.2090E+02  0.8345E+01
   42    1    2    1.83   10.24 -0.3798E+01  0.8783E+01
   31   22    9    1.84   18.68 -0.6931E+01  0.1336E+02
   29   21   11    1.84   19.89 -0.7380E+01  0.7896E+01
   20   18   16    1.85   13.93 -0.5168E+01  0.9235E+01
   30   21   10    1.85   23.15 -0.8589E+01  0.1319E+02
   24   19   14    1.86   18.43 -0.6838E+01  0.1176E+02
   27   20   12    1.87   33.76 -0.1253E+02  0.7549E+01
   40    1    5    1.87   29.61 -0.1098E+02  0.1125E+02
   25   19   13    1.90   13.73 -0.2058E+02  0.1032E+02
   30    1   13    1.95   10.57 -0.1584E+02  0.1256E+02
   37    1    7    1.95   13.36 -0.2002E+02  0.7966E+01
   32   20    3    2.02   15.42 -0.9971E+01  0.8152E+01
   28    1   13    2.02   13.79 -0.8917E+01  0.7647E+01
   29   19    8    2.03   28.21 -0.1824E+02  0.7978E+01
   32   20    1    2.05   39.94 -0.2582E+02  0.1509E+02
   34    1    8    2.06   16.65 -0.1076E+02  0.8958E+01
   21   16   13    2.09   38.17 -0.2468E+02  0.1008E+02
   35    2    5    2.12   38.41 -0.2483E+02  0.1148E+02
   22   16   12    2.13   16.54  0.1429E+02  0.1170E+02
   20   15   13    2.14   15.56  0.1344E+02  0.1011E+02
   30   18    4    2.15   12.79  0.1105E+02  0.1014E+02
   30   18    2    2.20   14.03  0.1212E+02  0.9147E+01
   30    1    9    2.21   14.83  0.1281E+02  0.1193E+02
   22   15   11    2.24   10.03  0.8666E+01  0.7346E+01
   28    1   10    2.25   28.84  0.2491E+02  0.1310E+02
   22    1   13    2.27   13.04  0.1126E+02  0.9030E+01
   29   17    2    2.28   20.14  0.1740E+02  0.8348E+01
   27   16    6    2.30   13.31 -0.1064E+02  0.9128E+01
   24   15    9    2.31   34.34 -0.2744E+02  0.1188E+02
   30    2    7    2.33   12.81 -0.1024E+02  0.7687E+01
   17    2   14    2.36   10.86 -0.8680E+01  0.7750E+01
   25   15    7    2.39   11.03 -0.8814E+01  0.7265E+01
   28   16    1    2.39   11.35 -0.9070E+01  0.8070E+01
   23   14    8    2.47   17.88 -0.1429E+02  0.7838E+01
   31    3    0    2.48   10.07 -0.1630E+02  0.1383E+02
   20   13   10    2.48   13.29 -0.1062E+02  0.9962E+01
   18   12   11    2.51   11.08 -0.8858E+01  0.8843E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        28
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1051
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)     987
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       64


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
