
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.326    77.326    38.682  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=       1       1
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       1
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000852  0.000719
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
 BEAM_DIVERGENCE_E.S.D.=   0.038
 REFLECTING_RANGE_E.S.D.=   0.064
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

    1571 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       1 ...       1
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1571 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               49
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
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    27        mC        500.0     159.4   38.7  109.4  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.7  109.4  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.4  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        999.0     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        999.0      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        999.0      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        999.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.5   86.5  116.0 116.6  63.4 101.5    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.5   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.6  63.4 101.5    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.4   86.5   86.5  78.5  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.5   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.5    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.7   86.5  247.6 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.5   86.5   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0  1  0  0  1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.2  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.5   86.5  177.2 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.5   86.5   86.5  78.5 101.5  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.2  109.4   38.7  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.7   77.3  177.2 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.5  177.2   77.3 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    4.8      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.703   -23.387     0.035
 UNIT_CELL_B-AXIS=    23.387    73.703     0.012
 UNIT_CELL_C-AXIS=    -0.019    -0.001    38.682
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       350 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.14
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.06
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.251    77.251    38.660  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.3E-01 1.3E-01 6.4E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012945  0.012945  0.025867  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.628   -23.379     0.025
 COORDINATES OF UNIT CELL B-AXIS    23.379    73.628     0.007
 COORDINATES OF UNIT CELL C-AXIS    -0.013     0.000    38.660
 CRYSTAL MOSAICITY (DEGREES)     0.064
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006546  0.002063
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.86   1232.90
 DETECTOR ORIGIN (PIXELS) AT                     1260.86   1242.43
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.21
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
 -0.05  0.05         4        955.976       0.376
  0.05  0.15       143        457.935       0.784
  0.15  0.25       184        305.083       1.177
  0.25  0.35       193        316.401       1.135
  0.35  0.45       154        189.615       1.894
  0.45  0.55       182        388.969       0.923
  0.55  0.65       196        368.323       0.975
  0.65  0.75       185        241.270       1.489
  0.75  0.85       170        336.635       1.067
  0.85  0.95       160        647.777       0.554
 
 -0.05  0.95      1571        359.135       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               240.6
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.037
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     0.0 XMAX=     1.0 NXBIN=    1
 YMIN= 0.00168 YMAX= 0.59763 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        524


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               242.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.034
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    11.0 XMAX=  2454.9 NXBIN=    3
 YMIN=     6.0 YMAX=  2519.7 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        524


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               242.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.030
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     0.0 XMAX=     1.0 NXBIN=    1
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
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        524


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
 1.192E+00  3.857E-03   14.75


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

   24.429   3.394    13.39   1.14      2.59      3.89      68      97      20
    3.394   2.412    11.59   1.54      5.03      4.13     134     176      12
    2.412   1.972     9.10   1.02      5.06      4.69     186     228       6
    1.972   1.709     0.02   0.70   3676.90   4278.03     202     256       0
    1.709   1.530    -0.03   0.52    -99.90    -99.90     228     300       0
    1.530   1.397     0.02   0.59    700.69    868.99     254     325       0
    1.397   1.294     0.05   0.60    672.84    819.56     100     148       2
   --------------------------------------------------------------------------
   24.429   1.294     3.54   0.80      5.11      5.58    1172    1530      40
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       1       1)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.679   0.1393E+02  0.1220E+02         1.14           1 
   47    0    0     1.644  -0.5921E+01  0.1187E+02        -0.50           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -1.5% OF MEAN INTENSITY


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

     3.84          73      48      1285        3.7%       2.5%      3.9%       50   16.67      3.5%    99.9*     0    0.000       0
     2.73         109      69      2098        3.3%       4.5%      4.0%       80   14.87      6.4%    98.8*     0    0.000       0
     2.23         175     105      2665        3.9%       4.9%      4.2%      140   14.34      6.9%    98.5*     0    0.000       0
     1.94         171     102      3101        3.3%       5.5%      5.3%      138    8.19      7.8%    99.2*     0    0.000       0
     1.73         197     119      3479        3.4%     -99.9%    -99.9%      156    0.00    -99.9%     9.8      0    0.000       0
     1.58         231     139      3850        3.6%     -99.9%    -99.9%      184    0.00    -99.9%    12.7      0    0.000       0
     1.47         260     163      4143        3.9%     689.2%    839.1%      194    0.06    974.7%     0.2      0    0.000       0
     1.37         225     136      4428        3.1%     798.3%   1015.2%      178   -0.00   1129.0%     2.1      0    0.000       0
     1.29          89      63      4739        1.3%    1233.4%   1492.3%       52    0.02   1744.4%   -30.2      0    0.000       0
    total        1530     944     29788        3.2%       5.1%      5.6%     1172    4.42      7.2%    99.7*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1571
 NUMBER OF REJECTED MISFITS                                40
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1530
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    944



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

   WILSON LINE (using all data) : A=   8.448 B=  38.403 CORRELATION=  0.80
      #      RES      SS        <I>       log(<I>)       BO
      19     7.383   0.005  3.8141E+03       8.246      21.9
      24     4.651   0.012  3.2057E+03       8.073      16.2
      32     3.507   0.020  2.0929E+03       7.646      19.7
      32     2.996   0.028  8.2617E+02       6.717      31.1
      42     2.641   0.036  9.5153E+02       6.858      22.2
      57     2.393   0.044  7.4104E+02       6.608      21.1
      45     2.204   0.051  5.0533E+02       6.225      21.6
      48     2.054   0.059  2.8533E+02       5.654      23.6
      53     1.930   0.067  1.6652E-01      -1.793      76.3
      65     1.825   0.075  2.2070E-01      -1.511      66.4
      56     1.736   0.083 -1.2427E-01    -999.999   -1000.0
      61     1.662   0.090  4.5278E-01      -0.792      51.1
      65     1.595   0.098 -9.5345E-01    -999.999   -1000.0
      85     1.534   0.106  9.7968E-01      -0.021      39.8
      83     1.479   0.114  7.5354E-01      -0.283      38.2
      61     1.430   0.122 -4.1299E-01    -999.999   -1000.0
      55     1.390   0.129 -2.1379E-02    -999.999   -1000.0
      41     1.350   0.137 -8.5256E-01    -999.999   -1000.0
      20     1.312   0.145  3.0704E+00       1.122      25.2


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       6     7.383        1.641        1.179        0.674
       3     4.651        1.141        0.653        0.279
       4     3.507        5.849        9.360       10.972
       4     2.996        1.261        0.607        0.222
       3     2.641        0.287        0.076        0.015
       5     2.393        0.656        0.277        0.087
       3     2.204        1.244        0.489        0.139
       8     2.054        0.230        0.072        0.017
       7     1.930     1014.093    13325.724   218033.446
      12     1.825      725.273    15147.318   296391.155
      11     1.736       17.684       10.600       47.085
       4     1.595       27.204       72.755      146.628
       6     1.479       51.939      -95.785      592.775
      96   overall      184.551     2841.729    53334.234


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      13     7.383        0.821        0.996        0.996
      21     4.651        1.078        1.211        1.187
      28     3.507        0.506        0.295        0.155
      28     2.996        0.735        0.494        0.273
      39     2.641        1.222        1.379        1.375
      52     2.393        0.954        0.788        0.547
      42     2.204        0.821        0.624        0.407
      40     2.054        1.309        1.424        1.402
      50     1.662      121.730      601.523     7561.069
      82     1.534       30.602       17.492      588.088
      77     1.479       58.403      103.661     1867.179
      54     1.390       21.717       17.849      222.032
      20     1.312       10.316       33.854      138.833
     848   overall      107.341      234.664    49758.469

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       6     7.383  0.00  0.00  0.00  0.00  0.32  0.30  0.56  0.53  0.51  0.49
       3     4.651  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       4     3.507  1.01  0.72  0.60  0.53  0.48  0.45  0.42  0.40  0.38  0.37
       4     2.996  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.641  1.34  0.96  0.80  0.70  0.64  0.59  1.12  1.06  1.01  0.98
       5     2.393  0.00  0.00  0.00  0.85  0.77  1.07  1.00  0.95  0.91  0.88
       3     2.204  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.054  1.51  1.45  1.20  1.32  1.20  1.11  1.05  1.19  1.14  1.10
       7     1.930  2.30  1.65  1.37  1.21  1.10  1.02  0.96  0.91  0.87  0.84
      12     1.825  1.68  1.21  1.00  0.88  0.80  0.74  0.70  0.66  0.63  0.61
      11     1.736  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.72  0.69  0.67
       4     1.595  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       6     1.479  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      96   overall  1.09  0.81  0.68  0.66  0.62  0.61  0.61  0.60  0.57  0.56


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      13     7.383  2.42  2.55  2.67  2.33  1.95  1.70  1.53  1.40  1.30  1.22
      21     4.651  1.50  1.58  1.10  1.16  1.21  1.16  1.14  1.04  1.12  1.05
      28     3.507  1.50  0.99  0.96  0.87  1.09  1.11  1.14  1.04  1.20  1.13
      28     2.996  1.88  1.38  1.10  1.19  1.09  1.03  0.99  1.10  1.08  1.02
      39     2.641  1.62  1.41  1.29  1.24  1.11  1.02  1.07  1.07  1.08  1.05
      52     2.393  1.41  1.17  1.04  1.17  1.03  1.07  1.11  1.05  1.07  1.03
      42     2.204  1.50  1.58  1.19  0.94  1.03  1.00  0.99  0.95  1.04  1.02
      40     2.054  2.63  1.52  1.35  1.14  1.02  0.94  0.89  0.95  0.93  0.95
      50     1.662  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      82     1.534  5.00  2.76  1.98  1.55  1.30  1.14  1.02  0.93  0.86  0.81
      77     1.479  4.91  2.58  1.80  1.54  1.29  1.12  1.01  0.92  0.85  0.80
      54     1.390  4.67  2.45  1.71  1.35  1.13  0.99  0.88  0.81  0.75  0.70
      20     1.312  4.73  2.48  1.74  1.36  1.14  1.00  0.89  0.91  0.84  0.79
     848   overall  2.19  1.32  1.00  0.86  0.76  0.70  0.66  0.62  0.61  0.58


 List of    171 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   48    6   14    1.38  733.32 -0.1568E+02  0.1204E+02 "alien"
   37   33   13    1.38  479.12 -0.1024E+02  0.9540E+01 "alien"
   54    4    6    1.39  428.54 -0.9162E+01  0.9198E+01 "alien"
   28   12   23    1.40  410.67 -0.8780E+01  0.9110E+01 "alien"
   38    9   20    1.38  406.78 -0.8697E+01  0.1067E+02 "alien"
   51    4   10    1.41  406.76 -0.8696E+01  0.1271E+02 "alien"
   39   34   11    1.37  318.10 -0.6801E+01  0.9781E+01 "alien"
   50    5   11    1.41  285.15 -0.6096E+01  0.1308E+02 "alien"
   31   30   17    1.41  268.01 -0.5730E+01  0.9409E+01 "alien"
   34   31   15    1.41  267.08 -0.5710E+01  0.8947E+01 "alien"
   38   33   12    1.39  260.13 -0.5561E+01  0.9406E+01 "alien"
   34   10   21    1.41  222.92 -0.4766E+01  0.1099E+02 "alien"
   36   32   14    1.39  220.44 -0.4713E+01  0.1562E+02 "alien"
   33   11   22    1.38  197.75 -0.4228E+01  0.7227E+01 "alien"
   32   23   10    1.75  185.20 -0.2301E+02  0.1341E+02 "alien"
   29   26   20    1.38  178.76 -0.3822E+01  0.1347E+02 "alien"
   46    6   15    1.40  167.37 -0.3578E+01  0.8433E+01 "alien"
   41    8   18    1.40  165.94 -0.3548E+01  0.1236E+02 "alien"
   54    4    7    1.38  165.60 -0.3540E+01  0.9557E+01 "alien"
   34    3   15    1.70  155.19 -0.1929E+02  0.1001E+02 "alien"
   28   19   10    1.97  145.74  0.2427E+02  0.1547E+02 "alien"
   28   21   22    1.37  139.50 -0.2982E+01  0.7480E+01 "alien"
   36   25    6    1.70  114.85 -0.1427E+02  0.1410E+02 "alien"
   27   20   22    1.40  103.40 -0.2211E+01  0.7286E+01 "alien"
   46    7   16    1.37   90.72 -0.1940E+01  0.8575E+01 "alien"
   38    1    7    1.91   87.31  0.1454E+02  0.8218E+01 "alien"
   30   28   19    1.38   82.32 -0.1760E+01  0.8348E+01 "alien"
   30   29   18    1.40   80.30 -0.1717E+01  0.1289E+02 "alien"
   37   26    2    1.70   78.67 -0.9776E+01  0.1556E+02 "alien"
   27   22   14    1.73   73.96 -0.9191E+01  0.7133E+01 "alien"
   29   20   10    1.91   73.16  0.1218E+02  0.1214E+02 "alien"
   24   19   14    1.86   72.57  0.1602E+02  0.1046E+02 "alien"
   35    2   13    1.77   72.05 -0.8954E+01  0.6873E+01 "alien"
   34   23    6    1.81   62.46  0.1378E+02  0.8933E+01 "alien"
   28   21   12    1.82   61.94  0.1367E+02  0.7931E+01 "alien"
   44    1    1    1.75   59.77 -0.7427E+01  0.8271E+01 "alien"
   25   21   15    1.74   58.81 -0.7308E+01  0.9556E+01 "alien"
   30    2   14    1.88   58.47  0.9736E+01  0.6516E+01 "alien"
   28    2   14    1.95   55.85  0.9300E+01  0.9264E+01 "alien"
   36   10   21    1.37   55.83 -0.1193E+01  0.7402E+01 "alien"
   26   21   14    1.77   55.47 -0.6893E+01  0.7165E+01 "alien"
   29   27   19    1.41   54.60 -0.1167E+01  0.8190E+01 "alien"
   31   30   18    1.38   52.92 -0.1131E+01  0.8889E+01 "alien"
   29    1   13    1.98   49.51  0.8244E+01  0.7230E+01 "alien"
   41    1    3    1.86   48.26  0.1065E+02  0.8909E+01 "alien"
   27   21   13    1.80   45.39  0.1002E+02  0.7555E+01 "alien"
   20   19   18    1.70   44.82 -0.5570E+01  0.6060E+01 "alien"
   37   25    1    1.73   44.03 -0.5472E+01  0.8832E+01 "alien"
   25   13   24    1.39   41.08 -0.8782E+00  0.5413E+01 "alien"
   37   26    5    1.67   40.29  0.1824E+02  0.9181E+01 "alien"
   31    5   17    1.67   39.00  0.1766E+02  0.1078E+02 "alien"
   38    2    4    1.99   38.41  0.6397E+01  0.9521E+01 "alien"
   53    3    3    1.45   36.84 -0.1522E+02  0.9515E+01 "alien"
   31   31   17    1.39   34.70 -0.1164E+02  0.1515E+02 "alien"
   51    3    6    1.47   34.41  0.2593E+02  0.1700E+02 "alien"
   40    1    4    1.89   34.38  0.5725E+01  0.1434E+02 "alien"
   43    2   10    1.63   32.94  0.1491E+02  0.9075E+01 "alien"
   21   18   15    1.89   32.01  0.5330E+01  0.6684E+01 "alien"
   36   24    2    1.78   30.76 -0.3822E+01  0.8709E+01 "alien"
   40   10   20    1.35   30.61 -0.2610E+02  0.1131E+02 "alien"
   35   25    7    1.71   28.72 -0.3570E+01  0.8827E+01 "alien"
   39    1   11    1.73   28.63 -0.3558E+01  0.1114E+02 "alien"
   34   22    4    1.87   27.90  0.6158E+01  0.9193E+01 "alien"
   38    8   19    1.42   27.41 -0.1132E+02  0.7420E+01 "alien"
   14   10   18    1.94   26.69  0.4445E+01  0.5343E+01 "alien"
   20   17   19    1.67   26.69  0.1208E+02  0.8642E+01 "alien"
   17   16   16    1.95   25.10  0.4180E+01  0.5715E+01 "alien"
   15    8   20    1.78   24.30 -0.3020E+01  0.4429E+01 "alien"
   33   23    7    1.81   24.21  0.5343E+01  0.1238E+02 "alien"
   40   32    8    1.44   23.96 -0.9896E+01  0.9552E+01 "alien"
   37    0   10    1.84   23.91  0.3056E+02  0.1445E+02 "alien"
   39    1    6    1.89   23.39  0.3894E+01  0.8283E+01 "alien"
   43    6   16    1.43   23.15 -0.9559E+01  0.7974E+01 "alien"
   43    1    1    1.79   23.12  0.5103E+01  0.8596E+01 "alien"
   13    6   19    1.90   22.84  0.3803E+01  0.3854E+01 "alien"
   49    5   12    1.41   22.77 -0.9403E+01  0.8851E+01 "alien"
   15   10   19    1.84   22.68  0.5006E+01  0.5738E+01 "alien"
   43    1    9    1.66   21.95  0.9940E+01  0.7807E+01 "alien"
   34   22    0    1.91   21.94  0.1789E+02  0.1297E+02 "alien"
   29   20   23    1.33   21.15 -0.1803E+02  0.1054E+02 "alien"
   42    1    1    1.84   21.10  0.4657E+01  0.8513E+01 "alien"
   46    4   12    1.48   21.09  0.1590E+02  0.8818E+01 "alien"
   31   27   15    1.52   21.03  0.2060E+02  0.1349E+02 "alien"
   29   21   11    1.84   20.85  0.4601E+01  0.7558E+01 "alien"
   35   26    8    1.66   20.61  0.9333E+01  0.1339E+02 "alien"
   36    1   11    1.83   20.38  0.4499E+01  0.1131E+02 "alien"
   30    3   15    1.82   20.22  0.4462E+01  0.8896E+01 "alien"
   21   19   19    1.63   20.20  0.9145E+01  0.8771E+01 "alien"
   29   28   17    1.47   20.02  0.1508E+02  0.8625E+01 "alien"
   37    1    7    1.95   19.92  0.3317E+01  0.8089E+01 "alien"
   27   24   16    1.60   19.72 -0.1880E+02  0.1179E+02 "alien"
   43    1    8    1.68   19.63  0.8886E+01  0.8495E+01 "alien"
   34   25    9    1.68   18.78  0.8503E+01  0.8603E+01 "alien"
   21   13   24    1.43   18.60 -0.7684E+01  0.6070E+01 "alien"
   44    0    5    1.71   18.58 -0.8584E+01  0.8162E+01 "alien"
   19   14   20    1.66   18.52  0.8385E+01  0.5180E+01 "alien"
   24   14   23    1.44   18.46 -0.7625E+01  0.5772E+01 "alien"
   53    2    1    1.46   18.44  0.1389E+02  0.1064E+02 "alien"
   42    3   13    1.56   18.30  0.1793E+02  0.1392E+02 "alien"
   34   26   10    1.64   18.24  0.8257E+01  0.8459E+01 "alien"
   31   29   16    1.45   18.10 -0.7475E+01  0.8714E+01 "alien"
   43    0    7    1.71   17.87 -0.7815E+01  0.9062E+01 "alien"
   42    0    6    1.77   17.78 -0.7716E+01  0.9435E+01 "alien"
   38    6   17    1.51   17.71  0.1334E+02  0.7744E+01 "alien"
   35   31   14    1.42   17.55 -0.7249E+01  0.1267E+02 "alien"
   34    0   10    1.96   17.38  0.9795E+01  0.7908E+01 "alien"
   17    9   21    1.67   17.11  0.7749E+01  0.5257E+01 "alien"
   37   30   11    1.47   16.99  0.1280E+02  0.1325E+02 "alien"
   39    2    2    1.97   16.81  0.2799E+01  0.8592E+01 "alien"
   54    3    3    1.42   16.39 -0.6770E+01  0.9355E+01 "alien"
   33    0   11    1.95   16.21  0.8367E+01  0.1122E+02 "alien"
   43    1    3    1.78   16.15 -0.2007E+01  0.8368E+01 "alien"
   41    1    5    1.83   15.98  0.3526E+01  0.8469E+01 "alien"
   53    7   13    1.30   15.92  0.4887E+02  0.1373E+02 "alien"
   51    3    8    1.44   15.77 -0.6512E+01  0.1000E+02 "alien"
   26   23   20    1.46   15.74  0.1186E+02  0.1124E+02 "alien"
   35   23    1    1.84   15.66  0.3457E+01  0.8976E+01 "alien"
   44    5   14    1.47   15.56  0.1173E+02  0.9466E+01 "alien"
   36   29   11    1.51   15.39  0.1508E+02  0.9219E+01 "alien"
   45    1    8    1.62   15.38 -0.1467E+02  0.8270E+01 "alien"
   28   28   18    1.44   15.16 -0.1795E+02  0.1182E+02 "alien"
   33   23    8    1.78   14.99  0.3308E+01  0.8238E+01 "alien"
   32   23    9    1.78   14.96  0.3301E+01  0.7950E+01 "alien"
   31    3   15    1.79   14.46  0.3190E+01  0.6509E+01 "alien"
   36   28    9    1.58   14.37 -0.1370E+02  0.1296E+02 "alien"
   44    1    0    1.76   14.34 -0.4812E+01  0.1114E+02 "alien"
   37   25    0    1.73   14.31 -0.4793E+01  0.1245E+02 "alien"
   30   23   12    1.73   14.26 -0.1772E+01  0.8098E+01 "alien"
   39   29    6    1.54   14.20  0.1391E+02  0.9684E+01 "alien"
   27   25   19    1.46   14.17  0.1068E+02  0.8561E+01 "alien"
   31   22   10    1.80   14.03  0.3096E+01  0.1228E+02 "alien"
   40   31    7    1.47   13.97  0.1053E+02  0.9969E+01 "alien"
   45    6   15    1.42   13.92 -0.5747E+01  0.8610E+01 "alien"
   30   28   16    1.48   13.75  0.1036E+02  0.8665E+01 "alien"
   33    3   15    1.73   13.71 -0.1704E+01  0.6575E+01 "alien"
   48    7   15    1.35   13.65 -0.1164E+02  0.1302E+02 "alien"
   38   29    8    1.53   13.59  0.1332E+02  0.9455E+01 "alien"
   44    0    6    1.69   13.52  0.1559E+02  0.8553E+01 "alien"
   39   28    4    1.59   13.48 -0.1285E+02  0.1075E+02 "alien"
   38    0    9    1.84   13.45  0.7522E+01  0.7837E+01 "alien"
   33   27   13    1.55   13.44  0.1317E+02  0.1254E+02 "alien"
   46    5   14    1.43   13.28 -0.5485E+01  0.1236E+02 "alien"
   18   18   17    1.82   13.01  0.7065E+01  0.6488E+01 "alien"
   29   27   16    1.52   13.01  0.1274E+02  0.1186E+02 "alien"
   38    1    6    1.94   12.99  0.2162E+01  0.8364E+01 "alien"
   41   31    1    1.50   12.89  0.9712E+01  0.9558E+01 "alien"
   35   10   21    1.39   12.86 -0.2750E+00  0.7604E+01 "alien"
   40    0    8    1.79   12.85  0.6896E+01  0.7954E+01 "alien"
   50    1    1    1.54   12.69  0.1243E+02  0.8942E+01 "alien"
   20   17   15    1.94   12.62  0.2101E+01  0.6538E+01 "alien"
   23   21   19    1.57   12.58 -0.1200E+02  0.6297E+01 "alien"
   42    7   17    1.42   12.50 -0.5161E+01  0.7864E+01 "alien"
   26   24   19    1.49   12.41  0.9351E+01  0.8690E+01 "alien"
   38   32   11    1.42   12.37 -0.5110E+01  0.9603E+01 "alien"
   36    0   10    1.88   11.91  0.5998E+01  0.7557E+01 "alien"
   19   18   18    1.74   11.74 -0.1459E+01  0.6361E+01 "alien"
   33   29   14    1.48   11.72  0.8830E+01  0.8851E+01 "alien"
   52    3    5    1.46   11.54  0.8696E+01  0.9160E+01 "alien"
   34    5   17    1.60   11.43 -0.1090E+02  0.6922E+01 "alien"
   47    1    5    1.61   11.36 -0.1083E+02  0.8203E+01 "alien"
   24   20   20    1.52   11.15  0.1093E+02  0.1055E+02 "alien"
   36   24    0    1.79   11.08  0.5288E+01  0.1251E+02 "alien"
   43    1    2    1.79   11.01  0.2431E+01  0.8192E+01 "alien"
   16   15   17    1.91   10.95  0.1823E+01  0.8207E+01 "alien"
   31   29   19    1.36   10.94 -0.9323E+01  0.8778E+01 "alien"
   37   27    6    1.63   10.81  0.4893E+01  0.1248E+02 "alien"
   26   21   21    1.44   10.50 -0.4335E+01  0.9730E+01 "alien"
   25   24   17    1.59   10.20 -0.9727E+01  0.7443E+01 "alien"
   17    8   20    1.75   10.13 -0.1259E+01  0.5028E+01 "alien"
   32   30   19    1.33   10.10 -0.8609E+01  0.1276E+02 "alien"
   33   33   17    1.34   10.08 -0.1741E+02  0.1283E+02 "alien"


 List of    171 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   15.92  0.4887E+02  0.1373E+02
   32   30   19    1.33   10.10 -0.8609E+01  0.1276E+02
   29   20   23    1.33   21.15 -0.1803E+02  0.1054E+02
   33   33   17    1.34   10.08 -0.1741E+02  0.1283E+02
   40   10   20    1.35   30.61 -0.2610E+02  0.1131E+02
   48    7   15    1.35   13.65 -0.1164E+02  0.1302E+02
   31   29   19    1.36   10.94 -0.9323E+01  0.8778E+01
   46    7   16    1.37   90.72 -0.1940E+01  0.8575E+01
   39   34   11    1.37  318.10 -0.6801E+01  0.9781E+01
   28   21   22    1.37  139.50 -0.2982E+01  0.7480E+01
   36   10   21    1.37   55.83 -0.1193E+01  0.7402E+01
   31   30   18    1.38   52.92 -0.1131E+01  0.8889E+01
   33   11   22    1.38  197.75 -0.4228E+01  0.7227E+01
   37   33   13    1.38  479.12 -0.1024E+02  0.9540E+01
   54    4    7    1.38  165.60 -0.3540E+01  0.9557E+01
   30   28   19    1.38   82.32 -0.1760E+01  0.8348E+01
   48    6   14    1.38  733.32 -0.1568E+02  0.1204E+02
   38    9   20    1.38  406.78 -0.8697E+01  0.1067E+02
   29   26   20    1.38  178.76 -0.3822E+01  0.1347E+02
   38   33   12    1.39  260.13 -0.5561E+01  0.9406E+01
   36   32   14    1.39  220.44 -0.4713E+01  0.1562E+02
   25   13   24    1.39   41.08 -0.8782E+00  0.5413E+01
   35   10   21    1.39   12.86 -0.2750E+00  0.7604E+01
   31   31   17    1.39   34.70 -0.1164E+02  0.1515E+02
   54    4    6    1.39  428.54 -0.9162E+01  0.9198E+01
   27   20   22    1.40  103.40 -0.2211E+01  0.7286E+01
   46    6   15    1.40  167.37 -0.3578E+01  0.8433E+01
   28   12   23    1.40  410.67 -0.8780E+01  0.9110E+01
   41    8   18    1.40  165.94 -0.3548E+01  0.1236E+02
   30   29   18    1.40   80.30 -0.1717E+01  0.1289E+02
   34   10   21    1.41  222.92 -0.4766E+01  0.1099E+02
   51    4   10    1.41  406.76 -0.8696E+01  0.1271E+02
   34   31   15    1.41  267.08 -0.5710E+01  0.8947E+01
   31   30   17    1.41  268.01 -0.5730E+01  0.9409E+01
   29   27   19    1.41   54.60 -0.1167E+01  0.8190E+01
   50    5   11    1.41  285.15 -0.6096E+01  0.1308E+02
   49    5   12    1.41   22.77 -0.9403E+01  0.8851E+01
   35   31   14    1.42   17.55 -0.7249E+01  0.1267E+02
   42    7   17    1.42   12.50 -0.5161E+01  0.7864E+01
   54    3    3    1.42   16.39 -0.6770E+01  0.9355E+01
   45    6   15    1.42   13.92 -0.5747E+01  0.8610E+01
   38   32   11    1.42   12.37 -0.5110E+01  0.9603E+01
   38    8   19    1.42   27.41 -0.1132E+02  0.7420E+01
   46    5   14    1.43   13.28 -0.5485E+01  0.1236E+02
   21   13   24    1.43   18.60 -0.7684E+01  0.6070E+01
   43    6   16    1.43   23.15 -0.9559E+01  0.7974E+01
   24   14   23    1.44   18.46 -0.7625E+01  0.5772E+01
   40   32    8    1.44   23.96 -0.9896E+01  0.9552E+01
   26   21   21    1.44   10.50 -0.4335E+01  0.9730E+01
   51    3    8    1.44   15.77 -0.6512E+01  0.1000E+02
   28   28   18    1.44   15.16 -0.1795E+02  0.1182E+02
   53    3    3    1.45   36.84 -0.1522E+02  0.9515E+01
   31   29   16    1.45   18.10 -0.7475E+01  0.8714E+01
   53    2    1    1.46   18.44  0.1389E+02  0.1064E+02
   52    3    5    1.46   11.54  0.8696E+01  0.9160E+01
   26   23   20    1.46   15.74  0.1186E+02  0.1124E+02
   27   25   19    1.46   14.17  0.1068E+02  0.8561E+01
   29   28   17    1.47   20.02  0.1508E+02  0.8625E+01
   40   31    7    1.47   13.97  0.1053E+02  0.9969E+01
   51    3    6    1.47   34.41  0.2593E+02  0.1700E+02
   37   30   11    1.47   16.99  0.1280E+02  0.1325E+02
   44    5   14    1.47   15.56  0.1173E+02  0.9466E+01
   33   29   14    1.48   11.72  0.8830E+01  0.8851E+01
   46    4   12    1.48   21.09  0.1590E+02  0.8818E+01
   30   28   16    1.48   13.75  0.1036E+02  0.8665E+01
   26   24   19    1.49   12.41  0.9351E+01  0.8690E+01
   41   31    1    1.50   12.89  0.9712E+01  0.9558E+01
   38    6   17    1.51   17.71  0.1334E+02  0.7744E+01
   36   29   11    1.51   15.39  0.1508E+02  0.9219E+01
   29   27   16    1.52   13.01  0.1274E+02  0.1186E+02
   31   27   15    1.52   21.03  0.2060E+02  0.1349E+02
   24   20   20    1.52   11.15  0.1093E+02  0.1055E+02
   38   29    8    1.53   13.59  0.1332E+02  0.9455E+01
   39   29    6    1.54   14.20  0.1391E+02  0.9684E+01
   50    1    1    1.54   12.69  0.1243E+02  0.8942E+01
   33   27   13    1.55   13.44  0.1317E+02  0.1254E+02
   42    3   13    1.56   18.30  0.1793E+02  0.1392E+02
   23   21   19    1.57   12.58 -0.1200E+02  0.6297E+01
   36   28    9    1.58   14.37 -0.1370E+02  0.1296E+02
   39   28    4    1.59   13.48 -0.1285E+02  0.1075E+02
   25   24   17    1.59   10.20 -0.9727E+01  0.7443E+01
   34    5   17    1.60   11.43 -0.1090E+02  0.6922E+01
   27   24   16    1.60   19.72 -0.1880E+02  0.1179E+02
   47    1    5    1.61   11.36 -0.1083E+02  0.8203E+01
   45    1    8    1.62   15.38 -0.1467E+02  0.8270E+01
   43    2   10    1.63   32.94  0.1491E+02  0.9075E+01
   21   19   19    1.63   20.20  0.9145E+01  0.8771E+01
   37   27    6    1.63   10.81  0.4893E+01  0.1248E+02
   34   26   10    1.64   18.24  0.8257E+01  0.8459E+01
   43    1    9    1.66   21.95  0.9940E+01  0.7807E+01
   35   26    8    1.66   20.61  0.9333E+01  0.1339E+02
   19   14   20    1.66   18.52  0.8385E+01  0.5180E+01
   37   26    5    1.67   40.29  0.1824E+02  0.9181E+01
   31    5   17    1.67   39.00  0.1766E+02  0.1078E+02
   17    9   21    1.67   17.11  0.7749E+01  0.5257E+01
   20   17   19    1.67   26.69  0.1208E+02  0.8642E+01
   43    1    8    1.68   19.63  0.8886E+01  0.8495E+01
   34   25    9    1.68   18.78  0.8503E+01  0.8603E+01
   44    0    6    1.69   13.52  0.1559E+02  0.8553E+01
   36   25    6    1.70  114.85 -0.1427E+02  0.1410E+02
   34    3   15    1.70  155.19 -0.1929E+02  0.1001E+02
   37   26    2    1.70   78.67 -0.9776E+01  0.1556E+02
   20   19   18    1.70   44.82 -0.5570E+01  0.6060E+01
   35   25    7    1.71   28.72 -0.3570E+01  0.8827E+01
   43    0    7    1.71   17.87 -0.7815E+01  0.9062E+01
   44    0    5    1.71   18.58 -0.8584E+01  0.8162E+01
   39    1   11    1.73   28.63 -0.3558E+01  0.1114E+02
   30   23   12    1.73   14.26 -0.1772E+01  0.8098E+01
   37   25    1    1.73   44.03 -0.5472E+01  0.8832E+01
   33    3   15    1.73   13.71 -0.1704E+01  0.6575E+01
   27   22   14    1.73   73.96 -0.9191E+01  0.7133E+01
   37   25    0    1.73   14.31 -0.4793E+01  0.1245E+02
   19   18   18    1.74   11.74 -0.1459E+01  0.6361E+01
   25   21   15    1.74   58.81 -0.7308E+01  0.9556E+01
   32   23   10    1.75  185.20 -0.2301E+02  0.1341E+02
   17    8   20    1.75   10.13 -0.1259E+01  0.5028E+01
   44    1    1    1.75   59.77 -0.7427E+01  0.8271E+01
   44    1    0    1.76   14.34 -0.4812E+01  0.1114E+02
   42    0    6    1.77   17.78 -0.7716E+01  0.9435E+01
   35    2   13    1.77   72.05 -0.8954E+01  0.6873E+01
   26   21   14    1.77   55.47 -0.6893E+01  0.7165E+01
   36   24    2    1.78   30.76 -0.3822E+01  0.8709E+01
   15    8   20    1.78   24.30 -0.3020E+01  0.4429E+01
   43    1    3    1.78   16.15 -0.2007E+01  0.8368E+01
   32   23    9    1.78   14.96  0.3301E+01  0.7950E+01
   33   23    8    1.78   14.99  0.3308E+01  0.8238E+01
   36   24    0    1.79   11.08  0.5288E+01  0.1251E+02
   31    3   15    1.79   14.46  0.3190E+01  0.6509E+01
   43    1    2    1.79   11.01  0.2431E+01  0.8192E+01
   40    0    8    1.79   12.85  0.6896E+01  0.7954E+01
   43    1    1    1.79   23.12  0.5103E+01  0.8596E+01
   27   21   13    1.80   45.39  0.1002E+02  0.7555E+01
   31   22   10    1.80   14.03  0.3096E+01  0.1228E+02
   34   23    6    1.81   62.46  0.1378E+02  0.8933E+01
   33   23    7    1.81   24.21  0.5343E+01  0.1238E+02
   30    3   15    1.82   20.22  0.4462E+01  0.8896E+01
   18   18   17    1.82   13.01  0.7065E+01  0.6488E+01
   28   21   12    1.82   61.94  0.1367E+02  0.7931E+01
   41    1    5    1.83   15.98  0.3526E+01  0.8469E+01
   36    1   11    1.83   20.38  0.4499E+01  0.1131E+02
   42    1    1    1.84   21.10  0.4657E+01  0.8513E+01
   37    0   10    1.84   23.91  0.3056E+02  0.1445E+02
   38    0    9    1.84   13.45  0.7522E+01  0.7837E+01
   15   10   19    1.84   22.68  0.5006E+01  0.5738E+01
   29   21   11    1.84   20.85  0.4601E+01  0.7558E+01
   35   23    1    1.84   15.66  0.3457E+01  0.8976E+01
   24   19   14    1.86   72.57  0.1602E+02  0.1046E+02
   41    1    3    1.86   48.26  0.1065E+02  0.8909E+01
   34   22    4    1.87   27.90  0.6158E+01  0.9193E+01
   36    0   10    1.88   11.91  0.5998E+01  0.7557E+01
   30    2   14    1.88   58.47  0.9736E+01  0.6516E+01
   39    1    6    1.89   23.39  0.3894E+01  0.8283E+01
   40    1    4    1.89   34.38  0.5725E+01  0.1434E+02
   21   18   15    1.89   32.01  0.5330E+01  0.6684E+01
   13    6   19    1.90   22.84  0.3803E+01  0.3854E+01
   38    1    7    1.91   87.31  0.1454E+02  0.8218E+01
   29   20   10    1.91   73.16  0.1218E+02  0.1214E+02
   34   22    0    1.91   21.94  0.1789E+02  0.1297E+02
   16   15   17    1.91   10.95  0.1823E+01  0.8207E+01
   14   10   18    1.94   26.69  0.4445E+01  0.5343E+01
   38    1    6    1.94   12.99  0.2162E+01  0.8364E+01
   20   17   15    1.94   12.62  0.2101E+01  0.6538E+01
   33    0   11    1.95   16.21  0.8367E+01  0.1122E+02
   28    2   14    1.95   55.85  0.9300E+01  0.9264E+01
   17   16   16    1.95   25.10  0.4180E+01  0.5715E+01
   37    1    7    1.95   19.92  0.3317E+01  0.8089E+01
   34    0   10    1.96   17.38  0.9795E+01  0.7908E+01
   28   19   10    1.97  145.74  0.2427E+02  0.1547E+02
   39    2    2    1.97   16.81  0.2799E+01  0.8592E+01
   29    1   13    1.98   49.51  0.8244E+01  0.7230E+01
   38    2    4    1.99   38.41  0.6397E+01  0.9521E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        89
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1571
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1392
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      179


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
