
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.673    77.673    38.876  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      49      49
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      49
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001063  0.000636
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
 BEAM_DIVERGENCE_E.S.D.=   0.050
 REFLECTING_RANGE_E.S.D.=   0.213
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

    1057 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      49 ...      49
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1057 REFLECTIONS ACCEPTED



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

 *  44        aP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  1  0  0  0  0 -1  0  0
 *  35        mP          0.0      77.7   38.9   77.7  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  34        mP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  33        mP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  32        oP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  23        oC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.7   77.7   38.9  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     160.1   38.9   77.7  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     160.1   38.9   77.7  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.9  160.1   77.7  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.9  160.1   77.7  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.9  160.1   77.7  90.0  90.0  75.9    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.9  160.1   77.7  90.0  90.0  75.9    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     160.1   38.9  109.8  90.0 133.3  75.9   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.9  109.8  116.5  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.9  160.1  160.1  86.6 104.1 104.1    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.8  116.5   38.9  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        997.9     173.7   77.7   38.9  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        997.9      77.7  173.7   38.9  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        997.9      77.7  173.7   38.9  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        997.9      77.7   77.7   38.9  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.5  116.5  116.5  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.9   86.9  116.5 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.9   86.9  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.9   86.9  116.5 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.8   86.9   86.9  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.9   86.9  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.9   86.9  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.9   86.9  248.7 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.9   86.9   77.7  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.9   86.9   77.7  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.9   77.7  178.0  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.9   86.9  178.0 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.9   86.9   77.7  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.9   86.9   86.9  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     178.0  109.8   38.9  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.9   77.7  178.0 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.9  178.0   77.7 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    4.0      10    0.01  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    74.027   -23.516     0.048
 UNIT_CELL_B-AXIS=    23.516    74.027     0.019
 UNIT_CELL_C-AXIS=    -0.026    -0.002    38.877
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM        88 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.22
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.682    77.682    38.871  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  5.3E-01 5.3E-01 2.4E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012873  0.012873  0.025726  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    74.042   -23.500    -0.014
 COORDINATES OF UNIT CELL B-AXIS    23.500    74.042    -0.008
 COORDINATES OF UNIT CELL C-AXIS     0.008     0.002    38.871
 CRYSTAL MOSAICITY (DEGREES)     0.213
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999979  0.006254  0.001964
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.85   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1260.99   1241.99
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.39
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.68    77.68    38.87  90.000  90.000  90.000


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
  0.05  0.15         6        562.666       0.155
  0.15  0.25       126         34.602       2.526
  0.25  0.35       152         93.967       0.930
  0.35  0.45       153         30.595       2.856
  0.45  0.55       169        144.265       0.606
  0.55  0.65       172         59.591       1.466
  0.65  0.75       171        134.146       0.651
  0.75  0.85       106         76.192       1.147
  0.85  0.95         2         12.041       7.257
 
 -0.05  0.95      1057         87.389       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT               112.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.718
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    48.1 XMAX=    48.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59623 NYBIN=    4
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        240


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT               112.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.716
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    79.5 XMAX=  2450.8 NXBIN=    2
 YMIN=     7.4 YMAX=  2517.3 NYBIN=    2
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        240


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               105.4
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.720
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    48.1 XMAX=    48.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1265    1262
 DETECTOR_SURFACE_POSITION=    1668    1689
 DETECTOR_SURFACE_POSITION=     862    1689
 DETECTOR_SURFACE_POSITION=     862     836
 DETECTOR_SURFACE_POSITION=    1668     836
 DETECTOR_SURFACE_POSITION=    2177    1662
 DETECTOR_SURFACE_POSITION=    1643    2227
 DETECTOR_SURFACE_POSITION=     888    2227
 DETECTOR_SURFACE_POSITION=     354    1662
 DETECTOR_SURFACE_POSITION=     354     863
 DETECTOR_SURFACE_POSITION=     888     298
 DETECTOR_SURFACE_POSITION=    1643     298
 DETECTOR_SURFACE_POSITION=    2177     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        240


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
 8.177E-01  8.991E-04   36.88


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

   24.565   2.880    16.77   1.05      2.06      2.03      66     100      16
    2.880   2.043     1.73   1.11     43.59     40.22     156     212       2
    2.043   1.670     0.06   0.82   1202.65   1274.80     200     256       2
    1.670   1.447    -0.07   1.12    -99.90    -99.90     196     278       4
    1.447   1.295     0.03   0.96   1031.28   1031.44     120     184       2
   --------------------------------------------------------------------------
   24.565   1.295     1.99   1.00     11.29     10.92     738    1030      26
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      49      49)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.689   0.1212E+02  0.1420E+02         0.85           1 
   47    0    0     1.653   0.4274E+01  0.1286E+02         0.33           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     4.9% OF MEAN INTENSITY


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

     3.84          50      36      1297        2.8%       1.7%      1.7%       28   27.52      2.5%    99.9*     0    0.000       0
     2.73          70      46      2129        2.2%       5.0%      5.8%       48   10.23      7.0%    99.5*     0    0.000       0
     2.24         121      76      2699        2.8%      31.0%     28.0%       90    2.53     43.9%    80.7*     0    0.000       0
     1.94         119      74      3139        2.4%     186.9%    181.6%       90    0.49    264.3%     3.3      0    0.000       0
     1.74         144      87      3535        2.5%     438.9%    444.7%      114    0.13    620.6%   -18.9      0    0.000       0
     1.59         149      93      3882        2.4%    1192.4%   1219.0%      112    0.03   1686.4%    -3.5      0    0.000       0
     1.47         160     103      4206        2.4%     -99.9%    -99.9%      114    0.00    -99.9%   -11.0      0    0.000       0
     1.37         145      94      4496        2.1%     -99.9%    -99.9%      102    0.00    -99.9%    12.7      0    0.000       0
     1.30          68      50      4758        1.1%     206.7%    160.7%       36    0.17    292.3%   -40.7      0    0.000       0
    total        1026     659     30141        2.2%      11.2%     10.8%      734    2.56     15.8%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1057
 NUMBER OF REJECTED MISFITS                                26
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1030
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    661



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

   WILSON LINE (using all data) : A=   5.738 B=  27.609 CORRELATION=  0.75
      #      RES      SS        <I>       log(<I>)       BO
      13     7.416   0.005  1.8431E+03       7.519    -195.9
      19     4.727   0.011  1.1977E+03       7.088     -60.3
      17     3.553   0.020  3.8654E+02       5.957      -5.5
      26     2.988   0.028  8.0675E+01       4.390      24.1
      32     2.637   0.036  5.5344E+01       4.014      24.0
      40     2.406   0.043  1.9442E+01       2.967      32.1
      33     2.208   0.051  9.7278E+00       2.275      33.8
      36     2.054   0.059  4.1501E+00       1.423      36.4
      38     1.931   0.067  1.8087E+00       0.593      38.4
      39     1.826   0.075  3.8143E-01      -0.964      44.7
      43     1.740   0.083  2.4415E-02      -3.713      57.2
      46     1.665   0.090 -6.6855E-01    -999.999   -1000.0
      46     1.597   0.098  1.6398E+00       0.495      26.7
      50     1.535   0.106 -1.9198E+00    -999.999   -1000.0
      50     1.481   0.114 -2.5899E+00    -999.999   -1000.0
      49     1.433   0.122 -2.9099E-01    -999.999   -1000.0
      39     1.389   0.130 -4.4850E-02    -999.999   -1000.0
      33     1.350   0.137  8.0302E-01      -0.219      21.7
       9     1.313   0.145  3.0793E+00       1.125      15.9


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       3     7.416        2.525        1.973        1.183
       3     4.727        1.559        1.055        0.531
       1     3.553        0.523        0.131        0.023
       3     2.988        2.698        1.636        0.732
       2     2.637        0.014        0.001        0.000
       2     2.406        0.043        0.001        0.000
       2     2.208        0.784        0.293        0.082
       5     2.054        1.565        0.965        0.512
       8     1.740    37470.323  2843305.266310795736.959
       4     1.597        3.693        3.065        2.424
       4     1.535        5.566        4.902        4.055
       3     1.481       35.350       16.928      141.021
       3     1.433       75.794       99.088      612.163
       1     1.389     1718.877    24686.394   253246.002
      70   overall     4308.468   324812.177 35520863.084


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.416        0.555        0.413        0.241
      16     4.727        1.430        2.054        2.337
      16     3.553        1.273        1.820        2.312
      23     2.988        0.649        0.482        0.297
      30     2.637        1.398        2.052        2.869
      38     2.406        1.025        0.795        0.532
      31     2.208        1.383        1.101        0.946
      31     2.054        2.786        1.845        3.375
      31     1.931       13.975       15.034       99.361
      32     1.826      211.492      762.147    20356.987
      42     1.597       28.538      -69.252      579.863
      32     1.350      101.151      398.994     4513.701
       9     1.313        8.225        4.628       23.589
     588   overall     5422.728   319065.642226119826.605

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       3     7.416  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.49
       3     4.727  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       1     3.553  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.988  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       2     2.637  2.01  1.45  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       2     2.406  2.01  1.45  1.20  1.06  1.92  1.78  1.67  1.59  1.52  1.47
       2     2.208  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.80  0.76  0.73
       5     2.054  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.59
       8     1.740  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       4     1.597  1.01  0.72  0.60  0.53  0.48  0.45  0.42  0.40  0.38  0.37
       4     1.535  3.02  2.17  1.80  1.58  1.44  1.34  1.26  1.19  1.14  1.10
       3     1.481  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.433  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       1     1.389  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      70   overall  1.04  0.74  0.65  0.57  0.55  0.51  0.48  0.48  0.46  0.46


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.416  4.20  3.31  2.31  2.12  1.78  1.55  1.59  1.45  1.35  1.27
      16     4.727  3.28  2.07  1.93  1.52  1.27  1.25  1.37  1.36  1.37  1.38
      16     3.553  1.97  1.03  0.72  1.14  1.11  0.97  0.99  1.13  1.16  1.09
      23     2.988  1.83  1.68  1.34  1.32  1.55  1.45  1.38  1.26  1.17  1.10
      30     2.637  1.05  1.29  1.16  1.31  1.10  1.03  0.99  1.03  1.01  1.11
      38     2.406  1.11  0.73  0.71  0.72  0.80  0.82  0.84  0.91  0.89  0.96
      31     2.208  2.71  1.42  1.24  1.08  0.90  0.79  0.77  0.70  0.76  0.77
      31     2.054  3.39  1.96  1.37  1.08  0.90  0.79  0.77  0.82  0.76  0.71
      31     1.931  4.41  2.31  1.62  1.27  1.07  0.93  0.90  0.82  0.76  0.77
      32     1.826  4.27  2.24  1.57  1.23  1.03  0.90  0.87  0.79  0.74  0.69
      42     1.597  2.75  1.44  1.01  0.79  0.67  0.58  0.52  0.52  0.52  0.53
      32     1.350  5.25  2.76  2.05  1.61  1.35  1.18  1.06  0.96  0.90  0.84
       9     1.313  4.67  2.45  1.71  1.35  1.13  0.99  0.88  0.81  0.75  0.70
     588   overall  1.75  1.03  0.79  0.69  0.61  0.55  0.53  0.52  0.50  0.50


 List of     93 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   36   24    3    1.78  959.09  0.2342E+02  0.1052E+02 "alien"
   36   25    6    1.71  860.77  0.2102E+02  0.1544E+02 "alien"
   25   21   15    1.75  770.42  0.1881E+02  0.1135E+02 "alien"
   34   32   16    1.37  748.48 -0.3357E+02  0.1724E+02 "alien"
   33   24    9    1.74  724.68  0.1769E+02  0.1036E+02 "alien"
   35    2   13    1.78  617.45  0.1507E+02  0.1209E+02 "alien"
   32   31   17    1.39  503.07 -0.2256E+02  0.1559E+02 "alien"
   35   32   15    1.38  502.23 -0.2253E+02  0.1586E+02 "alien"
   31   30   18    1.38  476.90 -0.2139E+02  0.1618E+02 "alien"
   34    3   15    1.71  416.62  0.1017E+02  0.1328E+02 "alien"
   53    4    8    1.40  367.59 -0.1649E+02  0.1018E+02 "alien"
   48    6   14    1.39  362.83 -0.1627E+02  0.1435E+02 "alien"
   29   26   20    1.39  314.94 -0.1413E+02  0.1399E+02 "alien"
   43    7   17    1.41  312.91 -0.1403E+02  0.1347E+02 "alien"
   42    8   18    1.39  300.14 -0.1346E+02  0.1372E+02 "alien"
   54    4    6    1.40  296.06 -0.1328E+02  0.1106E+02 "alien"
   45    7   16    1.40  280.69 -0.1259E+02  0.1567E+02 "alien"
   36   10   21    1.38  253.62 -0.1137E+02  0.1023E+02 "alien"
   27   22   14    1.74  235.15  0.5741E+01  0.8209E+01 "alien"
   35   24    5    1.78  210.90  0.5149E+01  0.9558E+01 "alien"
   37   25    1    1.74  207.74  0.5072E+01  0.8983E+01 "alien"
   39   33   11    1.40  191.27 -0.8579E+01  0.1104E+02 "alien"
   37   33   13    1.39  166.60 -0.7472E+01  0.1161E+02 "alien"
   33    3   15    1.74  157.05  0.3834E+01  0.8512E+01 "alien"
   29   22   12    1.78  142.08  0.3469E+01  0.8355E+01 "alien"
   54    4    7    1.39  118.64 -0.5321E+01  0.9969E+01 "alien"
   28   22   13    1.76  107.14  0.2616E+01  0.8466E+01 "alien"
   26   21   21    1.45  104.20 -0.3032E+02  0.1259E+02 "alien"
   50    5   12    1.40  100.25 -0.4496E+01  0.1112E+02 "alien"
   22   21   17    1.70   97.56  0.2382E+01  0.1102E+02 "alien"
   52    5   11    1.37   76.64 -0.3437E+01  0.1758E+02 "alien"
   37    9   20    1.41   64.12 -0.2876E+01  0.1286E+02 "alien"
   51    3    8    1.45   62.91 -0.1831E+02  0.1615E+02 "alien"
   27   22   21    1.42   56.65 -0.1649E+02  0.9840E+01 "alien"
   47    5   13    1.44   55.16 -0.1605E+02  0.1074E+02 "alien"
   49    6   14    1.37   54.50 -0.2444E+01  0.1741E+02 "alien"
   49    5   12    1.42   49.35 -0.1436E+02  0.1079E+02 "alien"
   20   18   16    1.86   47.17  0.1799E+02  0.1045E+02 "alien"
   31    2   14    1.86   44.59  0.1701E+02  0.8302E+01 "alien"
   38    9   20    1.39   42.88 -0.1923E+01  0.1340E+02 "alien"
   42    7   17    1.43   41.57 -0.1210E+02  0.1380E+02 "alien"
   37   25    2    1.73   38.88  0.9492E+00  0.1005E+02 "alien"
   29   27   19    1.42   37.71 -0.1097E+02  0.1381E+02 "alien"
   44    1    0    1.77   36.52  0.1671E+02  0.1439E+02 "alien"
   40   33    9    1.42   34.84 -0.1014E+02  0.9977E+01 "alien"
   19   18   18    1.75   34.59  0.8445E+00  0.1053E+02 "alien"
   38   32   11    1.43   32.38 -0.9422E+01  0.1009E+02 "alien"
   31   30   17    1.41   32.18 -0.9363E+01  0.1045E+02 "alien"
   37   25    0    1.74   32.18  0.9669E+01  0.1293E+02 "alien"
   43    1    1    1.80   31.87  0.1216E+02  0.8762E+01 "alien"
   33   32   18    1.33   31.50  0.2529E+02  0.1477E+02 "alien"
   27   21   13    1.81   31.38  0.1197E+02  0.8438E+01 "alien"
   41   32    6    1.46   31.35 -0.9123E+01  0.1007E+02 "alien"
   32   30   19    1.34   30.51  0.2450E+02  0.1470E+02 "alien"
   30   29   17    1.44   30.10 -0.8759E+01  0.9589E+01 "alien"
   43    0    7    1.72   28.10  0.5776E+01  0.1487E+02 "alien"
   34   26   10    1.64   27.45 -0.1835E+02  0.1454E+02 "alien"
   35   34   16    1.33   26.11  0.2097E+02  0.1531E+02 "alien"
   53    4    7    1.41   25.40 -0.7391E+01  0.1645E+02 "alien"
   46    1    6    1.63   22.79 -0.1524E+02  0.1369E+02 "alien"
   30   25   21    1.36   22.14  0.1778E+02  0.1056E+02 "alien"
   54    5    9    1.36   22.12  0.1776E+02  0.1024E+02 "alien"
   54    3    1    1.44   22.08 -0.6426E+01  0.1068E+02 "alien"
   36    3   14    1.70   21.51 -0.1438E+02  0.1371E+02 "alien"
   30   22   11    1.80   21.16  0.8072E+01  0.8259E+01 "alien"
   24   20   15    1.79   21.08  0.8039E+01  0.1260E+02 "alien"
   34    4   16    1.66   20.70 -0.1384E+02  0.1112E+02 "alien"
   24   23   17    1.63   20.21 -0.1351E+02  0.1270E+02 "alien"
   40    1    6    1.86   19.11  0.7291E+01  0.9876E+01 "alien"
   34   24    8    1.74   19.08  0.4658E+00  0.9214E+01 "alien"
   35   25    8    1.69   18.36 -0.1228E+02  0.1065E+02 "alien"
   40    2   12    1.66   17.89 -0.1196E+02  0.9039E+01 "alien"
   37   26    4    1.69   17.87 -0.1195E+02  0.9006E+01 "alien"
   32   31   18    1.36   16.52  0.1327E+02  0.1103E+02 "alien"
   54    3    0    1.44   16.42 -0.1503E+02  0.1534E+02 "alien"
   52    3    7    1.44   16.32 -0.4747E+01  0.1150E+02 "alien"
   37   31   12    1.44   16.27 -0.4735E+01  0.1589E+02 "alien"
   37    3   14    1.67   15.63 -0.1045E+02  0.8524E+01 "alien"
   52    4    8    1.42   15.53 -0.4518E+01  0.1532E+02 "alien"
   36   24    2    1.79   13.49  0.5144E+01  0.8798E+01 "alien"
   53    6   11    1.35   13.07  0.1050E+02  0.1197E+02 "alien"
   41    0    7    1.79   12.41  0.1118E+02  0.1493E+02 "alien"
   39    2   12    1.70   12.39 -0.8283E+01  0.8798E+01 "alien"
   33    1   12    1.90   12.34  0.2232E+02  0.1288E+02 "alien"
   45    1    8    1.63   12.27  0.2012E+02  0.1107E+02 "alien"
   39   32   10    1.43   11.80 -0.3432E+01  0.1035E+02 "alien"
   27   24   16    1.61   11.35  0.1861E+02  0.1281E+02 "alien"
   38   26    0    1.69   11.29 -0.1655E+02  0.1568E+02 "alien"
   35   23    2    1.85   10.87  0.4146E+01  0.8696E+01 "alien"
   31   26   21    1.33   10.57  0.8485E+01  0.1408E+02 "alien"
   27   26   17    1.54   10.27 -0.1972E+02  0.1473E+02 "alien"
   44    2   10    1.61   10.19  0.1671E+02  0.9098E+01 "alien"
   38   29    7    1.56   10.02 -0.1923E+02  0.1655E+02 "alien"


 List of     93 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   33   32   18    1.33   31.50  0.2529E+02  0.1477E+02
   35   34   16    1.33   26.11  0.2097E+02  0.1531E+02
   31   26   21    1.33   10.57  0.8485E+01  0.1408E+02
   32   30   19    1.34   30.51  0.2450E+02  0.1470E+02
   53    6   11    1.35   13.07  0.1050E+02  0.1197E+02
   30   25   21    1.36   22.14  0.1778E+02  0.1056E+02
   32   31   18    1.36   16.52  0.1327E+02  0.1103E+02
   54    5    9    1.36   22.12  0.1776E+02  0.1024E+02
   49    6   14    1.37   54.50 -0.2444E+01  0.1741E+02
   52    5   11    1.37   76.64 -0.3437E+01  0.1758E+02
   34   32   16    1.37  748.48 -0.3357E+02  0.1724E+02
   36   10   21    1.38  253.62 -0.1137E+02  0.1023E+02
   31   30   18    1.38  476.90 -0.2139E+02  0.1618E+02
   35   32   15    1.38  502.23 -0.2253E+02  0.1586E+02
   32   31   17    1.39  503.07 -0.2256E+02  0.1559E+02
   37   33   13    1.39  166.60 -0.7472E+01  0.1161E+02
   54    4    7    1.39  118.64 -0.5321E+01  0.9969E+01
   48    6   14    1.39  362.83 -0.1627E+02  0.1435E+02
   38    9   20    1.39   42.88 -0.1923E+01  0.1340E+02
   42    8   18    1.39  300.14 -0.1346E+02  0.1372E+02
   29   26   20    1.39  314.94 -0.1413E+02  0.1399E+02
   50    5   12    1.40  100.25 -0.4496E+01  0.1112E+02
   45    7   16    1.40  280.69 -0.1259E+02  0.1567E+02
   39   33   11    1.40  191.27 -0.8579E+01  0.1104E+02
   53    4    8    1.40  367.59 -0.1649E+02  0.1018E+02
   54    4    6    1.40  296.06 -0.1328E+02  0.1106E+02
   43    7   17    1.41  312.91 -0.1403E+02  0.1347E+02
   37    9   20    1.41   64.12 -0.2876E+01  0.1286E+02
   53    4    7    1.41   25.40 -0.7391E+01  0.1645E+02
   31   30   17    1.41   32.18 -0.9363E+01  0.1045E+02
   40   33    9    1.42   34.84 -0.1014E+02  0.9977E+01
   29   27   19    1.42   37.71 -0.1097E+02  0.1381E+02
   49    5   12    1.42   49.35 -0.1436E+02  0.1079E+02
   52    4    8    1.42   15.53 -0.4518E+01  0.1532E+02
   27   22   21    1.42   56.65 -0.1649E+02  0.9840E+01
   42    7   17    1.43   41.57 -0.1210E+02  0.1380E+02
   38   32   11    1.43   32.38 -0.9422E+01  0.1009E+02
   39   32   10    1.43   11.80 -0.3432E+01  0.1035E+02
   54    3    1    1.44   22.08 -0.6426E+01  0.1068E+02
   54    3    0    1.44   16.42 -0.1503E+02  0.1534E+02
   47    5   13    1.44   55.16 -0.1605E+02  0.1074E+02
   52    3    7    1.44   16.32 -0.4747E+01  0.1150E+02
   37   31   12    1.44   16.27 -0.4735E+01  0.1589E+02
   30   29   17    1.44   30.10 -0.8759E+01  0.9589E+01
   26   21   21    1.45  104.20 -0.3032E+02  0.1259E+02
   51    3    8    1.45   62.91 -0.1831E+02  0.1615E+02
   41   32    6    1.46   31.35 -0.9123E+01  0.1007E+02
   27   26   17    1.54   10.27 -0.1972E+02  0.1473E+02
   38   29    7    1.56   10.02 -0.1923E+02  0.1655E+02
   44    2   10    1.61   10.19  0.1671E+02  0.9098E+01
   27   24   16    1.61   11.35  0.1861E+02  0.1281E+02
   45    1    8    1.63   12.27  0.2012E+02  0.1107E+02
   46    1    6    1.63   22.79 -0.1524E+02  0.1369E+02
   24   23   17    1.63   20.21 -0.1351E+02  0.1270E+02
   34   26   10    1.64   27.45 -0.1835E+02  0.1454E+02
   34    4   16    1.66   20.70 -0.1384E+02  0.1112E+02
   40    2   12    1.66   17.89 -0.1196E+02  0.9039E+01
   37    3   14    1.67   15.63 -0.1045E+02  0.8524E+01
   38   26    0    1.69   11.29 -0.1655E+02  0.1568E+02
   37   26    4    1.69   17.87 -0.1195E+02  0.9006E+01
   35   25    8    1.69   18.36 -0.1228E+02  0.1065E+02
   39    2   12    1.70   12.39 -0.8283E+01  0.8798E+01
   36    3   14    1.70   21.51 -0.1438E+02  0.1371E+02
   22   21   17    1.70   97.56  0.2382E+01  0.1102E+02
   36   25    6    1.71  860.77  0.2102E+02  0.1544E+02
   34    3   15    1.71  416.62  0.1017E+02  0.1328E+02
   43    0    7    1.72   28.10  0.5776E+01  0.1487E+02
   37   25    2    1.73   38.88  0.9492E+00  0.1005E+02
   37   25    1    1.74  207.74  0.5072E+01  0.8983E+01
   33    3   15    1.74  157.05  0.3834E+01  0.8512E+01
   27   22   14    1.74  235.15  0.5741E+01  0.8209E+01
   37   25    0    1.74   32.18  0.9669E+01  0.1293E+02
   33   24    9    1.74  724.68  0.1769E+02  0.1036E+02
   34   24    8    1.74   19.08  0.4658E+00  0.9214E+01
   19   18   18    1.75   34.59  0.8445E+00  0.1053E+02
   25   21   15    1.75  770.42  0.1881E+02  0.1135E+02
   28   22   13    1.76  107.14  0.2616E+01  0.8466E+01
   44    1    0    1.77   36.52  0.1671E+02  0.1439E+02
   36   24    3    1.78  959.09  0.2342E+02  0.1052E+02
   35    2   13    1.78  617.45  0.1507E+02  0.1209E+02
   35   24    5    1.78  210.90  0.5149E+01  0.9558E+01
   29   22   12    1.78  142.08  0.3469E+01  0.8355E+01
   36   24    2    1.79   13.49  0.5144E+01  0.8798E+01
   41    0    7    1.79   12.41  0.1118E+02  0.1493E+02
   24   20   15    1.79   21.08  0.8039E+01  0.1260E+02
   30   22   11    1.80   21.16  0.8072E+01  0.8259E+01
   43    1    1    1.80   31.87  0.1216E+02  0.8762E+01
   27   21   13    1.81   31.38  0.1197E+02  0.8438E+01
   35   23    2    1.85   10.87  0.4146E+01  0.8696E+01
   31    2   14    1.86   44.59  0.1701E+02  0.8302E+01
   20   18   16    1.86   47.17  0.1799E+02  0.1045E+02
   40    1    6    1.86   19.11  0.7291E+01  0.9876E+01
   33    1   12    1.90   12.34  0.2232E+02  0.1288E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        51
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1057
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)     957
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      100


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
