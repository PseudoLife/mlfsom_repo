
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.669    77.669    38.886  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      50      50
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      50
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.001053  0.000641
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
 REFLECTING_RANGE_E.S.D.=   0.220
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

    1020 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      50 ...      50
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1020 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               31
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

 *  44        aP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  31        aP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  20        mC          0.0     109.8  109.8   38.9  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.8  109.8   38.9  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  35        mP          0.0      77.7   38.9   77.7  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  33        mP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  21        tP          0.0      77.7   77.7   38.9  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
 *  32        oP          0.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    37        mC        250.0     160.1   38.9   77.7  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     160.1   38.9   77.7  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.9  160.1   77.7  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.9  160.1   77.7  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.9  160.1   77.7  90.0  90.0  75.9    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.9  160.1   77.7  90.0  90.0  75.9    0  0 -1  0 -2  0 -1  0  0  1  0  0
    19        oI        500.0      38.9  109.8  116.5  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    27        mC        500.0     160.1   38.9  109.8  90.0 133.3  75.9   -2  0  1  0  0  0  1  0  1  1  0  0
    18        tI        625.0     109.8  116.5   38.9  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    26        oF        625.0      38.9  160.1  160.1  86.6 104.1 104.1    0  0 -1  0 -2  0  1  0  0  2  1  0
    22        hP        997.4      77.7   77.7   38.9  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    30        mC        997.4      77.7  173.7   38.9  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        997.4      77.7  173.7   38.9  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    41        mC        997.4     173.7   77.7   38.9  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.5  116.5  116.5  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.9   86.9  116.5 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     3        cP        999.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.9   86.9  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.9   86.9  116.5 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.8   86.9   86.9  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.9   86.9  109.8  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.9   86.9  109.8  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.9   86.9  248.7 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.9   86.9   77.7  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.9   77.7   77.7  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.9   86.9   77.7  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.9   77.7  178.0  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.9   86.9  178.0 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.9   86.9   77.7  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.9   86.9   86.9  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     178.0  109.8   38.9  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.9   77.7  178.0 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.9  178.0   77.7 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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
 UNIT_CELL_A-AXIS=    74.051   -23.428    -0.004
 UNIT_CELL_B-AXIS=    23.428    74.051     0.006
 UNIT_CELL_C-AXIS=     0.001    -0.004    38.886
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM        85 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.31
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.673    77.673    38.879  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.7E-01 2.7E-01 1.2E-01 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012874  0.012874  0.025721  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    74.054   -23.434     0.014
 COORDINATES OF UNIT CELL B-AXIS    23.434    74.054     0.002
 COORDINATES OF UNIT CELL C-AXIS    -0.007     0.001    38.879
 CRYSTAL MOSAICITY (DEGREES)     0.220
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999975  0.006595  0.002594
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.81   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1260.04   1242.48
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.33
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.67    77.67    38.88  90.000  90.000  90.000


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
  0.05  0.15         1          3.135      27.718
  0.15  0.25       122         83.263       1.044
  0.25  0.35       143         56.629       1.534
  0.35  0.45       158         99.940       0.869
  0.45  0.55       157         86.721       1.002
  0.55  0.65       178         63.857       1.361
  0.65  0.75       169        132.400       0.656
  0.75  0.85        92         78.481       1.107
 
 -0.05  0.95      1020         86.889       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        5
 DEGREES OF FREEDOM OF CHI^2 FIT               126.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.647
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    49.1 XMAX=    49.8 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59617 NYBIN=    5
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        266


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT               126.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.642
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    80.5 XMAX=  2399.7 NXBIN=    2
 YMIN=     6.4 YMAX=  2516.2 NYBIN=    2
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        266


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               123.8
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.635
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    49.1 XMAX=    49.8 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1240    1261
 DETECTOR_SURFACE_POSITION=    1634    1688
 DETECTOR_SURFACE_POSITION=     846    1688
 DETECTOR_SURFACE_POSITION=     846     835
 DETECTOR_SURFACE_POSITION=    1634     835
 DETECTOR_SURFACE_POSITION=    2132    1661
 DETECTOR_SURFACE_POSITION=    1609    2226
 DETECTOR_SURFACE_POSITION=     871    2226
 DETECTOR_SURFACE_POSITION=     349    1661
 DETECTOR_SURFACE_POSITION=     349     862
 DETECTOR_SURFACE_POSITION=     871     297
 DETECTOR_SURFACE_POSITION=    1609     297
 DETECTOR_SURFACE_POSITION=    2132     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        266


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
 6.633E-01  1.137E-03   36.41


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

   24.562   2.580    13.33   1.15      2.35      2.23      90     141      14
    2.580   1.829     0.80   0.89     72.66     77.41     206     275       4
    1.829   1.495    -0.00   0.97    857.02    878.01     236     329       2
    1.495   1.295     0.09   1.11   1695.39   1665.06     160     252       2
   --------------------------------------------------------------------------
   24.562   1.295     2.13   1.00      9.49      9.56     692     997      22
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      50      50)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.689  -0.3342E+01  0.1196E+02        -0.28           1 
   47    0    0     1.653  -0.3930E+01  0.1214E+02        -0.32           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -4.6% OF MEAN INTENSITY


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

     3.84          50      35      1297        2.7%       1.8%      1.7%       30   28.20      2.6%    99.9*     0    0.000       0
     2.73          68      46      2129        2.2%       4.9%      5.4%       44   10.72      6.9%    99.7*     0    0.000       0
     2.24         119      76      2697        2.8%      26.2%     27.7%       86    2.63     37.1%    89.1*     0    0.000       0
     1.94         113      71      3139        2.3%     123.1%    140.3%       84    0.44    174.1%    -4.6      0    0.000       0
     1.74         137      84      3537        2.4%     -99.9%    -99.9%      106    0.00    -99.9%     2.6      0    0.000       0
     1.59         143      89      3879        2.3%     498.0%    479.3%      108    0.14    704.2%     6.8      0    0.000       0
     1.47         155     104      4205        2.5%     359.5%    374.6%      102    0.02    508.5%    18.3      0    0.000       0
     1.37         140      93      4500        2.1%     -99.9%    -99.9%       94    0.00    -99.9%     7.0      0    0.000       0
     1.30          69      51      4751        1.1%     941.9%    962.0%       36    0.17   1332.1%    -8.6      0    0.000       0
    total         994     649     30134        2.2%       9.4%      9.5%      690    2.67     13.3%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1020
 NUMBER OF REJECTED MISFITS                                22
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                          997
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    651



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

   WILSON LINE (using all data) : A=   5.755 B=  24.626 CORRELATION=  0.81
      #      RES      SS        <I>       log(<I>)       BO
      13     7.667   0.004  1.6374E+03       7.401    -193.5
      18     4.727   0.011  1.2460E+03       7.128     -61.4
      18     3.551   0.020  3.6871E+02       5.910      -3.9
      25     2.978   0.028  7.5400E+01       4.323      25.4
      31     2.641   0.036  5.4446E+01       3.997      24.5
      40     2.415   0.043  2.0392E+01       3.015      32.0
      33     2.211   0.051  7.8008E+00       2.054      36.2
      34     2.052   0.059  3.0189E+00       1.105      39.2
      38     1.931   0.067  2.1654E+00       0.773      37.2
      38     1.825   0.075  3.1187E-01      -1.165      46.1
      40     1.739   0.083 -7.7921E-01    -999.999   -1000.0
      45     1.666   0.090 -1.0319E+00    -999.999   -1000.0
      46     1.596   0.098  9.4014E-01      -0.062      29.6
      50     1.536   0.106 -1.2670E+00    -999.999   -1000.0
      49     1.481   0.114  3.2301E+00       1.173      20.1
      47     1.433   0.122  6.7559E-01      -0.392      25.3
      40     1.390   0.129  2.5150E+00       0.922      18.7
      34     1.351   0.137  3.1162E+00       1.137      16.9
      10     1.311   0.146 -1.9962E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       2     7.667        4.256        4.057        2.840
       3     4.727        1.429        0.921        0.442
       1     3.551        0.551        0.142        0.026
       3     2.978        2.770        1.686        0.754
       3     2.415        0.151        0.029        0.004
       2     2.211        0.416        0.131        0.030
       4     2.052        5.597        3.921        5.805
       7     1.931       17.219       32.709       47.365
       4     1.536       67.395      129.973      801.392
       3     1.481        3.762        2.074        2.107
       1     1.390      228.266     1194.681     4466.166
      70   overall       35.963     -164.181     1338.537


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.667        0.654        0.556        0.371
      15     4.727        1.381        1.887        2.041
      17     3.551        1.329        2.009        2.698
      22     2.978        0.580        0.397        0.231
      30     2.641        1.481        2.447        3.854
      37     2.415        1.129        0.997        0.780
      31     2.211        1.393        1.093        0.858
      30     2.052        4.908        3.813       11.511
      30     1.825      435.595     -170.529   126322.366
      42     1.596       60.349      175.624     2658.531
      46     1.481        7.929        9.131       27.095
      44     1.433      130.760      578.324     9165.664
      39     1.390       17.743        4.880      126.996
      33     1.351        9.950       24.172       91.631
     579   overall       40.699       50.298     7455.882

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       2     7.667  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     4.727  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       1     3.551  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.978  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.415  0.00  0.96  0.80  0.70  1.28  1.19  1.12  1.06  1.01  0.98
       2     2.211  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       4     2.052  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       7     1.931  0.58  0.41  0.34  0.30  0.27  0.25  0.24  0.23  0.22  0.21
       4     1.536  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       3     1.481  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       1     1.390  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      70   overall  0.46  0.37  0.31  0.27  0.27  0.25  0.24  0.23  0.22  0.23


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.667  2.87  3.51  2.46  2.21  1.85  1.61  1.44  1.49  1.38  1.29
      15     4.727  3.50  2.21  1.80  1.42  1.19  1.18  1.46  1.45  1.46  1.37
      17     3.551  1.85  1.30  0.91  1.07  1.20  1.04  1.05  1.07  1.19  1.12
      22     2.978  1.43  1.76  1.40  1.52  1.27  1.21  1.17  1.24  1.15  1.22
      30     2.641  1.05  1.47  1.16  1.11  0.93  0.96  0.93  0.91  1.01  1.05
      37     2.415  2.27  1.49  1.04  0.98  0.89  0.84  0.81  0.79  0.82  0.86
      31     2.211  2.71  1.42  1.00  0.88  0.90  0.79  0.77  0.76  0.71  0.66
      30     2.052  3.50  1.84  1.54  1.21  1.02  0.89  0.79  0.73  0.73  0.69
      30     1.825  4.90  2.57  1.80  1.42  1.19  1.03  0.93  0.85  0.79  0.74
      42     1.596  4.75  2.50  1.75  1.37  1.15  1.00  0.90  0.82  0.76  0.72
      46     1.481  4.80  2.64  1.85  1.52  1.27  1.11  1.04  0.95  0.88  0.86
      44     1.433  5.73  3.01  2.19  1.72  1.50  1.31  1.17  1.07  1.00  0.93
      39     1.390  5.12  2.69  1.88  1.56  1.30  1.14  1.02  0.93  0.86  0.81
      33     1.351  4.78  2.51  1.75  1.47  1.31  1.14  1.14  1.05  0.97  0.91
     579   overall  2.81  1.65  1.19  1.01  0.88  0.78  0.74  0.70  0.68  0.66


 List of     57 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   40    1    6    1.86   75.27  0.2347E+02  0.9494E+01 "alien"
   35   23    2    1.85   49.46  0.1542E+02  0.8213E+01 "alien"
   46    5   14    1.44   45.04  0.3043E+02  0.1473E+02 "alien"
   40   29    3    1.56   35.97 -0.4558E+02  0.1536E+02 "alien"
   40    3   13    1.63   35.48  0.3335E+02  0.1187E+02 "alien"
   45    6   15    1.43   33.15  0.2239E+02  0.1339E+02 "alien"
   49    5   12    1.42   31.34  0.2118E+02  0.1027E+02 "alien"
   43    1    1    1.80   31.11  0.9703E+01  0.7672E+01 "alien"
   20   18   16    1.86   30.39  0.9478E+01  0.9564E+01 "alien"
   51    4   10    1.41   27.09  0.1830E+02  0.1583E+02 "alien"
   29   21   11    1.85   26.98  0.8414E+01  0.7793E+01 "alien"
   39   28    2    1.61   26.04  0.2448E+02  0.8542E+01 "alien"
   27   22   21    1.42   24.65  0.1666E+02  0.1296E+02 "alien"
   33   25   11    1.66   24.23 -0.2500E+02  0.1187E+02 "alien"
   41    1    4    1.86   23.42  0.7303E+01  0.8336E+01 "alien"
   43    1    2    1.80   23.21  0.7238E+01  0.8002E+01 "alien"
   25   20   14    1.83   23.14  0.7215E+01  0.7744E+01 "alien"
   33   23    8    1.79   22.66  0.7067E+01  0.7762E+01 "alien"
   35   26    9    1.65   21.87 -0.2256E+02  0.7792E+01 "alien"
   32   30   16    1.43   19.74  0.1334E+02  0.9286E+01 "alien"
   48    3    9    1.51   19.60 -0.2483E+02  0.1471E+02 "alien"
   32   23    9    1.79   19.58  0.6105E+01  0.1462E+02 "alien"
   41   33    7    1.43   19.51  0.1318E+02  0.9422E+01 "alien"
   25   21   15    1.75   19.03 -0.1483E+02  0.9873E+01 "alien"
   34    6   18    1.56   18.83 -0.2386E+02  0.1077E+02 "alien"
   42    7   17    1.43   18.68  0.1262E+02  0.1325E+02 "alien"
   31   27   15    1.53   17.27 -0.2188E+02  0.1241E+02 "alien"
   37    8   19    1.45   16.84  0.1138E+02  0.1311E+02 "alien"
   29   22   12    1.78   15.69 -0.1222E+02  0.7741E+01 "alien"
   49    7   15    1.34   15.61  0.4865E+02  0.1736E+02 "alien"
   37   27    6    1.64   15.46 -0.1595E+02  0.1280E+02 "alien"
   28   21   12    1.83   15.39  0.4801E+01  0.7685E+01 "alien"
   32   23   10    1.76   14.87 -0.1159E+02  0.1227E+02 "alien"
   36   26    7    1.67   14.80 -0.1527E+02  0.7760E+01 "alien"
   37   25    1    1.74   13.40 -0.1044E+02  0.8110E+01 "alien"
   46    3   11    1.52   13.26 -0.1680E+02  0.9158E+01 "alien"
   28   22   13    1.76   13.18 -0.1027E+02  0.7926E+01 "alien"
   53    4    7    1.41   13.14  0.8878E+01  0.1458E+02 "alien"
   53    3    5    1.44   13.01  0.8790E+01  0.1228E+02 "alien"
   47    1    5    1.62   12.61  0.1186E+02  0.9074E+01 "alien"
   50    1    2    1.55   12.58 -0.1593E+02  0.1012E+02 "alien"
   49    4   11    1.44   12.56  0.8486E+01  0.9895E+01 "alien"
   35   25    7    1.72   12.51 -0.9750E+01  0.8933E+01 "alien"
   40    1   10    1.74   12.43 -0.9685E+01  0.7861E+01 "alien"
   38    9   20    1.39   12.18  0.3063E+02  0.1288E+02 "alien"
   27   26   17    1.54   12.02 -0.1523E+02  0.1402E+02 "alien"
   34   11   22    1.37   12.02  0.3022E+02  0.1278E+02 "alien"
   55    3    0    1.41   11.66  0.6581E+02  0.1396E+02 "alien"
   53    3    4    1.45   11.62  0.7849E+01  0.1009E+02 "alien"
   27   24   16    1.61   11.38  0.1070E+02  0.1166E+02 "alien"
   27   22   14    1.74   11.36 -0.8854E+01  0.7713E+01 "alien"
   35    1   12    1.83   11.00  0.3431E+01  0.1085E+02 "alien"
   36   27    9    1.60   10.90  0.1025E+02  0.9224E+01 "alien"
   32   27   14    1.54   10.85 -0.1375E+02  0.1186E+02 "alien"
   36   24    2    1.79   10.75  0.3352E+01  0.8356E+01 "alien"
   34   27   12    1.57   10.33  0.9713E+01  0.1463E+02 "alien"
   24   23   17    1.63   10.32 -0.1065E+02  0.1186E+02 "alien"


 List of     57 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   49    7   15    1.34   15.61  0.4865E+02  0.1736E+02
   34   11   22    1.37   12.02  0.3022E+02  0.1278E+02
   38    9   20    1.39   12.18  0.3063E+02  0.1288E+02
   55    3    0    1.41   11.66  0.6581E+02  0.1396E+02
   53    4    7    1.41   13.14  0.8878E+01  0.1458E+02
   51    4   10    1.41   27.09  0.1830E+02  0.1583E+02
   49    5   12    1.42   31.34  0.2118E+02  0.1027E+02
   27   22   21    1.42   24.65  0.1666E+02  0.1296E+02
   42    7   17    1.43   18.68  0.1262E+02  0.1325E+02
   41   33    7    1.43   19.51  0.1318E+02  0.9422E+01
   45    6   15    1.43   33.15  0.2239E+02  0.1339E+02
   32   30   16    1.43   19.74  0.1334E+02  0.9286E+01
   46    5   14    1.44   45.04  0.3043E+02  0.1473E+02
   53    3    5    1.44   13.01  0.8790E+01  0.1228E+02
   49    4   11    1.44   12.56  0.8486E+01  0.9895E+01
   53    3    4    1.45   11.62  0.7849E+01  0.1009E+02
   37    8   19    1.45   16.84  0.1138E+02  0.1311E+02
   48    3    9    1.51   19.60 -0.2483E+02  0.1471E+02
   46    3   11    1.52   13.26 -0.1680E+02  0.9158E+01
   31   27   15    1.53   17.27 -0.2188E+02  0.1241E+02
   27   26   17    1.54   12.02 -0.1523E+02  0.1402E+02
   32   27   14    1.54   10.85 -0.1375E+02  0.1186E+02
   50    1    2    1.55   12.58 -0.1593E+02  0.1012E+02
   34    6   18    1.56   18.83 -0.2386E+02  0.1077E+02
   40   29    3    1.56   35.97 -0.4558E+02  0.1536E+02
   34   27   12    1.57   10.33  0.9713E+01  0.1463E+02
   36   27    9    1.60   10.90  0.1025E+02  0.9224E+01
   27   24   16    1.61   11.38  0.1070E+02  0.1166E+02
   39   28    2    1.61   26.04  0.2448E+02  0.8542E+01
   47    1    5    1.62   12.61  0.1186E+02  0.9074E+01
   40    3   13    1.63   35.48  0.3335E+02  0.1187E+02
   24   23   17    1.63   10.32 -0.1065E+02  0.1186E+02
   37   27    6    1.64   15.46 -0.1595E+02  0.1280E+02
   35   26    9    1.65   21.87 -0.2256E+02  0.7792E+01
   33   25   11    1.66   24.23 -0.2500E+02  0.1187E+02
   36   26    7    1.67   14.80 -0.1527E+02  0.7760E+01
   35   25    7    1.72   12.51 -0.9750E+01  0.8933E+01
   40    1   10    1.74   12.43 -0.9685E+01  0.7861E+01
   37   25    1    1.74   13.40 -0.1044E+02  0.8110E+01
   27   22   14    1.74   11.36 -0.8854E+01  0.7713E+01
   25   21   15    1.75   19.03 -0.1483E+02  0.9873E+01
   32   23   10    1.76   14.87 -0.1159E+02  0.1227E+02
   28   22   13    1.76   13.18 -0.1027E+02  0.7926E+01
   29   22   12    1.78   15.69 -0.1222E+02  0.7741E+01
   36   24    2    1.79   10.75  0.3352E+01  0.8356E+01
   32   23    9    1.79   19.58  0.6105E+01  0.1462E+02
   33   23    8    1.79   22.66  0.7067E+01  0.7762E+01
   43    1    2    1.80   23.21  0.7238E+01  0.8002E+01
   43    1    1    1.80   31.11  0.9703E+01  0.7672E+01
   25   20   14    1.83   23.14  0.7215E+01  0.7744E+01
   35    1   12    1.83   11.00  0.3431E+01  0.1085E+02
   28   21   12    1.83   15.39  0.4801E+01  0.7685E+01
   35   23    2    1.85   49.46  0.1542E+02  0.8213E+01
   29   21   11    1.85   26.98  0.8414E+01  0.7793E+01
   41    1    4    1.86   23.42  0.7303E+01  0.8336E+01
   20   18   16    1.86   30.39  0.9478E+01  0.9564E+01
   40    1    6    1.86   75.27  0.2347E+02  0.9494E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        18
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1020
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)     970
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       50


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
