
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.501    77.501    38.783  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      26      26
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      26
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000906  0.000695
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
 BEAM_DIVERGENCE_E.S.D.=   0.045
 REFLECTING_RANGE_E.S.D.=   0.129
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

    1349 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      26 ...      26
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1349 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               42
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

 *  44        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  1  0  0  0  0 -1  0  0
 *  35        mP          0.0      77.5   38.8   77.5  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  34        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  33        mP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.6  109.6   38.8  90.0  90.0  90.0    1 -1  0  0  1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.8   38.8   77.5  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.8   38.8   77.5  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  159.8   77.5  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  159.8   77.5  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  159.8   77.5  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  159.8   77.5  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     159.8   38.8  109.6  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.6  116.3  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  159.8  159.8  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.6  116.3   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        998.3     173.3   77.5   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.3      77.5  173.3   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.3      77.5  173.3   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.3      77.5   77.5   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.3  116.3  116.3  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.7   86.7  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.7   86.7  116.3 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.6   86.7   86.7  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.7   86.7  109.6  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.7   86.7  109.6  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.7  248.1 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.7   86.7   77.5  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.8   77.5   77.5  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.5  177.6  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.7   86.7  177.6 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.7   86.7   77.5  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.7   86.7   86.7  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.6  109.6   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.5  177.6 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.7  177.6   77.5 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    4.3      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.859   -23.478     0.048
 UNIT_CELL_B-AXIS=    23.478    73.859     0.020
 UNIT_CELL_C-AXIS=    -0.026    -0.002    38.784
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       208 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.30
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.10
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.443    77.443    38.749  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  9.1E-02 9.1E-02 4.4E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012913  0.012913  0.025807  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.808   -23.447     0.054
 COORDINATES OF UNIT CELL B-AXIS    23.447    73.808     0.020
 COORDINATES OF UNIT CELL C-AXIS    -0.029    -0.001    38.749
 CRYSTAL MOSAICITY (DEGREES)     0.129
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999974  0.006815  0.002184
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.88   1232.92
 DETECTOR ORIGIN (PIXELS) AT                     1260.70   1242.84
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.17
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.44    77.44    38.75  90.000  90.000  90.000


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
  0.05  0.15        65         97.091       1.413
  0.15  0.25       168        159.119       0.862
  0.25  0.35       185        116.771       1.175
  0.35  0.45       158         92.478       1.483
  0.45  0.55       167        196.277       0.699
  0.55  0.65       189        128.794       1.065
  0.65  0.75       193        136.847       1.002
  0.75  0.85       157        157.684       0.870
  0.85  0.95        67        112.264       1.222
 
 -0.05  0.95      1349        137.188       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               212.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.915
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    25.1 XMAX=    25.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59475 NYBIN=    9
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        472


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               212.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.910
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    10.3 XMAX=  2454.5 NXBIN=    3
 YMIN=     8.7 YMAX=  2519.9 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        472


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               212.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.902
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    25.1 XMAX=    25.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1648    1691
 DETECTOR_SURFACE_POSITION=     817    1691
 DETECTOR_SURFACE_POSITION=     817     838
 DETECTOR_SURFACE_POSITION=    1648     838
 DETECTOR_SURFACE_POSITION=    2172    1664
 DETECTOR_SURFACE_POSITION=    1622    2229
 DETECTOR_SURFACE_POSITION=     843    2229
 DETECTOR_SURFACE_POSITION=     293    1664
 DETECTOR_SURFACE_POSITION=     293     865
 DETECTOR_SURFACE_POSITION=     843     299
 DETECTOR_SURFACE_POSITION=    1622     299
 DETECTOR_SURFACE_POSITION=    2172     865
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        472


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
 7.120E-01  9.637E-04   38.17


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

   24.490   3.154    27.53   1.64      1.98      1.66      62      89      28
    3.154   2.240    10.21   1.13      5.05      4.55     178     209      16
    2.240   1.831     1.52   1.06     34.88     32.42     196     231       2
    1.831   1.587    -0.10   0.84    -99.90    -99.90     238     273       2
    1.587   1.420    -0.00   0.91    -99.90    -99.90     270     321       2
    1.420   1.297    -0.03   0.68   8369.04   9683.78     122     169       6
   --------------------------------------------------------------------------
   24.490   1.297     3.80   0.98      6.96      6.74    1066    1292      56
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      26      26)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.684  -0.1318E+02  0.1087E+02        -1.21           1 
   47    0    0     1.648   0.1952E+02  0.1122E+02         1.74           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    13.7% OF MEAN INTENSITY


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

     3.85          55      35      1285        2.7%       1.9%      1.6%       40   36.32      2.6%   100.0*     0    0.000       0
     2.74          86      53      2098        2.5%       3.3%      2.7%       66   24.25      4.6%    99.6*     0    0.000       0
     2.24         157      90      2665        3.4%       6.1%      5.7%      134   11.61      8.7%    99.1*     0    0.000       0
     1.94         150      86      3099        2.8%      24.4%     21.6%      128    3.04     34.4%    90.7*     0    0.000       0
     1.74         159      90      3482        2.6%     -99.9%    -99.9%      138    0.00    -99.9%    -0.1      0    0.000       0
     1.59         195     111      3850        2.9%     -99.9%    -99.9%      168    0.00    -99.9%    26.6      0    0.000       0
     1.47         207     121      4144        2.9%     -99.9%    -99.9%      172    0.00    -99.9%    12.7      0    0.000       0
     1.38         202     116      4427        2.6%     -99.9%    -99.9%      172    0.00    -99.9%    21.2      0    0.000       0
     1.30          80      56      4741        1.2%     375.6%    454.0%       48   -0.00    531.1%    18.0      0    0.000       0
    total        1291     758     29791        2.5%       7.0%      6.7%     1066    5.07      9.8%   100.0*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1349
 NUMBER OF REJECTED MISFITS                                56
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1292
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    759



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

   WILSON LINE (using all data) : A=   6.897 B=  30.251 CORRELATION=  0.77
      #      RES      SS        <I>       log(<I>)       BO
      13     7.281   0.005  2.0839E+03       7.642     -79.0
      18     4.606   0.012  2.0639E+03       7.632     -31.2
      25     3.510   0.020  8.2862E+02       6.720       4.4
      26     3.008   0.028  2.7568E+02       5.619      23.1
      33     2.639   0.036  2.3682E+02       5.467      19.9
      49     2.401   0.043  1.1511E+02       4.746      24.8
      38     2.210   0.051  6.1322E+01       4.116      27.2
      43     2.056   0.059  2.1081E+01       3.048      32.6
      42     1.932   0.067  5.8188E-01      -0.541      55.6
      48     1.830   0.075  3.2722E-02      -3.420      69.1
      42     1.740   0.083 -8.0595E-01    -999.999   -1000.0
      49     1.668   0.090 -1.6553E+00    -999.999   -1000.0
      52     1.598   0.098 -1.7375E-01    -999.999   -1000.0
      59     1.538   0.106 -7.4643E-01    -999.999   -1000.0
      69     1.481   0.114  1.0078E+00       0.008      30.2
      52     1.434   0.122 -1.0515E+00    -999.999   -1000.0
      44     1.392   0.129  1.5454E+00       0.435      25.0
      39     1.352   0.137 -1.2151E+00    -999.999   -1000.0
      17     1.315   0.144  2.2701E+00       0.820      21.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.281        1.692        1.157        0.619
       3     4.606        1.494        0.999        0.496
       2     3.510        0.401        0.088        0.014
       3     3.008        1.381        0.677        0.245
       2     2.639        0.299        0.080        0.015
       5     2.401        0.823        0.434        0.177
       2     2.210        1.496        0.740        0.275
       7     2.056        0.452        0.178        0.056
       6     1.932       54.467     -171.523      803.879
       3     1.598       79.891      277.481      884.965
       6     1.481       15.612       28.509       42.556
      78   overall     4773.300     1832.327 32242899.716


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

       9     7.281        0.526        0.361        0.195
      15     4.606        1.197        1.496        1.500
      23     3.510        0.877        0.764        0.629
      23     3.008        0.885        0.762        0.559
      31     2.639        1.257        1.516        1.657
      44     2.401        0.946        0.813        0.592
      36     2.210        0.749        0.530        0.331
      36     2.056        0.974        0.813        0.619
      36     1.932       82.163       60.367     3369.590
      38     1.830    42341.242  2011968.252*************
      63     1.481       58.998       52.069     1556.385
      44     1.392       22.582       38.947      303.071
      17     1.315        4.783        5.049        7.445
     680   overall     2377.827   112444.468 74734803.673

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.281  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.40  0.38  0.73
       3     4.606  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.510  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     3.008  0.00  0.00  0.00  0.00  0.00  0.00  0.56  0.53  0.51  0.49
       2     2.639  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       5     2.401  0.00  0.00  0.00  0.85  0.77  1.07  1.00  0.95  0.91  0.88
       2     2.210  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       7     2.056  0.00  0.00  0.34  0.91  0.82  1.02  0.96  0.91  1.09  1.05
       6     1.932  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       3     1.598  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       6     1.481  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
      78   overall  0.41  0.30  0.28  0.35  0.32  0.34  0.34  0.35  0.35  0.38


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       9     7.281  3.50  3.06  2.14  2.02  1.69  1.72  1.55  1.41  1.31  1.23
      15     4.606  1.40  1.84  1.54  1.21  1.02  1.48  1.32  1.33  1.46  1.37
      23     3.510  0.91  0.72  0.67  0.53  0.66  0.67  0.95  1.03  1.03  1.03
      23     3.008  1.37  2.16  1.51  1.32  1.22  1.25  1.30  1.18  1.17  1.10
      31     2.639  1.02  1.42  1.37  1.37  1.15  1.14  1.22  1.17  1.09  1.07
      44     2.401  0.72  0.75  0.79  1.10  1.16  1.01  1.17  1.16  1.07  1.04
      36     2.210  0.58  0.92  0.86  0.84  0.78  0.80  0.94  0.86  0.89  1.01
      36     2.056  1.17  0.77  0.86  0.67  0.85  0.80  0.77  0.86  0.98  0.97
      36     1.932  4.96  2.61  1.82  1.43  1.20  1.05  0.94  0.86  0.84  0.79
      38     1.830  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      63     1.481  4.84  2.54  1.78  1.40  1.17  1.02  0.91  0.86  0.80  0.75
      44     1.392  5.02  2.63  1.84  1.52  1.27  1.16  1.04  0.99  0.92  0.86
      17     1.315  5.56  2.92  2.04  1.61  1.35  1.17  1.05  0.96  0.89  0.84
     680   overall  1.81  1.15  0.88  0.76  0.68  0.64  0.63  0.61  0.59  0.57


 List of    105 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   31   22    9    1.84  825.19  0.2700E+02  0.1278E+02 "alien"
   24   19   14    1.87  656.83  0.2149E+02  0.1156E+02 "alien"
   28   21   12    1.83  624.76  0.2044E+02  0.6912E+01 "alien"
   35   23    2    1.84  259.55  0.8493E+01  0.7650E+01 "alien"
   31    2   14    1.85  241.74  0.7910E+01  0.6478E+01 "alien"
   34    1   12    1.86  193.74  0.6340E+01  0.6813E+01 "alien"
   34   23    6    1.81  170.66  0.5584E+01  0.7405E+01 "alien"
   26   20   13    1.85  165.85  0.5427E+01  0.6680E+01 "alien"
   33   23    8    1.79  163.37  0.5346E+01  0.7242E+01 "alien"
   43    1    2    1.79  158.53  0.5187E+01  0.7646E+01 "alien"
   25   20   14    1.82  127.49  0.4172E+01  0.6687E+01 "alien"
   37    1   11    1.80  126.10  0.4126E+01  0.7140E+01 "alien"
   43    2   11    1.60  123.88 -0.2152E+02  0.8447E+01 "alien"
   25   24   17    1.60  117.59 -0.2043E+02  0.1006E+02 "alien"
   30   22   11    1.79  116.43  0.3810E+01  0.7137E+01 "alien"
   45    1    8    1.62   98.50 -0.1711E+02  0.8308E+01 "alien"
   39   28    3    1.60   97.58 -0.1695E+02  0.7688E+01 "alien"
   41    3   13    1.59   88.82 -0.1543E+02  0.7549E+01 "alien"
   41    1    5    1.83   80.23  0.2625E+01  0.8123E+01 "alien"
   34   27   11    1.59   76.00 -0.1320E+02  0.8824E+01 "alien"
   46    1    6    1.63   62.67 -0.1089E+02  0.1134E+02 "alien"
   28   24   15    1.63   62.66 -0.1089E+02  0.6888E+01 "alien"
   48    1    3    1.60   55.18 -0.9587E+01  0.9179E+01 "alien"
   44    2   10    1.60   48.27 -0.8388E+01  0.7413E+01 "alien"
   14    8   19    1.88   47.80  0.1564E+01  0.5765E+01 "alien"
   17   15   18    1.82   46.07  0.1507E+01  0.5340E+01 "alien"
   35    1   12    1.82   46.05  0.1507E+01  0.7288E+01 "alien"
   21   17   20    1.61   45.47 -0.7900E+01  0.7918E+01 "alien"
   39   28    2    1.61   45.22 -0.7857E+01  0.7581E+01 "alien"
   22   18   20    1.58   37.25 -0.6471E+01  0.5850E+01 "alien"
   21   15   21    1.57   36.42 -0.6328E+01  0.8793E+01 "alien"
   43    1    0    1.80   35.12  0.1878E+02  0.1089E+02 "alien"
   24   23   17    1.63   34.63 -0.6017E+01  0.7306E+01 "alien"
   38    4   15    1.59   32.14 -0.5583E+01  0.6888E+01 "alien"
   40    0    8    1.80   30.83  0.1093E+02  0.7246E+01 "alien"
   39    0    8    1.84   30.72  0.1078E+02  0.1336E+02 "alien"
   48    2    7    1.55   29.44 -0.2198E+02  0.7990E+01 "alien"
   39   30    6    1.53   28.26 -0.2110E+02  0.1417E+02 "alien"
   34    7   19    1.51   27.37  0.2759E+02  0.1069E+02 "alien"
   36    1    8    1.97   26.97  0.1570E+02  0.7920E+01 "alien"
   48    1    5    1.58   26.13 -0.4540E+01  0.7750E+01 "alien"
   33    6   18    1.57   25.52 -0.4435E+01  0.6200E+01 "alien"
   37   30   11    1.48   25.36  0.2556E+02  0.1014E+02 "alien"
   26   19   12    1.93   25.27  0.1471E+02  0.7063E+01 "alien"
   33   26   12    1.60   25.09 -0.4359E+01  0.8498E+01 "alien"
   33    1   12    1.90   24.21  0.1409E+02  0.7573E+01 "alien"
   47    1    6    1.60   23.72 -0.4121E+01  0.1133E+02 "alien"
   30    8   20    1.53   23.24 -0.1735E+02  0.8762E+01 "alien"
   34   22    4    1.88   22.17  0.7254E+00  0.8709E+01 "alien"
   42    4   14    1.53   21.76 -0.1624E+02  0.7443E+01 "alien"
   35   32   15    1.38   21.30  0.3291E+02  0.1657E+02 "alien"
   37    8   19    1.44   21.02 -0.2210E+02  0.1145E+02 "alien"
   54    3    4    1.42   20.76 -0.2183E+02  0.8627E+01 "alien"
   46    5   14    1.43   20.68 -0.2174E+02  0.9822E+01 "alien"
   33   22    5    1.89   20.65  0.1202E+02  0.1438E+02 "alien"
   36   27    8    1.62   20.50 -0.3562E+01  0.8091E+01 "alien"
   40   32    7    1.46   19.66 -0.2067E+02  0.9342E+01 "alien"
   25   10   22    1.50   18.38  0.1852E+02  0.8251E+01 "alien"
   27   20   12    1.88   17.77  0.5816E+00  0.6729E+01 "alien"
   37    2   13    1.71   17.32 -0.1396E+02  0.7571E+01 "alien"
   52    2    3    1.48   16.96  0.1709E+02  0.8378E+01 "alien"
   51    2    2    1.51   16.76 -0.1251E+02  0.1051E+02 "alien"
   31   27   15    1.52   16.51 -0.1233E+02  0.1367E+02 "alien"
   33   27   13    1.55   16.50 -0.1231E+02  0.1271E+02 "alien"
   30   26   15    1.56   16.41 -0.1225E+02  0.1114E+02 "alien"
   48    2    8    1.53   16.23 -0.1212E+02  0.8331E+01 "alien"
   23   23   18    1.60   15.83 -0.8286E+01  0.9464E+01 "alien"
   32   27   14    1.54   15.82 -0.1181E+02  0.1172E+02 "alien"
   35   31   14    1.42   15.65 -0.1646E+02  0.8266E+01 "alien"
   21   17   14    1.99   14.98  0.8718E+01  0.9432E+01 "alien"
   36   34   15    1.34   14.85 -0.1804E+02  0.1228E+02 "alien"
   35    2   13    1.77   14.62 -0.1178E+02  0.7433E+01 "alien"
   37   31   11    1.46   14.16  0.1427E+02  0.1046E+02 "alien"
   23   15   22    1.49   14.04  0.1415E+02  0.7866E+01 "alien"
   30   23   12    1.73   13.70 -0.1104E+02  0.6995E+01 "alien"
   39   32    9    1.45   13.67 -0.1438E+02  0.9568E+01 "alien"
   52    2    4    1.47   13.52  0.1363E+02  0.1077E+02 "alien"
   34    3   15    1.70   13.18 -0.1062E+02  0.1126E+02 "alien"
   41   31    2    1.50   13.12  0.1322E+02  0.8241E+01 "alien"
   43    1    8    1.69   13.08 -0.2165E+02  0.1402E+02 "alien"
   35   25    7    1.71   12.77 -0.1029E+02  0.7190E+01 "alien"
   27   24   20    1.44   12.28 -0.1291E+02  0.1151E+02 "alien"
   26   24   19    1.49   12.26  0.1236E+02  0.1099E+02 "alien"
   30   12   23    1.38   12.20  0.1886E+02  0.9681E+01 "alien"
   33    8   20    1.48   11.90  0.1200E+02  0.1072E+02 "alien"
   35   24    6    1.76   11.82 -0.9523E+01  0.7721E+01 "alien"
   37   25    2    1.73   11.65 -0.9386E+01  0.8586E+01 "alien"
   52    5   11    1.37   11.58 -0.1407E+02  0.1377E+02 "alien"
   38   29    8    1.54   11.52 -0.8595E+01  0.7952E+01 "alien"
   34   28   12    1.54   11.48 -0.8566E+01  0.9217E+01 "alien"
   28   20   11    1.90   11.38  0.6624E+01  0.7149E+01 "alien"
   29   20   23    1.34   11.24 -0.1365E+02  0.1054E+02 "alien"
   42   32    1    1.47   11.11  0.1119E+02  0.8637E+01 "alien"
   44    1    1    1.76   10.98 -0.8847E+01  0.8152E+01 "alien"
   38   34   13    1.35   10.87 -0.1321E+02  0.8709E+01 "alien"
   50    3    8    1.47   10.75  0.1083E+02  0.8298E+01 "alien"
   44    5   14    1.48   10.72  0.1081E+02  0.1457E+02 "alien"
   51    6   12    1.37   10.64 -0.1292E+02  0.1527E+02 "alien"
   49    7   15    1.34   10.63 -0.1292E+02  0.1317E+02 "alien"
   39    2    1    1.98   10.53  0.6126E+01  0.7986E+01 "alien"
   37    1    7    1.96   10.40  0.6054E+01  0.7543E+01 "alien"
   24   10   22    1.52   10.39 -0.7757E+01  0.7365E+01 "alien"
   30   29   17    1.44   10.27 -0.1080E+02  0.7756E+01 "alien"
   41   32    6    1.45   10.23 -0.1075E+02  0.8318E+01 "alien"
   36   28   10    1.56   10.18 -0.7595E+01  0.7899E+01 "alien"


 List of    105 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   29   20   23    1.34   11.24 -0.1365E+02  0.1054E+02
   36   34   15    1.34   14.85 -0.1804E+02  0.1228E+02
   49    7   15    1.34   10.63 -0.1292E+02  0.1317E+02
   38   34   13    1.35   10.87 -0.1321E+02  0.8709E+01
   52    5   11    1.37   11.58 -0.1407E+02  0.1377E+02
   51    6   12    1.37   10.64 -0.1292E+02  0.1527E+02
   30   12   23    1.38   12.20  0.1886E+02  0.9681E+01
   35   32   15    1.38   21.30  0.3291E+02  0.1657E+02
   54    3    4    1.42   20.76 -0.2183E+02  0.8627E+01
   35   31   14    1.42   15.65 -0.1646E+02  0.8266E+01
   46    5   14    1.43   20.68 -0.2174E+02  0.9822E+01
   27   24   20    1.44   12.28 -0.1291E+02  0.1151E+02
   30   29   17    1.44   10.27 -0.1080E+02  0.7756E+01
   37    8   19    1.44   21.02 -0.2210E+02  0.1145E+02
   39   32    9    1.45   13.67 -0.1438E+02  0.9568E+01
   41   32    6    1.45   10.23 -0.1075E+02  0.8318E+01
   40   32    7    1.46   19.66 -0.2067E+02  0.9342E+01
   37   31   11    1.46   14.16  0.1427E+02  0.1046E+02
   42   32    1    1.47   11.11  0.1119E+02  0.8637E+01
   52    2    4    1.47   13.52  0.1363E+02  0.1077E+02
   50    3    8    1.47   10.75  0.1083E+02  0.8298E+01
   37   30   11    1.48   25.36  0.2556E+02  0.1014E+02
   33    8   20    1.48   11.90  0.1200E+02  0.1072E+02
   52    2    3    1.48   16.96  0.1709E+02  0.8378E+01
   44    5   14    1.48   10.72  0.1081E+02  0.1457E+02
   26   24   19    1.49   12.26  0.1236E+02  0.1099E+02
   23   15   22    1.49   14.04  0.1415E+02  0.7866E+01
   25   10   22    1.50   18.38  0.1852E+02  0.8251E+01
   41   31    2    1.50   13.12  0.1322E+02  0.8241E+01
   34    7   19    1.51   27.37  0.2759E+02  0.1069E+02
   51    2    2    1.51   16.76 -0.1251E+02  0.1051E+02
   24   10   22    1.52   10.39 -0.7757E+01  0.7365E+01
   31   27   15    1.52   16.51 -0.1233E+02  0.1367E+02
   39   30    6    1.53   28.26 -0.2110E+02  0.1417E+02
   48    2    8    1.53   16.23 -0.1212E+02  0.8331E+01
   42    4   14    1.53   21.76 -0.1624E+02  0.7443E+01
   30    8   20    1.53   23.24 -0.1735E+02  0.8762E+01
   38   29    8    1.54   11.52 -0.8595E+01  0.7952E+01
   32   27   14    1.54   15.82 -0.1181E+02  0.1172E+02
   34   28   12    1.54   11.48 -0.8566E+01  0.9217E+01
   48    2    7    1.55   29.44 -0.2198E+02  0.7990E+01
   33   27   13    1.55   16.50 -0.1231E+02  0.1271E+02
   36   28   10    1.56   10.18 -0.7595E+01  0.7899E+01
   30   26   15    1.56   16.41 -0.1225E+02  0.1114E+02
   21   15   21    1.57   36.42 -0.6328E+01  0.8793E+01
   33    6   18    1.57   25.52 -0.4435E+01  0.6200E+01
   22   18   20    1.58   37.25 -0.6471E+01  0.5850E+01
   48    1    5    1.58   26.13 -0.4540E+01  0.7750E+01
   34   27   11    1.59   76.00 -0.1320E+02  0.8824E+01
   41    3   13    1.59   88.82 -0.1543E+02  0.7549E+01
   38    4   15    1.59   32.14 -0.5583E+01  0.6888E+01
   25   24   17    1.60  117.59 -0.2043E+02  0.1006E+02
   47    1    6    1.60   23.72 -0.4121E+01  0.1133E+02
   23   23   18    1.60   15.83 -0.8286E+01  0.9464E+01
   48    1    3    1.60   55.18 -0.9587E+01  0.9179E+01
   39   28    3    1.60   97.58 -0.1695E+02  0.7688E+01
   33   26   12    1.60   25.09 -0.4359E+01  0.8498E+01
   44    2   10    1.60   48.27 -0.8388E+01  0.7413E+01
   43    2   11    1.60  123.88 -0.2152E+02  0.8447E+01
   21   17   20    1.61   45.47 -0.7900E+01  0.7918E+01
   39   28    2    1.61   45.22 -0.7857E+01  0.7581E+01
   45    1    8    1.62   98.50 -0.1711E+02  0.8308E+01
   36   27    8    1.62   20.50 -0.3562E+01  0.8091E+01
   46    1    6    1.63   62.67 -0.1089E+02  0.1134E+02
   24   23   17    1.63   34.63 -0.6017E+01  0.7306E+01
   28   24   15    1.63   62.66 -0.1089E+02  0.6888E+01
   43    1    8    1.69   13.08 -0.2165E+02  0.1402E+02
   34    3   15    1.70   13.18 -0.1062E+02  0.1126E+02
   37    2   13    1.71   17.32 -0.1396E+02  0.7571E+01
   35   25    7    1.71   12.77 -0.1029E+02  0.7190E+01
   37   25    2    1.73   11.65 -0.9386E+01  0.8586E+01
   30   23   12    1.73   13.70 -0.1104E+02  0.6995E+01
   35   24    6    1.76   11.82 -0.9523E+01  0.7721E+01
   44    1    1    1.76   10.98 -0.8847E+01  0.8152E+01
   35    2   13    1.77   14.62 -0.1178E+02  0.7433E+01
   33   23    8    1.79  163.37  0.5346E+01  0.7242E+01
   30   22   11    1.79  116.43  0.3810E+01  0.7137E+01
   43    1    2    1.79  158.53  0.5187E+01  0.7646E+01
   40    0    8    1.80   30.83  0.1093E+02  0.7246E+01
   37    1   11    1.80  126.10  0.4126E+01  0.7140E+01
   43    1    0    1.80   35.12  0.1878E+02  0.1089E+02
   34   23    6    1.81  170.66  0.5584E+01  0.7405E+01
   17   15   18    1.82   46.07  0.1507E+01  0.5340E+01
   25   20   14    1.82  127.49  0.4172E+01  0.6687E+01
   35    1   12    1.82   46.05  0.1507E+01  0.7288E+01
   28   21   12    1.83  624.76  0.2044E+02  0.6912E+01
   41    1    5    1.83   80.23  0.2625E+01  0.8123E+01
   39    0    8    1.84   30.72  0.1078E+02  0.1336E+02
   35   23    2    1.84  259.55  0.8493E+01  0.7650E+01
   31   22    9    1.84  825.19  0.2700E+02  0.1278E+02
   26   20   13    1.85  165.85  0.5427E+01  0.6680E+01
   31    2   14    1.85  241.74  0.7910E+01  0.6478E+01
   34    1   12    1.86  193.74  0.6340E+01  0.6813E+01
   24   19   14    1.87  656.83  0.2149E+02  0.1156E+02
   27   20   12    1.88   17.77  0.5816E+00  0.6729E+01
   34   22    4    1.88   22.17  0.7254E+00  0.8709E+01
   14    8   19    1.88   47.80  0.1564E+01  0.5765E+01
   33   22    5    1.89   20.65  0.1202E+02  0.1438E+02
   28   20   11    1.90   11.38  0.6624E+01  0.7149E+01
   33    1   12    1.90   24.21  0.1409E+02  0.7573E+01
   26   19   12    1.93   25.27  0.1471E+02  0.7063E+01
   37    1    7    1.96   10.40  0.6054E+01  0.7543E+01
   36    1    8    1.97   26.97  0.1570E+02  0.7920E+01
   39    2    1    1.98   10.53  0.6126E+01  0.7986E+01
   21   17   14    1.99   14.98  0.8718E+01  0.9432E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        56
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1349
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1197
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      152


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
