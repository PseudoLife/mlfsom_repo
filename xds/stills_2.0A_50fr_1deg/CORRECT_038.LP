
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.588    77.588    38.830  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      38      38
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      38
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000981  0.000646
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
 BEAM_DIVERGENCE_E.S.D.=   0.048
 REFLECTING_RANGE_E.S.D.=   0.178
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

    1180 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      38 ...      38
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1180 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               39
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

 *  44        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.6   38.8   77.6  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  33        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  32        oP          0.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.7  109.7   38.8  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.7  109.7   38.8  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.6   77.6   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     160.0   38.8   77.6  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     160.0   38.8   77.6  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.8  160.0   77.6  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.8  160.0   77.6  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.8  160.0   77.6  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.8  160.0   77.6  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     160.0   38.8  109.7  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.8  109.7  116.4  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.8  160.0  160.0  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.7  116.4   38.8  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    30        mC        998.1      77.6  173.5   38.8  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.1      77.6  173.5   38.8  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    41        mC        998.1     173.5   77.6   38.8  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    22        hP        998.1      77.6   77.6   38.8  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.4  116.4  116.4  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.8   86.8  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.8   86.8  116.4 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.7   86.8   86.8  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.8   86.8  109.7  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.8   86.8  109.7  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.8   86.8  248.4 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.8   86.8   77.6  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0 -1  0  0  0  0 -1  0  0
    12        hP        999.0      38.8   77.6   77.6  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.8   77.6  177.8  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.8   86.8  177.8 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.8   86.8   77.6  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.8   86.8   86.8  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.8  109.7   38.8  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.8   77.6  177.8 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.8  177.8   77.6 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    4.7      15    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.948   -23.486     0.049
 UNIT_CELL_B-AXIS=    23.486    73.948     0.014
 UNIT_CELL_C-AXIS=    -0.025     0.001    38.830
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       132 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.36
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.12
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.634    77.634    38.837  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  2.1E-01 2.1E-01 7.7E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012881  0.012881  0.025749  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    74.012   -23.437     0.084
 COORDINATES OF UNIT CELL B-AXIS    23.437    74.012     0.024
 COORDINATES OF UNIT CELL C-AXIS    -0.044     0.001    38.837
 CRYSTAL MOSAICITY (DEGREES)     0.178
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999970  0.007462  0.002124
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.88   1232.95
 DETECTOR ORIGIN (PIXELS) AT                     1260.78   1243.82
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.48
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.63    77.63    38.84  90.000  90.000  90.000


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
  0.05  0.15        26         55.055       1.909
  0.15  0.25       149        102.305       1.027
  0.25  0.35       165        118.481       0.887
  0.35  0.45       159         42.712       2.461
  0.45  0.55       164        183.636       0.572
  0.55  0.65       175         67.437       1.558
  0.65  0.75       193        127.034       0.827
  0.75  0.85       136         99.656       1.055
  0.85  0.95        13         78.022       1.347
 
 -0.05  0.95      1180        105.100       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        7
 DEGREES OF FREEDOM OF CHI^2 FIT               172.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.804
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    37.1 XMAX=    37.9 NXBIN=    1
 YMIN= 0.00166 YMAX= 0.59192 NYBIN=    7
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        358


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               172.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.794
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    67.0 XMAX=  2451.6 NXBIN=    3
 YMIN=     5.8 YMAX=  2518.5 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        358


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               172.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.787
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    37.1 XMAX=    37.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1259    1262
 DETECTOR_SURFACE_POSITION=    1664    1689
 DETECTOR_SURFACE_POSITION=     854    1689
 DETECTOR_SURFACE_POSITION=     854     835
 DETECTOR_SURFACE_POSITION=    1664     835
 DETECTOR_SURFACE_POSITION=    2176    1662
 DETECTOR_SURFACE_POSITION=    1639    2228
 DETECTOR_SURFACE_POSITION=     880    2228
 DETECTOR_SURFACE_POSITION=     343    1662
 DETECTOR_SURFACE_POSITION=     343     862
 DETECTOR_SURFACE_POSITION=     880     296
 DETECTOR_SURFACE_POSITION=    1639     296
 DETECTOR_SURFACE_POSITION=    2176     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        358


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
 8.850E-01  1.846E-03   24.74


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

   24.550   2.890    15.59   1.52      2.88      2.65      88     114      14
    2.890   2.051     3.38   0.86     13.55     14.24     210     236       0
    2.051   1.676     0.06   0.87    864.07    883.98     242     277       0
    1.676   1.453     0.00   0.92   2432.50   2389.66     258     313       2
    1.453   1.300     0.01   0.73   1663.70   1893.65     172     217       6
   --------------------------------------------------------------------------
   24.550   1.300     2.24   0.91      9.89      9.98     970    1157      22
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      38      38)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.688  -0.7482E+01  0.1297E+02        -0.58           1 
   47    0    0     1.652   0.1221E+02  0.1286E+02         0.95           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    10.7% OF MEAN INTENSITY


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

     3.86          59      37      1285        2.9%       2.6%      2.4%       44   24.41      3.6%    99.9*     0    0.000       0
     2.74          77      45      2097        2.1%       5.7%      4.6%       64   13.57      8.0%    99.1*     0    0.000       0
     2.24         143      80      2667        3.0%      11.0%     13.5%      126    5.11     15.6%    97.0*     0    0.000       0
     1.95         129      72      3098        2.3%      76.2%     71.1%      114    0.92    107.7%    47.7*     0    0.000       0
     1.74         144      82      3481        2.4%    1204.5%   1392.1%      124    0.08   1703.3%     7.5      0    0.000       0
     1.59         178     102      3853        2.6%     -99.9%    -99.9%      152    0.00    -99.9%   -19.4      0    0.000       0
     1.47         173     102      4150        2.5%    6975.1%   6676.3%      142   -0.02   9864.2%     1.5      0    0.000       0
     1.38         174     100      4423        2.3%   16776.8%  19683.2%      148   -0.01  23726.0%     3.2      0    0.000       0
     1.30          80      52      4735        1.1%     375.2%    388.7%       56    0.07    530.6%    14.4      0    0.000       0
    total        1157     672     29789        2.3%       9.9%     10.0%      970    2.97     14.0%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1180
 NUMBER OF REJECTED MISFITS                                22
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1157
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    672



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

   WILSON LINE (using all data) : A=   6.684 B=  32.245 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      16     7.335   0.005  1.7132E+03       7.446     -81.9
      17     4.744   0.011  1.6909E+03       7.433     -33.7
      19     3.587   0.019  5.8752E+02       6.376       7.9
      25     2.990   0.028  1.3135E+02       4.878      32.3
      31     2.641   0.036  1.1329E+02       4.730      27.3
      42     2.410   0.043  4.0215E+01       3.694      34.7
      32     2.215   0.051  1.9618E+01       2.976      36.4
      35     2.062   0.059  6.8072E+00       1.918      40.5
      37     1.939   0.066  3.3170E-01      -1.104      58.6
      44     1.832   0.074  2.9058E+00       1.067      37.7
      40     1.744   0.082 -3.5431E-01    -999.999   -1000.0
      43     1.672   0.089 -1.0123E+00    -999.999   -1000.0
      49     1.602   0.097 -2.4887E-01    -999.999   -1000.0
      48     1.541   0.105  1.3958E-01      -1.969      41.1
      58     1.487   0.113  7.6978E-01      -0.262      30.7
      45     1.437   0.121 -1.2741E+00    -999.999   -1000.0
      40     1.396   0.128  4.9764E-01      -0.698      28.8
      35     1.355   0.136  1.8734E+00       0.628      22.2
      16     1.318   0.144 -1.7942E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       4     7.335        2.322        1.862        1.167
       3     4.744        1.317        0.825        0.384
       2     3.587        0.402        0.089        0.014
       3     2.990        2.717        1.632        0.720
       1     2.641        0.000        0.000        0.000
       3     2.410        0.095        0.014        0.002
       1     2.215        0.642        0.178        0.035
       6     2.062        1.084        0.264        0.168
       6     1.939      254.401     1320.558    10264.947
       9     1.832        4.228        2.561        2.990
       7     1.744       55.698      144.643      633.284
       3     1.602       90.860      443.181     1624.721
       3     1.541      680.698     5993.636    52983.242
       3     1.437        5.487        3.794        4.678
       1     1.396     7409.818   220954.037  4706177.480
       1     1.355       19.784       30.482       33.548
      72   overall      163.314     3366.408    68795.083


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.335        0.580        0.474        0.307
      14     4.744        1.210        1.474        1.431
      17     3.587        1.015        1.075        1.038
      22     2.990        0.667        0.542        0.377
      30     2.641        1.260        1.648        2.031
      39     2.410        1.070        1.025        0.883
      31     2.215        1.045        0.863        0.631
      29     2.062        1.637        1.836        1.984
      35     1.832        5.589        5.367       13.738
      53     1.487       96.705      295.216     4079.751
      34     1.355       24.761       20.772      263.395
     600   overall      233.597     -989.248   253488.954

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       4     7.335  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.38  0.37
       3     4.744  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       2     3.587  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       3     2.990  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     2.641  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       3     2.410  0.00  0.96  0.80  0.70  1.28  1.19  1.12  1.06  1.52  1.47
       1     2.215  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     2.062  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  1.01  0.98
       6     1.939  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       9     1.832  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       7     1.744  3.45  2.48  2.06  1.81  1.65  1.53  1.44  1.36  1.30  1.26
       3     1.602  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       3     1.541  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       3     1.437  1.34  0.96  0.80  0.70  1.28  1.19  1.12  1.06  1.01  0.98
       1     1.396  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     1.355  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      72   overall  1.17  0.88  0.73  0.65  0.64  0.59  0.56  0.53  0.57  0.57


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.335  1.75  3.22  2.25  1.77  1.91  1.66  1.66  1.51  1.40  1.32
      14     4.744  2.25  1.18  1.38  1.08  1.09  1.42  1.42  1.30  1.44  1.36
      17     3.587  1.24  0.65  0.68  0.71  0.75  0.78  1.29  1.28  1.19  1.12
      22     2.990  0.96  2.01  1.40  1.38  1.39  1.41  1.35  1.32  1.30  1.22
      30     2.641  0.70  1.10  1.03  1.11  0.93  1.03  0.99  1.03  1.01  1.16
      39     2.410  1.08  0.85  0.79  0.78  0.98  0.97  0.87  0.88  1.04  1.01
      31     2.215  1.36  1.07  1.00  1.17  0.98  0.93  0.90  0.94  0.87  0.87
      29     2.062  2.17  1.33  1.20  0.94  0.88  0.92  0.82  0.81  0.87  0.82
      35     1.832  4.50  2.36  1.76  1.39  1.23  1.20  1.08  0.99  0.91  0.86
      53     1.487  4.96  2.60  1.82  1.43  1.20  1.05  0.94  0.86  0.83  0.78
      34     1.355  4.33  2.27  1.59  1.34  1.12  0.98  0.93  0.91  0.84  0.79
     600   overall  1.38  0.91  0.71  0.63  0.58  0.57  0.54  0.53  0.52  0.51


 List of    103 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   30   26   15    1.56  151.55  0.2115E+02  0.1322E+02 "alien"
   38   28    5    1.61  142.10 -0.3536E+02  0.1601E+02 "alien"
   27   26   17    1.53  130.88  0.1827E+02  0.1452E+02 "alien"
   40   30    3    1.54  127.81  0.1784E+02  0.1117E+02 "alien"
   48    2    7    1.55  124.06  0.1732E+02  0.1029E+02 "alien"
   33   27   13    1.55   90.58  0.1264E+02  0.1480E+02 "alien"
   35    5   17    1.58   85.99 -0.2140E+02  0.1443E+02 "alien"
   37   29   10    1.52   68.03  0.9496E+01  0.1587E+02 "alien"
   22   18   20    1.58   65.62 -0.1633E+02  0.1027E+02 "alien"
   24   10   22    1.52   62.27  0.8692E+01  0.1052E+02 "alien"
   39   28    1    1.62   61.64 -0.1534E+02  0.9269E+01 "alien"
   27    9   21    1.53   61.46  0.8579E+01  0.1128E+02 "alien"
   35   25    7    1.72   61.08 -0.2164E+02  0.9026E+01 "alien"
   21   17   20    1.61   59.76 -0.1487E+02  0.1036E+02 "alien"
   50    1    1    1.55   59.27  0.8274E+01  0.9716E+01 "alien"
   35   24    6    1.76   57.89 -0.2051E+02  0.9013E+01 "alien"
   39    1   11    1.73   56.49 -0.2001E+02  0.1606E+02 "alien"
   39    4   15    1.57   55.13 -0.1372E+02  0.1352E+02 "alien"
   45    1    8    1.63   49.97 -0.1244E+02  0.1040E+02 "alien"
   47    1    5    1.62   48.57 -0.1209E+02  0.9516E+01 "alien"
   40    2    1    1.94   44.83  0.1487E+02  0.1187E+02 "alien"
   29    1   13    1.99   43.67  0.1449E+02  0.8865E+01 "alien"
   45    2   10    1.58   43.59 -0.1085E+02  0.1052E+02 "alien"
   30    8   20    1.53   43.55  0.6079E+01  0.1182E+02 "alien"
   31   25   13    1.63   42.57 -0.1059E+02  0.8202E+01 "alien"
   48    1    3    1.60   41.00 -0.1020E+02  0.1163E+02 "alien"
   39   29    5    1.56   40.77  0.5690E+01  0.9036E+01 "alien"
   33   11   22    1.38   38.32  0.1907E+02  0.1358E+02 "alien"
   40   30    5    1.52   37.68  0.5259E+01  0.9529E+01 "alien"
   37   29    9    1.54   37.21  0.5194E+01  0.1016E+02 "alien"
   35   28   11    1.56   36.86  0.5145E+01  0.9381E+01 "alien"
   50    2    5    1.52   33.64  0.4695E+01  0.9539E+01 "alien"
   52    2    3    1.48   31.58  0.2431E+02  0.1007E+02 "alien"
   42    8   18    1.39   31.57  0.1571E+02  0.1049E+02 "alien"
   31    9   21    1.47   30.22  0.2326E+02  0.1425E+02 "alien"
   30    1   13    1.96   29.94  0.9931E+01  0.1445E+02 "alien"
   34    6   18    1.56   28.53  0.3982E+01  0.8296E+01 "alien"
   25   18   12    1.99   27.68  0.9181E+01  0.1312E+02 "alien"
   39   30    7    1.52   27.12  0.3785E+01  0.9259E+01 "alien"
   40   29    1    1.57   26.93  0.3758E+01  0.9341E+01 "alien"
   37   32   12    1.42   26.28 -0.3348E+02  0.1686E+02 "alien"
   31   30   18    1.38   26.07  0.1297E+02  0.1628E+02 "alien"
   51    2    5    1.49   25.82  0.1988E+02  0.1050E+02 "alien"
   48    3    9    1.51   25.21  0.1940E+02  0.1706E+02 "alien"
   29   28   17    1.47   24.79  0.1908E+02  0.9384E+01 "alien"
   25   21   15    1.75   24.59 -0.8713E+01  0.8918E+01 "alien"
   55    3    0    1.41   24.48  0.7420E+02  0.1503E+02 "alien"
   38   29    7    1.56   23.89  0.3334E+01  0.1116E+02 "alien"
   28   24   15    1.63   23.48 -0.5844E+01  0.8296E+01 "alien"
   39   33   11    1.40   23.33  0.1161E+02  0.1094E+02 "alien"
   34   29   13    1.50   22.92  0.1764E+02  0.1146E+02 "alien"
   47    4   12    1.47   22.82  0.1757E+02  0.1006E+02 "alien"
   26   19   12    1.93   22.25  0.7381E+01  0.8271E+01 "alien"
   52    4    9    1.41   22.08  0.1099E+02  0.1067E+02 "alien"
   34   24    7    1.77   21.91 -0.7764E+01  0.1503E+02 "alien"
   44    4   13    1.51   20.99  0.2930E+01  0.9603E+01 "alien"
   30   28   19    1.39   20.19  0.1005E+02  0.9265E+01 "alien"
   42   32    1    1.47   19.95  0.1536E+02  0.9620E+01 "alien"
   27   24   16    1.61   19.94 -0.4962E+01  0.1428E+02 "alien"
   27   22   14    1.74   19.47 -0.6898E+01  0.7735E+01 "alien"
   34   22    1    1.91   19.44  0.6448E+01  0.9893E+01 "alien"
   30   25   14    1.62   18.60 -0.4628E+01  0.8484E+01 "alien"
   46    1    6    1.63   18.13 -0.4512E+01  0.1337E+02 "alien"
   33    3   15    1.74   17.83 -0.6318E+01  0.7954E+01 "alien"
   38    1   11    1.77   17.68 -0.6265E+01  0.8536E+01 "alien"
   25   25   18    1.54   17.53  0.8384E+01  0.8133E+01 "alien"
   54    4    6    1.40   16.95  0.8435E+01  0.1076E+02 "alien"
   37   10   21    1.37   16.80  0.3148E+02  0.1527E+02 "alien"
   37   30   11    1.48   16.59  0.1277E+02  0.1136E+02 "alien"
   36   10   21    1.38   16.58  0.8249E+01  0.9187E+01 "alien"
   32   26   13    1.59   15.67 -0.3899E+01  0.9087E+01 "alien"
   37    0    9    1.89   15.64  0.1541E+02  0.1492E+02 "alien"
   33    7   19    1.53   15.27  0.2132E+01  0.8064E+01 "alien"
   46    4   12    1.49   15.08  0.1161E+02  0.1088E+02 "alien"
   45    2    9    1.60   14.96 -0.3722E+01  0.1546E+02 "alien"
   32   21    6    1.94   14.94  0.4957E+01  0.8732E+01 "alien"
   50    1    0    1.55   14.91  0.5857E+01  0.1381E+02 "alien"
   42    0    7    1.75   14.90 -0.1484E+02  0.8364E+01 "alien"
   38    2    4    2.00   14.59  0.4838E+01  0.1177E+02 "alien"
   46    7   16    1.37   14.20  0.7067E+01  0.9868E+01 "alien"
   25   10   22    1.51   13.99  0.1077E+02  0.1090E+02 "alien"
   37    1    7    1.96   13.82  0.4585E+01  0.8993E+01 "alien"
   30   20    8    1.97   13.70  0.4546E+01  0.8587E+01 "alien"
   41   11   21    1.30   13.58 -0.2437E+02  0.1601E+02 "alien"
   28   22   13    1.76   13.12 -0.4648E+01  0.7732E+01 "alien"
   22   14   22    1.52   12.93  0.1805E+01  0.1011E+02 "alien"
   32   28   20    1.33   12.85 -0.2305E+02  0.1653E+02 "alien"
   44    2   10    1.61   12.41 -0.3089E+01  0.8841E+01 "alien"
   34    0   10    1.97   12.06  0.9225E+01  0.1389E+02 "alien"
   54    5    9    1.36   12.02  0.2252E+02  0.1107E+02 "alien"
   40    2    0    1.94   11.94  0.9059E+01  0.1509E+02 "alien"
   47    2    8    1.56   11.41  0.1592E+01  0.9519E+01 "alien"
   19   14   20    1.67   11.21 -0.1135E+02  0.8848E+01 "alien"
   33   25   11    1.66   11.20 -0.1134E+02  0.1408E+02 "alien"
   36   27    8    1.63   11.19 -0.2786E+01  0.1037E+02 "alien"
   39    2    2    1.98   11.02  0.3654E+01  0.9399E+01 "alien"
   50    5   12    1.39   10.59  0.5272E+01  0.1091E+02 "alien"
   40   31    7    1.48   10.58  0.8141E+01  0.1054E+02 "alien"
   53    3    5    1.44   10.34 -0.1317E+02  0.9619E+01 "alien"
   51    6   13    1.35   10.23  0.1916E+02  0.1883E+02 "alien"
   43    1    3    1.79   10.16 -0.3601E+01  0.9467E+01 "alien"
   48    7   15    1.36   10.16  0.1904E+02  0.1101E+02 "alien"
   19   18   18    1.74   10.09 -0.3577E+01  0.6779E+01 "alien"


 List of    103 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   41   11   21    1.30   13.58 -0.2437E+02  0.1601E+02
   32   28   20    1.33   12.85 -0.2305E+02  0.1653E+02
   51    6   13    1.35   10.23  0.1916E+02  0.1883E+02
   54    5    9    1.36   12.02  0.2252E+02  0.1107E+02
   48    7   15    1.36   10.16  0.1904E+02  0.1101E+02
   37   10   21    1.37   16.80  0.3148E+02  0.1527E+02
   46    7   16    1.37   14.20  0.7067E+01  0.9868E+01
   36   10   21    1.38   16.58  0.8249E+01  0.9187E+01
   31   30   18    1.38   26.07  0.1297E+02  0.1628E+02
   33   11   22    1.38   38.32  0.1907E+02  0.1358E+02
   30   28   19    1.39   20.19  0.1005E+02  0.9265E+01
   42    8   18    1.39   31.57  0.1571E+02  0.1049E+02
   50    5   12    1.39   10.59  0.5272E+01  0.1091E+02
   39   33   11    1.40   23.33  0.1161E+02  0.1094E+02
   54    4    6    1.40   16.95  0.8435E+01  0.1076E+02
   52    4    9    1.41   22.08  0.1099E+02  0.1067E+02
   55    3    0    1.41   24.48  0.7420E+02  0.1503E+02
   37   32   12    1.42   26.28 -0.3348E+02  0.1686E+02
   53    3    5    1.44   10.34 -0.1317E+02  0.9619E+01
   31    9   21    1.47   30.22  0.2326E+02  0.1425E+02
   47    4   12    1.47   22.82  0.1757E+02  0.1006E+02
   42   32    1    1.47   19.95  0.1536E+02  0.9620E+01
   29   28   17    1.47   24.79  0.1908E+02  0.9384E+01
   40   31    7    1.48   10.58  0.8141E+01  0.1054E+02
   37   30   11    1.48   16.59  0.1277E+02  0.1136E+02
   52    2    3    1.48   31.58  0.2431E+02  0.1007E+02
   46    4   12    1.49   15.08  0.1161E+02  0.1088E+02
   51    2    5    1.49   25.82  0.1988E+02  0.1050E+02
   34   29   13    1.50   22.92  0.1764E+02  0.1146E+02
   25   10   22    1.51   13.99  0.1077E+02  0.1090E+02
   48    3    9    1.51   25.21  0.1940E+02  0.1706E+02
   44    4   13    1.51   20.99  0.2930E+01  0.9603E+01
   39   30    7    1.52   27.12  0.3785E+01  0.9259E+01
   22   14   22    1.52   12.93  0.1805E+01  0.1011E+02
   24   10   22    1.52   62.27  0.8692E+01  0.1052E+02
   37   29   10    1.52   68.03  0.9496E+01  0.1587E+02
   50    2    5    1.52   33.64  0.4695E+01  0.9539E+01
   40   30    5    1.52   37.68  0.5259E+01  0.9529E+01
   33    7   19    1.53   15.27  0.2132E+01  0.8064E+01
   27    9   21    1.53   61.46  0.8579E+01  0.1128E+02
   30    8   20    1.53   43.55  0.6079E+01  0.1182E+02
   27   26   17    1.53  130.88  0.1827E+02  0.1452E+02
   25   25   18    1.54   17.53  0.8384E+01  0.8133E+01
   40   30    3    1.54  127.81  0.1784E+02  0.1117E+02
   37   29    9    1.54   37.21  0.5194E+01  0.1016E+02
   50    1    1    1.55   59.27  0.8274E+01  0.9716E+01
   48    2    7    1.55  124.06  0.1732E+02  0.1029E+02
   50    1    0    1.55   14.91  0.5857E+01  0.1381E+02
   33   27   13    1.55   90.58  0.1264E+02  0.1480E+02
   35   28   11    1.56   36.86  0.5145E+01  0.9381E+01
   34    6   18    1.56   28.53  0.3982E+01  0.8296E+01
   38   29    7    1.56   23.89  0.3334E+01  0.1116E+02
   30   26   15    1.56  151.55  0.2115E+02  0.1322E+02
   47    2    8    1.56   11.41  0.1592E+01  0.9519E+01
   39   29    5    1.56   40.77  0.5690E+01  0.9036E+01
   40   29    1    1.57   26.93  0.3758E+01  0.9341E+01
   39    4   15    1.57   55.13 -0.1372E+02  0.1352E+02
   45    2   10    1.58   43.59 -0.1085E+02  0.1052E+02
   22   18   20    1.58   65.62 -0.1633E+02  0.1027E+02
   35    5   17    1.58   85.99 -0.2140E+02  0.1443E+02
   32   26   13    1.59   15.67 -0.3899E+01  0.9087E+01
   45    2    9    1.60   14.96 -0.3722E+01  0.1546E+02
   48    1    3    1.60   41.00 -0.1020E+02  0.1163E+02
   44    2   10    1.61   12.41 -0.3089E+01  0.8841E+01
   21   17   20    1.61   59.76 -0.1487E+02  0.1036E+02
   27   24   16    1.61   19.94 -0.4962E+01  0.1428E+02
   38   28    5    1.61  142.10 -0.3536E+02  0.1601E+02
   47    1    5    1.62   48.57 -0.1209E+02  0.9516E+01
   39   28    1    1.62   61.64 -0.1534E+02  0.9269E+01
   30   25   14    1.62   18.60 -0.4628E+01  0.8484E+01
   45    1    8    1.63   49.97 -0.1244E+02  0.1040E+02
   36   27    8    1.63   11.19 -0.2786E+01  0.1037E+02
   31   25   13    1.63   42.57 -0.1059E+02  0.8202E+01
   46    1    6    1.63   18.13 -0.4512E+01  0.1337E+02
   28   24   15    1.63   23.48 -0.5844E+01  0.8296E+01
   33   25   11    1.66   11.20 -0.1134E+02  0.1408E+02
   19   14   20    1.67   11.21 -0.1135E+02  0.8848E+01
   35   25    7    1.72   61.08 -0.2164E+02  0.9026E+01
   39    1   11    1.73   56.49 -0.2001E+02  0.1606E+02
   33    3   15    1.74   17.83 -0.6318E+01  0.7954E+01
   27   22   14    1.74   19.47 -0.6898E+01  0.7735E+01
   19   18   18    1.74   10.09 -0.3577E+01  0.6779E+01
   25   21   15    1.75   24.59 -0.8713E+01  0.8918E+01
   42    0    7    1.75   14.90 -0.1484E+02  0.8364E+01
   35   24    6    1.76   57.89 -0.2051E+02  0.9013E+01
   28   22   13    1.76   13.12 -0.4648E+01  0.7732E+01
   38    1   11    1.77   17.68 -0.6265E+01  0.8536E+01
   34   24    7    1.77   21.91 -0.7764E+01  0.1503E+02
   43    1    3    1.79   10.16 -0.3601E+01  0.9467E+01
   37    0    9    1.89   15.64  0.1541E+02  0.1492E+02
   34   22    1    1.91   19.44  0.6448E+01  0.9893E+01
   26   19   12    1.93   22.25  0.7381E+01  0.8271E+01
   32   21    6    1.94   14.94  0.4957E+01  0.8732E+01
   40    2    1    1.94   44.83  0.1487E+02  0.1187E+02
   40    2    0    1.94   11.94  0.9059E+01  0.1509E+02
   30    1   13    1.96   29.94  0.9931E+01  0.1445E+02
   37    1    7    1.96   13.82  0.4585E+01  0.8993E+01
   30   20    8    1.97   13.70  0.4546E+01  0.8587E+01
   34    0   10    1.97   12.06  0.9225E+01  0.1389E+02
   39    2    2    1.98   11.02  0.3654E+01  0.9399E+01
   25   18   12    1.99   27.68  0.9181E+01  0.1312E+02
   29    1   13    1.99   43.67  0.1449E+02  0.8865E+01
   38    2    4    2.00   14.59  0.4838E+01  0.1177E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        57
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1180
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1066
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      114


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
