
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.350    77.350    38.712  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      10      10
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      10
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  1.000000 -0.000741  0.000692
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
 REFLECTING_RANGE_E.S.D.=   0.088
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

    1488 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      10 ...      10
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1488 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               47
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

 *  44        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0 -1  0  0  0  0  1  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.5   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.5   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.7  159.5   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.7  159.5   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.7  159.5   77.3  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.7  159.5   77.3  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     159.5   38.7  109.4  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.7  109.4  116.0  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.7  159.5  159.5  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.4  116.0   38.7  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    41        mC        998.1     173.0   77.3   38.7  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    30        mC        998.1      77.3  173.0   38.7  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.1      77.3  173.0   38.7  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    22        hP        998.1      77.3   77.3   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    1  0 -1  0  0  1  1  0  1 -1  1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.4   86.5   86.5  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.5  247.6 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.5   86.5   77.3  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.7   77.3  177.2  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.5   86.5  177.2 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.2  109.4   38.7  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.7   77.3  177.2 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.5  177.2   77.3 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    6.0      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.719   -23.419     0.030
 UNIT_CELL_B-AXIS=    23.419    73.719     0.012
 UNIT_CELL_C-AXIS=    -0.016    -0.001    38.712
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       322 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.20
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.326    77.326    38.699  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  6.7E-02 6.7E-02 3.5E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012932  0.012932  0.025841  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.690   -23.432     0.025
 COORDINATES OF UNIT CELL B-AXIS    23.432    73.690     0.007
 COORDINATES OF UNIT CELL C-AXIS    -0.013     0.000    38.699
 CRYSTAL MOSAICITY (DEGREES)     0.088
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006295  0.001890
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.87
 DETECTOR ORIGIN (PIXELS) AT                     1261.12   1242.03
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.23
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.33    77.33    38.70  90.000  90.000  90.000


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
  0.05  0.15       120        249.500       0.990
  0.15  0.25       170        227.027       1.088
  0.25  0.35       195        193.173       1.279
  0.35  0.45       146        171.562       1.440
  0.45  0.55       172        324.782       0.761
  0.55  0.65       204        230.243       1.073
  0.65  0.75       182        169.698       1.456
  0.75  0.85       167        248.925       0.992
  0.85  0.95       132        462.681       0.534
 
 -0.05  0.95      1488        247.052       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               241.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.921
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     9.1 XMAX=     9.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59645 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        522


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               241.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.920
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     7.1 XMAX=  2437.7 NXBIN=    3
 YMIN=     7.3 YMAX=  2518.5 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        522


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               241.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.917
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     9.1 XMAX=     9.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1222    1263
 DETECTOR_SURFACE_POSITION=    1635    1689
 DETECTOR_SURFACE_POSITION=     810    1689
 DETECTOR_SURFACE_POSITION=     810     836
 DETECTOR_SURFACE_POSITION=    1635     836
 DETECTOR_SURFACE_POSITION=    2157    1663
 DETECTOR_SURFACE_POSITION=    1609    2228
 DETECTOR_SURFACE_POSITION=     835    2228
 DETECTOR_SURFACE_POSITION=     288    1663
 DETECTOR_SURFACE_POSITION=     288     863
 DETECTOR_SURFACE_POSITION=     835     298
 DETECTOR_SURFACE_POSITION=    1609     298
 DETECTOR_SURFACE_POSITION=    2157     863
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
 1.177E+00  2.236E-03   19.49


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

   24.453   3.397    17.44   1.48      2.61      2.98      58      84      26
    3.397   2.414    13.11   1.43      4.65      3.60     126     167      14
    2.414   1.974     7.66   0.89      5.56      5.71     188     220       2
    1.974   1.711     0.09   0.60    433.87    564.41     190     240       0
    1.711   1.531    -0.02   0.65   2111.40   2665.00     218     280       2
    1.531   1.398    -0.03   0.69    -99.90    -99.90     242     306       0
    1.398   1.295     0.15   0.76    421.73    470.70     100     146       0
   --------------------------------------------------------------------------
   24.453   1.295     3.72   0.83      5.58      5.91    1122    1443      44
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      10      10)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681  -0.1064E+02  0.1245E+02        -0.85           1 
   47    0    0     1.645   0.1363E+02  0.1271E+02         1.07           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     5.2% OF MEAN INTENSITY


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

     3.84          64      42      1285        3.3%       2.5%      3.0%       44   21.29      3.6%    99.9*     0    0.000       0
     2.73         105      67      2098        3.2%       4.0%      3.3%       76   18.42      5.7%    99.3*     0    0.000       0
     2.24         163      98      2665        3.7%       4.9%      4.0%      130   14.48      7.0%    98.9*     0    0.000       0
     1.94         164      93      3101        3.0%       6.4%      8.0%      142    6.56      9.0%    99.2*     0    0.000       0
     1.74         185     112      3481        3.2%     451.3%    584.5%      146    0.10    638.3%    -3.3      0    0.000       0
     1.58         218     131      3848        3.4%    4287.3%   5771.4%      174   -0.01   6063.2%    -6.1      0    0.000       0
     1.47         237     146      4143        3.5%    7941.6%   9026.4%      182   -0.10  11231.1%     4.2      0    0.000       0
     1.37         220     131      4428        3.0%     -99.9%    -99.9%      178    0.00    -99.9%    -1.2      0    0.000       0
     1.29          87      62      4743        1.3%     168.0%    187.8%       50    0.32    237.5%    15.3      0    0.000       0
    total        1443     882     29792        3.0%       5.6%      5.9%     1122    4.73      7.9%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1488
 NUMBER OF REJECTED MISFITS                                44
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1443
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    882



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

   WILSON LINE (using all data) : A=   7.532 B=  29.639 CORRELATION=  0.79
      #      RES      SS        <I>       log(<I>)       BO
      16     7.466   0.004  3.2573E+03       8.089     -62.1
      22     4.648   0.012  2.8647E+03       7.960     -18.5
      29     3.490   0.021  1.2957E+03       7.167       8.9
      32     2.993   0.028  6.0374E+02       6.403      20.2
      40     2.647   0.036  6.2864E+02       6.444      15.3
      53     2.395   0.044  3.7771E+02       5.934      18.3
      42     2.206   0.051  2.3479E+02       5.459      20.2
      45     2.056   0.059  1.0189E+02       4.624      24.6
      47     1.930   0.067  3.6922E-01      -0.996      63.5
      61     1.828   0.075  1.6288E+00       0.488      47.1
      52     1.739   0.083  5.2070E-01      -0.653      49.5
      57     1.664   0.090 -1.1277E+00    -999.999   -1000.0
      62     1.597   0.098  1.4747E+00       0.388      36.4
      75     1.535   0.106 -1.0776E+00    -999.999   -1000.0
      78     1.480   0.114 -1.2694E+00    -999.999   -1000.0
      58     1.431   0.122  1.4950E+00       0.402      29.2
      53     1.391   0.129 -3.7172E-01    -999.999   -1000.0
      42     1.351   0.137  3.7993E+00       1.335      22.6
      18     1.314   0.145  6.0621E+00       1.802      19.8


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.466        0.812        0.355        0.123
       3     4.648        1.138        0.657        0.282
       3     3.490        0.336        0.082        0.014
       4     2.993        1.060        0.461        0.152
       3     2.647        0.243        0.055        0.009
       5     2.395        0.734        0.339        0.118
       3     2.206        1.302        0.533        0.160
       8     2.056        0.287        0.097        0.025
       7     1.930      258.752      344.762    15332.858
      11     1.828       16.012       26.768       77.057
       4     1.597        4.030        4.116        4.112
       3     1.535       10.616        3.855       14.871
      90   overall       29.728        9.627     1268.978


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      11     7.466        1.279        1.822        2.128
      19     4.648        1.049        1.118        1.029
      26     3.490        0.836        0.631        0.443
      28     2.993        0.884        0.804        0.662
      37     2.647        1.123        1.187        1.137
      48     2.395        0.953        0.791        0.551
      39     2.206        0.789        0.598        0.394
      37     2.056        1.499        2.098        2.680
      40     1.930      170.787        2.551    19724.563
      50     1.828        8.181        6.733       35.364
      43     1.739       93.679      528.087     5552.400
      58     1.597       13.766       23.151      128.853
      55     1.431       13.647       17.485       83.595
      40     1.351        4.702        6.598       11.457
      18     1.314        3.887       12.545       34.295
     792   overall       16.831       33.085     1316.802

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.466  0.00  0.00  0.00  0.00  0.38  0.36  0.67  0.64  0.61  0.59
       3     4.648  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.490  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.993  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.647  1.34  0.96  0.80  0.70  0.64  0.59  0.56  1.06  1.01  0.98
       5     2.395  0.00  0.00  0.00  0.85  0.77  1.07  1.00  0.95  0.91  0.88
       3     2.206  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.056  1.51  1.81  1.50  1.32  1.20  1.11  1.05  0.99  0.95  1.10
       7     1.930  1.73  1.24  1.03  0.91  0.82  0.76  0.72  0.68  0.65  0.84
      11     1.828  2.20  1.58  1.31  1.15  1.05  0.97  0.91  0.87  0.83  0.80
       4     1.597  1.01  0.72  0.60  0.53  0.96  0.89  0.84  0.80  1.14  1.10
       3     1.535  1.34  0.96  0.80  0.70  1.28  1.19  1.12  1.06  1.01  0.98
      90   overall  0.76  0.61  0.51  0.49  0.51  0.51  0.50  0.49  0.49  0.52


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      11     7.466  1.91  2.01  2.10  2.21  1.85  1.61  1.44  1.32  1.23  1.15
      19     4.648  1.11  1.74  1.22  0.96  1.20  1.17  1.05  1.15  1.15  1.08
      26     3.490  0.40  0.42  0.45  0.58  0.59  0.77  0.92  0.98  0.91  0.85
      28     2.993  1.50  1.58  1.10  1.08  1.18  1.03  1.21  1.17  1.08  1.07
      37     2.647  0.85  1.04  1.36  1.15  1.03  1.02  1.02  1.03  1.05  1.03
      48     2.395  1.53  1.26  1.04  1.07  1.06  1.11  1.12  1.10  1.09  1.02
      39     2.206  1.35  1.13  0.99  0.93  0.91  0.97  0.92  1.02  0.99  1.01
      37     2.056  1.14  1.49  1.36  1.31  1.10  1.02  1.02  1.08  1.05  0.98
      40     1.930  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
      50     1.828  4.41  2.43  1.70  1.33  1.17  1.02  0.91  0.84  0.78  0.73
      43     1.739  5.13  2.69  1.88  1.48  1.30  1.13  1.02  0.97  0.90  0.85
      58     1.597  4.53  2.47  1.73  1.41  1.18  1.03  0.92  0.88  0.81  0.76
      55     1.431  4.20  2.21  1.54  1.21  1.11  1.01  0.94  0.86  0.80  0.81
      40     1.351  4.20  2.48  1.93  1.59  1.33  1.22  1.09  1.00  1.01  0.95
      18     1.314  2.92  1.53  1.29  1.35  1.41  1.35  1.32  1.21  1.12  1.05
     792   overall  2.10  1.32  1.02  0.88  0.80  0.74  0.71  0.69  0.66  0.63


 List of     76 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   41    8   18    1.40   63.67 -0.2367E+02  0.1579E+02 "alien"
   34   22    3    1.89   54.81  0.2024E+02  0.9492E+01 "alien"
   51    5   11    1.39   41.22 -0.1532E+02  0.9763E+01 "alien"
   29   27   19    1.41   38.79 -0.1442E+02  0.9054E+01 "alien"
   25   21   15    1.74   38.42  0.2000E+02  0.1132E+02 "alien"
   32   23   10    1.75   35.06  0.1826E+02  0.1572E+02 "alien"
   33   32   16    1.38   34.95 -0.1299E+02  0.1568E+02 "alien"
   45    7   16    1.39   33.43 -0.1243E+02  0.9784E+01 "alien"
   50    5   12    1.39   32.31 -0.1201E+02  0.9940E+01 "alien"
   38   33   12    1.39   30.70 -0.1141E+02  0.9985E+01 "alien"
   40    1   10    1.73   29.35  0.1528E+02  0.8506E+01 "alien"
   36   32   14    1.39   27.87 -0.1036E+02  0.1871E+02 "alien"
   25   19   13    1.90   27.53  0.1017E+02  0.8172E+01 "alien"
   30   20    8    1.96   23.62  0.8721E+01  0.9053E+01 "alien"
   34   24    7    1.76   23.45  0.1221E+02  0.8972E+01 "alien"
   39   33   11    1.39   21.97 -0.8167E+01  0.1063E+02 "alien"
   30   12   23    1.38   20.80 -0.7733E+01  0.7787E+01 "alien"
   38    9   20    1.38   20.76 -0.7717E+01  0.9804E+01 "alien"
   27   18   23    1.37   20.36 -0.7567E+01  0.1040E+02 "alien"
   52    2    3    1.48   19.96 -0.2534E+02  0.9304E+01 "alien"
   37   25    2    1.72   19.58  0.1020E+02  0.9965E+01 "alien"
   52    5   10    1.38   19.50 -0.7248E+01  0.9822E+01 "alien"
   30   28   19    1.38   19.18 -0.7129E+01  0.9196E+01 "alien"
   35   26    8    1.67   18.89 -0.2130E+02  0.1531E+02 "alien"
   24   13   24    1.40   18.30 -0.6803E+01  0.8128E+01 "alien"
   42    3   13    1.56   18.21 -0.1963E+02  0.1531E+02 "alien"
   47    2    9    1.54   18.14 -0.1955E+02  0.1556E+02 "alien"
   43    4   13    1.53   17.80 -0.1918E+02  0.1531E+02 "alien"
   36   25    4    1.74   17.72  0.9225E+01  0.9112E+01 "alien"
   36   28    9    1.58   17.55  0.2589E+02  0.1670E+02 "alien"
   43    7   17    1.40   17.14 -0.6371E+01  0.9394E+01 "alien"
   38    2   12    1.72   17.06  0.8882E+01  0.1235E+02 "alien"
   20   17   15    1.94   16.66  0.6150E+01  0.6915E+01 "alien"
   33    0   11    1.95   16.17  0.1844E+02  0.1347E+02 "alien"
   33   11   22    1.38   15.71 -0.5840E+01  0.8530E+01 "alien"
   26   19   12    1.93   15.59  0.5756E+01  0.8146E+01 "alien"
   28   23   14    1.69   15.53 -0.1752E+02  0.7926E+01 "alien"
   14   10   18    1.94   15.51  0.5726E+01  0.6110E+01 "alien"
   31    8   20    1.51   14.91 -0.1606E+02  0.6895E+01 "alien"
   32    1   12    1.93   14.83  0.5475E+01  0.1327E+02 "alien"
   26   23   16    1.64   14.77 -0.1665E+02  0.1409E+02 "alien"
   30   26   20    1.37   14.38 -0.5344E+01  0.1505E+02 "alien"
   35   24    5    1.77   14.00  0.7291E+01  0.8890E+01 "alien"
   27   19   11    1.95   13.73  0.5069E+01  0.8765E+01 "alien"
   31   30   18    1.38   13.63 -0.5068E+01  0.9601E+01 "alien"
   26   22   15    1.70   13.00  0.6767E+01  0.8074E+01 "alien"
   39    2    3    1.96   12.86  0.4748E+01  0.1077E+02 "alien"
   15   12   18    1.90   12.85  0.4743E+01  0.6465E+01 "alien"
   28   20   11    1.89   12.71  0.4693E+01  0.8330E+01 "alien"
   34   24    8    1.73   12.58  0.6549E+01  0.8663E+01 "alien"
   42    6   16    1.46   12.34  0.1845E+02  0.9153E+01 "alien"
   37    1    7    1.95   12.21  0.4509E+01  0.8719E+01 "alien"
   29   22   12    1.77   12.20  0.6351E+01  0.8204E+01 "alien"
   29   20   10    1.91   11.95  0.4411E+01  0.1337E+02 "alien"
   29   20    9    1.96   11.67  0.4308E+01  0.1351E+02 "alien"
   36   33   14    1.37   11.61 -0.4317E+01  0.1026E+02 "alien"
   46    7   16    1.37   11.36 -0.4221E+01  0.9294E+01 "alien"
   39   30    6    1.53   11.34 -0.1222E+02  0.1052E+02 "alien"
   30    3   15    1.82   11.30  0.1841E+02  0.1097E+02 "alien"
   40    7   18    1.43   11.23  0.1679E+02  0.1562E+02 "alien"
   49    4   11    1.44   11.22  0.1677E+02  0.9977E+01 "alien"
   26   13   24    1.38   11.11 -0.4129E+01  0.9532E+01 "alien"
   34   25    9    1.69   11.04 -0.1245E+02  0.8898E+01 "alien"
   53    7   13    1.30   11.02  0.6682E+02  0.1522E+02 "alien"
   35   30   13    1.46   10.90 -0.1384E+02  0.9471E+01 "alien"
   33   29   14    1.48   10.86 -0.1378E+02  0.9273E+01 "alien"
   43    2   10    1.63   10.85 -0.1224E+02  0.1311E+02 "alien"
   33   25   10    1.68   10.84 -0.1222E+02  0.9933E+01 "alien"
   38   31   10    1.46   10.64 -0.1351E+02  0.9669E+01 "alien"
   43    1    9    1.66   10.54 -0.1189E+02  0.8493E+01 "alien"
   39    2    0    1.98   10.49  0.8058E+01  0.1306E+02 "alien"
   37   31   12    1.43   10.43  0.1560E+02  0.1633E+02 "alien"
   35    4   16    1.63   10.37  0.1529E+02  0.1398E+02 "alien"
   42    8   18    1.38   10.36 -0.3853E+01  0.9115E+01 "alien"
   38   28    6    1.59   10.26  0.1513E+02  0.1373E+02 "alien"
   27   24   16    1.60   10.25  0.1511E+02  0.1401E+02 "alien"


 List of     76 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   11.02  0.6682E+02  0.1522E+02
   46    7   16    1.37   11.36 -0.4221E+01  0.9294E+01
   30   26   20    1.37   14.38 -0.5344E+01  0.1505E+02
   36   33   14    1.37   11.61 -0.4317E+01  0.1026E+02
   27   18   23    1.37   20.36 -0.7567E+01  0.1040E+02
   30   12   23    1.38   20.80 -0.7733E+01  0.7787E+01
   31   30   18    1.38   13.63 -0.5068E+01  0.9601E+01
   26   13   24    1.38   11.11 -0.4129E+01  0.9532E+01
   33   11   22    1.38   15.71 -0.5840E+01  0.8530E+01
   33   32   16    1.38   34.95 -0.1299E+02  0.1568E+02
   52    5   10    1.38   19.50 -0.7248E+01  0.9822E+01
   30   28   19    1.38   19.18 -0.7129E+01  0.9196E+01
   38    9   20    1.38   20.76 -0.7717E+01  0.9804E+01
   42    8   18    1.38   10.36 -0.3853E+01  0.9115E+01
   51    5   11    1.39   41.22 -0.1532E+02  0.9763E+01
   38   33   12    1.39   30.70 -0.1141E+02  0.9985E+01
   36   32   14    1.39   27.87 -0.1036E+02  0.1871E+02
   50    5   12    1.39   32.31 -0.1201E+02  0.9940E+01
   45    7   16    1.39   33.43 -0.1243E+02  0.9784E+01
   39   33   11    1.39   21.97 -0.8167E+01  0.1063E+02
   43    7   17    1.40   17.14 -0.6371E+01  0.9394E+01
   24   13   24    1.40   18.30 -0.6803E+01  0.8128E+01
   41    8   18    1.40   63.67 -0.2367E+02  0.1579E+02
   29   27   19    1.41   38.79 -0.1442E+02  0.9054E+01
   40    7   18    1.43   11.23  0.1679E+02  0.1562E+02
   37   31   12    1.43   10.43  0.1560E+02  0.1633E+02
   49    4   11    1.44   11.22  0.1677E+02  0.9977E+01
   42    6   16    1.46   12.34  0.1845E+02  0.9153E+01
   38   31   10    1.46   10.64 -0.1351E+02  0.9669E+01
   35   30   13    1.46   10.90 -0.1384E+02  0.9471E+01
   52    2    3    1.48   19.96 -0.2534E+02  0.9304E+01
   33   29   14    1.48   10.86 -0.1378E+02  0.9273E+01
   31    8   20    1.51   14.91 -0.1606E+02  0.6895E+01
   39   30    6    1.53   11.34 -0.1222E+02  0.1052E+02
   43    4   13    1.53   17.80 -0.1918E+02  0.1531E+02
   47    2    9    1.54   18.14 -0.1955E+02  0.1556E+02
   42    3   13    1.56   18.21 -0.1963E+02  0.1531E+02
   36   28    9    1.58   17.55  0.2589E+02  0.1670E+02
   38   28    6    1.59   10.26  0.1513E+02  0.1373E+02
   27   24   16    1.60   10.25  0.1511E+02  0.1401E+02
   35    4   16    1.63   10.37  0.1529E+02  0.1398E+02
   43    2   10    1.63   10.85 -0.1224E+02  0.1311E+02
   26   23   16    1.64   14.77 -0.1665E+02  0.1409E+02
   43    1    9    1.66   10.54 -0.1189E+02  0.8493E+01
   35   26    8    1.67   18.89 -0.2130E+02  0.1531E+02
   33   25   10    1.68   10.84 -0.1222E+02  0.9933E+01
   34   25    9    1.69   11.04 -0.1245E+02  0.8898E+01
   28   23   14    1.69   15.53 -0.1752E+02  0.7926E+01
   26   22   15    1.70   13.00  0.6767E+01  0.8074E+01
   38    2   12    1.72   17.06  0.8882E+01  0.1235E+02
   37   25    2    1.72   19.58  0.1020E+02  0.9965E+01
   40    1   10    1.73   29.35  0.1528E+02  0.8506E+01
   34   24    8    1.73   12.58  0.6549E+01  0.8663E+01
   36   25    4    1.74   17.72  0.9225E+01  0.9112E+01
   25   21   15    1.74   38.42  0.2000E+02  0.1132E+02
   32   23   10    1.75   35.06  0.1826E+02  0.1572E+02
   34   24    7    1.76   23.45  0.1221E+02  0.8972E+01
   35   24    5    1.77   14.00  0.7291E+01  0.8890E+01
   29   22   12    1.77   12.20  0.6351E+01  0.8204E+01
   30    3   15    1.82   11.30  0.1841E+02  0.1097E+02
   34   22    3    1.89   54.81  0.2024E+02  0.9492E+01
   28   20   11    1.89   12.71  0.4693E+01  0.8330E+01
   15   12   18    1.90   12.85  0.4743E+01  0.6465E+01
   25   19   13    1.90   27.53  0.1017E+02  0.8172E+01
   29   20   10    1.91   11.95  0.4411E+01  0.1337E+02
   26   19   12    1.93   15.59  0.5756E+01  0.8146E+01
   32    1   12    1.93   14.83  0.5475E+01  0.1327E+02
   14   10   18    1.94   15.51  0.5726E+01  0.6110E+01
   20   17   15    1.94   16.66  0.6150E+01  0.6915E+01
   27   19   11    1.95   13.73  0.5069E+01  0.8765E+01
   33    0   11    1.95   16.17  0.1844E+02  0.1347E+02
   37    1    7    1.95   12.21  0.4509E+01  0.8719E+01
   29   20    9    1.96   11.67  0.4308E+01  0.1351E+02
   39    2    3    1.96   12.86  0.4748E+01  0.1077E+02
   30   20    8    1.96   23.62  0.8721E+01  0.9053E+01
   39    2    0    1.98   10.49  0.8058E+01  0.1306E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        19
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1488
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1412
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       76


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
