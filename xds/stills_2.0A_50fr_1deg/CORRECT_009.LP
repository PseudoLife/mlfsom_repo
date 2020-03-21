
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.339    77.339    38.707  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=       9       9
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       9
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  1.000000 -0.000714  0.000699
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
 REFLECTING_RANGE_E.S.D.=   0.089
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

    1479 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       9 ...       9
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1479 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               46
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
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0  1  0  0  1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  2  0  1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    27        mC        500.0     159.4   38.7  109.4  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.7  109.4  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.4  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    30        mC        998.1      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.1      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    41        mC        998.1     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    22        hP        998.1      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0    1 -1 -1  0 -1  1 -1  0  1  1  1  0
     2        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.4   86.5   86.5  78.4  50.8  50.8    1  1  0  0  1  0 -1  0  0  1 -1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    1  0 -1  0  0  1 -1  0  1  1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0 -1  1  0 -1  0  1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.5  247.6 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.5   86.5   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.2  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.5   86.5  177.2 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
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

  *   0.00        0    6.7      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.710   -23.409     0.034
 UNIT_CELL_B-AXIS=    23.409    73.710     0.015
 UNIT_CELL_C-AXIS=    -0.019    -0.002    38.707
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       329 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.21
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.316    77.316    38.695  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.1E-01 1.1E-01 4.8E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012934  0.012934  0.025843  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.683   -23.421     0.027
 COORDINATES OF UNIT CELL B-AXIS    23.421    73.683     0.010
 COORDINATES OF UNIT CELL C-AXIS    -0.014    -0.001    38.695
 CRYSTAL MOSAICITY (DEGREES)     0.089
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006359  0.001947
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.88   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1261.05   1242.13
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.22
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.32    77.32    38.69  90.000  90.000  90.000


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
  0.05  0.15       114        233.907       1.109
  0.15  0.25       179        236.147       1.098
  0.25  0.35       191        228.879       1.133
  0.35  0.45       146        146.263       1.773
  0.45  0.55       180        326.386       0.795
  0.55  0.65       196        249.267       1.040
  0.65  0.75       181        218.629       1.186
  0.75  0.85       166        220.694       1.175
  0.85  0.95       126        521.602       0.497
 
 -0.05  0.95      1479        259.324       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               244.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.864
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     8.1 XMAX=     8.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59659 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        532


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               244.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.860
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    30.9 XMAX=  2437.8 NXBIN=    3
 YMIN=     7.2 YMAX=  2518.7 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        532


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               244.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.848
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     8.1 XMAX=     8.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1234    1263
 DETECTOR_SURFACE_POSITION=    1643    1690
 DETECTOR_SURFACE_POSITION=     826    1690
 DETECTOR_SURFACE_POSITION=     826     836
 DETECTOR_SURFACE_POSITION=    1643     836
 DETECTOR_SURFACE_POSITION=    2159    1663
 DETECTOR_SURFACE_POSITION=    1618    2228
 DETECTOR_SURFACE_POSITION=     851    2228
 DETECTOR_SURFACE_POSITION=     309    1663
 DETECTOR_SURFACE_POSITION=     309     863
 DETECTOR_SURFACE_POSITION=     851     298
 DETECTOR_SURFACE_POSITION=    1618     298
 DETECTOR_SURFACE_POSITION=    2159     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        532


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
 1.155E+00  2.412E-03   18.94


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

   24.449   3.397    16.77   1.50      2.76      3.06      64      89      20
    3.397   2.414    13.18   1.49      4.53      3.59     128     168      14
    2.414   1.974     8.03   0.69      4.94      5.52     176     209      12
    1.974   1.711    -0.07   0.69    -99.90    -99.90     190     240       0
    1.711   1.531     0.02   0.76    684.42    790.08     218     276       0
    1.531   1.398     0.04   0.70   1406.82   1688.13     238     303       2
    1.398   1.295     0.03   0.75   1777.81   2049.63     100     145       0
   --------------------------------------------------------------------------
   24.449   1.295     3.77   0.85      5.49      5.78    1114    1430      48
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       9       9)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681  -0.2442E+01  0.1255E+02        -0.19           1 
   47    0    0     1.645   0.1270E+00  0.1270E+02         0.01           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     0.0% OF MEAN INTENSITY


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

     3.84          67      43      1285        3.3%       2.8%      3.1%       48   20.88      3.9%    99.9*     0    0.000       0
     2.73         107      68      2098        3.2%       3.8%      3.4%       78   18.03      5.4%    99.2*     0    0.000       0
     2.24         157      95      2665        3.6%       4.5%      4.0%      124   14.83      6.3%    99.3*     0    0.000       0
     1.94         160      91      3101        2.9%       6.4%      7.1%      138    7.13      9.0%    99.3*     0    0.000       0
     1.74         185     112      3481        3.2%     -99.9%    -99.9%      146    0.00    -99.9%     1.7      0    0.000       0
     1.58         213     126      3848        3.3%    1824.7%   2171.5%      174   -0.07   2580.5%    -7.0      0    0.000       0
     1.47         237     147      4143        3.5%     878.5%   1012.4%      180    0.10   1242.4%    -3.8      0    0.000       0
     1.37         218     130      4429        2.9%    1161.4%   1396.9%      176    0.00   1642.5%     9.3      0    0.000       0
     1.29          86      61      4743        1.3%    2342.2%   2500.7%       50    0.07   3312.4%    25.9      0    0.000       0
    total        1430     873     29793        2.9%       5.5%      5.8%     1114    4.79      7.8%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1479
 NUMBER OF REJECTED MISFITS                                48
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1430
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    873



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

   WILSON LINE (using all data) : A=   8.489 B=  35.350 CORRELATION=  0.92
      #      RES      SS        <I>       log(<I>)       BO
      18     7.398   0.005  2.9969E+03       8.005      53.0
      21     4.621   0.012  2.9040E+03       7.974      22.0
      30     3.495   0.020  1.3169E+03       7.183      31.9
      32     2.987   0.028  6.3444E+02       6.453      36.3
      40     2.647   0.036  6.8514E+02       6.530      27.5
      51     2.394   0.044  3.8638E+02       5.957      29.0
      41     2.204   0.051  2.6800E+02       5.591      28.2
      43     2.054   0.059  1.1814E+02       4.772      31.4
      47     1.930   0.067 -1.3837E+00    -999.999   -1000.0
      61     1.828   0.075 -1.4485E+00    -999.999   -1000.0
      52     1.739   0.083  9.1203E-01      -0.092      51.9
      54     1.665   0.090 -9.7606E-01    -999.999   -1000.0
      60     1.595   0.098 -2.4655E-01    -999.999   -1000.0
      76     1.535   0.106  1.5566E+00       0.443      37.9
      78     1.480   0.114  2.8056E-01      -1.271      42.7
      57     1.432   0.122  4.6888E-01      -0.757      37.9
      53     1.390   0.129 -6.5634E-02    -999.999   -1000.0
      41     1.351   0.137 -4.4269E-01    -999.999   -1000.0
      18     1.314   0.145  3.5801E+00       1.275      24.9


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.398        0.983        0.468        0.176
       3     4.621        1.181        0.692        0.302
       3     3.495        0.360        0.092        0.017
       4     2.987        1.038        0.448        0.147
       3     2.647        0.250        0.058        0.010
       5     2.394        0.856        0.421        0.155
       3     2.204        1.318        0.528        0.152
       8     2.054        0.303        0.103        0.027
      11     1.828       17.439       19.561       43.614
       9     1.739       25.803       13.660      128.477
      10     1.665       24.491       29.205      121.160
       3     1.595        8.669        9.375        7.473
       3     1.535        0.541        0.105        0.029
       6     1.480      199.480     1367.286     7108.489
      88   overall       25.632       33.930      805.272


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      13     7.398        1.393        2.275        3.038
      18     4.621        1.082        1.201        1.135
      27     3.495        0.851        0.652        0.460
      28     2.987        0.905        0.837        0.700
      37     2.647        1.111        1.158        1.086
      46     2.394        0.893        0.701        0.462
      38     2.204        0.775        0.577        0.372
      35     2.054        1.462        1.920        2.286
      43     1.739       23.252       14.520      265.258
      73     1.535       14.892       26.088      142.352
      54     1.432      220.065      328.300    20111.294
      18     1.314       11.813       51.509      233.981
     785   overall       61.442     -113.888    13151.967

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.398  0.00  0.00  0.00  0.00  0.00  0.36  0.33  0.64  0.61  0.59
       3     4.621  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.495  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.987  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.647  1.34  0.96  0.80  0.70  0.64  0.59  0.56  1.06  1.01  0.98
       5     2.394  0.00  0.00  0.00  0.00  1.15  1.07  1.00  0.95  0.91  0.88
       3     2.204  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.054  2.01  1.45  1.20  1.32  1.20  1.11  1.05  0.99  0.95  0.92
      11     1.828  2.56  1.84  1.53  1.34  1.22  1.13  1.07  1.01  1.11  1.07
       9     1.739  1.34  0.96  0.80  0.70  0.64  0.79  0.74  0.71  0.68  0.65
      10     1.665  1.61  1.16  0.96  0.85  0.77  0.71  0.67  0.64  0.61  0.59
       3     1.595  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       3     1.535  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       6     1.480  1.34  0.96  0.80  0.70  0.64  0.89  0.84  0.80  0.76  0.73
      88   overall  1.24  0.89  0.74  0.67  0.68  0.71  0.67  0.67  0.66  0.65


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      13     7.398  2.42  2.55  1.78  2.10  1.95  1.70  1.53  1.40  1.30  1.22
      18     4.621  1.75  1.84  1.29  1.01  1.27  1.23  1.10  1.01  1.22  1.14
      27     3.495  0.78  0.61  0.57  0.67  0.66  0.74  0.88  1.01  0.94  0.88
      28     2.987  1.50  1.77  1.24  1.19  1.09  1.03  1.21  1.17  1.08  1.07
      37     2.647  0.85  1.19  1.25  1.07  0.96  0.90  1.02  1.08  1.05  1.03
      46     2.394  1.37  1.08  1.01  1.12  1.05  1.11  1.08  1.07  1.06  1.03
      38     2.204  1.66  1.16  0.91  0.80  0.94  0.99  0.94  1.00  1.11  1.08
      35     2.054  1.50  1.42  1.21  1.21  1.09  1.14  1.08  1.04  1.06  0.99
      43     1.739  4.40  2.31  1.70  1.34  1.12  0.98  0.92  0.84  0.78  0.74
      73     1.535  4.32  2.34  1.69  1.33  1.11  0.97  0.87  0.82  0.78  0.74
      54     1.432  4.86  2.55  1.79  1.40  1.22  1.07  0.99  0.91  0.84  0.79
      18     1.314  4.09  2.76  2.14  1.69  1.41  1.23  1.10  1.01  1.03  0.97
     785   overall  1.50  0.99  0.76  0.66  0.61  0.57  0.55  0.54  0.53  0.51


 List of    131 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   37   33   13    1.38  308.38 -0.2024E+02  0.1009E+02 "alien"
   41    8   18    1.40  302.18 -0.1983E+02  0.1606E+02 "alien"
   36   32   14    1.39  241.02 -0.1582E+02  0.1741E+02 "alien"
   39   33   11    1.39  222.28 -0.1459E+02  0.1028E+02 "alien"
   34   10   21    1.41  191.15 -0.1255E+02  0.1349E+02 "alien"
   51    4   10    1.41  191.07 -0.1254E+02  0.1509E+02 "alien"
   46    7   16    1.37  166.43 -0.1092E+02  0.9188E+01 "alien"
   50    5   11    1.41  148.06 -0.9718E+01  0.1663E+02 "alien"
   27   18   23    1.37  142.54 -0.9356E+01  0.1063E+02 "alien"
   36   10   21    1.38  124.73 -0.8187E+01  0.8059E+01 "alien"
   27   20   22    1.40  114.70 -0.7528E+01  0.7794E+01 "alien"
   52    5   10    1.38   89.58 -0.5879E+01  0.9629E+01 "alien"
   33   32   16    1.38   87.22 -0.5725E+01  0.1616E+02 "alien"
   43    7   17    1.40   83.78 -0.5499E+01  0.9027E+01 "alien"
   53    5    9    1.38   82.75 -0.5432E+01  0.1049E+02 "alien"
   42    3   12    1.60   79.42 -0.1958E+02  0.1442E+02 "alien"
   28   21   22    1.38   76.02 -0.4989E+01  0.8002E+01 "alien"
   29   27   19    1.41   75.08 -0.4928E+01  0.8820E+01 "alien"
   29   12   23    1.39   68.94 -0.4525E+01  0.9323E+01 "alien"
   45    7   16    1.39   68.89 -0.4522E+01  0.9495E+01 "alien"
   37   10   21    1.36   65.92 -0.2918E+02  0.1192E+02 "alien"
   39   31    8    1.48   65.05  0.1825E+02  0.9857E+01 "alien"
   51    3    7    1.46   60.31  0.1692E+02  0.9527E+01 "alien"
   39   28    4    1.59   60.24 -0.1485E+02  0.1604E+02 "alien"
   45    1    8    1.62   57.77 -0.1424E+02  0.9309E+01 "alien"
   50    2    6    1.50   57.49  0.1613E+02  0.9511E+01 "alien"
   28   13   24    1.36   54.65 -0.2419E+02  0.1156E+02 "alien"
   48    6   14    1.38   51.68 -0.3392E+01  0.1305E+02 "alien"
   36   33   14    1.37   50.39 -0.3307E+01  0.1011E+02 "alien"
   45    5   14    1.45   48.77  0.2287E+02  0.9502E+01 "alien"
   37    4   15    1.62   47.15 -0.1162E+02  0.8595E+01 "alien"
   48    3    9    1.51   46.77  0.1312E+02  0.1427E+02 "alien"
   38   28    5    1.60   46.29 -0.1141E+02  0.9140E+01 "alien"
   50    3    8    1.47   44.18  0.1240E+02  0.9472E+01 "alien"
   46    5   14    1.43   42.90  0.2012E+02  0.1407E+02 "alien"
   37   31   11    1.46   41.90  0.1176E+02  0.1024E+02 "alien"
   33    8   20    1.47   41.58  0.1166E+02  0.1147E+02 "alien"
   31    6   18    1.62   41.05 -0.1012E+02  0.1138E+02 "alien"
   34   33   16    1.35   39.17 -0.1734E+02  0.9794E+01 "alien"
   46    4   12    1.49   37.14  0.1042E+02  0.9142E+01 "alien"
   31   30   18    1.38   37.14 -0.2437E+01  0.9687E+01 "alien"
   25   20   21    1.46   36.64  0.1028E+02  0.1106E+02 "alien"
   54    5    8    1.37   35.36 -0.1565E+02  0.1150E+02 "alien"
   35   31   14    1.42   35.05  0.1643E+02  0.1391E+02 "alien"
   34    7   19    1.50   35.03  0.9828E+01  0.1159E+02 "alien"
   52    2    4    1.47   34.54  0.9691E+01  0.1585E+02 "alien"
   21    9   21    1.62   34.48 -0.8502E+01  0.7731E+01 "alien"
   23   13   23    1.46   33.66  0.9442E+01  0.8825E+01 "alien"
   50    6   13    1.36   33.20 -0.1470E+02  0.9816E+01 "alien"
   32   30   16    1.42   32.89  0.1542E+02  0.9512E+01 "alien"
   21   10   23    1.50   32.08  0.9000E+01  0.8012E+01 "alien"
   30   29   18    1.40   31.61 -0.2075E+01  0.1410E+02 "alien"
   49    6   14    1.36   30.79 -0.1363E+02  0.1586E+02 "alien"
   38    9   20    1.38   30.58 -0.2007E+01  0.9752E+01 "alien"
   48    1    3    1.60   29.96 -0.7386E+01  0.8953E+01 "alien"
   41   31    4    1.49   29.58  0.8300E+01  0.1028E+02 "alien"
   43    9   19    1.33   29.09 -0.1288E+02  0.1366E+02 "alien"
   39    4   15    1.57   28.31 -0.6980E+01  0.1098E+02 "alien"
   41   31    3    1.49   27.82  0.7805E+01  0.9904E+01 "alien"
   22   13   24    1.42   27.79  0.1303E+02  0.7962E+01 "alien"
   33   26   12    1.60   27.53 -0.6787E+01  0.1301E+02 "alien"
   37    5   16    1.57   27.44 -0.6766E+01  0.1343E+02 "alien"
   26   23   20    1.46   26.84  0.7531E+01  0.1247E+02 "alien"
   36    9   20    1.42   25.98  0.1218E+02  0.1157E+02 "alien"
   29   20   23    1.34   24.45 -0.1083E+02  0.1101E+02 "alien"
   20   13   21    1.60   24.41 -0.6018E+01  0.7097E+01 "alien"
   39   32    9    1.44   23.63  0.1108E+02  0.1027E+02 "alien"
   39   29    4    1.57   23.31 -0.5747E+01  0.9416E+01 "alien"
   52    3    6    1.45   23.01  0.1079E+02  0.9906E+01 "alien"
   48    5   13    1.41   22.82  0.1070E+02  0.9783E+01 "alien"
   41    9   19    1.37   22.65 -0.1003E+02  0.1237E+02 "alien"
   27   25   19    1.46   22.61  0.6344E+01  0.9358E+01 "alien"
   33   12   23    1.34   22.61 -0.1001E+02  0.1436E+02 "alien"
   35   27   10    1.59   22.01 -0.5426E+01  0.8807E+01 "alien"
   54    3    1    1.43   22.01  0.1032E+02  0.9887E+01 "alien"
   49    5   12    1.41   21.95  0.1029E+02  0.9585E+01 "alien"
   29   24   14    1.65   21.50 -0.2099E+02  0.1163E+02 "alien"
   21    7   24    1.46   21.04  0.5904E+01  0.6803E+01 "alien"
   53    3    5    1.43   20.99  0.9842E+01  0.9673E+01 "alien"
   41   32    5    1.46   20.92  0.5869E+01  0.1009E+02 "alien"
   43    5   15    1.47   20.86  0.5853E+01  0.1013E+02 "alien"
   49    4   11    1.44   20.53  0.9628E+01  0.9619E+01 "alien"
   35   26    8    1.66   20.37 -0.1988E+02  0.1603E+02 "alien"
   32    8   20    1.49   19.76  0.5543E+01  0.8192E+01 "alien"
   37   32   12    1.42   19.20  0.9001E+01  0.1009E+02 "alien"
   48    1    4    1.59   18.40 -0.4537E+01  0.9096E+01 "alien"
   32   29   15    1.47   18.33  0.5144E+01  0.9467E+01 "alien"
   54    3    3    1.42   18.31  0.8587E+01  0.9806E+01 "alien"
   31   27   15    1.52   18.15  0.2825E+02  0.1565E+02 "alien"
   43    2   10    1.63   18.12 -0.1769E+02  0.1340E+02 "alien"
   26   21   14    1.77   17.93  0.1635E+02  0.7804E+01 "alien"
   49    3    8    1.50   17.91  0.5024E+01  0.1700E+02 "alien"
   42    5   15    1.49   17.84  0.5004E+01  0.8736E+01 "alien"
   53    7   13    1.30   17.75  0.6354E+02  0.1537E+02 "alien"
   41    3   13    1.59   16.90 -0.4167E+01  0.8342E+01 "alien"
   27   22   21    1.42   16.81  0.7883E+01  0.8484E+01 "alien"
   30   29   17    1.44   16.76  0.7857E+01  0.9111E+01 "alien"
   51    6   12    1.36   16.44 -0.7277E+01  0.1102E+02 "alien"
   46    1    6    1.63   15.85 -0.3909E+01  0.1277E+02 "alien"
   46    1    7    1.61   15.79 -0.3893E+01  0.8839E+01 "alien"
   38   32   11    1.42   15.70  0.7363E+01  0.9818E+01 "alien"
   22    9   21    1.60   15.15 -0.3734E+01  0.8387E+01 "alien"
   28   26   19    1.44   14.92  0.6995E+01  0.9049E+01 "alien"
   34   25    9    1.69   14.87 -0.1452E+02  0.8876E+01 "alien"
   27   27   18    1.47   14.72  0.1146E+02  0.9552E+01 "alien"
   47    1    6    1.59   14.66 -0.3616E+01  0.1277E+02 "alien"
   40   30    2    1.54   14.54  0.2263E+02  0.1109E+02 "alien"
   30   12   23    1.38   13.76 -0.9032E+00  0.7479E+01 "alien"
   21   13   24    1.43   13.75  0.6447E+01  0.7213E+01 "alien"
   33   33   17    1.34   13.09 -0.1434E+02  0.1405E+02 "alien"
   51    2    3    1.50   13.00  0.3648E+01  0.9453E+01 "alien"
   32   30   19    1.33   12.81 -0.5672E+01  0.1424E+02 "alien"
   53    2    0    1.46   12.79  0.8687E+01  0.1537E+02 "alien"
   31   25   13    1.63   12.56 -0.3096E+01  0.8296E+01 "alien"
   41   31    0    1.50   12.47  0.8298E+01  0.1419E+02 "alien"
   23   21   16    1.73   12.44  0.1134E+02  0.1206E+02 "alien"
   37   27    6    1.63   12.44 -0.1214E+02  0.1262E+02 "alien"
   43    8   18    1.37   12.35 -0.5467E+01  0.9523E+01 "alien"
   37   26    5    1.67   11.96 -0.1167E+02  0.8921E+01 "alien"
   41    5   15    1.52   11.53  0.1795E+02  0.1374E+02 "alien"
   32   28   14    1.52   11.53  0.1794E+02  0.1768E+02 "alien"
   37   25    1    1.73   11.26  0.1027E+02  0.9404E+01 "alien"
   44    1    1    1.75   11.03  0.1006E+02  0.9414E+01 "alien"
   34   22    4    1.87   10.89 -0.1577E+02  0.1015E+02 "alien"
   37    3   14    1.66   10.69 -0.1044E+02  0.7636E+01 "alien"
   30   21   10    1.85   10.51 -0.1523E+02  0.1209E+02 "alien"
   39   28    1    1.61   10.31 -0.2543E+01  0.9297E+01 "alien"
   51    2    4    1.50   10.15  0.2848E+01  0.9282E+01 "alien"
   37    7   18    1.48   10.11  0.2836E+01  0.1357E+02 "alien"
   39    2    1    1.98   10.10 -0.1398E+02  0.9136E+01 "alien"
   21   15   21    1.57   10.09 -0.2488E+01  0.8940E+01 "alien"


 List of    131 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   17.75  0.6354E+02  0.1537E+02
   43    9   19    1.33   29.09 -0.1288E+02  0.1366E+02
   32   30   19    1.33   12.81 -0.5672E+01  0.1424E+02
   29   20   23    1.34   24.45 -0.1083E+02  0.1101E+02
   33   12   23    1.34   22.61 -0.1001E+02  0.1436E+02
   33   33   17    1.34   13.09 -0.1434E+02  0.1405E+02
   34   33   16    1.35   39.17 -0.1734E+02  0.9794E+01
   28   13   24    1.36   54.65 -0.2419E+02  0.1156E+02
   37   10   21    1.36   65.92 -0.2918E+02  0.1192E+02
   49    6   14    1.36   30.79 -0.1363E+02  0.1586E+02
   51    6   12    1.36   16.44 -0.7277E+01  0.1102E+02
   50    6   13    1.36   33.20 -0.1470E+02  0.9816E+01
   43    8   18    1.37   12.35 -0.5467E+01  0.9523E+01
   41    9   19    1.37   22.65 -0.1003E+02  0.1237E+02
   54    5    8    1.37   35.36 -0.1565E+02  0.1150E+02
   46    7   16    1.37  166.43 -0.1092E+02  0.9188E+01
   36   33   14    1.37   50.39 -0.3307E+01  0.1011E+02
   27   18   23    1.37  142.54 -0.9356E+01  0.1063E+02
   53    5    9    1.38   82.75 -0.5432E+01  0.1049E+02
   28   21   22    1.38   76.02 -0.4989E+01  0.8002E+01
   36   10   21    1.38  124.73 -0.8187E+01  0.8059E+01
   30   12   23    1.38   13.76 -0.9032E+00  0.7479E+01
   31   30   18    1.38   37.14 -0.2437E+01  0.9687E+01
   33   32   16    1.38   87.22 -0.5725E+01  0.1616E+02
   37   33   13    1.38  308.38 -0.2024E+02  0.1009E+02
   52    5   10    1.38   89.58 -0.5879E+01  0.9629E+01
   48    6   14    1.38   51.68 -0.3392E+01  0.1305E+02
   38    9   20    1.38   30.58 -0.2007E+01  0.9752E+01
   36   32   14    1.39  241.02 -0.1582E+02  0.1741E+02
   29   12   23    1.39   68.94 -0.4525E+01  0.9323E+01
   45    7   16    1.39   68.89 -0.4522E+01  0.9495E+01
   39   33   11    1.39  222.28 -0.1459E+02  0.1028E+02
   27   20   22    1.40  114.70 -0.7528E+01  0.7794E+01
   43    7   17    1.40   83.78 -0.5499E+01  0.9027E+01
   41    8   18    1.40  302.18 -0.1983E+02  0.1606E+02
   30   29   18    1.40   31.61 -0.2075E+01  0.1410E+02
   34   10   21    1.41  191.15 -0.1255E+02  0.1349E+02
   51    4   10    1.41  191.07 -0.1254E+02  0.1509E+02
   29   27   19    1.41   75.08 -0.4928E+01  0.8820E+01
   50    5   11    1.41  148.06 -0.9718E+01  0.1663E+02
   48    5   13    1.41   22.82  0.1070E+02  0.9783E+01
   49    5   12    1.41   21.95  0.1029E+02  0.9585E+01
   36    9   20    1.42   25.98  0.1218E+02  0.1157E+02
   27   22   21    1.42   16.81  0.7883E+01  0.8484E+01
   35   31   14    1.42   35.05  0.1643E+02  0.1391E+02
   37   32   12    1.42   19.20  0.9001E+01  0.1009E+02
   54    3    3    1.42   18.31  0.8587E+01  0.9806E+01
   22   13   24    1.42   27.79  0.1303E+02  0.7962E+01
   38   32   11    1.42   15.70  0.7363E+01  0.9818E+01
   32   30   16    1.42   32.89  0.1542E+02  0.9512E+01
   54    3    1    1.43   22.01  0.1032E+02  0.9887E+01
   46    5   14    1.43   42.90  0.2012E+02  0.1407E+02
   53    3    5    1.43   20.99  0.9842E+01  0.9673E+01
   21   13   24    1.43   13.75  0.6447E+01  0.7213E+01
   28   26   19    1.44   14.92  0.6995E+01  0.9049E+01
   49    4   11    1.44   20.53  0.9628E+01  0.9619E+01
   30   29   17    1.44   16.76  0.7857E+01  0.9111E+01
   39   32    9    1.44   23.63  0.1108E+02  0.1027E+02
   52    3    6    1.45   23.01  0.1079E+02  0.9906E+01
   45    5   14    1.45   48.77  0.2287E+02  0.9502E+01
   37   31   11    1.46   41.90  0.1176E+02  0.1024E+02
   53    2    0    1.46   12.79  0.8687E+01  0.1537E+02
   23   13   23    1.46   33.66  0.9442E+01  0.8825E+01
   51    3    7    1.46   60.31  0.1692E+02  0.9527E+01
   41   32    5    1.46   20.92  0.5869E+01  0.1009E+02
   26   23   20    1.46   26.84  0.7531E+01  0.1247E+02
   27   25   19    1.46   22.61  0.6344E+01  0.9358E+01
   21    7   24    1.46   21.04  0.5904E+01  0.6803E+01
   25   20   21    1.46   36.64  0.1028E+02  0.1106E+02
   43    5   15    1.47   20.86  0.5853E+01  0.1013E+02
   52    2    4    1.47   34.54  0.9691E+01  0.1585E+02
   50    3    8    1.47   44.18  0.1240E+02  0.9472E+01
   32   29   15    1.47   18.33  0.5144E+01  0.9467E+01
   27   27   18    1.47   14.72  0.1146E+02  0.9552E+01
   33    8   20    1.47   41.58  0.1166E+02  0.1147E+02
   39   31    8    1.48   65.05  0.1825E+02  0.9857E+01
   37    7   18    1.48   10.11  0.2836E+01  0.1357E+02
   46    4   12    1.49   37.14  0.1042E+02  0.9142E+01
   41   31    4    1.49   29.58  0.8300E+01  0.1028E+02
   42    5   15    1.49   17.84  0.5004E+01  0.8736E+01
   32    8   20    1.49   19.76  0.5543E+01  0.8192E+01
   41   31    3    1.49   27.82  0.7805E+01  0.9904E+01
   51    2    4    1.50   10.15  0.2848E+01  0.9282E+01
   49    3    8    1.50   17.91  0.5024E+01  0.1700E+02
   21   10   23    1.50   32.08  0.9000E+01  0.8012E+01
   50    2    6    1.50   57.49  0.1613E+02  0.9511E+01
   34    7   19    1.50   35.03  0.9828E+01  0.1159E+02
   41   31    0    1.50   12.47  0.8298E+01  0.1419E+02
   51    2    3    1.50   13.00  0.3648E+01  0.9453E+01
   48    3    9    1.51   46.77  0.1312E+02  0.1427E+02
   41    5   15    1.52   11.53  0.1795E+02  0.1374E+02
   32   28   14    1.52   11.53  0.1794E+02  0.1768E+02
   31   27   15    1.52   18.15  0.2825E+02  0.1565E+02
   40   30    2    1.54   14.54  0.2263E+02  0.1109E+02
   39    4   15    1.57   28.31 -0.6980E+01  0.1098E+02
   21   15   21    1.57   10.09 -0.2488E+01  0.8940E+01
   39   29    4    1.57   23.31 -0.5747E+01  0.9416E+01
   37    5   16    1.57   27.44 -0.6766E+01  0.1343E+02
   48    1    4    1.59   18.40 -0.4537E+01  0.9096E+01
   39   28    4    1.59   60.24 -0.1485E+02  0.1604E+02
   41    3   13    1.59   16.90 -0.4167E+01  0.8342E+01
   47    1    6    1.59   14.66 -0.3616E+01  0.1277E+02
   35   27   10    1.59   22.01 -0.5426E+01  0.8807E+01
   42    3   12    1.60   79.42 -0.1958E+02  0.1442E+02
   48    1    3    1.60   29.96 -0.7386E+01  0.8953E+01
   33   26   12    1.60   27.53 -0.6787E+01  0.1301E+02
   20   13   21    1.60   24.41 -0.6018E+01  0.7097E+01
   38   28    5    1.60   46.29 -0.1141E+02  0.9140E+01
   22    9   21    1.60   15.15 -0.3734E+01  0.8387E+01
   46    1    7    1.61   15.79 -0.3893E+01  0.8839E+01
   39   28    1    1.61   10.31 -0.2543E+01  0.9297E+01
   31    6   18    1.62   41.05 -0.1012E+02  0.1138E+02
   37    4   15    1.62   47.15 -0.1162E+02  0.8595E+01
   21    9   21    1.62   34.48 -0.8502E+01  0.7731E+01
   45    1    8    1.62   57.77 -0.1424E+02  0.9309E+01
   46    1    6    1.63   15.85 -0.3909E+01  0.1277E+02
   31   25   13    1.63   12.56 -0.3096E+01  0.8296E+01
   43    2   10    1.63   18.12 -0.1769E+02  0.1340E+02
   37   27    6    1.63   12.44 -0.1214E+02  0.1262E+02
   29   24   14    1.65   21.50 -0.2099E+02  0.1163E+02
   37    3   14    1.66   10.69 -0.1044E+02  0.7636E+01
   35   26    8    1.66   20.37 -0.1988E+02  0.1603E+02
   37   26    5    1.67   11.96 -0.1167E+02  0.8921E+01
   34   25    9    1.69   14.87 -0.1452E+02  0.8876E+01
   37   25    1    1.73   11.26  0.1027E+02  0.9404E+01
   23   21   16    1.73   12.44  0.1134E+02  0.1206E+02
   44    1    1    1.75   11.03  0.1006E+02  0.9414E+01
   26   21   14    1.77   17.93  0.1635E+02  0.7804E+01
   30   21   10    1.85   10.51 -0.1523E+02  0.1209E+02
   34   22    4    1.87   10.89 -0.1577E+02  0.1015E+02
   39    2    1    1.98   10.10 -0.1398E+02  0.9136E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        83
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1479
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1305
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      174


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
