
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.384    77.384    38.726  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      13      13
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      13
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000769  0.000692
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
 REFLECTING_RANGE_E.S.D.=   0.095
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

    1464 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      13 ...      13
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1464 REFLECTIONS ACCEPTED



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

 *  44        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  31        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  35        mP          0.0      77.4   38.7   77.4  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  33        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.4   77.4   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.5   38.7   77.4  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.5   38.7   77.4  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.7  159.5   77.4  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.5   77.4  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.7  159.5   77.4  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.7  159.5   77.4  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    27        mC        500.0     159.5   38.7  109.4  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.7  109.4  116.1  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.7  159.5  159.5  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.4  116.1   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        998.3     173.0   77.4   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.3      77.4  173.0   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.3      77.4  173.0   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.3      77.4   77.4   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.1  116.1  116.1  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.5   86.5  116.1 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.5   86.5  116.1 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.4   86.5   86.5  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.5   86.5  109.4  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.5   86.5  109.4  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.7   86.5  247.8 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.5   86.5   77.4  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.5   86.5   77.4  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.4  177.3  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.5   86.5  177.3 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.5   86.5   77.4  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.3  109.4   38.7  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.7   77.4  177.3 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.5  177.3   77.4 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    6.7      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.745   -23.457     0.024
 UNIT_CELL_B-AXIS=    23.457    73.745     0.013
 UNIT_CELL_C-AXIS=    -0.013    -0.003    38.725
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       304 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.23
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.340    77.340    38.704  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  4.5E-02 4.5E-02 2.3E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012930  0.012930  0.025837  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.705   -23.434     0.031
 COORDINATES OF UNIT CELL B-AXIS    23.434    73.705     0.011
 COORDINATES OF UNIT CELL C-AXIS    -0.016    -0.001    38.704
 CRYSTAL MOSAICITY (DEGREES)     0.095
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006416  0.002011
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.88   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1260.95   1242.22
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.19
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.34    77.34    38.70  90.000  90.000  90.000


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
  0.05  0.15       110        176.994       1.195
  0.15  0.25       168        218.406       0.968
  0.25  0.35       197        173.949       1.216
  0.35  0.45       147        152.971       1.382
  0.45  0.55       170        308.827       0.685
  0.55  0.65       198        194.564       1.087
  0.65  0.75       186        188.898       1.119
  0.75  0.85       170        180.338       1.172
  0.85  0.95       118        337.392       0.627
 
 -0.05  0.95      1464        211.438       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               241.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.987
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    12.1 XMAX=    12.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59625 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        518


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               241.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.983
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     7.6 XMAX=  2436.9 NXBIN=    3
 YMIN=     8.1 YMAX=  2517.8 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        518


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               241.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.973
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    12.1 XMAX=    12.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1222    1263
 DETECTOR_SURFACE_POSITION=    1635    1689
 DETECTOR_SURFACE_POSITION=     810    1689
 DETECTOR_SURFACE_POSITION=     810     837
 DETECTOR_SURFACE_POSITION=    1635     837
 DETECTOR_SURFACE_POSITION=    2156    1663
 DETECTOR_SURFACE_POSITION=    1609    2228
 DETECTOR_SURFACE_POSITION=     835    2228
 DETECTOR_SURFACE_POSITION=     289    1663
 DETECTOR_SURFACE_POSITION=     289     863
 DETECTOR_SURFACE_POSITION=     835     298
 DETECTOR_SURFACE_POSITION=    1609     298
 DETECTOR_SURFACE_POSITION=    2156     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        518


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
 1.101E+00  2.044E-03   21.08


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

   24.457   3.398    17.95   1.65      2.68      2.79      62      87      20
    3.398   2.414    13.20   1.50      4.14      3.55     126     167      14
    2.414   1.975     7.04   0.95      6.74      6.38     190     220       0
    1.975   1.711    -0.04   0.66    -99.90    -99.90     188     236       0
    1.711   1.531    -0.02   0.68    -99.90    -99.90     220     277       0
    1.531   1.398     0.03   0.69   2299.48   2690.42     240     298       0
    1.398   1.295     0.05   0.51    823.06   1169.42      94     142       2
   --------------------------------------------------------------------------
   24.457   1.295     3.73   0.86      5.80      6.15    1120    1427      36
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      13      13)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681  -0.2054E+01  0.1200E+02        -0.17           1 
   47    0    0     1.646  -0.1350E+02  0.1209E+02        -1.12           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -6.1% OF MEAN INTENSITY


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

     3.84          64      42      1285        3.3%       2.9%      2.8%       44   22.39      4.1%    99.8*     0    0.000       0
     2.73         104      66      2098        3.1%       3.0%      3.2%       76   19.18      4.2%    99.6*     0    0.000       0
     2.24         166      99      2665        3.7%       5.1%      4.2%      134   14.08      7.2%    99.4*     0    0.000       0
     1.94         165      93      3101        3.0%       8.6%      9.3%      144    5.74     12.1%    98.7*     0    0.000       0
     1.74         181     108      3479        3.1%     -99.9%    -99.9%      146    0.00    -99.9%     2.0      0    0.000       0
     1.58         214     128      3850        3.3%     776.9%   1020.3%      172    0.02   1098.7%   -10.1      0    0.000       0
     1.47         233     142      4143        3.4%    1046.3%   1292.4%      182    0.04   1479.7%     3.1      0    0.000       0
     1.37         217     128      4428        2.9%     -99.9%    -99.9%      178    0.00    -99.9%   -16.0      0    0.000       0
     1.30          83      61      4739        1.3%     419.5%    579.6%       44    0.11    593.3%    23.2      0    0.000       0
    total        1427     867     29788        2.9%       5.8%      6.2%     1120    4.78      8.2%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1464
 NUMBER OF REJECTED MISFITS                                36
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1427
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    867



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

   WILSON LINE (using all data) : A=   8.375 B=  33.709 CORRELATION=  0.97
      #      RES      SS        <I>       log(<I>)       BO
      17     7.253   0.005  2.1326E+03       7.665      74.7
      21     4.622   0.012  2.6815E+03       7.894      20.6
      30     3.496   0.020  1.1501E+03       7.048      32.4
      30     2.989   0.028  5.6338E+02       6.334      36.5
      41     2.645   0.036  5.1863E+02       6.251      29.7
      54     2.395   0.044  3.0313E+02       5.714      30.5
      41     2.204   0.051  1.7811E+02       5.182      31.0
      45     2.056   0.059  7.3793E+01       4.301      34.4
      46     1.931   0.067 -4.7705E-01    -999.999   -1000.0
      59     1.828   0.075 -3.2759E-01    -999.999   -1000.0
      52     1.739   0.083 -4.7191E-01    -999.999   -1000.0
      56     1.665   0.090 -3.6666E-01    -999.999   -1000.0
      59     1.596   0.098 -5.0388E-01    -999.999   -1000.0
      73     1.536   0.106  1.7620E+00       0.566      36.8
      76     1.480   0.114  8.0251E-01      -0.220      37.7
      58     1.432   0.122 -7.5512E-01    -999.999   -1000.0
      50     1.390   0.129  2.9208E-01      -1.231      37.1
      41     1.351   0.137  2.7682E+00       1.018      26.9
      18     1.314   0.145 -1.3289E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.253        1.839        1.206        0.625
       3     4.622        1.220        0.728        0.324
       3     3.496        0.341        0.084        0.015
       4     2.989        0.955        0.396        0.125
       3     2.645        0.258        0.058        0.010
       5     2.395        0.739        0.335        0.117
       2     2.204        1.365        0.557        0.163
       8     2.056        0.295        0.112        0.033
       7     1.931       29.394       38.440      146.028
      11     1.828       23.935       61.109      183.738
       9     1.739       19.542      -23.249       92.315
       6     1.480       36.396      108.023      240.560
       2     1.351        1.936        1.274        0.610
      88   overall        2.941        7.534       16.511


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.253        0.441        0.301        0.166
      18     4.622        1.089        1.218        1.156
      27     3.496        0.860        0.680        0.504
      26     2.989        0.900        0.830        0.690
      38     2.645        1.153        1.253        1.231
      49     2.395        0.943        0.783        0.545
      39     2.204        0.835        0.653        0.443
      37     2.056        1.479        1.992        2.419
      56     1.596       34.479       11.008      843.855
      70     1.536       10.929       15.466       64.375
      70     1.480       44.830       38.889      900.440
      50     1.390      320.384      543.500    39260.861
      39     1.351        4.863        4.356       11.083
     779   overall       26.137       40.309     2607.506

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.253  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.32  0.30  0.29
       3     4.622  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.496  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.989  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.645  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       5     2.395  0.00  0.00  0.00  0.00  0.38  1.07  1.00  0.95  0.91  0.88
       2     2.204  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.056  2.01  1.45  1.50  1.32  1.20  1.11  1.05  0.99  1.14  1.10
       7     1.931  2.88  2.07  1.71  1.51  1.37  1.27  1.20  1.14  1.09  1.05
      11     1.828  2.20  1.58  1.31  1.15  1.05  0.97  0.91  0.87  0.83  0.80
       9     1.739  3.13  2.25  1.87  1.64  1.49  1.39  1.30  1.24  1.18  1.14
       6     1.480  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.27  0.25  0.24
       2     1.351  0.00  0.00  0.00  0.00  0.00  0.89  0.84  0.80  0.76  0.73
      88   overall  1.14  0.82  0.71  0.62  0.59  0.63  0.59  0.60  0.61  0.60


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.253  2.63  2.76  1.93  1.52  1.91  1.85  1.66  1.51  1.40  1.32
      18     4.622  1.17  1.84  1.29  1.01  1.13  1.23  1.10  1.11  1.22  1.14
      27     3.496  0.78  0.61  0.57  0.67  0.66  0.66  0.88  0.94  0.94  0.88
      26     2.989  1.21  1.49  1.19  1.17  1.17  1.28  1.22  1.12  1.10  1.03
      38     2.645  0.83  1.16  1.42  1.20  1.07  0.99  1.05  1.10  1.06  1.04
      49     2.395  1.50  1.13  1.10  0.99  1.09  1.13  1.09  1.07  1.10  1.03
      39     2.204  1.62  1.27  0.89  1.01  0.98  0.91  0.97  1.02  1.04  1.01
      37     2.056  1.42  0.89  1.04  1.23  1.03  1.02  1.02  0.93  1.00  1.07
      56     1.596  5.25  2.86  2.00  1.57  1.32  1.15  1.06  0.97  0.90  0.85
      70     1.536  4.80  2.52  1.76  1.39  1.20  1.04  0.94  0.88  0.82  0.77
      70     1.480  4.95  2.60  1.82  1.43  1.23  1.08  0.96  0.88  0.82  0.77
      50     1.390  5.04  2.65  1.85  1.46  1.22  1.11  0.99  0.94  0.88  0.82
      39     1.351  4.58  2.40  1.68  1.40  1.17  1.08  1.02  1.02  0.95  0.89
     779   overall  2.23  1.35  1.02  0.87  0.79  0.73  0.70  0.68  0.66  0.63


 List of    116 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   30    3   15    1.82   83.19 -0.2725E+02  0.1097E+02 "alien"
   43    1    8    1.69   63.94 -0.2344E+02  0.1383E+02 "alien"
   30    1   13    1.95   49.32 -0.2353E+02  0.1154E+02 "alien"
   52    4    9    1.40   49.28  0.1439E+02  0.9988E+01 "alien"
   46    6   15    1.40   48.75  0.1424E+02  0.9530E+01 "alien"
   36   28    9    1.58   47.88 -0.2413E+02  0.1705E+02 "alien"
   54    4    7    1.38   44.62  0.1303E+02  0.9984E+01 "alien"
   29   12   23    1.39   44.26  0.1293E+02  0.1041E+02 "alien"
   28    2   14    1.95   41.23 -0.1967E+02  0.1094E+02 "alien"
   50    5   11    1.41   40.12  0.1172E+02  0.1848E+02 "alien"
   43    1    1    1.80   36.75 -0.1204E+02  0.8637E+01 "alien"
   38    9   20    1.38   36.27  0.1059E+02  0.9785E+01 "alien"
   33   23    8    1.79   34.80 -0.1140E+02  0.8365E+01 "alien"
   36   26    6    1.68   34.19 -0.1253E+02  0.1404E+02 "alien"
   39    2    3    1.96   33.59 -0.1603E+02  0.1046E+02 "alien"
   35   24    6    1.75   33.55 -0.1583E+02  0.8771E+01 "alien"
   37   31   12    1.43   33.10 -0.2500E+02  0.1498E+02 "alien"
   54    4    6    1.39   32.10  0.9375E+01  0.1004E+02 "alien"
   45    1    7    1.64   30.94 -0.1135E+02  0.8489E+01 "alien"
   37   26    5    1.67   30.52 -0.1119E+02  0.8843E+01 "alien"
   33   11   22    1.38   29.28  0.8554E+01  0.8665E+01 "alien"
   42    1    9    1.69   28.75 -0.1054E+02  0.8377E+01 "alien"
   46    1    6    1.63   28.34 -0.1428E+02  0.1251E+02 "alien"
   38    2   12    1.72   27.62 -0.1304E+02  0.1317E+02 "alien"
   29   20    9    1.96   27.28 -0.1301E+02  0.1351E+02 "alien"
   34   23    5    1.83   26.84 -0.8794E+01  0.8861E+01 "alien"
   17   11   20    1.73   26.83 -0.1266E+02  0.6024E+01 "alien"
   21   17   20    1.60   26.69 -0.1345E+02  0.8132E+01 "alien"
   35    1   12    1.82   26.18 -0.8575E+01  0.8149E+01 "alien"
   25   19   13    1.90   25.72 -0.1227E+02  0.7854E+01 "alien"
   36   30   12    1.47   25.46  0.2043E+02  0.9942E+01 "alien"
   39    2    1    1.98   25.33 -0.1208E+02  0.8970E+01 "alien"
   36   33   14    1.37   25.27  0.7381E+01  0.9939E+01 "alien"
   43    2   11    1.60   24.84 -0.1252E+02  0.9497E+01 "alien"
   35   27   10    1.59   24.73 -0.1246E+02  0.8664E+01 "alien"
   34   26   10    1.64   23.97 -0.8789E+01  0.8588E+01 "alien"
   31   22    9    1.84   23.79 -0.7793E+01  0.9640E+01 "alien"
   29   22   12    1.77   23.78 -0.1122E+02  0.7972E+01 "alien"
   49    6   13    1.39   23.71  0.6924E+01  0.1718E+02 "alien"
   51    5   11    1.39   23.01  0.6720E+01  0.9846E+01 "alien"
   44    1    8    1.65   22.98 -0.8425E+01  0.8391E+01 "alien"
   35    2   13    1.77   22.69 -0.1071E+02  0.7906E+01 "alien"
   32   11   22    1.39   21.70  0.6339E+01  0.1300E+02 "alien"
   38   28    5    1.60   20.79 -0.1048E+02  0.8925E+01 "alien"
   26   23   20    1.46   20.10  0.1613E+02  0.1236E+02 "alien"
   36   32   13    1.41   20.09 -0.1517E+02  0.1626E+02 "alien"
   41   32    4    1.47   19.11  0.1534E+02  0.1061E+02 "alien"
   52    5   10    1.38   19.10  0.5579E+01  0.9769E+01 "alien"
   51    8   15    1.30   18.80 -0.2498E+02  0.1368E+02 "alien"
   21   20   18    1.67   18.76 -0.6879E+01  0.1189E+02 "alien"
   25   24   17    1.59   18.44 -0.9291E+01  0.8462E+01 "alien"
   35   24    4    1.79   18.41 -0.6030E+01  0.1618E+02 "alien"
   47    5   13    1.43   18.38 -0.1388E+02  0.9102E+01 "alien"
   40   33    9    1.41   18.28  0.5338E+01  0.1036E+02 "alien"
   51    4    9    1.43   18.16 -0.1372E+02  0.9397E+01 "alien"
   54    3    3    1.42   17.99 -0.1359E+02  0.9840E+01 "alien"
   39   33   10    1.41   17.92  0.5235E+01  0.1095E+02 "alien"
   14   10   18    1.94   17.87 -0.8524E+01  0.6211E+01 "alien"
   42    1    1    1.84   17.78 -0.5826E+01  0.9438E+01 "alien"
   34   23    6    1.81   17.59 -0.5762E+01  0.8734E+01 "alien"
   52    2    4    1.47   17.42  0.1398E+02  0.1664E+02 "alien"
   36   24    1    1.79   16.98 -0.5563E+01  0.8793E+01 "alien"
   19   18   18    1.74   16.95 -0.7998E+01  0.6570E+01 "alien"
   45    9   18    1.33   16.82 -0.2235E+02  0.9657E+01 "alien"
   28   24   15    1.63   16.78 -0.8454E+01  0.8199E+01 "alien"
   38   26    0    1.68   16.54 -0.1925E+02  0.1693E+02 "alien"
   19    9   22    1.59   16.26 -0.8192E+01  0.6335E+01 "alien"
   46    0    3    1.67   15.91 -0.1766E+02  0.9051E+01 "alien"
   46    7   16    1.37   15.78  0.4610E+01  0.9519E+01 "alien"
   18   10   21    1.65   15.78 -0.5784E+01  0.6340E+01 "alien"
   36   24    3    1.77   15.70 -0.7411E+01  0.9633E+01 "alien"
   31   25   13    1.63   15.60 -0.7862E+01  0.8260E+01 "alien"
   40   29    0    1.57   15.46 -0.2282E+02  0.1301E+02 "alien"
   51    4   10    1.41   15.28  0.4462E+01  0.1500E+02 "alien"
   31    2   14    1.85   14.95 -0.4896E+01  0.7062E+01 "alien"
   36   25    5    1.72   14.85 -0.7008E+01  0.8649E+01 "alien"
   39   31    8    1.48   14.70  0.1180E+02  0.9747E+01 "alien"
   40   31    5    1.50   14.61  0.1173E+02  0.1207E+02 "alien"
   19   16   19    1.70   14.52 -0.6852E+01  0.5719E+01 "alien"
   39    1   11    1.73   13.71 -0.6471E+01  0.1232E+02 "alien"
   27   22   21    1.42   13.68 -0.1033E+02  0.8083E+01 "alien"
   42    0    6    1.77   13.64 -0.1654E+02  0.1558E+02 "alien"
   15    7   19    1.87   13.52 -0.4428E+01  0.5863E+01 "alien"
   39    1   10    1.76   13.32 -0.6287E+01  0.1337E+02 "alien"
   20   19   18    1.71   13.15 -0.6204E+01  0.6640E+01 "alien"
   37   26    4    1.68   12.97 -0.4757E+01  0.8855E+01 "alien"
   39   29    4    1.57   12.82 -0.6461E+01  0.9625E+01 "alien"
   51    3    8    1.44   12.51 -0.9446E+01  0.1064E+02 "alien"
   31   22   10    1.80   12.49 -0.4091E+01  0.1286E+02 "alien"
   19    9   21    1.65   12.46 -0.4569E+01  0.6098E+01 "alien"
   44    4   13    1.51   12.39  0.2183E+02  0.9182E+01 "alien"
   32   24   11    1.69   12.35 -0.4528E+01  0.8551E+01 "alien"
   35   23    1    1.84   12.33 -0.4038E+01  0.8897E+01 "alien"
   53    3    3    1.45   12.21 -0.9219E+01  0.1201E+02 "alien"
   33   25   10    1.68   12.06 -0.4423E+01  0.1021E+02 "alien"
   25   13   24    1.39   11.94  0.3488E+01  0.5952E+01 "alien"
   28   28   18    1.45   11.94 -0.2061E+02  0.1213E+02 "alien"
   20   18   16    1.85   11.93 -0.3907E+01  0.9011E+01 "alien"
   22   18   20    1.58   11.91 -0.5999E+01  0.6439E+01 "alien"
   47    2    9    1.54   11.74  0.2068E+02  0.1532E+02 "alien"
   49    3    9    1.48   11.70  0.9391E+01  0.9376E+01 "alien"
   33    0   11    1.95   11.55 -0.1229E+02  0.1308E+02 "alien"
   32   21    6    1.93   11.42 -0.5446E+01  0.9168E+01 "alien"
   39    0    8    1.83   11.41 -0.8258E+01  0.9888E+01 "alien"
   54    3    0    1.43   11.40 -0.1900E+02  0.1453E+02 "alien"
   29   24   14    1.65   11.22 -0.4115E+01  0.8398E+01 "alien"
   27   24   16    1.60   11.22 -0.5653E+01  0.1293E+02 "alien"
   35   23    2    1.84   10.83 -0.3548E+01  0.8844E+01 "alien"
   50    2    4    1.53   10.72  0.1888E+02  0.1732E+02 "alien"
   25   21   15    1.74   10.62 -0.5012E+01  0.1036E+02 "alien"
   41    1    4    1.85   10.42 -0.3415E+01  0.8808E+01 "alien"
   34   22    4    1.87   10.32 -0.3380E+01  0.1004E+02 "alien"
   48    3    9    1.51   10.17  0.8162E+01  0.1014E+02 "alien"
   43    7   17    1.40   10.15  0.2966E+01  0.9266E+01 "alien"
   48    5   13    1.41   10.14 -0.7659E+01  0.9805E+01 "alien"
   37   30   10    1.50   10.09  0.8101E+01  0.9769E+01 "alien"


 List of    116 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   51    8   15    1.30   18.80 -0.2498E+02  0.1368E+02
   45    9   18    1.33   16.82 -0.2235E+02  0.9657E+01
   46    7   16    1.37   15.78  0.4610E+01  0.9519E+01
   36   33   14    1.37   25.27  0.7381E+01  0.9939E+01
   33   11   22    1.38   29.28  0.8554E+01  0.8665E+01
   52    5   10    1.38   19.10  0.5579E+01  0.9769E+01
   54    4    7    1.38   44.62  0.1303E+02  0.9984E+01
   38    9   20    1.38   36.27  0.1059E+02  0.9785E+01
   49    6   13    1.39   23.71  0.6924E+01  0.1718E+02
   51    5   11    1.39   23.01  0.6720E+01  0.9846E+01
   29   12   23    1.39   44.26  0.1293E+02  0.1041E+02
   25   13   24    1.39   11.94  0.3488E+01  0.5952E+01
   32   11   22    1.39   21.70  0.6339E+01  0.1300E+02
   54    4    6    1.39   32.10  0.9375E+01  0.1004E+02
   43    7   17    1.40   10.15  0.2966E+01  0.9266E+01
   46    6   15    1.40   48.75  0.1424E+02  0.9530E+01
   52    4    9    1.40   49.28  0.1439E+02  0.9988E+01
   51    4   10    1.41   15.28  0.4462E+01  0.1500E+02
   40   33    9    1.41   18.28  0.5338E+01  0.1036E+02
   39   33   10    1.41   17.92  0.5235E+01  0.1095E+02
   50    5   11    1.41   40.12  0.1172E+02  0.1848E+02
   48    5   13    1.41   10.14 -0.7659E+01  0.9805E+01
   36   32   13    1.41   20.09 -0.1517E+02  0.1626E+02
   27   22   21    1.42   13.68 -0.1033E+02  0.8083E+01
   54    3    3    1.42   17.99 -0.1359E+02  0.9840E+01
   51    4    9    1.43   18.16 -0.1372E+02  0.9397E+01
   54    3    0    1.43   11.40 -0.1900E+02  0.1453E+02
   47    5   13    1.43   18.38 -0.1388E+02  0.9102E+01
   37   31   12    1.43   33.10 -0.2500E+02  0.1498E+02
   51    3    8    1.44   12.51 -0.9446E+01  0.1064E+02
   28   28   18    1.45   11.94 -0.2061E+02  0.1213E+02
   53    3    3    1.45   12.21 -0.9219E+01  0.1201E+02
   26   23   20    1.46   20.10  0.1613E+02  0.1236E+02
   52    2    4    1.47   17.42  0.1398E+02  0.1664E+02
   36   30   12    1.47   25.46  0.2043E+02  0.9942E+01
   41   32    4    1.47   19.11  0.1534E+02  0.1061E+02
   39   31    8    1.48   14.70  0.1180E+02  0.9747E+01
   49    3    9    1.48   11.70  0.9391E+01  0.9376E+01
   37   30   10    1.50   10.09  0.8101E+01  0.9769E+01
   40   31    5    1.50   14.61  0.1173E+02  0.1207E+02
   48    3    9    1.51   10.17  0.8162E+01  0.1014E+02
   44    4   13    1.51   12.39  0.2183E+02  0.9182E+01
   50    2    4    1.53   10.72  0.1888E+02  0.1732E+02
   47    2    9    1.54   11.74  0.2068E+02  0.1532E+02
   40   29    0    1.57   15.46 -0.2282E+02  0.1301E+02
   39   29    4    1.57   12.82 -0.6461E+01  0.9625E+01
   22   18   20    1.58   11.91 -0.5999E+01  0.6439E+01
   36   28    9    1.58   47.88 -0.2413E+02  0.1705E+02
   19    9   22    1.59   16.26 -0.8192E+01  0.6335E+01
   25   24   17    1.59   18.44 -0.9291E+01  0.8462E+01
   35   27   10    1.59   24.73 -0.1246E+02  0.8664E+01
   43    2   11    1.60   24.84 -0.1252E+02  0.9497E+01
   38   28    5    1.60   20.79 -0.1048E+02  0.8925E+01
   27   24   16    1.60   11.22 -0.5653E+01  0.1293E+02
   21   17   20    1.60   26.69 -0.1345E+02  0.8132E+01
   31   25   13    1.63   15.60 -0.7862E+01  0.8260E+01
   46    1    6    1.63   28.34 -0.1428E+02  0.1251E+02
   28   24   15    1.63   16.78 -0.8454E+01  0.8199E+01
   34   26   10    1.64   23.97 -0.8789E+01  0.8588E+01
   45    1    7    1.64   30.94 -0.1135E+02  0.8489E+01
   19    9   21    1.65   12.46 -0.4569E+01  0.6098E+01
   29   24   14    1.65   11.22 -0.4115E+01  0.8398E+01
   44    1    8    1.65   22.98 -0.8425E+01  0.8391E+01
   18   10   21    1.65   15.78 -0.5784E+01  0.6340E+01
   46    0    3    1.67   15.91 -0.1766E+02  0.9051E+01
   37   26    5    1.67   30.52 -0.1119E+02  0.8843E+01
   21   20   18    1.67   18.76 -0.6879E+01  0.1189E+02
   38   26    0    1.68   16.54 -0.1925E+02  0.1693E+02
   36   26    6    1.68   34.19 -0.1253E+02  0.1404E+02
   33   25   10    1.68   12.06 -0.4423E+01  0.1021E+02
   37   26    4    1.68   12.97 -0.4757E+01  0.8855E+01
   43    1    8    1.69   63.94 -0.2344E+02  0.1383E+02
   42    1    9    1.69   28.75 -0.1054E+02  0.8377E+01
   32   24   11    1.69   12.35 -0.4528E+01  0.8551E+01
   19   16   19    1.70   14.52 -0.6852E+01  0.5719E+01
   20   19   18    1.71   13.15 -0.6204E+01  0.6640E+01
   38    2   12    1.72   27.62 -0.1304E+02  0.1317E+02
   36   25    5    1.72   14.85 -0.7008E+01  0.8649E+01
   17   11   20    1.73   26.83 -0.1266E+02  0.6024E+01
   39    1   11    1.73   13.71 -0.6471E+01  0.1232E+02
   19   18   18    1.74   16.95 -0.7998E+01  0.6570E+01
   25   21   15    1.74   10.62 -0.5012E+01  0.1036E+02
   35   24    6    1.75   33.55 -0.1583E+02  0.8771E+01
   39    1   10    1.76   13.32 -0.6287E+01  0.1337E+02
   36   24    3    1.77   15.70 -0.7411E+01  0.9633E+01
   42    0    6    1.77   13.64 -0.1654E+02  0.1558E+02
   35    2   13    1.77   22.69 -0.1071E+02  0.7906E+01
   29   22   12    1.77   23.78 -0.1122E+02  0.7972E+01
   36   24    1    1.79   16.98 -0.5563E+01  0.8793E+01
   33   23    8    1.79   34.80 -0.1140E+02  0.8365E+01
   35   24    4    1.79   18.41 -0.6030E+01  0.1618E+02
   43    1    1    1.80   36.75 -0.1204E+02  0.8637E+01
   31   22   10    1.80   12.49 -0.4091E+01  0.1286E+02
   34   23    6    1.81   17.59 -0.5762E+01  0.8734E+01
   30    3   15    1.82   83.19 -0.2725E+02  0.1097E+02
   35    1   12    1.82   26.18 -0.8575E+01  0.8149E+01
   34   23    5    1.83   26.84 -0.8794E+01  0.8861E+01
   39    0    8    1.83   11.41 -0.8258E+01  0.9888E+01
   35   23    2    1.84   10.83 -0.3548E+01  0.8844E+01
   42    1    1    1.84   17.78 -0.5826E+01  0.9438E+01
   31   22    9    1.84   23.79 -0.7793E+01  0.9640E+01
   35   23    1    1.84   12.33 -0.4038E+01  0.8897E+01
   31    2   14    1.85   14.95 -0.4896E+01  0.7062E+01
   20   18   16    1.85   11.93 -0.3907E+01  0.9011E+01
   41    1    4    1.85   10.42 -0.3415E+01  0.8808E+01
   15    7   19    1.87   13.52 -0.4428E+01  0.5863E+01
   34   22    4    1.87   10.32 -0.3380E+01  0.1004E+02
   25   19   13    1.90   25.72 -0.1227E+02  0.7854E+01
   32   21    6    1.93   11.42 -0.5446E+01  0.9168E+01
   14   10   18    1.94   17.87 -0.8524E+01  0.6211E+01
   30    1   13    1.95   49.32 -0.2353E+02  0.1154E+02
   33    0   11    1.95   11.55 -0.1229E+02  0.1308E+02
   28    2   14    1.95   41.23 -0.1967E+02  0.1094E+02
   29   20    9    1.96   27.28 -0.1301E+02  0.1351E+02
   39    2    3    1.96   33.59 -0.1603E+02  0.1046E+02
   39    2    1    1.98   25.33 -0.1208E+02  0.8970E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        46
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1464
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1354
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      110


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
