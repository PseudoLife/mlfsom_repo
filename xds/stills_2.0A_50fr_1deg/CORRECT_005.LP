
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.305    77.305    38.689  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=       5       5
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       5
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  1.000000 -0.000631  0.000704
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
 BEAM_DIVERGENCE_E.S.D.=   0.040
 REFLECTING_RANGE_E.S.D.=   0.082
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

    1518 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       5 ...       5
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1518 REFLECTIONS ACCEPTED



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
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas               48
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
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  1  0  0  0  0 -1  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  20        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    27        mC        500.0     159.4   38.7  109.3  90.0 133.3  76.0   -2  0  1  0  0  0  1  0  1  1  0  0
    19        oI        500.0      38.7  109.3  116.0  90.0  70.5  90.0    0  0  1  0  1  1  0  0 -1  1  1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0 -1  0 -2  0  1  0  0  2  1  0
    18        tI        625.0     109.3  116.0   38.7  70.5  90.0  90.0    1  1  0  0  1 -1 -1  0  0  0 -1  0
    30        mC        998.1      77.3  172.9   38.7  90.0  90.0  63.4   -1  0  0  0 -1 -2  0  0  0  0  1  0
    40        oC        998.1      77.3  172.9   38.7  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    41        mC        998.1     172.9   77.3   38.7  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    22        hP        998.1      77.3   77.3   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0    1  1 -1  0 -1 -1 -1  0 -1  1  1  0
     2        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.6    1  0  1  0  0 -1 -1  0  1  1 -1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.3   86.4   86.4  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.4   86.4  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.4  247.5 102.1  99.0 116.6    0  0 -1  0 -1  0  1  0  1  3  1  0
    10        mC        999.0      86.4   86.4   77.3  90.0  90.0 126.8   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.7   77.3  177.1  64.1  77.4  90.0    0  0  1  0  1  0  0  0  1  2  1  0
    16        oF        999.0      86.4   86.4  177.1 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    17        mC        999.0      86.4   86.4   86.4  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.1  109.3   38.7  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.7   77.3  177.1 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.4  177.1   77.3 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    5.8      19    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.681   -23.387     0.045
 UNIT_CELL_B-AXIS=    23.387    73.681     0.014
 UNIT_CELL_C-AXIS=    -0.023     0.000    38.688
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       344 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.19
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.07
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.266    77.266    38.667  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  8.0E-02 8.0E-02 4.5E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012942  0.012942  0.025862  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.640   -23.391     0.020
 COORDINATES OF UNIT CELL B-AXIS    23.391    73.640     0.005
 COORDINATES OF UNIT CELL C-AXIS    -0.010     0.001    38.667
 CRYSTAL MOSAICITY (DEGREES)     0.082
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999977  0.006492  0.002016
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.86   1232.90
 DETECTOR ORIGIN (PIXELS) AT                     1260.92   1242.34
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.14
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.27    77.27    38.67  90.000  90.000  90.000


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
  0.05  0.15       126        306.885       0.970
  0.15  0.25       185        259.987       1.145
  0.25  0.35       183        268.612       1.108
  0.35  0.45       154        157.567       1.889
  0.45  0.55       183        435.448       0.684
  0.55  0.65       199        252.375       1.180
  0.65  0.75       179        224.160       1.328
  0.75  0.85       169        267.087       1.115
  0.85  0.95       140        546.857       0.544
 
 -0.05  0.95      1518        297.706       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               238.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.954
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     4.1 XMAX=     4.9 NXBIN=    1
 YMIN= 0.00168 YMAX= 0.59740 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        528


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               238.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.948
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    12.3 XMAX=  2453.9 NXBIN=    3
 YMIN=     6.4 YMAX=  2519.7 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        528


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               238.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.943
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     4.1 XMAX=     4.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1233    1263
 DETECTOR_SURFACE_POSITION=    1648    1690
 DETECTOR_SURFACE_POSITION=     818    1690
 DETECTOR_SURFACE_POSITION=     818     836
 DETECTOR_SURFACE_POSITION=    1648     836
 DETECTOR_SURFACE_POSITION=    2172    1663
 DETECTOR_SURFACE_POSITION=    1622    2229
 DETECTOR_SURFACE_POSITION=     844    2229
 DETECTOR_SURFACE_POSITION=     295    1663
 DETECTOR_SURFACE_POSITION=     295     863
 DETECTOR_SURFACE_POSITION=     844     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2172     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        528


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
 1.245E+00  2.233E-03   18.96


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

   24.434   3.395    17.16   1.11      2.30      3.05      56      84      26
    3.395   2.412    13.70   1.53      4.67      3.48     128     170      14
    2.412   1.973     9.24   0.98      5.18      4.71     176     217      12
    1.973   1.710    -0.03   0.62    -99.90    -99.90     192     246       0
    1.710   1.530     0.02   0.63    841.67   1011.31     212     284       0
    1.530   1.397    -0.08   0.58    -99.90    -99.90     244     315       0
    1.397   1.294     0.03   0.71    -99.90    -99.90     100     147       2
   --------------------------------------------------------------------------
   24.434   1.294     3.93   0.81      5.23      5.44    1108    1463      54
 


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

   46    0    0     1.680  -0.8336E+01  0.1237E+02        -0.67           1 
   47    0    0     1.644  -0.1242E+02  0.1267E+02        -0.98           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -4.0% OF MEAN INTENSITY


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

     3.84          63      42      1285        3.3%       2.4%      3.0%       42   21.07      3.3%    99.9*     0    0.000       0
     2.73         104      67      2098        3.2%       4.1%      3.3%       74   18.44      5.7%    98.9*     0    0.000       0
     2.23         162      99      2665        3.7%       4.4%      3.8%      126   15.73      6.2%    99.4*     0    0.000       0
     1.94         167      98      3101        3.2%       6.3%      5.7%      138    8.52      9.0%    99.2*     0    0.000       0
     1.73         190     116      3479        3.3%    8695.6%  10381.9%      148   -0.06  12297.4%    10.5      0    0.000       0
     1.58         218     134      3850        3.5%     682.3%    845.9%      168    0.03    964.9%     4.8      0    0.000       0
     1.47         250     157      4143        3.8%     -99.9%    -99.9%      186    0.00    -99.9%    -4.9      0    0.000       0
     1.37         221     134      4428        3.0%     -99.9%    -99.9%      174    0.00    -99.9%    -1.5      0    0.000       0
     1.29          87      61      4739        1.3%     -99.9%    -99.9%       52    0.00    -99.9%   -19.5      0    0.000       0
    total        1463     909     29788        3.1%       5.2%      5.4%     1108    4.95      7.4%    99.8*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1518
 NUMBER OF REJECTED MISFITS                                54
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1463
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    909



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

   WILSON LINE (using all data) : A=   8.547 B=  37.077 CORRELATION=  0.88
      #      RES      SS        <I>       log(<I>)       BO
      17     7.549   0.004  3.2485E+03       8.086      52.5
      21     4.661   0.012  3.3968E+03       8.131      18.1
      30     3.493   0.020  1.4788E+03       7.299      30.4
      31     2.988   0.028  7.7094E+02       6.648      33.9
      40     2.645   0.036  8.4324E+02       6.737      25.3
      54     2.395   0.044  5.1698E+02       6.248      26.4
      41     2.204   0.051  3.7968E+02       5.939      25.3
      49     2.055   0.059  1.8977E+02       5.246      27.9
      50     1.930   0.067 -5.1057E-01    -999.999   -1000.0
      63     1.825   0.075 -6.2238E-01    -999.999   -1000.0
      54     1.738   0.083  9.8817E-01      -0.012      51.7
      60     1.663   0.090  2.7928E-01      -1.276      54.3
      62     1.595   0.098  1.4649E-01      -1.921      53.2
      81     1.534   0.106  8.4350E-01      -0.170      41.0
      82     1.479   0.114 -1.7339E+00    -999.999   -1000.0
      59     1.431   0.122 -6.0984E-01    -999.999   -1000.0
      55     1.390   0.129 -5.2007E-01    -999.999   -1000.0
      40     1.350   0.137  1.2846E+00       0.250      30.2
      20     1.313   0.145  1.8946E+00       0.639      27.2


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.549        0.877        0.393        0.139
       3     4.661        0.944        0.494        0.193
       3     3.493        0.367        0.094        0.017
       4     2.988        1.016        0.434        0.141
       3     2.645        0.236        0.055        0.010
       5     2.395        0.810        0.375        0.130
       3     2.204        1.206        0.464        0.128
       8     2.055        0.246        0.079        0.020
      11     1.825       14.766       10.669       49.583
      11     1.738       32.362       99.658      297.718
       4     1.595      343.705     1846.034    18061.393
       3     1.534        6.041        5.540        3.726
       6     1.479       10.894        8.917       29.084
       1     1.390       60.528      163.127      314.028
       2     1.350       28.740       52.284      102.266
      93   overall       34.626       81.652     1020.686


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.549        1.330        2.026        2.527
      18     4.661        0.988        0.982        0.849
      27     3.493        0.837        0.618        0.418
      27     2.988        0.859        0.741        0.574
      37     2.645        1.131        1.188        1.118
      49     2.395        0.919        0.745        0.506
      38     2.204        0.792        0.593        0.381
      41     2.055        1.282        1.437        1.476
      43     1.738       30.952       39.734      446.510
      49     1.663      378.648      648.764    76012.847
      58     1.595     1743.246    10319.399  1756852.539
      78     1.534       43.970      -89.887     2172.153
      38     1.350       29.330       28.290      308.926
      20     1.313       29.913      127.181      851.691
     816   overall      154.883      770.675   129705.504

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.549  0.00  0.00  0.00  0.00  0.38  0.36  0.33  0.64  0.61  0.59
       3     4.661  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       3     3.493  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.988  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.645  1.34  0.96  0.80  0.70  0.64  0.59  1.12  1.06  1.01  0.98
       5     2.395  0.00  0.00  0.00  0.00  0.77  0.71  1.00  0.95  0.91  0.88
       3     2.204  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.055  1.51  1.45  1.50  1.32  1.20  1.11  1.05  0.99  1.14  1.10
      11     1.825  2.93  2.10  1.75  1.54  1.40  1.30  1.37  1.30  1.25  1.20
      11     1.738  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.72  0.69  0.67
       4     1.595  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
       3     1.534  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       6     1.479  2.68  1.93  1.60  1.41  1.28  1.19  1.40  1.33  1.27  1.22
       1     1.390  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       2     1.350  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
      93   overall  1.17  0.87  0.75  0.66  0.66  0.63  0.67  0.65  0.65  0.63


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.549  1.75  2.30  1.61  2.27  1.91  1.66  1.49  1.36  1.26  1.19
      18     4.661  1.17  1.23  0.86  0.84  1.13  1.11  0.99  1.11  1.12  1.05
      27     3.493  0.78  0.61  0.57  0.67  0.66  0.74  0.88  0.94  0.94  0.88
      27     2.988  1.56  1.63  1.14  1.24  1.04  1.07  1.10  1.14  1.06  1.05
      37     2.645  0.85  1.34  1.25  1.15  1.03  0.90  1.07  1.08  1.05  1.03
      49     2.395  1.50  1.24  1.02  1.18  1.09  1.04  1.09  1.07  1.10  1.03
      38     2.204  0.83  1.16  1.02  0.80  0.94  0.93  0.94  1.00  1.11  1.08
      41     2.055  1.79  1.21  1.22  1.11  0.99  0.86  0.92  0.97  0.90  0.93
      43     1.738  5.38  2.95  2.06  1.62  1.36  1.19  1.06  1.01  0.94  0.92
      49     1.663  4.93  2.59  1.89  1.49  1.24  1.09  0.97  0.89  0.83  0.77
      58     1.595  5.07  2.66  1.86  1.46  1.23  1.07  0.96  0.88  0.81  0.76
      78     1.534  4.58  2.40  1.83  1.48  1.24  1.08  0.97  0.91  0.84  0.81
      38     1.350  4.98  2.61  1.83  1.44  1.20  1.05  0.94  0.86  0.80  0.75
      20     1.313  5.25  2.76  2.12  1.82  1.52  1.33  1.19  1.09  1.01  0.95
     816   overall  2.12  1.30  0.99  0.86  0.76  0.69  0.66  0.65  0.62  0.60


 List of    132 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   37    5   16    1.57  164.51  0.2410E+02  0.1291E+02 "alien"
   36   28    9    1.58  151.73  0.2223E+02  0.1634E+02 "alien"
   38    4   15    1.59  110.75  0.1622E+02  0.8460E+01 "alien"
   31    6   18    1.61   89.07  0.1305E+02  0.1100E+02 "alien"
   45    2    9    1.59   81.13  0.1188E+02  0.9012E+01 "alien"
   45    2   10    1.57   70.64  0.1035E+02  0.1027E+02 "alien"
   33   25   10    1.68   68.48  0.1913E+02  0.1037E+02 "alien"
   46    1    6    1.63   65.71  0.9626E+01  0.1313E+02 "alien"
   36   32   14    1.39   54.55 -0.2837E+02  0.1711E+02 "alien"
   31    5   17    1.67   54.32  0.1517E+02  0.1208E+02 "alien"
   33   26   12    1.60   52.15  0.7640E+01  0.1277E+02 "alien"
   18   10   21    1.65   51.33  0.1434E+02  0.6841E+01 "alien"
   41    3   13    1.59   49.49  0.7250E+01  0.8531E+01 "alien"
   34    4   16    1.65   44.94  0.1255E+02  0.8138E+01 "alien"
   39   28    2    1.60   44.28  0.6486E+01  0.9569E+01 "alien"
   21   17   20    1.60   44.00  0.6445E+01  0.8778E+01 "alien"
   22   21   17    1.69   43.94  0.1227E+02  0.1040E+02 "alien"
   41    8   18    1.40   43.71 -0.2273E+02  0.1407E+02 "alien"
   44    2   10    1.60   40.06  0.5868E+01  0.8792E+01 "alien"
   42    3   12    1.59   38.75  0.5676E+01  0.1347E+02 "alien"
   46    1    7    1.61   38.67  0.5664E+01  0.9140E+01 "alien"
   30    7   19    1.58   38.46  0.5633E+01  0.1096E+02 "alien"
   27   23   15    1.66   35.24  0.9843E+01  0.8203E+01 "alien"
   35   27   10    1.59   34.48  0.5052E+01  0.9074E+01 "alien"
   38    9   20    1.38   31.36 -0.1631E+02  0.1211E+02 "alien"
   28   12   23    1.40   29.12 -0.1515E+02  0.1046E+02 "alien"
   39    2    2    1.97   27.93 -0.1426E+02  0.9207E+01 "alien"
   34   25    9    1.68   27.61  0.7711E+01  0.8936E+01 "alien"
   51    4   10    1.41   27.55 -0.1433E+02  0.1520E+02 "alien"
   45    1    8    1.62   27.39  0.4012E+01  0.9847E+01 "alien"
   36   31   13    1.43   26.80 -0.1634E+02  0.1006E+02 "alien"
   43    2   10    1.63   26.66  0.7445E+01  0.1033E+02 "alien"
   40   30    4    1.53   26.47  0.2233E+02  0.1053E+02 "alien"
   30   25   14    1.61   26.34  0.3858E+01  0.8326E+01 "alien"
   25   20   14    1.82   25.86 -0.1610E+02  0.7979E+01 "alien"
   23   21   19    1.57   25.44  0.3727E+01  0.7214E+01 "alien"
   53    7   13    1.30   24.58  0.4657E+02  0.1474E+02 "alien"
   28   21   12    1.82   23.74 -0.1478E+02  0.8152E+01 "alien"
   34    5   17    1.60   23.00  0.3369E+01  0.7734E+01 "alien"
   29   26   20    1.38   22.65 -0.1178E+02  0.1488E+02 "alien"
   49    4   11    1.43   22.55 -0.1375E+02  0.9353E+01 "alien"
   46    6   15    1.40   22.28 -0.1159E+02  0.9361E+01 "alien"
   19   14   20    1.66   22.09  0.6170E+01  0.5641E+01 "alien"
   30   20    8    1.96   22.00 -0.1123E+02  0.9020E+01 "alien"
   34   26   10    1.64   21.52  0.6010E+01  0.8742E+01 "alien"
   38   32   11    1.42   21.48 -0.1310E+02  0.9841E+01 "alien"
   38    3   14    1.63   21.33  0.5956E+01  0.8412E+01 "alien"
   37   26    3    1.69   21.17  0.5913E+01  0.9306E+01 "alien"
   31    4   16    1.73   21.08  0.2083E+02  0.1059E+02 "alien"
   28   23   14    1.69   20.81  0.5812E+01  0.8068E+01 "alien"
   37    3   14    1.66   20.68  0.5777E+01  0.7873E+01 "alien"
   48    1    4    1.59   20.55  0.3010E+01  0.9113E+01 "alien"
   30   21    9    1.89   20.30 -0.1037E+02  0.1314E+02 "alien"
   54    4    6    1.39   20.13 -0.1047E+02  0.9989E+01 "alien"
   33    5   17    1.62   20.07  0.2941E+01  0.7340E+01 "alien"
   26   19   12    1.92   19.15 -0.9776E+01  0.8324E+01 "alien"
   38    1    6    1.94   18.85 -0.9622E+01  0.8564E+01 "alien"
   42    7   17    1.42   18.83 -0.1149E+02  0.8922E+01 "alien"
   35    5   17    1.58   18.49  0.2709E+01  0.1189E+02 "alien"
   35   23    3    1.83   18.34 -0.1142E+02  0.9057E+01 "alien"
   31   22   10    1.80   18.28 -0.1138E+02  0.1239E+02 "alien"
   34   22    4    1.87   18.12 -0.1128E+02  0.1030E+02 "alien"
   36    1   11    1.83   18.10 -0.1127E+02  0.1255E+02 "alien"
   29   27   19    1.41   17.58 -0.9142E+01  0.8649E+01 "alien"
   33   32   16    1.38   17.38 -0.9041E+01  0.1567E+02 "alien"
   15   12   18    1.89   17.19 -0.8777E+01  0.6362E+01 "alien"
   25   24   17    1.59   17.13  0.2509E+01  0.8510E+01 "alien"
   39   28    0    1.61   16.28  0.7424E+01  0.1373E+02 "alien"
   29   20   10    1.91   16.25 -0.8296E+01  0.1329E+02 "alien"
   29    2   14    1.92   16.14 -0.8238E+01  0.7261E+01 "alien"
   36   27    8    1.62   15.87  0.2325E+01  0.9247E+01 "alien"
   33   12   23    1.34   15.81  0.2031E+02  0.1456E+02 "alien"
   45    5   14    1.45   15.77 -0.9620E+01  0.9153E+01 "alien"
   51    6   13    1.34   15.65  0.2011E+02  0.1673E+02 "alien"
   49    4   10    1.46   15.59 -0.2702E+02  0.1450E+02 "alien"
   54    3    3    1.42   15.41 -0.9398E+01  0.9782E+01 "alien"
   31   25   13    1.63   14.98  0.2195E+01  0.8676E+01 "alien"
   53    5    9    1.38   14.89 -0.7741E+01  0.1037E+02 "alien"
   51    4    9    1.42   14.82 -0.9041E+01  0.1008E+02 "alien"
   39   29    5    1.56   14.76  0.1245E+02  0.9935E+01 "alien"
   43    1    3    1.78   14.61  0.1443E+02  0.9162E+01 "alien"
   24   20   15    1.78   14.51 -0.9028E+01  0.1324E+02 "alien"
   28   21   22    1.38   14.14 -0.7352E+01  0.8075E+01 "alien"
   47    1    5    1.61   14.07  0.2061E+01  0.9222E+01 "alien"
   31   27   15    1.52   14.05  0.1185E+02  0.1459E+02 "alien"
   40    7   18    1.42   13.96 -0.8516E+01  0.1643E+02 "alien"
   52    2    4    1.47   13.84 -0.2400E+02  0.1790E+02 "alien"
   36   10   21    1.38   13.82 -0.7185E+01  0.8145E+01 "alien"
   28   28   18    1.44   13.81 -0.2191E+02  0.1225E+02 "alien"
   27   19   11    1.95   13.74 -0.7013E+01  0.9039E+01 "alien"
   21   12   22    1.54   13.67  0.1153E+02  0.7967E+01 "alien"
   39   33   10    1.41   13.52 -0.7031E+01  0.1097E+02 "alien"
   54    4    7    1.38   13.50 -0.7019E+01  0.1009E+02 "alien"
   54    3    1    1.43   13.22 -0.8063E+01  0.1001E+02 "alien"
   34   22    3    1.89   13.22 -0.6748E+01  0.9382E+01 "alien"
   22   19   15    1.85   13.14 -0.8181E+01  0.1206E+02 "alien"
   53    4    8    1.39   13.14 -0.6832E+01  0.9849E+01 "alien"
   31    3   15    1.79   13.08 -0.8142E+01  0.7045E+01 "alien"
   35   28   11    1.55   13.02  0.1098E+02  0.9487E+01 "alien"
   54    3    4    1.41   12.99 -0.7920E+01  0.9985E+01 "alien"
   31   24   12    1.68   12.96  0.3621E+01  0.8645E+01 "alien"
   37   32   12    1.42   12.91 -0.7873E+01  0.9979E+01 "alien"
   47    2    9    1.53   12.88  0.1087E+02  0.1665E+02 "alien"
   16   15   17    1.91   12.86 -0.6564E+01  0.9299E+01 "alien"
   33    7   19    1.52   12.65  0.1067E+02  0.7592E+01 "alien"
   33   24    9    1.73   12.27  0.1212E+02  0.8738E+01 "alien"
   25   23   19    1.52   12.16  0.1026E+02  0.1118E+02 "alien"
   49    6   14    1.36   12.05  0.1548E+02  0.1724E+02 "alien"
   17   15   18    1.82   12.04 -0.7496E+01  0.5577E+01 "alien"
   21   17   14    1.99   11.74 -0.5993E+01  0.7908E+01 "alien"
   35   23    1    1.84   11.74 -0.7304E+01  0.9238E+01 "alien"
   44    0    6    1.69   11.60  0.7244E+01  0.9116E+01 "alien"
   13    6   19    1.90   11.54 -0.5894E+01  0.4632E+01 "alien"
   36   25    4    1.73   11.49  0.1136E+02  0.9151E+01 "alien"
   32   25   12    1.64   11.27  0.3147E+01  0.8786E+01 "alien"
   37   29   10    1.51   11.15  0.9404E+01  0.1484E+02 "alien"
   32   27   14    1.53   11.15  0.9402E+01  0.1296E+02 "alien"
   23   17   21    1.52   11.03  0.9300E+01  0.9609E+01 "alien"
   25   21   15    1.74   10.90  0.1078E+02  0.1085E+02 "alien"
   49    2    6    1.53   10.81  0.9119E+01  0.1449E+02 "alien"
   30    3   15    1.82   10.80 -0.6719E+01  0.1032E+02 "alien"
   31    2   14    1.85   10.72 -0.6673E+01  0.7136E+01 "alien"
   32    1   12    1.93   10.71 -0.5469E+01  0.1251E+02 "alien"
   26   24   19    1.49   10.65 -0.1847E+02  0.1139E+02 "alien"
   32   30   19    1.33   10.58  0.1358E+02  0.1391E+02 "alien"
   22   22   18    1.62   10.30  0.3102E+01  0.1239E+02 "alien"
   33    4   16    1.68   10.29  0.2872E+01  0.6973E+01 "alien"
   44    1    8    1.65   10.28  0.2870E+01  0.8746E+01 "alien"
   42    0    6    1.77   10.24  0.2071E+02  0.1339E+02 "alien"
   34    0   10    1.96   10.12 -0.1050E+02  0.8494E+01 "alien"
   36   34   15    1.33   10.05  0.1291E+02  0.1423E+02 "alien"
   34   28   12    1.54   10.04  0.8469E+01  0.9670E+01 "alien"


 List of    132 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   24.58  0.4657E+02  0.1474E+02
   32   30   19    1.33   10.58  0.1358E+02  0.1391E+02
   36   34   15    1.33   10.05  0.1291E+02  0.1423E+02
   33   12   23    1.34   15.81  0.2031E+02  0.1456E+02
   51    6   13    1.34   15.65  0.2011E+02  0.1673E+02
   49    6   14    1.36   12.05  0.1548E+02  0.1724E+02
   28   21   22    1.38   14.14 -0.7352E+01  0.8075E+01
   53    5    9    1.38   14.89 -0.7741E+01  0.1037E+02
   36   10   21    1.38   13.82 -0.7185E+01  0.8145E+01
   33   32   16    1.38   17.38 -0.9041E+01  0.1567E+02
   54    4    7    1.38   13.50 -0.7019E+01  0.1009E+02
   38    9   20    1.38   31.36 -0.1631E+02  0.1211E+02
   29   26   20    1.38   22.65 -0.1178E+02  0.1488E+02
   36   32   14    1.39   54.55 -0.2837E+02  0.1711E+02
   53    4    8    1.39   13.14 -0.6832E+01  0.9849E+01
   54    4    6    1.39   20.13 -0.1047E+02  0.9989E+01
   46    6   15    1.40   22.28 -0.1159E+02  0.9361E+01
   28   12   23    1.40   29.12 -0.1515E+02  0.1046E+02
   41    8   18    1.40   43.71 -0.2273E+02  0.1407E+02
   51    4   10    1.41   27.55 -0.1433E+02  0.1520E+02
   29   27   19    1.41   17.58 -0.9142E+01  0.8649E+01
   39   33   10    1.41   13.52 -0.7031E+01  0.1097E+02
   54    3    4    1.41   12.99 -0.7920E+01  0.9985E+01
   37   32   12    1.42   12.91 -0.7873E+01  0.9979E+01
   42    7   17    1.42   18.83 -0.1149E+02  0.8922E+01
   54    3    3    1.42   15.41 -0.9398E+01  0.9782E+01
   38   32   11    1.42   21.48 -0.1310E+02  0.9841E+01
   40    7   18    1.42   13.96 -0.8516E+01  0.1643E+02
   51    4    9    1.42   14.82 -0.9041E+01  0.1008E+02
   36   31   13    1.43   26.80 -0.1634E+02  0.1006E+02
   54    3    1    1.43   13.22 -0.8063E+01  0.1001E+02
   49    4   11    1.43   22.55 -0.1375E+02  0.9353E+01
   28   28   18    1.44   13.81 -0.2191E+02  0.1225E+02
   45    5   14    1.45   15.77 -0.9620E+01  0.9153E+01
   49    4   10    1.46   15.59 -0.2702E+02  0.1450E+02
   52    2    4    1.47   13.84 -0.2400E+02  0.1790E+02
   26   24   19    1.49   10.65 -0.1847E+02  0.1139E+02
   37   29   10    1.51   11.15  0.9404E+01  0.1484E+02
   25   23   19    1.52   12.16  0.1026E+02  0.1118E+02
   31   27   15    1.52   14.05  0.1185E+02  0.1459E+02
   33    7   19    1.52   12.65  0.1067E+02  0.7592E+01
   23   17   21    1.52   11.03  0.9300E+01  0.9609E+01
   40   30    4    1.53   26.47  0.2233E+02  0.1053E+02
   49    2    6    1.53   10.81  0.9119E+01  0.1449E+02
   47    2    9    1.53   12.88  0.1087E+02  0.1665E+02
   32   27   14    1.53   11.15  0.9402E+01  0.1296E+02
   21   12   22    1.54   13.67  0.1153E+02  0.7967E+01
   34   28   12    1.54   10.04  0.8469E+01  0.9670E+01
   35   28   11    1.55   13.02  0.1098E+02  0.9487E+01
   39   29    5    1.56   14.76  0.1245E+02  0.9935E+01
   45    2   10    1.57   70.64  0.1035E+02  0.1027E+02
   37    5   16    1.57  164.51  0.2410E+02  0.1291E+02
   23   21   19    1.57   25.44  0.3727E+01  0.7214E+01
   35    5   17    1.58   18.49  0.2709E+01  0.1189E+02
   36   28    9    1.58  151.73  0.2223E+02  0.1634E+02
   30    7   19    1.58   38.46  0.5633E+01  0.1096E+02
   48    1    4    1.59   20.55  0.3010E+01  0.9113E+01
   41    3   13    1.59   49.49  0.7250E+01  0.8531E+01
   38    4   15    1.59  110.75  0.1622E+02  0.8460E+01
   25   24   17    1.59   17.13  0.2509E+01  0.8510E+01
   35   27   10    1.59   34.48  0.5052E+01  0.9074E+01
   45    2    9    1.59   81.13  0.1188E+02  0.9012E+01
   42    3   12    1.59   38.75  0.5676E+01  0.1347E+02
   33   26   12    1.60   52.15  0.7640E+01  0.1277E+02
   44    2   10    1.60   40.06  0.5868E+01  0.8792E+01
   34    5   17    1.60   23.00  0.3369E+01  0.7734E+01
   21   17   20    1.60   44.00  0.6445E+01  0.8778E+01
   39   28    2    1.60   44.28  0.6486E+01  0.9569E+01
   46    1    7    1.61   38.67  0.5664E+01  0.9140E+01
   47    1    5    1.61   14.07  0.2061E+01  0.9222E+01
   30   25   14    1.61   26.34  0.3858E+01  0.8326E+01
   39   28    0    1.61   16.28  0.7424E+01  0.1373E+02
   31    6   18    1.61   89.07  0.1305E+02  0.1100E+02
   45    1    8    1.62   27.39  0.4012E+01  0.9847E+01
   36   27    8    1.62   15.87  0.2325E+01  0.9247E+01
   33    5   17    1.62   20.07  0.2941E+01  0.7340E+01
   22   22   18    1.62   10.30  0.3102E+01  0.1239E+02
   31   25   13    1.63   14.98  0.2195E+01  0.8676E+01
   46    1    6    1.63   65.71  0.9626E+01  0.1313E+02
   43    2   10    1.63   26.66  0.7445E+01  0.1033E+02
   38    3   14    1.63   21.33  0.5956E+01  0.8412E+01
   34   26   10    1.64   21.52  0.6010E+01  0.8742E+01
   32   25   12    1.64   11.27  0.3147E+01  0.8786E+01
   34    4   16    1.65   44.94  0.1255E+02  0.8138E+01
   44    1    8    1.65   10.28  0.2870E+01  0.8746E+01
   18   10   21    1.65   51.33  0.1434E+02  0.6841E+01
   37    3   14    1.66   20.68  0.5777E+01  0.7873E+01
   27   23   15    1.66   35.24  0.9843E+01  0.8203E+01
   19   14   20    1.66   22.09  0.6170E+01  0.5641E+01
   31    5   17    1.67   54.32  0.1517E+02  0.1208E+02
   33    4   16    1.68   10.29  0.2872E+01  0.6973E+01
   33   25   10    1.68   68.48  0.1913E+02  0.1037E+02
   31   24   12    1.68   12.96  0.3621E+01  0.8645E+01
   34   25    9    1.68   27.61  0.7711E+01  0.8936E+01
   28   23   14    1.69   20.81  0.5812E+01  0.8068E+01
   37   26    3    1.69   21.17  0.5913E+01  0.9306E+01
   44    0    6    1.69   11.60  0.7244E+01  0.9116E+01
   22   21   17    1.69   43.94  0.1227E+02  0.1040E+02
   31    4   16    1.73   21.08  0.2083E+02  0.1059E+02
   33   24    9    1.73   12.27  0.1212E+02  0.8738E+01
   36   25    4    1.73   11.49  0.1136E+02  0.9151E+01
   25   21   15    1.74   10.90  0.1078E+02  0.1085E+02
   42    0    6    1.77   10.24  0.2071E+02  0.1339E+02
   43    1    3    1.78   14.61  0.1443E+02  0.9162E+01
   24   20   15    1.78   14.51 -0.9028E+01  0.1324E+02
   31    3   15    1.79   13.08 -0.8142E+01  0.7045E+01
   31   22   10    1.80   18.28 -0.1138E+02  0.1239E+02
   17   15   18    1.82   12.04 -0.7496E+01  0.5577E+01
   25   20   14    1.82   25.86 -0.1610E+02  0.7979E+01
   30    3   15    1.82   10.80 -0.6719E+01  0.1032E+02
   28   21   12    1.82   23.74 -0.1478E+02  0.8152E+01
   35   23    3    1.83   18.34 -0.1142E+02  0.9057E+01
   36    1   11    1.83   18.10 -0.1127E+02  0.1255E+02
   35   23    1    1.84   11.74 -0.7304E+01  0.9238E+01
   31    2   14    1.85   10.72 -0.6673E+01  0.7136E+01
   22   19   15    1.85   13.14 -0.8181E+01  0.1206E+02
   34   22    4    1.87   18.12 -0.1128E+02  0.1030E+02
   34   22    3    1.89   13.22 -0.6748E+01  0.9382E+01
   30   21    9    1.89   20.30 -0.1037E+02  0.1314E+02
   15   12   18    1.89   17.19 -0.8777E+01  0.6362E+01
   13    6   19    1.90   11.54 -0.5894E+01  0.4632E+01
   29   20   10    1.91   16.25 -0.8296E+01  0.1329E+02
   16   15   17    1.91   12.86 -0.6564E+01  0.9299E+01
   29    2   14    1.92   16.14 -0.8238E+01  0.7261E+01
   26   19   12    1.92   19.15 -0.9776E+01  0.8324E+01
   32    1   12    1.93   10.71 -0.5469E+01  0.1251E+02
   38    1    6    1.94   18.85 -0.9622E+01  0.8564E+01
   27   19   11    1.95   13.74 -0.7013E+01  0.9039E+01
   30   20    8    1.96   22.00 -0.1123E+02  0.9020E+01
   34    0   10    1.96   10.12 -0.1050E+02  0.8494E+01
   39    2    2    1.97   27.93 -0.1426E+02  0.9207E+01
   21   17   14    1.99   11.74 -0.5993E+01  0.7908E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        55
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1518
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1374
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      144


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
