
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.314    77.314    38.695  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_3.0A_5fr_1deg/stills_3.0A_5fr_1deg_???_001.img SMV
 DATA_RANGE=       2       2
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       2
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000988  0.000641
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
 BEAM_DIVERGENCE_E.S.D.=   0.039
 REFLECTING_RANGE_E.S.D.=   0.098
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

    1473 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       2 ...       2
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1473 REFLECTIONS ACCEPTED



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
 *  31        aP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  0  1  0  0 -1  0  0  0
 *  35        mP          0.0      77.3   38.7   77.3  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  25        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  23        oC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  20        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  21        tP          0.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    39        mC        250.0     159.4   38.7   77.3  90.0  90.0  76.0    0  2  1  0  0  0  1  0  1  0  0  0
    36        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0 -1  0  2  0  1  0  0 -1  0  0
    38        oC        250.0      38.7  159.4   77.3  90.0  90.0 104.0    0  0  1  0  0 -2 -1  0  1  0  0  0
    28        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0 -1  0 -2  0 -1  0  0  1  0  0
    29        mC        250.0      38.7  159.4   77.3  90.0  90.0  76.0    0  0  1  0  0  2  1  0 -1  0  0  0
    27        mC        500.0     159.4   38.7  109.3  90.0 133.3  76.0    0  2 -1  0  0  0 -1  0 -1 -1  0  0
    19        oI        500.0      38.7  109.3  116.0  90.0  70.5  90.0    0  0 -1  0 -1 -1  0  0 -1  1 -1  0
    26        oF        625.0      38.7  159.4  159.4  86.6 104.0 104.0    0  0  1  0  0  2 -1  0 -2  0 -1  0
    18        tI        625.0     109.3  116.0   38.7  70.5  90.0  90.0   -1 -1  0  0  1 -1  1  0  0  0  1  0
    41        mC        998.0     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.0      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.0      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.0      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0 -1 -1  0  1  0  1  0 -1 -1  1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.5   86.5  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.5   86.5  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.3   86.5   86.5  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.5   86.5  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.5   86.5  109.3  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.7   86.5  247.5 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.5   86.5   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.2  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.5   86.5  177.2 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.5   86.5   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.2  109.3   38.7  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
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
 UNIT_CELL_A-AXIS=    73.599   -23.682     0.117
 UNIT_CELL_B-AXIS=    23.682    73.599     0.046
 UNIT_CELL_C-AXIS=    -0.063    -0.004    38.695
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       142 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.10
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.222    77.222    38.648  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  4.1E-02 4.1E-02 2.7E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012950  0.012950  0.025874  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.595   -23.389     0.034
 COORDINATES OF UNIT CELL B-AXIS    23.389    73.595     0.012
 COORDINATES OF UNIT CELL C-AXIS    -0.018    -0.001    38.648
 CRYSTAL MOSAICITY (DEGREES)     0.098
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006332  0.002113
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1260.80   1242.10
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.04
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.22    77.22    38.65  90.000  90.000  90.000


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
  0.05  0.15       111        111.991       1.479
  0.15  0.25       175        123.278       1.344
  0.25  0.35       189        151.734       1.092
  0.35  0.45       149         73.171       2.264
  0.45  0.55       180        234.548       0.706
  0.55  0.65       200        180.962       0.915
  0.65  0.75       181        109.939       1.507
  0.75  0.85       169        139.880       1.184
  0.85  0.95       119        407.154       0.407
 
 -0.05  0.95      1473        165.639       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        5
 DEGREES OF FREEDOM OF CHI^2 FIT               121.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.641
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     1.1 XMAX=     1.9 NXBIN=    1
 YMIN= 0.00168 YMAX= 0.59803 NYBIN=    5
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        284


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        4
 DEGREES OF FREEDOM OF CHI^2 FIT               121.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.640
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     9.1 XMAX=  2442.0 NXBIN=    2
 YMIN=     6.7 YMAX=  2519.5 NYBIN=    2
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        284


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               119.1
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.644
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     1.1 XMAX=     1.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1226    1263
 DETECTOR_SURFACE_POSITION=    1639    1690
 DETECTOR_SURFACE_POSITION=     812    1690
 DETECTOR_SURFACE_POSITION=     812     836
 DETECTOR_SURFACE_POSITION=    1639     836
 DETECTOR_SURFACE_POSITION=    2161    1663
 DETECTOR_SURFACE_POSITION=    1613    2229
 DETECTOR_SURFACE_POSITION=     838    2229
 DETECTOR_SURFACE_POSITION=     290    1663
 DETECTOR_SURFACE_POSITION=     290     863
 DETECTOR_SURFACE_POSITION=     838     297
 DETECTOR_SURFACE_POSITION=    1613     297
 DETECTOR_SURFACE_POSITION=    2161     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        284


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
 8.397E-01  1.864E-03   25.28


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

   24.420   3.393    22.23   1.75      2.45      2.31      50      76      32
    3.393   2.411     3.72   0.99      6.33      6.04     128     170      10
    2.411   1.972    -0.01   0.86    -99.90    -99.90     186     221       2
    1.972   1.709    -0.01   0.69    789.93    963.85     182     234       6
    1.709   1.529     0.03   0.84   1485.03   1649.44     216     276       2
    1.529   1.396     0.01   0.87   3372.90   3524.57     230     297       4
    1.396   1.293    -0.01   0.87    738.25    782.32      94     138       4
   --------------------------------------------------------------------------
   24.420   1.293     1.65   0.89      6.57      6.78    1086    1412      60
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       2       2)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.679  -0.1173E+02  0.1030E+02        -1.14           1 
   47    0    0     1.643   0.6194E+01  0.1048E+02         0.59           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     3.4% OF MEAN INTENSITY


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

     3.84          55      37      1285        2.9%       2.4%      2.3%       36   27.38      3.3%    99.8*     0    0.000       0
     2.73         100      66      2098        3.1%       3.4%      3.0%       68   13.31      4.8%    99.7*     0    0.000       0
     2.23         171     102      2665        3.8%     -99.9%    -99.9%      138    0.00    -99.9%    -2.1      0    0.000       0
     1.94         164      94      3101        3.0%    3592.5%   3773.4%      140   -0.00   5080.5%   -28.4      0    0.000       0
     1.73         181     110      3481        3.2%     401.7%    520.6%      142   -0.02    568.1%    11.4      0    0.000       0
     1.58         209     127      3848        3.3%     879.8%    988.3%      164    0.11   1244.3%    -3.1      0    0.000       0
     1.47         236     145      4143        3.5%     -99.9%    -99.9%      182    0.00    -99.9%    -0.4      0    0.000       0
     1.37         214     129      4429        2.9%     469.3%    502.9%      170    0.14    663.7%     2.2      0    0.000       0
     1.29          82      59      4739        1.2%   22090.1%  19900.2%       46   -0.09  31240.1%     3.4      0    0.000       0
    total        1412     869     29789        2.9%       6.6%      6.8%     1086    2.17      9.3%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1473
 NUMBER OF REJECTED MISFITS                                60
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1412
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    869



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

   WILSON LINE (using all data) : A=   6.222 B=  32.285 CORRELATION=  0.75
      #      RES      SS        <I>       log(<I>)       BO
      15     7.601   0.004  3.6328E+03       8.198    -228.3
      18     4.612   0.012  3.1186E+03       8.045     -77.5
      30     3.491   0.021  1.2192E+03       7.106     -21.5
      30     2.960   0.029  3.1151E+02       5.741       8.4
      43     2.639   0.036 -7.5833E-01    -999.999   -1000.0
      54     2.390   0.044 -1.4737E-01    -999.999   -1000.0
      42     2.203   0.052  1.8393E-01      -1.693      76.8
      47     2.054   0.059  6.5782E-02      -2.721      75.4
      45     1.927   0.067 -1.8237E+00    -999.999   -1000.0
      61     1.825   0.075  4.6948E-01      -0.756      46.5
      52     1.736   0.083  3.7885E-01      -0.971      43.4
      55     1.663   0.090  9.2214E-01      -0.081      34.8
      60     1.593   0.098 -3.1178E-01    -999.999   -1000.0
      74     1.533   0.106 -2.1863E+00    -999.999   -1000.0
      77     1.478   0.115  3.6185E-01      -1.017      31.6
      59     1.430   0.122  1.7008E+00       0.531      23.3
      50     1.388   0.130  1.8226E-01      -1.702      30.6
      40     1.349   0.137 -1.7072E+00    -999.999   -1000.0
      16     1.313   0.145  8.0676E-01      -0.215      22.2


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.601        0.724        0.298        0.097
       3     4.612        0.984        0.526        0.209
       3     3.491        0.379        0.099        0.018
       4     2.960        1.004        0.696        0.345
       3     2.639       18.995       22.974       49.390
       8     2.054     1818.936    27205.021   835505.336
      11     1.736      292.615       74.947    14973.843
       3     1.430       14.026       27.679       43.992
       2     1.349        8.749        9.372        7.350
      89   overall      288.304     1930.109    93599.430


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      10     7.601        1.275        1.749        1.959
      15     4.612        1.040        1.057        0.903
      27     3.491        0.858        0.664        0.476
      26     2.960        2.475        4.902        8.274
      49     2.390       18.895       -7.522      204.082
      39     2.203      565.152      732.805   114941.361
      50     1.825       99.263      347.921     5274.604
      41     1.736      164.403      500.854    11843.296
      46     1.663       27.676      -17.280      482.250
      57     1.593       29.332       -8.363      618.617
      71     1.478      223.004      290.116    25300.680
      56     1.430        9.125       11.151       54.809
      50     1.388     1021.538     5002.753   439606.637
      16     1.313       70.986      182.110     2002.827
     779   overall      314.765     1018.816   525253.858

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.601  0.00  0.00  0.00  0.00  0.38  0.36  0.67  0.64  0.61  0.59
       3     4.612  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       3     3.491  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.960  3.02  2.17  1.80  1.58  1.44  1.34  1.26  1.19  1.14  1.10
       3     2.639  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       8     2.054  2.01  1.45  1.20  1.06  0.96  0.89  0.84  0.80  0.76  0.73
      11     1.736  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.72  0.69  0.67
       3     1.430  1.34  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       2     1.349  4.03  2.89  2.40  2.11  1.92  1.78  1.67  1.59  1.52  1.47
      89   overall  0.86  0.65  0.54  0.47  0.45  0.42  0.41  0.39  0.39  0.38


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      10     7.601  2.10  1.65  2.31  2.12  1.78  1.55  1.39  1.27  1.18  1.11
      15     4.612  1.40  1.84  1.29  1.21  1.02  1.03  0.93  1.09  1.12  1.16
      27     3.491  1.17  0.61  0.57  0.67  0.66  0.90  0.96  0.94  0.94  0.94
      26     2.960  7.28  4.03  2.82  2.22  1.86  1.62  1.45  1.33  1.23  1.16
      49     2.390  5.58  2.93  2.05  1.61  1.35  1.18  1.05  0.96  0.89  0.84
      39     2.203  5.66  2.97  2.08  1.63  1.37  1.19  1.07  0.98  0.91  0.85
      50     1.825  4.83  2.54  1.77  1.40  1.17  1.11  0.99  0.91  0.84  0.79
      41     1.736  5.64  2.96  2.07  1.70  1.43  1.24  1.11  1.02  0.95  0.89
      46     1.663  3.20  1.68  1.26  1.06  0.88  0.77  0.69  0.67  0.62  0.62
      57     1.593  5.53  2.90  2.03  1.60  1.34  1.17  1.05  0.96  0.89  0.83
      71     1.478  4.88  2.56  1.79  1.41  1.18  1.03  0.92  0.84  0.78  0.74
      56     1.430  4.32  2.36  1.72  1.35  1.23  1.15  1.06  0.97  0.93  0.88
      50     1.388  5.04  2.65  1.85  1.46  1.22  1.06  0.95  0.87  0.81  0.76
      16     1.313  5.91  3.10  2.17  1.71  1.43  1.25  1.12  1.02  0.95  0.89
     779   overall  3.37  1.81  1.30  1.04  0.88  0.80  0.72  0.67  0.63  0.60


 List of    172 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   26   18   11    2.01  226.23  0.1488E+02  0.1207E+02 "alien"
   31   20    6    1.99  204.60  0.1346E+02  0.7611E+01 "alien"
   36    2    4    2.09  112.53  0.7403E+01  0.7719E+01 "alien"
   31    3    2    2.46  110.73 -0.1632E+02  0.6866E+01 "alien"
   31   19    4    2.07  105.04  0.6909E+01  0.7909E+01 "alien"
   32   11   22    1.39  100.40  0.1830E+02  0.1021E+02 "alien"
   32   20    3    2.02   98.68  0.6492E+01  0.7754E+01 "alien"
   36   33   14    1.37   94.84  0.1729E+02  0.8379E+01 "alien"
   53    4    8    1.39   91.41  0.1666E+02  0.8551E+01 "alien"
   30    2    7    2.33   89.61 -0.1321E+02  0.6587E+01 "alien"
   42    3   12    1.59   74.92 -0.2336E+02  0.1250E+02 "alien"
   27   15    2    2.48   73.12 -0.1078E+02  0.7101E+01 "alien"
   31   19    1    2.12   71.48  0.4702E+01  0.8680E+01 "alien"
   27   17    7    2.22   69.92  0.1286E+02  0.1308E+02 "alien"
   40   33    9    1.41   68.70  0.1252E+02  0.8338E+01 "alien"
   30   12   23    1.37   67.69  0.1234E+02  0.6494E+01 "alien"
   27   18    9    2.08   67.24  0.4423E+01  0.7253E+01 "alien"
   31    3    4    2.40   59.85 -0.8821E+01  0.1142E+02 "alien"
   27    2    8    2.46   59.12 -0.8713E+01  0.9073E+01 "alien"
   35    5   17    1.58   58.45 -0.1822E+02  0.1048E+02 "alien"
   32    2    7    2.21   56.76  0.1044E+02  0.8384E+01 "alien"
   31    9   21    1.46   56.03  0.2028E+02  0.9841E+01 "alien"
   54    4    6    1.39   54.42  0.9918E+01  0.8354E+01 "alien"
   32    2    6    2.26   54.24  0.9977E+01  0.7305E+01 "alien"
   14   11   13    2.45   53.26 -0.7849E+01  0.6685E+01 "alien"
   30   29   18    1.40   52.88  0.9638E+01  0.1123E+02 "alien"
   21   14   11    2.31   52.29 -0.7706E+01  0.1053E+02 "alien"
   37    2    4    2.04   48.18  0.3169E+01  0.7609E+01 "alien"
   45    7   16    1.39   47.57  0.8670E+01  0.7932E+01 "alien"
   26   13   24    1.38   47.15  0.8593E+01  0.7290E+01 "alien"
   52    4    9    1.40   46.83  0.8535E+01  0.8521E+01 "alien"
   46    7   16    1.37   46.77  0.8525E+01  0.7650E+01 "alien"
   37    2    3    2.06   46.28  0.3044E+01  0.7574E+01 "alien"
   29   22   12    1.77   43.96  0.1666E+02  0.6724E+01 "alien"
   26   16    8    2.24   43.51  0.8003E+01  0.1218E+02 "alien"
   50    5   12    1.39   43.45  0.7919E+01  0.8610E+01 "alien"
   17   15   15    2.05   42.92  0.2824E+01  0.5361E+01 "alien"
   42   32    1    1.46   42.75  0.1547E+02  0.8476E+01 "alien"
   25   17   10    2.13   41.98  0.7721E+01  0.1092E+02 "alien"
   52    2    1    1.48   41.78  0.1512E+02  0.8035E+01 "alien"
   38   31   10    1.46   41.30  0.1495E+02  0.7956E+01 "alien"
   19    1   13    2.40   40.38 -0.5952E+01  0.4518E+01 "alien"
   42    8   18    1.38   40.07  0.7302E+01  0.7785E+01 "alien"
   30    1    9    2.21   39.19  0.7209E+01  0.7022E+01 "alien"
   26   17    9    2.15   38.85  0.7147E+01  0.6880E+01 "alien"
   25    1   13    2.14   38.76  0.7130E+01  0.5975E+01 "alien"
   35   23    1    1.84   38.43  0.1804E+02  0.7611E+01 "alien"
   34    1    8    2.05   36.93  0.2429E+01  0.7150E+01 "alien"
   34   29   13    1.49   35.89  0.1299E+02  0.7569E+01 "alien"
   34    3   15    1.70   35.88  0.1359E+02  0.1052E+02 "alien"
   28   16    3    2.35   35.85 -0.5283E+01  0.7338E+01 "alien"
   45    1    8    1.62   35.42 -0.1104E+02  0.8070E+01 "alien"
   36   24    2    1.78   33.58  0.1272E+02  0.7487E+01 "alien"
   24    9   21    1.57   32.67 -0.1018E+02  0.7247E+01 "alien"
   31    4   16    1.73   32.51  0.1232E+02  0.9175E+01 "alien"
   40    1    6    1.85   32.31  0.1517E+02  0.7474E+01 "alien"
   15    2   14    2.43   32.26 -0.4755E+01  0.5914E+01 "alien"
   38   28    5    1.60   30.70 -0.9572E+01  0.7481E+01 "alien"
   30   28   19    1.38   30.64  0.5585E+01  0.7407E+01 "alien"
   32   20    2    2.04   30.20  0.1987E+01  0.7873E+01 "alien"
   31   26   14    1.57   29.78 -0.9286E+01  0.7265E+01 "alien"
   26   10   22    1.48   29.18  0.1056E+02  0.8009E+01 "alien"
   48    1    3    1.60   29.16 -0.9091E+01  0.7465E+01 "alien"
   20   18   16    1.85   29.02  0.1363E+02  0.7833E+01 "alien"
   29   28   17    1.46   27.83  0.1007E+02  0.7063E+01 "alien"
   16    2   14    2.39   26.73 -0.3939E+01  0.4053E+01 "alien"
   36   24    1    1.78   26.48  0.1243E+02  0.7472E+01 "alien"
   28   16    0    2.39   25.67 -0.2558E+02  0.1006E+02 "alien"
   45    4   13    1.48   25.45  0.9210E+01  0.7948E+01 "alien"
   30   18    3    2.18   25.42  0.4676E+01  0.7654E+01 "alien"
   37   31   11    1.46   25.24  0.9134E+01  0.8424E+01 "alien"
   27   22   14    1.73   25.24  0.9562E+01  0.6293E+01 "alien"
   24    1   13    2.18   25.02  0.4602E+01  0.5924E+01 "alien"
   38    2    0    2.03   24.93  0.1039E+02  0.1115E+02 "alien"
   21   15   12    2.19   24.85  0.4570E+01  0.6048E+01 "alien"
   23   15   10    2.27   24.42  0.4491E+01  0.6412E+01 "alien"
   39   28    2    1.60   24.19 -0.7543E+01  0.7485E+01 "alien"
   38    2    3    2.00   23.93  0.1574E+01  0.7391E+01 "alien"
   30    7   19    1.58   23.89 -0.7448E+01  0.9724E+01 "alien"
   34    2    5    2.18   23.77  0.4371E+01  0.8031E+01 "alien"
   25   10   22    1.50   23.49  0.8500E+01  0.7008E+01 "alien"
   24   13   24    1.40   23.48  0.4279E+01  0.6947E+01 "alien"
   48    9   17    1.30   22.69  0.1830E+02  0.1226E+02 "alien"
   41   11   21    1.29   22.12  0.1785E+02  0.1159E+02 "alien"
   29   18    6    2.13   21.88  0.4024E+01  0.7728E+01 "alien"
   23   21   16    1.73   20.94  0.7934E+01  0.9697E+01 "alien"
   37    5   16    1.57   20.87 -0.6508E+01  0.1097E+02 "alien"
   34    7   19    1.50   20.77  0.7515E+01  0.9207E+01 "alien"
   27    2    7    2.53   20.53 -0.1557E+02  0.1038E+02 "alien"
   31    3    3    2.43   19.95 -0.2940E+01  0.6842E+01 "alien"
   52    3    5    1.46   19.65  0.7110E+01  0.8774E+01 "alien"
   29   19    8    2.02   19.62  0.1291E+01  0.7243E+01 "alien"
   30   18    0    2.21   19.47  0.1431E+02  0.1295E+02 "alien"
   41   31    4    1.48   19.38  0.7013E+01  0.8303E+01 "alien"
   49    3    8    1.50   19.18  0.6940E+01  0.1413E+02 "alien"
   31   25   13    1.62   18.90 -0.5892E+01  0.6565E+01 "alien"
   29   23   13    1.71   18.83  0.7133E+01  0.6524E+01 "alien"
   41   31    2    1.50   18.60  0.6729E+01  0.8089E+01 "alien"
   17   11   20    1.72   18.44  0.6985E+01  0.5160E+01 "alien"
   14   12   14    2.30   18.22 -0.2685E+01  0.4447E+01 "alien"
   29    1    9    2.26   17.92  0.3296E+01  0.1168E+02 "alien"
   41    1    4    1.85   17.88  0.8393E+01  0.7520E+01 "alien"
   27   15    0    2.50   17.86 -0.9246E+01  0.9773E+01 "alien"
   14    5   17    2.08   17.51  0.1152E+01  0.4903E+01 "alien"
   28   24   15    1.63   17.27 -0.5386E+01  0.6670E+01 "alien"
   21   15   21    1.57   17.01 -0.5304E+01  0.7826E+01 "alien"
   33   11   22    1.38   16.96  0.3091E+01  0.6900E+01 "alien"
   21   12    7    2.76   16.76 -0.1271E+02  0.8418E+01 "alien"
   42    1    4    1.81   16.73  0.7856E+01  0.7364E+01 "alien"
   50    8   15    1.31   16.70  0.1347E+02  0.9574E+01 "alien"
   22    1   13    2.27   16.62  0.3056E+01  0.5234E+01 "alien"
   52    2    2    1.48   16.39  0.5930E+01  0.7937E+01 "alien"
   52    5   10    1.38   16.31  0.2973E+01  0.8291E+01 "alien"
   43    2   10    1.63   16.25  0.1499E+02  0.1127E+02 "alien"
   30    0   11    2.08   16.23  0.3312E+01  0.6547E+01 "alien"
   29    2    7    2.39   16.21 -0.2389E+01  0.6324E+01 "alien"
   12   11   14    2.39   16.02 -0.2361E+01  0.5745E+01 "alien"
   34   31   15    1.41   15.98  0.2913E+01  0.7962E+01 "alien"
   42    1    1    1.84   15.97  0.7496E+01  0.7652E+01 "alien"
   28   27   17    1.50   15.87  0.5741E+01  0.6786E+01 "alien"
   31   22    9    1.84   15.82  0.7425E+01  0.7882E+01 "alien"
   32   29   15    1.47   15.68  0.5675E+01  0.7717E+01 "alien"
   27   20   12    1.87   15.58  0.7313E+01  0.6550E+01 "alien"
   45    0    2    1.71   15.46  0.1715E+02  0.1449E+02 "alien"
   26   14    2    2.59   15.43 -0.1170E+02  0.8484E+01 "alien"
   43    0    7    1.71   15.37  0.1695E+02  0.1127E+02 "alien"
   14    4   16    2.20   15.01  0.2761E+01  0.4769E+01 "alien"
   39   31    8    1.48   14.87  0.5380E+01  0.8101E+01 "alien"
   21   17   20    1.60   14.79 -0.4610E+01  0.7108E+01 "alien"
   31   19    0    2.12   14.67  0.2670E+01  0.1331E+02 "alien"
   20   17   19    1.67   14.62  0.1348E+02  0.8328E+01 "alien"
   23   23   18    1.59   14.44 -0.1224E+02  0.8644E+01 "alien"
   14   14   16    2.05   14.22  0.2506E+01  0.4639E+01 "alien"
   25   20   14    1.82   14.07  0.6606E+01  0.6580E+01 "alien"
   30    3    3    2.51   14.05 -0.1066E+02  0.6844E+01 "alien"
   17   13   13    2.29   13.99 -0.2061E+01  0.5169E+01 "alien"
   40   29    0    1.56   13.72 -0.1106E+02  0.1096E+02 "alien"
   31    2    6    2.32   13.66 -0.2013E+01  0.1321E+02 "alien"
   18   11    9    2.79   13.38 -0.1015E+02  0.5770E+01 "alien"
   29   21   11    1.84   13.26  0.6224E+01  0.6531E+01 "alien"
   37   27    7    1.61   12.99 -0.4049E+01  0.7343E+01 "alien"
   45    1    7    1.64   12.93  0.1192E+02  0.7429E+01 "alien"
   35   11   22    1.35   12.71 -0.2170E+02  0.9528E+01 "alien"
   35   30   13    1.46   12.68  0.4589E+01  0.7591E+01 "alien"
   32    3    0    2.40   12.61 -0.4447E+01  0.9877E+01 "alien"
   46    6   15    1.40   12.56  0.2289E+01  0.7978E+01 "alien"
   36   24    0    1.78   12.27  0.1347E+02  0.1053E+02 "alien"
   15   10   19    1.84   12.18  0.5716E+01  0.5491E+01 "alien"
   20   19   18    1.70   12.17  0.4612E+01  0.5432E+01 "alien"
   42   32    3    1.45   12.17  0.2070E+02  0.1503E+02 "alien"
   26   21   14    1.77   12.16  0.4607E+01  0.6171E+01 "alien"
   23   15   22    1.49   11.88  0.4300E+01  0.5550E+01 "alien"
   44    0    5    1.71   11.85  0.1021E+02  0.7298E+01 "alien"
   42   32    0    1.46   11.82  0.9704E+01  0.1158E+02 "alien"
   35    2    5    2.12   11.64  0.7654E+00  0.1061E+02 "alien"
   20   13   21    1.60   11.62 -0.3622E+01  0.6069E+01 "alien"
   49    6   14    1.36   11.46 -0.1957E+02  0.1377E+02 "alien"
   22   19   15    1.85   11.43  0.5365E+01  0.1003E+02 "alien"
   25   14    5    2.54   11.43 -0.8665E+01  0.6441E+01 "alien"
   34    2   13    1.80   11.14  0.5231E+01  0.1112E+02 "alien"
   32    8   20    1.49   11.12  0.4025E+01  0.6826E+01 "alien"
   40   10   20    1.34   11.10 -0.1896E+02  0.1139E+02 "alien"
   39   28    0    1.61   10.98 -0.7350E+01  0.1098E+02 "alien"
   20   12    9    2.62   10.93 -0.8292E+01  0.9276E+01 "alien"
   32   31   18    1.35   10.91 -0.1862E+02  0.1115E+02 "alien"
   54    3    4    1.41   10.85  0.1846E+02  0.8442E+01 "alien"
   51    7   14    1.32   10.81  0.8719E+01  0.8620E+01 "alien"
   27   24   20    1.43   10.63  0.1809E+02  0.1057E+02 "alien"
   29   20   10    1.91   10.58 -0.1930E+02  0.1120E+02 "alien"
   30   23   12    1.73   10.51  0.3983E+01  0.6714E+01 "alien"
   21    7   24    1.46   10.46  0.3786E+01  0.5477E+01 "alien"
   25    2    8    2.60   10.44 -0.7920E+01  0.5922E+01 "alien"


 List of    172 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   41   11   21    1.29   22.12  0.1785E+02  0.1159E+02
   48    9   17    1.30   22.69  0.1830E+02  0.1226E+02
   50    8   15    1.31   16.70  0.1347E+02  0.9574E+01
   51    7   14    1.32   10.81  0.8719E+01  0.8620E+01
   40   10   20    1.34   11.10 -0.1896E+02  0.1139E+02
   32   31   18    1.35   10.91 -0.1862E+02  0.1115E+02
   35   11   22    1.35   12.71 -0.2170E+02  0.9528E+01
   49    6   14    1.36   11.46 -0.1957E+02  0.1377E+02
   46    7   16    1.37   46.77  0.8525E+01  0.7650E+01
   36   33   14    1.37   94.84  0.1729E+02  0.8379E+01
   30   12   23    1.37   67.69  0.1234E+02  0.6494E+01
   26   13   24    1.38   47.15  0.8593E+01  0.7290E+01
   33   11   22    1.38   16.96  0.3091E+01  0.6900E+01
   52    5   10    1.38   16.31  0.2973E+01  0.8291E+01
   30   28   19    1.38   30.64  0.5585E+01  0.7407E+01
   42    8   18    1.38   40.07  0.7302E+01  0.7785E+01
   50    5   12    1.39   43.45  0.7919E+01  0.8610E+01
   45    7   16    1.39   47.57  0.8670E+01  0.7932E+01
   53    4    8    1.39   91.41  0.1666E+02  0.8551E+01
   32   11   22    1.39  100.40  0.1830E+02  0.1021E+02
   54    4    6    1.39   54.42  0.9918E+01  0.8354E+01
   46    6   15    1.40   12.56  0.2289E+01  0.7978E+01
   24   13   24    1.40   23.48  0.4279E+01  0.6947E+01
   52    4    9    1.40   46.83  0.8535E+01  0.8521E+01
   30   29   18    1.40   52.88  0.9638E+01  0.1123E+02
   34   31   15    1.41   15.98  0.2913E+01  0.7962E+01
   40   33    9    1.41   68.70  0.1252E+02  0.8338E+01
   54    3    4    1.41   10.85  0.1846E+02  0.8442E+01
   27   24   20    1.43   10.63  0.1809E+02  0.1057E+02
   42   32    3    1.45   12.17  0.2070E+02  0.1503E+02
   37   31   11    1.46   25.24  0.9134E+01  0.8424E+01
   52    3    5    1.46   19.65  0.7110E+01  0.8774E+01
   38   31   10    1.46   41.30  0.1495E+02  0.7956E+01
   31    9   21    1.46   56.03  0.2028E+02  0.9841E+01
   35   30   13    1.46   12.68  0.4589E+01  0.7591E+01
   42   32    1    1.46   42.75  0.1547E+02  0.8476E+01
   21    7   24    1.46   10.46  0.3786E+01  0.5477E+01
   42   32    0    1.46   11.82  0.9704E+01  0.1158E+02
   29   28   17    1.46   27.83  0.1007E+02  0.7063E+01
   32   29   15    1.47   15.68  0.5675E+01  0.7717E+01
   39   31    8    1.48   14.87  0.5380E+01  0.8101E+01
   52    2    2    1.48   16.39  0.5930E+01  0.7937E+01
   45    4   13    1.48   25.45  0.9210E+01  0.7948E+01
   52    2    1    1.48   41.78  0.1512E+02  0.8035E+01
   26   10   22    1.48   29.18  0.1056E+02  0.8009E+01
   41   31    4    1.48   19.38  0.7013E+01  0.8303E+01
   23   15   22    1.49   11.88  0.4300E+01  0.5550E+01
   32    8   20    1.49   11.12  0.4025E+01  0.6826E+01
   34   29   13    1.49   35.89  0.1299E+02  0.7569E+01
   28   27   17    1.50   15.87  0.5741E+01  0.6786E+01
   49    3    8    1.50   19.18  0.6940E+01  0.1413E+02
   41   31    2    1.50   18.60  0.6729E+01  0.8089E+01
   25   10   22    1.50   23.49  0.8500E+01  0.7008E+01
   34    7   19    1.50   20.77  0.7515E+01  0.9207E+01
   40   29    0    1.56   13.72 -0.1106E+02  0.1096E+02
   21   15   21    1.57   17.01 -0.5304E+01  0.7826E+01
   31   26   14    1.57   29.78 -0.9286E+01  0.7265E+01
   24    9   21    1.57   32.67 -0.1018E+02  0.7247E+01
   37    5   16    1.57   20.87 -0.6508E+01  0.1097E+02
   35    5   17    1.58   58.45 -0.1822E+02  0.1048E+02
   30    7   19    1.58   23.89 -0.7448E+01  0.9724E+01
   23   23   18    1.59   14.44 -0.1224E+02  0.8644E+01
   42    3   12    1.59   74.92 -0.2336E+02  0.1250E+02
   48    1    3    1.60   29.16 -0.9091E+01  0.7465E+01
   20   13   21    1.60   11.62 -0.3622E+01  0.6069E+01
   38   28    5    1.60   30.70 -0.9572E+01  0.7481E+01
   21   17   20    1.60   14.79 -0.4610E+01  0.7108E+01
   39   28    2    1.60   24.19 -0.7543E+01  0.7485E+01
   39   28    0    1.61   10.98 -0.7350E+01  0.1098E+02
   37   27    7    1.61   12.99 -0.4049E+01  0.7343E+01
   45    1    8    1.62   35.42 -0.1104E+02  0.8070E+01
   31   25   13    1.62   18.90 -0.5892E+01  0.6565E+01
   28   24   15    1.63   17.27 -0.5386E+01  0.6670E+01
   43    2   10    1.63   16.25  0.1499E+02  0.1127E+02
   45    1    7    1.64   12.93  0.1192E+02  0.7429E+01
   20   17   19    1.67   14.62  0.1348E+02  0.8328E+01
   34    3   15    1.70   35.88  0.1359E+02  0.1052E+02
   20   19   18    1.70   12.17  0.4612E+01  0.5432E+01
   29   23   13    1.71   18.83  0.7133E+01  0.6524E+01
   43    0    7    1.71   15.37  0.1695E+02  0.1127E+02
   45    0    2    1.71   15.46  0.1715E+02  0.1449E+02
   44    0    5    1.71   11.85  0.1021E+02  0.7298E+01
   17   11   20    1.72   18.44  0.6985E+01  0.5160E+01
   30   23   12    1.73   10.51  0.3983E+01  0.6714E+01
   31    4   16    1.73   32.51  0.1232E+02  0.9175E+01
   27   22   14    1.73   25.24  0.9562E+01  0.6293E+01
   23   21   16    1.73   20.94  0.7934E+01  0.9697E+01
   29   22   12    1.77   43.96  0.1666E+02  0.6724E+01
   26   21   14    1.77   12.16  0.4607E+01  0.6171E+01
   36   24    2    1.78   33.58  0.1272E+02  0.7487E+01
   36   24    1    1.78   26.48  0.1243E+02  0.7472E+01
   36   24    0    1.78   12.27  0.1347E+02  0.1053E+02
   34    2   13    1.80   11.14  0.5231E+01  0.1112E+02
   42    1    4    1.81   16.73  0.7856E+01  0.7364E+01
   25   20   14    1.82   14.07  0.6606E+01  0.6580E+01
   42    1    1    1.84   15.97  0.7496E+01  0.7652E+01
   31   22    9    1.84   15.82  0.7425E+01  0.7882E+01
   15   10   19    1.84   12.18  0.5716E+01  0.5491E+01
   29   21   11    1.84   13.26  0.6224E+01  0.6531E+01
   35   23    1    1.84   38.43  0.1804E+02  0.7611E+01
   20   18   16    1.85   29.02  0.1363E+02  0.7833E+01
   41    1    4    1.85   17.88  0.8393E+01  0.7520E+01
   40    1    6    1.85   32.31  0.1517E+02  0.7474E+01
   22   19   15    1.85   11.43  0.5365E+01  0.1003E+02
   27   20   12    1.87   15.58  0.7313E+01  0.6550E+01
   29   20   10    1.91   10.58 -0.1930E+02  0.1120E+02
   31   20    6    1.99  204.60  0.1346E+02  0.7611E+01
   38    2    3    2.00   23.93  0.1574E+01  0.7391E+01
   26   18   11    2.01  226.23  0.1488E+02  0.1207E+02
   32   20    3    2.02   98.68  0.6492E+01  0.7754E+01
   29   19    8    2.02   19.62  0.1291E+01  0.7243E+01
   38    2    0    2.03   24.93  0.1039E+02  0.1115E+02
   32   20    2    2.04   30.20  0.1987E+01  0.7873E+01
   37    2    4    2.04   48.18  0.3169E+01  0.7609E+01
   14   14   16    2.05   14.22  0.2506E+01  0.4639E+01
   34    1    8    2.05   36.93  0.2429E+01  0.7150E+01
   17   15   15    2.05   42.92  0.2824E+01  0.5361E+01
   37    2    3    2.06   46.28  0.3044E+01  0.7574E+01
   31   19    4    2.07  105.04  0.6909E+01  0.7909E+01
   30    0   11    2.08   16.23  0.3312E+01  0.6547E+01
   27   18    9    2.08   67.24  0.4423E+01  0.7253E+01
   14    5   17    2.08   17.51  0.1152E+01  0.4903E+01
   36    2    4    2.09  112.53  0.7403E+01  0.7719E+01
   35    2    5    2.12   11.64  0.7654E+00  0.1061E+02
   31   19    1    2.12   71.48  0.4702E+01  0.8680E+01
   31   19    0    2.12   14.67  0.2670E+01  0.1331E+02
   25   17   10    2.13   41.98  0.7721E+01  0.1092E+02
   29   18    6    2.13   21.88  0.4024E+01  0.7728E+01
   25    1   13    2.14   38.76  0.7130E+01  0.5975E+01
   26   17    9    2.15   38.85  0.7147E+01  0.6880E+01
   30   18    3    2.18   25.42  0.4676E+01  0.7654E+01
   34    2    5    2.18   23.77  0.4371E+01  0.8031E+01
   24    1   13    2.18   25.02  0.4602E+01  0.5924E+01
   21   15   12    2.19   24.85  0.4570E+01  0.6048E+01
   14    4   16    2.20   15.01  0.2761E+01  0.4769E+01
   30    1    9    2.21   39.19  0.7209E+01  0.7022E+01
   30   18    0    2.21   19.47  0.1431E+02  0.1295E+02
   32    2    7    2.21   56.76  0.1044E+02  0.8384E+01
   27   17    7    2.22   69.92  0.1286E+02  0.1308E+02
   26   16    8    2.24   43.51  0.8003E+01  0.1218E+02
   32    2    6    2.26   54.24  0.9977E+01  0.7305E+01
   29    1    9    2.26   17.92  0.3296E+01  0.1168E+02
   22    1   13    2.27   16.62  0.3056E+01  0.5234E+01
   23   15   10    2.27   24.42  0.4491E+01  0.6412E+01
   17   13   13    2.29   13.99 -0.2061E+01  0.5169E+01
   14   12   14    2.30   18.22 -0.2685E+01  0.4447E+01
   21   14   11    2.31   52.29 -0.7706E+01  0.1053E+02
   31    2    6    2.32   13.66 -0.2013E+01  0.1321E+02
   30    2    7    2.33   89.61 -0.1321E+02  0.6587E+01
   28   16    3    2.35   35.85 -0.5283E+01  0.7338E+01
   12   11   14    2.39   16.02 -0.2361E+01  0.5745E+01
   16    2   14    2.39   26.73 -0.3939E+01  0.4053E+01
   29    2    7    2.39   16.21 -0.2389E+01  0.6324E+01
   28   16    0    2.39   25.67 -0.2558E+02  0.1006E+02
   19    1   13    2.40   40.38 -0.5952E+01  0.4518E+01
   31    3    4    2.40   59.85 -0.8821E+01  0.1142E+02
   32    3    0    2.40   12.61 -0.4447E+01  0.9877E+01
   15    2   14    2.43   32.26 -0.4755E+01  0.5914E+01
   31    3    3    2.43   19.95 -0.2940E+01  0.6842E+01
   14   11   13    2.45   53.26 -0.7849E+01  0.6685E+01
   27    2    8    2.46   59.12 -0.8713E+01  0.9073E+01
   31    3    2    2.46  110.73 -0.1632E+02  0.6866E+01
   27   15    2    2.48   73.12 -0.1078E+02  0.7101E+01
   27   15    0    2.50   17.86 -0.9246E+01  0.9773E+01
   30    3    3    2.51   14.05 -0.1066E+02  0.6844E+01
   27    2    7    2.53   20.53 -0.1557E+02  0.1038E+02
   25   14    5    2.54   11.43 -0.8665E+01  0.6441E+01
   26   14    2    2.59   15.43 -0.1170E+02  0.8484E+01
   25    2    8    2.60   10.44 -0.7920E+01  0.5922E+01
   20   12    9    2.62   10.93 -0.8292E+01  0.9276E+01
   21   12    7    2.76   16.76 -0.1271E+02  0.8418E+01
   18   11    9    2.79   13.38 -0.1015E+02  0.5770E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        87
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1473
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1268
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      205


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
