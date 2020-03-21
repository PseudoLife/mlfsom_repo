
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.396    77.396    38.740  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      15      15
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      15
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000793  0.000678
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
 BEAM_DIVERGENCE_E.S.D.=   0.043
 REFLECTING_RANGE_E.S.D.=   0.097
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

    1458 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      15 ...      15
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1458 REFLECTIONS ACCEPTED



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

 *  31        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0 -1  0  0  0  0 -1  0  0
 *  20        mC          0.0     109.5  109.5   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  44        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  25        mC          0.0     109.5  109.5   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  33        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  23        oC          0.0     109.5  109.5   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  35        mP          0.0      77.4   38.7   77.4  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  21        tP          0.0      77.4   77.4   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
 *  32        oP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    28        mC        250.0      38.7  159.6   77.4  90.0  90.0  75.9    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    29        mC        250.0      38.7  159.6   77.4  90.0  90.0  75.9    0  0  1  0  2  0  1  0  0  1  0  0
    37        mC        250.0     159.6   38.7   77.4  90.0  90.0  75.9    0  2  1  0  0  0  1  0  1  0  0  0
    39        mC        250.0     159.6   38.7   77.4  90.0  90.0  75.9   -2  0 -1  0  0  0 -1  0  0 -1  0  0
    36        oC        250.0      38.7  159.6   77.4  90.0  90.0 104.1    0  0  1  0  0 -2 -1  0  1  0  0  0
    38        oC        250.0      38.7  159.6   77.4  90.0  90.0 104.1    0  0 -1  0  2  0  1  0  0 -1  0  0
    19        oI        500.0      38.7  109.5  116.1  90.0  70.5  90.0    0  0 -1  0  1 -1  0  0 -1 -1 -1  0
    27        mC        500.0     159.6   38.7  109.5  90.0 133.3  75.9   -2  0 -1  0  0  0 -1  0  1 -1  0  0
    18        tI        625.0     109.5  116.1   38.7  70.5  90.0  90.0    1 -1  0  0  1  1  1  0  0  0  1  0
    26        oF        625.0      38.7  159.6  159.6  86.6 104.1 104.1    0  0  1  0 -2  0 -1  0  0 -2 -1  0
    30        mC        997.9      77.4  173.1   38.7  90.0  90.0  63.4   -1  0  0  0 -1  2  0  0  0  0 -1  0
    41        mC        997.9     173.1   77.4   38.7  90.0  90.0  63.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    22        hP        997.9      77.4   77.4   38.7  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
    40        oC        997.9      77.4  173.1   38.7  90.0  90.0 116.6   -1  0  0  0  1  2  0  0  0  0 -1  0
     1        cF        999.0     116.1  116.1  116.1  96.4  96.4 141.0    1 -1  1  0 -1  1  1  0 -1 -1 -1  0
     2        hR        999.0      86.5   86.5  116.1 116.5  63.5 101.6    1  0  1  0  0 -1 -1  0  1  1 -1  0
     3        cP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
     5        cI        999.0      86.5   86.5  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     4        hR        999.0      86.5   86.5  116.1 116.5  63.5 101.6   -1  0  1  0  0  1 -1  0 -1 -1 -1  0
     6        tI        999.0     109.5   86.5   86.5  78.4  50.8  50.8    1  1  0  0  0  1  1  0  1  0  1  0
     7        tI        999.0      86.5   86.5  109.5  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     8        oI        999.0      86.5   86.5  109.5  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     9        hR        999.0      38.7   86.5  247.8 102.1  99.0 116.6    0  0  1  0 -1  0 -1  0  1 -3 -1  0
    10        mC        999.0      86.5   86.5   77.4  90.0  90.0 126.8   -1  0  1  0  1  0  1  0  0  1  0  0
    11        tP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    12        hP        999.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
    13        oC        999.0      86.5   86.5   77.4  90.0  90.0  53.2    1  0  1  0  1  0 -1  0  0  1  0  0
    15        tI        999.0      38.7   77.4  177.3  64.1  77.4  90.0    0  0  1  0  0 -1  0  0  2 -1  1  0
    16        oF        999.0      86.5   86.5  177.3 107.0 119.2  53.2   -1  0 -1  0 -1  0  1  0  1  2  1  0
    14        mC        999.0      86.5   86.5   77.4  90.0  90.0  53.2    0 -1  1  0  0 -1 -1  0  1  0  0  0
    17        mC        999.0      86.5   86.5   86.5  78.4 101.6  53.2   -1  0  1  0 -1  0 -1  0  0 -1 -1  0
    24        hR        999.0     177.3  109.5   38.7  90.0  77.4 108.0    2  1  1  0 -1  1  0  0  0  0  1  0
    42        oI        999.0      38.7   77.4  177.3 115.9 102.6  90.0    0  0 -1  0 -1  0  0  0  1  2  1  0
    43        mI        999.0      86.5  177.3   77.4 115.9 153.4  60.8    1  0  1  0  1  2  1  0 -1  0  0  0

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

  *   0.00        0    6.0      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.761   -23.440     0.027
 UNIT_CELL_B-AXIS=    23.440    73.761     0.016
 UNIT_CELL_C-AXIS=    -0.016    -0.004    38.739
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       298 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.26
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.342    77.342    38.706  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  7.2E-02 7.2E-02 3.6E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012930  0.012930  0.025836  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.704   -23.440     0.035
 COORDINATES OF UNIT CELL B-AXIS    23.440    73.704     0.010
 COORDINATES OF UNIT CELL C-AXIS    -0.018     0.000    38.706
 CRYSTAL MOSAICITY (DEGREES)     0.097
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006281  0.002062
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1260.87   1242.02
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.13
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.34    77.34    38.71  90.000  90.000  90.000


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
  0.05  0.15       108        153.774       1.282
  0.15  0.25       170        200.488       0.983
  0.25  0.35       199        182.330       1.081
  0.35  0.45       144        115.282       1.709
  0.45  0.55       165        287.944       0.684
  0.55  0.65       206        185.932       1.060
  0.65  0.75       178        181.830       1.084
  0.75  0.85       170        174.779       1.128
  0.85  0.95       118        303.842       0.649
 
 -0.05  0.95      1458        197.064       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               250.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.953
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    14.1 XMAX=    14.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59621 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        540


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               250.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.952
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     8.3 XMAX=  2436.2 NXBIN=    3
 YMIN=     6.3 YMAX=  2517.7 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        540


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               250.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.949
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    14.1 XMAX=    14.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1222    1262
 DETECTOR_SURFACE_POSITION=    1635    1689
 DETECTOR_SURFACE_POSITION=     810    1689
 DETECTOR_SURFACE_POSITION=     810     835
 DETECTOR_SURFACE_POSITION=    1635     835
 DETECTOR_SURFACE_POSITION=    2155    1662
 DETECTOR_SURFACE_POSITION=    1609    2227
 DETECTOR_SURFACE_POSITION=     836    2227
 DETECTOR_SURFACE_POSITION=     289    1662
 DETECTOR_SURFACE_POSITION=     289     862
 DETECTOR_SURFACE_POSITION=     836     297
 DETECTOR_SURFACE_POSITION=    1609     297
 DETECTOR_SURFACE_POSITION=    2155     862
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        540


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
 1.091E+00  1.468E-03   24.98


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

   24.458   3.398    20.86   1.60      2.37      2.39      60      85      22
    3.398   2.414    13.85   1.53      4.36      3.36     128     169      12
    2.414   1.975     6.50   1.05      7.57      7.10     186     216       2
    1.975   1.711     0.02   0.84   1028.17   1105.92     188     234       0
    1.711   1.532     0.13   0.50    275.17    370.64     220     272       0
    1.532   1.399    -0.07   0.64    -99.90    -99.90     240     299       0
    1.399   1.295     0.06   0.71    969.05   1178.34      94     144       2
   --------------------------------------------------------------------------
   24.458   1.295     3.91   0.87      6.04      6.21    1116    1419      38
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      15      15)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681   0.1232E+02  0.1271E+02         0.97           1 
   47    0    0     1.646   0.1041E+02  0.1286E+02         0.81           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     5.0% OF MEAN INTENSITY


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

     3.84          63      41      1285        3.2%       2.4%      2.3%       44   26.30      3.3%    99.9*     0    0.000       0
     2.73         107      68      2098        3.2%       3.5%      2.9%       78   20.59      5.0%    99.7*     0    0.000       0
     2.24         163      97      2665        3.6%       5.1%      4.1%      132   14.13      7.3%    98.9*     0    0.000       0
     1.94         162      92      3101        3.0%      11.4%     11.4%      140    5.03     16.1%    98.0*     0    0.000       0
     1.74         179     106      3479        3.0%    6406.2%   6663.3%      146    0.00   9059.7%   -22.7      0    0.000       0
     1.59         211     125      3850        3.2%     241.0%    338.0%      172    0.14    340.9%     3.4      0    0.000       0
     1.47         230     140      4143        3.4%     -99.9%    -99.9%      180    0.00    -99.9%     7.9      0    0.000       0
     1.37         220     130      4428        2.9%    1023.9%   1249.0%      180    0.01   1448.0%     2.5      0    0.000       0
     1.30          84      62      4744        1.3%     -99.9%    -99.9%       44    0.00    -99.9%   -13.6      0    0.000       0
    total        1419     861     29793        2.9%       6.0%      6.2%     1116    5.03      8.5%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1458
 NUMBER OF REJECTED MISFITS                                38
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1419
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    861



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

   WILSON LINE (using all data) : A=   7.630 B=  32.720 CORRELATION=  0.84
      #      RES      SS        <I>       log(<I>)       BO
      17     7.253   0.005  2.0968E+03       7.648      -2.0
      20     4.647   0.012  2.7070E+03       7.904     -11.8
      29     3.503   0.020  1.0988E+03       7.002      15.4
      33     2.993   0.028  4.8802E+02       6.190      25.8
      40     2.648   0.036  4.7693E+02       6.167      20.5
      53     2.395   0.044  2.5051E+02       5.523      24.2
      40     2.204   0.051  1.4087E+02       4.948      26.1
      45     2.056   0.059  6.0509E+01       4.103      29.8
      46     1.931   0.067 -3.7623E-01    -999.999   -1000.0
      57     1.829   0.075  1.0527E+00       0.051      50.7
      52     1.739   0.083  4.3507E-01      -0.832      51.2
      55     1.665   0.090  1.3770E+00       0.320      40.5
      56     1.596   0.098  2.7182E-01      -1.303      45.5
      72     1.536   0.106  1.0397E+00       0.039      35.8
      76     1.480   0.114 -1.7618E+00    -999.999   -1000.0
      59     1.432   0.122 -9.4671E-01    -999.999   -1000.0
      51     1.390   0.129  1.2602E+00       0.231      28.6
      42     1.351   0.137 -2.1173E+00    -999.999   -1000.0
      18     1.314   0.145  9.9887E+00       2.301      18.4


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.253        1.851        1.230        0.648
       3     4.647        1.144        0.663        0.286
       3     3.503        0.333        0.081        0.014
       4     2.993        1.042        0.449        0.146
       3     2.648        0.276        0.065        0.011
       5     2.395        0.851        0.435        0.170
       2     2.204        1.768        0.846        0.295
       8     2.056        0.282        0.099        0.027
      11     1.665       12.952        2.879       23.412
      87   overall       44.089     -191.964     3078.042


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.253        0.443        0.304        0.168
      17     4.647        1.052        1.132        1.030
      26     3.503        0.875        0.703        0.529
      29     2.993        0.947        0.927        0.817
      37     2.648        1.130        1.189        1.127
      48     2.395        0.953        0.814        0.589
      38     2.204        0.774        0.582        0.392
      37     2.056        1.509        2.144        2.829
      39     1.931       29.446      -13.288      317.768
      47     1.829       29.072       59.321      404.718
      43     1.739      146.504      143.990    12505.969
      44     1.665       12.044       -3.290       69.106
      53     1.596      492.454    -6325.793   263631.009
      69     1.536       37.100      -19.228      832.915
      51     1.390       22.116       38.746      332.532
      18     1.314        1.173        0.921        0.746
     774   overall       49.420     -420.560    18872.032

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.253  0.00  0.00  0.00  0.00  0.00  0.00  0.33  0.32  0.30  0.29
       3     4.647  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.503  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.993  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.648  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       5     2.395  0.00  0.00  0.00  0.42  1.15  1.07  1.00  0.95  0.91  0.88
       2     2.204  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.056  1.51  1.45  1.20  1.06  1.20  1.11  1.05  0.99  0.95  1.10
      11     1.665  1.10  0.79  0.65  0.58  0.52  0.49  0.46  0.58  0.55  0.53
      87   overall  0.42  0.33  0.28  0.27  0.31  0.31  0.31  0.31  0.31  0.34


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.253  2.63  2.76  1.93  1.52  1.91  1.85  1.66  1.51  1.40  1.32
      17     4.647  1.24  1.62  1.13  0.89  1.20  1.17  1.05  1.18  1.19  1.12
      26     3.503  0.81  0.64  0.59  0.70  0.68  0.77  0.92  0.98  0.91  0.85
      29     2.993  1.09  1.71  1.20  1.05  1.23  1.30  1.23  1.19  1.10  1.09
      37     2.648  0.85  1.34  1.25  1.15  0.96  1.02  1.02  1.08  1.09  1.03
      48     2.395  1.31  1.15  1.04  1.07  1.11  1.11  1.08  1.06  1.05  1.05
      38     2.204  0.55  1.02  1.02  0.88  0.94  0.93  0.89  0.96  0.98  1.04
      37     2.056  1.14  1.19  1.15  1.23  1.10  1.08  1.07  1.03  1.05  1.03
      39     1.931  5.12  2.83  1.98  1.56  1.30  1.14  1.02  0.93  0.86  0.81
      47     1.829  4.92  2.58  1.81  1.42  1.19  1.04  0.93  0.85  0.82  0.81
      43     1.739  5.62  2.95  2.06  1.69  1.42  1.24  1.11  1.01  0.94  0.88
      44     1.665  4.30  2.38  1.67  1.31  1.10  0.96  0.86  0.83  0.77  0.72
      53     1.596  4.56  2.39  1.67  1.32  1.10  1.00  0.90  0.82  0.76  0.72
      69     1.536  5.03  2.64  1.85  1.45  1.25  1.09  0.98  0.89  0.83  0.78
      51     1.390  4.95  2.70  1.97  1.55  1.30  1.13  1.01  0.93  0.86  0.84
      18     1.314  2.34  1.53  1.29  1.01  0.85  0.74  0.77  0.71  0.84  0.79
     774   overall  2.59  1.62  1.21  1.01  0.91  0.84  0.79  0.76  0.73  0.70


 List of     80 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   47    1    5    1.61   52.02  0.1414E+02  0.9242E+01 "alien"
   38   28    6    1.59   48.51  0.1319E+02  0.1347E+02 "alien"
   30    1   13    1.95   44.38 -0.1670E+02  0.1156E+02 "alien"
   45    2   10    1.57   42.48  0.1155E+02  0.1008E+02 "alien"
   21    9   21    1.62   41.86  0.1138E+02  0.8831E+01 "alien"
   31   25   13    1.63   39.59  0.1076E+02  0.8716E+01 "alien"
   21   15   21    1.57   36.59  0.9946E+01  0.9024E+01 "alien"
   39    1   11    1.73   34.39  0.1496E+02  0.1323E+02 "alien"
   34    5   17    1.60   34.07  0.9261E+01  0.7644E+01 "alien"
   37   25    2    1.73   32.71  0.1423E+02  0.1017E+02 "alien"
   26   21   14    1.77   32.63  0.1420E+02  0.7853E+01 "alien"
   27   22   14    1.73   31.59  0.1375E+02  0.7948E+01 "alien"
   28    2   14    1.95   30.48 -0.1147E+02  0.1119E+02 "alien"
   25   18   12    1.98   30.33 -0.1141E+02  0.1173E+02 "alien"
   33   22    5    1.89   29.75 -0.1119E+02  0.1470E+02 "alien"
   45    2    9    1.59   29.34  0.7976E+01  0.9653E+01 "alien"
   41    3   13    1.59   28.30  0.7693E+01  0.8765E+01 "alien"
   39   28    3    1.60   27.82  0.7561E+01  0.9369E+01 "alien"
   36    1    8    1.96   26.65 -0.1003E+02  0.8834E+01 "alien"
   33    1   12    1.90   26.40 -0.9934E+01  0.8153E+01 "alien"
   48    1    5    1.58   26.39  0.7174E+01  0.8996E+01 "alien"
   35   27   10    1.59   26.25  0.7134E+01  0.8929E+01 "alien"
   20   19   18    1.71   25.93  0.1128E+02  0.6977E+01 "alien"
   38    2    4    1.99   25.49 -0.9590E+01  0.1058E+02 "alien"
   36   24    3    1.77   25.30  0.1101E+02  0.1002E+02 "alien"
   20   13   21    1.60   25.09  0.6820E+01  0.8182E+01 "alien"
   23   21   16    1.73   24.83  0.1080E+02  0.1325E+02 "alien"
   46    1    7    1.61   24.50  0.6658E+01  0.8690E+01 "alien"
   14   10   18    1.94   24.25 -0.9125E+01  0.6245E+01 "alien"
   23   21   19    1.58   21.76  0.5914E+01  0.7020E+01 "alien"
   34   27   11    1.59   21.13  0.5743E+01  0.9635E+01 "alien"
   31   21    8    1.90   20.95 -0.7881E+01  0.8597E+01 "alien"
   54    3    2    1.43   20.51 -0.1942E+02  0.9902E+01 "alien"
   51    4   10    1.41   19.79  0.2494E+02  0.1584E+02 "alien"
   43    2   11    1.60   18.76  0.5099E+01  0.1035E+02 "alien"
   21   19   16    1.81   18.34  0.1931E+02  0.9829E+01 "alien"
   38    1    7    1.91   17.05 -0.6415E+01  0.8541E+01 "alien"
   35    5   17    1.58   17.04  0.4632E+01  0.1201E+02 "alien"
   32   27   14    1.54   16.86  0.1753E+02  0.1349E+02 "alien"
   40    3   13    1.62   16.83  0.4574E+01  0.8271E+01 "alien"
   29    2   14    1.92   16.59 -0.6242E+01  0.7660E+01 "alien"
   31   30   17    1.41   16.31  0.2056E+02  0.1033E+02 "alien"
   52    3    7    1.43   16.30 -0.1543E+02  0.1008E+02 "alien"
   19   16   19    1.70   16.23  0.7060E+01  0.5862E+01 "alien"
   40    1    5    1.88   16.19  0.1704E+02  0.1266E+02 "alien"
   26   21   21    1.44   15.96 -0.1511E+02  0.8389E+01 "alien"
   37    2   13    1.71   15.84  0.6894E+01  0.8895E+01 "alien"
   42   32    3    1.46   15.47 -0.1465E+02  0.1710E+02 "alien"
   40   30    4    1.53   15.35  0.1596E+02  0.1012E+02 "alien"
   46    3   11    1.51   15.08  0.1568E+02  0.9276E+01 "alien"
   21   18   15    1.90   14.88 -0.5599E+01  0.7460E+01 "alien"
   27   20   12    1.87   14.86  0.1565E+02  0.8234E+01 "alien"
   39   28    0    1.61   14.72  0.1111E+02  0.1352E+02 "alien"
   32   31   17    1.38   14.58  0.1837E+02  0.1804E+02 "alien"
   17   12   23    1.53   14.36  0.1494E+02  0.6924E+01 "alien"
   38   32   11    1.42   14.16 -0.1341E+02  0.9787E+01 "alien"
   47    2    9    1.54   14.15  0.1471E+02  0.1567E+02 "alien"
   20   18   16    1.85   14.05  0.1479E+02  0.9451E+01 "alien"
   30   21   10    1.85   13.84  0.1457E+02  0.1275E+02 "alien"
   44    1    0    1.76   13.76  0.1551E+02  0.1311E+02 "alien"
   40   32    7    1.46   13.57 -0.1285E+02  0.1040E+02 "alien"
   35   31   14    1.42   13.51 -0.1279E+02  0.1324E+02 "alien"
   25    9   21    1.56   13.37  0.1390E+02  0.9678E+01 "alien"
   36   25    5    1.72   13.18  0.5736E+01  0.8950E+01 "alien"
   48    7   15    1.36   12.95 -0.2742E+02  0.1441E+02 "alien"
   53    3    3    1.45   12.88 -0.1219E+02  0.1213E+02 "alien"
   44    5   14    1.48   12.83 -0.2260E+02  0.1432E+02 "alien"
   47    2    8    1.56   12.80  0.1331E+02  0.9396E+01 "alien"
   28   28   18    1.45   12.36 -0.2751E+02  0.1237E+02 "alien"
   38   29    8    1.53   12.10  0.1258E+02  0.9607E+01 "alien"
   34   22    0    1.91   11.81 -0.1008E+02  0.1285E+02 "alien"
   38    4   15    1.59   11.77  0.3200E+01  0.8101E+01 "alien"
   48    1    3    1.60   11.77  0.3200E+01  0.9219E+01 "alien"
   37   33   13    1.38   11.69  0.1473E+02  0.1032E+02 "alien"
   35   28   11    1.55   11.66  0.1212E+02  0.9337E+01 "alien"
   36    1   11    1.83   11.63  0.1224E+02  0.1580E+02 "alien"
   39   32   10    1.43   11.47 -0.1086E+02  0.1034E+02 "alien"
   42    0    6    1.77   10.96  0.1023E+02  0.1571E+02 "alien"
   36   25    4    1.74   10.42  0.4534E+01  0.9039E+01 "alien"
   36    2   13    1.74   10.40  0.4526E+01  0.7891E+01 "alien"


 List of     80 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   48    7   15    1.36   12.95 -0.2742E+02  0.1441E+02
   32   31   17    1.38   14.58  0.1837E+02  0.1804E+02
   37   33   13    1.38   11.69  0.1473E+02  0.1032E+02
   51    4   10    1.41   19.79  0.2494E+02  0.1584E+02
   31   30   17    1.41   16.31  0.2056E+02  0.1033E+02
   35   31   14    1.42   13.51 -0.1279E+02  0.1324E+02
   38   32   11    1.42   14.16 -0.1341E+02  0.9787E+01
   39   32   10    1.43   11.47 -0.1086E+02  0.1034E+02
   54    3    2    1.43   20.51 -0.1942E+02  0.9902E+01
   52    3    7    1.43   16.30 -0.1543E+02  0.1008E+02
   26   21   21    1.44   15.96 -0.1511E+02  0.8389E+01
   28   28   18    1.45   12.36 -0.2751E+02  0.1237E+02
   53    3    3    1.45   12.88 -0.1219E+02  0.1213E+02
   42   32    3    1.46   15.47 -0.1465E+02  0.1710E+02
   40   32    7    1.46   13.57 -0.1285E+02  0.1040E+02
   44    5   14    1.48   12.83 -0.2260E+02  0.1432E+02
   46    3   11    1.51   15.08  0.1568E+02  0.9276E+01
   40   30    4    1.53   15.35  0.1596E+02  0.1012E+02
   17   12   23    1.53   14.36  0.1494E+02  0.6924E+01
   38   29    8    1.53   12.10  0.1258E+02  0.9607E+01
   47    2    9    1.54   14.15  0.1471E+02  0.1567E+02
   32   27   14    1.54   16.86  0.1753E+02  0.1349E+02
   35   28   11    1.55   11.66  0.1212E+02  0.9337E+01
   47    2    8    1.56   12.80  0.1331E+02  0.9396E+01
   25    9   21    1.56   13.37  0.1390E+02  0.9678E+01
   45    2   10    1.57   42.48  0.1155E+02  0.1008E+02
   21   15   21    1.57   36.59  0.9946E+01  0.9024E+01
   23   21   19    1.58   21.76  0.5914E+01  0.7020E+01
   48    1    5    1.58   26.39  0.7174E+01  0.8996E+01
   35    5   17    1.58   17.04  0.4632E+01  0.1201E+02
   38   28    6    1.59   48.51  0.1319E+02  0.1347E+02
   34   27   11    1.59   21.13  0.5743E+01  0.9635E+01
   41    3   13    1.59   28.30  0.7693E+01  0.8765E+01
   38    4   15    1.59   11.77  0.3200E+01  0.8101E+01
   35   27   10    1.59   26.25  0.7134E+01  0.8929E+01
   45    2    9    1.59   29.34  0.7976E+01  0.9653E+01
   39   28    3    1.60   27.82  0.7561E+01  0.9369E+01
   48    1    3    1.60   11.77  0.3200E+01  0.9219E+01
   43    2   11    1.60   18.76  0.5099E+01  0.1035E+02
   34    5   17    1.60   34.07  0.9261E+01  0.7644E+01
   20   13   21    1.60   25.09  0.6820E+01  0.8182E+01
   46    1    7    1.61   24.50  0.6658E+01  0.8690E+01
   47    1    5    1.61   52.02  0.1414E+02  0.9242E+01
   39   28    0    1.61   14.72  0.1111E+02  0.1352E+02
   40    3   13    1.62   16.83  0.4574E+01  0.8271E+01
   21    9   21    1.62   41.86  0.1138E+02  0.8831E+01
   31   25   13    1.63   39.59  0.1076E+02  0.8716E+01
   19   16   19    1.70   16.23  0.7060E+01  0.5862E+01
   20   19   18    1.71   25.93  0.1128E+02  0.6977E+01
   37    2   13    1.71   15.84  0.6894E+01  0.8895E+01
   36   25    5    1.72   13.18  0.5736E+01  0.8950E+01
   37   25    2    1.73   32.71  0.1423E+02  0.1017E+02
   39    1   11    1.73   34.39  0.1496E+02  0.1323E+02
   27   22   14    1.73   31.59  0.1375E+02  0.7948E+01
   23   21   16    1.73   24.83  0.1080E+02  0.1325E+02
   36   25    4    1.74   10.42  0.4534E+01  0.9039E+01
   36    2   13    1.74   10.40  0.4526E+01  0.7891E+01
   44    1    0    1.76   13.76  0.1551E+02  0.1311E+02
   36   24    3    1.77   25.30  0.1101E+02  0.1002E+02
   42    0    6    1.77   10.96  0.1023E+02  0.1571E+02
   26   21   14    1.77   32.63  0.1420E+02  0.7853E+01
   21   19   16    1.81   18.34  0.1931E+02  0.9829E+01
   36    1   11    1.83   11.63  0.1224E+02  0.1580E+02
   20   18   16    1.85   14.05  0.1479E+02  0.9451E+01
   30   21   10    1.85   13.84  0.1457E+02  0.1275E+02
   27   20   12    1.87   14.86  0.1565E+02  0.8234E+01
   40    1    5    1.88   16.19  0.1704E+02  0.1266E+02
   33   22    5    1.89   29.75 -0.1119E+02  0.1470E+02
   33    1   12    1.90   26.40 -0.9934E+01  0.8153E+01
   21   18   15    1.90   14.88 -0.5599E+01  0.7460E+01
   31   21    8    1.90   20.95 -0.7881E+01  0.8597E+01
   38    1    7    1.91   17.05 -0.6415E+01  0.8541E+01
   34   22    0    1.91   11.81 -0.1008E+02  0.1285E+02
   29    2   14    1.92   16.59 -0.6242E+01  0.7660E+01
   14   10   18    1.94   24.25 -0.9125E+01  0.6245E+01
   30    1   13    1.95   44.38 -0.1670E+02  0.1156E+02
   28    2   14    1.95   30.48 -0.1147E+02  0.1119E+02
   36    1    8    1.96   26.65 -0.1003E+02  0.8834E+01
   25   18   12    1.98   30.33 -0.1141E+02  0.1173E+02
   38    2    4    1.99   25.49 -0.9590E+01  0.1058E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        33
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1458
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1365
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))       93


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
