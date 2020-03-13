
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)   7-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.307    77.307    38.685  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/demo_damaged/myfake_???.img SMV
 DATA_RANGE=       1      50
 ROTATION_AXIS=  0.000963  0.000022  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       1
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  0.999999 -0.000806  0.000651
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
 REFLECTING_RANGE_E.S.D.=   0.090
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

   91391 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       1 ...      50
       0 OVERLOADED REFLECTIONS (IGNORED)
      79 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
   91312 REFLECTIONS ACCEPTED



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
 UNIT_CELL_CONSTANTS=    77.00    77.00    38.00  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_B-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_C-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 REIDX=  0  0  0  0  0  0  0  0  0  0  0  0 (all 0 : not specified)
 TEST_RESOLUTION_RANGE=    8.00    4.50 (Angstrom) for space group determination
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas             2178
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
 *  33        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  34        mP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  32        oP          0.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.3  109.3   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.3  109.3   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
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
    41        mC        998.3     172.9   77.3   38.7  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    30        mC        998.3      77.3  172.9   38.7  90.0  90.0  63.4    0  1  0  0  2  1  0  0  0  0 -1  0
    40        oC        998.3      77.3  172.9   38.7  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    22        hP        998.3      77.3   77.3   38.7  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     1        cF        999.0     116.0  116.0  116.0  96.4  96.4 141.0   -1 -1  1  0  1  1  1  0 -1  1 -1  0
     2        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.6    0 -1  1  0 -1  0 -1  0  1 -1 -1  0
     3        cP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
     5        cI        999.0      86.4   86.4  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     4        hR        999.0      86.4   86.4  116.0 116.5  63.5 101.6    0  1 -1  0 -1  0  1  0  1  1  1  0
     6        tI        999.0     109.3   86.4   86.4  78.4  50.8  50.8   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     7        tI        999.0      86.4   86.4  109.3  50.8  50.8  78.4   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     8        oI        999.0      86.4   86.4  109.3  50.8  50.8  78.4    0  1  1  0  1  0  1  0  1  1  0  0
     9        hR        999.0      38.7   86.4  247.5 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
    10        mC        999.0      86.4   86.4   77.3  90.0  90.0 126.8    0  1  1  0  0 -1  1  0  1  0  0  0
    11        tP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    12        hP        999.0      38.7   77.3   77.3  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    13        oC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    15        tI        999.0      38.7   77.3  177.1  64.1  77.4  90.0    0  0 -1  0  0 -1  0  0 -2 -1 -1  0
    16        oF        999.0      86.4   86.4  177.1 107.0 119.2  53.2    0  1  1  0  0  1 -1  0 -2 -1 -1  0
    14        mC        999.0      86.4   86.4   77.3  90.0  90.0  53.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    17        mC        999.0      86.4   86.4   86.4  78.4 101.6  53.2    0  1 -1  0  0  1  1  0  1  0  1  0
    24        hR        999.0     177.1  109.3   38.7  90.0  77.4 108.0   -1 -2 -1  0 -1  1  0  0  0  0 -1  0
    42        oI        999.0      38.7   77.3  177.1 115.9 102.6  90.0    0  0  1  0  0  1  0  0 -2 -1 -1  0
    43        mI        999.0      86.4  177.1   77.3 115.9 153.4  60.8    0 -1 -1  0 -2 -1 -1  0  0  1  0  0

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

  *   0.00        0    5.0    1648    0.02  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.00    77.00    38.00  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.670   -23.432     0.023
 UNIT_CELL_B-AXIS=    23.432    73.670    -0.034
 UNIT_CELL_C-AXIS=    -0.006     0.020    38.685
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM     13515 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM AXIS ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.13
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.203    77.203    38.638  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  1.2E-02 1.2E-02 5.5E-03 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012953  0.012953  0.025881  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.572   -23.395     0.026
 COORDINATES OF UNIT CELL B-AXIS    23.395    73.572     0.007
 COORDINATES OF UNIT CELL C-AXIS    -0.013     0.000    38.638
 CRYSTAL MOSAICITY (DEGREES)     0.090
 LAB COORDINATES OF ROTATION AXIS -0.000330  0.000014  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006290  0.002014
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.88
 DETECTOR ORIGIN (PIXELS) AT                     1260.94   1242.03
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.03
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.20    77.20    38.64  90.000  90.000  90.000


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
 -0.05  0.05      8879        520.240       0.915
  0.05  0.15      9153        498.081       0.955
  0.15  0.25      9215        460.675       1.033
  0.25  0.35      9114        443.489       1.073
  0.35  0.45      9114        469.560       1.013
  0.45  0.55      9233        472.315       1.007
  0.55  0.65      9031        458.637       1.037
  0.65  0.75      9197        470.464       1.011
  0.75  0.85      9257        467.419       1.018
  0.85  0.95      9119        498.094       0.955
 
 -0.05  0.95     91312        475.770       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED      200
 DEGREES OF FREEDOM OF CHI^2 FIT             24292.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.630
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     0.1 XMAX=    49.9 NXBIN=   10
 YMIN= 0.00067 YMAX= 0.60124 NYBIN=   20
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS      39040


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED      784
 DEGREES OF FREEDOM OF CHI^2 FIT             24267.2
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.627
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     6.0 XMAX=  2457.3 NXBIN=   28
 YMIN=     5.7 YMAX=  2521.3 NYBIN=   28
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS      39040


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED      130
 DEGREES OF FREEDOM OF CHI^2 FIT             24292.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.625
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     0.1 XMAX=    49.9 NXBIN=   10
 DETECTOR_SURFACE_POSITION=    1232    1264
 DETECTOR_SURFACE_POSITION=    1648    1691
 DETECTOR_SURFACE_POSITION=     815    1691
 DETECTOR_SURFACE_POSITION=     815     836
 DETECTOR_SURFACE_POSITION=    1648     836
 DETECTOR_SURFACE_POSITION=    2174    1664
 DETECTOR_SURFACE_POSITION=    1622    2230
 DETECTOR_SURFACE_POSITION=     841    2230
 DETECTOR_SURFACE_POSITION=     289    1664
 DETECTOR_SURFACE_POSITION=     289     863
 DETECTOR_SURFACE_POSITION=     841     297
 DETECTOR_SURFACE_POSITION=    1622     297
 DETECTOR_SURFACE_POSITION=    2174     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS      39040


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
 7.982E-01  1.812E-03   26.30


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

   38.601  12.238    26.11   1.11      2.88      2.54      88      91      32
   12.238   8.880    24.88   1.65      3.30      2.66     199     201      40
    8.880   7.315    24.73   1.63      3.20      2.57     219     225      31
    7.315   6.364    24.33   1.24      2.85      2.64     266     268      55
    6.364   5.707    24.22   1.42      3.04      2.63     293     294      53
    5.707   5.220    24.29   1.09      2.30      2.68     327     331      60
    5.220   4.839    24.23   1.22      2.64      2.68     364     369      60
    4.839   4.531    24.46   1.17      2.52      2.70     411     415      66
    4.531   4.275    24.93   0.90      2.10      2.70     475     480      50
    4.275   4.058    24.52   1.08      2.43      2.68     493     497      54
    4.058   3.871    24.15   1.01      2.31      2.69     497     502      64
    3.871   3.708    23.74   1.15      2.67      2.63     457     470      62
    3.708   3.564    24.63   0.95      2.18      2.64     564     569      39
    3.564   3.435    22.52   1.00      2.27      2.66     500     506      64
    3.435   3.319    23.94   0.98      2.22      2.68     523     530      66
    3.319   3.215    23.51   1.12      2.31      2.64     482     489      56
    3.215   3.119    21.98   1.13      2.56      2.74     597     614      72
    3.119   3.032    22.39   1.16      2.75      2.75     591     601     106
    3.032   2.952    22.37   1.00      2.32      2.70     552     562      85
    2.952   2.877    20.16   1.11      2.71      2.83     621     632     101
    2.877   2.808    21.78   1.07      2.75      2.81     672     679     104
    2.808   2.744    20.95   1.13      2.84      2.84     626     634     124
    2.744   2.684    20.79   1.18      2.72      2.84     704     715      89
    2.684   2.628    20.57   1.14      2.84      2.90     726     736      88
    2.628   2.575    19.83   1.26      2.89      2.92     692     700     125
    2.575   2.525    20.24   1.22      3.03      2.85     655     669     125
    2.525   2.478    10.21   1.10      3.51      3.48     788     792      55
    2.478   2.434    -0.00   1.01  40380.42  39342.81     908     912       5
    2.434   2.392    -0.01   0.97    -99.90    -99.90     836     842       2
    2.392   2.352    -0.04   0.95    -99.90    -99.90     936     942       3
    2.352   2.313     0.01   0.86   2511.75   2719.70     832     836       2
    2.313   2.277    -0.01   0.97    -99.90    -99.90     840     842       1
    2.277   2.242     0.07   0.96    739.31    753.45     966     971       3
    2.242   2.209    -0.03   1.03    -99.90    -99.90     862     869       2
    2.209   2.178    -0.01   1.00    -99.90    -99.90    1034    1038       5
    2.178   2.147     0.06   0.91    757.96    811.95     965     970       1
    2.147   2.118    -0.01   1.01    -99.90    -99.90     971     983       2
    2.118   2.090    -0.05   0.98    -99.90    -99.90     989     999       6
    2.090   2.063     0.00   1.01   2966.50   2977.72    1008    1013       6
    2.063   2.037     0.02   1.05   4956.42   4761.83    1029    1037       6
    2.037   2.012    -0.04   0.95    -99.90    -99.90     966     983       4
    2.012   1.988     0.01   0.91   6746.96   6996.53    1030    1038       0
    1.988   1.965    -0.04   0.98    -99.90    -99.90     918     939       6
    1.965   1.943    -0.05   0.92    -99.90    -99.90    1018    1032       4
    1.943   1.921    -0.05   0.99    -99.90    -99.90    1071    1086       6
    1.921   1.900    -0.06   0.96    -99.90    -99.90    1134    1148       1
    1.900   1.880     0.02   0.91   3270.69   3459.05    1064    1072       3
    1.880   1.860    -0.02   0.87    -99.90    -99.90    1114    1121       3
    1.860   1.841    -0.07   0.89    -99.90    -99.90     947     958       2
    1.841   1.823    -0.01   0.95    -99.90    -99.90    1161    1163       3
    1.823   1.805     0.00   0.97   6653.05   6602.75    1099    1106       3
    1.805   1.788    -0.01   0.92    -99.90    -99.90    1129    1131       8
    1.788   1.771    -0.01   0.93    -99.90    -99.90    1189    1192       0
    1.771   1.754    -0.01   0.88  27050.41  28698.85    1122    1129       3
    1.754   1.738     0.04   1.00   1383.11   1365.08    1289    1292       4
    1.738   1.723     0.01   0.98   3614.12   3604.34    1300    1303       3
    1.723   1.707     0.00   0.91   9295.38   9584.07    1111    1118       2
    1.707   1.693     0.01   0.94   2815.55   2937.84    1171    1180       7
    1.693   1.678    -0.03   0.91    -99.90    -99.90    1244    1246       3
    1.678   1.664     0.01   0.91   3590.86   3704.99    1322    1333       1
    1.664   1.651     0.01   0.91   3138.65   3263.52    1193    1206       2
    1.651   1.637    -0.01   0.87    -99.90    -99.90    1155    1163       6
    1.637   1.624    -0.06   0.96    -99.90    -99.90    1222    1225       4
    1.624   1.612     0.00   0.91   9442.07   9781.90    1295    1313       0
    1.612   1.599    -0.01   0.98    -99.90    -99.90    1261    1275       5
    1.599   1.587    -0.03   1.02    -99.90    -99.90    1185    1190       1
    1.587   1.575     0.08   0.92    750.07    775.74    1229    1236       1
    1.575   1.564     0.01   0.94   6033.25   6026.70    1198    1206       4
    1.564   1.552    -0.03   0.95    -99.90    -99.90    1290    1293       0
    1.552   1.541    -0.02   0.92    -99.90    -99.90    1237    1241       0
    1.541   1.530     0.00   0.96   3528.15   3581.64    1383    1383       1
    1.530   1.520    -0.02   0.92    -99.90    -99.90    1342    1345       4
    1.520   1.509    -0.03   0.89    -99.90    -99.90    1418    1421       2
    1.509   1.499     0.04   0.88   1359.98   1461.08    1294    1308       4
    1.499   1.489     0.04   0.92   1422.59   1487.68    1317    1323       5
    1.489   1.479     0.02   0.97   2744.38   2731.09    1449    1453       4
    1.479   1.469    -0.01   0.98    -99.90    -99.90    1466    1467       1
    1.469   1.460     0.05   0.93   1176.41   1206.81    1375    1375       3
    1.460   1.451     0.00   0.91   4509.10   4674.54    1340    1352       1
    1.451   1.442     0.01   1.01   2941.50   2876.35    1306    1328       1
    1.442   1.433     0.01   0.99   5960.30   5884.47    1236    1255       6
    1.433   1.424     0.00   0.94   5759.72   5842.52    1230    1247       2
    1.424   1.415    -0.00   0.99   4995.48   4899.65    1177    1196       4
    1.415   1.407     0.01   0.98   3753.52   3727.01    1232    1255       2
    1.407   1.399    -0.08   0.99    -99.90    -99.90    1039    1059      19
    1.399   1.391    -0.07   0.86    -99.90    -99.90     905     909      10
    1.391   1.383    -0.05   0.95    -99.90    -99.90     982     983      11
    1.383   1.375    -0.00   0.94  46565.11  48136.45     945     951      17
    1.375   1.367    -0.04   0.91    -99.90    -99.90     803     812       9
    1.367   1.359    -0.03   0.94    -99.90    -99.90     779     783      12
    1.359   1.352     0.07   0.99    856.13    852.99     718     725       6
    1.352   1.345    -0.03   0.90    -99.90    -99.90     655     680       8
    1.345   1.337    -0.04   0.96    -99.90    -99.90     540     574       7
    1.337   1.330     0.04   1.01   1983.69   1946.47     494     536       0
    1.330   1.323     0.08   0.97    539.40    537.53     426     468      15
    1.323   1.316    -0.02   0.80    -99.90    -99.90     377     392       2
    1.316   1.309     0.27   0.90    177.44    178.33     233     262       9
    1.309   1.303     0.06   0.89    743.46    780.14     134     181       2
    1.303   1.296     0.21   0.80    336.82    384.84     100     151       2
    1.296   1.290     0.00   1.01    416.19    381.52      74     123       0
   --------------------------------------------------------------------------
   38.601   1.290     3.33   0.97      4.11      4.27   87991   89040    2228
 


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
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1      11
 ******************************************************************************


 AVERAGE INTENSITY FOR      13 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     1.1% OF MEAN INTENSITY


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

     3.85         837     489      1260       38.8%       2.4%      2.2%      683   31.40      3.3%    99.7*     0    1.566       5
     2.73        1406     784      2111       37.1%       2.5%      2.3%     1169   28.90      3.4%    99.7*    68*   1.220      33
     2.23        2032     974      2693       36.2%       4.0%      3.5%     1916    9.94      5.4%    99.8*    11    0.781      95
     1.93        2443    1216      3113       39.1%     -99.9%    -99.9%     2309    0.00    -99.9%    -0.9     14    0.485      69
     1.73        2742    1395      3521       39.6%     -99.9%    -99.9%     2607    0.00    -99.9%    -0.3     16    0.560      40
     1.58        3014    1569      3881       40.4%    3055.1%   3104.4%     2827    0.01   4291.5%    -3.4    -61    0.360      27
     1.46        3268    1676      4179       40.1%    1262.4%   1294.1%     3092    0.04   1767.2%    -1.6     -5    0.416      43
     1.37        2754    1396      4492       31.1%    7766.4%   7910.3%     2629    0.00  10858.2%     1.1    -12    0.411      41
     1.29        1077     674      4733       14.2%     893.4%    938.3%      794    0.07   1253.5%    41.3*   -81    0.476       6
    total       19573   10173     29983       33.9%       3.8%      3.6%    18026    4.70      5.3%    99.9*     5    0.628     359


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     20053
 NUMBER OF REJECTED MISFITS                               429
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                   13
 NUMBER OF ACCEPTED OBSERVATIONS                        19611
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                  10193



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1      21
 ******************************************************************************


 AVERAGE INTENSITY FOR      23 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     0.8% OF MEAN INTENSITY


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

     3.85        1593     658      1260       52.2%       2.4%      2.5%     1468   36.47      3.0%    99.8*    72*   1.505     119
     2.73        2691    1209      2111       57.3%       2.5%      2.5%     2425   32.14      3.2%    99.8*    79*   1.453     221
     2.23        3881    1661      2693       61.7%       4.0%      3.7%     3729   10.42      5.3%    99.8*    16    0.740     324
     1.93        4640    2031      3113       65.2%     -99.9%    -99.9%     4460    0.00    -99.9%     0.2      6    0.502     354
     1.73        5262    2309      3521       65.6%     -99.9%    -99.9%     5094    0.00    -99.9%     0.9     -1    0.505     367
     1.58        5767    2547      3881       65.6%    8576.8%   8828.6%     5534   -0.01  11362.1%    -0.8      2    0.541     407
     1.46        6247    2716      4179       65.0%    1656.7%   1697.9%     6037    0.03   2182.5%     0.0      2    0.531     451
     1.37        5274    2310      4492       51.4%   39322.1%  39383.5%     5112   -0.01  52089.6%    -0.1     -1    0.512     360
     1.29        2082    1154      4733       24.4%    1263.9%   1344.0%     1665    0.06   1699.5%    32.4*    -7    0.520      54
    total       37437   16595     29983       55.3%       3.9%      3.9%    35524    4.83      5.0%    99.9*    11*   0.668    2657


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     38391
 NUMBER OF REJECTED MISFITS                               861
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                   23
 NUMBER OF ACCEPTED OBSERVATIONS                        37507
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                  16629



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1      30
 ******************************************************************************


 AVERAGE INTENSITY FOR      35 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     0.6% OF MEAN INTENSITY


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

     3.85        2281     683      1260       54.2%       2.5%      2.6%     2221   42.82      3.0%    99.9*    75*   1.531     220
     2.73        3845    1276      2111       60.4%       2.4%      2.6%     3708   37.35      2.9%    99.9*    82*   1.407     474
     2.23        5549    1754      2693       65.1%       4.1%      3.9%     5462   12.08      4.9%    99.8*    15*   0.770     891
     1.93        6611    2116      3113       68.0%     -99.9%    -99.9%     6528    0.00    -99.9%    -2.6      3    0.514    1140
     1.73        7538    2415      3521       68.6%     -99.9%    -99.9%     7453    0.00    -99.9%     2.6     -2    0.495    1267
     1.58        8227    2657      3881       68.5%   14910.2%  15369.7%     8091   -0.01  18093.7%     1.2      0    0.509    1362
     1.46        8912    2831      4179       67.7%    2127.9%   2173.9%     8799    0.03   2574.9%     1.3      3    0.537    1475
     1.37        7558    2553      4492       56.8%   19925.8%  19959.6%     7421   -0.02  24464.3%    -1.5     -2    0.520    1121
     1.29        2966    1353      4733       28.6%    1494.0%   1551.8%     2559    0.08   1891.1%    28.7*    -5    0.536     206
    total       53487   17638     29983       58.8%       4.0%      4.1%    52242    5.56      4.8%    99.9*     7    0.623    8156


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     54879
 NUMBER OF REJECTED MISFITS                              1241
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                   35
 NUMBER OF ACCEPTED OBSERVATIONS                        53603
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                  17680



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1      40
 ******************************************************************************


 AVERAGE INTENSITY FOR      44 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     0.2% OF MEAN INTENSITY


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

     3.85        3024     778      1260       61.7%       2.6%      2.7%     2992   46.45      2.9%    99.9*    78*   1.600     323
     2.73        5114    1443      2111       68.4%       2.5%      2.7%     5052   40.63      2.9%    99.9*    80*   1.485     656
     2.23        7373    2014      2693       74.8%       4.1%      4.0%     7322   12.98      4.8%    99.9*    16*   0.794    1195
     1.93        8826    2460      3113       79.0%     -99.9%    -99.9%     8758    0.00    -99.9%    -0.9      2    0.507    1569
     1.73       10061    2821      3521       80.1%     -99.9%    -99.9%     9995    0.00    -99.9%     3.5     -4    0.496    1804
     1.58       11018    3137      3881       80.8%   42146.4%  43361.1%    10916   -0.02  49842.2%     3.2      0    0.509    2023
     1.46       11939    3393      4179       81.2%    2206.7%   2264.4%    11867    0.03   2616.5%     2.2      2    0.520    2188
     1.37       10046    3119      4492       69.4%     -99.9%    -99.9%     9879    0.00    -99.9%    -2.1     -1    0.508    1625
     1.29        3971    1594      4733       33.7%    1825.8%   1881.9%     3602    0.09   2257.5%    23.1*    -3    0.540     350
    total       71372   20759     29983       69.2%       4.1%      4.2%    70383    5.82      4.8%    99.9*     6    0.623   11733


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     73277
 NUMBER OF REJECTED MISFITS                              1710
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                   44
 NUMBER OF ACCEPTED OBSERVATIONS                        71523
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                  20811



 ******************************************************************************
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       1      50)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

    2    0    0    38.601   0.2254E+04  0.8600E+02        26.21           1 
    3    0    0    25.734  -0.2913E+01  0.9998E+00        -2.91           2*
    4    0    0    19.301   0.1176E+05  0.3169E+03        37.12           2 
    5    0    0    15.441   0.2455E+01  0.1702E+01         1.44           2*
    7    0    0    11.029   0.2356E+00  0.1724E+01         0.14           2*
    8    0    0     9.650   0.4550E+03  0.1342E+02        33.92           2 
    9    0    0     8.578   0.2154E+01  0.2263E+01         0.95           2*
   10    0    0     7.720   0.2386E+04  0.9224E+02        25.87           1 
   11    0    0     7.018   0.8604E+01  0.4043E+01         2.13           2*
   12    0    0     6.434   0.1127E+03  0.5350E+01        21.06           2 
   13    0    0     5.939  -0.2099E+01  0.6700E+01        -0.31           1*
   14    0    0     5.514   0.3937E+04  0.1074E+03        36.65           2 
   15    0    0     5.147   0.3689E+01  0.4758E+01         0.78           2*
   16    0    0     4.825   0.3375E+04  0.9276E+02        36.38           2 
   17    0    0     4.541  -0.8473E+01  0.5441E+01        -1.56           2*
   18    0    0     4.289   0.1922E+05  0.5186E+03        37.06           2 
   19    0    0     4.063  -0.5940E+01  0.7926E+01        -0.75           2*
   20    0    0     3.860   0.1015E+05  0.3888E+03        26.10           1 
   21    0    0     3.676  -0.2622E+01  0.7068E+01        -0.37           2*
   22    0    0     3.509   0.2708E+04  0.7547E+02        35.88           2 
   23    0    0     3.357  -0.3857E+01  0.8874E+01        -0.43           1*
   24    0    0     3.217   0.2049E+04  0.5837E+02        35.10           2 
   25    0    0     3.088   0.4149E+01  0.6703E+01         0.62           2*
   26    0    0     2.969   0.5927E+04  0.1621E+03        36.56           2 
   27    0    0     2.859   0.8370E+01  0.9882E+01         0.85           2*
   28    0    0     2.757   0.1349E+04  0.3980E+02        33.90           2 
   29    0    0     2.662  -0.5603E+01  0.9014E+01        -0.62           2*
   30    0    0     2.573   0.9431E+01  0.1367E+02         0.69           1 
   31    0    0     2.490   0.6262E+01  0.7366E+01         0.85           2*
   32    0    0     2.413  -0.1815E+02  0.1041E+02        -1.74           2 
   33    0    0     2.339  -0.9603E+01  0.1046E+02        -0.92           1*
   34    0    0     2.271  -0.5479E+01  0.7694E+01        -0.71           2 
   35    0    0     2.206   0.1756E+01  0.1116E+02         0.16           2*
   36    0    0     2.145   0.3443E+01  0.8018E+01         0.43           2 
   37    0    0     2.087   0.2767E+01  0.1108E+02         0.25           2*
   38    0    0     2.032   0.3964E+01  0.9864E+01         0.40           2 
   39    0    0     1.980   0.2853E+01  0.1118E+02         0.26           1*
   40    0    0     1.930  -0.1097E+01  0.1104E+02        -0.10           2 
   41    0    0     1.883  -0.5780E+01  0.7667E+01        -0.75           2*
   42    0    0     1.838   0.4168E+01  0.1022E+02         0.41           2 
   43    0    0     1.795   0.4548E+01  0.1037E+02         0.44           2*
   44    0    0     1.755  -0.6958E+01  0.7488E+01        -0.93           2 
   45    0    0     1.716   0.1413E+02  0.1064E+02         1.33           2*
   46    0    0     1.678   0.2572E+01  0.7723E+01         0.33           2 
   47    0    0     1.643   0.2861E+01  0.7923E+01         0.36           2*
   48    0    0     1.608   0.2843E+00  0.1605E+02         0.02           1 
   49    0    0     1.576  -0.1069E+02  0.1156E+02        -0.93           1*
   50    0    0     1.544  -0.5700E+01  0.1594E+02        -0.36           1 
   51    0    0     1.514  -0.1109E+02  0.1461E+02        -0.76           1*
   52    0    0     1.485  -0.3290E+01  0.1227E+02        -0.27           1 

 AVERAGE INTENSITY FOR      44 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS     0.1% OF MEAN INTENSITY


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

     3.85        3735     967      1260       76.7%       2.6%      2.7%     3692   46.27      3.0%    99.9*    78*   1.627     369
     2.73        6356    1759      2111       83.3%       2.5%      2.7%     6256   40.99      2.9%    99.9*    81*   1.576     706
     2.23        9135    2404      2693       89.3%       4.0%      4.0%     9060   13.15      4.7%    99.9*    16*   0.810    1314
     1.93       11053    2923      3113       93.9%     -99.9%    -99.9%    10934    0.00    -99.9%     2.1      2    0.519    1763
     1.73       12552    3219      3521       91.4%     -99.9%    -99.9%    12483    0.00    -99.9%     4.0     -1    0.500    2141
     1.58       13721    3524      3881       90.8%   16632.3%  17103.7%    13621   -0.02  19272.8%     4.7      0    0.522    2421
     1.46       14856    3747      4179       89.7%    3883.6%   3977.0%    14810    0.00   4500.9%     2.2      1    0.520    2678
     1.37       12518    3467      4492       77.2%   26465.0%  26745.6%    12369   -0.02  31109.3%     1.4      0    0.510    2058
     1.29        4942    1758      4733       37.1%    1229.3%   1265.0%     4603    0.08   1485.2%    45.1*    -3    0.537     580
    total       88868   23768     29983       79.3%       4.1%      4.3%    87828    6.24      4.8%    99.9*     6    0.626   14030


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES     91312
 NUMBER OF REJECTED MISFITS                              2228
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                   44
 NUMBER OF ACCEPTED OBSERVATIONS                        89040
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                  23821



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

   WILSON LINE (using all data) : A=   7.465 B=  35.956 CORRELATION=  0.75
      #      RES      SS        <I>       log(<I>)       BO
     320     7.296   0.005  5.7563E+03       8.658    -127.0
     571     4.518   0.012  4.7162E+03       8.459     -40.6
     763     3.530   0.020  3.6221E+03       8.195     -18.2
     867     2.990   0.028  2.0844E+03       7.642      -3.2
    1017     2.642   0.036  1.3115E+03       7.179       4.0
    1186     2.392   0.044  3.3871E+01       3.523      45.1
    1281     2.200   0.052  4.0271E-01      -0.910      81.1
    1390     2.050   0.059 -1.4066E-02    -999.999   -1000.0
    1471     1.927   0.067 -2.7665E-01    -999.999   -1000.0
    1527     1.823   0.075  1.2510E-02      -4.381      78.7
    1587     1.734   0.083  2.0254E-01      -1.597      54.5
    1696     1.657   0.091 -1.9222E-01    -999.999   -1000.0
    1704     1.589   0.099  1.9201E-01      -1.650      46.0
    1749     1.529   0.107 -5.6220E-02    -999.999   -1000.0
    1819     1.476   0.115  3.2518E-01      -1.123      37.4
    1859     1.428   0.123  2.2671E-01      -1.484      36.5
    1322     1.385   0.130 -3.6305E-01    -999.999   -1000.0
    1029     1.345   0.138  2.2680E+00       0.819      24.0
     614     1.309   0.146  2.8822E+00       1.059      22.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

     117     7.296        0.941        0.792        0.653
     141     4.518        0.539        0.293        0.131
     135     3.530        1.382        1.987        2.849
     141     2.990        0.715        0.607        0.482
     152     2.642        0.779        0.743        0.653
     153     2.392        2.571       13.810       57.429
     164     2.200      142.248       84.915     7486.948
     159     2.050       28.954       34.941      258.497
     170     1.927       23.095        4.818      215.569
     153     1.823   145331.720  9491958.686*************
     151     1.734      490.296        3.376    65895.660
     137     1.428      750.061     1530.367   168714.537
      46     1.385       34.039       57.584      401.889
    2483   overall     9109.532   584933.773392775749.257


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

     203     7.296        0.978        0.984        0.954
     430     4.518        1.137        1.160        1.102
     628     3.530        0.896        0.750        0.580
     726     2.990        1.116        1.355        1.910
     865     2.642        1.058        1.132        1.317
    1033     2.392       36.816     1259.396    40801.944
    1117     2.200      114.704      230.565     7734.659
    1436     1.734      439.621      674.274   118971.928
    1558     1.589      621.600      408.649   281373.694
    1667     1.476      216.114      441.490    29964.311
    1722     1.428      580.267      578.311   223436.771
     991     1.345      179.423    22111.674  2198437.432
     599     1.309       83.920     5968.762   394774.487
   21289   overall     8148.524   -22907.603737640672.273

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

     117     7.296  1.10  0.96  0.84  0.81  0.85  0.87  0.87  0.88  0.90  0.91
     141     4.518  1.11  1.15  1.06  1.12  1.12  1.09  1.09  1.07  1.08  1.07
     135     3.530  0.98  1.03  0.96  1.03  0.97  0.95  0.93  0.97  0.98  0.99
     141     2.990  1.11  1.07  1.02  1.05  1.08  1.07  1.07  1.06  1.07  1.05
     152     2.642  1.06  1.08  1.14  1.14  1.10  1.11  1.11  1.10  1.10  1.09
     153     2.392  3.00  2.42  2.20  1.98  1.83  1.71  1.62  1.55  1.48  1.43
     164     2.200  2.09  1.50  1.24  1.10  1.00  0.92  0.87  0.82  0.79  0.76
     159     2.050  2.08  1.49  1.24  1.09  0.99  0.92  0.86  0.82  0.79  0.76
     170     1.927  2.04  1.46  1.21  1.07  0.98  0.91  0.86  0.81  0.78  0.75
     153     1.823  1.97  1.42  1.18  1.04  0.94  0.87  0.82  0.78  0.75  0.72
     151     1.734  1.89  1.36  1.14  1.01  0.92  0.85  0.81  0.77  0.75  0.72
     137     1.428  1.73  1.25  1.03  0.91  0.83  0.77  0.72  0.69  0.66  0.63
      46     1.385  1.75  1.26  1.04  0.92  0.84  0.77  0.73  0.69  0.66  0.64
    2483   overall  1.26  1.00  0.88  0.82  0.77  0.73  0.71  0.69  0.67  0.66


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

     203     7.296  1.14  1.03  1.10  1.14  1.08  1.09  1.09  1.04  1.04  1.00
     430     4.518  0.93  0.86  1.01  0.95  0.96  0.96  0.98  0.98  0.98  0.99
     628     3.530  1.10  1.07  1.02  0.97  0.96  0.98  0.96  0.98  0.99  0.99
     726     2.990  1.24  1.12  1.02  0.99  0.97  0.97  0.96  0.97  0.97  0.96
     865     2.642  0.96  0.92  0.96  0.99  0.99  0.96  0.96  0.96  0.96  0.97
    1033     2.392  7.35  4.70  3.53  2.88  2.43  2.13  1.91  1.75  1.62  1.53
    1117     2.200  5.04  2.67  1.89  1.49  1.26  1.10  1.00  0.92  0.86  0.81
    1436     1.734  5.20  2.74  1.93  1.52  1.28  1.12  1.00  0.92  0.86  0.81
    1558     1.589  5.00  2.63  1.84  1.45  1.22  1.07  0.96  0.89  0.83  0.78
    1667     1.476  4.93  2.59  1.82  1.44  1.21  1.06  0.95  0.87  0.81  0.77
    1722     1.428  5.14  2.70  1.89  1.49  1.25  1.09  0.99  0.90  0.84  0.79
     991     1.345  5.15  2.78  2.01  1.62  1.40  1.24  1.14  1.06  1.01  0.97
     599     1.309  5.00  2.71  1.93  1.59  1.38  1.24  1.13  1.05  1.00  0.95
   21289   overall  2.66  1.51  1.11  0.91  0.79  0.71  0.65  0.61  0.58  0.55


 List of    500 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   32   24    5    1.87 2426.75  0.3036E+02  0.1582E+02 "alien"
   36   20    4    1.84 2210.61  0.2765E+02  0.1128E+02 "alien"
   30   21    3    2.08 2039.50 -0.2869E+02  0.1123E+02 "alien"
   30   21    6    2.00 2007.14 -0.2823E+02  0.1028E+02 "alien"
   35   24    3    1.80 1908.61  0.2388E+02  0.1021E+02 "alien"
   31    9   13    1.86 1634.08  0.2044E+02  0.9028E+01 "alien"
   28   22   11    1.84 1559.52  0.1951E+02  0.8034E+01 "alien"
   33   14    5    2.07 1553.56 -0.2185E+02  0.7975E+01 "alien"
   33   20    8    1.85 1525.54  0.1908E+02  0.9715E+01 "alien"
   38    7    8    1.85 1517.00  0.1898E+02  0.7528E+01 "alien"
   26   14   12    2.03 1511.82 -0.2127E+02  0.8900E+01 "alien"
   36   15    7    1.86 1462.57  0.1830E+02  0.1057E+02 "alien"
   34   13    3    2.09 1430.39 -0.2012E+02  0.8101E+01 "alien"
   30   23    2    2.03 1399.92 -0.1969E+02  0.1110E+02 "alien"
   33   13    5    2.10 1393.33 -0.1960E+02  0.1069E+02 "alien"
   31   21    2    2.05 1381.76 -0.1944E+02  0.1140E+02 "alien"
   37    4   10    1.83 1366.96  0.1710E+02  0.5623E+01 "alien"
   23    5   15    2.03 1356.19 -0.1908E+02  0.7080E+01 "alien"
   31   21   10    1.82 1347.13  0.1685E+02  0.7110E+01 "alien"
   34   12    5    2.06 1341.69 -0.1887E+02  0.1129E+02 "alien"
   29    6   15    1.83 1334.22  0.1669E+02  0.8253E+01 "alien"
   34    3   12    1.85 1307.70  0.1636E+02  0.4912E+01 "alien"
   32   18    5    2.03 1278.74 -0.1799E+02  0.7873E+01 "alien"
   31   25    7    1.83 1269.55  0.1588E+02  0.1072E+02 "alien"
   26   23    7    2.06 1268.85 -0.1785E+02  0.7643E+01 "alien"
   38   13    5    1.87 1265.51  0.1583E+02  0.8963E+01 "alien"
   38   17    4    1.82 1260.29  0.1577E+02  0.1101E+02 "alien"
   28   27    8    1.84 1237.64  0.1548E+02  0.7022E+01 "alien"
   25   17   11    2.07 1235.01 -0.1737E+02  0.6549E+01 "alien"
   32   23    9    1.78 1228.23  0.1536E+02  0.6096E+01 "alien"
   34   25    4    1.80 1228.18  0.1536E+02  0.1067E+02 "alien"
   27   24    5    2.06 1189.21 -0.1673E+02  0.1557E+02 "alien"
   32   25    5    1.85 1187.68  0.1486E+02  0.1062E+02 "alien"
   40   12    1    1.85 1185.09  0.1483E+02  0.1090E+02 "alien"
   32   21    2    2.01 1162.37 -0.1635E+02  0.1048E+02 "alien"
   28    9   16    1.78 1146.66  0.1434E+02  0.7656E+01 "alien"
   30   23   10    1.81 1137.08  0.1422E+02  0.9079E+01 "alien"
   37   14    7    1.84 1131.20  0.1415E+02  0.7689E+01 "alien"
   30   16   12    1.86 1130.22  0.1414E+02  0.6468E+01 "alien"
   28   24    5    2.02 1129.28 -0.1588E+02  0.1478E+02 "alien"
   29    5   11    2.10 1125.39 -0.1583E+02  0.6802E+01 "alien"
   26   24    6    2.07 1124.42 -0.1582E+02  0.1078E+02 "alien"
   40   12    4    1.82 1115.56  0.1396E+02  0.7924E+01 "alien"
   35    9   11    1.83 1110.12  0.1389E+02  0.6428E+01 "alien"
   29   20    6    2.07 1107.05 -0.1557E+02  0.1084E+02 "alien"
   24   20   11    2.02 1096.79 -0.1543E+02  0.9157E+01 "alien"
   18    4   19    1.83 1077.50  0.1348E+02  0.5032E+01 "alien"
   34   14    3    2.07 1072.00 -0.1508E+02  0.7306E+01 "alien"
   28   21    8    2.01 1063.23 -0.1496E+02  0.1047E+02 "alien"
   30   22    1    2.07 1053.81 -0.1482E+02  0.8301E+01 "alien"
   38   11    8    1.81 1043.67  0.1306E+02  0.7019E+01 "alien"
   29   21    7    2.01 1032.65 -0.1453E+02  0.7617E+01 "alien"
   31   17   11    1.85 1024.25  0.1281E+02  0.6442E+01 "alien"
   26   18    9    2.12 1017.51 -0.1431E+02  0.1016E+02 "alien"
   23   15   12    2.12 1005.75 -0.1415E+02  0.7252E+01 "alien"
   31   12    9    2.04 1004.90 -0.1413E+02  0.7419E+01 "alien"
   31   19    2    2.11  999.65 -0.1406E+02  0.6571E+01 "alien"
   33   14   12    1.79  998.76  0.1249E+02  0.6111E+01 "alien"
   39   16    4    1.80  998.73  0.1249E+02  0.7716E+01 "alien"
   36   19    7    1.79  993.55  0.1243E+02  0.6215E+01 "alien"
   29   16    9    2.05  974.64 -0.1371E+02  0.6854E+01 "alien"
   31   20   10    1.84  970.99  0.1215E+02  0.6829E+01 "alien"
   30   10    9    2.12  965.37 -0.1358E+02  0.5918E+01 "alien"
   30   20    7    2.00  962.11 -0.1353E+02  0.1056E+02 "alien"
   32   21    9    1.83  948.89  0.1187E+02  0.8283E+01 "alien"
   32   12    8    2.05  945.78 -0.1330E+02  0.7486E+01 "alien"
   30   14    9    2.05  934.70 -0.1315E+02  0.7295E+01 "alien"
   23   22   14    1.82  932.09  0.1166E+02  0.5804E+01 "alien"
   39   14    1    1.86  927.64  0.1160E+02  0.1094E+02 "alien"
   27    7   12    2.10  925.77 -0.1302E+02  0.6076E+01 "alien"
   28   17   10    2.01  925.25 -0.1301E+02  0.8149E+01 "alien"
   27   12   12    2.03  924.82 -0.1301E+02  0.6584E+01 "alien"
   38    1    3    2.01  924.20 -0.1300E+02  0.5436E+01 "alien"
   38   10    8    1.82  923.36  0.1155E+02  0.5366E+01 "alien"
   40    7    4    1.87  922.29  0.1154E+02  0.6091E+01 "alien"
   32    8    8    2.11  921.82 -0.1297E+02  0.6345E+01 "alien"
   32   17    7    1.99  921.23 -0.1296E+02  0.5406E+01 "alien"
   38   19    4    1.79  919.60  0.1150E+02  0.7920E+01 "alien"
   29   14   13    1.87  916.84  0.1147E+02  0.7621E+01 "alien"
   32   23    8    1.82  915.94  0.1146E+02  0.7340E+01 "alien"
   35   20    7    1.81  912.10  0.1141E+02  0.6846E+01 "alien"
   32   13    8    2.03  911.58 -0.1282E+02  0.6207E+01 "alien"
   37   21    3    1.80  905.69  0.1133E+02  0.6402E+01 "alien"
   33   27    2    1.80  904.44  0.1131E+02  0.1086E+02 "alien"
   26   12   16    1.80  903.74  0.1131E+02  0.6269E+01 "alien"
   35    1    6    2.09  902.72 -0.1270E+02  0.4884E+01 "alien"
   34   20    8    1.81  899.93  0.1126E+02  0.7049E+01 "alien"
   34    1   13    1.80  899.20  0.1125E+02  0.5972E+01 "alien"
   28   24    4    2.05  892.50 -0.1255E+02  0.1122E+02 "alien"
   30   18   11    1.87  891.16  0.1115E+02  0.5029E+01 "alien"
   32    7    8    2.12  888.12 -0.1249E+02  0.6303E+01 "alien"
   33   18    3    2.03  880.97 -0.1239E+02  0.7617E+01 "alien"
   25    8   17    1.80  878.87  0.1099E+02  0.4929E+01 "alien"
   25   15   16    1.78  874.22  0.1094E+02  0.4956E+01 "alien"
   36   21    4    1.82  872.81  0.1092E+02  0.1120E+02 "alien"
   27   19   13    1.84  872.25  0.1091E+02  0.6426E+01 "alien"
   24   14   12    2.10  869.85 -0.1224E+02  0.6435E+01 "alien"
   29   23    5    2.01  868.44 -0.1222E+02  0.1122E+02 "alien"
   31   19    6    2.02  868.38 -0.1221E+02  0.5585E+01 "alien"
   22    9   17    1.86  867.09  0.1085E+02  0.6883E+01 "alien"
   30   23    9    1.84  861.67  0.1078E+02  0.1003E+02 "alien"
   27   25    5    2.02  861.00 -0.1211E+02  0.1566E+02 "alien"
   28   22    5    2.09  860.53 -0.1210E+02  0.1501E+02 "alien"
   35   21    7    1.79  857.59  0.1073E+02  0.7121E+01 "alien"
   20   19   12    2.11  854.53 -0.1202E+02  0.5715E+01 "alien"
   21    7   15    2.07  852.50 -0.1199E+02  0.4645E+01 "alien"
   28   23   10    1.87  844.89  0.1057E+02  0.7180E+01 "alien"
   32   15    5    2.10  844.25 -0.1188E+02  0.1099E+02 "alien"
   33   16   10    1.85  840.43  0.1051E+02  0.7046E+01 "alien"
   31   11    8    2.11  838.82 -0.1180E+02  0.8664E+01 "alien"
   31   10    9    2.07  838.61 -0.1180E+02  0.1028E+02 "alien"
   30   14   13    1.83  837.09  0.1047E+02  0.6172E+01 "alien"
   35   10    6    2.01  832.21 -0.1171E+02  0.6118E+01 "alien"
   35   15    3    2.00  831.82 -0.1170E+02  0.5744E+01 "alien"
   27   23    5    2.10  828.56 -0.1165E+02  0.1581E+02 "alien"
   39   18    1    1.80  827.05  0.1035E+02  0.5645E+01 "alien"
   36   20    6    1.80  826.73  0.1034E+02  0.9130E+01 "alien"
   35   13   11    1.78  826.20  0.1034E+02  0.5006E+01 "alien"
   36    3   11    1.83  825.79  0.1033E+02  0.5158E+01 "alien"
   33   22    7    1.84  824.79  0.1032E+02  0.7412E+01 "alien"
   30   17    7    2.07  820.81 -0.1155E+02  0.5415E+01 "alien"
   37    2   11    1.79  817.25  0.1022E+02  0.5484E+01 "alien"
   28   21   12    1.82  816.48  0.1021E+02  0.4869E+01 "alien"
   42    2    1    1.83  815.61  0.1020E+02  0.5663E+01 "alien"
   38   13    7    1.82  814.19  0.1019E+02  0.7573E+01 "alien"
   21    8   15    2.06  814.03 -0.1145E+02  0.5196E+01 "alien"
   40   14    1    1.82  809.20  0.1012E+02  0.5671E+01 "alien"
   35   22    4    1.83  808.51  0.1011E+02  0.1085E+02 "alien"
   35    1    5    2.12  806.82 -0.1135E+02  0.6275E+01 "alien"
   37   18    3    1.86  806.58  0.1009E+02  0.7507E+01 "alien"
   26   22    9    2.00  806.55 -0.1134E+02  0.1007E+02 "alien"
   41    6    6    1.79  806.49  0.1009E+02  0.6205E+01 "alien"
   37   12    7    1.87  805.82  0.1008E+02  0.7698E+01 "alien"
   24    5   14    2.08  796.73 -0.1121E+02  0.4718E+01 "alien"
   41    1    7    1.78  796.51  0.9964E+01  0.6142E+01 "alien"
   30    5   15    1.81  795.26  0.9948E+01  0.4964E+01 "alien"
   39   14    4    1.83  791.00  0.9895E+01  0.7839E+01 "alien"
   32    4    9    2.09  789.16 -0.1110E+02  0.6948E+01 "alien"
   19   10   15    2.09  781.72 -0.1100E+02  0.4704E+01 "alien"
   30   19    6    2.06  778.21 -0.1095E+02  0.5539E+01 "alien"
   32   13    7    2.07  776.88 -0.1093E+02  0.7768E+01 "alien"
   38   10    9    1.79  774.15  0.9684E+01  0.6887E+01 "alien"
   28   20    6    2.12  774.14 -0.1089E+02  0.1090E+02 "alien"
   34   16   10    1.81  771.43  0.9650E+01  0.7036E+01 "alien"
   24   10   17    1.80  770.09  0.9634E+01  0.4711E+01 "alien"
   41    4    2    1.87  769.00  0.9620E+01  0.6822E+01 "alien"
   27   25    1    2.09  768.84 -0.1081E+02  0.1141E+02 "alien"
   32   19   10    1.83  767.16  0.9597E+01  0.5497E+01 "alien"
   34    5    8    2.04  766.25 -0.1078E+02  0.5344E+01 "alien"
   27   16   14    1.84  762.77  0.9542E+01  0.5285E+01 "alien"
   31   15    7    2.08  762.34 -0.1072E+02  0.7177E+01 "alien"
   34   16    1    2.05  761.15 -0.1071E+02  0.8185E+01 "alien"
   30    8   15    1.79  755.14  0.9447E+01  0.5831E+01 "alien"
   31   11    9    2.06  753.37 -0.1060E+02  0.7564E+01 "alien"
   28   10   12    2.02  751.15 -0.1057E+02  0.9109E+01 "alien"
   27   24   12    1.78  745.95  0.9332E+01  0.6260E+01 "alien"
   26   22    7    2.10  743.99 -0.1046E+02  0.8744E+01 "alien"
   35   25    1    1.79  741.59  0.9277E+01  0.7771E+01 "alien"
   37   16    6    1.84  734.42  0.9187E+01  0.7645E+01 "alien"
   20   17   14    2.01  734.41 -0.1033E+02  0.5332E+01 "alien"
   28   11   11    2.07  734.36 -0.1033E+02  0.9720E+01 "alien"
   29   12   10    2.08  733.68 -0.1032E+02  0.7083E+01 "alien"
   35    7    5    2.08  727.99 -0.1024E+02  0.8321E+01 "alien"
   25    8   16    1.87  727.94  0.9106E+01  0.6052E+01 "alien"
   33    3   10    2.00  724.39 -0.1019E+02  0.5463E+01 "alien"
   41    3    7    1.78  724.30  0.9061E+01  0.6373E+01 "alien"
   24    3   18    1.78  724.23  0.9060E+01  0.6490E+01 "alien"
   24    8   14    2.05  720.77 -0.1014E+02  0.5695E+01 "alien"
   28   16   14    1.81  718.70  0.8991E+01  0.7461E+01 "alien"
   40   10    3    1.85  717.72  0.8978E+01  0.7178E+01 "alien"
   34   16    9    1.85  717.47  0.8975E+01  0.6763E+01 "alien"
   34    9    8    2.00  713.69 -0.1004E+02  0.5840E+01 "alien"
   41    3    6    1.80  711.07  0.8895E+01  0.6509E+01 "alien"
   29   27    6    1.86  707.59  0.8852E+01  0.6697E+01 "alien"
   29   19    6    2.10  707.01 -0.9945E+01  0.9139E+01 "alien"
   41    7    4    1.82  705.91  0.8831E+01  0.6645E+01 "alien"
   30   19    4    2.12  705.20 -0.9919E+01  0.8093E+01 "alien"
   29    5   15    1.84  703.47  0.8800E+01  0.8073E+01 "alien"
   28   15    9    2.12  700.55 -0.9854E+01  0.5111E+01 "alien"
   28   13   11    2.04  699.02 -0.9832E+01  0.5599E+01 "alien"
   21    5   18    1.84  694.11  0.8683E+01  0.5875E+01 "alien"
   26   25    3    2.11  693.97 -0.9761E+01  0.1125E+02 "alien"
   27    4   16    1.84  689.75  0.8629E+01  0.7444E+01 "alien"
   23   22   10    2.05  688.79 -0.9688E+01  0.5977E+01 "alien"
   24   15   13    2.01  688.40 -0.9683E+01  0.5334E+01 "alien"
   39   11    7    1.80  686.77  0.8591E+01  0.6187E+01 "alien"
   31   22    1    2.03  685.53 -0.9643E+01  0.1150E+02 "alien"
   27   15   15    1.79  679.16  0.8496E+01  0.5722E+01 "alien"
   42    7    4    1.78  677.70  0.8478E+01  0.6157E+01 "alien"
   32   29    1    1.79  675.24  0.8447E+01  0.6712E+01 "alien"
   35   23    2    1.84  668.84  0.8367E+01  0.7759E+01 "alien"
   34    3    7    2.09  667.44 -0.9388E+01  0.6172E+01 "alien"
   36    7    5    2.03  664.98 -0.9353E+01  0.8219E+01 "alien"
   28   19   12    1.86  663.81  0.8304E+01  0.6690E+01 "alien"
   34   17    3    2.01  663.80 -0.9337E+01  0.8147E+01 "alien"
   29    3   11    2.11  663.69 -0.9335E+01  0.4884E+01 "alien"
   22   18   13    2.00  661.16 -0.9300E+01  0.5595E+01 "alien"
   41    5    3    1.85  660.19  0.8259E+01  0.7220E+01 "alien"
   28   14   14    1.84  655.43  0.8199E+01  0.5952E+01 "alien"
   42    5    2    1.82  650.61  0.8139E+01  0.5722E+01 "alien"
   31   19    3    2.10  649.98 -0.9142E+01  0.5727E+01 "alien"
   30    4   15    1.81  649.79  0.8129E+01  0.4327E+01 "alien"
   27   17   10    2.05  649.71 -0.9139E+01  0.5032E+01 "alien"
   33    3   13    1.83  648.70  0.8115E+01  0.4776E+01 "alien"
   35    3   11    1.86  646.33  0.8085E+01  0.4192E+01 "alien"
   23   18   12    2.04  645.39 -0.9078E+01  0.7464E+01 "alien"
   32   28    3    1.80  642.84  0.8042E+01  0.1068E+02 "alien"
   35    9    6    2.03  641.80 -0.9028E+01  0.4518E+01 "alien"
   26   17   11    2.03  641.24 -0.9020E+01  0.5520E+01 "alien"
   32   11   12    1.86  639.07  0.7995E+01  0.6331E+01 "alien"
   30   13   10    2.01  638.68 -0.8984E+01  0.7193E+01 "alien"
   14   13   19    1.82  637.96  0.7981E+01  0.2631E+01 "alien"
   34   25    2    1.82  634.61  0.7939E+01  0.7791E+01 "alien"
   24   15   12    2.08  633.76 -0.8914E+01  0.5318E+01 "alien"
   20   12   15    2.03  633.63 -0.8913E+01  0.6415E+01 "alien"
   25    7   14    2.02  632.81 -0.8901E+01  0.4876E+01 "alien"
   27   25    6    1.99  630.89 -0.8874E+01  0.6867E+01 "alien"
   31   15   13    1.79  630.44  0.7887E+01  0.5313E+01 "alien"
   24   20   15    1.78  630.34  0.7885E+01  0.4666E+01 "alien"
   12    4   18    2.02  626.51 -0.8812E+01  0.5349E+01 "alien"
   38   16    1    1.87  624.46  0.7812E+01  0.7979E+01 "alien"
   30   25    8    1.83  623.47  0.7799E+01  0.7305E+01 "alien"
   38    1    9    1.84  622.49  0.7787E+01  0.5849E+01 "alien"
   36   20    2    1.87  621.63  0.7776E+01  0.5870E+01 "alien"
   31    4   14    1.84  619.66  0.7752E+01  0.4775E+01 "alien"
   36    8    5    2.02  615.54 -0.8658E+01  0.6295E+01 "alien"
   32   10    8    2.08  614.85 -0.8648E+01  0.7456E+01 "alien"
   28   12   11    2.06  614.74 -0.8647E+01  0.6478E+01 "alien"
   35   23    5    1.79  614.15  0.7683E+01  0.5588E+01 "alien"
   30    9   15    1.78  612.27  0.7659E+01  0.5880E+01 "alien"
   31   18    6    2.04  611.74 -0.8605E+01  0.7326E+01 "alien"
   33    1    9    2.05  610.62 -0.8589E+01  0.6314E+01 "alien"
   32    8   14    1.79  610.54  0.7638E+01  0.4487E+01 "alien"
   38   12    6    1.86  609.99  0.7631E+01  0.7245E+01 "alien"
   18    1   17    2.01  608.26 -0.8556E+01  0.5652E+01 "alien"
   38   17    1    1.85  603.07  0.7544E+01  0.8037E+01 "alien"
   15    8   17    2.03  601.94 -0.8467E+01  0.3075E+01 "alien"
   37   15    8    1.80  600.70  0.7515E+01  0.5326E+01 "alien"
   40    2    6    1.85  600.08  0.7507E+01  0.4579E+01 "alien"
   21    6   18    1.83  598.93  0.7492E+01  0.4339E+01 "alien"
   21   14   14    2.05  597.69 -0.8407E+01  0.5366E+01 "alien"
   30    9   11    2.02  597.44 -0.8404E+01  0.6805E+01 "alien"
   23   16   16    1.82  596.33  0.7460E+01  0.5877E+01 "alien"
   18    7   16    2.07  596.02 -0.8383E+01  0.4308E+01 "alien"
   26   17   15    1.79  594.84  0.7441E+01  0.5111E+01 "alien"
   19   14   18    1.79  594.33  0.7435E+01  0.6268E+01 "alien"
   20    5   19    1.79  592.59  0.7413E+01  0.5680E+01 "alien"
   38    9    8    1.83  592.26  0.7409E+01  0.5223E+01 "alien"
   29   27    7    1.84  592.15  0.7408E+01  0.6459E+01 "alien"
   35    6   12    1.80  591.18  0.7395E+01  0.5007E+01 "alien"
   30    6   10    2.11  590.63 -0.8308E+01  0.5387E+01 "alien"
   25    6   17    1.81  590.29  0.7384E+01  0.5605E+01 "alien"
   34   11    6    2.05  589.81 -0.8296E+01  0.7956E+01 "alien"
   35   14    2    2.04  589.01 -0.8285E+01  0.5890E+01 "alien"
   35    5   12    1.81  589.00  0.7368E+01  0.4444E+01 "alien"
   41    9    1    1.84  588.76  0.7365E+01  0.4997E+01 "alien"
   42    6    3    1.80  588.57  0.7363E+01  0.5637E+01 "alien"
   39   18    3    1.78  586.66  0.7339E+01  0.5571E+01 "alien"
   28    1   12    2.09  585.56 -0.8236E+01  0.4642E+01 "alien"
   27   22   13    1.78  584.78  0.7315E+01  0.5311E+01 "alien"
   24    6   15    1.99  584.09 -0.8216E+01  0.5217E+01 "alien"
   29    7   12    2.02  582.01 -0.8186E+01  0.9216E+01 "alien"
   35    7    7    2.01  581.87 -0.8184E+01  0.5879E+01 "alien"
   25   12   16    1.82  580.63  0.7263E+01  0.7877E+01 "alien"
   31   19    1    2.12  579.46 -0.8151E+01  0.7192E+01 "alien"
   34    2   13    1.80  577.26  0.7221E+01  0.5185E+01 "alien"
   33   15    2    2.12  575.77 -0.8099E+01  0.6298E+01 "alien"
   32   24    8    1.79  574.40  0.7186E+01  0.7207E+01 "alien"
   33   24    3    1.87  574.15  0.7182E+01  0.7859E+01 "alien"
   33   24    6    1.82  573.07  0.7169E+01  0.7556E+01 "alien"
   32   18    4    2.05  572.04 -0.8046E+01  0.8111E+01 "alien"
   16    6   17    2.03  570.07 -0.8018E+01  0.3831E+01 "alien"
   33   15    7    1.99  567.64 -0.7984E+01  0.7185E+01 "alien"
   37   16    5    1.86  567.01  0.7093E+01  0.5610E+01 "alien"
   41    7    6    1.78  566.84  0.7091E+01  0.6015E+01 "alien"
   33   11   12    1.83  564.99  0.7068E+01  0.6776E+01 "alien"
   34   17    1    2.03  564.07 -0.7934E+01  0.5788E+01 "alien"
   32   20    2    2.03  561.98 -0.7905E+01  0.6213E+01 "alien"
   31   21    4    2.02  560.06 -0.7878E+01  0.8207E+01 "alien"
   39   12    3    1.87  559.30  0.6997E+01  0.5651E+01 "alien"
   33    6   12    1.87  558.75  0.6990E+01  0.4734E+01 "alien"
   41    2    2    1.87  556.79  0.6965E+01  0.5410E+01 "alien"
   30   26    7    1.83  555.00  0.6943E+01  0.1026E+02 "alien"
   41    4    4    1.84  554.52  0.6937E+01  0.5362E+01 "alien"
   33   12    6    2.08  552.63 -0.7773E+01  0.6531E+01 "alien"
   40    6    6    1.83  551.96  0.6905E+01  0.6098E+01 "alien"
   29   24    9    1.85  547.76  0.6852E+01  0.9360E+01 "alien"
   27    9   12    2.07  547.72 -0.7704E+01  0.7533E+01 "alien"
   17    8   19    1.82  545.86  0.6829E+01  0.2949E+01 "alien"
   37    8    2    2.03  545.24 -0.7669E+01  0.5944E+01 "alien"
   29   26    8    1.83  544.63  0.6813E+01  0.6952E+01 "alien"
   31   24    9    1.79  544.28  0.6809E+01  0.9727E+01 "alien"
   37   20    3    1.82  543.43  0.6798E+01  0.6380E+01 "alien"
   33   27    4    1.78  541.51  0.6774E+01  0.7886E+01 "alien"
   37   11    9    1.81  541.09  0.6769E+01  0.9945E+01 "alien"
   33   12    5    2.11  540.90 -0.7608E+01  0.1031E+02 "alien"
   13    3   18    2.01  540.89 -0.7608E+01  0.5269E+01 "alien"
   33   19    2    2.02  540.18 -0.7598E+01  0.5670E+01 "alien"
   40    5    5    1.86  540.12  0.6757E+01  0.6175E+01 "alien"
   34   19    8    1.83  538.59  0.6738E+01  0.7678E+01 "alien"
   25    9   17    1.79  538.46  0.6736E+01  0.5125E+01 "alien"
   23    4   14    2.12  538.29 -0.7572E+01  0.4761E+01 "alien"
   28   17    8    2.12  537.57 -0.7561E+01  0.1037E+02 "alien"
   34    7    7    2.06  537.19 -0.7556E+01  0.4882E+01 "alien"
   14    3   20    1.82  536.94  0.6717E+01  0.3600E+01 "alien"
   33   10   12    1.84  536.80  0.6715E+01  0.4801E+01 "alien"
   36   17    7    1.83  536.43  0.6711E+01  0.5415E+01 "alien"
   36    9    5    2.01  535.90 -0.7538E+01  0.8924E+01 "alien"
   43    1    1    1.79  535.13  0.6694E+01  0.4605E+01 "alien"
   32   28    4    1.78  534.39  0.6685E+01  0.7776E+01 "alien"
   31    8   14    1.82  533.55  0.6675E+01  0.8754E+01 "alien"
   34   17   10    1.80  532.99  0.6668E+01  0.6442E+01 "alien"
   42    5    4    1.79  531.10  0.6644E+01  0.4964E+01 "alien"
   35   20    8    1.78  530.65  0.6638E+01  0.6984E+01 "alien"
   19   14   14    2.11  530.43 -0.7461E+01  0.5517E+01 "alien"
   42    1    4    1.81  530.22  0.6633E+01  0.4794E+01 "alien"
   28   23   11    1.82  529.42  0.6623E+01  0.6231E+01 "alien"
   33   11   13    1.78  526.48  0.6586E+01  0.6558E+01 "alien"
   28   25    4    2.01  525.70 -0.7394E+01  0.1109E+02 "alien"
   39   13    6    1.80  524.96  0.6567E+01  0.7713E+01 "alien"
   26    8   16    1.84  524.44  0.6561E+01  0.5493E+01 "alien"
   13    4   17    2.11  522.43 -0.7348E+01  0.4430E+01 "alien"
   29   18    9    2.00  521.40 -0.7334E+01  0.5084E+01 "alien"
   24   11   14    2.01  521.15 -0.7330E+01  0.4903E+01 "alien"
   27    8   13    2.02  521.05 -0.7329E+01  0.5930E+01 "alien"
   24    4   15    2.00  516.51 -0.7265E+01  0.5176E+01 "alien"
   38    8    7    1.87  515.40  0.6448E+01  0.7557E+01 "alien"
   32    1   11    1.99  512.30 -0.7206E+01  0.5089E+01 "alien"
   28   21    5    2.12  510.93 -0.7187E+01  0.1464E+02 "alien"
   35   17    8    1.84  510.58  0.6387E+01  0.9373E+01 "alien"
   18   11   15    2.11  509.96 -0.7173E+01  0.4912E+01 "alien"
   30   23    1    2.04  508.35 -0.7150E+01  0.8219E+01 "alien"
   39   14    5    1.81  508.17  0.6357E+01  0.8911E+01 "alien"
   32   25    6    1.82  505.66  0.6326E+01  0.7563E+01 "alien"
   27   22    6    2.10  503.98 -0.7089E+01  0.7866E+01 "alien"
   39   17    4    1.78  503.42  0.6298E+01  0.7764E+01 "alien"
   26   22   13    1.80  501.43  0.6273E+01  0.5791E+01 "alien"
   36    9    9    1.87  500.93  0.6267E+01  0.5238E+01 "alien"
   31   15    8    2.03  499.75 -0.7029E+01  0.7385E+01 "alien"
   39    9    8    1.79  499.15  0.6244E+01  0.4204E+01 "alien"
   24   19   11    2.05  497.86 -0.7003E+01  0.9604E+01 "alien"
   30    4   11    2.06  496.79 -0.6988E+01  0.5470E+01 "alien"
   27   22    8    2.01  495.36 -0.6968E+01  0.7419E+01 "alien"
   27   12   11    2.10  491.47 -0.6913E+01  0.7874E+01 "alien"
   35    3    7    2.04  491.47 -0.6913E+01  0.4922E+01 "alien"
   26    9   12    2.12  491.37 -0.6912E+01  0.6243E+01 "alien"
   39   24   10    1.55  491.03 -0.2761E+02  0.1095E+02 "alien"
   41    5    1    1.87  490.30  0.6134E+01  0.5419E+01 "alien"
   29    8   11    2.07  490.14 -0.6894E+01  0.9796E+01 "alien"
   27   18    9    2.08  489.19 -0.6881E+01  0.7485E+01 "alien"
   28   24   10    1.84  488.79  0.6115E+01  0.4984E+01 "alien"
   18    2   16    2.10  487.96 -0.6864E+01  0.5669E+01 "alien"
   40   16    2    1.78  487.44  0.6098E+01  0.5545E+01 "alien"
   29   15    9    2.07  486.62 -0.6845E+01  0.7330E+01 "alien"
   35   14    1    2.05  485.57 -0.6830E+01  0.5827E+01 "alien"
   14    2   17    2.10  484.82 -0.6819E+01  0.4580E+01 "alien"
   37    7    2    2.04  484.65 -0.6817E+01  0.5399E+01 "alien"
   30   10   14    1.83  484.30  0.6058E+01  0.6193E+01 "alien"
   27   26    9    1.86  482.86  0.6040E+01  0.5408E+01 "alien"
   32   14   12    1.82  482.15  0.6032E+01  0.6375E+01 "alien"
   29   11   10    2.09  480.46 -0.6758E+01  0.5899E+01 "alien"
   31   30    1    1.79  479.32  0.5996E+01  0.6526E+01 "alien"
   19   11   15    2.08  479.17 -0.6740E+01  0.5101E+01 "alien"
   18    9   16    2.04  478.33 -0.6728E+01  0.5177E+01 "alien"
   29   24   10    1.81  477.77  0.5977E+01  0.6915E+01 "alien"
   31   16    7    2.05  476.05 -0.6696E+01  0.7645E+01 "alien"
   26    3   14    2.02  476.05 -0.6696E+01  0.4163E+01 "alien"
   41   13    2    1.79  475.79  0.5952E+01  0.5609E+01 "alien"
   33    2    8    2.10  475.69 -0.6691E+01  0.6913E+01 "alien"
   34   16    2    2.04  475.26 -0.6685E+01  0.8226E+01 "alien"
   34    5    7    2.08  475.07 -0.6682E+01  0.4436E+01 "alien"
   36    2    3    2.11  474.63 -0.6676E+01  0.4854E+01 "alien"
   41    8    4    1.82  472.89  0.5916E+01  0.6166E+01 "alien"
   30   18    6    2.09  472.26 -0.6643E+01  0.5928E+01 "alien"
   18   17   14    2.07  471.83 -0.6637E+01  0.5535E+01 "alien"
   28   16    9    2.09  471.12 -0.6627E+01  0.9986E+01 "alien"
   12    8   17    2.09  470.89 -0.6624E+01  0.4769E+01 "alien"
   40   15    2    1.80  470.58  0.5887E+01  0.5542E+01 "alien"
   16    3   19    1.87  470.41  0.5885E+01  0.6194E+01 "alien"
   37    1    2    2.07  468.84 -0.6595E+01  0.3954E+01 "alien"
   41   12    3    1.79  468.26  0.5858E+01  0.7799E+01 "alien"
   40    7    7    1.80  467.90  0.5853E+01  0.5841E+01 "alien"
   39    2    7    1.86  466.29  0.5833E+01  0.4419E+01 "alien"
   37    1    4    2.04  465.88 -0.6553E+01  0.5486E+01 "alien"
   22    6   18    1.81  465.41  0.5822E+01  0.4483E+01 "alien"
   39   16    5    1.78  464.46  0.5810E+01  0.5455E+01 "alien"
   23   18   16    1.78  463.50  0.5798E+01  0.4420E+01 "alien"
   32   21   10    1.79  463.14  0.5794E+01  0.9880E+01 "alien"
   28   25   10    1.82  461.84  0.5777E+01  0.6079E+01 "alien"
   42    7    1    1.81  457.52  0.5723E+01  0.4460E+01 "alien"
   25   14   13    2.00  457.12 -0.6430E+01  0.4491E+01 "alien"
   27   23    7    2.02  456.49 -0.6421E+01  0.1092E+02 "alien"
   35   16    9    1.82  456.39  0.5709E+01  0.1001E+02 "alien"
   20   12   18    1.80  456.27  0.5708E+01  0.4892E+01 "alien"
   26   15   15    1.82  454.34  0.5684E+01  0.5721E+01 "alien"
   37    8   10    1.80  454.32  0.5683E+01  0.5655E+01 "alien"
   31   26    4    1.87  454.20  0.5682E+01  0.1023E+02 "alien"
   28   25    5    1.99  453.56 -0.6380E+01  0.1429E+02 "alien"
   25   15   12    2.05  453.48 -0.6379E+01  0.4729E+01 "alien"
   36   11   10    1.81  453.04  0.5667E+01  0.6656E+01 "alien"
   36    1    6    2.03  449.38 -0.6321E+01  0.3822E+01 "alien"
   30   12   10    2.03  447.45 -0.6294E+01  0.1023E+02 "alien"
   28   12   15    1.81  445.72  0.5576E+01  0.4823E+01 "alien"
   42    2    4    1.80  445.24  0.5570E+01  0.5027E+01 "alien"
   20   18   17    1.78  444.12  0.5556E+01  0.4203E+01 "alien"
   24   23   13    1.83  443.05  0.5542E+01  0.5624E+01 "alien"
   43    3    2    1.78  442.76  0.5539E+01  0.4416E+01 "alien"
   36   12    3    2.01  440.86 -0.6201E+01  0.7869E+01 "alien"
   34   12   11    1.83  440.69  0.5513E+01  0.6511E+01 "alien"
   36    9   11    1.79  439.89  0.5503E+01  0.6828E+01 "alien"
   32   22    8    1.84  439.31  0.5496E+01  0.7225E+01 "alien"
   22   15   13    2.08  439.17 -0.6177E+01  0.5072E+01 "alien"
   33   13    7    2.02  439.09 -0.6176E+01  0.7716E+01 "alien"
   34   11    5    2.08  438.57 -0.6169E+01  0.7970E+01 "alien"
   15    2   17    2.08  437.75 -0.6157E+01  0.4615E+01 "alien"
   26   14   15    1.83  436.95  0.5466E+01  0.6820E+01 "alien"
   25   10   17    1.78  436.40  0.5459E+01  0.5201E+01 "alien"
   26   12   13    2.00  435.71 -0.6129E+01  0.6289E+01 "alien"
   29   24    1    2.05  435.23 -0.6122E+01  0.8075E+01 "alien"
   40   13    1    1.83  434.57  0.5436E+01  0.7330E+01 "alien"
   38   15    3    1.87  434.26  0.5432E+01  0.7917E+01 "alien"
   25   19   15    1.78  433.93  0.5428E+01  0.4344E+01 "alien"
   31   19    5    2.05  433.83 -0.6102E+01  0.1088E+02 "alien"
   36   11    9    1.85  433.53  0.5423E+01  0.5193E+01 "alien"
   30   20    6    2.03  433.21 -0.6094E+01  0.9081E+01 "alien"
   26   24   12    1.81  433.14  0.5418E+01  0.6282E+01 "alien"
   15    4   20    1.80  433.07  0.5418E+01  0.3948E+01 "alien"
   32    9    8    2.09  431.88 -0.6075E+01  0.4934E+01 "alien"
   35    7    4    2.11  431.79 -0.6074E+01  0.5167E+01 "alien"
   42    8    2    1.80  431.35  0.5396E+01  0.5036E+01 "alien"
   21   13   15    1.99  431.17 -0.6065E+01  0.4936E+01 "alien"
   23   15   16    1.83  430.65  0.5387E+01  0.5442E+01 "alien"
   33   10    7    2.07  428.79 -0.6031E+01  0.6252E+01 "alien"
   39   16    2    1.82  427.99  0.5354E+01  0.5564E+01 "alien"
   14    3   18    1.99  425.65 -0.5987E+01  0.5761E+01 "alien"
   29   22    5    2.05  425.47 -0.5985E+01  0.1127E+02 "alien"
   29   15   13    1.85  424.49  0.5310E+01  0.8781E+01 "alien"
   40   10    1    1.87  424.14  0.5306E+01  0.5482E+01 "alien"
   32   19   11    1.79  421.61  0.5274E+01  0.6508E+01 "alien"
   36    4   10    1.87  420.10  0.5255E+01  0.4647E+01 "alien"
   29   19   13    1.78  419.66  0.5250E+01  0.6191E+01 "alien"
   37   14    8    1.81  418.39  0.5234E+01  0.5255E+01 "alien"
   28   20    8    2.03  418.10 -0.5881E+01  0.1002E+02 "alien"
   30    6   11    2.05  415.46 -0.5844E+01  0.5321E+01 "alien"
   24   10   13    2.10  414.64 -0.5832E+01  0.6153E+01 "alien"
   34   12    7    2.00  414.05 -0.5824E+01  0.7491E+01 "alien"
   25   24    8    2.02  413.58 -0.5817E+01  0.6287E+01 "alien"
   23   12   17    1.81  413.18  0.5169E+01  0.4719E+01 "alien"
   28    9   11    2.10  412.52 -0.5802E+01  0.6868E+01 "alien"
   26   16   14    1.86  412.38  0.5159E+01  0.5168E+01 "alien"
   33   19    1    2.02  412.27 -0.5799E+01  0.5741E+01 "alien"
   31   17    5    2.10  412.23 -0.5798E+01  0.7995E+01 "alien"
   38    3    3    2.00  411.99 -0.5795E+01  0.6245E+01 "alien"
   25    7   13    2.10  411.97 -0.5795E+01  0.6129E+01 "alien"
   26   13   11    2.12  411.76 -0.5792E+01  0.5537E+01 "alien"
   36    3    3    2.11  410.44 -0.5773E+01  0.5816E+01 "alien"
   29   11   11    2.03  410.32 -0.5771E+01  0.6528E+01 "alien"
   25    7   17    1.81  409.75  0.5126E+01  0.4686E+01 "alien"
   48   25   10    1.34  409.63  0.9290E+03  0.4785E+02 "alien"
   36   13   10    1.79  408.94  0.5116E+01  0.7133E+01 "alien"
   23   16   13    2.02  407.98 -0.5739E+01  0.6070E+01 "alien"
   27   20    8    2.07  407.90 -0.5737E+01  0.1039E+02 "alien"
   36   24    2    1.78  407.11  0.5093E+01  0.5530E+01 "alien"
   25    6   14    2.03  406.46 -0.5717E+01  0.5515E+01 "alien"
   39    9    7    1.82  406.36  0.5083E+01  0.5966E+01 "alien"
   30   17    6    2.11  405.16 -0.5699E+01  0.7827E+01 "alien"
   34    1   12    1.86  401.85  0.5027E+01  0.4390E+01 "alien"
   27    8   16    1.81  400.60  0.5011E+01  0.7887E+01 "alien"
   32   16    6    2.05  400.43 -0.5632E+01  0.5617E+01 "alien"
   31    2   11    2.03  400.35 -0.5631E+01  0.5313E+01 "alien"
   36    8    2    2.08  399.38 -0.5618E+01  0.4696E+01 "alien"
   42    3    1    1.83  399.34  0.4996E+01  0.5832E+01 "alien"
   23   11   17    1.82  398.71  0.4988E+01  0.4777E+01 "alien"
   35    8    3    2.12  398.56 -0.5606E+01  0.6081E+01 "alien"
   23    7   18    1.78  398.06  0.4980E+01  0.6578E+01 "alien"
   27   17   14    1.82  396.07  0.4955E+01  0.6049E+01 "alien"
   21   10   18    1.80  396.04  0.4954E+01  0.4301E+01 "alien"
   17   12   15    2.12  395.84 -0.5568E+01  0.4353E+01 "alien"
   40    8    7    1.79  393.68  0.4925E+01  0.5308E+01 "alien"
   53   10   11    1.33  393.14  0.8916E+03  0.4533E+02 "alien"
   41   11    2    1.81  392.26  0.4907E+01  0.6335E+01 "alien"
   26   25    7    2.00  391.30 -0.5504E+01  0.5670E+01 "alien"
   29   24    4    2.01  391.25 -0.5503E+01  0.1122E+02 "alien"
   28    6   15    1.86  390.61  0.4886E+01  0.5415E+01 "alien"
   36   18    7    1.81  390.10  0.4880E+01  0.5348E+01 "alien"
   29   17    9    2.03  389.84 -0.5483E+01  0.7075E+01 "alien"
   37   11    2    1.99  387.02 -0.5444E+01  0.1169E+02 "alien"
   34   10    7    2.03  386.06 -0.5430E+01  0.5814E+01 "alien"
   24    4   17    1.85  386.00  0.4829E+01  0.6558E+01 "alien"
   24   21   11    1.99  385.85 -0.5427E+01  0.5797E+01 "alien"
   31    3   14    1.84  385.76  0.4826E+01  0.4769E+01 "alien"
   35    6    5    2.09  385.57 -0.5423E+01  0.7344E+01 "alien"
   21   11   17    1.86  384.47  0.4810E+01  0.6555E+01 "alien"
   38   16    6    1.80  384.11  0.4805E+01  0.1062E+02 "alien"
   38    1    8    1.87  383.81  0.4801E+01  0.4815E+01 "alien"
   32   10   13    1.82  383.21  0.4794E+01  0.6463E+01 "alien"
   34   23    4    1.85  380.64  0.4762E+01  0.1107E+02 "alien"
   30   27    7    1.81  380.19  0.4756E+01  0.7417E+01 "alien"
   29   21   11    1.84  378.96  0.4741E+01  0.5588E+01 "alien"
   30   13   14    1.79  377.83  0.4727E+01  0.7460E+01 "alien"


 List of    500 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53   10   11    1.33  393.14  0.8916E+03  0.4533E+02
   48   25   10    1.34  409.63  0.9290E+03  0.4785E+02
   39   24   10    1.55  491.03 -0.2761E+02  0.1095E+02
   36   24    2    1.78  407.11  0.5093E+01  0.5530E+01
   27   22   13    1.78  584.78  0.7315E+01  0.5311E+01
   28    9   16    1.78 1146.66  0.1434E+02  0.7656E+01
   41    3    7    1.78  724.30  0.9061E+01  0.6373E+01
   33   11   13    1.78  526.48  0.6586E+01  0.6558E+01
   25   19   15    1.78  433.93  0.5428E+01  0.4344E+01
   33   27    4    1.78  541.51  0.6774E+01  0.7886E+01
   39   18    3    1.78  586.66  0.7339E+01  0.5571E+01
   35   20    8    1.78  530.65  0.6638E+01  0.6984E+01
   27   24   12    1.78  745.95  0.9332E+01  0.6260E+01
   30    9   15    1.78  612.27  0.7659E+01  0.5880E+01
   25   10   17    1.78  436.40  0.5459E+01  0.5201E+01
   24    3   18    1.78  724.23  0.9060E+01  0.6490E+01
   20   18   17    1.78  444.12  0.5556E+01  0.4203E+01
   41    1    7    1.78  796.51  0.9964E+01  0.6142E+01
   35   13   11    1.78  826.20  0.1034E+02  0.5006E+01
   42    7    4    1.78  677.70  0.8478E+01  0.6157E+01
   39   16    5    1.78  464.46  0.5810E+01  0.5455E+01
   29   19   13    1.78  419.66  0.5250E+01  0.6191E+01
   32   23    9    1.78 1228.23  0.1536E+02  0.6096E+01
   23   18   16    1.78  463.50  0.5798E+01  0.4420E+01
   24   20   15    1.78  630.34  0.7885E+01  0.4666E+01
   43    3    2    1.78  442.76  0.5539E+01  0.4416E+01
   39   17    4    1.78  503.42  0.6298E+01  0.7764E+01
   41    7    6    1.78  566.84  0.7091E+01  0.6015E+01
   25   15   16    1.78  874.22  0.1094E+02  0.4956E+01
   40   16    2    1.78  487.44  0.6098E+01  0.5545E+01
   32   28    4    1.78  534.39  0.6685E+01  0.7776E+01
   23    7   18    1.78  398.06  0.4980E+01  0.6578E+01
   32    8   14    1.79  610.54  0.7638E+01  0.4487E+01
   32   29    1    1.79  675.24  0.8447E+01  0.6712E+01
   38   19    4    1.79  919.60  0.1150E+02  0.7920E+01
   32   19   11    1.79  421.61  0.5274E+01  0.6508E+01
   38   10    9    1.79  774.15  0.9684E+01  0.6887E+01
   20    5   19    1.79  592.59  0.7413E+01  0.5680E+01
   41   13    2    1.79  475.79  0.5952E+01  0.5609E+01
   31   30    1    1.79  479.32  0.5996E+01  0.6526E+01
   32   21   10    1.79  463.14  0.5794E+01  0.9880E+01
   36   13   10    1.79  408.94  0.5116E+01  0.7133E+01
   26   17   15    1.79  594.84  0.7441E+01  0.5111E+01
   30    8   15    1.79  755.14  0.9447E+01  0.5831E+01
   35   21    7    1.79  857.59  0.1073E+02  0.7121E+01
   41   12    3    1.79  468.26  0.5858E+01  0.7799E+01
   41    6    6    1.79  806.49  0.1009E+02  0.6205E+01
   31   15   13    1.79  630.44  0.7887E+01  0.5313E+01
   31   24    9    1.79  544.28  0.6809E+01  0.9727E+01
   36    9   11    1.79  439.89  0.5503E+01  0.6828E+01
   33   14   12    1.79  998.76  0.1249E+02  0.6111E+01
   25    9   17    1.79  538.46  0.6736E+01  0.5125E+01
   40    8    7    1.79  393.68  0.4925E+01  0.5308E+01
   39    9    8    1.79  499.15  0.6244E+01  0.4204E+01
   37    2   11    1.79  817.25  0.1022E+02  0.5484E+01
   32   24    8    1.79  574.40  0.7186E+01  0.7207E+01
   35   25    1    1.79  741.59  0.9277E+01  0.7771E+01
   43    1    1    1.79  535.13  0.6694E+01  0.4605E+01
   35   23    5    1.79  614.15  0.7683E+01  0.5588E+01
   42    5    4    1.79  531.10  0.6644E+01  0.4964E+01
   36   19    7    1.79  993.55  0.1243E+02  0.6215E+01
   27   15   15    1.79  679.16  0.8496E+01  0.5722E+01
   30   13   14    1.79  377.83  0.4727E+01  0.7460E+01
   19   14   18    1.79  594.33  0.7435E+01  0.6268E+01
   37   15    8    1.80  600.70  0.7515E+01  0.5326E+01
   39   18    1    1.80  827.05  0.1035E+02  0.5645E+01
   37   21    3    1.80  905.69  0.1133E+02  0.6402E+01
   34   25    4    1.80 1228.18  0.1536E+02  0.1067E+02
   40    7    7    1.80  467.90  0.5853E+01  0.5841E+01
   34   17   10    1.80  532.99  0.6668E+01  0.6442E+01
   42    8    2    1.80  431.35  0.5396E+01  0.5036E+01
   32   28    3    1.80  642.84  0.8042E+01  0.1068E+02
   38   16    6    1.80  384.11  0.4805E+01  0.1062E+02
   25    8   17    1.80  878.87  0.1099E+02  0.4929E+01
   26   12   16    1.80  903.74  0.1131E+02  0.6269E+01
   40   15    2    1.80  470.58  0.5887E+01  0.5542E+01
   39   16    4    1.80  998.73  0.1249E+02  0.7716E+01
   36   20    6    1.80  826.73  0.1034E+02  0.9130E+01
   15    4   20    1.80  433.07  0.5418E+01  0.3948E+01
   39   11    7    1.80  686.77  0.8591E+01  0.6187E+01
   20   12   18    1.80  456.27  0.5708E+01  0.4892E+01
   35   24    3    1.80 1908.61  0.2388E+02  0.1021E+02
   42    6    3    1.80  588.57  0.7363E+01  0.5637E+01
   35    6   12    1.80  591.18  0.7395E+01  0.5007E+01
   26   22   13    1.80  501.43  0.6273E+01  0.5791E+01
   34    2   13    1.80  577.26  0.7221E+01  0.5185E+01
   21   10   18    1.80  396.04  0.4954E+01  0.4301E+01
   33   27    2    1.80  904.44  0.1131E+02  0.1086E+02
   41    3    6    1.80  711.07  0.8895E+01  0.6509E+01
   39   13    6    1.80  524.96  0.6567E+01  0.7713E+01
   37    8   10    1.80  454.32  0.5683E+01  0.5655E+01
   42    2    4    1.80  445.24  0.5570E+01  0.5027E+01
   34    1   13    1.80  899.20  0.1125E+02  0.5972E+01
   24   10   17    1.80  770.09  0.9634E+01  0.4711E+01
   42    1    4    1.81  530.22  0.6633E+01  0.4794E+01
   30   23   10    1.81 1137.08  0.1422E+02  0.9079E+01
   25    7   17    1.81  409.75  0.5126E+01  0.4686E+01
   26   24   12    1.81  433.14  0.5418E+01  0.6282E+01
   23   12   17    1.81  413.18  0.5169E+01  0.4719E+01
   28   12   15    1.81  445.72  0.5576E+01  0.4823E+01
   35    5   12    1.81  589.00  0.7368E+01  0.4444E+01
   30   27    7    1.81  380.19  0.4756E+01  0.7417E+01
   30    5   15    1.81  795.26  0.9948E+01  0.4964E+01
   28   16   14    1.81  718.70  0.8991E+01  0.7461E+01
   35   20    7    1.81  912.10  0.1141E+02  0.6846E+01
   37   14    8    1.81  418.39  0.5234E+01  0.5255E+01
   38   11    8    1.81 1043.67  0.1306E+02  0.7019E+01
   41   11    2    1.81  392.26  0.4907E+01  0.6335E+01
   42    7    1    1.81  457.52  0.5723E+01  0.4460E+01
   39   14    5    1.81  508.17  0.6357E+01  0.8911E+01
   36   11   10    1.81  453.04  0.5667E+01  0.6656E+01
   29   24   10    1.81  477.77  0.5977E+01  0.6915E+01
   36   18    7    1.81  390.10  0.4880E+01  0.5348E+01
   27    8   16    1.81  400.60  0.5011E+01  0.7887E+01
   25    6   17    1.81  590.29  0.7384E+01  0.5605E+01
   30    4   15    1.81  649.79  0.8129E+01  0.4327E+01
   22    6   18    1.81  465.41  0.5822E+01  0.4483E+01
   37   11    9    1.81  541.09  0.6769E+01  0.9945E+01
   34   20    8    1.81  899.93  0.1126E+02  0.7049E+01
   34   16   10    1.81  771.43  0.9650E+01  0.7036E+01
   41    8    4    1.82  472.89  0.5916E+01  0.6166E+01
   33   24    6    1.82  573.07  0.7169E+01  0.7556E+01
   38   13    7    1.82  814.19  0.1019E+02  0.7573E+01
   32   23    8    1.82  915.94  0.1146E+02  0.7340E+01
   28   25   10    1.82  461.84  0.5777E+01  0.6079E+01
   40   12    4    1.82 1115.56  0.1396E+02  0.7924E+01
   31    8   14    1.82  533.55  0.6675E+01  0.8754E+01
   23   16   16    1.82  596.33  0.7460E+01  0.5877E+01
   14   13   19    1.82  637.96  0.7981E+01  0.2631E+01
   42    5    2    1.82  650.61  0.8139E+01  0.5722E+01
   37   20    3    1.82  543.43  0.6798E+01  0.6380E+01
   35   16    9    1.82  456.39  0.5709E+01  0.1001E+02
   23   11   17    1.82  398.71  0.4988E+01  0.4777E+01
   14    3   20    1.82  536.94  0.6717E+01  0.3600E+01
   31   21   10    1.82 1347.13  0.1685E+02  0.7110E+01
   36   21    4    1.82  872.81  0.1092E+02  0.1120E+02
   27   17   14    1.82  396.07  0.4955E+01  0.6049E+01
   40   14    1    1.82  809.20  0.1012E+02  0.5671E+01
   28   21   12    1.82  816.48  0.1021E+02  0.4869E+01
   38   10    8    1.82  923.36  0.1155E+02  0.5366E+01
   26   15   15    1.82  454.34  0.5684E+01  0.5721E+01
   32   10   13    1.82  383.21  0.4794E+01  0.6463E+01
   39    9    7    1.82  406.36  0.5083E+01  0.5966E+01
   34   25    2    1.82  634.61  0.7939E+01  0.7791E+01
   38   17    4    1.82 1260.29  0.1577E+02  0.1101E+02
   28   23   11    1.82  529.42  0.6623E+01  0.6231E+01
   23   22   14    1.82  932.09  0.1166E+02  0.5804E+01
   32   14   12    1.82  482.15  0.6032E+01  0.6375E+01
   17    8   19    1.82  545.86  0.6829E+01  0.2949E+01
   41    7    4    1.82  705.91  0.8831E+01  0.6645E+01
   39   16    2    1.82  427.99  0.5354E+01  0.5564E+01
   32   25    6    1.82  505.66  0.6326E+01  0.7563E+01
   25   12   16    1.82  580.63  0.7263E+01  0.7877E+01
   35    9   11    1.83 1110.12  0.1389E+02  0.6428E+01
   32   21    9    1.83  948.89  0.1187E+02  0.8283E+01
   36    3   11    1.83  825.79  0.1033E+02  0.5158E+01
   33   11   12    1.83  564.99  0.7068E+01  0.6776E+01
   32   19   10    1.83  767.16  0.9597E+01  0.5497E+01
   37    4   10    1.83 1366.96  0.1710E+02  0.5623E+01
   34   12   11    1.83  440.69  0.5513E+01  0.6511E+01
   30   10   14    1.83  484.30  0.6058E+01  0.6193E+01
   31   25    7    1.83 1269.55  0.1588E+02  0.1072E+02
   18    4   19    1.83 1077.50  0.1348E+02  0.5032E+01
   39   14    4    1.83  791.00  0.9895E+01  0.7839E+01
   36   17    7    1.83  536.43  0.6711E+01  0.5415E+01
   30   25    8    1.83  623.47  0.7799E+01  0.7305E+01
   38    9    8    1.83  592.26  0.7409E+01  0.5223E+01
   40    6    6    1.83  551.96  0.6905E+01  0.6098E+01
   24   23   13    1.83  443.05  0.5542E+01  0.5624E+01
   42    3    1    1.83  399.34  0.4996E+01  0.5832E+01
   23   15   16    1.83  430.65  0.5387E+01  0.5442E+01
   29    6   15    1.83 1334.22  0.1669E+02  0.8253E+01
   40   13    1    1.83  434.57  0.5436E+01  0.7330E+01
   35   22    4    1.83  808.51  0.1011E+02  0.1085E+02
   33    3   13    1.83  648.70  0.8115E+01  0.4776E+01
   29   26    8    1.83  544.63  0.6813E+01  0.6952E+01
   34   19    8    1.83  538.59  0.6738E+01  0.7678E+01
   42    2    1    1.83  815.61  0.1020E+02  0.5663E+01
   30   26    7    1.83  555.00  0.6943E+01  0.1026E+02
   30   14   13    1.83  837.09  0.1047E+02  0.6172E+01
   21    6   18    1.83  598.93  0.7492E+01  0.4339E+01
   26   14   15    1.83  436.95  0.5466E+01  0.6820E+01
   35   23    2    1.84  668.84  0.8367E+01  0.7759E+01
   35   17    8    1.84  510.58  0.6387E+01  0.9373E+01
   37   16    6    1.84  734.42  0.9187E+01  0.7645E+01
   33   22    7    1.84  824.79  0.1032E+02  0.7412E+01
   28   27    8    1.84 1237.64  0.1548E+02  0.7022E+01
   38    1    9    1.84  622.49  0.7787E+01  0.5849E+01
   27   16   14    1.84  762.77  0.9542E+01  0.5285E+01
   27    4   16    1.84  689.75  0.8629E+01  0.7444E+01
   41    9    1    1.84  588.76  0.7365E+01  0.4997E+01
   29   27    7    1.84  592.15  0.7408E+01  0.6459E+01
   29   21   11    1.84  378.96  0.4741E+01  0.5588E+01
   27   19   13    1.84  872.25  0.1091E+02  0.6426E+01
   29    5   15    1.84  703.47  0.8800E+01  0.8073E+01
   33   10   12    1.84  536.80  0.6715E+01  0.4801E+01
   32   22    8    1.84  439.31  0.5496E+01  0.7225E+01
   28   14   14    1.84  655.43  0.8199E+01  0.5952E+01
   26    8   16    1.84  524.44  0.6561E+01  0.5493E+01
   41    4    4    1.84  554.52  0.6937E+01  0.5362E+01
   37   14    7    1.84 1131.20  0.1415E+02  0.7689E+01
   31   20   10    1.84  970.99  0.1215E+02  0.6829E+01
   36   20    4    1.84 2210.61  0.2765E+02  0.1128E+02
   21    5   18    1.84  694.11  0.8683E+01  0.5875E+01
   31    4   14    1.84  619.66  0.7752E+01  0.4775E+01
   28   24   10    1.84  488.79  0.6115E+01  0.4984E+01
   31    3   14    1.84  385.76  0.4826E+01  0.4769E+01
   30   23    9    1.84  861.67  0.1078E+02  0.1003E+02
   28   22   11    1.84 1559.52  0.1951E+02  0.8034E+01
   34   23    4    1.85  380.64  0.4762E+01  0.1107E+02
   32   25    5    1.85 1187.68  0.1486E+02  0.1062E+02
   38    7    8    1.85 1517.00  0.1898E+02  0.7528E+01
   40   12    1    1.85 1185.09  0.1483E+02  0.1090E+02
   40    2    6    1.85  600.08  0.7507E+01  0.4579E+01
   24    4   17    1.85  386.00  0.4829E+01  0.6558E+01
   33   20    8    1.85 1525.54  0.1908E+02  0.9715E+01
   33   16   10    1.85  840.43  0.1051E+02  0.7046E+01
   41    5    3    1.85  660.19  0.8259E+01  0.7220E+01
   29   15   13    1.85  424.49  0.5310E+01  0.8781E+01
   29   24    9    1.85  547.76  0.6852E+01  0.9360E+01
   36   11    9    1.85  433.53  0.5423E+01  0.5193E+01
   34    3   12    1.85 1307.70  0.1636E+02  0.4912E+01
   38   17    1    1.85  603.07  0.7544E+01  0.8037E+01
   40   10    3    1.85  717.72  0.8978E+01  0.7178E+01
   34   16    9    1.85  717.47  0.8975E+01  0.6763E+01
   31   17   11    1.85 1024.25  0.1281E+02  0.6442E+01
   34    1   12    1.86  401.85  0.5027E+01  0.4390E+01
   38   12    6    1.86  609.99  0.7631E+01  0.7245E+01
   30   16   12    1.86 1130.22  0.1414E+02  0.6468E+01
   37   18    3    1.86  806.58  0.1009E+02  0.7507E+01
   27   26    9    1.86  482.86  0.6040E+01  0.5408E+01
   37   16    5    1.86  567.01  0.7093E+01  0.5610E+01
   40    5    5    1.86  540.12  0.6757E+01  0.6175E+01
   39   14    1    1.86  927.64  0.1160E+02  0.1094E+02
   39    2    7    1.86  466.29  0.5833E+01  0.4419E+01
   28   19   12    1.86  663.81  0.8304E+01  0.6690E+01
   32   11   12    1.86  639.07  0.7995E+01  0.6331E+01
   22    9   17    1.86  867.09  0.1085E+02  0.6883E+01
   28    6   15    1.86  390.61  0.4886E+01  0.5415E+01
   35    3   11    1.86  646.33  0.8085E+01  0.4192E+01
   31    9   13    1.86 1634.08  0.2044E+02  0.9028E+01
   36   15    7    1.86 1462.57  0.1830E+02  0.1057E+02
   21   11   17    1.86  384.47  0.4810E+01  0.6555E+01
   26   16   14    1.86  412.38  0.5159E+01  0.5168E+01
   29   27    6    1.86  707.59  0.8852E+01  0.6697E+01
   41    4    2    1.87  769.00  0.9620E+01  0.6822E+01
   40    7    4    1.87  922.29  0.1154E+02  0.6091E+01
   38   13    5    1.87 1265.51  0.1583E+02  0.8963E+01
   28   23   10    1.87  844.89  0.1057E+02  0.7180E+01
   36   20    2    1.87  621.63  0.7776E+01  0.5870E+01
   29   14   13    1.87  916.84  0.1147E+02  0.7621E+01
   36    4   10    1.87  420.10  0.5255E+01  0.4647E+01
   25    8   16    1.87  727.94  0.9106E+01  0.6052E+01
   41    5    1    1.87  490.30  0.6134E+01  0.5419E+01
   37   12    7    1.87  805.82  0.1008E+02  0.7698E+01
   30   18   11    1.87  891.16  0.1115E+02  0.5029E+01
   16    3   19    1.87  470.41  0.5885E+01  0.6194E+01
   38   15    3    1.87  434.26  0.5432E+01  0.7917E+01
   38   16    1    1.87  624.46  0.7812E+01  0.7979E+01
   40   10    1    1.87  424.14  0.5306E+01  0.5482E+01
   38    8    7    1.87  515.40  0.6448E+01  0.7557E+01
   41    2    2    1.87  556.79  0.6965E+01  0.5410E+01
   33   24    3    1.87  574.15  0.7182E+01  0.7859E+01
   39   12    3    1.87  559.30  0.6997E+01  0.5651E+01
   31   26    4    1.87  454.20  0.5682E+01  0.1023E+02
   38    1    8    1.87  383.81  0.4801E+01  0.4815E+01
   36    9    9    1.87  500.93  0.6267E+01  0.5238E+01
   33    6   12    1.87  558.75  0.6990E+01  0.4734E+01
   32   24    5    1.87 2426.75  0.3036E+02  0.1582E+02
   24    6   15    1.99  584.09 -0.8216E+01  0.5217E+01
   33   15    7    1.99  567.64 -0.7984E+01  0.7185E+01
   28   25    5    1.99  453.56 -0.6380E+01  0.1429E+02
   32   17    7    1.99  921.23 -0.1296E+02  0.5406E+01
   21   13   15    1.99  431.17 -0.6065E+01  0.4936E+01
   32    1   11    1.99  512.30 -0.7206E+01  0.5089E+01
   37   11    2    1.99  387.02 -0.5444E+01  0.1169E+02
   24   21   11    1.99  385.85 -0.5427E+01  0.5797E+01
   14    3   18    1.99  425.65 -0.5987E+01  0.5761E+01
   27   25    6    1.99  630.89 -0.8874E+01  0.6867E+01
   33    3   10    2.00  724.39 -0.1019E+02  0.5463E+01
   26   25    7    2.00  391.30 -0.5504E+01  0.5670E+01
   25   14   13    2.00  457.12 -0.6430E+01  0.4491E+01
   34   12    7    2.00  414.05 -0.5824E+01  0.7491E+01
   30   20    7    2.00  962.11 -0.1353E+02  0.1056E+02
   26   12   13    2.00  435.71 -0.6129E+01  0.6289E+01
   34    9    8    2.00  713.69 -0.1004E+02  0.5840E+01
   24    4   15    2.00  516.51 -0.7265E+01  0.5176E+01
   38    3    3    2.00  411.99 -0.5795E+01  0.6245E+01
   29   18    9    2.00  521.40 -0.7334E+01  0.5084E+01
   35   15    3    2.00  831.82 -0.1170E+02  0.5744E+01
   30   21    6    2.00 2007.14 -0.2823E+02  0.1028E+02
   26   22    9    2.00  806.55 -0.1134E+02  0.1007E+02
   22   18   13    2.00  661.16 -0.9300E+01  0.5595E+01
   32   21    2    2.01 1162.37 -0.1635E+02  0.1048E+02
   34   17    3    2.01  663.80 -0.9337E+01  0.8147E+01
   38    1    3    2.01  924.20 -0.1300E+02  0.5436E+01
   29   24    4    2.01  391.25 -0.5503E+01  0.1122E+02
   28   21    8    2.01 1063.23 -0.1496E+02  0.1047E+02
   24   11   14    2.01  521.15 -0.7330E+01  0.4903E+01
   18    1   17    2.01  608.26 -0.8556E+01  0.5652E+01
   29   21    7    2.01 1032.65 -0.1453E+02  0.7617E+01
   36    9    5    2.01  535.90 -0.7538E+01  0.8924E+01
   36   12    3    2.01  440.86 -0.6201E+01  0.7869E+01
   24   15   13    2.01  688.40 -0.9683E+01  0.5334E+01
   28   25    4    2.01  525.70 -0.7394E+01  0.1109E+02
   28   17   10    2.01  925.25 -0.1301E+02  0.8149E+01
   13    3   18    2.01  540.89 -0.7608E+01  0.5269E+01
   20   17   14    2.01  734.41 -0.1033E+02  0.5332E+01
   29   23    5    2.01  868.44 -0.1222E+02  0.1122E+02
   35    7    7    2.01  581.87 -0.8184E+01  0.5879E+01
   35   10    6    2.01  832.21 -0.1171E+02  0.6118E+01
   27   22    8    2.01  495.36 -0.6968E+01  0.7419E+01
   30   13   10    2.01  638.68 -0.8984E+01  0.7193E+01
   27    8   13    2.02  521.05 -0.7329E+01  0.5930E+01
   26    3   14    2.02  476.05 -0.6696E+01  0.4163E+01
   33   19    2    2.02  540.18 -0.7598E+01  0.5670E+01
   31   21    4    2.02  560.06 -0.7878E+01  0.8207E+01
   31   19    6    2.02  868.38 -0.1221E+02  0.5585E+01
   29    7   12    2.02  582.01 -0.8186E+01  0.9216E+01
   30    9   11    2.02  597.44 -0.8404E+01  0.6805E+01
   28   24    5    2.02 1129.28 -0.1588E+02  0.1478E+02
   36    8    5    2.02  615.54 -0.8658E+01  0.6295E+01
   23   16   13    2.02  407.98 -0.5739E+01  0.6070E+01
   24   20   11    2.02 1096.79 -0.1543E+02  0.9157E+01
   28   10   12    2.02  751.15 -0.1057E+02  0.9109E+01
   25    7   14    2.02  632.81 -0.8901E+01  0.4876E+01
   25   24    8    2.02  413.58 -0.5817E+01  0.6287E+01
   33   19    1    2.02  412.27 -0.5799E+01  0.5741E+01
   27   25    5    2.02  861.00 -0.1211E+02  0.1566E+02
   27   23    7    2.02  456.49 -0.6421E+01  0.1092E+02
   33   13    7    2.02  439.09 -0.6176E+01  0.7716E+01
   12    4   18    2.02  626.51 -0.8812E+01  0.5349E+01
   29   17    9    2.03  389.84 -0.5483E+01  0.7075E+01
   23    5   15    2.03 1356.19 -0.1908E+02  0.7080E+01
   34   10    7    2.03  386.06 -0.5430E+01  0.5814E+01
   35    9    6    2.03  641.80 -0.9028E+01  0.4518E+01
   34   17    1    2.03  564.07 -0.7934E+01  0.5788E+01
   31   22    1    2.03  685.53 -0.9643E+01  0.1150E+02
   37    8    2    2.03  545.24 -0.7669E+01  0.5944E+01
   33   18    3    2.03  880.97 -0.1239E+02  0.7617E+01
   32   13    8    2.03  911.58 -0.1282E+02  0.6207E+01
   26   17   11    2.03  641.24 -0.9020E+01  0.5520E+01
   31    2   11    2.03  400.35 -0.5631E+01  0.5313E+01
   27   12   12    2.03  924.82 -0.1301E+02  0.6584E+01
   32   18    5    2.03 1278.74 -0.1799E+02  0.7873E+01
   26   14   12    2.03 1511.82 -0.2127E+02  0.8900E+01
   16    6   17    2.03  570.07 -0.8018E+01  0.3831E+01
   29   11   11    2.03  410.32 -0.5771E+01  0.6528E+01
   30   23    2    2.03 1399.92 -0.1969E+02  0.1110E+02
   36    7    5    2.03  664.98 -0.9353E+01  0.8219E+01
   30   20    6    2.03  433.21 -0.6094E+01  0.9081E+01
   25    6   14    2.03  406.46 -0.5717E+01  0.5515E+01
   30   12   10    2.03  447.45 -0.6294E+01  0.1023E+02
   15    8   17    2.03  601.94 -0.8467E+01  0.3075E+01
   20   12   15    2.03  633.63 -0.8913E+01  0.6415E+01
   31   15    8    2.03  499.75 -0.7029E+01  0.7385E+01
   36    1    6    2.03  449.38 -0.6321E+01  0.3822E+01
   32   20    2    2.03  561.98 -0.7905E+01  0.6213E+01
   28   20    8    2.03  418.10 -0.5881E+01  0.1002E+02
   35   14    2    2.04  589.01 -0.8285E+01  0.5890E+01
   34    5    8    2.04  766.25 -0.1078E+02  0.5344E+01
   28   13   11    2.04  699.02 -0.9832E+01  0.5599E+01
   37    7    2    2.04  484.65 -0.6817E+01  0.5399E+01
   37    1    4    2.04  465.88 -0.6553E+01  0.5486E+01
   30   23    1    2.04  508.35 -0.7150E+01  0.8219E+01
   35    3    7    2.04  491.47 -0.6913E+01  0.4922E+01
   31   18    6    2.04  611.74 -0.8605E+01  0.7326E+01
   31   12    9    2.04 1004.90 -0.1413E+02  0.7419E+01
   34   16    2    2.04  475.26 -0.6685E+01  0.8226E+01
   23   18   12    2.04  645.39 -0.9078E+01  0.7464E+01
   18    9   16    2.04  478.33 -0.6728E+01  0.5177E+01
   35   14    1    2.05  485.57 -0.6830E+01  0.5827E+01
   25   15   12    2.05  453.48 -0.6379E+01  0.4729E+01
   29   22    5    2.05  425.47 -0.5985E+01  0.1127E+02
   28   24    4    2.05  892.50 -0.1255E+02  0.1122E+02
   32   16    6    2.05  400.43 -0.5632E+01  0.5617E+01
   32   12    8    2.05  945.78 -0.1330E+02  0.7486E+01
   24    8   14    2.05  720.77 -0.1014E+02  0.5695E+01
   31   19    5    2.05  433.83 -0.6102E+01  0.1088E+02
   29   24    1    2.05  435.23 -0.6122E+01  0.8075E+01
   34   11    6    2.05  589.81 -0.8296E+01  0.7956E+01
   29   16    9    2.05  974.64 -0.1371E+02  0.6854E+01
   24   19   11    2.05  497.86 -0.7003E+01  0.9604E+01
   21   14   14    2.05  597.69 -0.8407E+01  0.5366E+01
   30   14    9    2.05  934.70 -0.1315E+02  0.7295E+01
   30    6   11    2.05  415.46 -0.5844E+01  0.5321E+01
   31   21    2    2.05 1381.76 -0.1944E+02  0.1140E+02
   27   17   10    2.05  649.71 -0.9139E+01  0.5032E+01
   34   16    1    2.05  761.15 -0.1071E+02  0.8185E+01
   33    1    9    2.05  610.62 -0.8589E+01  0.6314E+01
   31   16    7    2.05  476.05 -0.6696E+01  0.7645E+01
   23   22   10    2.05  688.79 -0.9688E+01  0.5977E+01
   32   18    4    2.05  572.04 -0.8046E+01  0.8111E+01
   28   12   11    2.06  614.74 -0.8647E+01  0.6478E+01
   31   11    9    2.06  753.37 -0.1060E+02  0.7564E+01
   27   24    5    2.06 1189.21 -0.1673E+02  0.1557E+02
   30   19    6    2.06  778.21 -0.1095E+02  0.5539E+01
   21    8   15    2.06  814.03 -0.1145E+02  0.5196E+01
   26   23    7    2.06 1268.85 -0.1785E+02  0.7643E+01
   34    7    7    2.06  537.19 -0.7556E+01  0.4882E+01
   34   12    5    2.06 1341.69 -0.1887E+02  0.1129E+02
   30    4   11    2.06  496.79 -0.6988E+01  0.5470E+01
   25   17   11    2.07 1235.01 -0.1737E+02  0.6549E+01
   26   24    6    2.07 1124.42 -0.1582E+02  0.1078E+02
   18   17   14    2.07  471.83 -0.6637E+01  0.5535E+01
   18    7   16    2.07  596.02 -0.8383E+01  0.4308E+01
   29   15    9    2.07  486.62 -0.6845E+01  0.7330E+01
   32   13    7    2.07  776.88 -0.1093E+02  0.7768E+01
   21    7   15    2.07  852.50 -0.1199E+02  0.4645E+01
   29    8   11    2.07  490.14 -0.6894E+01  0.9796E+01
   28   11   11    2.07  734.36 -0.1033E+02  0.9720E+01
   30   22    1    2.07 1053.81 -0.1482E+02  0.8301E+01
   34   14    3    2.07 1072.00 -0.1508E+02  0.7306E+01
   37    1    2    2.07  468.84 -0.6595E+01  0.3954E+01
   27    9   12    2.07  547.72 -0.7704E+01  0.7533E+01
   33   14    5    2.07 1553.56 -0.2185E+02  0.7975E+01
   29   20    6    2.07 1107.05 -0.1557E+02  0.1084E+02
   30   17    7    2.07  820.81 -0.1155E+02  0.5415E+01
   33   10    7    2.07  428.79 -0.6031E+01  0.6252E+01
   27   20    8    2.07  407.90 -0.5737E+01  0.1039E+02
   31   10    9    2.07  838.61 -0.1180E+02  0.1028E+02
   29   12   10    2.08  733.68 -0.1032E+02  0.7083E+01
   22   15   13    2.08  439.17 -0.6177E+01  0.5072E+01
   24    5   14    2.08  796.73 -0.1121E+02  0.4718E+01
   15    2   17    2.08  437.75 -0.6157E+01  0.4615E+01
   31   15    7    2.08  762.34 -0.1072E+02  0.7177E+01
   19   11   15    2.08  479.17 -0.6740E+01  0.5101E+01
   32   10    8    2.08  614.85 -0.8648E+01  0.7456E+01
   30   21    3    2.08 2039.50 -0.2869E+02  0.1123E+02
   34   11    5    2.08  438.57 -0.6169E+01  0.7970E+01
   33   12    6    2.08  552.63 -0.7773E+01  0.6531E+01
   34    5    7    2.08  475.07 -0.6682E+01  0.4436E+01
   27   18    9    2.08  489.19 -0.6881E+01  0.7485E+01
   36    8    2    2.08  399.38 -0.5618E+01  0.4696E+01
   24   15   12    2.08  633.76 -0.8914E+01  0.5318E+01
   35    7    5    2.08  727.99 -0.1024E+02  0.8321E+01
   35    1    6    2.09  902.72 -0.1270E+02  0.4884E+01
   28   22    5    2.09  860.53 -0.1210E+02  0.1501E+02
   30   18    6    2.09  472.26 -0.6643E+01  0.5928E+01
   28   16    9    2.09  471.12 -0.6627E+01  0.9986E+01
   32    4    9    2.09  789.16 -0.1110E+02  0.6948E+01
   12    8   17    2.09  470.89 -0.6624E+01  0.4769E+01
   29   11   10    2.09  480.46 -0.6758E+01  0.5899E+01
   34   13    3    2.09 1430.39 -0.2012E+02  0.8101E+01
   35    6    5    2.09  385.57 -0.5423E+01  0.7344E+01
   34    3    7    2.09  667.44 -0.9388E+01  0.6172E+01
   32    9    8    2.09  431.88 -0.6075E+01  0.4934E+01
   28    1   12    2.09  585.56 -0.8236E+01  0.4642E+01
   19   10   15    2.09  781.72 -0.1100E+02  0.4704E+01
   27   25    1    2.09  768.84 -0.1081E+02  0.1141E+02
   31   19    3    2.10  649.98 -0.9142E+01  0.5727E+01
   27   23    5    2.10  828.56 -0.1165E+02  0.1581E+02
   33   13    5    2.10 1393.33 -0.1960E+02  0.1069E+02
   27   22    6    2.10  503.98 -0.7089E+01  0.7866E+01
   27   12   11    2.10  491.47 -0.6913E+01  0.7874E+01
   26   22    7    2.10  743.99 -0.1046E+02  0.8744E+01
   14    2   17    2.10  484.82 -0.6819E+01  0.4580E+01
   27    7   12    2.10  925.77 -0.1302E+02  0.6076E+01
   24   10   13    2.10  414.64 -0.5832E+01  0.6153E+01
   18    2   16    2.10  487.96 -0.6864E+01  0.5669E+01
   31   17    5    2.10  412.23 -0.5798E+01  0.7995E+01
   29    5   11    2.10 1125.39 -0.1583E+02  0.6802E+01
   32   15    5    2.10  844.25 -0.1188E+02  0.1099E+02
   25    7   13    2.10  411.97 -0.5795E+01  0.6129E+01
   33    2    8    2.10  475.69 -0.6691E+01  0.6913E+01
   28    9   11    2.10  412.52 -0.5802E+01  0.6868E+01
   24   14   12    2.10  869.85 -0.1224E+02  0.6435E+01
   29   19    6    2.10  707.01 -0.9945E+01  0.9139E+01
   32    8    8    2.11  921.82 -0.1297E+02  0.6345E+01
   18   11   15    2.11  509.96 -0.7173E+01  0.4912E+01
   36    3    3    2.11  410.44 -0.5773E+01  0.5816E+01
   19   14   14    2.11  530.43 -0.7461E+01  0.5517E+01
   13    4   17    2.11  522.43 -0.7348E+01  0.4430E+01
   31   19    2    2.11  999.65 -0.1406E+02  0.6571E+01
   35    7    4    2.11  431.79 -0.6074E+01  0.5167E+01
   31   11    8    2.11  838.82 -0.1180E+02  0.8664E+01
   26   25    3    2.11  693.97 -0.9761E+01  0.1125E+02
   36    2    3    2.11  474.63 -0.6676E+01  0.4854E+01
   20   19   12    2.11  854.53 -0.1202E+02  0.5715E+01
   30    6   10    2.11  590.63 -0.8308E+01  0.5387E+01
   29    3   11    2.11  663.69 -0.9335E+01  0.4884E+01
   33   12    5    2.11  540.90 -0.7608E+01  0.1031E+02
   30   17    6    2.11  405.16 -0.5699E+01  0.7827E+01
   28   15    9    2.12  700.55 -0.9854E+01  0.5111E+01
   26    9   12    2.12  491.37 -0.6912E+01  0.6243E+01
   17   12   15    2.12  395.84 -0.5568E+01  0.4353E+01
   33   15    2    2.12  575.77 -0.8099E+01  0.6298E+01
   23   15   12    2.12 1005.75 -0.1415E+02  0.7252E+01
   28   17    8    2.12  537.57 -0.7561E+01  0.1037E+02
   32    7    8    2.12  888.12 -0.1249E+02  0.6303E+01
   26   13   11    2.12  411.76 -0.5792E+01  0.5537E+01
   28   20    6    2.12  774.14 -0.1089E+02  0.1090E+02
   23    4   14    2.12  538.29 -0.7572E+01  0.4761E+01
   31   19    1    2.12  579.46 -0.8151E+01  0.7192E+01
   35    1    5    2.12  806.82 -0.1135E+02  0.6275E+01
   35    8    3    2.12  398.56 -0.5606E+01  0.6081E+01
   30   19    4    2.12  705.20 -0.9919E+01  0.8093E+01
   28   21    5    2.12  510.93 -0.7187E+01  0.1464E+02
   26   18    9    2.12 1017.51 -0.1431E+02  0.1016E+02
   30   10    9    2.12  965.37 -0.1358E+02  0.5918E+01

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT      4696
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"       91312
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)   71460
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))    19852


 cpu time used                  4.9 sec
 elapsed wall-clock time        2.6 sec
