
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  20-Mar-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=   96 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    77.372    77.372    38.720  90.000  90.000  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.100
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/stills_2.0A_50fr_1deg/stills_2.0A_50fr_1deg_???_001.img SMV
 DATA_RANGE=      12      12
 ROTATION_AXIS=  0.000000  0.000000  1.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=      12
 X-RAY_WAVELENGTH=  1.000000
 INCIDENT_BEAM_DIRECTION=  1.000000 -0.000756  0.000692
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
 REFLECTING_RANGE_E.S.D.=   0.093
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

    1468 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE      12 ...      12
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
    1468 REFLECTIONS ACCEPTED



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
 *  31        aP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  0 -1  0  0  1  0  0  0
 *  35        mP          0.0      77.4   38.7   77.4  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  34        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  32        oP          0.0      38.7   77.4   77.4  90.0  90.0  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
 *  25        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  20        mC          0.0     109.4  109.4   38.7  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  23        oC          0.0     109.4  109.4   38.7  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
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
     9        hR        999.0      38.7   86.5  247.7 102.1  99.0 116.6    0  0  1  0  0  1 -1  0 -3 -1 -1  0
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

  *   0.00        0    6.4      18    0.00  1  0  0  0  0  1  0  0  0  0  1  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.25    77.25    38.66  90.000  90.000  90.000
 UNIT_CELL_A-AXIS=    73.732   -23.454     0.040
 UNIT_CELL_B-AXIS=    23.454    73.732     0.015
 UNIT_CELL_C-AXIS=    -0.021    -0.001    38.719
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


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM       308 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.18
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.08
 SPACE GROUP NUMBER     96
 UNIT CELL PARAMETERS     77.305    77.305    38.689  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  7.2E-02 7.2E-02 3.2E-02 0.0E+00 0.0E+00 0.0E+00
 REC. CELL PARAMETERS   0.012936  0.012936  0.025847  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS    73.669   -23.427     0.029
 COORDINATES OF UNIT CELL B-AXIS    23.427    73.669     0.012
 COORDINATES OF UNIT CELL C-AXIS    -0.016    -0.001    38.689
 CRYSTAL MOSAICITY (DEGREES)     0.093
 LAB COORDINATES OF ROTATION AXIS  0.000000  0.000000  1.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.999978  0.006317  0.001929
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1263.87   1232.89
 DETECTOR ORIGIN (PIXELS) AT                     1261.07   1242.08
 CRYSTAL TO DETECTOR DISTANCE (mm)       250.09
 LAB COORDINATES OF DETECTOR X-AXIS  0.000000  0.000000  1.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000 -1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=   96
 UNIT_CELL_CONSTANTS=    77.30    77.30    38.69  90.000  90.000  90.000


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
  0.05  0.15       109        185.936       1.172
  0.15  0.25       169        221.122       0.985
  0.25  0.35       199        182.781       1.192
  0.35  0.45       146        158.374       1.376
  0.45  0.55       170        311.043       0.700
  0.55  0.65       196        198.098       1.100
  0.65  0.75       186        201.614       1.081
  0.75  0.85       172        190.619       1.143
  0.85  0.95       121        336.153       0.648
 
 -0.05  0.95      1468        217.846       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       10
 DEGREES OF FREEDOM OF CHI^2 FIT               237.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.029
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=    11.1 XMAX=    11.9 NXBIN=    1
 YMIN= 0.00167 YMAX= 0.59677 NYBIN=   10
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        516


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED        9
 DEGREES OF FREEDOM OF CHI^2 FIT               237.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.026
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=    31.4 XMAX=  2437.2 NXBIN=    3
 YMIN=     7.9 YMAX=  2518.1 NYBIN=    3
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        516


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED       13
 DEGREES OF FREEDOM OF CHI^2 FIT               237.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      1.023
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=    11.1 XMAX=    11.9 NXBIN=    1
 DETECTOR_SURFACE_POSITION=    1234    1263
 DETECTOR_SURFACE_POSITION=    1643    1689
 DETECTOR_SURFACE_POSITION=     826    1689
 DETECTOR_SURFACE_POSITION=     826     837
 DETECTOR_SURFACE_POSITION=    1643     837
 DETECTOR_SURFACE_POSITION=    2159    1663
 DETECTOR_SURFACE_POSITION=    1617    2228
 DETECTOR_SURFACE_POSITION=     851    2228
 DETECTOR_SURFACE_POSITION=     310    1663
 DETECTOR_SURFACE_POSITION=     310     863
 DETECTOR_SURFACE_POSITION=     851     298
 DETECTOR_SURFACE_POSITION=    1617     298
 DETECTOR_SURFACE_POSITION=    2159     863
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS        516


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
 1.014E+00  1.794E-03   23.45


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

   24.446   3.396    19.59   1.87      2.53      2.52      64      88      20
    3.396   2.413    14.62   1.36      3.69      3.23     122     163      18
    2.413   1.974     7.73   0.97      5.92      5.92     188     218       2
    1.974   1.711    -0.01   0.69    -99.90    -99.90     186     235       2
    1.711   1.531     0.00   0.76   2610.80   3024.84     222     276       2
    1.531   1.398    -0.05   0.62    -99.90    -99.90     242     300       0
    1.398   1.294     0.06   0.76    -99.90    -99.90      94     141       2
   --------------------------------------------------------------------------
   24.446   1.294     4.07   0.88      5.44      5.72    1118    1421      46
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=      12      12)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

   46    0    0     1.681   0.1553E+01  0.1153E+02         0.13           1 
   47    0    0     1.645  -0.1002E+02  0.1140E+02        -0.88           1*

 AVERAGE INTENSITY FOR       1 REFLECTIONS WHICH SHOULD
 BE SYSTEMATICALLY ABSENT IS    -4.4% OF MEAN INTENSITY


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

     3.84          65      42      1285        3.3%       2.6%      2.5%       46   24.54      3.7%   100.0*     0    0.000       0
     2.73         106      67      2098        3.2%       3.0%      2.9%       78   20.94      4.2%    99.6*     0    0.000       0
     2.24         158      95      2665        3.6%       4.2%      3.9%      126   15.43      6.0%    99.4*     0    0.000       0
     1.94         165      93      3101        3.0%       7.9%      8.1%      144    6.59     11.1%    98.8*     0    0.000       0
     1.73         180     108      3481        3.1%     -99.9%    -99.9%      144    0.00    -99.9%    -1.7      0    0.000       0
     1.58         215     127      3848        3.3%     857.9%    966.6%      176    0.03   1213.2%    -8.7      0    0.000       0
     1.47         232     142      4143        3.4%     -99.9%    -99.9%      180    0.00    -99.9%    17.1      0    0.000       0
     1.37         218     128      4429        2.9%    2300.7%   2794.0%      180    0.03   3253.6%    -1.0      0    0.000       0
     1.29          82      60      4740        1.3%     -99.9%    -99.9%       44    0.00    -99.9%     1.7      0    0.000       0
    total        1421     862     29790        2.9%       5.4%      5.7%     1118    5.23      7.7%    99.9*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES      1468
 NUMBER OF REJECTED MISFITS                                46
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    1
 NUMBER OF ACCEPTED OBSERVATIONS                         1421
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    862



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

   WILSON LINE (using all data) : A=   7.607 B=  31.243 CORRELATION=  0.81
      #      RES      SS        <I>       log(<I>)       BO
      17     7.250   0.005  2.1316E+03       7.665      -6.1
      21     4.617   0.012  2.7075E+03       7.904     -12.7
      29     3.502   0.020  1.1918E+03       7.083      12.8
      32     2.995   0.028  5.5974E+02       6.327      22.9
      40     2.648   0.036  5.2671E+02       6.267      18.8
      51     2.392   0.044  3.1806E+02       5.762      21.1
      41     2.203   0.051  1.9702E+02       5.283      22.6
      45     2.055   0.059  8.2382E+01       4.411      27.0
      46     1.930   0.067  3.5052E-01      -1.048      64.5
      59     1.828   0.075 -1.8242E-01    -999.999   -1000.0
      52     1.739   0.083  8.0875E-01      -0.212      47.3
      55     1.665   0.090  5.2450E-01      -0.645      45.7
      59     1.595   0.098 -9.7304E-02    -999.999   -1000.0
      72     1.535   0.106 -7.1982E-01    -999.999   -1000.0
      77     1.479   0.114 -8.3273E-01    -999.999   -1000.0
      58     1.431   0.122 -1.5265E-01    -999.999   -1000.0
      50     1.389   0.130 -1.6804E-01    -999.999   -1000.0
      40     1.351   0.137  9.2373E-01      -0.079      28.1
      18     1.314   0.145  5.1612E+00       1.641      20.6


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       5     7.250        1.842        1.211        0.629
       3     4.617        1.243        0.753        0.339
       3     3.502        0.334        0.082        0.014
       4     2.995        1.021        0.439        0.143
       3     2.648        0.291        0.071        0.013
       5     2.392        0.761        0.350        0.123
       2     2.203        1.453        0.621        0.192
       8     2.055        0.292        0.107        0.031
      11     1.828       15.480       27.908       52.974
      11     1.665       56.389       55.152      771.140
       3     1.595       17.679       25.780       33.902
       3     1.431       32.852       48.509      142.165
       2     1.351      339.170     2933.519    18515.066
      88   overall       17.589       71.735      523.714


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

      12     7.250        0.438        0.298        0.163
      18     4.617        1.096        1.225        1.162
      26     3.502        0.865        0.675        0.490
      28     2.995        0.914        0.865        0.745
      37     2.648        1.161        1.339        1.428
      46     2.392        0.926        0.750        0.508
      39     2.203        0.820        0.628        0.418
      37     2.055        1.510        2.110        2.709
      39     1.930      273.632      529.136    31233.063
      43     1.739       43.228       57.926      837.680
      44     1.665      114.745     -119.907     9686.324
      18     1.314        6.993       22.812       79.136
     774   overall       25.360       24.477     2208.069

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       5     7.250  0.00  0.00  0.00  0.00  0.00  0.00  0.33  0.32  0.30  0.29
       3     4.617  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.98
       3     3.502  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       4     2.995  0.00  0.00  0.00  0.00  0.00  0.45  0.42  0.40  0.38  0.37
       3     2.648  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  1.01  0.98
       5     2.392  0.00  0.00  0.00  0.00  0.38  1.07  1.00  0.95  0.91  0.88
       2     2.203  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8     2.055  2.01  1.45  1.50  1.32  1.20  1.11  1.05  0.99  1.14  1.10
      11     1.828  2.56  1.84  1.53  1.34  1.22  1.13  1.07  1.01  0.97  0.93
      11     1.665  1.83  1.31  1.09  0.96  0.87  0.81  0.76  0.72  0.69  0.67
       3     1.595  2.68  1.93  1.60  1.41  1.28  1.19  1.12  1.06  1.01  0.98
       3     1.431  1.34  0.96  0.80  0.70  0.64  0.59  0.56  0.53  0.51  0.49
       2     1.351  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
      88   overall  1.01  0.72  0.63  0.55  0.52  0.55  0.53  0.51  0.52  0.52


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

      12     7.250  2.63  2.76  1.93  1.52  1.91  1.85  1.66  1.51  1.40  1.32
      18     4.617  1.75  1.84  1.29  1.01  1.13  1.23  1.10  1.01  1.22  1.14
      26     3.502  0.81  0.85  0.59  0.70  0.68  0.77  0.92  0.98  0.91  0.85
      28     2.995  1.13  1.58  1.10  1.30  1.18  1.27  1.21  1.17  1.08  1.07
      37     2.648  0.85  1.19  1.15  1.15  1.03  0.96  1.07  1.13  1.05  1.07
      46     2.392  1.60  1.32  1.09  0.99  0.99  1.11  1.04  1.07  1.10  1.03
      39     2.203  1.62  1.27  0.99  0.86  1.04  0.97  0.87  0.93  0.99  1.10
      37     2.055  1.42  1.04  1.25  1.23  1.10  1.02  0.97  0.98  0.96  1.03
      39     1.930  4.31  2.26  1.58  1.24  1.04  0.91  0.81  0.75  0.73  0.69
      43     1.739  5.13  2.69  1.88  1.48  1.24  1.08  0.97  0.93  0.86  0.81
      44     1.665  4.78  2.51  1.75  1.45  1.21  1.06  0.95  0.87  0.80  0.76
      18     1.314  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.97
     774   overall  1.33  0.89  0.68  0.60  0.55  0.53  0.50  0.49  0.48  0.48


 List of    168 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   36   28    9    1.58  242.79 -0.2362E+02  0.1587E+02 "alien"
   36   32   13    1.41  205.95 -0.3144E+02  0.1482E+02 "alien"
   35    5   17    1.58  195.68 -0.1904E+02  0.1110E+02 "alien"
   20   13   21    1.60  130.85 -0.1273E+02  0.7133E+01 "alien"
   37   31   12    1.43  128.65 -0.1964E+02  0.1507E+02 "alien"
   48    1    4    1.59  124.26 -0.1209E+02  0.8464E+01 "alien"
   40    1    5    1.87  121.63 -0.2219E+02  0.1127E+02 "alien"
   19    9   22    1.59  102.77 -0.1000E+02  0.5884E+01 "alien"
   24    9   21    1.57  101.95 -0.9920E+01  0.8812E+01 "alien"
   23   21   19    1.57   92.15 -0.8966E+01  0.6551E+01 "alien"
   25   15   23    1.42   91.40 -0.1395E+02  0.9813E+01 "alien"
   34   22    4    1.87   90.59 -0.1653E+02  0.9371E+01 "alien"
   48    5   13    1.41   87.24 -0.1332E+02  0.9473E+01 "alien"
   33   30   15    1.44   84.79 -0.1294E+02  0.8695E+01 "alien"
   33   23    7    1.82   83.67 -0.1526E+02  0.1362E+02 "alien"
   43    3   12    1.57   81.84 -0.7964E+01  0.8051E+01 "alien"
   37    4   15    1.62   78.10 -0.7600E+01  0.8464E+01 "alien"
   39    4   15    1.57   77.54 -0.7545E+01  0.1126E+02 "alien"
   31   30   18    1.38   76.89 -0.1292E+02  0.9175E+01 "alien"
   53    3    3    1.45   75.47 -0.1152E+02  0.1131E+02 "alien"
   22    9   21    1.60   73.31 -0.7134E+01  0.8345E+01 "alien"
   33   32   16    1.38   72.74 -0.1222E+02  0.1536E+02 "alien"
   50    5   11    1.41   71.97 -0.1209E+02  0.1652E+02 "alien"
   37    5   16    1.57   68.74 -0.6689E+01  0.1324E+02 "alien"
   42    1    4    1.81   68.67 -0.1253E+02  0.8218E+01 "alien"
   36   32   14    1.39   63.51 -0.1067E+02  0.1671E+02 "alien"
   35   24    4    1.79   61.43 -0.1121E+02  0.1485E+02 "alien"
   46    7   16    1.37   60.93 -0.1024E+02  0.8833E+01 "alien"
   40    7   18    1.43   60.86 -0.9289E+01  0.1525E+02 "alien"
   32   26   13    1.59   60.73 -0.5909E+01  0.8370E+01 "alien"
   35   10   21    1.39   60.51 -0.1017E+02  0.8204E+01 "alien"
   50    4   10    1.43   58.66 -0.8954E+01  0.9135E+01 "alien"
   37   27    7    1.61   58.35 -0.5678E+01  0.8819E+01 "alien"
   25   18   12    1.98   58.29  0.2043E+02  0.1164E+02 "alien"
   28   21   12    1.82   56.85 -0.1037E+02  0.7590E+01 "alien"
   43    2   11    1.60   55.42 -0.5393E+01  0.9244E+01 "alien"
   40    1    6    1.85   52.91 -0.9652E+01  0.8074E+01 "alien"
   30   26   20    1.37   52.59 -0.8837E+01  0.1427E+02 "alien"
   41    3   13    1.59   52.41 -0.5100E+01  0.7986E+01 "alien"
   53    3    4    1.44   52.14 -0.7959E+01  0.9217E+01 "alien"
   54    4    7    1.38   48.69 -0.8182E+01  0.9403E+01 "alien"
   30   22   11    1.79   48.45 -0.8839E+01  0.7886E+01 "alien"
   43    1    2    1.79   47.97 -0.8750E+01  0.8128E+01 "alien"
   52    3    6    1.45   47.68 -0.7279E+01  0.8851E+01 "alien"
   52    4    8    1.42   47.53 -0.7255E+01  0.9168E+01 "alien"
   22   18   20    1.58   46.15 -0.4491E+01  0.6020E+01 "alien"
   30   21    9    1.89   44.58  0.1562E+02  0.1274E+02 "alien"
   24   14   23    1.44   43.15 -0.6587E+01  0.6266E+01 "alien"
   36    1   11    1.83   42.81 -0.7810E+01  0.1382E+02 "alien"
   47    1    6    1.59   41.05 -0.3994E+01  0.1198E+02 "alien"
   31    3   15    1.79   40.32 -0.7355E+01  0.6873E+01 "alien"
   38   33   12    1.39   40.13 -0.6744E+01  0.9393E+01 "alien"
   32   11   22    1.39   40.03 -0.6727E+01  0.1284E+02 "alien"
   38   28    6    1.59   40.01 -0.3893E+01  0.1248E+02 "alien"
   40   32    8    1.44   38.72 -0.5910E+01  0.9158E+01 "alien"
   36   31   13    1.43   37.62 -0.5743E+01  0.9328E+01 "alien"
   33   25   10    1.68   34.80  0.1825E+02  0.9981E+01 "alien"
   42    1    1    1.84   33.69 -0.6146E+01  0.9210E+01 "alien"
   38    9   20    1.38   33.51 -0.5632E+01  0.9164E+01 "alien"
   26   19   22    1.42   33.29 -0.5082E+01  0.9681E+01 "alien"
   34   22    1    1.91   32.68  0.1146E+02  0.8689E+01 "alien"
   34   31   15    1.41   31.07 -0.5221E+01  0.8950E+01 "alien"
   27   24   20    1.44   29.60 -0.4518E+01  0.1203E+02 "alien"
   50    5   12    1.39   28.61 -0.4808E+01  0.9534E+01 "alien"
   32    1   12    1.93   28.59  0.1002E+02  0.1319E+02 "alien"
   30    3   15    1.82   28.37 -0.5176E+01  0.1054E+02 "alien"
   37   26    5    1.67   28.26  0.1482E+02  0.8688E+01 "alien"
   45    7   16    1.39   28.15 -0.4731E+01  0.9409E+01 "alien"
   51    4    9    1.43   28.01 -0.4275E+01  0.9026E+01 "alien"
   36    9   20    1.42   27.57 -0.4209E+01  0.1139E+02 "alien"
   37   30   10    1.50   27.39 -0.2281E+02  0.9357E+01 "alien"
   36    1    8    1.96   26.97  0.9454E+01  0.8312E+01 "alien"
   48    1    3    1.60   26.91 -0.2618E+01  0.8743E+01 "alien"
   45    5   14    1.45   26.72 -0.4079E+01  0.8483E+01 "alien"
   33   22    5    1.89   26.71  0.9362E+01  0.9747E+01 "alien"
   41   32    6    1.45   26.60 -0.4060E+01  0.9405E+01 "alien"
   22   21   17    1.70   26.50  0.1390E+02  0.9201E+01 "alien"
   28   26   19    1.44   26.47 -0.4041E+01  0.8491E+01 "alien"
   32   30   16    1.42   26.29 -0.4013E+01  0.8663E+01 "alien"
   46    5   14    1.43   25.77 -0.3934E+01  0.1345E+02 "alien"
   27   21   13    1.80   25.60 -0.4671E+01  0.7358E+01 "alien"
   52    3    5    1.46   25.48 -0.2122E+02  0.1004E+02 "alien"
   46    6   15    1.40   24.58 -0.4131E+01  0.9336E+01 "alien"
   21   15   21    1.57   24.36 -0.2371E+01  0.8765E+01 "alien"
   20   13   24    1.44   23.68 -0.3615E+01  0.7085E+01 "alien"
   14    8   19    1.87   23.60 -0.4305E+01  0.4967E+01 "alien"
   30   29   17    1.44   23.56 -0.3596E+01  0.8430E+01 "alien"
   37   27    6    1.63   22.87  0.1199E+02  0.1209E+02 "alien"
   38   28    5    1.60   22.66 -0.2205E+01  0.8776E+01 "alien"
   25   24   17    1.59   22.20 -0.2160E+01  0.7867E+01 "alien"
   23   21   16    1.73   22.18  0.1794E+02  0.1192E+02 "alien"
   38   29    7    1.55   22.17 -0.1596E+02  0.9010E+01 "alien"
   29   20   10    1.91   22.07  0.7736E+01  0.1300E+02 "alien"
   29   24   14    1.65   22.02  0.1155E+02  0.8052E+01 "alien"
   41    1    5    1.83   21.97 -0.4009E+01  0.8445E+01 "alien"
   39   28    3    1.60   21.39 -0.2081E+01  0.8737E+01 "alien"
   27   22   21    1.42   20.77 -0.3171E+01  0.7839E+01 "alien"
   34    1   12    1.86   20.76 -0.3788E+01  0.7282E+01 "alien"
   43    1    9    1.66   20.42  0.1071E+02  0.8319E+01 "alien"
   30   20    8    1.96   20.30  0.7117E+01  0.8487E+01 "alien"
   50    2    4    1.53   20.08 -0.1445E+02  0.1551E+02 "alien"
   34    5   17    1.60   19.83 -0.1929E+01  0.6944E+01 "alien"
   48    2    8    1.53   19.19 -0.1381E+02  0.8743E+01 "alien"
   54    3    4    1.41   19.17 -0.2926E+01  0.9291E+01 "alien"
   29   27   19    1.41   19.10 -0.3209E+01  0.8398E+01 "alien"
   40   31    6    1.49   19.07 -0.1588E+02  0.9284E+01 "alien"
   20   17   15    1.94   19.01  0.6665E+01  0.6440E+01 "alien"
   39   30    6    1.53   18.55 -0.1335E+02  0.9920E+01 "alien"
   37   29   10    1.51   18.46 -0.1329E+02  0.1465E+02 "alien"
   28   27   17    1.50   18.42 -0.1534E+02  0.7914E+01 "alien"
   20    9   21    1.63   18.17  0.9529E+01  0.6774E+01 "alien"
   34    3   15    1.70   18.04  0.1459E+02  0.1174E+02 "alien"
   34   11   22    1.36   17.45  0.1612E+02  0.7694E+01 "alien"
   32   29   15    1.47   17.17 -0.1430E+02  0.8826E+01 "alien"
   39    1    6    1.89   17.11  0.5998E+01  0.8068E+01 "alien"
   50    2    6    1.50   17.07 -0.1421E+02  0.8673E+01 "alien"
   34    4   16    1.65   17.05  0.8943E+01  0.7608E+01 "alien"
   30   25   21    1.35   16.97  0.1567E+02  0.1213E+02 "alien"
   35   27   10    1.59   16.80 -0.1634E+01  0.8653E+01 "alien"
   31   27   15    1.52   16.57 -0.1193E+02  0.1401E+02 "alien"
   37    1    7    1.95   16.36  0.5734E+01  0.8260E+01 "alien"
   31    6   18    1.62   16.12 -0.1569E+01  0.1186E+02 "alien"
   24   22   19    1.55   16.04 -0.1155E+02  0.7066E+01 "alien"
   40   30    5    1.52   15.94 -0.1147E+02  0.9152E+01 "alien"
   41   32    4    1.47   15.75 -0.1312E+02  0.9628E+01 "alien"
   29   29   18    1.42   15.74 -0.7185E+01  0.1426E+02 "alien"
   38    4   15    1.59   15.69 -0.1527E+01  0.7693E+01 "alien"
   26   21   21    1.44   15.32 -0.2339E+01  0.8185E+01 "alien"
   33   33   17    1.34   15.32  0.4101E+02  0.1444E+02 "alien"
   54    3    0    1.43   15.12 -0.6593E+01  0.1327E+02 "alien"
   35    2   13    1.77   14.88  0.1204E+02  0.7642E+01 "alien"
   39    2    3    1.96   14.84  0.5200E+01  0.1000E+02 "alien"
   33   11   22    1.38   14.83 -0.2492E+01  0.8083E+01 "alien"
   33   27   13    1.55   14.71 -0.1059E+02  0.1321E+02 "alien"
   23   23   18    1.59   14.67 -0.3950E+01  0.1027E+02 "alien"
   39   31    9    1.46   14.60 -0.1216E+02  0.1016E+02 "alien"
   50    7   14    1.34   14.44  0.1334E+02  0.9759E+01 "alien"
   34   24    8    1.73   14.24  0.1151E+02  0.8425E+01 "alien"
   32   23   10    1.75   14.21  0.1149E+02  0.1523E+02 "alien"
   52    2    2    1.48   13.85 -0.1153E+02  0.8806E+01 "alien"
   35   11   22    1.35   13.75  0.1271E+02  0.1099E+02 "alien"
   32   12   23    1.35   13.66  0.1262E+02  0.1032E+02 "alien"
   53    7   13    1.30   13.53  0.6984E+02  0.1455E+02 "alien"
   40   30    2    1.54   13.47 -0.9694E+01  0.1014E+02 "alien"
   33   25   11    1.65   13.34  0.6995E+01  0.1330E+02 "alien"
   43    1    0    1.80   13.30 -0.6091E+01  0.1184E+02 "alien"
   39   33   10    1.41   13.14 -0.2208E+01  0.1044E+02 "alien"
   31   21    8    1.90   12.49  0.4380E+01  0.8098E+01 "alien"
   29   20   23    1.34   11.94  0.1103E+02  0.1091E+02 "alien"
   53    6   11    1.34   11.93  0.1102E+02  0.9655E+01 "alien"
   38    1    7    1.91   11.87  0.4161E+01  0.8189E+01 "alien"
   39    2    1    1.98   11.58  0.4057E+01  0.8690E+01 "alien"
   46    0    1    1.68   11.49  0.1338E+02  0.8486E+01 "alien"
   17   17   17    1.86   11.42 -0.4606E+01  0.7697E+01 "alien"
   35   24    6    1.75   11.37  0.9198E+01  0.8640E+01 "alien"
   40    1   10    1.73   11.24  0.9092E+01  0.7889E+01 "alien"
   35    6   18    1.53   11.13 -0.8015E+01  0.1197E+02 "alien"
   49    3    9    1.48   11.11 -0.9253E+01  0.8664E+01 "alien"
   30    1   13    1.95   11.10  0.3892E+01  0.1080E+02 "alien"
   33    5   17    1.62   10.98 -0.1069E+01  0.6764E+01 "alien"
   24   10   22    1.51   10.82 -0.7787E+01  0.7283E+01 "alien"
   31   22   10    1.80   10.53 -0.1921E+01  0.1283E+02 "alien"
   36   24    3    1.77   10.38  0.8392E+01  0.9383E+01 "alien"
   40   29    1    1.56   10.37 -0.7466E+01  0.8595E+01 "alien"
   49    7   15    1.34   10.17  0.9399E+01  0.1454E+02 "alien"
   14   10   18    1.94   10.07  0.3531E+01  0.5798E+01 "alien"
   24   13   24    1.40   10.06 -0.1691E+01  0.7705E+01 "alien"
   26   22   15    1.70   10.06  0.8134E+01  0.7823E+01 "alien"


 List of    168 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   53    7   13    1.30   13.53  0.6984E+02  0.1455E+02
   29   20   23    1.34   11.94  0.1103E+02  0.1091E+02
   49    7   15    1.34   10.17  0.9399E+01  0.1454E+02
   33   33   17    1.34   15.32  0.4101E+02  0.1444E+02
   50    7   14    1.34   14.44  0.1334E+02  0.9759E+01
   53    6   11    1.34   11.93  0.1102E+02  0.9655E+01
   30   25   21    1.35   16.97  0.1567E+02  0.1213E+02
   32   12   23    1.35   13.66  0.1262E+02  0.1032E+02
   35   11   22    1.35   13.75  0.1271E+02  0.1099E+02
   34   11   22    1.36   17.45  0.1612E+02  0.7694E+01
   46    7   16    1.37   60.93 -0.1024E+02  0.8833E+01
   30   26   20    1.37   52.59 -0.8837E+01  0.1427E+02
   31   30   18    1.38   76.89 -0.1292E+02  0.9175E+01
   33   11   22    1.38   14.83 -0.2492E+01  0.8083E+01
   33   32   16    1.38   72.74 -0.1222E+02  0.1536E+02
   54    4    7    1.38   48.69 -0.8182E+01  0.9403E+01
   38    9   20    1.38   33.51 -0.5632E+01  0.9164E+01
   38   33   12    1.39   40.13 -0.6744E+01  0.9393E+01
   36   32   14    1.39   63.51 -0.1067E+02  0.1671E+02
   50    5   12    1.39   28.61 -0.4808E+01  0.9534E+01
   45    7   16    1.39   28.15 -0.4731E+01  0.9409E+01
   35   10   21    1.39   60.51 -0.1017E+02  0.8204E+01
   32   11   22    1.39   40.03 -0.6727E+01  0.1284E+02
   46    6   15    1.40   24.58 -0.4131E+01  0.9336E+01
   24   13   24    1.40   10.06 -0.1691E+01  0.7705E+01
   34   31   15    1.41   31.07 -0.5221E+01  0.8950E+01
   29   27   19    1.41   19.10 -0.3209E+01  0.8398E+01
   39   33   10    1.41   13.14 -0.2208E+01  0.1044E+02
   50    5   11    1.41   71.97 -0.1209E+02  0.1652E+02
   48    5   13    1.41   87.24 -0.1332E+02  0.9473E+01
   36   32   13    1.41  205.95 -0.3144E+02  0.1482E+02
   54    3    4    1.41   19.17 -0.2926E+01  0.9291E+01
   29   29   18    1.42   15.74 -0.7185E+01  0.1426E+02
   52    4    8    1.42   47.53 -0.7255E+01  0.9168E+01
   36    9   20    1.42   27.57 -0.4209E+01  0.1139E+02
   27   22   21    1.42   20.77 -0.3171E+01  0.7839E+01
   26   19   22    1.42   33.29 -0.5082E+01  0.9681E+01
   25   15   23    1.42   91.40 -0.1395E+02  0.9813E+01
   32   30   16    1.42   26.29 -0.4013E+01  0.8663E+01
   40    7   18    1.43   60.86 -0.9289E+01  0.1525E+02
   51    4    9    1.43   28.01 -0.4275E+01  0.9026E+01
   36   31   13    1.43   37.62 -0.5743E+01  0.9328E+01
   54    3    0    1.43   15.12 -0.6593E+01  0.1327E+02
   46    5   14    1.43   25.77 -0.3934E+01  0.1345E+02
   50    4   10    1.43   58.66 -0.8954E+01  0.9135E+01
   37   31   12    1.43  128.65 -0.1964E+02  0.1507E+02
   27   24   20    1.44   29.60 -0.4518E+01  0.1203E+02
   28   26   19    1.44   26.47 -0.4041E+01  0.8491E+01
   30   29   17    1.44   23.56 -0.3596E+01  0.8430E+01
   33   30   15    1.44   84.79 -0.1294E+02  0.8695E+01
   24   14   23    1.44   43.15 -0.6587E+01  0.6266E+01
   53    3    4    1.44   52.14 -0.7959E+01  0.9217E+01
   40   32    8    1.44   38.72 -0.5910E+01  0.9158E+01
   26   21   21    1.44   15.32 -0.2339E+01  0.8185E+01
   20   13   24    1.44   23.68 -0.3615E+01  0.7085E+01
   52    3    6    1.45   47.68 -0.7279E+01  0.8851E+01
   53    3    3    1.45   75.47 -0.1152E+02  0.1131E+02
   41   32    6    1.45   26.60 -0.4060E+01  0.9405E+01
   45    5   14    1.45   26.72 -0.4079E+01  0.8483E+01
   52    3    5    1.46   25.48 -0.2122E+02  0.1004E+02
   39   31    9    1.46   14.60 -0.1216E+02  0.1016E+02
   41   32    4    1.47   15.75 -0.1312E+02  0.9628E+01
   32   29   15    1.47   17.17 -0.1430E+02  0.8826E+01
   49    3    9    1.48   11.11 -0.9253E+01  0.8664E+01
   52    2    2    1.48   13.85 -0.1153E+02  0.8806E+01
   40   31    6    1.49   19.07 -0.1588E+02  0.9284E+01
   37   30   10    1.50   27.39 -0.2281E+02  0.9357E+01
   28   27   17    1.50   18.42 -0.1534E+02  0.7914E+01
   50    2    6    1.50   17.07 -0.1421E+02  0.8673E+01
   37   29   10    1.51   18.46 -0.1329E+02  0.1465E+02
   24   10   22    1.51   10.82 -0.7787E+01  0.7283E+01
   40   30    5    1.52   15.94 -0.1147E+02  0.9152E+01
   31   27   15    1.52   16.57 -0.1193E+02  0.1401E+02
   50    2    4    1.53   20.08 -0.1445E+02  0.1551E+02
   39   30    6    1.53   18.55 -0.1335E+02  0.9920E+01
   48    2    8    1.53   19.19 -0.1381E+02  0.8743E+01
   35    6   18    1.53   11.13 -0.8015E+01  0.1197E+02
   40   30    2    1.54   13.47 -0.9694E+01  0.1014E+02
   24   22   19    1.55   16.04 -0.1155E+02  0.7066E+01
   33   27   13    1.55   14.71 -0.1059E+02  0.1321E+02
   38   29    7    1.55   22.17 -0.1596E+02  0.9010E+01
   40   29    1    1.56   10.37 -0.7466E+01  0.8595E+01
   39    4   15    1.57   77.54 -0.7545E+01  0.1126E+02
   43    3   12    1.57   81.84 -0.7964E+01  0.8051E+01
   21   15   21    1.57   24.36 -0.2371E+01  0.8765E+01
   24    9   21    1.57  101.95 -0.9920E+01  0.8812E+01
   37    5   16    1.57   68.74 -0.6689E+01  0.1324E+02
   23   21   19    1.57   92.15 -0.8966E+01  0.6551E+01
   22   18   20    1.58   46.15 -0.4491E+01  0.6020E+01
   35    5   17    1.58  195.68 -0.1904E+02  0.1110E+02
   36   28    9    1.58  242.79 -0.2362E+02  0.1587E+02
   32   26   13    1.59   60.73 -0.5909E+01  0.8370E+01
   19    9   22    1.59  102.77 -0.1000E+02  0.5884E+01
   38   28    6    1.59   40.01 -0.3893E+01  0.1248E+02
   48    1    4    1.59  124.26 -0.1209E+02  0.8464E+01
   41    3   13    1.59   52.41 -0.5100E+01  0.7986E+01
   38    4   15    1.59   15.69 -0.1527E+01  0.7693E+01
   25   24   17    1.59   22.20 -0.2160E+01  0.7867E+01
   47    1    6    1.59   41.05 -0.3994E+01  0.1198E+02
   35   27   10    1.59   16.80 -0.1634E+01  0.8653E+01
   23   23   18    1.59   14.67 -0.3950E+01  0.1027E+02
   39   28    3    1.60   21.39 -0.2081E+01  0.8737E+01
   48    1    3    1.60   26.91 -0.2618E+01  0.8743E+01
   43    2   11    1.60   55.42 -0.5393E+01  0.9244E+01
   34    5   17    1.60   19.83 -0.1929E+01  0.6944E+01
   20   13   21    1.60  130.85 -0.1273E+02  0.7133E+01
   38   28    5    1.60   22.66 -0.2205E+01  0.8776E+01
   22    9   21    1.60   73.31 -0.7134E+01  0.8345E+01
   37   27    7    1.61   58.35 -0.5678E+01  0.8819E+01
   31    6   18    1.62   16.12 -0.1569E+01  0.1186E+02
   37    4   15    1.62   78.10 -0.7600E+01  0.8464E+01
   33    5   17    1.62   10.98 -0.1069E+01  0.6764E+01
   37   27    6    1.63   22.87  0.1199E+02  0.1209E+02
   20    9   21    1.63   18.17  0.9529E+01  0.6774E+01
   29   24   14    1.65   22.02  0.1155E+02  0.8052E+01
   33   25   11    1.65   13.34  0.6995E+01  0.1330E+02
   34    4   16    1.65   17.05  0.8943E+01  0.7608E+01
   43    1    9    1.66   20.42  0.1071E+02  0.8319E+01
   37   26    5    1.67   28.26  0.1482E+02  0.8688E+01
   46    0    1    1.68   11.49  0.1338E+02  0.8486E+01
   33   25   10    1.68   34.80  0.1825E+02  0.9981E+01
   22   21   17    1.70   26.50  0.1390E+02  0.9201E+01
   34    3   15    1.70   18.04  0.1459E+02  0.1174E+02
   26   22   15    1.70   10.06  0.8134E+01  0.7823E+01
   40    1   10    1.73   11.24  0.9092E+01  0.7889E+01
   23   21   16    1.73   22.18  0.1794E+02  0.1192E+02
   34   24    8    1.73   14.24  0.1151E+02  0.8425E+01
   32   23   10    1.75   14.21  0.1149E+02  0.1523E+02
   35   24    6    1.75   11.37  0.9198E+01  0.8640E+01
   36   24    3    1.77   10.38  0.8392E+01  0.9383E+01
   35    2   13    1.77   14.88  0.1204E+02  0.7642E+01
   31    3   15    1.79   40.32 -0.7355E+01  0.6873E+01
   30   22   11    1.79   48.45 -0.8839E+01  0.7886E+01
   43    1    2    1.79   47.97 -0.8750E+01  0.8128E+01
   35   24    4    1.79   61.43 -0.1121E+02  0.1485E+02
   43    1    0    1.80   13.30 -0.6091E+01  0.1184E+02
   27   21   13    1.80   25.60 -0.4671E+01  0.7358E+01
   31   22   10    1.80   10.53 -0.1921E+01  0.1283E+02
   42    1    4    1.81   68.67 -0.1253E+02  0.8218E+01
   33   23    7    1.82   83.67 -0.1526E+02  0.1362E+02
   30    3   15    1.82   28.37 -0.5176E+01  0.1054E+02
   28   21   12    1.82   56.85 -0.1037E+02  0.7590E+01
   41    1    5    1.83   21.97 -0.4009E+01  0.8445E+01
   36    1   11    1.83   42.81 -0.7810E+01  0.1382E+02
   42    1    1    1.84   33.69 -0.6146E+01  0.9210E+01
   40    1    6    1.85   52.91 -0.9652E+01  0.8074E+01
   34    1   12    1.86   20.76 -0.3788E+01  0.7282E+01
   17   17   17    1.86   11.42 -0.4606E+01  0.7697E+01
   34   22    4    1.87   90.59 -0.1653E+02  0.9371E+01
   14    8   19    1.87   23.60 -0.4305E+01  0.4967E+01
   40    1    5    1.87  121.63 -0.2219E+02  0.1127E+02
   33   22    5    1.89   26.71  0.9362E+01  0.9747E+01
   39    1    6    1.89   17.11  0.5998E+01  0.8068E+01
   30   21    9    1.89   44.58  0.1562E+02  0.1274E+02
   31   21    8    1.90   12.49  0.4380E+01  0.8098E+01
   34   22    1    1.91   32.68  0.1146E+02  0.8689E+01
   38    1    7    1.91   11.87  0.4161E+01  0.8189E+01
   29   20   10    1.91   22.07  0.7736E+01  0.1300E+02
   32    1   12    1.93   28.59  0.1002E+02  0.1319E+02
   14   10   18    1.94   10.07  0.3531E+01  0.5798E+01
   20   17   15    1.94   19.01  0.6665E+01  0.6440E+01
   30    1   13    1.95   11.10  0.3892E+01  0.1080E+02
   37    1    7    1.95   16.36  0.5734E+01  0.8260E+01
   39    2    3    1.96   14.84  0.5200E+01  0.1000E+02
   30   20    8    1.96   20.30  0.7117E+01  0.8487E+01
   36    1    8    1.96   26.97  0.9454E+01  0.8312E+01
   39    2    1    1.98   11.58  0.4057E+01  0.8690E+01
   25   18   12    1.98   58.29  0.2043E+02  0.1164E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT       101
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"        1468
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)    1259
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))      209


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.2 sec
