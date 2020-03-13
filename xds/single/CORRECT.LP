
 ***** CORRECT ***** (VERSION Jan 31, 2020  BUILT=20200131)  29-Feb-2020

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=    4 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    27.985    62.713    59.812  90.000  90.923  90.000 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=2000000000     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    50.000     1.000
 NAME_TEMPLATE_OF_DATA_FRAMES=/Users/atakisi/MLFSOM/xds/single/lm_h97_c01_rt_10_???.cbf CBF
 DATA_RANGE=       1       1
 ROTATION_AXIS=  1.000000  0.000000  0.000000
 OSCILLATION_RANGE=  1.000000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       1
 X-RAY_WAVELENGTH=  0.630900
 INCIDENT_BEAM_DIRECTION=  0.000974  0.001408  1.585036
 FRACTION_OF_POLARIZATION=   0.900
 POLARIZATION_PLANE_NORMAL=  0.000000  1.000000  0.000000
 AIR=   0.00008
 DETECTOR=PILATUS
 SILICON=  0.981077 SENSOR_THICKNESS=  0.320000
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  2463  NY=  2527    QX=  0.172000  QY=  0.172000
 ORGX=   1229.00  ORGY=   1207.00
 DETECTOR_DISTANCE=   400.000
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 BEAM_DIVERGENCE_E.S.D.=   0.031
 REFLECTING_RANGE_E.S.D.=   0.414
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

     170 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       1 ...       1
       0 OVERLOADED REFLECTIONS (IGNORED)
       0 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
     170 REFLECTIONS ACCEPTED



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
 SPACE_GROUP_NUMBER=    4 (0: unknown symmetry and unit cell)
 UNIT_CELL_CONSTANTS=    27.00    63.00    60.00  90.000  90.500  90.000
 UNIT_CELL_A-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_B-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_C-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 REIDX=  0  0  0  0  0  0  0  0  0  0  0  0 (all 0 : not specified)
 TEST_RESOLUTION_RANGE=    8.00    4.50 (Angstrom) for space group determination
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas                0
 NUMBER OF ACCEPTED UNIQUE REFLECTIONS FROM REFERENCE DATA SET       0


 !!! WARNING !!! NO REFLECTIONS TO COMPUTE Rmeas.
                 INCORRECT CHOICE OF TEST_RESOLUTION_RANGE= ?


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

 *  44        aP          0.0      28.0   59.8   62.7  90.0  90.0  90.9   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  31        aP          0.0      28.0   59.8   62.7  90.0  90.0  89.1    1  0  0  0  0  0 -1  0  0  1  0  0
 *  34        mP          0.0      28.0   62.7   59.8  90.0  90.9  90.0    1  0  0  0  0  1  0  0  0  0  1  0
 *  35        mP         16.0      59.8   28.0   62.7  90.0  90.0  90.9    0  0  1  0  1  0  0  0  0  1  0  0
 *  33        mP         16.0      28.0   59.8   62.7  90.0  90.0  90.9   -1  0  0  0  0  0 -1  0  0 -1  0  0
 *  32        oP         16.0      28.0   59.8   62.7  90.0  90.0  90.9   -1  0  0  0  0  0 -1  0  0 -1  0  0
    39        mC        234.0     122.4   28.0   62.7  90.0  90.0  77.6    1  0  2  0  1  0  0  0  0  1  0  0
    38        oC        234.0      28.0  122.4   62.7  90.0  90.0 102.4    1  0  0  0 -1  0 -2  0  0  1  0  0
    29        mC        234.0      28.0  122.4   62.7  90.0  90.0  77.6    1  0  0  0  1  0  2  0  0 -1  0  0
    25        mC        246.2      86.6   86.6   28.0  89.4  90.6  87.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    20        mC        246.2      86.6   86.6   28.0  90.6  90.6  92.8    0  1  1  0  0 -1  1  0  1  0  0  0
    23        oC        246.2      86.6   86.6   28.0  89.4  90.6  87.2    0 -1 -1  0  0 -1  1  0 -1  0  0  0
    21        tP        246.2      59.8   62.7   28.0  90.0  90.9  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    28        mC        258.0      28.0  128.5   59.8  90.2  90.9  77.4   -1  0  0  0 -1  2  0  0  0  0 -1  0
    37        mC        266.0     128.5   28.0   59.8  90.9  90.2  77.4   -1 -2  0  0 -1  0  0  0  0  0 -1  0
    36        oC        266.0      28.0  128.5   59.8  89.8  90.9 102.6   -1  0  0  0  1  2  0  0  0  0 -1  0
    27        mC        500.0     122.4   28.0   86.6  89.4 132.2  77.6   -1  0 -2  0 -1  0  0  0  0  1  1  0
    26        oF        625.0      28.0  122.4  128.5  87.3 102.6 102.4    1  0  0  0 -1  0 -2  0 -1  2  0  0
    19        oI        730.2      28.0   86.6   90.8  87.2  72.6  89.4   -1  0  0  0  0  1  1  0 -1  1 -1  0
    18        tI        855.2      86.6   90.8   28.0  72.6  90.6  92.8    0  1  1  0  1 -1  1  0  1  0  0  0
     1        cF        999.0      90.8   91.3   90.8  97.9  92.6 144.2    1  1  1  0  1 -1 -1  0 -1  1 -1  0
     2        hR        999.0      65.6   68.7   91.3 119.9  62.2  99.7   -1  0 -1  0  1  1  0  0  1 -1 -1  0
     3        cP        999.0      28.0   59.8   62.7  90.0  90.0  90.9   -1  0  0  0  0  0 -1  0  0 -1  0  0
     5        cI        999.0      68.7   65.6   86.6  51.4  48.9  80.3   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     4        hR        999.0      66.4   68.7   90.8 120.6  62.4 100.2   -1  0  1  0  1 -1  0  0  1  1  1  0
     6        tI        999.0      86.6   68.7   65.6  80.3  51.4  48.9    0 -1 -1  0 -1 -1  0  0 -1  0 -1  0
     7        tI        999.0      68.7   65.6   86.6  51.4  48.9  80.3   -1 -1  0  0 -1  0 -1  0  0 -1 -1  0
     8        oI        999.0      65.6   68.7   86.6  48.9  51.4  80.3    1  0  1  0  1  1  0  0  0  1  1  0
     9        hR        999.0      28.0   65.6  199.5 102.3  98.3 114.4    1  0  0  0 -1  0 -1  0 -1  3  1  0
    10        mC        999.0      66.4   65.6   62.7  90.0  90.0 129.8   -1  0  1  0 -1  0 -1  0  0 -1  0  0
    11        tP        999.0      28.0   59.8   62.7  90.0  90.0  90.9   -1  0  0  0  0  0 -1  0  0 -1  0  0
    12        hP        999.0      28.0   59.8   62.7  90.0  90.0  90.9   -1  0  0  0  0  0 -1  0  0 -1  0  0
    13        oC        999.0      65.6   66.4   62.7  90.0  90.0  50.2   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    15        tI        999.0      28.0   59.8  141.5  65.2  79.0  90.9   -1  0  0  0  0  0 -1  0 -1 -2 -1  0
    16        oF        999.0      66.4   65.6  141.5 117.6 107.3  50.2    1  0 -1  0 -1  0 -1  0  1  2  1  0
    14        mC        999.0      65.6   66.4   62.7  90.0  90.0  50.2   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    17        mC        999.0      65.6   66.4   62.7  90.0  90.0  50.2   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    22        hP        999.0      59.8   62.7   28.0  90.0  90.9  90.0    0  0 -1  0  0 -1  0  0 -1  0  0  0
    24        hR        999.0     137.5   86.6   28.0  89.4  79.0 105.5   -1 -1 -2  0  0 -1  1  0 -1  0  0  0
    30        mC        999.0      59.8  138.9   28.0  90.4  90.9  64.5    0  0 -1  0  0 -2 -1  0 -1  0  0  0
    40        oC        999.0      59.8  138.9   28.0  89.6  90.9 115.5    0  0  1  0  0 -2 -1  0  1  0  0  0
    42        oI        999.0      28.0   59.8  141.5 114.8 101.0  90.9    1  0  0  0  0  0  1  0 -1 -2 -1  0
    41        mC        999.0     138.9   59.8   28.0  90.9  90.4  64.5    0  2  1  0  0  0  1  0  1  0  0  0
    43        mI        999.0      66.4  141.5   59.8 114.8 155.1  72.7   -1  0  1  0  1  2  1  0  0  0 -1  0

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

  *   0.00        0    0.0       0    0.03  1  0  0  0  0  1  0  0  0  0  1  0
      0.00        0    0.0       0    0.04 -1  0  0  0  0 -1  0  0  0  0  1  0
      0.00        0    0.0       0    0.05  1  0  0  0  0  0 -1  0  0  1  0  0
      0.00        0    0.0       0    0.05 -1  0  0  0  0  0  1  0  0  1  0  0

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=    4
 UNIT_CELL_CONSTANTS=    27.00    63.00    60.00  90.000  90.500  90.000
 UNIT_CELL_A-AXIS=    -0.620    21.738   -17.597
 UNIT_CELL_B-AXIS=   -62.668     0.115     2.351
 UNIT_CELL_C-AXIS=     1.829    36.920    46.969
 REFLECTIONS H,K,L FROM "INTEGRATE.HKL" WILL BE REINDEXED BY
 H' = (  1*H +  0*K +  0*L)/ 1 +  0
 K' = (  0*H +  1*K +  0*L)/ 1 +  0
 L' = (  0*H +  0*K +  1*L)/ 1 +  0

 *****  2 EQUIVALENT POSITIONS IN SPACE GROUP #  4 *****

    If x',y',z' is an equivalent position to x,y,z, then
        x'=x*ML(1)+y*ML( 2)+z*ML( 3)+ML( 4)/12.0
        y'=x*ML(5)+y*ML( 6)+z*ML( 7)+ML( 8)/12.0
        z'=x*ML(9)+y*ML(10)+z*ML(11)+ML(12)/12.0

    #    1  2  3  4    5  6  7  8    9 10 11 12 
    1    1  0  0  0    0  1  0  0    0  0  1  0
    2   -1  0  0  0    0  1  0  6    0  0 -1  0

 ******************************************************************************
      MEAN DISCREPANCIES BETWEEN OBSERVED AND CALCULATED SPOT LOCATIONS
 ******************************************************************************

 The discrepancies in X- and Y-coordinates of the spots are depicted in the
 two images DX-CORRECTIONS.cbf and DY-CORRECTIONS.cbf for inspection with
 the XDS-Viewer.



 ******************************************************************************
  REFINEMENT OF DIFFRACTION PARAMETERS USING ALL IMAGES
 ******************************************************************************


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM        36 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM  ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.42
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.09
 SPACE GROUP NUMBER      4
 UNIT CELL PARAMETERS     27.855    62.569    59.404  90.000  90.504  90.000
 E.S.D. OF CELL PARAMETERS -1.0E+00-1.0E+00-1.0E+00-1.0E+00-1.0E+00-1.0E+00
 REC. CELL PARAMETERS   0.035901  0.015982  0.016834  90.000  89.496  90.000
 COORDINATES OF UNIT CELL A-AXIS    -0.566    21.650   -17.517
 COORDINATES OF UNIT CELL B-AXIS   -62.530     0.133     2.185
 COORDINATES OF UNIT CELL C-AXIS     1.703    36.969    46.468
 CRYSTAL MOSAICITY (DEGREES)     0.414
 LAB COORDINATES OF ROTATION AXIS  1.000000  0.000000  0.000000
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)  -0.018699  0.014380  1.584862
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    1231.34   1208.46
 DETECTOR ORIGIN (PIXELS) AT                     1258.73   1187.39
 CRYSTAL TO DETECTOR DISTANCE (mm)       399.31
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=    4
 UNIT_CELL_CONSTANTS=    27.86    62.57    59.40  90.000  90.504  90.000


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
  0.15  0.25         5          0.199     242.326
  0.25  0.35        18         39.160       1.229
  0.35  0.45        40         45.697       1.053
  0.45  0.55        51         69.627       0.691
  0.55  0.65        40         45.123       1.066
  0.65  0.75        11         22.304       2.158
  0.75  0.85         5          9.143       5.263
 
 -0.05  0.95       170         48.122       1.000
 

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
 4.000E+00  1.000E-04   50.00


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

    3.371   0.959     2.28   0.17     12.47     25.70      12     170       0
   --------------------------------------------------------------------------
    3.371   1.134     2.28   0.17     12.47     25.70      12     170       0
 


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
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       1       1)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED


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

     2.46           6       6      7451        0.1%     -99.9%    -99.9%        0    0.00    -99.9%     0.0      0    0.000       0
     2.04          28      27      5719        0.5%      15.0%     16.6%        2    4.96     21.3%     0.0      0    0.000       0
     1.77          26      25      6678        0.4%       8.3%     10.3%        2    3.21     11.7%     0.0      0    0.000       0
     1.59          25      23      7533        0.3%      14.3%     84.1%        4    1.08     20.3%     0.0      0    0.000       0
     1.46          25      24      8277        0.3%      88.8%    298.7%        2    0.49    125.6%     0.0      0    0.000       0
     1.35          20      20      9028        0.2%     -99.9%    -99.9%        0    0.00    -99.9%     0.0      0    0.000       0
     1.27          26      25      9591        0.3%     -99.9%    -99.9%        2    0.00    -99.9%     0.0      0    0.000       0
     1.19          11      11     10284        0.1%     -99.9%    -99.9%        0    0.00    -99.9%     0.0      0    0.000       0
     1.13           3       3     10815        0.0%     -99.9%    -99.9%        0    0.00    -99.9%     0.0      0    0.000       0
    total         170     164     75376        0.2%      12.5%     25.7%       12    2.33     17.6%    97.5*     0    0.000       0


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES       170
 NUMBER OF REJECTED MISFITS                                 0
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                          170
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                    164



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

   WILSON LINE (using all data) : A=   6.517 B=  19.511 CORRELATION=  0.89
      #      RES      SS        <I>       log(<I>)       BO
       3     3.229   0.024  5.4001E+02       6.292       4.7
       3     2.643   0.036  7.7389E+02       6.651      -1.9
      14     2.387   0.044  1.0145E+02       4.620      21.6
      10     2.117   0.056  5.7785E+01       4.057      22.0
      14     2.008   0.062  8.2672E+01       4.415      16.9
      14     1.858   0.072  2.9090E+01       3.370      21.7
       3     1.743   0.082  3.0454E+01       3.416      18.8
      16     1.652   0.092  1.1774E+01       2.466      22.1
      11     1.588   0.099  8.2897E+00       2.115      22.2
       7     1.515   0.109  5.7556E+00       1.750      21.9
      14     1.464   0.117  2.6022E+00       0.956      23.8
      12     1.413   0.125  3.6917E+00       1.306      20.8
       6     1.360   0.135  1.9207E+00       0.653      21.7
      17     1.325   0.142  6.7069E-01      -0.399      24.3
       7     1.275   0.154 -9.7052E-01    -999.999   -1000.0
       9     1.238   0.163  5.4915E-01      -0.599      21.8
       1     1.195   0.175  4.4369E+00       1.490      14.3
       2     1.179   0.180  3.9650E+00       1.378      14.3
       1     1.134   0.194 -4.3002E+00    -999.999   -1000.0


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

       1     3.229        0.302        0.057        0.008
       1     2.387        0.389        0.084        0.013
       1     2.117        0.322        0.063        0.009
       1     2.008        3.358        2.132        0.967
       1     1.858        0.081        0.008        0.001
       1     1.652        0.039        0.003        0.000
       1     1.464        1.173        0.440        0.118
       8   overall        1.087        0.120        0.238


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

       2     3.229        0.576        0.232        0.074
       3     2.643        0.700        0.394        0.179
      13     2.387        1.096        1.208        1.189
       9     2.117        0.898        0.661        0.396
      13     2.008        0.569        0.315        0.147
      13     1.858        1.290        1.414        1.310
       3     1.743        0.598        0.270        0.099
      15     1.652        1.050        0.809        0.510
      11     1.588        1.295        1.735        1.973
       7     1.515        1.982        2.668        3.034
      13     1.464        2.488        2.322        2.689
      12     1.413        2.161        2.415        2.416
       6     1.360        3.678        1.276        4.732
      16     1.325       17.631       34.593      148.565
       9     1.238       65.377       55.878     1182.284
       1     1.195        0.500        0.167        0.042
       2     1.179        1.240        0.907        0.503
     156   overall        6.726        7.831       84.619

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       1     3.229  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  1.47
       1     2.387  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     2.117  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  1.47
       1     2.008  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       1     1.858  0.00  0.00  0.00  0.00  1.92  1.78  1.67  1.59  1.52  1.47
       1     1.652  0.00  0.00  0.00  2.11  1.92  1.78  1.67  1.59  1.52  1.47
       1     1.464  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
       8   overall  0.00  0.00  0.00  0.26  0.48  0.45  0.42  0.40  0.38  0.73


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

       2     3.229  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.91  0.84  0.79
       3     2.643  0.00  0.00  0.00  0.00  0.85  0.74  1.32  1.21  1.12  1.05
      13     2.387  0.81  1.27  1.19  0.93  1.17  1.19  1.07  0.98  0.91  0.85
       9     2.117  1.17  0.61  1.29  1.01  1.13  0.99  1.10  1.01  1.12  1.05
      13     2.008  0.00  2.12  1.48  1.17  0.98  0.85  1.07  0.98  0.91  0.85
      13     1.858  1.62  1.70  1.19  1.17  1.37  1.19  1.22  1.12  1.04  1.10
       3     1.743  0.00  0.00  0.00  0.00  0.00  0.74  0.66  0.61  1.12  1.05
      15     1.652  2.10  1.47  1.03  1.01  1.02  0.89  0.93  0.85  1.01  0.95
      11     1.588  0.96  0.50  1.05  1.10  0.92  0.81  1.08  1.16  1.38  1.29
       7     1.515  3.00  2.36  1.65  1.73  1.45  1.27  1.14  1.04  1.20  1.13
      13     1.464  3.23  1.70  1.48  1.40  1.17  1.02  1.07  0.98  0.91  0.85
      12     1.413  4.38  2.30  1.61  1.52  1.27  1.11  1.16  1.06  0.98  0.92
       6     1.360  1.75  0.92  0.64  0.51  0.42  0.37  0.66  0.61  0.56  0.53
      16     1.325  5.91  3.10  2.17  1.71  1.43  1.25  1.12  1.02  0.95  0.89
       9     1.238  4.67  2.45  1.71  1.35  1.13  0.99  0.88  0.81  0.75  0.70
       1     1.195  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  1.58
       2     1.179  5.25  2.76  1.93  1.52  1.27  1.11  0.99  0.91  0.84  0.79
     156   overall  2.29  1.59  1.26  1.11  1.04  0.94  0.98  0.92  0.93  0.89


 List of      3 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

   20   14   18    1.23   18.13  0.9954E+01  0.1235E+02 "alien"
   20   12   17    1.25   14.82  0.8138E+01  0.1193E+02 "alien"
   12   18   32    1.33   13.94  0.9351E+01  0.1163E+02 "alien"


 List of      3 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   20   14   18    1.23   18.13  0.9954E+01  0.1235E+02
   20   12   17    1.25   14.82  0.8138E+01  0.1193E+02
   12   18   32    1.33   13.94  0.9351E+01  0.1163E+02

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT         0
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"         170
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES)     170
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))        0


 cpu time used                  0.1 sec
 elapsed wall-clock time        0.1 sec
