
 ***** DEFPIX ***** (VERSION Jan 31, 2020  BUILT=20200131)   6-Mar-2020

 modifies the initial background image BKGINIT.cbf from the INIT step by
 labelling detector pixels that should be ignored in the integration step
 and saves the results as a control image BKGPIX.cbf. Excluded pixels are
 marked by -3. The control images should be inspected with the XDS-Viewer.

 The pixels are labelled based on parameters in the input file XDS.INP:
 UNTRUSTED_ELLIPSE=, UNTRUSTED_RECTANGLE=, UNTRUSTED_QUADRILATERAL=,
 EXCLUDE_RESOLUTION_RANGE=, INCLUDE_RESOLUTION_RANGE=.

 In addition, DEFPIX attempts to recognize and exclude shaded regions
 of the detector area that are caused by intruding hardware obscuring
 x-ray diffraction. For this purpose each pixel value in the table is
 compared with the mean of all pixels at the same resolution; the result
 is saved in the control image ABS.cbf. Unshaded pixels should have
 values around 10000 while shaded regions are recognized by their lower
 values, typically below 7000. The acceptable value range for unshaded
 pixels is defined in the input file XDS.INP by the parameter
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=.

 UNTRUSTED_RECTANGLE=    487.   495.       1.  2527.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=    981.   989.       1.  2527.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=   1475.  1483.       1.  2527.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=   1969.  1977.       1.  2527.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.     195.   213.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.     407.   425.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.     619.   637.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.     831.   849.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.    1043.  1061.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.    1255.  1273.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.    1467.  1485.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.    1679.  1697.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.    1891.  1909.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.    2103.  2121.   (X1,X2, Y1,Y2)
 UNTRUSTED_RECTANGLE=      1.  2463.    2315.  2333.   (X1,X2, Y1,Y2)
 INCLUDE_RESOLUTION_RANGE=     50.000     1.100
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=      6000.    30000.
 RESOLUTION RANGE RECORDED BY DETECTOR (ANGSTROM)    49.996     1.290
 NUMBER OF PIXELS REMOVED FROM   "BKGINIT.cbf"         2544
 NUMBER OF TRUSTED PIXELS IN      "BKGPIX.cbf"      5224574
 MEAN VALUE OF TRUSTED PIXELS IN  "BKGPIX.cbf"     2801.303


 cpu time used                  0.6 sec
 elapsed wall-clock time        0.6 sec
