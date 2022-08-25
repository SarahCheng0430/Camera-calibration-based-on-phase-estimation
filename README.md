# Camera-calibration-based-on-phase-estimation
+Camera_calibration_toolbox for test.

1. Verify that MATLAB is installed. 

2. Open 'Camera_calibration_toolbox\ app.mlapp'. 
    (Open MATLAB--> Open the Folder 'Camera_calibration_toolbox' --> Open ' app.mlapp',  and run it. )

3. In the GUI, click 'open files', select all the images from the folder 'samples_for_test'.

4. Fill in Three parameters as follows, which are parameters of the captured samples in  'samples_for_test'.
                    T(Pattern Pitch)--4000;  
                    t(Pixel Size)--10; 
                    f(Norminal Focal length)--9000

5. Click 'Calibration Starts !'  then the Intrinsic parameters and radial distortion coefficients would be displayed.

6. Click 'Export Data' to export all the intrinsic and extrinsic parameters.

7. Click 'Calculate Reprojection Error',
    then reprojection error of  each image would be displayed in a Bar Figure 
    and root mean square errors displayed in a Scatter Figure.
