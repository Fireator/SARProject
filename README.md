# SARProject

Contains the CMake files used to generate sample image of the altitudes calculated from the point cloud values of the ranges taken from txt files.


Things completed:

- Import data from file
- Calculate altitudes
- Calculate heading
- Coordinate mapping
- Generate test file in bmp format for viewing
- Resize using bilinear interpolation to 480 x 480


  Unfinished:

- Mirror final image
- Use the coordinate mapping to generate an image (using bilinear interpolation for missing data)
- Or use OpenCV and apply inverse perspective warp to skew 2 corners of the image + Crop for final square image




