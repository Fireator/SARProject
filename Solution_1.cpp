//Solution_1.cpp : Defines the entry point for the application.
//

#include "Solution_1.h"
#include <vector>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cstring>
#include <sstream>
#include <cmath>



using namespace std;

//Structure for function arguments
struct SarMapArgs{

    size_t depthMapRowCount;
    size_t depthMapColCount;

    double depthMapMinAzD;
    double depthMapMinElD;
    double depthMapMaxAzD;
    double depthMapMaxElD;

    double srcLatD;
    double srcLonD;
    double srcAltM;

    double sarMapCorner1LatD;
    double sarMapCorner1LonD;
    double sarMapCorner2LatD;
    double sarMapCorner2LonD;
    double sarMapCorner3LatD;
    double sarMapCorner3LonD;
    double sarMapCorner4LatD;
    double sarMapCorner4LonD;



};

struct latLong {

    double y;
    double x;

};

struct Point {
    double x, y;
};

//Function prototypes
vector <uint8_t> generateSarMapTexture(const vector<float>& depthMapM, const SarMapArgs& args);
double altitudeCalc(double elevation, float range, double currentAltitude);
double slope(double x1, double y1, double x2, double y2);
vector <uint8_t> resize(int in_width, int in_height, int out_width, int out_height, const vector <uint8_t>& input);
void coordinateMapping(int x, int y, vector <latLong>& position, const SarMapArgs args);
double calculateLineLength(const latLong& p1, const latLong& p2);




//Test prototypes
void generateBMP(vector <uint8_t>& input, const char* filename, int w, int h);
vector <float> importList(const char* fileName);


int main()
{
    SarMapArgs info;

    //720p sample
    info.depthMapColCount = 1024;
    info.depthMapRowCount = 1024;


    info.depthMapMaxElD = -3.01042;
    info.depthMapMinElD = -4.41941;
    info.depthMapMaxAzD = 42.4402;
    info.depthMapMinAzD = 16.0127;

    info.srcAltM = 3048;
    info.srcLatD = 35.6328;
    info.srcLonD = -115;

    info.sarMapCorner1LatD = 36.0109;
    info.sarMapCorner1LonD = -114.866;

    info.sarMapCorner2LatD = 35.9229;
    info.sarMapCorner2LonD = -114.898;

    info.sarMapCorner3LatD = 35.825;
    info.sarMapCorner3LonD = -114.784;

    info.sarMapCorner4LatD = 35.9508;
    info.sarMapCorner4LonD = -114.642;
    
    


    //generate sample vector
    vector <float> sample = importList("DepthMapRangeValuesM1024.txt");    

    //run tests
    // 
    //build image
    vector <uint8_t> result = generateSarMapTexture(sample, info);
    //generate file
    generateBMP(result, "resize_sample.bmp", 480, 480);

	return 0;
}

//Returns 480x480 SAR map image generated
vector <uint8_t> generateSarMapTexture(const vector<float>& depthMapM, const SarMapArgs& args) {
    


    //Variables
    vector <uint8_t> output;
    vector <double> altitudes;
    vector <latLong> position;
    vector <latLong> cornerPts = { {args.sarMapCorner1LatD, args.sarMapCorner1LonD},{args.sarMapCorner2LatD, args.sarMapCorner2LonD},{args.sarMapCorner3LatD, args.sarMapCorner3LonD},{args.sarMapCorner4LatD, args.sarMapCorner4LonD} };
    int x = 0;
    int y = 0;

    size_t vectorMax = args.depthMapColCount * args.depthMapRowCount - 1;
    double maxAltitude = 0;
    double minAltitude = 0;

    // calculate elevation ratio
    double elevationRatio = (args.depthMapMaxElD - args.depthMapMinElD) / args.depthMapRowCount;

    //Calculate image bearings
    double heading = (args.depthMapMaxAzD - args.depthMapMinAzD) / 2 + args.depthMapMinAzD;
    //cout << heading;

    //find shortest side to calculate desired captured area
    vector<double> lineLengths;
    for (int i = 0; i < 4; i++) {
        int j = (i + 1) % 4;  //index of the next point

        latLong p1 = cornerPts[i];
        latLong p2 = cornerPts[j];
        
        lineLengths.push_back(calculateLineLength(p1, p2));
        
    }

    //start at the end of the array to invert image
    for (int i = vectorMax; i >= 0; i--) {

        //Add correct elevation ratio to current elevation for every row / maxElevation because of inversion
        double elevation = args.depthMapMaxElD - floor((vectorMax - i) / args.depthMapRowCount) * elevationRatio;


        //calculate altitude
        double altitude = altitudeCalc(elevation, depthMapM[i], args.srcAltM);
        
        //increment x and y coordinates
        if (x > 1023) {
            y++;
            x = 0;
        }
        
        //populate coordinate mapping
        coordinateMapping(x, y, position, args);

        //add the altitudes to the vector
        altitudes.push_back(altitude);
       
        //find max and min altitudes
        if (i == vectorMax) {
            minAltitude = altitude;
            maxAltitude = altitude;
        }

        //skewed data
        if (altitude > args.srcAltM || altitude < 0) {
            altitude = minAltitude;
        }
        
        if (altitude < minAltitude)
            minAltitude = altitude;

        if (altitude > maxAltitude)
            maxAltitude = altitude;
       
        //advance x
        x++;
        

    }

    //Calculate gradient scale for convering altitude to grayscale image
    double gradient = maxAltitude - minAltitude;
    double scale = 255.0 / gradient;
    //cout << endl << "||" << maxAltitude << endl;
    

    //define the stretching factors for each side
    double leftFactor = lineLengths[0] / lineLengths[1];
    double rightFactor = lineLengths[2] / lineLengths[1];
    double topFactor = lineLengths[3] / lineLengths[1];
    double bottomFactor = 1.0;

    
    //define the four points bounding the desired cropped region in the stretched image
    double startPoint = (1024.0 - (1024.0 * leftFactor)) / 2.0;
    Point topLeft = { startPoint, 1024 * leftFactor };
    Point topRight = { startPoint + 1024*leftFactor, 1024 * leftFactor };
    Point bottomLeft = { startPoint, 0};
    Point bottomRight = { startPoint + 1024 * leftFactor, 0 };




    //image generation
    for (int i = 0; i <= vectorMax; i++) {
        //calculate pixel value
        uint8_t pixel = (uint8_t)((altitudes[i] - minAltitude) * scale);
        
        //push the 4 bytes to the image buffer
        output.push_back(pixel);
        output.push_back(pixel);
        output.push_back(pixel);
        output.push_back(255);
      
    }

    //Perform prespective warping using 


    //480p output
    if (args.depthMapColCount == 480 && args.depthMapRowCount == 480) {
        return output; 
    }

    
    //Resize the current image to 480p
    else {
       return resize(args.depthMapRowCount, args.depthMapColCount, 480, 480, output);
    }

}

//Calculates altitude for the given point.
double altitudeCalc(double elevation, float range, double currentAltitude) {

    //Range elevation angorithm
    double rad_elevation = abs(elevation * 3.14159 / 180);
    double altitude = currentAltitude - (range * sin(rad_elevation));

    return altitude;

}

//Useful tools for mapping image vector to coordinates
void coordinateMapping(int x, int y, vector <latLong>& position, const SarMapArgs args)
{

    //for each point
    double x1, x2, x3, x4, y1, y2, y3, y4;

    y1 = args.sarMapCorner2LatD;
    y2 = args.sarMapCorner3LatD;
    y3 = args.sarMapCorner4LatD;
    y4 = args.sarMapCorner1LatD;

    x1 = args.sarMapCorner2LonD;
    x2 = args.sarMapCorner3LonD;
    x3 = args.sarMapCorner4LonD;
    x4 = args.sarMapCorner1LonD;

    double u = x / static_cast<double>((1024 - 1));
    double v = y / static_cast<double>((1024 - 1));

    double alpha = (1 - u) * (1 - v);
    double beta = u * (1 - v);
    double gamma = (1 - u) * v;
    double delta = u * v;

    double mappedX = alpha * x1 + beta * x2 + gamma * x3 + delta * x4;
    double mappedY = alpha * y1 + beta * y2 + gamma * y3 + delta * y4;

    //map points to the grid
    latLong curPosition = { mappedX, mappedY };
    position.push_back(curPosition);
}


double calculateLineLength(const latLong& p1, const latLong& p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    return sqrt(dx * dx + dy * dy);
}


//Bilinear Scaling Grayscale Algorithm adapted from tech-algorithm.com
vector <uint8_t> resize(int in_width, int in_height, int out_width, int out_height, const vector <uint8_t>& input) {
    vector <uint8_t> output;
    int A, B, C, D, x, y, index, gray;
    float x_ratio = ((float)(in_width - 1)) / out_width;
    float y_ratio = ((float)(in_height - 1)) / out_height;
    float x_diff, y_diff;
    
    for (int i = 0; i < out_height; i++) {
        for (int j = 0; j < out_width; j++) {

            //keep track of widths and columns
            x = (int)(x_ratio * j);
            y = (int)(y_ratio * i);
            x_diff = (x_ratio * j) - x;
            y_diff = (y_ratio * i) - y;
            index = y * in_width + x;

            // Get reference pixelf from main image
            A = input[index * 4] & 0xff;
            B = input[(index + 1) * 4] & 0xff;
            C = input[(index + in_width) * 4] & 0xff;
            D = input[(index + in_width + 1) * 4] & 0xff;

            // Bilinear scaling algorithm
            gray = (uint8_t)(
                A * (1 - x_diff) * (1 - y_diff) + B * (x_diff) * (1 - y_diff) +
                C * (y_diff) * (1 - x_diff) + D * (x_diff * y_diff)
                );

            //Push pixels into output image vector
            output.push_back(gray);
            output.push_back(gray);
            output.push_back(gray);
            output.push_back(255);

        }
    }

    //cout << output.size();
    return output;
}


//test generate image from output
void generateBMP(vector <uint8_t> &input, const char * filename, int w, int h) {
    FILE* f;
    unsigned char* img = NULL;
    int filesize = 54 + 3 * w * h;  //w is your image width, h is image height

    img = (unsigned char*)malloc(3 * w * h);
    memset(img, 0, 3 * w * h);

    for (int i = 0; i < w * h ; i++)
    {
        int pixel = input[i * 4];

        if (pixel < 0) {
            pixel = 0;
            
        }
        img[(i) * 3 + 2] = (unsigned char)(pixel);
        img[(i) * 3 + 1] = (unsigned char)(pixel);
        img[(i) * 3 + 0] = (unsigned char)(pixel);
       
    }

    unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
    unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
    unsigned char bmppad[3] = { 0,0,0 };

    bmpfileheader[2] = (unsigned char)(filesize);
    bmpfileheader[3] = (unsigned char)(filesize >> 8);
    bmpfileheader[4] = (unsigned char)(filesize >> 16);
    bmpfileheader[5] = (unsigned char)(filesize >> 24);

    bmpinfoheader[4] = (unsigned char)(w);
    bmpinfoheader[5] = (unsigned char)(w >> 8);
    bmpinfoheader[6] = (unsigned char)(w >> 16);
    bmpinfoheader[7] = (unsigned char)(w >> 24);
    bmpinfoheader[8] = (unsigned char)(h);
    bmpinfoheader[9] = (unsigned char)(h >> 8);
    bmpinfoheader[10] = (unsigned char)(h >> 16);
    bmpinfoheader[11] = (unsigned char)(h >> 24);

    f = fopen(filename, "wb");
    fwrite(bmpfileheader, 1, 14, f);
    fwrite(bmpinfoheader, 1, 40, f);
    for (int i = 0; i < h; i++)
    {
        fwrite(img + (w * (h - i - 1) * 3), 3, w, f);
        fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
    }

    free(img);
    fclose(f);
}


//For testing opening .txt files containing point coud range data
vector<float> importList(const char* fileName)
{
    vector<float> numbers;
    ifstream inputFile(fileName);

    if (inputFile.is_open()) {
        string line;
        while (std::getline(inputFile, line)) {
            istringstream iss(line);
            float num;
            while (iss >> num) {
                numbers.push_back(num);
            }
        }
        inputFile.close();
    }
    else {
        std::cout << "Unable to open the file." << std::endl;
    }

    return numbers;
}


