#include "Class_PPM_Image.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>

#include<chrono> / Timer 

using namespace std;
using namespace chrono;


//read the PPM image from fname
//write the PPM image in fname

//Split "mem" into "parts", e.g. if mem = 10 and parts = 4 you will have: 0,2,4,6,10
//if possible the function will split mem into equal chuncks, if not 
//the last chunck will be slightly larger

std::vector<int> bounds(int parts, int mem) {
    std::vector<int>bnd;
    int delta = mem / parts;
    int reminder = mem % parts;
    int N1 = 0, N2 = 0;
    bnd.push_back(N1);
    for (int i = 0; i < parts; ++i) {
        N2 = N1 + delta;
        if (i == parts - 1)
            N2 += reminder;
        bnd.push_back(N2);
        N1 = N2;
    }
    return bnd;
}

//Test if a given position (ii,jj) is "inside" the limits 0..nr_lines and 0..nr_columns

bool border(int ii, int jj, int nr_lines, int nr_columns) {
    if (ii >= 0 && ii < nr_lines && jj >= 0 && jj < nr_columns)
        return true;
    else
        return false;
}

// New Filter Function Here!!

// PPM Image Filter processing
void invertColors(ppm& image, ppm& image2, int i, int j) {
    int ii, jj, nr_lines, nr_columns, max_color , indx;
    unsigned int r, g, b;
   // float r_sum, g_sum, b_sum;

    nr_lines = image.height;
    nr_columns = image.width;
    max_color = image.max_col_val;

    //Apply the filter:
    r = 0;
    g = 0;
    b = 0;

    //check all  center  PIXEL Image
    ii = i ;
    jj = j;
    if (border(ii, jj, nr_lines, nr_columns)) {
        indx = ii * image.width + jj;

        r = (unsigned int)image.r[indx];
        g = (unsigned int)image.g[indx];
        b = (unsigned int)image.b[indx];

        r = max_color - r;
        g = max_color - g;
        b = max_color - b;

    }
    //Save the modifed pixel value in image2 
    indx = i * image.width + j;
    image2.r[indx] = (unsigned char)r;
    image2.g[indx] = (unsigned char)g;
    image2.b[indx] = (unsigned char)b;
 

};

//Blur the pixel at (i,j) using information from the neighbour pixels

void process(ppm& image, ppm& image2, int i, int j) {
    int ii, jj, nr_lines, nr_columns, indx;
    unsigned int r, g, b;
    float r_sum, g_sum, b_sum;
    //Filter used for bluring an image
    float filter[] = {
        0.10179640718562874, 0.11377245508982035, 0.10179640718562874,
        0.11377245508982035, 0.1377245508982036, 0.11377245508982035,
        0.10179640718562874, 0.11377245508982035, 0.10179640718562874
    };

    nr_lines = image.height;
    nr_columns = image.width;

    //Apply the filter:
    r_sum = 0;
    g_sum = 0;
    b_sum = 0;

    //check North-West
    ii = i - 1;
    jj = j - 1;
    if (border(ii, jj, nr_lines, nr_columns)) {
        indx = ii * image.width + jj;

        r = (unsigned int)image.r[indx];
        g = (unsigned int)image.g[indx];
        b = (unsigned int)image.b[indx];

        r_sum += r * filter[0];
        g_sum += g * filter[0];
        b_sum += b * filter[0];
    }

    //check North
    ii = i - 1;
    jj = j;
    if (border(ii, jj, nr_lines, nr_columns)) {
        indx = ii * image.width + jj;

        r = (unsigned int)image.r[indx];
        g = (unsigned int)image.g[indx];
        b = (unsigned int)image.b[indx];

        r_sum += r * filter[1];
        g_sum += g * filter[1];
        b_sum += b * filter[1];
    }

    //check North-East
    ii = i - 1;
    jj = j + 1;
    if (border(ii, jj, nr_lines, nr_columns)) {
        indx = ii * image.width + jj;

        r = (unsigned int)image.r[indx];
        g = (unsigned int)image.g[indx];
        b = (unsigned int)image.b[indx];

        r_sum += r * filter[2];
        g_sum += g * filter[2];
        b_sum += b * filter[2];
    }

    //check West
    ii = i;
    jj = j - 1;
    if (border(ii, jj, nr_lines, nr_columns)) {
        indx = ii * image.width + jj;

        r = (unsigned int)image.r[indx];
        g = (unsigned int)image.g[indx];
        b = (unsigned int)image.b[indx];

        r_sum += r * filter[3];
        g_sum += g * filter[3];
        b_sum += b * filter[3];
    }

    //center
    ii = i;
    jj = j;
    indx = ii * image.width + jj;

    r = (unsigned int)image.r[indx];
    g = (unsigned int)image.g[indx];
    b = (unsigned int)image.b[indx];

    r_sum += r * filter[4];
    g_sum += g * filter[4];
    b_sum += b * filter[4];


    //check East
    ii = i;
    jj = j + 1;
    if (border(ii, jj, nr_lines, nr_columns)) {
        indx = ii * image.width + jj;

        r = (unsigned int)image.r[indx];
        g = (unsigned int)image.g[indx];
        b = (unsigned int)image.b[indx];

        r_sum += r * filter[5];
        g_sum += g * filter[5];
        b_sum += b * filter[5];
    }
    //check South-West
    ii = i + 1;
    jj = j - 1;
    if (border(ii, jj, nr_lines, nr_columns)) {
        indx = ii * image.width + jj;

        r = (unsigned int)image.r[indx];
        g = (unsigned int)image.g[indx];
        b = (unsigned int)image.b[indx];

        r_sum += r * filter[6];
        g_sum += g * filter[6];
        b_sum += b * filter[6];
    }
    //check South
    ii = i + 1;
    jj = j;
    if (border(ii, jj, nr_lines, nr_columns)) {
        indx = ii * image.width + jj;

        r = (unsigned int)image.r[indx];
        g = (unsigned int)image.g[indx];
        b = (unsigned int)image.b[indx];

        r_sum += r * filter[7];
        g_sum += g * filter[7];
        b_sum += b * filter[7];
    }
    //check South-East
    ii = i + 1;
    jj = j + 1;
    if (border(ii, jj, nr_lines, nr_columns)) {
        indx = ii * image.width + jj;

        r = (unsigned int)image.r[indx];
        g = (unsigned int)image.g[indx];
        b = (unsigned int)image.b[indx];

        r_sum += r * filter[8];
        g_sum += g * filter[8];
        b_sum += b * filter[8];
    }

    //Save the modifed pixel value in image2 
    indx = i * image.width + j;
    image2.r[indx] = (unsigned char)r_sum;
    image2.g[indx] = (unsigned char)g_sum;
    image2.b[indx] = (unsigned char)b_sum;
}

//Blur a chunck of an image

void tst(ppm& image, ppm& image2, int left, int right) {
    for (int i = left; i < right; ++i) {
        int ii = i / image.width;
        int jj = i - ii * image.width;
        process(image, image2, ii, jj);
        // add new filter process Here!!
        invertColors(image2, image, ii, jj);
        process(image, image2, ii, jj);
    }
}


int main() {
    //*******************************************************
    // Imput images: fname
    // std::string fname = std::string("your_file_name.ppm");
    std::string fname = std::string("Images\\sample_5184×3456.ppm");
   // std::string fname = std::string("Images/sample_1280×853.ppm");
    //std::string fname = std::string("Images\\sample_640×426.ppm");
    //std::string fname = std::string("Images\\Image.ppm");
    //*********************************************************
    // Output Image fnameout
    // std::string fnameout = std::string("your_file_name.ppm");
    //std::string fnameout = std::string("Images\\Output\\INVERTCOLOROutput_sample_640×426.ppm");
   std::string fnameout = std::string("Images/Output/INVERTEO11utput_sample_1280×853.ppm");
   // std::string fnameout = std::string("Images/Output/OutputinvertInvertcolor_sample_5184×3456.ppm");
     // std::string fnameout = std::string("Images/Output/O2utputImage.ppm");

    ppm image(fname);
    ppm image2(image.width, image.height);

    //Number of threads to use (the image will be divided between threads)
    const int num_threads = thread::hardware_concurrency();
    int parts = num_threads;
    cout << "num_threads: " << num_threads << endl;

    std::vector<int>bnd = bounds(parts, image.size);

    std::vector<std::thread> tt(parts - 1);
   // std::vector<std::thread> tt;

    // Timer starts now
    auto start1 = high_resolution_clock::now();

    //Lauch parts-1 threads
    for (int i = 0; i < parts - 1; ++i) {
        //tt.push_back(std::thread(tst, std::ref(image), std::ref(image2), bnd[i], bnd[i + 1]));
        tt[i] = std::thread(tst, std::ref(image), std::ref(image2), bnd[i], bnd[i + 1]);
    }

    //Use the main thread to do part of the work !!!
    for (int i = parts - 1; i < parts; ++i) {
        tst(image, image2, bnd[i], bnd[i + 1]);
    }

    //Join parts-1 threads
 /*   for (auto& e : tt) {
        e.join();
    }*/
    for (int i = 0; i < parts - 1; ++i) {
        tt[i].join();
    }
    // Timer ends now
    auto stop1 = high_resolution_clock::now();

    // Difference is calculated in microseconds
    auto durationmicro = duration_cast<microseconds>(stop1 - start1);

    // Difference is calculated in nanoseconds
    auto durationnano = duration_cast<nanoseconds>(stop1 - start1);

    // Displays time in output file in microseconds
    cerr << "Time taken in microseconds : "<< (double)(durationmicro.count() / 1000.0) << " microseconds" << endl;
    // Displays time in output file in nanoseconds
    cout << "Time taken in nanoseconds : " << (double)(durationnano.count() / 1000.0) << " nanoseconds" << endl;


   // std::cout << "Time Diff in Total Micro Seconds " << duration.total_microseconds() << std::endl;

    //Save the result
    // Output Image fnameout
    image2.write(fnameout);

    return 0;
}