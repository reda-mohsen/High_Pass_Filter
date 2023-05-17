#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include "omp.h"
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#pragma once

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;


// function to take in the path of an image file, reads the image, and converts it to a grayscale image stored as an array of integers.
int* inputImage(int* w, int* h, System::String^ imagePath, int filterSize) //put the size of image in w & h
{
	int* input;
	int OriginalImageWidth, OriginalImageHeight;
	//Read Image and save it to local arrayss
	System::Drawing::Bitmap BM(imagePath); // create a Bitmap object from the image file

	OriginalImageWidth = BM.Width; // set the width of the image
	OriginalImageHeight = BM.Height;  // set the height of the image

	// Calculate padded width and height
	*w = BM.Width + 2 * (filterSize / 2);
	*h = BM.Height + 2 * (filterSize / 2);

	int* Red = new int[*w * *h]; // create an array for the red channel
	int* Green = new int[*w * *h];  // create an array for the green channel
	int* Blue = new int[*w * *h];  // create an array for the blue channel
	input = new int[*w * *h]; // create an array for the grayscale image

	// Calculate padding offset
	int paddingOffset = filterSize / 2;

	for (int i = 0; i < *h; i++) {  // loop through the rows of the image
		for (int j = 0; j < *w; j++) { // loop through the columns of the image
			int originalRow = i - paddingOffset;
			int originalCol = j - paddingOffset;

			// Check if the current pixel is within the original image boundaries
			if (originalRow >= 0 && originalRow < OriginalImageHeight && originalCol >= 0 && originalCol < OriginalImageWidth) {
				System::Drawing::Color c = BM.GetPixel(originalCol, originalRow); // get the color of the pixel at (j,i)

				Red[i * BM.Width + j] = c.R; // store the red channel value
				Blue[i * BM.Width + j] = c.B; // store the blue channel value
				Green[i * BM.Width + j] = c.G; // store the green channel value

				input[i * *w + j] = ((c.R + c.B + c.G) / 3); // gray scale value equals the average of RGB values
			}
			else {
				// Set the values to 0 for the padded pixels
				Red[i * *w + j] = 0;
				Blue[i * *w + j] = 0;
				Green[i * *w + j] = 0;
				input[i * *w + j] = 0;
			}
		}
	}

	return input; // return the grayscale image array
}


// This function takes an input image array, its width and height, and the desired filter size as arguments
// and applies a parallel high pass filter to the input image array to obtain an output image array
int* openMp_highPassFilter(int* input, int width, int height, int filterSize) {

	// Serial Section
	// Allocate memory dynamically for output image array to ensure enough memory 
	int* output = new int[width * height];
	// Set the filter kernel dynamically based on the filter size
	std::vector<int> filter(filterSize * filterSize);
	int filterOffset = filterSize / 2;

	// Parallel Section
#pragma omp parallel for collapse(2) schedule(dynamic) shared(filterOffset, filter) num_threads(8)
	for (int j = 0; j < filterSize; j++) {
		for (int i = 0; i < filterSize; i++) {
			if (i == filterOffset && j == filterOffset) {
				// Set the center pixel of the filter kernel to the sum of all other pixels
				filter[j * filterSize + i] = (filterSize * filterSize - 1);
			}
			else {
				// Set all other pixels of the filter kernel to -1
				filter[j * filterSize + i] = -1;
			}
		}
	}

	// Apply the high pass filter to each pixel in the input image array reduction(+:sum)
	int sum = 0;
	// Parallel Section
#pragma omp parallel for collapse(2) schedule(dynamic) shared(input,output,filter,filterOffset) reduction(+:sum) num_threads(8)

	for (int y = filterOffset; y < height - filterOffset; y++) {
		for (int x = filterOffset; x < width - filterOffset; x++) {
			sum = 0;
			for (int j = -filterOffset; j <= filterOffset; j++) {
				for (int i = -filterOffset; i <= filterOffset; i++) {
					// The filtered value of each pixel is stored in the corresponding position in the output image array
					sum += input[(y + j) * width + (x + i)] * filter[(j + filterOffset) * filterSize + (i + filterOffset)];
				}
			}
			// Set the output pixel value
			//output[y * width + x] = sum;
			//output[y * width + x] = 128 - sum;		//  Scale the output pixel for display purpose (visualization)
			output[y * width + x] = sum + input[y * width + x];		// Add the filtered image array to the input image array to get final sharp image
		}
	}

	// Return the filtered image array
	return output;
}

void createImage(int* image, int width, int height, int index)
{
	// create a new bitmap with the specified width and height
	System::Drawing::Bitmap MyNewImage(width, height);

	// iterate over each pixel in the image
	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			// ensure that the pixel value is within the range of 0-255
			//i * OriginalImageWidth + j
			if (image[i * width + j] < 0)
			{
				image[i * width + j] = 0;
			}
			if (image[i * width + j] > 255)
			{
				image[i * width + j] = 255;
			}

			// create a new color using the pixel value and set it as the pixel color in the new bitmap
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}

	// save the new bitmap to a file with a specified name (using the index parameter)
	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");

	// print a message indicating that the image has been saved
	cout << "result Image Saved " << index << endl;
}




int main()
{
	int ImageWidth = 4, ImageHeight = 4, FilterSize = 3; // Initialize the dimensions of the input image and the filter size(odd value).

	do {
		cout << "Enter the filter size you want (should be odd value): ";
		cin >> FilterSize;
	} while (FilterSize % 2 != 1 || FilterSize < 1);

	cout << endl << "Processing...." << endl;

	System::String^ imagePath; // Declare the path of the input image as a System::String^ variable
	std::string img = "..//Data//Input//dog.png"; // Initialize the path of the input image as a std::string variable

	imagePath = marshal_as<System::String^>(img); // Convert the std::string variable to a System::String^ variable using the marshal_as function

	int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath, FilterSize); // Call the inputImage function to retrieve the grayscale values of the pixels of the input image and store them in the integer array imageData

	int* filteredImage = new int[ImageWidth * ImageHeight];

	int start_s, stop_s, openmp_TotalTime = 0; // Declare variables to measure the time taken to process the image

	start_s = clock(); // Start the timer
	// Perform any image processing tasks or calculations
	// Apply the high pass filter using openmp
	filteredImage = openMp_highPassFilter(imageData, ImageWidth, ImageHeight, FilterSize);

	stop_s = clock(); // Stop the timer
	openmp_TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000; // Calculate the total time taken to process the image

	createImage(filteredImage, ImageWidth, ImageHeight, 9); // Call the createImage function to save the resulting image

	cout << "Time using OpenMp: " << openmp_TotalTime << " ms" << endl; // Print the total time taken using openmp/parallel

	delete[] filteredImage;

	free(imageData); // Free the memory allocated for the imageData array

	system("pause");
	return 0;
}