//#include <iostream>
//#include <math.h>
//#include <stdlib.h>
//#include<string.h>
//#include "mpi.h"
//#include<msclr\marshal_cppstd.h>
//#include <ctime>// include this header 
//#pragma once
//
//#using <mscorlib.dll>
//#using <System.dll>
//#using <System.Drawing.dll>
//#using <System.Windows.Forms.dll>
//using namespace std;
//using namespace msclr::interop;
//
//
//// function to take in the path of an image file, reads the image, and converts it to a grayscale image stored as an array of integers.
//int* inputImage(int* w, int* h, System::String^ imagePath, int filterSize) //put the size of image in w & h
//{
//	int* input;
//	int OriginalImageWidth, OriginalImageHeight;
//	//Read Image and save it to local arrayss
//	System::Drawing::Bitmap BM(imagePath); // create a Bitmap object from the image file
//
//	OriginalImageWidth = BM.Width; // set the width of the image
//	OriginalImageHeight = BM.Height;  // set the height of the image
//
//	// Calculate padded width and height
//	*w = BM.Width + 2 * (filterSize / 2);
//	*h = BM.Height + 2 * (filterSize / 2);
//
//	int* Red = new int[*w * *h]; // create an array for the red channel
//	int* Green = new int[*w * *h];  // create an array for the green channel
//	int* Blue = new int[*w * *h];  // create an array for the blue channel
//	input = new int[*w * *h]; // create an array for the grayscale image
//
//	// Calculate padding offset
//	int paddingOffset = filterSize / 2;
//
//	for (int i = 0; i < *h; i++) {  // loop through the rows of the image
//		for (int j = 0; j < *w; j++) { // loop through the columns of the image
//			int originalRow = i - paddingOffset;
//			int originalCol = j - paddingOffset;
//
//			// Check if the current pixel is within the original image boundaries
//			if (originalRow >= 0 && originalRow < OriginalImageHeight && originalCol >= 0 && originalCol < OriginalImageWidth) {
//				System::Drawing::Color c = BM.GetPixel(originalCol, originalRow); // get the color of the pixel at (j,i)
//
//				Red[i * BM.Width + j] = c.R; // store the red channel value
//				Blue[i * BM.Width + j] = c.B; // store the blue channel value
//				Green[i * BM.Width + j] = c.G; // store the green channel value
//
//				input[i * *w + j] = ((c.R + c.B + c.G) / 3); // gray scale value equals the average of RGB values
//			}
//			else {
//				// Set the values to 0 for the padded pixels
//				Red[i * *w + j] = 0;
//				Blue[i * *w + j] = 0;
//				Green[i * *w + j] = 0;
//				input[i * *w + j] = 0;
//			}
//		}
//	}
//
//	return input; // return the grayscale image array
//}
//
//
//void createImage(int* image, int width, int height, int index)
//{
//	// create a new bitmap with the specified width and height
//	System::Drawing::Bitmap MyNewImage(width, height);
//
//	// iterate over each pixel in the image
//	for (int i = 0; i < MyNewImage.Height; i++)
//	{
//		for (int j = 0; j < MyNewImage.Width; j++)
//		{
//			// ensure that the pixel value is within the range of 0-255
//			//i * OriginalImageWidth + j
//			if (image[i * width + j] < 0)
//			{
//				image[i * width + j] = 0;
//			}
//			if (image[i * width + j] > 255)
//			{
//				image[i * width + j] = 255;
//			}
//
//			// create a new color using the pixel value and set it as the pixel color in the new bitmap
//			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j], image[i * MyNewImage.Width + j]);
//			MyNewImage.SetPixel(j, i, c);
//		}
//	}
//
//	// save the new bitmap to a file with a specified name (using the index parameter)
//	MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
//
//	// print a message indicating that the image has been saved
//	cout << "result Image Saved " << index << endl;
//}
//
//
//int main()
//{
//	int width = 4, height = 4, filterSize = 3; // Initialize the dimensions of the input image and the filter size(odd value).
//
//	MPI_Init(NULL, NULL);  // Initialize MPI
//	int size, rank;
//	MPI_Comm_size(MPI_COMM_WORLD, &size);  // Get the number of processes
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get the rank of the current process
//
//	do {
//		if (rank == 0) {
//			cout << "Enter the filter size you want (should be an odd value): ";
//			cin >> filterSize;
//			cout << endl << "Processing...." << endl;
//		}
//
//		// Broadcast the filter size to all other processes
//		MPI_Bcast(&filterSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//		// Check if the filter size is valid
//		if (filterSize % 2 != 1 || filterSize < 1) {
//			if (rank == 0) {
//				cout << "Invalid filter size. Please enter an odd value." << endl;
//			}
//		}
//	} while (filterSize % 2 != 1 || filterSize < 1);
//
//
//	System::String^ imagePath; // Declare the path of the input image as a System::String^ variable
//	std::string img = "..//Data//Input//moon.png"; // Initialize the path of the input image as a std::string variable
//
//	imagePath = marshal_as<System::String^>(img); // Convert the std::string variable to a System::String^ variable using the marshal_as function
//
//	int* imageData = inputImage(&width, &height, imagePath, filterSize); // Call the inputImage function to retrieve the grayscale values of the pixels of the input image and store them in the integer array imageData
//
//	int* filteredImage = new int[width * height];
//
//	int start_s, stop_s, mpi_TotalTime = 0; // Declare variables to measure the time taken to process the image
//
//	start_s = clock(); // Start the timer
//	// Perform any image processing tasks or calculations
//	// Apply the high pass filter using mpi
//	int local_height = height / size; // Calculate the local height for each process
//	int local_input_size = (local_height + 2) * width;
//	int local_output_size = local_height * width;
//	int* local_input = new int[local_input_size];
//	int* local_output = new int[local_output_size];
//
//	// Scatter the input image array to all processes
//	MPI_Scatter(imageData, local_height * width, MPI_INT, local_input + width, local_height * width, MPI_INT, 0, MPI_COMM_WORLD);
//
//	// Set the filter kernel dynamically based on the filter size
//	std::vector<int> filter(filterSize * filterSize);
//	int filterOffset = filterSize / 2;
//	for (int j = 0; j < filterSize; j++) {
//		for (int i = 0; i < filterSize; i++) {
//			if (i == filterOffset && j == filterOffset) {
//				// Set the center pixel of the filter kernel to the sum of all other pixels
//				filter[j * filterSize + i] = (filterSize * filterSize - 1);
//			}
//			else {
//				// Set all other pixels of the filter kernel to -1
//				filter[j * filterSize + i] = -1;
//			}
//		}
//	}
//
//	// Add border padding to local_input
//	if (rank == 0) {
//		// Copy the top border pixels to the extra top row in the local_input array
//		memcpy(local_input, imageData, width * sizeof(int));
//		// Copy the bottom border pixels to the extra bottom row in the local_input array
//		memcpy(local_input + (local_height + 1) * width, imageData + (height - 1) * width, width * sizeof(int));
//	}
//
//	// Copy the top border pixels to the extra top row in the local_input array of each process
//	if (rank > 0) {
//		memcpy(local_input, local_input + width, width * sizeof(int));
//	}
//
//	// Copy the bottom border pixels to the extra bottom row in the local_input array of each process
//	if (rank < size - 1) {
//		memcpy(local_input + (local_height + 1) * width, local_input + local_height * width, width * sizeof(int));
//	}
//
//	// Apply the high pass filter to each pixel in the input image array
//	for (int y = 1; y <= local_height; y++) {
//		for (int x = filterOffset; x < width - filterOffset; x++) {
//			int sum = 0;
//			for (int j = -filterOffset; j <= filterOffset; j++) {
//				for (int i = -filterOffset; i <= filterOffset; i++) {
//					// The filtered value of each pixel is stored in the corresponding position in the output image array
//					sum += local_input[((y + j) * width) + (x + i)] * filter[(j + filterOffset) * filterSize + (i + filterOffset)];
//				}
//			}
//			// Set the output pixel value
//			//local_output[(y - 1) * width + x] = sum;
//			//local_output[(y - 1) * width + x] = 128-sum;
//			local_output[(y - 1) * width + x] = sum + local_input[(y - 1) * width + x];
//		}
//	}
//
//	// Gather the filtered output image arrays from all processes
//	MPI_Gather(local_output, local_output_size, MPI_INT, filteredImage, local_output_size, MPI_INT, 0, MPI_COMM_WORLD);
//
//	if (rank == 0) {
//		stop_s = clock(); // Stop the timer
//		mpi_TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000; // Calculate the total time taken to process the image
//		createImage(filteredImage, width, height, 3); // Call the createImage function to save the resulting image
//		cout << "Time using MPI: " << mpi_TotalTime << " ms" << endl; // Print the total time taken using openmp/parallel
//		delete[] filteredImage;
//		MPI_Finalize();  // Finalize MPI
//	}
//	else {
//		MPI_Finalize();  // Finalize MPI
//	}
//
//	delete[] imageData; // Free the memory allocated for the input image array
//
//	return 0;
//}