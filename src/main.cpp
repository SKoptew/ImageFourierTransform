#include "BufferConversions.h"
#include "FourierTransform.h"

#include <ctime>
#include <iostream>
#include <iomanip>

int main(int argc, char* argv[])
{
	auto img = new BMP("data/tank_2048.bmp");
	//auto img = new BMP("data/tank_512.bmp");
	//auto img = new BMP("data/stripes.bmp");

	Complex* buffer0 = CreateComplexBuffer(img->Width(), img->Height());
	//Complex* buffer1 = CreateComplexBuffer(img->Width(), img->Height());

	ImageToComplexArray(img, buffer0);

	auto c_start = std::clock();
	{
		//FourierTransform2D(buffer0, buffer1, img->Width(), img->Height());
		//FourierTransform2D(buffer0, buffer1, img->Width(), img->Height(), true);

		//FFT2DRecursive(buffer0, img->Width(), img->Height());
		//FFT2DRecursive(buffer0, img->Width(), img->Height(), true);

		FFT2D(buffer0, img->Width(), img->Height());
		FFT2D(buffer0, img->Width(), img->Height(), true);
	}
	auto c_end = std::clock();

	std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
		<< 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n";

	ComplexArrayToImage(buffer0, img);
	//ComplexArrayLogToImage(buffer0, img);
	
	DisposeComplexBuffer(buffer0);
	//DisposeComplexBuffer(buffer1);
	img->write("out.bmp");

	return 0;
}

