#include "BufferConversions.h"
#include "FourierTransform.h"

int main(int argc, char* argv[])
{
	//auto img = new BMP("tank.bmp");
	auto img = new BMP("tank_512.bmp");
	//auto img = new BMP("stripes.bmp");


	Complex* buffer0 = CreateComplexBuffer(img->Width(), img->Height());
	Complex* buffer1 = CreateComplexBuffer(img->Width(), img->Height());

	ImageToComplexArray(img, buffer0);

	//FourierTransform2D(buffer0, buffer1, img->Width(), img->Height());
	//FourierTransform2D(buffer0, buffer1, img->Width(), img->Height(), true);

	FFT2D(buffer0, img->Width(), img->Height());
	FFT2D(buffer0, img->Width(), img->Height(), true);

	ComplexArrayToImage(buffer0, img);
	//ComplexArrayLogToImage(buffer0, img);
	
	DisposeComplexBuffer(buffer0);
	DisposeComplexBuffer(buffer1);
	img->write("out.bmp");

	return 0;
}

