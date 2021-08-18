#include "BufferConversions.h"

int main(int argc, char* argv[])
{
	auto img = new BMP("tank.bmp");


	Complex* data = CreateComplexBuffer(img->Width(), img->Height());
	ImageToComplexArray(img,  data);
	ComplexArrayToImage(data, img);
	
	DisposeComplexBuffer(data);
	img->write("out.bmp");

	return 0;
}

