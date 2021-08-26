#include "BufferConversions.h"

Complex* CreateComplexBuffer(int w, int h)
{
	return new Complex[w * h];
}

void DisposeComplexBuffer(Complex* data)
{
	delete[] data;
	data = nullptr;
}

void ImageToComplexArray(BMP* img, Complex* data0, Complex* data1)
{
	const int w = img->Width();
	const int h = img->Height();

	float r, g, b;

	for (uint32_t y = 0; y < h; ++y)
	for (uint32_t x = 0; x < w; ++x)
	{
		img->get_pixel(x, y, r, g, b);

		data0[y*w + x] = Complex(r, g);
		data1[y*w + x] = Complex(b, 0.f);
	}
}

void ComplexArrayToImage(BMP* img, Complex* data0, Complex* data1)
{
	const int w = img->Width();
	const int h = img->Height();

	for (uint32_t y = 0; y < h; ++y)
	for (uint32_t x = 0; x < w; ++x)
	{
		//const float magnitude = data0[y * w + x].abs();
		const Complex v0 = data0[y * w + x];
		const Complex v1 = data1[y * w + x];

		img->set_pixel(x, y, v0.re, v0.im, v1.re);
	}
}
