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

void ImageToComplexArray(BMP* img, Complex* data)
{
	const int w = img->Width();
	const int h = img->Height();

	float r, g, b;

	for (uint32_t y = 0; y < h; ++y)
	for (uint32_t x = 0; x < w; ++x)
	{
		img->get_pixel(x, y, r, g, b);
		float brightness = 0.2126729f * r + 0.7151522f * g + 0.0721750f * b;

		data[y*w + x] = Complex(brightness, 0.f);
	}
}

// writes magnitude of complex numbers as image value
void ComplexArrayToImage(Complex* data, BMP* img)
{
	const int w = img->Width();
	const int h = img->Height();

	for (uint32_t y = 0; y < h; ++y)
	for (uint32_t x = 0; x < w; ++x)
	{
		const float magnitude = data[y * w + x].abs();

		img->set_pixel(x, y, magnitude, magnitude, magnitude);
	}
}