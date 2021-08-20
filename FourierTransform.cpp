#include "FourierTransform.h"

void FourierTransform1D(Complex* src, Complex* dst, int start, int stride, int N, bool inverse)
{
	for (int k = 0; k < N; ++k) // calc column of FFT indexes
	{
		Complex sum = { 0.f, 0.f };
		const float phase_mul = (inverse ? PI2 : -PI2) * k / N;

		for (int a = 0, idx = start; a < N; ++a, idx += stride)
		{
			float phase = a * phase_mul;
			Complex twiddle = Complex(cos(phase), sin(phase));

			sum = sum + src[idx] * twiddle;
		}

		dst[start + k * stride] = sum;
	}
}

void FourierTransform2D(Complex* src, Complex* tmp, int w, int h, bool inverse)
{
	//-- horizontal pass
	for (int y = 0; y < h; ++y)
		FourierTransform1D(src, tmp, y*w, 1, w, inverse);

	// vertical pass
	for (int x = 0; x < w; ++x)
		FourierTransform1D(tmp, src, x, w, h, inverse);

	// normalization
	if (!inverse)
	{
		const float inv_N2 = 1.f / (w * h);

		for (int i = 0; i < w * h; ++i)
			src[i] *= inv_N2;
	}
}

//--------------------------------------------------------------------------------------------------------------------------

void FFT1D(Complex* src, int start, int stride, int N, bool inverse)
{
	if (N > 1)
	{
		FFT1D(src, start,          stride * 2, N/2, inverse);
		FFT1D(src, start + stride, stride * 2, N/2, inverse);
		
		Complex twiddle = { 1.f, 0.f };
		float   angle   = (inverse ? PI2 : -PI2) / N;
		Complex wn      = Complex(cos(angle), sin(angle));

		for (int k = 0; k < N / 2; ++k)
		{
			Complex p = src[start +  k        * stride] * twiddle;
			Complex q = src[start + (k + N/2) * stride] * twiddle;

			src[start +  k        * stride] = p + q;
			src[start + (k + N/2) * stride] = p - q;

			twiddle = twiddle * wn;
		}
	}
}

void FFT2D(Complex* src, int w, int h, bool inverse)
{
	// horizontal pass
	for (int y = 0; y < h; ++y)
		FFT1D(src, y * w, 1, w, inverse);

	// vertical pass
	for (int x = 0; x < w; ++x)
		FFT1D(src, x, w, h, inverse);

	// normalization
	if (!inverse)
	{
		const float inv_N2 = 1.f / (w * h);

		for (int i = 0; i < w * h; ++i)
			src[i] *= inv_N2;
	}
}
