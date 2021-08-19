#include "FourierTransform.h"

void FourierTransform2D(Complex* src, Complex* tmp, int w, int h, bool inverse)
{
	//-- horizontal pass
	for (int y = 0; y < h; ++y)
		FourierTransform1DHoriz(&(src[y * w]), &(tmp[y * w]), w, inverse);

	// vertical pass
	for (int x = 0; x < w; ++x)
		FourierTransform1DVert(tmp, src, x, w, h, inverse);
}

void FourierTransform1DHoriz(Complex* src, Complex* dst, int N, bool inverse)
{
	for (int k = 0; k < N; ++k) // calc row of FFT indexes
	{
		Complex sum = { 0.f, 0.f };

		const float phase_mul = (inverse ? PI2 : -PI2) * k / N;

		for (int a = 0; a < N; ++a)
		{
			float phase = a * phase_mul;
			auto e = Complex(cos(phase), sin(phase));

			auto res = src[a] * e;
			sum = sum + res;
		}

		if (!inverse)
		{
			sum.re /= N;
			sum.im /= N;
		}

		dst[k] = sum;
	}
}

void FourierTransform1DVert(Complex* src, Complex* dst, int start, int stride, int N, bool inverse)
{
	for (int l = 0; l < N; ++l) // calc column of FFT indexes
	{		
		Complex sum = { 0.f, 0.f };
		const float phase_mul = (inverse ? PI2 : -PI2) * l / N;

		for (int b = 0, idx = start; b < N; ++b, idx += stride)
		{
			float phase = b * phase_mul;
			auto e = Complex(cos(phase), sin(phase));

			auto res = src[idx] * e;
			sum = sum + res;
		}

		if (!inverse)
		{
			sum.re /= N;
			sum.im /= N;
		}

		dst[start + l*stride] = sum;
	}
}
