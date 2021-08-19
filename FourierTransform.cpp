#include "FourierTransform.h"

void FourierTransform2D(Complex* src, Complex* tmp, int w, int h, bool inverse)
{
	//-- horizontal pass
	for (int y = 0; y < h; ++y)
		FourierTransform1D(src, tmp, y*w, 1, w, inverse);

	// vertical pass
	for (int x = 0; x < w; ++x)
		FourierTransform1D(tmp, src, x, w, h, inverse);
}

void FourierTransform1D(Complex* src, Complex* dst, int start, int stride, int N, bool inverse)
{
	for (int k = 0; k < N; ++k) // calc column of FFT indexes
	{		
		Complex sum = { 0.f, 0.f };
		const float phase_mul = (inverse ? PI2 : -PI2) * k / N;

		for (int a = 0, idx = start; a < N; ++a, idx += stride)
		{
			float phase = a * phase_mul;
			Complex e = Complex(cos(phase), sin(phase));

			sum = sum + src[idx] * e;
		}

		if (!inverse)
		{
			sum.re /= N;
			sum.im /= N;
		}

		dst[start + k*stride] = sum;
	}
}
