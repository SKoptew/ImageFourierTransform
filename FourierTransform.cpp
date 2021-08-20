#include "FourierTransform.h"
#include <utility>

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

void FFT1DRecursive(Complex* src, int start, int stride, int N, bool inverse)
{
	if (N > 1)
	{
		FFT1DRecursive(src, start,          stride * 2, N/2, inverse);
		FFT1DRecursive(src, start + stride, stride * 2, N/2, inverse);
		
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

void FFT2DRecursive(Complex* src, int w, int h, bool inverse)
{
	// horizontal pass
	for (int y = 0; y < h; ++y)
		FFT1DRecursive(src, y * w, 1, w, inverse);

	// vertical pass
	for (int x = 0; x < w; ++x)
		FFT1DRecursive(src, x, w, h, inverse);

	// normalization
	if (!inverse)
	{
		const float inv_N2 = 1.f / (w * h);

		for (int i = 0; i < w * h; ++i)
			src[i] *= inv_N2;
	}
}
//--------------------------------------------------------------------------------------------------------------------------

int reverse_last_bits(int num, int lg_n)
{
	int res = 0;

	for (int i = 0; i < lg_n; ++i)
	if (num & (1 << i))
	{
		res |= 1 << (lg_n - 1 - i);
	}
	return res;
}

void FFT1D(Complex* src, int start, int stride, int N, bool inverse)
{
	//-- get lg(N)
	int lg_n = 0;
	while ((1 << lg_n) < N)
		++lg_n;

	//-- reorder elements
	for (int i = 0; i < N; ++i)
	{
		int i_rev = reverse_last_bits(i, lg_n);

		if (i < i_rev)
		{
			std::swap(src[start + i * stride], src[start + i_rev * stride]);
		}
	}

	//-- convolve: pairs, quadruples, ...
	for (int len = 2; len <= N; len <<= 1)
	{
		const float   angle = (inverse ? PI2 : -PI2) / len;
		const Complex wlen  = Complex(cos(angle), sin(angle));

		for (int i = 0; i < N; i += len)
		{
			Complex twiddle = { 1.f, 0.f };

			for (int j = 0; j < len / 2; ++j)
			{
				int idx_u = start + (i + j) * stride;
				int idx_v = start + (i + j + len / 2) * stride;

				Complex u = src[idx_u];
				Complex v = src[idx_v] * twiddle;

				src[idx_u] = u + v;
				src[idx_v] = u - v;

				twiddle = twiddle * wlen;
			}
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