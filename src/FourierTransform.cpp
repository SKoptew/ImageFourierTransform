#include "FourierTransform.h"
#include <utility>

//-- naive implementation
void FT1D_Bruteforce(Complex* src, Complex* dst, int start, int stride, int N, bool inverse)
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

void FT2D_Bruteforce(Complex* src, Complex* tmp, int w, int h, bool inverse)
{
	//-- horizontal pass
	for (int y = 0; y < h; ++y)
		FT1D_Bruteforce(src, tmp, y*w, 1, w, inverse);

	// vertical pass
	for (int x = 0; x < w; ++x)
		FT1D_Bruteforce(tmp, src, x, w, h, inverse);

	// normalization
	if (!inverse)
	{
		const float inv_N2 = 1.f / (w * h);

		for (int i = 0; i < w * h; ++i)
			src[i] *= inv_N2;
	}
}
//--------------------------------------------------------------------------------------------------------------------------

void FFT1D_CT_Recursive(Complex* src, int start, int stride, int N, bool inverse)
{
	if (N > 1)
	{
		FFT1D_CT_Recursive(src, start,          stride * 2, N/2, inverse);
		FFT1D_CT_Recursive(src, start + stride, stride * 2, N/2, inverse);
		
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

void FFT2D_CT_Recursive(Complex* src, int w, int h, bool inverse)
{
	for (int y = 0; y < h; ++y)
		FFT1D_CT_Recursive(src, y * w, 1, w, inverse);

	for (int x = 0; x < w; ++x)
		FFT1D_CT_Recursive(src, x, w, h, inverse);

	if (!inverse)
	{
		const float inv_N2 = 1.f / (w * h);

		for (int i = 0; i < w * h; ++i)
			src[i] *= inv_N2;
	}
}
//--------------------------------------------------------------------------------------------------------------------------
void FFT1D_CT_Bitreversal(Complex* src, int start, int stride, int N, bool inverse)
{
	//-- bit-reversal permutation of elements
	for (int i = 1, i_rev = 0; i < N; ++i)
	{
		int bit = N >> 1;
		for (; i_rev >= bit; bit >>= 1)
			i_rev -= bit;
	
		i_rev += bit;
	
		if (i < i_rev)
			std::swap(src[start + i * stride], src[start + i_rev * stride]);
	}

	//-- convolve: pairs
	int len = 2;
	{
		for (int i = 0; i < N; i += len)
		{
			int idx_u = start +  i      * stride;
			int idx_v = start + (i + 1) * stride;
	
			Complex u = src[idx_u];
			Complex v = src[idx_v];
	
			src[idx_u] = u + v;
			src[idx_v] = u - v;
		}
	}
	len <<= 1;

	//-- convolve: quadruples, ...
	for (; len <= N; len <<= 1)
	{
		const float   angle = (inverse ? PI2 : -PI2) / len;
		const Complex wlen  = Complex(cos(angle), sin(angle));

		for (int i = 0; i < N; i += len)
		{
			Complex twiddle = { 1.f, 0.f };

			for (int j = 0; j < len / 2; ++j)
			{
				int idx_u = start + (i + j)           * stride;
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

void FFT2D_CT_Bitreversal(Complex* src, int w, int h, bool inverse)
{
	for (int y = 0; y < h; ++y)
		FFT1D_CT_Bitreversal(src, y * w, 1, w, inverse);

	for (int x = 0; x < w; ++x)
		FFT1D_CT_Bitreversal(src, x, w, h, inverse);

	if (!inverse)
	{
		const float inv_N2 = 1.f / (w * h);

		for (int i = 0; i < w * h; ++i)
			src[i] *= inv_N2;
	}
}
//--------------------------------------------------------------------------------------------------------------------------

// http://wwwa.pikara.ne.jp/okojisan/otfft-en/stockham2.html
// N: sequence length
// s: idx stride
// eo : x is output if eo == true, y is output if eo == false
// x, y - swappable buffers. input <-> output

void FFT1D_Stockham(Complex* x, Complex *y, int N, int s, bool eo, bool inverse)
{
	const int M = N / 2;
	const double theta0 = (inverse ? PI2 : -PI2) / N;

	if (N == 1)
	{
		if (eo) // y is output
		{
			for (int q = 0; q < s; ++q)
				y[q] = x[q];
		}
	}
	else
	{
		for (int p = 0; p < M; ++p)
		{
			const float   phase = p * theta0;
			const Complex wp = Complex(cos(phase), sin(phase));

			for (int q = 0; q < s; ++q)
			{
				Complex a = x[q + s*(p + 0)];
				Complex b = x[q + s*(p + M)];

				y[q + s*(2*p + 0)] =  a + b;
				y[q + s*(2*p + 1)] = (a - b) * wp;
			}			
		}
		FFT1D_Stockham(y, x, N/2, 2*s, !eo, inverse); // $$$ tail recursion!
	}
}

void FT2D_Stockham(Complex* src, int w, int h, bool inverse)
{
	//-- horizontal pass
	{
		Complex* row_x = new Complex[w];
		Complex* row_y = new Complex[w];
	
		for (int y = 0; y < h; ++y)
		{
			for (int k = 0; k < w; ++k)
				row_x[k] = src[y*w + k];
	
			FFT1D_Stockham(row_x, row_y, w, 1, false, inverse);
	
			for (int k = 0; k < w; ++k)
				src[y*w + k] = row_x[k];
		}
	
		delete[] row_x;
		delete[] row_y;
	}

	//-- vertical pass
	{
		Complex* col_x = new Complex[h];
		Complex* col_y = new Complex[h];
	
		for (int x = 0; x < w; ++x)
		{
			for (int k = 0; k < h; ++k)
				col_x[k] = src[k*w + x];
	
			FFT1D_Stockham(col_x, col_y, h, 1, false, inverse);
	
			for (int k = 0; k < h; ++k)
				src[k*w + x] = col_x[k];
		}
	
		delete[] col_x;
		delete[] col_y;
	}

	// normalization
	if (!inverse)
	{
		const float inv_N2 = 1.f / (w * h);
	
		for (int i = 0; i < w * h; ++i)
			src[i] *= inv_N2;
	}
}