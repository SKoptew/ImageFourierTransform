#pragma once

#include "Complex.h"

void FourierTransform2D(Complex* src, Complex* dst, int w, int h, bool inverse = false);

void FFT2D(Complex* src, int w, int h, bool inverse = false);