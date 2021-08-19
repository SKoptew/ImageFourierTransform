#pragma once

#include "Complex.h"

void FourierTransform2D(Complex* src, Complex* dst, int w, int h, bool inverse = false);
void FourierTransform1D(Complex* src, Complex* dst, int start, int stride, int N, bool inverse);