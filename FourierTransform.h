#pragma once

#include "Complex.h"

void FourierTransform2D(Complex* src, Complex* dst, int w, int h, bool inverse = false);

void FourierTransform1DHoriz(Complex* src, Complex* dst, int N, bool inverse);
void FourierTransform1DVert (Complex* src, Complex* dst, int start, int stride, int N, bool inverse);