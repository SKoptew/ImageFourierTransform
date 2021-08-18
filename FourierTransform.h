#pragma once

#include "Complex.h"

void FourierTransform2D(Complex* src, Complex* dst, int w, int h);

void FourierTransform1DHoriz(Complex* src, Complex* dst, int N);
void FourierTransform1DVert(Complex* src, Complex* dst, int start, int stride, int N);