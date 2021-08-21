#pragma once

#include "Complex.h"

void FT2D_Bruteforce(Complex* src, Complex* tmp, int w, int h, bool inverse = false);

void FFT2D_CT_Recursive  (Complex* src, int w, int h, bool inverse = false); // Cooley-Turkey
void FFT2D_CT_Bitreversal(Complex* src, int w, int h, bool inverse = false);