#pragma once

#include "BMP.h"
#include "Complex.h"

Complex* CreateComplexBuffer(int w, int h);
void     DisposeComplexBuffer(Complex* data);

void ImageToComplexArray(BMP* img, Complex* data);
void ComplexArrayToImage(Complex* data, BMP* img);