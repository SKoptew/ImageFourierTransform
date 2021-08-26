#pragma once

#include "BMP.h"
#include "Complex.h"

Complex* CreateComplexBuffer(int w, int h);
void     DisposeComplexBuffer(Complex* data);

void ImageToComplexArray(BMP* img, Complex* data0, Complex* data1);
void ComplexArrayToImage(BMP* img, Complex* data0, Complex* data1);
