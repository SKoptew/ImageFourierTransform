#pragma once

#include <cmath>
#include <float.h>

struct Complex
{
	float re;
	float im;

	Complex() {}
	Complex(float _re, float _im) : re(_re), im(_im) {}

	//-- unary
	Complex operator - () const { return Complex(-re, -im); }

	//-- binary
	Complex operator + (const Complex& rhs) const { return Complex(re + rhs.re, im + rhs.im); }
	Complex operator - (const Complex& rhs) const { return Complex(re - rhs.re, im - rhs.im); }
	Complex operator * (const Complex& rhs) const { return Complex(re*rhs.re - im*rhs.im, im*rhs.re + re*rhs.im); }

	Complex& operator *= (float rhs) { re *= rhs; im *= rhs; return *this; }

	float abs()   const { return sqrt(re * re + im * im); }
	float angle() const { return atan2(im, re); }
};

const float   PI2 = 3.14159265358979323846 * 2.f;
const Complex I   = Complex(-1.f, 0.f);