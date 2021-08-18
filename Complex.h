#pragma once

struct Complex
{
	float re;
	float im;

	Complex() {}

	Complex(float _re, float _im) :
		re(_re),
		im(_im)
	{}

	float abs() const { return sqrtf(re * re + im * im); }
};
