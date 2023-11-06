#pragma once

// Standard C++ library headers
#include <cassert>
#include <string>

// Aleph0 library headers
#include "vector2d.h"

namespace aleph0
{

class Matrix2x2
{
public:
	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2()
	{
		// The default constructor of the union Data will do the job
		// of initializing the matrix entries to zero.
	}

	Matrix2x2(const Matrix2x2& other)
	{
		// Use the Vector2d operator=
		data_.rows_[0] = other.data_.rows_[0];
		data_.rows_[1] = other.data_.rows_[1];
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2(float a, float b, float c, float d)
	{
		(*this)[0][0] = a;
		(*this)[0][1] = b;
		(*this)[1][0] = c;
		(*this)[1][1] = d;
	}

public:
	///////////////////////////////////////////////////////////////////////////////
	void SetZeroMatrix()
	{
		memset(data_.elements_, 0, 4 * sizeof(float));
	}

	///////////////////////////////////////////////////////////////////////////////
	void SetIdentityMatrix()
	{
		(*this)[0][0] = 1.f;
		(*this)[0][1] = 0.f;
		(*this)[1][0] = 0.f;
		(*this)[1][1] = 1.f;
	}

	///////////////////////////////////////////////////////////////////////////////
	const float* operator[](const size_t row_index) const
	{
		assert(row_index < 2);

		return data_.array_[row_index];
	}

	///////////////////////////////////////////////////////////////////////////////
	float* operator[](const size_t row_index)
	{
		assert(row_index < 2);

		return data_.array_[row_index];
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2& operator=(const Matrix2x2& other)
	{
		(*this)[0][0] = other[0][0];
		(*this)[0][1] = other[0][1];
		(*this)[1][0] = other[1][0];
		(*this)[1][1] = other[1][1];

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2 operator+(const Matrix2x2& other) const
	{
		Matrix2x2 result;

		result[0][0] = (*this)[0][0] + other[0][0];
		result[0][1] = (*this)[0][1] + other[0][1];
		result[1][0] = (*this)[1][0] + other[1][0];
		result[1][1] = (*this)[1][1] + other[1][1];

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2 operator+=(const Matrix2x2& other)
	{
		(*this) = (*this) + other;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2 operator-(const Matrix2x2& other) const
	{
		Matrix2x2 result;

		result[0][0] = (*this)[0][0] - other[0][0];
		result[0][1] = (*this)[0][1] - other[0][1];
		result[1][0] = (*this)[1][0] - other[1][0];
		result[1][1] = (*this)[1][1] - other[1][1];

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2& operator-=(const Matrix2x2& other)
	{
		(*this) = (*this) - other;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2 operator*(const float scalar) const
	{
		Matrix2x2 result;

		result[0][0] = (*this)[0][0] * scalar;
		result[0][1] = (*this)[0][1] * scalar;
		result[1][0] = (*this)[1][0] * scalar;
		result[1][1] = (*this)[1][1] * scalar;

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2& operator*=(const float scalar)
	{
		(*this) = (*this) * scalar;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2 operator*(const Matrix2x2& other) const
	{
		return NaiveMatrixProduct(*this, other);
	}

	///////////////////////////////////////////////////////////////////////////////
	static Matrix2x2 NaiveMatrixProduct(const Matrix2x2& A, const Matrix2x2& B)
	{
		Matrix2x2 result;

		for (size_t i = 0; i < 2; ++i)
		{
			for (size_t j = 0; j < 2; ++j)
			{
				for (size_t k = 0; k < 2; ++k)
				{
					result[i][j] += A[i][k] * B[k][j];
				}
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	static Matrix2x2 UnrolledMatrixProduct(const Matrix2x2& A, const Matrix2x2& B)
	{
		Matrix2x2 result;

		// Same as above, but without the three for loops, i.e. with loops unrolled.
		// Notice that there are 8 multiplications in total.
		result[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0];
		result[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1];
		result[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0];
		result[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1];

		return result;
	}

	static Matrix2x2 StrassenMatrixProduct(const Matrix2x2& A, const Matrix2x2& B)
	{
		// Strassen's algorithm for 2x2 matrices. Note that we have decreased the number
		// of multiplications to 7, at the cost of more addition and substraction operations.
		const float A_11 = A[0][0];
		const float A_12 = A[0][1];
		const float A_21 = A[1][0];
		const float A_22 = A[1][1];
		const float B_11 = B[0][0];
		const float B_12 = B[0][1];
		const float B_21 = B[1][0];
		const float B_22 = B[1][1];

		const float M_1 = (A_11 + A_22) * (B_11 + B_22);
		const float M_2 = (A_21 + A_22) * B_11;
		const float M_3 = A_11 * (B_12 - B_22);
		const float M_4 = A_22 * (B_21 - B_11);
		const float M_5 = (A_11 + A_12) * B_22;
		const float M_6 = (A_21 - A_11) * (B_11 + B_12);
		const float M_7 = (A_12 - A_22) * (B_21 + B_22);

		const float C_11 = (M_1 + M_4 - M_5 + M_7);
		const float C_12 = (M_3 + M_5);
		const float C_21 = (M_2 + M_4);;
		const float C_22 = (M_1 - M_2 + M_3 + M_6);

		const Matrix2x2 result(C_11, C_12, C_21, C_22);

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	std::string ToString() const
	{
		std::string result =
			"\n" + data_.rows_[0].ToString() +
			"\n" + data_.rows_[1].ToString();

		return result;
	}

private:
	///////////////////////////////////////////////////////////////////////////////
	union Data
	{
		Data()
		{
			memset(elements_, 0, 4 * sizeof(float));
		}

		// We can represent our 2x2 matrix in a number of ways:

		// 1) As structure of four separate float variables
		struct
		{
			float a_;
			float b_;
			float c_;
			float d_;
		};

		// 2) As an array of four floats.
		float elements_[4];

		// 3) As a 2x2 array of floats.
		float array_[2][2];

		// 4) Even as an array of two Vector2d objects.
		Vector2d rows_[2];
	};

	Data data_;
}; // class Matrix2x2

} // namespace aleph0