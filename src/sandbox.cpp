// Standard C++ library headers
#include <cassert>
#include <iostream>
#include <random>

// Aleph0 library headers
#include "vector2d.h"
#include "vector3d.h"
#include "matrix2x2.h"
#include "matrix3x3.h"
#include "matrix4x4.h"

using namespace std;
using namespace aleph0;

int main()
{
	cout << "This is the Aleph_0 mathematical library.\n";
	cout << "--------------------------------------------------------------------------------\n";

	///////////////////////////////////////////////////////////////////////////////
	//
	// Check the basic 2d and 3d vector operations.
	{
		cout << "\nLet's see if our Vector3d class works as expected:\n";
		cout << "--------------------------------------------------------------------------------\n";

		const Vector3d v_1;
		const Vector3d v_2(1.f, 2.f, 3.f);
		const Vector3d v_3(v_2 * (-1));
		const Vector3d v_4(Vector2d(0.f, 1.f), 1.f);

		cout << "v_1 = " << v_1.ToString() << endl;
		cout << "v_2 = " << v_2.ToString() << endl;
		cout << "v_3 = " << v_3.ToString() << endl;
		cout << "v_3 - v_2 = " << (v_3 + v_2).ToString() << endl;
		cout << "v_3 - v_2 == v_1: " << ((v_3 + v_2) == v_1 ? "yes" : "no") << endl;
		cout << "v_4 = " << v_4.ToString() << endl;

		cout << "\nPerforming some tests with the Matrix2x2 and Matrix4x4 classes:\n";
		cout << "--------------------------------------------------------------------------------\n";
	}

	///////////////////////////////////////////////////////////////////////////////
	//
	// All three ways of calculating the 2x2 matrix product should produce the same result.
	{
		const Matrix2x2 A(1, 2, 3, 4);
		const Matrix2x2 B(5, 6, 7, 8);

		cout << "A = " << A.ToString() << endl;
		cout << "B = " << B.ToString() << endl;
		cout << "A *_naive B = " << (Matrix2x2::NaiveMatrixProduct(A, B)).ToString() << endl;
		cout << "A *_unrolled B = " << (Matrix2x2::UnrolledMatrixProduct(A, B)).ToString() << endl;
		cout << "A *_strassen B = " << (Matrix2x2::StrassenMatrixProduct(A, B)).ToString() << endl;

		assert(Matrix2x2::NaiveMatrixProduct(A, B) == Matrix2x2::UnrolledMatrixProduct(A, B));
		assert(Matrix2x2::UnrolledMatrixProduct(A, B) == Matrix2x2::StrassenMatrixProduct(A, B));
	}

	///////////////////////////////////////////////////////////////////////////////
	//
	// The two ways of calculating the 4x4 matrix product should produce the same result.
	{
		const Matrix4x4 A(Matrix2x2(1, 2, 3, 4), Matrix2x2(5, 6, 7, 8), Matrix2x2(9, 10, 11, 12), Matrix2x2(13, 14, 15, 16));
		const Matrix4x4 B(Matrix2x2(9, 10, 11, 12), Matrix2x2(13, 14, 15, 16), Matrix2x2(1, 2, 3, 4), Matrix2x2(5, 6, 7, 8));

		cout << "A = " << A.ToString() << endl;
		cout << "B = " << B.ToString() << endl;
		cout << "A *_naive B = " << (Matrix4x4::NaiveMatrixProduct(A, B)).ToString() << endl;
		cout << "A *_strassen B = " << (Matrix4x4::StrassenMatrixProduct(A, B)).ToString() << endl;

		assert(Matrix4x4::NaiveMatrixProduct(A, B) == Matrix4x4::StrassenMatrixProduct(A, B));
	}

	///////////////////////////////////////////////////////////////////////////////
	//
	// The determinant of the identity matrix should be 1.
	{
		const Matrix2x2 I(1, 0, 0, 1);

		cout << "I = " << I.ToString() << endl;
		cout << "det(I) = " << I.Determinant() << endl;

		assert(I.Determinant() == 1.0f);
	}

	///////////////////////////////////////////////////////////////////////////////
	//
	// The determinant of a matrix and the determinant of its transpose should be the same.
	{
		const Matrix2x2 A(1, 2, 3, 4);

		cout << "A = " << A.ToString() << endl;
		cout << "det(A) = " << A.Determinant() << endl;
		cout << "det(A') = " << A.Transposed().Determinant() << endl;

		assert(A.Determinant() == A.Transposed().Determinant());
	}

	///////////////////////////////////////////////////////////////////////////////
	//
	// Adding a multiple of some row/column to another row/column of a matrix should not change the determinant.
	// Here, we have added 3 times the first row of A to the second row of A to obtain B.
	{
		const Matrix2x2 A(1, 2, 3, 4);
		const Matrix2x2 B(1, 2, 6, 10);

		cout << "B = " << B.ToString() << endl;
		cout << "det(B) = " << B.Determinant() << endl;

		assert(B.Determinant() == B.Determinant());
	}

	///////////////////////////////////////////////////////////////////////////////
	//
	// The determinant of a matrix with two identical rows/columns should be zero.
	{
		const Matrix2x2 L(1, 2, 1, 2);

		cout << "L = " << L.ToString() << endl;
		cout << "det(L) = " << L.Determinant() << endl;

		assert(L.Determinant() == 0.f);
	}

	///////////////////////////////////////////////////////////////////////////////
	//
	// Check if inverting 3x3 matrices works as expected.
	{
		random_device seed_provider;
		mt19937 random_generator(seed_provider());
		uniform_real_distribution<float> distribution(0, 1);

		Matrix3x3 M;

		for (int i = 0; i < 9; ++i)
		{
			M(i) = distribution(random_generator);
		}

		Matrix3x3 I; I.SetIdentityMatrix();

		cout << "M = " << M.ToString() << endl;
		cout << "M^-1 = " << M.Inverse().ToString() << endl;
		cout << "M * M^-1" << (M * M.Inverse()).ToString() << endl;

		assert((M * M.Inverse()).IsApproximatelyEqual(I, 0.00001f));
	}

	///////////////////////////////////////////////////////////////////////////////
	//
	// Check if inverting 4x4 matrices works as expected.
	{
		random_device seed_provider;
		mt19937 random_generator(seed_provider());
		uniform_real_distribution<float> distribution(0, 1);
		Matrix4x4 M;

		for (int i = 0; i < 16; ++i)
		{
			M(i) = distribution(random_generator);
		}

		Matrix4x4 I; I.SetIdentityMatrix();

		cout << "M = " << M.ToString() << endl;
		cout << "M^-1 = " << M.Inverse().ToString() << endl;
		cout << "M * M^-1" << (M * M.Inverse()).ToString() << endl;

		assert((M * M.Inverse()).IsApproximatelyEqual(I, 0.00001f));
	}

	return 0;
}