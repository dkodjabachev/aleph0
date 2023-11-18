// Standard C++ library headers
#include <iostream>

// Aleph0 library headers
#include "vector2d.h"
#include "vector3d.h"
#include "matrix2x2.h"
#include "matrix4x4.h"

using namespace std;
using namespace aleph0;

int main()
{
	cout << "This is the Aleph_0 mathematical library.\n";
	cout << "--------------------------------------------------------------------------------\n";

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

	///////////////////////////////////////////////////////////////////////////////
	//
	// All three ways of calculating the 2x2 matrix product should produce the same result.
	const Matrix2x2 A(1, 2, 3, 4);
	const Matrix2x2 B(4, 5, 6, 7);

	cout << "A = " << A.ToString() << endl;
	cout << "B = " << B.ToString() << endl;
	cout << "A *_naive B = " << (Matrix2x2::NaiveMatrixProduct(A, B)).ToString() << endl;
	cout << "A *_unrolled B = " << (Matrix2x2::UnrolledMatrixProduct(A, B)).ToString() << endl;
	cout << "A *_strassen B = " << (Matrix2x2::StrassenMatrixProduct(A, B)).ToString() << endl;

	assert(Matrix2x2::NaiveMatrixProduct(A, B) == Matrix2x2::UnrolledMatrixProduct(A, B));
	assert(Matrix2x2::UnrolledMatrixProduct(A, B) == Matrix2x2::StrassenMatrixProduct(A, B));

	///////////////////////////////////////////////////////////////////////////////
	//
	// The two ways of calculating the 4x4 matrix product should produce the same result.
	const Matrix4x4 C(Matrix2x2(1, 2, 3, 4), Matrix2x2(5, 6, 7, 8), Matrix2x2(9, 10, 11, 12), Matrix2x2(13, 14, 15, 16));
	const Matrix4x4 D(Matrix2x2(9, 10, 11, 12), Matrix2x2(13, 14, 15, 16), Matrix2x2(1, 2, 3, 4), Matrix2x2(5, 6, 7, 8));

	cout << "C = " << C.ToString() << endl;
	cout << "D = " << D.ToString() << endl;
	cout << "C *_naive D = " << (Matrix4x4::NaiveMatrixProduct(C, D)).ToString() << endl;
	cout << "C *_strassen D = " << (Matrix4x4::StrassenMatrixProduct(C, D)).ToString() << endl;

	assert(Matrix4x4::NaiveMatrixProduct(C, D) == Matrix4x4::StrassenMatrixProduct(C, D));

	///////////////////////////////////////////////////////////////////////////////
	//
	// The determinant of the identity matrix should be 1.
	const Matrix2x2 I(1, 0, 0, 1);

	cout << "I = " << I.ToString() << endl;
	cout << "det(I) = " << I.Determinant() << endl;

	assert(I.Determinant() == 1.0f);

	///////////////////////////////////////////////////////////////////////////////
	//
	// The determinant of a matrix and the determinant of its transpose should be the same.
	const Matrix2x2 M(1, 2, 3, 4);

	cout << "M = " << M.ToString() << endl;
	cout << "det(M) = " << M.Determinant() << endl;
	cout << "det(M') = " << M.Transposed().Determinant() << endl;

	assert(M.Determinant() == M.Transposed().Determinant());

	///////////////////////////////////////////////////////////////////////////////
	//
	// If we add a multiple of some row/column to another row/column, the determinant should not change.
	// Here, we have added 3 times the first row of M to the second row of M to obtain N.
	const Matrix2x2 N(1, 2, 6, 10);

	cout << "N = " << N.ToString() << endl;
	cout << "det(N) = " << N.Determinant() << endl;

	assert(N.Determinant() == M.Determinant());

	///////////////////////////////////////////////////////////////////////////////
	//
	// The determinant of a matrix with two identical rows/columns should be zero.
	const Matrix2x2 L(1, 2, 1, 2);

	cout << "L = " << L.ToString() << endl;
	cout << "det(L) = " << L.Determinant() << endl;

	assert(L.Determinant() == 0.f);

	return 0;
}