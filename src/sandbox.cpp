#include <iostream>
#include "vector2d.h"
#include "vector3d.h"

using namespace std;
using namespace aleph0;

int main()
{
	cout << "This is the Aleph_0 mathematical library." << endl << endl;

	Vector3d v_1;
	Vector3d v_2(1.f, 2.f, 3.f);
	Vector3d v_3(v_2 * (-1));
	Vector3d v_4(Vector2d(0.f, 1.f), 1.f);

	cout << "v_1 = " << v_1.ToString() << endl;
	cout << "v_2 = " << v_2.ToString() << endl;
	cout << "v_3 = " << v_3.ToString() << endl;
	cout << "v_3 - v_2 = " << (v_3 + v_2).ToString() << endl;
	cout << "v_3 - v_2 == v_1: " << ((v_3 + v_2) == v_1 ? "yes" : "no") << endl;
	cout << "v_4 = " << v_4.ToString() << endl;

	return 0;
}