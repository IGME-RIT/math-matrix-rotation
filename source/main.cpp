/*
Title: Matrix Mathematics
File Name: main.cpp
Copyright © 2016
Author: Andrew Litfin
Written under the supervision of David I. Schwartz, Ph.D., and
supported by a professional development seed grant from the B. Thomas
Golisano College of Computing & Information Sciences
(https://www.rit.edu/gccis) at the Rochester Institute of Technology.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// The primary objects of study in linear algebra are matrices.
// This tutorial series will explore the applications of matrices to computer games and simulation,
//  especially in the realm of physical transformations.
// The exposition follows that of Eric Lengyel in "Foundations of Game Engine Development" (Volume 1).
// We have included the Vector structs from the previous series, and introduced Matrix structs that act similarly.
// These structs are based upon and largely follow code samples given in FGED.
//  As before, Matrix2D is heavily annotated, with other structs being annotated in places of difference.

// This tutorial explores rotation matrices.

#include "../header/Matrix4D.h"
#include "../header/tests.h"
#include "../header/helpers.h"

#include <iostream>
#include <ctime>

int main()
{
	// Required for some helper functions
	srand((unsigned)time(0));

	// Rotations are one of the most common types of transforms used in games.
	// Whenever an object spins, a hinged door is opened, etc., a rotation transformation is being applied.
	// A rotation matrix is just that--a matrix that encodes a transformation of a vector being rotated about an axis.
	// We most often rotate around coordinate axis but it is possible to rotate about an arbitrary axis.
	// To find what a rotation does, we can consider what happens to the component in the direction of the axis,
	//  and what happens to the component orthogonal to the axis.
	// As per convention, we consider a positive angle of rotation as counter-clockwise, and a negative angle of rotation as clockwise
	//  when viewed down the axis of rotation.

	// First, let's consider what happens in 2D space.
	// As mentioned in a previous tutorial, we can construct a matrix by its columns by considering what happens to the standard basis vectors.
	// So, when undergoing a rotation by theta degrees, the vector (1, 0) goes to (cos(theta), sin(theta)).
	// Similarly, the vector (0, 1) goes to (-sin(theta), cos(theta)).
	// Hence for an arbitrary vector (x, y), the matrix
	// [ cos(theta) -sin(theta) ]
	// [ sin(theta)  cos(theta) ]
	// encodes a rotation by theta degrees.
	// (To see this, draw a coordinate grid yourself, and consider where a vector should go under rotation.)

	float pi4 = 3.14159265f / 4.0f;
	Matrix2D rot45_2d = MakeRotation(pi4);
	std::cout << "In 2D\n-----\n"
		<< "e1 = " << Vector2D(1, 0) << "\nrot45_2d * e1 = " << rot45_2d * Vector2D(1, 0) << "\n"
		<< "e2 = " << Vector2D(0, 1) << "\nrot45_2d * e2 = " << rot45_2d * Vector2D(0, 1) << "\n";

	// We can embed this in 3D space by letting the z axis be the axis of rotation.
	// Then the component of any vector in the direction of the z axis should be unchanged under rotation,
	//  while the component perpendicular to the z-axis (i.e. the part in the xy plane) should undergo the rotation above.
	// Hence we have the rotation
	// [ cos -sin 0 ]
	// [ sin  cos 0 ]
	// [  0    0  1 ]

	Matrix3D rot45_Z = MakeRotationZ(pi4);
	std::cout << "In 3D\n-----\n"
		<< "e1 = " << Vector3D(1, 0, 0) << "\nrot45_Z * e1 = " << rot45_Z * Vector3D(1, 0, 0) << "\n"
		<< "e2 = " << Vector3D(0, 1, 0) << "\nrot45_Z * e2 = " << rot45_Z * Vector3D(0, 1, 0) << "\n"
		<< "e3 = " << Vector3D(0, 0, 1) << "\nrot45_Z * e3 = " << rot45_Z * Vector3D(0, 0, 1) << "\n";

	// The same argument can be used for the other two coordinate axis as well.

	Vector3D v(1, 1, 1);
	Matrix3D rot45_X = MakeRotationX(pi4);
	Matrix3D rot45_Y = MakeRotationY(pi4);
	std::cout << "v = " << v << "\n"
		<< "rot45_X * v = " << rot45_X * v << "\n"
		<< "rot45_Y * v = " << rot45_Y * v << "\n"
		<< "rot45_Z * v = " << rot45_Z * v << "\n";

	// At this point we can consider rotations about an arbitrary axis.
	// Let a be the unit length axis of rotation, v be our vector, and theta be the angle of rotation.
	// Then the component of v onto a remains unchanged during a rotation, but the rejection of v from a is what needs to be rotated.
	// It can be expressed as a linear combination of Reject(v, a) and another vector orthogonal to both that lies in the plane.
	// We choose Cross(a, v) to be that vector. Note that Magnitude(Reject(v, a)) == Cross(a, v) == Magnitude(v) * sin(alpha),
	//  where alpha is the angle between a and v.
	// Then as before if we equate Reject(v, a) to e1 and Cross(a, v) to e2 then the result after rotation is
	//  cos(theta) * Reject(v, a) + Cross(a, v) * sin(theta).
	// Hence our final rotation is
	//  v' = Project(v, a) + Reject(v, a) * cos(theta) + Cross(a, v) * sin(theta).
	//     = MakeProjection(a) * v + MakeRejection(a) * v * cos(theta) + MatCross(a) * v * sin(theta)
	//     = (MakeProjection(a) + MakeRejection(a) * cos(theta) + MatCross(a) * sin(theta)) * v
	//     = (Matrix3D() * cos(theta) + Outer(a) * (1 - cos(theta)) + MatCross(a) * sin(theta)) * v
	// This is one of a few equivalent statements of what is called `Rodrigues' Rotation Formula.'
	// Another equivalent expression is
	//  v' = (Matrix3D() + sin(theta) * MatCross(a) + (1 - cos(theta)) * MatCross(a) * MatCross(a)) * v
	// (To see that they are equivalent, note that MatCross(a) * MatCross(a) == Outer(a) + Matrix3D() and play around with that).

	// For example,

	Vector3D a(1, 1, 0);
	float theta = pi4;
	Matrix3D R = MakeRotation(theta, a);

	std::cout << "a = " << a << ", theta = " << theta << "\n"
		<< "R = MakeRotation(theta, a) =\n" << R
		<< "R * v = " << R * v << "\n";

	// Rotation matrices themselves have several interesting properties.

	// First, they always have determinant +1.
	// This means that they are volume and orientation preserving transformations.
	std::cout << "Determinant(R) = " << Determinant(R) << "\n";

	// Second, they are orthogonal matrices.
	// This means that their inverse is their transpose.
	// e.g.
	std::cout << "Inverse(R) =\n" << Inverse(R)
		<< "Transpose(R) =\n" << Transpose(R);

	std::cout << "\nPress Enter to exit . . . ";
	std::cin.get();
	return 0;
}
