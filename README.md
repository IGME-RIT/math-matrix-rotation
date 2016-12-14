# Matrix Rotation

This tutorial explores rotation matrices and uses several results and properties from linear algebra to show how nicely they behave.

This tutorial also derives one variant of Rodrigues' rotation formula.

Rodrigues' Rotation Formula is a matrix formula for a rotation in 3D space centered at the origin about an arbitrary axis of rotation and arbitrary angle.

Given an arbitrary unit vector __k__ and angle &theta;, Rodrigues' Rotation formula gives us that the matrix for the rotation about __k__ by &theta; radians is

R = I<sub>3</sub> + sin(&theta;)\[__k__\]<sub>&times;</sub> + (1-cos(&theta;))\[__k__\]<sub>&times;</sub><sup>2</sup>,

or

R = cos(&theta;)I<sub>3</sub> + sin(&theta;)\[__k__\]<sub>&times;</sub> + (1-cos(&theta;))__k__&otimes;__k__.

where I<sub>3</sub> is the 3D identity matrix, \[__k__\]<sub>&times;</sub> is the [skew-symmetric cross product matrix](<https://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication>), and __k__&otimes;__k__ is the [outer product](<https://en.wikipedia.org/wiki/Outer_product>) of __k__ with itself.

To see that the two are equivalent, note that if ||__k__|| = 1, then \[__k__\]<sub>&times;</sub><sup>2</sup> = __k__&otimes;__k__ - I<sub>3</sub>

# Setup

You will need to have CMake installed on your computer, and properly added to your path. In order to setup, run the following in a shell, then open the project in your preferred editor. Windows setup has been configured for use with Visual Studio.

Windows:
```
cd path/to/folder
setup.cmd
```
Linux:
```
cd path/to/folder
./setup
```
