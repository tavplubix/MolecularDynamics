#include <iostream>
#include <conio.h>
#include "Vector.h"

int main() {
	Vector v(0, 0, 0, 1, 1, 1);
	Vector v1(v);
	v1 = v + v1;

	std::cout << v << " " << v1;
	_getch(); // Ожидание нажатия клавиши.
	return 0;
}