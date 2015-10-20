#include <iostream>
#include <conio.h>
#include "Vector.h"

int main() {
	Vector v(0,0,0,1,2,3);
	Vector v1(v);

	std::cout << v1.GetABS();
	_getch(); // Ожидание нажатия клавиши.
	return 0;
}