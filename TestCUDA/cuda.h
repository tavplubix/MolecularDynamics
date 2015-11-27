
#ifdef QT
#include <QDebug>
#endif

#define checkCudaErrors(val) check( (val), #val, __FILE__, __LINE__)

template<typename T>
void check(T err, const char* const func, const char* const file, const int line)
{
	if (err != cudaSuccess) {
#ifdef QT
		qDebug << "CUDA error at: " << file << ":" << line;
		qDebug << cudaGetErrorString(err) << " " << func;
		qApp->quit(1);
#else
		std::cerr << "CUDA error at: " << file << ":" << line << std::endl;
		std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
		exit(1);
#endif
	}
}













