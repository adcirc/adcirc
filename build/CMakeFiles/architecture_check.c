#include <stdint.h>
int main()
{
#if INTPTR_MAX == INT64_MAX
        return 0;
#elif INTPTR_MAX == INT32_MAX
#error 32-bit max integer pointer
#else
#error Unknown pointer size or missing size macros
#endif
}