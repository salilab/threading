Place the private header files in this directory. They will be
available to your code with

     #include <IMP/threading/internal/myheader.h>

All headers should include `IMP/threading/threading_config.h` as their
first include and surround all code with `IMPTHREADING_BEGIN_INTERNAL_NAMESPACE`
and `IMPTHREADING_END_INTERNAL_NAMESPACE` to put it in the
IMP::threading::internal namespace and manage compiler warnings.
