#include "bn.h"
#include "bn_ext.h"


typedef struct {
    struct bn ndim;
    struct bn nx;
    struct bn nz;
    double **** data; //(**) pointers are
}