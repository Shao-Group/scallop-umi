#include "RO_read.h"

RO_read::RO_read()
{

}

RO_read::~RO_read()
{
}

int RO_read::print()
{
    printf("read:%s, chrm:%s, BSJ:%d-%d\n",read_name.c_str(),chrm.c_str(),BSJ_pos,BSJ_rpos);
    return 0;
}
