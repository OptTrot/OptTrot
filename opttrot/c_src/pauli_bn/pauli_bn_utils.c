#include "pauli_bn_utils.h"


char xz_to_pstr(bool x, bool z){
    if (x||z){
        if (x && z){return ASCII_Y;}
        else if(x){return ASCII_X;}
        else{return ASCII_Z;}
    }
    else{return ASCII_I;}
}

void _ints_to_pstr(DTYPE nx, DTYPE nz, size_t type_size, char * buff){
    int bit_size = 8*type_size;
    unsigned int mask = 1 << (bit_size- 1);

    for(int i=0; i<bit_size; i++){
        buff[i] = (char)xz_to_pstr((bool)(nx&mask), (bool)(nz&mask));
        mask >>=1;
    }
    
}

bool commute_test(struct bn * nx1, struct bn * nz1, struct bn * nx2, struct bn * nz2)
{
    struct bn tmp1, tmp2;
    bignum_and(nx1, nz2, &tmp1);
    bignum_and(nz1, nx2, &tmp2);
    /* Naive implementation
    i = bignum_bit_count(&tmp1)&1;
    j = bignum_bit_count(&tmp2)&1;
    if(i==j){Py_RETURN_TRUE;}
    else{Py_RETURN_FALSE;}
    */

    // Using a simpler method 
    bignum_xor(&tmp1, &tmp2, &tmp2);
    return (bool)bignum_bit_count(&tmp2)^1;
}

size_t bignum_tuple_to_pstr(
    struct bn * nx, struct bn * nz, 
    size_t qubits,
    char * buff, size_t buff_size)
    {
        // Error bit <=31 work well, but not for >32.

    // Calculate access units
    int bit_unit = sizeof(DTYPE) * 8;
    int max_index_arr = (int)qubits/bit_unit;
    int empty_str_len = qubits%bit_unit;

    int j = 0, k=0;
    int i=BN_ARRAY_SIZE-max_index_arr;
    for(; i < BN_ARRAY_SIZE+1; i++)
    {
        j = BN_ARRAY_SIZE -i;
        _ints_to_pstr(nx->array[j], nz->array[j], sizeof(DTYPE), buff+k*sizeof(DTYPE));
        k++;
    }

    //int null_posiiton = (k)*bit_unit ;
    //(buff+null_posiiton)[0] = '\0';

    memmove(buff, buff + bit_unit - empty_str_len, qubits+1);
    (buff+qubits)[0]= '\0';
    size_t length = strlen(buff);
    return length;
}


/* 
Convert the x, z pauli code to  
Pauli string in ASCII code.
*/

//char _xz_to_pstr(bool x, bool z){
//    int i_term =  ASCII_I*((int)(!(x||z)));
//    int x_term = 0, z_term = 0;
//    int r = 0;
//
//    if(x)
//    {
//        x_term = ASCII_X; 
//        z_term = ASCII_X;
//        r++;
//
//    };
//    if(z)
//    {
//        z_term =  ASCII_Z;
//        r++;
//    };
//
//    return (char)(i_term+((x_term+z_term)/r));
//}
