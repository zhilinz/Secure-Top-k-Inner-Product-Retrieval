#ifndef _IPE_H_INCLUDED
#define _IPE_H_INCLUDED

//#include <getopt.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <unistd.h>
//#include <gmp.h>
#include <assert.h>
#include <omp.h>
#include <libhcs.h>


#include <iostream>
#include <fstream>
#include <string>

#include "matrices.h"
#include "monitor.h"


using namespace std;


// public parameters
typedef struct
{
    fmpz_t N;
    fmpz_t N2;  //n^2
    fmpz_t g;   // g=1+N
    int  dim;

}PublicParameter;

// master secret key
typedef struct
{
    djcs_public_key *pk;
    djcs_private_key *sk;

    fmpz_t h;
    fmpz_t lambda;

    fmpz_mat_t  M1;
    fmpz_mat_t  invM1;

    fmpz_mat_t  M2;
    fmpz_mat_t  invM2;

    fmpz *  S;
    fmpz *  hS; // h^{s_i}

} MasterSecretKey;


/***********************************************************************
 *               initialize pp and msk
 ***********************************************************************/

void initialize(PublicParameter *pp, MasterSecretKey *msk, int dim);

/***********************************************************************
 *               generate pp and msk
 ***********************************************************************/

void setup(PublicParameter *pp, MasterSecretKey *msk, int bits);

// free the space

void clearup (PublicParameter *pp, MasterSecretKey *msk);

// random split q into qa and qb such that q=qa+qb mod N

void random_split(fmpz * q, fmpz * qa, fmpz * qb, PublicParameter *pp, flint_rand_t state);

// given the plaintext query q, generate the ciphertext query K

void keyGen(fmpz * q, fmpz * K, PublicParameter *pp, MasterSecretKey *msk);

// given a plaintext data p, output the ciphertext C

void encrypt(fmpz * p, fmpz * C, PublicParameter *pp, MasterSecretKey *msk);


void decrypt(fmpz * C, fmpz * K, fmpz_t res, PublicParameter *pp);


// given the last \ell elements [p], reveal the original data vector p

void decrypt2(fmpz * C, PublicParameter *pp,  MasterSecretKey *msk);


#endif // IPE_H_INCLUDED
