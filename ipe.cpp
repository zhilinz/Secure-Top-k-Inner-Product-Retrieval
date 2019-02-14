#include "ipe.h"
#include "monitor.h"

/***********************************************************************
 *               initialize pp and msk
 ***********************************************************************/

void initialize(PublicParameter *pp, MasterSecretKey *msk, int dim)
{
    fmpz_init(pp->N);
    fmpz_init(pp->N2);
    fmpz_init(pp->g);
    pp->dim=dim;

    msk->pk  = djcs_init_public_key();
    msk->sk  = djcs_init_private_key();

    fmpz_init(msk->h);
    fmpz_init(msk->lambda);

    fmpz_mat_init(msk->M1, dim, dim);
    fmpz_mat_init(msk->invM1, dim, dim);

    fmpz_mat_init(msk->M2, dim, dim);
    fmpz_mat_init(msk->invM2, dim, dim);


    msk->S  = _fmpz_vec_init(2*dim);
    msk->hS = _fmpz_vec_init(2*dim);

}

/***********************************************************************
 *               generate pp and msk
 ***********************************************************************/

void setup(PublicParameter *pp, MasterSecretKey *msk, int bits)
{
    /***********************1. generate msk ****************/

    // generate msk.pk, msk.sk

    hcs_random *hr = hcs_init_random();

    hcs_reseed_random(hr);

    djcs_generate_key_pair(msk->pk, msk->sk, hr, 1, bits);

    hcs_free_random(hr);

    char *N = mpz_get_str(NULL,10, msk->pk->n[0]);
    fmpz_set_str(pp->N, N, 10);
    char *NSquare = mpz_get_str(NULL,10, msk->pk->n[1]);
    fmpz_set_str(pp->N2, NSquare, 10);
    fmpz_add_ui(pp->g, pp->N, 1);


    // generate msk.M, msk.invM

    flint_rand_t state;

    flint_randinit(state);

    // this function can quickly generate M1, M2 and invM1, invM2 with arbitrary large dim and N, while M1 and M2 is a diagonal matrix

    fmpz_mat_inv(pp->dim, state, pp->N, msk->M1, msk->invM1);

    fmpz_mat_transpose(msk->invM1, msk->invM1);   // invM1 ={invM1}^T

    fmpz_mat_inv(pp->dim, state, pp->N, msk->M2, msk->invM2);

    fmpz_mat_transpose(msk->invM2, msk->invM2);   // invM2 ={invM2}^T


//    fmpz_mat_swap(msk->M1, msk->invM1);         // it makes invM1 and invM2 to be a full matrix
//    fmpz_mat_swap(msk->M2, msk->invM2);


    // h=h0^{2N}

    fmpz_random_in_mult_group(msk->h, state, pp->N2);  // h \in Z_{N^2}^*

    fmpz_powm(msk->h, msk->h, pp->N, pp->N2);

    fmpz_powm_ui(msk->h, msk->h, 2, pp->N2);


    // generate S that is randomly selected over Z_N, each s_i is a random value in [1, ord(G)]

    fmpz_t ord;
    fmpz_init(ord);

    char *ordString = mpz_get_str(NULL,10,msk->sk->d);

    fmpz_set_str(ord, ordString, 10);               // currently, ord = lambda(N)=lcm(p-1. q-1)

    fmpz_set(msk->lambda, ord);

    fmpz_mul(ord, ord, pp->N);

    fmpz_divexact_ui(ord, ord, 2);                  // currently, ord =N*lambda(N)/2

    fmpz_sub_ui(ord, ord, 1);

    // s \in [1, ord(G)]

    for(int i = 0; i < 2*pp->dim; i++)
    {
        random(msk->S+i, state, ord);

        fmpz_add_ui(ord, ord, 1);

        fmpz_powm(msk->hS+i, msk->h, msk->S+i, pp->N2);
    }


    //////////////////////////////////////////////////////////////

    void (*freefunc)(void *, size_t);
    mp_get_memory_functions (NULL, NULL, &freefunc);

    freefunc(N, strlen(N) + 1);
    freefunc(NSquare, strlen(NSquare) + 1);
    freefunc(ordString, strlen(ordString) + 1);

    fmpz_clear(ord);

    flint_randclear(state);

}


// free the space

void clearup (PublicParameter *pp, MasterSecretKey *msk)
{
    fmpz_clear(pp->N);
    fmpz_clear(pp->N2);
    fmpz_clear(pp->g);

    djcs_free_public_key(msk->pk);
    djcs_free_private_key(msk->sk);

    fmpz_clear(msk->h);
    fmpz_clear(msk->lambda);

    fmpz_mat_clear(msk->M1);
    fmpz_mat_clear(msk->invM1);
    fmpz_mat_clear(msk->M2);
    fmpz_mat_clear(msk->invM2);

    _fmpz_vec_clear(msk->S, 2*pp->dim);
    _fmpz_vec_clear(msk->hS, 2*pp->dim);

}


void random_split(fmpz * q, fmpz * qa, fmpz * qb, PublicParameter *pp, flint_rand_t state)
{

    long i;

    for (i=0; i<pp->dim; i++)
    {
        random(qa+i, state, pp->N);

        fmpz_sub (qb+i, q+i, qa+i);

        fmpz_mod(qb+i, qb+i, pp->N);
    }

}

// given the plaintext query q, generate the ciphertext query K

void keyGen(fmpz * q, fmpz * K, PublicParameter *pp, MasterSecretKey *msk)
{
    fmpz * qa=_fmpz_vec_init(pp->dim);
    fmpz * qb=_fmpz_vec_init(pp->dim);
    fmpz * qc=_fmpz_vec_init(2*pp->dim);  // qc=(qa, qb)

    // random split a into qa and qb

    flint_rand_t state;

    flint_randinit(state);

    random_split(q, qa, qb, pp, state);

    flint_randclear(state);

    // compute qa * inv1, qb*inv2

    fmpz_vec_mat_mul_mod_fmpz(qa, msk->invM1, pp->N);

    fmpz_vec_mat_mul_mod_fmpz(qb, msk->invM2, pp->N);


    _fmpz_vec_set(qc, qa, pp->dim);

    for(long i=pp->dim; i<2*pp->dim; i++)
        fmpz_set(qc+i, qb+i-pp->dim);

    // generate the final ciphertext K

    _fmpz_vec_dot (K, qc , msk->S, 2*pp->dim);

    fmpz_mod(K, K, msk->lambda);                // compute b=dot(q1, S) mod lamada

    // K[i]=q'[i-1], i=1, ..., 2*\ell

    for(long i=1; i<2*pp->dim+1; i++)
        fmpz_set(K+i, qc+i-1);


    // free the space

    _fmpz_vec_clear(qa, pp->dim);
    _fmpz_vec_clear(qb, pp->dim);
    _fmpz_vec_clear(qc, 2*pp->dim);

}


// given a plaintext data p, output the ciphertext C

void encrypt(fmpz * p, fmpz * C, PublicParameter *pp, MasterSecretKey *msk)
{
    flint_rand_t state;

    flint_randinit(state);

    fmpz * pa=_fmpz_vec_init(pp->dim);
    fmpz * pb=_fmpz_vec_init(pp->dim);
    fmpz * pc=_fmpz_vec_init(2*pp->dim);

    // compute pa=p*M1 and pb=p*M2, pc=(pa, pb)

     _fmpz_vec_set(pa, p, pp->dim);

     _fmpz_vec_set(pb, p, pp->dim);

    fmpz_vec_mat_mul_mod_fmpz(pa, msk->M1, pp->N);

    fmpz_vec_mat_mul_mod_fmpz(pb, msk->M2, pp->N);

    // concat pc=(pa, pb)

    _fmpz_vec_set(pc, pa, pp->dim);

    for(long i=pp->dim; i<2*pp->dim; i++)
        fmpz_set(pc+i, pb+i-pp->dim);

    // generate the final ciphertext C

    fmpz_t a;
    fmpz_t s;

    fmpz_init(a);
    fmpz_init(s);

    random(a, state, pp->N);

    fmpz_powm(C, msk->h, a, pp->N2);   // C0=h^a

    int i;

    for(i=1;i<2*pp->dim+1; i++)
    {
        fmpz_mul(C+i, pc+i-1, pp->N);

        fmpz_add_ui(C+i, C+i, 1);           // 1+pi*N

        fmpz_mod(C+i, C+i, pp->N2);

        fmpz_powm(s, msk->hS+i-1, a, pp->N2);

        fmpz_mul(C+i, C+i, s);

        fmpz_mod(C+i, C+i,pp->N2);

    }

    _fmpz_vec_clear(pa, pp->dim);
    _fmpz_vec_clear(pb, pp->dim);
    _fmpz_vec_clear(pc, 2*pp->dim);

    fmpz_clear(a);
    fmpz_clear(s);

    flint_randclear(state);

}



void decrypt(fmpz * C, fmpz * K, fmpz_t res, PublicParameter *pp)
{
    fmpz_t temp;

    fmpz_init(temp);

    fmpz_powm(temp, C, K, pp->N2);     // obtain h^{b*(a+b)}

    fmpz_invmod(res, temp, pp->N2);    // obtain h^{-ab mod lamada}


    for(int i=1; i<2*pp->dim+1; i++)
    {

        fmpz_powm(temp, C+i, K+i, pp->N2);
        fmpz_mul(res, res, temp);
        fmpz_mod(res, res, pp->N2);
    }

    fmpz_sub_ui(res, res, 1);

    fmpz_mod(res, res, pp->N2);

    fmpz_divexact(res, res, pp->N);

    fmpz_clear(temp);

}


void decrypt2(fmpz * C, PublicParameter *pp,  MasterSecretKey *msk)
{
    int dim=pp->dim;

    fmpz_t temp;

    fmpz_init(temp);

    fmpz_mat_t invM2;

    fmpz_mat_init(invM2, dim, dim);

    for(int i=1; i<=dim; i++)
    {
        fmpz_powm(temp, C, msk->S+dim+i-1, pp->N2);     // obtain C0^{s_i}

        fmpz_invmod(temp, temp, pp->N2);                // obtain C0^{-s_i}

        fmpz_mul(C+dim+i, C+dim+i, temp);

        fmpz_sub_ui(C+dim+i, C+dim+i, 1);

        fmpz_mod(C+dim+i, C+dim+i, pp->N2);

        fmpz_divexact(C+dim+i, C+dim+i, pp->N);

    }

    fmpz_mat_transpose(invM2, msk->invM2);

    fmpz_vec_mat_mul_mod_fmpz(C+dim+1, invM2, pp->N);

    fmpz_clear(temp);

    fmpz_mat_clear(invM2);
}
