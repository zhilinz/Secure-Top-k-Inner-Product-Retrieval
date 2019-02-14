#ifndef _GEN_MATRICES_H_
#define _GEN_MATRICES_H_

//#include <getopt.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <unistd.h>
//#include <gmp.h>
#include <assert.h>
//#include <string>
//#include <omp.h>
#include <flint/fmpz_mat.h>
#include <flint/perm.h>

#include "cryptorand.h"
#include "base.h"

using namespace std;

/* random select rop from [0, rop] such that gcd(rop, op)==1 */

void fmpz_random_in_mult_group(fmpz_t rop, flint_rand_t state, fmpz_t op);

// M = M mod p
void fmpz_mat_mod_fmpz(fmpz_mat_t M, fmpz_t p);

// compute v=v * M mod p, v is a vector and M is a matrix
void fmpz_vec_mat_mul_mod_fmpz(fmpz* v, fmpz_mat_t M, fmpz_t p);

// set res as the dot product of v1 and v2 mod p
void fmpz_dot_mod_fmpz(fmpz_mat_t v1, fmpz_mat_t v2, fmpz_t res, fmpz_t p);

void print_random_matrices_with_adj(char *n_str, char *p_str, char *simulated,
    char *seed);

/**
 * Wrapper around fmpz_mat_mul, which simply multiplies two matrices b and c and
 * stores the result in a. Then, each entry of a is modded by p.
 *
 * All matrices are of dimension n x n.
 *
 * @param a The product of the two matrices
 * @param b The first matrix
 * @param c The second matrix
 * @param n The dimension of all matrices
 * @param p The modulus
 * @return void
 */
void fmpz_mat_mul_modp(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, int n,
    fmpz_t p);

/**
 * Computes the determinant of the matrix a, of dimension n x n, mod p. The
 * result is stored in det.
 *
 * This function uses Gaussian elimination to compute the determinant, and has
 * complexity O(n^3).
 *
 * @param det The determinant, to be stored
 * @param a The matrix to compute the determinant of
 * @param n The dimension of the matrix
 * @param p The modulus
 * @return void
 */
void fmpz_modp_matrix_det(fmpz_t det, fmpz_mat_t a, int n, fmpz_t p);

/**
 * Computes the adjugate of the matrix a, of dimension n x n, mod p. The result
 * is stored in b.
 *
 * @param det The adjugate of the matrix, to be stored
 * @param a The matrix to compute the determinant of
 * @param n The dimension of the matrix
 * @param p The modulus
 * @return void
 */
void fmpz_modp_matrix_adjugate(fmpz_mat_t b, fmpz_mat_t a, int n, fmpz_t p);


// a <- inv(x) mod p, b is the adj matrix of x
void fmpz_mat_inverse_modp(fmpz_t det, fmpz_mat_t b, fmpz_mat_t a, int n,fmpz_t p);


// randomly generate a within Z_N

void random(fmpz_t a, flint_rand_t state, fmpz_t N);

// compute dot(v1, v2) over Z
void fmpz_dot(fmpz_mat_t v1, fmpz_mat_t v2, fmpz_t res);

// random generate a inveretible triangle matrix over Z_N and compute its inverse

void fmpz_tri_inv(int n, flint_rand_t state, fmpz_t modp, fmpz_mat_t a, fmpz_mat_t inv);

// random generate a inveretible matrix over Z_N and compute its inverse

void fmpz_mat_inv(int n, flint_rand_t state, fmpz_t modp, fmpz_mat_t a, fmpz_mat_t inv);

// compute the eculidian norm of a vector, norm =|v|^2

void fmpz_norm(fmpz_mat_t v, fmpz_t norm);

// read a file into a matrix

void fmpz_mat_readFile(fmpz_mat_t mat, string filepath);

// random permute the rows of a matrix

void fmpz_mat_randomperm(fmpz_mat_t mat, flint_rand_t state);



#endif /* _GEN_MATRICES_H_ */
