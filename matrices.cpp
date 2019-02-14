#include "matrices.h"

#include <libhcs.h>


/**
 * Main function
 *
 * Must be called with the parameters n, p, and a seed, separated by spaces
 *
 */

void fmpz_mat_mul_modp(fmpz_mat_t a, fmpz_mat_t b, fmpz_mat_t c, int n,
    fmpz_t p) {
  fmpz_mat_mul(a, b, c);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      fmpz_mod(fmpz_mat_entry(a, i, j), fmpz_mat_entry(a, i, j), p);
    }
  }
}

void fmpz_mat_inverse_modp(fmpz_t det, fmpz_mat_t b, fmpz_mat_t a, int n,fmpz_t p)
{
    fmpz_invmod(det,det, p);

    fmpz_mat_scalar_mul_fmpz(b, a, det);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            fmpz_mod(fmpz_mat_entry(a, i, j), fmpz_mat_entry(a, i, j), p);
    }
  }

}

void print_random_matrices_with_adj(char *n_str, char *p_str, char *simulated,
    char *seed) {
  int n = atoi(n_str);
  int is_simulated_setup = atoi(simulated);
  cryptorand_t randstate;
  cryptorand_initseed(randstate, seed ? seed : "", NULL);

  fmpz_t modp;
  fmpz_init(modp);
  fmpz_set_str(modp, p_str, 10);

  fmpz_mat_t a;
  fmpz_mat_init(a, n, n);

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      fmpz_randm_crypto(fmpz_mat_entry(a, i, j), randstate, modp);
    }
  }

  fmpz_t det;
  fmpz_init(det);

  fmpz_mat_t adjugate;
  fmpz_mat_init(adjugate, n, n);

  fmpz_mat_t prod;
  fmpz_mat_init(prod, n, n);

  fmpz_mat_t check;
  fmpz_mat_init(check, n, n);

  if(is_simulated_setup) {
    /* set det and adj randomly */
    fmpz_randm_crypto(det, randstate, modp);

    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        fmpz_randm_crypto(fmpz_mat_entry(adjugate, i, j), randstate, modp);
      }
    }
  } else {
    fmpz_modp_matrix_det(det, a, n, modp);
    if (fmpz_is_zero(det)) {
      fprintf(stderr, "ERROR: Random matrix was not invertible.\n");
      goto exit_det;
    }

    fmpz_modp_matrix_adjugate(adjugate, a, n, modp);
    fmpz_mat_transpose(adjugate, adjugate);

    fmpz_mat_mul_modp(prod, a, adjugate, n, modp);

    /* check that the adjugate and determinant were computed correctly */
    fmpz_mat_one(check);
    fmpz_mat_scalar_mul_fmpz(check, check, det);

    int status = fmpz_mat_equal(prod, check);
    if (status == 0) {
      fprintf(stderr, "ERROR: Failed to produce the proper matrices.\n");
      goto exit;
    }
  }

  /* print the resulting values */
//  fmpz_fprint(stdout, det);
//  printf("\n");
//  fmpz_mat_fprint(stdout, a);
//  printf("\n");
//  fmpz_mat_transpose(adjugate, adjugate);
//  fmpz_mat_fprint(stdout, adjugate);
//  printf("\n");

  // print the matrix and its inverse

  fmpz_fprint(stdout, det);
  printf("\n");
  fmpz_mat_fprint(stdout, a);
  printf("\n\n");

  fmpz_mat_inverse_modp(det, adjugate, adjugate, n, modp);
  fmpz_mat_fprint(stdout, adjugate);
  printf("\n\n");

  fmpz_mat_mul_modp(prod, a, adjugate, n, modp);
  fmpz_mat_fprint(stdout, prod);
  printf("\n\n");

exit:
  fmpz_mat_clear(a);
  fmpz_mat_clear(prod);
  fmpz_mat_clear(check);

exit_det:
  fmpz_mat_clear(adjugate);
  fmpz_clear(det);

  cryptorand_clear(randstate);
}


void fmpz_modp_matrix_det(fmpz_t det, fmpz_mat_t a, int n, fmpz_t p) {
  assert(n >= 1);

  if(n == 1) {
    fmpz_set(det, fmpz_mat_entry(a, 0, 0));
    return;
  }

  if (n == 2) {
    fmpz_t tmp1;
    fmpz_init(tmp1);
    fmpz_mul(tmp1, fmpz_mat_entry(a,0,0), fmpz_mat_entry(a,1,1));
    fmpz_mod(tmp1, tmp1, p);
    fmpz_t tmp2;
    fmpz_init(tmp2);
    fmpz_mul(tmp2, fmpz_mat_entry(a,1,0), fmpz_mat_entry(a,0,1));
    fmpz_mod(tmp2, tmp2, p);
    fmpz_sub(det, tmp1, tmp2);
    fmpz_mod(det, det, p);
    fmpz_clear(tmp1);
    fmpz_clear(tmp2);
    return;
  }

  fmpz_mat_t m;
  fmpz_mat_init_set(m, a);

  fmpz_t tmp;
  fmpz_init(tmp);
  fmpz_t multfactor;
  fmpz_init(multfactor);

  int num_swaps = 0;

  for(int j = 0; j < n; j++) {
    for(int i = j+1; i < n; i++) {

      if(fmpz_is_zero(fmpz_mat_entry(m, j, j))) {
        // find first row that isn't a zero, and swap
        int h;
        for(h = j+1; h < n; h++) {
          if(!fmpz_is_zero(fmpz_mat_entry(m, h, j))) {
            // found the row
            break;
          }
        }

        if(h == n) {
          // matrix is not invertible
          fmpz_set_ui(det, 0);
          fmpz_clear(multfactor);
          fmpz_clear(tmp);
          fmpz_mat_clear(m);
          return;
        }

        // swap row h with row j
        for(int k = 0; k < n; k++) {
          fmpz_set(tmp, fmpz_mat_entry(m, h, k));
          fmpz_set(fmpz_mat_entry(m, h, k), fmpz_mat_entry(m, j, k));
          fmpz_set(fmpz_mat_entry(m, j, k), tmp);
        }

        num_swaps++;
      }

      fmpz_invmod(multfactor, fmpz_mat_entry(m, j, j), p);
      fmpz_mul(multfactor, multfactor, fmpz_mat_entry(m, i, j));
      fmpz_mod(multfactor, multfactor, p);

#pragma omp parallel for
      for(int k = j; k < n; k++) {
        fmpz_t tmp2;
        fmpz_init(tmp2);
        fmpz_mul(tmp2, fmpz_mat_entry(m, j, k), multfactor);
        fmpz_sub(fmpz_mat_entry(m, i, k), fmpz_mat_entry(m, i, k), tmp2);
        fmpz_mod(fmpz_mat_entry(m, i, k), fmpz_mat_entry(m, i, k), p);
        fmpz_clear(tmp2);
      }
    }
  }

  fmpz_clear(multfactor);
  fmpz_clear(tmp);

  fmpz_set_ui(det, 1);

  for(int j = 0; j < n; j++) {
    fmpz_mul(det, det, fmpz_mat_entry(m, j, j));
  }
  if(num_swaps % 2 == 1) {
    fmpz_neg(det, det);
  }
  fmpz_mod(det, det, p);
  fmpz_mat_clear(m);
}

void fmpz_modp_matrix_adjugate(fmpz_mat_t b, fmpz_mat_t a, int n, fmpz_t p) {
  if(n == 1) {
    fmpz_set_ui(fmpz_mat_entry(b, 0, 0), 1);
    return;
  }

  fmpz_t det;
  fmpz_init(det);

  fmpz_mat_t c;
  fmpz_mat_init(c, n-1, n-1);

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      /* Form the adjoint a_ij */
      for (int i_iter = 0, i1 = 0; i_iter < n; i_iter++, i1++) {
        if (i_iter == i) {
          i1--;
          continue;
        }
        for (int j_iter = 0, j1 = 0; j_iter < n; j_iter++, j1++) {
          if (j_iter == j) {
            j1--;
            continue;
          }
          fmpz_set(fmpz_mat_entry(c, i1, j1), fmpz_mat_entry(a, i_iter, j_iter));
        }
      }

      /* Calculate the determinant */
      fmpz_modp_matrix_det(det, c, n-1, p);

      /* Fill in the elements of the adjugate */
      if((i+j) % 2 == 1) {
        fmpz_negmod(det, det, p);
      }
      fmpz_mod(det, det, p);
      fmpz_set(fmpz_mat_entry(b, i, j), det);
    }
  }

  fmpz_clear(det);
  fmpz_mat_clear(c);
}



/* random select rop from [0, rop] such that gcd(rop, op)==1 */

void fmpz_random_in_mult_group(fmpz_t rop, flint_rand_t state, fmpz_t op)
{
    fmpz_t t1;
    fmpz_init(t1);

    do {
        random(rop, state, op);
        fmpz_gcd(t1, rop, op);
    } while (fmpz_cmp_ui(t1, 1) != 0);

    fmpz_clear(t1);
}

// M = M mod p
void fmpz_mat_mod_fmpz(fmpz_mat_t M, fmpz_t p)
{
    int i, j;

    for(i = 0; i < M->r; i++)
    {
        for(j = 0; j < M->c; j++)
        {
            fmpz_mod(fmpz_mat_entry(M, i, j), fmpz_mat_entry(M, i, j), p);
        }
    }

}

// compute v=v * M mod p, v is a vector and M is a squared matrix

void fmpz_vec_mat_mul_mod_fmpz(fmpz* v, fmpz_mat_t M, fmpz_t p)
{

    fmpz* v1=_fmpz_vec_init(M->c);

    for(long i=0; i<M->c; i++)
    {
        fmpz_zero(v1+i);

        for(long j=0; j<M->r; j++)
            fmpz_addmul(v1+i, v+j, fmpz_mat_entry(M, j, i));
    }

    _fmpz_vec_scalar_mod_fmpz(v1, v1, M->c, p);

    _fmpz_vec_set(v, v1, M->c);

    _fmpz_vec_clear(v1, M->c);

}


void fmpz_dot_mod_fmpz(fmpz_mat_t v1, fmpz_mat_t v2, fmpz_t res, fmpz_t p)
{
    if(v1->r!=v2->r || v1->c !=v2->c)
        perror("the dimension is un-matching");

    fmpz_zero(res);

    for (int i=0; i<v1->c; i++)
    {
        fmpz_addmul(res, fmpz_mat_entry(v1,0,i), fmpz_mat_entry(v2, 0, i));
        fmpz_mod(res, res, p);
    }

}


void fmpz_dot(fmpz_mat_t v1, fmpz_mat_t v2, fmpz_t res)
{
    if(v1->r!=v2->r || v1->c !=v2->c)
        perror("the dimension is un-matching");

    fmpz_zero(res);

    for (int i=0; i<v1->c; i++)
    {
        fmpz_addmul(res, fmpz_mat_entry(v1,0,i), fmpz_mat_entry(v2, 0, i));
    }

}

// randomly generate a within Z_N

void random(fmpz_t a, flint_rand_t state, fmpz_t N)
{
    fmpz_randm(a, state, N);
}


void fmpz_tri_inv(int n, flint_rand_t state, fmpz_t modp, fmpz_mat_t a, fmpz_mat_t inv)
{
    int i;

    fmpz_mat_zero(a);       // initialize a to all zeros
    fmpz_mat_zero(inv);     // initialize a to all zeros


    #pragma omp parallel for private(i)
    for(i=1;i<=n;i++)
    {
        int k;

        for(k=1;k<=i;k++)
        {
            if (k==i)
                fmpz_random_in_mult_group(fmpz_mat_entry(a, i-1, k-1),  state, modp); // from Z_N*
            else
                random(fmpz_mat_entry(a, i-1, k-1), state, modp);  // from Z_N

        }
    }

    // compute inv column by column

    #pragma omp parallel for private(i)
    for(i=1;i<=n;i++)
    {
        int k;

        fmpz_invmod(fmpz_mat_entry(inv, i-1, i-1), fmpz_mat_entry(a, i-1, i-1), modp);

        // consider inv is also a triangle matrix, for the i-column of Inv, its first i-1 elements are zero, and the i element is computed as above
        // we only consider to compute the remaining elements Inv[i+1][i],..., Inv[n][i]

        fmpz_t sum;
        fmpz_t current;

        fmpz_init(sum);
        fmpz_init(current);

        // compute the value at [k,i] of Inv

        for (k=i+1;k<=n;k++)
        {
            int l;

            fmpz_set_ui(sum, 0);
            fmpz_set_ui(current, 0);

            // compute a[k][i].Inv[i][i] + ... +  a[k][k-1].Inv[k-1][i]

            for (l=i;l<k;l++)
            {
                fmpz_mul(current, fmpz_mat_entry(a, k-1, l-1),fmpz_mat_entry(inv, l-1, i-1));

                fmpz_mod(current, current, modp);

                fmpz_add(sum, sum, current);

                fmpz_mod(sum, sum, modp);

            }

            // compute Inv[k][i] by sum + a[k][k]*Inv[k][i]=0 mod N

            fmpz_sub(sum, modp, sum);

            fmpz_invmod(current, fmpz_mat_entry(a, k-1, k-1), modp);

            fmpz_mul(fmpz_mat_entry(inv, k-1, i-1), sum, current);

            fmpz_mod(fmpz_mat_entry(inv, k-1, i-1), fmpz_mat_entry(inv, k-1, i-1), modp);
        }

        fmpz_clear(sum);
        fmpz_clear(current);

    }

}

// random generate a inveretible matrix over Z_N and compute its inverse

void fmpz_mat_inv(int n, flint_rand_t state, fmpz_t modp, fmpz_mat_t a, fmpz_mat_t inv)
{

    // we compute a with LU decomposition, i.e., a=L*U and a^{-1}=U^{-1}*L^{-1}, where L is a lower triangle matrix and U is a upper triangle matrix

    fmpz_mat_t L;
    fmpz_mat_t invL;
    fmpz_mat_t U;
    fmpz_mat_t invU;

    fmpz_mat_init(L, n, n);
    fmpz_mat_init(invL, n, n);
    fmpz_mat_init(U, n, n);
    fmpz_mat_init(invU, n, n);


    // random generate L, U, InvL, InvU

    fmpz_tri_inv(n, state, modp,  L, invL);
    fmpz_tri_inv(n, state, modp,  U, invU);

    // compute M=L* (U^T); M^{-1}= (U^T)^{-1} * L^{-1}

    fmpz_mat_transpose(U, U);

    fmpz_mat_mul_modp(a, L, U, n, modp);

    fmpz_mat_transpose(invU, invU);

    fmpz_mat_mul_modp(inv, invU, invL, n, modp);

    // free the space

    fmpz_mat_clear(L);
    fmpz_mat_clear(invL);
    fmpz_mat_clear(U);
    fmpz_mat_clear(invU);

}


void fmpz_norm(fmpz_mat_t v, fmpz_t norm)
{
    fmpz_dot(v, v, norm);
}

void fmpz_mat_readFile(fmpz_mat_t mat, string filepath)
{
    string line;

    long tempValue;

    ifstream fin(filepath.c_str(), ifstream::in);

    long rows=0, cols=0;

    while (getline(fin, line))
    {
        if (line.length() == 0)
            continue;

        rows++;

        // check the column size

        if(rows==1)
        {
            stringstream test(line);

            while (test >> tempValue)
            {
                if (test.peek() == ' ')
                {
                    test.ignore();
                    cols++;
                }
            }
        }
    }

    fin.close();

    // read the file into the matrix

    fmpz_mat_init(mat, rows, cols);

    fin.open(filepath.c_str(), ifstream::in);

    rows=0;

    while (getline(fin, line))
    {
        if (line.length() == 0)
            continue;

        rows++;

        stringstream test(line);

        cols=0;

        while (test >> tempValue)
        {
            fmpz_set_si(fmpz_mat_entry(mat, rows-1, cols), tempValue);

            if (test.peek() == ' ')
            {
                test.ignore();
                cols++;
            }
        }

    }

    fin.close();

}


void fmpz_mat_randomperm(fmpz_mat_t mat, flint_rand_t state)
{
    slong i, j;
    slong * rows;

    rows = _perm_init(mat->r);
    _perm_randtest(rows, mat->r, state);

    fmpz_mat_t temp;
    fmpz_mat_init(temp, mat->r, mat->c);
    fmpz_mat_zero(temp);

    for (i = 0; i < mat->r; i++)
    {
        for(j=0; j<mat->c; j++)
        {
            fmpz_set(fmpz_mat_entry(temp, rows[i], j), fmpz_mat_entry(mat, i, j));
        }
    }

    fmpz_mat_swap(mat, temp);

    _perm_clear(rows);

    fmpz_mat_clear(temp);
}
