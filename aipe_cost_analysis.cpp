#include <iostream>
#include <fstream>
#include <sys/time.h>


#include "skip.h"

using namespace std;

int main()
{
    int n=10;
    int bits=1024;
    int msgbits=30;

    PublicParameter pp;
    MasterSecretKey msk;

    srand(time(NULL));

    flint_rand_t state;
    flint_randinit(state);

    Monitor tt;

    ofstream fout;


    fmpz_t res;
    fmpz_init(res);

    for(int i=1; i<=5; i++)
    {
        int dim=i*50;

        cout<<"dim="<<dim<<"......."<<endl;

        initialize(&pp, &msk, dim);
        setup(&pp, &msk, bits);

//        // test flint
//
//        fmpz_t a, b, c;
//
//        fmpz_init(a);
//        fmpz_init(b);
//        fmpz_init(c);
//
//        fmpz_set_ui(a, 11);
//        fmpz_set_ui(b, 5);
//
//        fmpz_pow_ui(a, a, 1023);
//        fmpz_pow_ui(b, b, 2047);
//
//
//        cout<<"size b="<<fmpz_sizeinbase(b, 2)<<endl;
//        cout<<"size a="<<fmpz_sizeinbase(a, 2)<<endl;
//
//        tt.start();
//        fmpz_powm(c, b, a, pp.N2);
//        tt.stop();
//
//        fmpz_print(pp.N2);
//
//        cout<<"flint Exp:"<<tt.getElapsedTime()<<endl;
//
//        fmpz_clear(a);
//        fmpz_clear(b);
//        fmpz_clear(c);



        fmpz_mat_t P;
        fmpz_mat_t EncP;
        fmpz_mat_t Q;
        fmpz_mat_t EncQ;

        fmpz_mat_init(P, n, dim);
        fmpz_mat_init(EncP, n, 2*dim+1);

        fmpz_mat_init(Q, n, dim);
        fmpz_mat_init(EncQ, n, 2*dim+1);

        // random generate P and Q

        for(int j=1; j<=n; j++)
        {
            _fmpz_vec_randtest_unsigned(P->rows[j-1], state, dim, msgbits);
            _fmpz_vec_randtest_unsigned(Q->rows[j-1], state, dim, msgbits);
        }

        // test KeyGen

        tt.start();

        for(int j=1; j<=n; j++)
            keyGen(Q->rows[j-1], EncQ->rows[j-1], &pp, &msk);

        tt.stop();

        fout.open("data/analysis.txt", ofstream::out | ofstream::app);
        fout << fixed<<setprecision(5)<<tt.getElapsedTime()/n<< " ";
        fout.close();


        // test Encrypt

        tt.start();

        for(int j=1; j<=n; j++)
            encrypt(P->rows[j-1], EncP->rows[j-1], &pp, &msk);

        tt.stop();

        fout.open("data/analysis.txt", ofstream::out | ofstream::app);
        fout << fixed<<setprecision(5)<<tt.getElapsedTime()/n<< " ";
        fout.close();


        // test Decrypt

        tt.start();

        for(int j=1; j<=n; j++)
        {
            for (int k=1; k<=n; k++)
                decrypt(EncP->rows[j-1], EncQ->rows[k-1], res, &pp);
        }

        tt.stop();

        fout.open("data/analysis.txt", ofstream::out | ofstream::app);
        fout << fixed<<setprecision(5)<< tt.getElapsedTime()/pow(n,2)<<endl;
        fout.close();


        // free the memory

        clearup(&pp, &msk);
        fmpz_mat_clear(P);
        fmpz_mat_clear(EncP);
        fmpz_mat_clear(Q);
        fmpz_mat_clear(EncQ);

    }

    return 0;
}
