#ifndef PACK_H_INCLUDED
#define PACK_H_INCLUDED

#include <math.h>
#include <vector>
#include <queue>
#include <iomanip>

#include "skip.h"


// the pure IP packing algorithm without sequential scan

class IP_Packing
{
private:

    vector<vector<IPTuple> > results;            // store the top-k answers of each query
    Monitor tt;                                 // use it to track the running time of each stage

    long u;
    long d;         // # of data vectors packed in one auxiliary vector
    long dim;       // dim of data and query vectors

    PublicParameter ppl;
    MasterSecretKey mskl;

    fmpz_mat_t EncP;           // encrypted auxiliary vectors
    fmpz_mat_t EncQ;           // encrypted queries

    fmpz* U;                   // [1, u, u^2, ...,u^{d-1}]

    void pack(fmpz_mat_t bucket, fmpz* auxilary);

    void unpack(fmpz_t ip, long j, fmpz_t ip1);

    void unpack2(fmpz* auxilary, long j, fmpz* p);     // restore the j-th data vector in auxiliary of size len

public:

    void initialization(int bits, unsigned long maxIP, fmpz_mat_t P);

    void phase1(fmpz_mat_t Q, int queryNum=1);

    void phase2(int k);

    void phase3();

    void printTopK();                                   // print the results

    void clearResult();

    ~IP_Packing();
};


/***************************************************************************/

// auxilary = U * bucket

void IP_Packing::pack(fmpz_mat_t bucket, fmpz * auxilary)
{
    for(long i=0; i<bucket->c; i++)
    {
        fmpz_zero(auxilary+i);

        for(long j=0; j< bucket->r; j++)
            fmpz_addmul(auxilary+i, U+j, fmpz_mat_entry(bucket, j, i));
    }

}

void IP_Packing::unpack(fmpz_t ip, long j, fmpz_t ip1)
{
    fmpz_t v1, v2;
    fmpz_init(v1);
    fmpz_init(v2);

    // v1= ip % u^j

    if(j< d)
        fmpz_mod(v1, ip, U+j);
    else
        fmpz_set(v1, ip);

    // v2= ip % u^{j-1}

    fmpz_mod(v2, ip, U+j-1);

    // ip1=(v1-v2)/u^{j-1} % u

    fmpz_sub(ip1, v1, v2);

    fmpz_divexact(ip1, ip1, U+j-1);

    fmpz_mod(ip1, ip1, U+1);

    // free the memory

    fmpz_clear(v1);
    fmpz_clear(v2);


}

void IP_Packing::unpack2(fmpz* auxilary, long j, fmpz* p)
{
    fmpz* u1 = _fmpz_vec_init(dim);
    fmpz* u2 = _fmpz_vec_init(dim);

    // u1= auxilary % u^j

    if(j< d)
        _fmpz_vec_scalar_mod_fmpz(u1, auxilary, dim, U+j);
    else
        _fmpz_vec_set(u1, auxilary, dim);

    // u2= auxilary % u^{j-1}

    _fmpz_vec_scalar_mod_fmpz(u2, auxilary, dim, U+j-1);

    // p=(u1-u2)/u^{j-1} % u

    _fmpz_vec_sub(p, u1, u2, dim);

    _fmpz_vec_scalar_divexact_fmpz(p, p, dim, U+j-1);

    _fmpz_vec_scalar_mod_fmpz(p, p, dim, U+1);


    // free the memory

    _fmpz_vec_clear(u1, dim);

    _fmpz_vec_clear(u2, dim);

}


/***************************************************************************/

void IP_Packing::initialization(int bits, unsigned long maxIP, fmpz_mat_t P)
{

    /**********************************************************
     *                      generate key and parameters
     *********************************************************/

    initialize(&ppl, &mskl, P->c);
    setup(&ppl, &mskl, bits);

    u=maxIP+1;

    d=floor(fmpz_dlog(ppl.N)/log(u))-1;

    if(d<1)
    {
         perror("d is too small\n");
         exit(0);
    }

    dim=P->c;

    // pre-compute U

    U = _fmpz_vec_init(d);

    fmpz_set_ui(U, 1);

    for(long i=1; i<d; i++)
        fmpz_mul_si(U+i, U+i-1, u);


    /**********************************************************
     *                    packing and encryption
     *********************************************************/
    tt.start();

    long rows, cols;

    rows=P->r;
    cols=P->c;

    long n1=ceil(double(rows)/d);

    fmpz_mat_init(EncP, n1, 2*cols+1);

    fmpz_mat_t bucket;
    fmpz * auxilary;

    fmpz_mat_init(bucket, d, cols);
    auxilary   = _fmpz_vec_init(cols);

    fmpz_mat_zero(bucket);

    flint_rand_t state;
    flint_randinit(state);

    // iterate on each bucket

    for (long i=0; i<n1; i++)
    {
        // fill the i-th bucket

        if(i==n1-1 && d*n1>rows)
        {
            // the last bucket might be filled with d*n1-row dummy all-zero vectors

            for(long j=0; j<rows-i*d; j++)
                for(long k=0; k<P->c; k++)
                    fmpz_set(fmpz_mat_entry(bucket, j, k), fmpz_mat_entry(P, i*d+j, k));

            for(long j=rows; j<d*n1; j++)
                for(long k=0; k<P->c; k++)
                    fmpz_set_ui(fmpz_mat_entry(bucket, j-i*d, k), 0);
        }
        else
        {
            for(long j=0; j<d; j++)
                for(long k=0; k<P->c; k++)
                    fmpz_set(fmpz_mat_entry(bucket, j, k), fmpz_mat_entry(P, i*d+j, k));
        }


        // random permute all vectors in bucket

        fmpz_mat_randomperm(bucket, state);

        // pack d vectors in bucket into an auxiliary vector

        pack(bucket, auxilary);


        // encrypt [a]

        encrypt(auxilary, EncP->rows[i], &ppl, &mskl);

    }

    tt.stop();

    // collect the the communication cost

    FILE *pfile;

    pfile=fopen("data/EncP.txt", "w+");

    fmpz_mat_fprint_pretty (pfile, EncP);

    fclose(pfile);

    ofstream fout;
    fout.open("data/log.txt", ofstream::out | ofstream::app);
    fout<<endl;
    fout <<"IP-Packing: d="<<d <<" ";
    fout <<"init time="<<tt.getElapsedTime() << " ";
    fout <<"outsourcing data size="<<evaluateFileSize("data/EncP.txt") << endl;
    fout.close();

    remove("data/EncP.txt");

    // free the memory

    fmpz_mat_clear(bucket);
    _fmpz_vec_clear(auxilary, dim);

    flint_randclear(state);
}


void IP_Packing::phase1(fmpz_mat_t Q, int queryNum)
{
    if(Q->r<queryNum)
    {
        perror("query number exceeds\n");
        exit(0);
    }

    tt.start();

    fmpz_mat_init(EncQ, queryNum, 2*Q->c+1);

    for(long i=0; i<queryNum; i++)
    {
        keyGen(Q->rows[i], EncQ->rows[i], &ppl, &mskl);
    }


    tt.stop();

    ofstream fout;
    fout.open("data/log.txt", ofstream::out | ofstream::app);
    fout << fixed<<setprecision(5)<<tt.getElapsedTime()/EncQ->r << " ";
    fout.close();
}


void IP_Packing::phase2(int k)
{
    tt.start();

    fmpz_t ip, ip1;

    fmpz_init(ip);

    fmpz_init(ip1);

    for(long i=0; i<EncQ->r; i++)
    {
        priority_queue<IPTuple, vector<IPTuple>, comparator > topK;

        // push k empty tuples into the heap

        for(int j=0; j<k; j++)
        {
            IPTuple *t= new IPTuple;

            fmpz_init(t->ip);

            t->auxilary=_fmpz_vec_init(2*dim+1);

            fmpz_set_ui(t->ip, 0);

            topK.push(*t);
        }

        // linear scan

        for(long j=0; j<EncP->r; j++)
        {

            decrypt(EncP->rows[j], EncQ->rows[i], ip, &ppl);                     // compute q^T a

            IPTuple result;

            for(long l=1; l<=d; l++)
            {
                unpack(ip, l, ip1);

                if(fmpz_cmp(ip1, topK.top().ip) >0)
                {
                    result=topK.top();

                    topK.pop();

                    fmpz_set(result.ip, ip1);                                  // set ip

                    _fmpz_vec_set(result.auxilary, EncP->rows[j], 2*dim+1);

                    result.index=l;

                    topK.push(result);
                }

            }

        }

        // store the results into answer sets

        vector<IPTuple> answers;

        IPTuple result;

        for(int j=0; j<k; j++)
        {
            result=topK.top();

            topK.pop();

            answers.push_back(result);

        }

        results.push_back(answers);

    }

    fmpz_clear(ip);

    fmpz_clear(ip1);

    tt.stop();

    ofstream fout;
    fout.open("data/log.txt", ofstream::out | ofstream::app);
    fout << fixed<<setprecision(5)<<tt.getElapsedTime()/EncQ->r << " ";
    fout << fixed<<setprecision(3)<<0.0<< " ";
    fout.close();
}


void IP_Packing::phase3()
{
    tt.start();

    for (ulong i=0; i<results.size(); i++)
    {
        // deal with the answer set of i-th query

        for(ulong j=0; j<results[i].size(); j++)
        {
            decrypt2(results[i][j].auxilary, &ppl, &mskl);    // the plain auxilary vector is stored in the last dim items

            unpack2(results[i][j].auxilary+dim+1, results[i][j].index, results[i][j].auxilary+dim+1);
        }

    }

    tt.stop();


    ofstream fout;
    fout.open("data/log.txt", ofstream::out | ofstream::app);
    fout << fixed <<setprecision(5)<< tt.getElapsedTime()/EncQ->r <<endl;
    fout.close();

}

void IP_Packing::printTopK()
{
    for (ulong i=0; i<results.size(); i++)
    {
        // deal with the answer set of i-th query

        for(ulong j=0; j<results[i].size(); j++)
        {
            fmpz_print(results[i][j].ip);

            cout<<" ";

            if(j==results[i].size()-1)
                cout<<endl;

//            fmpz_print(results[i][j].ip);
//
//            cout<<"->";
//
//            _fmpz_vec_print(results[i][j].auxilary+dim+1, dim);
//
//            cout<<endl;
//
//            if(j==results[i].size()-1)
//                cout<<"--------------------------------------------------"<<endl;
        }

    }

    cout<<"****************************************************************************"<<endl;

}


void IP_Packing::clearResult()
{
    fmpz_mat_clear(EncQ);

    for (ulong i=0; i<results.size(); i++)
    {
        for(ulong j=0; j<results[i].size(); j++)
        {
            fmpz_clear(results[i][j].ip);

            _fmpz_vec_clear(results[i][j].auxilary, 2*dim+1);
        }

        results[i].clear();
    }
    results.clear();
}

IP_Packing::~IP_Packing()
{
    fmpz_mat_clear(EncP);

    _fmpz_vec_clear(U, dim);

}

#endif // PACK_H_INCLUDED
