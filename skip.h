#ifndef SKIP_H_INCLUDED
#define SKIP_H_INCLUDED

#include <math.h>
#include <vector>
#include <queue>
#include <iomanip>

#include "ipe.h"


// IPTuple is used to store the top-k results

typedef struct
{
    fmpz_t  ip;         // inner product

    long    index;
    fmpz *  auxilary;  // the auxiliary vector

}IPTuple;



struct comparator
{
    bool operator()(IPTuple a, IPTuple b)
    {
        if(fmpz_cmp(a.ip, b.ip)>0)
            return true;
        else
            return false;
    }
};



class SKIP
{
private:

    vector<vector<IPTuple> > results;            // store the top-k answers of each query
    Monitor tt;                                 // use it to track the running time of each stage

    long u;
    long d;         // # of data vectors packed in one auxiliary vector
    long dim;       // dim of data and query vectors

    PublicParameter pp1;
    PublicParameter ppl;

    MasterSecretKey msk1;
    MasterSecretKey mskl;

    fmpz_mat_t EncP;           // encrypted auxiliary vectors
    fmpz_mat_t EncQ;           // encrypted queries

    fmpz_mat_t EncIndexNorm;   // encrypted index norm of each auxiliary vectors
    fmpz_mat_t EncL2Norm;      // encrypted query norm

    fmpz_t partialPNorm;       // the first value of all P (the same to all rows)
    fmpz*  partialQNorm;


    fmpz* U;                   // [1, u, u^2, ...,u^{d-1}]

    void pack(fmpz_mat_t bucket, fmpz* auxilary);

    void unpack(fmpz_t ip, long j, fmpz_t ip1);

    void unpack2(fmpz* auxilary, long j, fmpz* p);     // restore the j-th data vector in auxiliary of size len

public:

    void initialization(int bits, unsigned long maxIP, fmpz_mat_t P);

    void phase1(fmpz_mat_t Q, int queryNum=1);

    void phase2(int k);

    void phase3();

    void clearResult();

    void printTopK();                                   // print the results

    ~SKIP();
};


/***************************************************************************/

// auxilary = U * bucket

void SKIP::pack(fmpz_mat_t bucket, fmpz * auxilary)
{
    for(long i=0; i<bucket->c; i++)
    {
        fmpz_zero(auxilary+i);

        for(long j=0; j< bucket->r; j++)
            fmpz_addmul(auxilary+i, U+j, fmpz_mat_entry(bucket, j, i));
    }

}

void SKIP::unpack(fmpz_t ip, long j, fmpz_t ip1)
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

void SKIP::unpack2(fmpz* auxilary, long j, fmpz* p)
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

void SKIP::initialization(int bits, unsigned long maxIP, fmpz_mat_t P)
{

    /**********************************************************
     *                      generate key and parameters
     *********************************************************/

    initialize(&pp1, &msk1, 1);
    setup(&pp1, &msk1, bits);

    initialize(&ppl, &mskl, P->c);
    setup(&ppl, &mskl, bits);

    fmpz_init(partialPNorm);
    fmpz_set(partialPNorm, fmpz_mat_entry(P, 0, 0));

    u=maxIP+1;

    d=floor(fmpz_dlog(ppl.N)/log(u))-1;


    cout<<d<<endl;
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
    fmpz_mat_init(EncIndexNorm, n1, 3);

    fmpz_mat_t bucket;
    fmpz * auxilary;
    fmpz * index_norm;
    fmpz_t v;

    fmpz_mat_init(bucket, d, cols);
    auxilary   = _fmpz_vec_init(cols);
    index_norm = _fmpz_vec_init(1);
    fmpz_init(v);

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

        // set index norm of auxiliary

        fmpz_set_ui(index_norm, 0);

        for(long j=0; j<d; j++)
        {
           _fmpz_vec_dot (v, bucket->rows[j]+1, bucket->rows[j]+1, bucket->c-1);

           if(fmpz_cmp(v, index_norm)>0)
                fmpz_set(index_norm, v);
        }

        // encrypt [a], [||a||]

        encrypt(index_norm, EncIndexNorm->rows[i], &pp1, &msk1);

        encrypt(auxilary, EncP->rows[i], &ppl, &mskl);

    }

    tt.stop();

    // collect the the communication cost

    FILE *pfile;

    pfile=fopen("data/EncP.txt", "w+");

    fmpz_mat_fprint_pretty (pfile, EncP);

    fclose(pfile);


    pfile=fopen("data/EncIndexNorm.txt", "w+");

    fmpz_mat_fprint_pretty (pfile, EncIndexNorm);

    fclose(pfile);


    ofstream fout;
    fout.open("data/log.txt", ofstream::out | ofstream::app);
    fout<<endl;
    fout <<"SKIP: d="<<d <<" ";
    fout <<"init time="<<tt.getElapsedTime() << " ";
    fout << "outsourcing data size="<<evaluateFileSize("data/EncP.txt") + evaluateFileSize("data/EncIndexNorm.txt")<<endl;
    fout.close();


    remove("data/EncP.txt");
    remove("data/EncIndexNorm.txt");

    // free the memory

    fmpz_mat_clear(bucket);
    _fmpz_vec_clear(auxilary, dim);
    _fmpz_vec_clear(index_norm, 1);
    fmpz_clear(v);
    flint_randclear(state);

}


void SKIP::phase1(fmpz_mat_t Q, int queryNum)
{
    if(Q->r<queryNum)
    {
        perror("query number exceeds\n");
        exit(0);
    }

    tt.start();

    fmpz_mat_init(EncQ, queryNum, 2*Q->c+1);
    fmpz_mat_init(EncL2Norm, queryNum, 3);

    partialQNorm=_fmpz_vec_init(queryNum);

    fmpz * norm = _fmpz_vec_init(1);


    for(long i=0; i<queryNum; i++)
    {
        _fmpz_vec_dot (norm, Q->rows[i]+1, Q->rows[i]+1, Q->c-1);

        fmpz_set(partialQNorm+i, fmpz_mat_entry(Q, i, 0));

        keyGen(norm,  EncL2Norm->rows[i], &pp1, &msk1);

        keyGen(Q->rows[i], EncQ->rows[i], &ppl, &mskl);

    }

    _fmpz_vec_clear(norm, 1);

    tt.stop();

    ofstream fout;
    fout.open("data/log.txt", ofstream::out | ofstream::app);
    fout <<fixed<<setprecision(5)<<tt.getElapsedTime()/EncQ->r << " ";
    fout.close();
}


void SKIP::phase2(int k)
{
    tt.start();

    fmpz_t ip, ip1;

    fmpz_init(ip);

    fmpz_init(ip1);

    double pruning=0;

    for(long i=0; i<EncQ->r; i++)
    {
        priority_queue<IPTuple, vector<IPTuple>, comparator > topK;

        fmpz_t partialIP;

        fmpz_init(partialIP);

        fmpz_mul(partialIP, partialPNorm, partialQNorm+i);

        // push k empty tuples into the heap

        for(int j=0; j<k; j++)
        {
            IPTuple *t= new IPTuple;

            fmpz_init(t->ip);

            t->auxilary=_fmpz_vec_init(2*dim+1);

            fmpz_set_ui(t->ip, 0);

            topK.push(*t);
        }

        // sequential scan over IP Packing

        long j;


        for(j=0; j<EncP->r; j++)
        {
            // first try to prune with sequential scan

            decrypt(EncIndexNorm->rows[j], EncL2Norm->rows[i], ip, &pp1);        // compute (|a| * |q|)^2

//            fmpz_sub(ip1, topK.top().ip, partialIP);
//
//            fmpz_mul(ip1, ip1, ip1);

            fmpz_sqrt(ip, ip);

            fmpz_add(ip, ip, partialIP);

            fmpz_set(ip1, topK.top().ip);


            if(fmpz_cmp(ip, ip1) <=0)
                break;

            // if fails, regular scan

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

        pruning+=double(EncP->r-j)/EncP->r;

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

        fmpz_clear(partialIP);

    }

    fmpz_clear(ip);

    fmpz_clear(ip1);

    tt.stop();

    ofstream fout;
    fout.open("data/log.txt", ofstream::out | ofstream::app);
    fout << fixed<<setprecision(5)<<tt.getElapsedTime()/EncQ->r << " ";
    fout << fixed<<setprecision(3)<<pruning/EncQ->r << " ";
    fout.close();
}


void SKIP::phase3()
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

void SKIP::printTopK()
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

//            cout<<"->";
//
//            _fmpz_vec_print(results[i][j].auxilary+dim+1, dim);
//
//            cout<<endl;

//            if(j==results[i].size()-1)
//                   cout<<"--------------------------------------------------"<<endl;
        }

    }

     cout<<"****************************************************************************"<<endl;

}


void SKIP::clearResult()
{
    _fmpz_vec_clear(partialQNorm, EncQ->r);
    fmpz_mat_clear(EncQ);
    fmpz_mat_clear(EncL2Norm);

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


SKIP::~SKIP()
{
    clearup(&pp1, &msk1);
    clearup(&ppl, &mskl);

    fmpz_mat_clear(EncP);
    fmpz_mat_clear(EncIndexNorm);

    _fmpz_vec_clear(U, dim);
    fmpz_clear(partialPNorm);

}

#endif // SKIP_H_INCLUDED
