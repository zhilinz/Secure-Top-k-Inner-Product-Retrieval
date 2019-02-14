#ifndef AIPE_H_INCLUDED
#define AIPE_H_INCLUDED

#include <math.h>
#include <vector>
#include <queue>
#include <iomanip>

#include "skip.h"


// the pure AIPE without ip packing and sequential scan

class AIPE
{
private:

    vector<vector<IPTuple> > results;            // store the top-k answers of each query
    Monitor tt;                                 // use it to track the running time of each stage

    long dim;                    // dim of data and query vectors

    PublicParameter ppl;
    MasterSecretKey mskl;

    fmpz_mat_t EncP;           // encrypted data vectors
    fmpz_mat_t EncQ;           // encrypted query vectors

public:

    void initialization(int bits, fmpz_mat_t P);

    void phase1(fmpz_mat_t Q, int queryNum=1);

    void phase2(int k);

    void phase3();

    void printTopK();                                   // print the results

    void clearResult();

    ~AIPE();
};



/***************************************************************************/

void AIPE::initialization(int bits, fmpz_mat_t P)
{

    /**********************************************************
     *                      generate key and parameters
     *********************************************************/

    initialize(&ppl, &mskl, P->c);
    setup(&ppl, &mskl, bits);


    dim=P->c;

    /**********************************************************
     *                          encryption
     *********************************************************/
    tt.start();

    long rows, cols;

    rows=P->r;
    cols=P->c;

    fmpz_mat_init(EncP, rows, 2*cols+1);

    // iterate on each record

    for (long i=0; i<rows; i++)
    {
        encrypt(P->rows[i], EncP->rows[i], &ppl, &mskl);
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
    fout <<"AIPE: init time="<<tt.getElapsedTime() << " ";
    fout <<"outsourcing data size=" << evaluateFileSize("data/EncP.txt") << endl;
    fout.close();

    remove("data/EncP.txt");
}


void AIPE::phase1(fmpz_mat_t Q, int queryNum)
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


void AIPE::phase2(int k)
{
    tt.start();

    fmpz_t ip;

    fmpz_init(ip);

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

        // linear scan on every encrypted data

        for(long j=0; j<EncP->r; j++)
        {

            decrypt(EncP->rows[j], EncQ->rows[i], ip, &ppl);              // compute q^T p

            IPTuple result;

            if(fmpz_cmp(ip, topK.top().ip) >0)
            {
                result=topK.top();

                topK.pop();

                fmpz_set(result.ip, ip);                                  // set ip

                _fmpz_vec_set(result.auxilary, EncP->rows[j], 2*dim+1);

                topK.push(result);
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

    tt.stop();

    ofstream fout;
    fout.open("data/log.txt", ofstream::out | ofstream::app);
    fout << fixed<<setprecision(5)<<tt.getElapsedTime()/EncQ->r << " ";
    fout << fixed<<setprecision(3)<<0.0<< " ";
    fout.close();
}


void AIPE::phase3()
{
    tt.start();

    for (ulong i=0; i<results.size(); i++)
    {
        // deal with the answer set of i-th query

        for(ulong j=0; j<results[i].size(); j++)
        {
            decrypt2(results[i][j].auxilary, &ppl, &mskl);    // the plain auxilary vector is stored in the last dim items
        }

    }

    tt.stop();

    ofstream fout;
    fout.open("data/log.txt", ofstream::out | ofstream::app);
    fout << fixed <<setprecision(5)<< tt.getElapsedTime()/EncQ->r <<endl;
    fout.close();

}

void AIPE::printTopK()
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
//            cout<<" -> ";
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


void AIPE::clearResult()
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

AIPE::~AIPE()
{
    clearup(&ppl, &mskl);

    fmpz_mat_clear(EncP);
}

#endif // AIPE_H_INCLUDED
