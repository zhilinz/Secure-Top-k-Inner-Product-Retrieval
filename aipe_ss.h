#ifndef AIPE_SS_H_INCLUDED
#define AIPE_SS_H_INCLUDED

#include <math.h>
#include <vector>
#include <queue>
#include <iomanip>

#include "skip.h"

// AIPE with sequential scan

class AIPE_SequentialScan
{
private:

    vector<vector<IPTuple> > results;            // store the top-k answers of each query
    Monitor tt;                                 // use it to track the running time of each stage

    long dim;                  // dim of data and query vectors

    PublicParameter pp1;
    PublicParameter ppl;

    MasterSecretKey msk1;
    MasterSecretKey mskl;

    fmpz_mat_t EncP;           // encrypted data vectors
    fmpz_mat_t EncQ;           // encrypted query vectors

    fmpz_mat_t EncIndexNorm;   // encrypted norm of each data vectors
    fmpz_mat_t EncL2Norm;      // encrypted query norm

    fmpz_t partialPNorm;       // the first value of all P (the same to all rows)
    fmpz*  partialQNorm;

public:

    void initialization(int bits, fmpz_mat_t P);

    void phase1(fmpz_mat_t Q, int queryNum=1);

    void phase2(int k);

    void phase3();

    void printTopK();                                   // print the results

    void clearResult();

    ~AIPE_SequentialScan();
};


/***************************************************************************/

void AIPE_SequentialScan::initialization(int bits, fmpz_mat_t P)
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

    dim=P->c;


    /**********************************************************
     *                          encryption
     *********************************************************/
    tt.start();

    long rows, cols;

    fmpz * norm = _fmpz_vec_init(1);

    rows=P->r;
    cols=P->c;

    fmpz_mat_init(EncP, rows, 2*cols+1);
    fmpz_mat_init(EncIndexNorm, rows, 3);

    // iterate on each record

    for (long i=0; i<rows; i++)
    {
        _fmpz_vec_dot (norm, P->rows[i]+1, P->rows[i]+1, P->c-1);

        // encrypt [p], [||p||]

        encrypt(norm, EncIndexNorm->rows[i], &pp1, &msk1);

        encrypt(P->rows[i], EncP->rows[i], &ppl, &mskl);

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
    fout <<"AIPE-SS: init time="<<tt.getElapsedTime() << " ";
    fout <<"outsourcing data size=" << evaluateFileSize("data/EncP.txt") + evaluateFileSize("data/EncIndexNorm.txt") <<endl;
    fout.close();

    remove("data/EncP.txt");
    remove("data/EncIndexNorm.txt");

    // free the memory

    _fmpz_vec_clear(norm, 1);

}


void AIPE_SequentialScan::phase1(fmpz_mat_t Q, int queryNum)
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


void AIPE_SequentialScan::phase2(int k)
{
    tt.start();

    fmpz_t ip, ip1, partialIP;

    fmpz_init(ip);
    fmpz_init(ip1);
    fmpz_init(partialIP);

    double pruning=0;

    for(long i=0; i<EncQ->r; i++)
    {
        fmpz_mul(partialIP, partialPNorm, partialQNorm+i);

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

        // sequential scan over IP Packing

        long j;

        for(j=0; j<EncP->r; j++)
        {
            // first try to prune with sequential scan

            decrypt(EncIndexNorm->rows[j], EncL2Norm->rows[i], ip, &pp1);        // compute (|p| * |q|)^2

            fmpz_sqrt(ip, ip);

            fmpz_add(ip, ip, partialIP);

            fmpz_set(ip1, topK.top().ip);

            if(fmpz_cmp(ip, ip1) <=0)
                break;


            // if fails, regular scan

            decrypt(EncP->rows[j], EncQ->rows[i], ip, &ppl);                     // compute q^T p

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

    }

    fmpz_clear(ip);

    fmpz_clear(ip1);

    fmpz_clear(partialIP);

    tt.stop();

    ofstream fout;
    fout.open("data/log.txt", ofstream::out | ofstream::app);
    fout << fixed<<setprecision(5)<<tt.getElapsedTime()/EncQ->r << " ";
    fout << fixed<<setprecision(3)<<pruning/EncQ->r << " ";
    fout.close();
}


void AIPE_SequentialScan::phase3()
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

void AIPE_SequentialScan::printTopK()
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
//                   cout<<"--------------------------------------------------"<<endl;
        }

    }

     cout<<"****************************************************************************"<<endl;

}


void AIPE_SequentialScan::clearResult()
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

AIPE_SequentialScan::~AIPE_SequentialScan()
{
    clearup(&pp1, &msk1);
    clearup(&ppl, &mskl);

    fmpz_mat_clear(EncP);
    fmpz_mat_clear(EncIndexNorm);
    fmpz_clear(partialPNorm);

}


#endif // AIPE_SS_H_INCLUDED
