#include <iostream>
#include <fstream>
#include <sys/time.h>


#include "skip.h"
#include "pack.h"
#include "aipe.h"
#include "aipe_ss.h"

using namespace std;


// this function is used to generate the initial paillier encryption key,
// which determine the value of N used in generating the invertible matrix in Z_N

int main(int argc, char* argv[])
{
    srand(time(NULL));

    int bits;
    vector<int> k_list;
    long maxip;
    string dataset;
    string alg;

    if(argc==6)
    {
        alg=argv[1];


        stringstream ss(argv[2]);
        int k;

        while(ss>>k)
        {
            k_list.push_back(k);

            if(ss.peek()==',')
                ss.ignore();
        }

        dataset=argv[3];
        maxip=atoi(argv[4]);
        bits=atoi(argv[5]);
    }
    else
    {
        alg="SKIP";
        k_list.push_back(1);
        dataset="./data/MovieLens/";
        maxip=100000;
        bits=1024;
    }


    fmpz_mat_t P;
    fmpz_mat_t Q;

    fmpz_mat_readFile(P, dataset+"IntP.txt");
    fmpz_mat_readFile(Q, dataset+"IntQ.txt");


    if(alg=="SKIP")
    {
        SKIP skip;

        for (unsigned i=0; i<k_list.size();i++)
        {
            cout<<"alg="<<alg<<", k="<<k_list[i]<<", dataset=\""<<dataset<<"\", maxip="<<maxip<<", bits="<<bits<<endl;

            if(i==0)
                skip.initialization(bits, maxip, P);

            skip.phase1(Q, 10);

            skip.phase2(k_list[i]);

            skip.phase3();

            skip.printTopK();

            skip.clearResult();
        }

    }
    else if(alg=="IP-PACKING")
    {
        IP_Packing ippack;

        for (unsigned i=0; i<k_list.size();i++)
        {
            cout<<"alg="<<alg<<", k="<<k_list[i]<<", dataset=\""<<dataset<<"\", maxip="<<maxip<<", bits="<<bits<<endl;

            if(i==0)
                ippack.initialization(bits, maxip, P);

            ippack.phase1(Q, 1);

            ippack.phase2(k_list[i]);

            ippack.phase3();

            ippack.printTopK();

            ippack.clearResult();
        }

    }
    else if(alg=="AIPE")
    {
        AIPE aipe;

        for (unsigned i=0; i<k_list.size();i++)
        {
            cout<<"alg="<<alg<<", k="<<k_list[i]<<", dataset=\""<<dataset<<"\", maxip="<<maxip<<", bits="<<bits<<endl;

            if(i==0)
                  aipe.initialization(bits, P);

            aipe.phase1(Q, 1);

            aipe.phase2(k_list[i]);

            aipe.phase3();

            aipe.printTopK();

            aipe.clearResult();
        }

    }
    else if(alg=="AIPE-SS")
    {
        AIPE_SequentialScan aipe_ss;

        for (unsigned i=0; i<k_list.size();i++)
        {
            cout<<"alg="<<alg<<", k="<<k_list[i]<<", dataset=\""<<dataset<<"\", maxip="<<maxip<<", bits="<<bits<<endl;

            if(i==0)
                  aipe_ss.initialization(bits, P);

            aipe_ss.phase1(Q, 10);

            aipe_ss.phase2(k_list[i]);

            aipe_ss.phase3();

            aipe_ss.printTopK();

            aipe_ss.clearResult();
        }
    }

}
