#include <cstdio>
#include <random>
#include <cmath>
#include <ctime>

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// the code uses variable length array for convenience (and performance)

using namespace std;

#define L  32
#define NR 4

int s0[L][L];
int s1[L/2][L/2];
int s2[L/2/2][L/2/2];
int s3[L/2/2/2][L/2/2/2];

void *ptrs[] = {s0, s1, s2, s3};

int (*s)[L] = (int (*)[L])ptrs[0];

#define K 0.44068679350977147

#define TO_STRING(STR) #STR
#define EXPAND_TO_STRING(STR) TO_STRING(STR)
const char filename[] = "plots/ops_l" EXPAND_TO_STRING(L) "_maj.dat";
const char filename_callen[] = "plots/callen_l" EXPAND_TO_STRING(L) "_maj.dat";


vector<vector<vector<vector<int>>>> pts;
vector<vector<vector<vector<int>>>> extend;

#define N (L*L)


// global variables for wolff algorithm
int cluster; // cluster size
double padd = 1.0 - exp(-2.0*K);
int istack[N];
int jstack[N];


// random generator
random_device rd;
//mt19937_64 gen(rd());
mt19937_64 gen(6432524);
uniform_real_distribution<> dis(0.0, 1.0);


void initialize()
{
    // initialize lattice of spins
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            //if (dis(gen) < 0.5) s[i][j] = 1;
            //else s[i][j] = -1;
            s[i][j] = 1;
        }
    }

    // parse operator files to vector pts
    ifstream infile("ops.txt");
    string line;
    while (getline(infile, line)) {
        istringstream iss(line);

        vector<vector<vector<int>>> inds;

        char c;
        iss >> c;

        while (true) {
            iss >> c;

            vector<vector<int>> ind;
            int i;
            while ((iss >> i)) {
                ind.push_back({i / 3, i % 3});
                iss >> c;
                if (c != ',') break;
            }

            inds.push_back(ind);

            iss >> c;
            if (c != ',') break;
        }

        pts.push_back(inds);
    }

    // extend for duplicate used in the coupling analysis
    for (size_t c = 0; c < pts.size(); ++c) {
        vector<vector<vector<int>>> inds;
        for (auto & vs : pts[c]) {
            for (auto & v : vs) {
                // double looping - prepare original and duplicate
                vector<vector<int>> ind;
                for (auto & vv : vs) {
                    ind.push_back({vv[0]-v[0], vv[1]-v[1]});
                }
                inds.push_back(ind);
            }
        }
        extend.push_back(inds);
    }
}

void step()
{
    int i, j; // current site popped from stack
    int in, jn; // neighbor sites

    int sp; // stack size
    int oldspin, newspin; // old and new spin alignment

    cluster = 1; // reset 
    
    // choose the seed spin for the cluster
    i = L * dis(gen);
    j = L * dis(gen);

    sp = 1;
    istack[0] = i;
    jstack[0] = j;

    oldspin =  s[i][j];
    newspin = -oldspin;
    s[i][j] = newspin; // flip the seed site

    // check the neighbors
    while (sp) {
        sp--;
        i = istack[sp];
        j = jstack[sp];

        if ((in = i + 1) >= L) in -= L;
        if (s[in][j] == oldspin)
            if (dis(gen) < padd) {
                istack[sp] = in;
                jstack[sp] = j;
                sp++;
                s[in][j] = newspin;
                cluster++;
            }

        if ((in = i - 1) <  0) in += L;
        if (s[in][j] == oldspin)
            if (dis(gen) < padd) {
                istack[sp] = in;
                jstack[sp] = j;
                sp++;
                s[in][j] = newspin;
                cluster++;
            }

        if ((jn = j + 1) >= L) jn -= L;
        if (s[i][jn] == oldspin)
            if (dis(gen) < padd) {
                istack[sp] = i;
                jstack[sp] = jn;
                sp++;
                s[i][jn] = newspin;
                cluster++;
            }

        if ((jn = j - 1) <  0) jn += L;
        if (s[i][jn] == oldspin)
            if (dis(gen) < padd) {
                istack[sp] = i;
                jstack[sp] = jn;
                sp++;
                s[i][jn] = newspin;
                cluster++;
            }
    }
}

void renormalize_majority()
{
    int l = L;

    for (int r = 0; r < NR-1; ++r) {

        int (*s)[l] = (int (*)[l])ptrs[r];

        l /= 2;
        int (*sb)[l] = (int (*)[l])ptrs[r+1];

        for (int i = 0; i < l; ++i)
        for (int j = 0; j < l; ++j) {
            // each block
            int sum = 0;
            for (int ib = 2*i; ib < 2*i+2; ++ib)    
                for (int jb = 2*j; jb < 2*j+2; ++jb)    
                    sum += s[ib][jb];

            if (sum > 0) sb[i][j] = 1;
            else if (sum < 0) sb[i][j] = -1;
            else {
                if (dis(gen) < 0.5) sb[i][j] = 1;
                else sb[i][j] = -1;
            }
        }
    }
}


void operate(FILE *pf)
{
    int l = L;

    for (int r = 0; r < NR; ++r) {

        int (*s)[l] = (int (*)[l])ptrs[r];

        int ops[pts.size()]; // variable length array convenient
        for (size_t c = 0; c < pts.size(); ++c)
            ops[c] = 0;

        for (int i = 0; i < l; ++i)
        for (int j = 0; j < l; ++j) {
            for (size_t c = 0; c < pts.size(); ++c) {
                for (auto & vs : pts[c]) {
                    int prod = 1;
                    for (auto & v : vs) {
                        int in, jn;
                        if ((in = i + v[0]) >= l) in -= l;
                        if ((jn = j + v[1]) >= l) jn -= l;
                        prod *= s[in][jn];
                    }
                    ops[c] += prod;
                }
            }
        }

        fwrite(ops, sizeof(int), pts.size(), pf);

        // reduce block size
        l /= 2;
    }
}


void callen(FILE *pf)
{
    int l = L;

    for (int r = 0; r < NR; ++r) {

        int (*s)[l] = (int (*)[l])ptrs[r];

        int hats[extend.size()];

        for (size_t c = 0; c < extend.size(); ++c)
            hats[c] = 0;

        for (int i = 0; i < l; ++i)
        for (int j = 0; j < l; ++j) {

            // compute hat ops
            for (size_t c = 0; c < extend.size(); ++c) {

                hats[c] = 0;

                for (auto & vs : extend[c]) {

                    int prod = 1;

                    for (auto & v : vs) {
                        int in = i + v[0]; 
                        int jn = j + v[1];

                        if (in >= l) in -= l;
                        if (jn >= l) jn -= l;
                        if (in < 0) in += l;
                        if (jn < 0) jn += l;

                        prod *= s[in][jn];
                    }

                    hats[c] += prod / s[i][j];
                }
            }

            fwrite(hats, sizeof(int), extend.size(), pf);

        }


        // reduce block size
        l /= 2;
    }
}

int main()
{
    initialize();

#define EQUIL 100000
#define EVERY 5
#define NCONF 100000

    struct timespec start, finish;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);


    FILE *pf = fopen(filename, "wb");
    FILE *pfc = fopen(filename_callen, "wb");

    for (int k = 0; k <= EQUIL+EVERY*NCONF; ++k) {

        step();

        if (k > EQUIL && k % EVERY == 0) {
            renormalize_majority();
            operate(pf);
            callen(pfc);

            // flush progress
            if ((((k-EQUIL)/EVERY) % (NCONF/100)) == 0) {
                clock_gettime(CLOCK_MONOTONIC, &finish);
                elapsed = (finish.tv_sec - start.tv_sec);
                elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
                fprintf(stderr, "%4d%%   %fs", (k-EQUIL)/EVERY/(NCONF/100), elapsed);
                //fprintf(stderr, "\r");
            }
        }
    }

    fclose(pf);
    fclose(pfc);
}


