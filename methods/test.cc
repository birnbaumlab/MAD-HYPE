
/*

MAD-HYPE Algorithm

This is the core scorer
It is written in C++ so that execution speed can be improved 20-fold

*/


// Library importation
#include <iterator>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <array>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <unordered_set>

using namespace std;

struct chainset_t {
    int alphas[2];
    int betas[2];
};

struct chain_t {
    bool is_alpha;
    int id;
};

///////////////////////
/// DATA PROCESSING ///
///      METHODS    ///
///////////////////////

// Index a 3D matrix using specified values
int Index(int a, int b, int c, int total){
    int ind = a*(total+1)*(total+1) + b*(total+1) + c;
    return ind;
}


// Print method for vectors
void Print(const vector<int>& v){
    for(unsigned i = 0; i< v.size(); ++i) {
        cout << v[i] << " ";
    }
    cout << endl;
}


void intersection(char *a, char *b, char *ab, int size) {
    for (int i = 0; i < size; i++)
        ab[i] = a[i]&b[i];
}
int sum(char *a, int size) {
    int r = 0;
    for (int i = 0; i < size; i++)
        r += a[i];
    return r;
}
int intersection_sum(char *a, char *b, int size)
{
    int r = 0; 
    for (int i = 0; i < size; i++)
    {
        //cout << a[i] << "  " << b[i] << "  " << r << endl;
        r += a[i]&b[i];
    }
    return r;
}

// Extract vector from string
vector<int> ParseCSString(string str)
{
    vector<int> vect;
    stringstream ss(str);
    int i;
    while (ss >> i)
    {
        vect.push_back(i);

        if (ss.peek() == ',')
            ss.ignore();
    }
    return vect;
}


// Finds the unique indices in data
vector<int> LoadUniques(string fname)
{
    // Initialize values 
    vector<int> uniques;
    
    // Process uniques file and pass to vector
    ifstream file(fname);
    string str;
    while (getline(file,str))
    {
        uniques.push_back(stoi(str));
    }

    // Return value
    return uniques;
}


// Returns counts of the number of wells per chain  
int* LoadChainCount(string fname, int unique_count)
{
    // Initialize values 
    vector<int> multilevel; 
    
    // Initialize data array
    int *data = new int[unique_count]; 
    
    // Process well file and pass to vector
    ifstream file(fname);
    int index = 0; string str;

    // Parse lines in file
    while (getline(file,str))
    {
        // Temporary storage of line as vector
        multilevel = ParseCSString(str);
        
        // Assign well #'s as 1's in matrix
        data[index] = (int)(multilevel.size());

        // Increase index
        index += 1;

    }

    return data;
}


// Well data load
char** LoadChainData(string fname,int unique_count,int w_tot) 
{
    // Initialize values 
    vector<int> multilevel; 
    
    // Initialize data array
    char** data = new char*[unique_count];
    for (int i = 0; i < unique_count; ++i)
    {
       data[i] = new char[w_tot];
       memset(data[i],0,sizeof(char)*w_tot);
    }
    
    // Process well file and pass to vector
    ifstream file(fname);
    int index = 0; string str;

    // Parse lines in file
    while (getline(file,str))
    {
        // Temporary storage of line as vector
        multilevel = ParseCSString(str);
        
        // Assign well #'s as 1's in matrix
        for (int j = 0; j < multilevel.size(); ++j)
        {
            data[index][multilevel[j]] = 1;
            //cout << data[index][multilevel[j]] << endl;
        }

        
        // Increase index
        index += 1;

    }


    // Return value
    return data;
}

void LoadInitChainsets(string &fname, vector<chainset_t> &init_chainsets) {
    ifstream file(fname);

    vector<int> splitline;
    string line;

    while (getline(file, line)) {
        // split string along commas
        splitline = ParseCSString(line);

        // assign vals to new chainset
        chainset_t cs;
        cs.alphas[0] = splitline[0];
        cs.alphas[1] = splitline[1];
        cs.betas[0] = splitline[2];
        cs.betas[1] = splitline[3];

        // append to list of chainsets
        init_chainsets.push_back(cs);
    }

    file.close();
}


/////////////////////
/// COMPUTATIONAL ///
///    METHODS    ///
/////////////////////

// Returns N choose K integer

double nCk(int n, int k)
{
    double res = 1;
    if ( k > n-k )
        k = n - k;
    for (int i = 0; i < k; ++i )
    {
        res *= (n-i);
        res/= (i+1);
    }
    return res;
}

double multinomial_prob4(int n1, int n2, int n3, int n4, double p1, double p2, double p3, double p4) {
    //double val = nCk(n1+n2, n2) * nCk(n1+n2+n3, n3) * nCk(n1+n2+n3+n4, n4) * pow(p1, n1) * pow(n2, p2) * pow(n3, p3) * pow(n4, p4);
    if ((n1>0 && p1==0) || (n2>0 && p2==0) || (n3>0 && p3==0) || (n4>0 && p4==0))
        return 0.0;
    double coeff, t1, t2, t3, t4;

    coeff = lgamma(n1+n2+n3+n4+1) - lgamma(n1+1) - lgamma(n2+1) - lgamma(n3+1) - lgamma(n4+1);
    t1 = p1!=0 ? n1*log(p1) : 0.0;
    t2 = p2!=0 ? n2*log(p2) : 0.0;
    t3 = p3!=0 ? n3*log(p3) : 0.0;
    t4 = p4!=0 ? n4*log(p4) : 0.0;
    return exp(coeff+t1+t2+t3+t4);
}
double multinomial_prob8(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8) {
    if ((n1>0 && p1==0) || (n2>0 && p2==0) || (n3>0 && p3==0) || (n4>0 && p4==0) || (n5>0 && p5==0) || (n6>0 && p6==0) || (n7>0 && p7==0) || (n8>0 && p8==0))
        return 0.0;
    double coeff, t1, t2, t3, t4, t5, t6, t7, t8;

    coeff = lgamma(n1+n2+n3+n4+1) - lgamma(n1+1) - lgamma(n2+1) - lgamma(n3+1) - lgamma(n4+1) - lgamma(n5+1) - lgamma(n6+1) - lgamma(n7+1) - lgamma(n8+1);
    t1 = p1!=0 ? n1*log(p1) : 0.0;
    t2 = p2!=0 ? n2*log(p2) : 0.0;
    t3 = p3!=0 ? n3*log(p3) : 0.0;
    t4 = p4!=0 ? n4*log(p4) : 0.0;
    t5 = p5!=0 ? n5*log(p5) : 0.0;
    t6 = p6!=0 ? n6*log(p6) : 0.0;
    t7 = p7!=0 ? n7*log(p7) : 0.0;
    t8 = p8!=0 ? n8*log(p8) : 0.0;
    return exp(coeff+t1+t2+t3+t4+t5+t6+t7+t8);
}


// Non-match MLE estimator for f_ab,f_a,f_b

void nonmatch_frequency(int w_ab,int w_a,int w_b,int w_tot,double& f_a,double& f_b)
{
    f_a = (double)(w_a + w_ab)/(double)(w_tot);
    f_b = (double)(w_b + w_ab)/(double)(w_tot);
}
void nonmatch_frequency_dual(int w_abc, int w_ab, int w_ac, int w_bc, int w_a, int w_b, int w_c, int w_tot, double &f_ab, double &f_ac, double &f_bc, double &f_a, double &f_b, double &f_c) {

    // This is complicated. Basically, we reformulate the well counts in terms of the following well counts:
    //   w_o = # wells with none of the chains
    //   w_nbc = # wells without b or c = # wells with only a or nothing
    //   w_nac = # wells without a or c
    //   w_nab = # wells without a or b
    //   w_na = # wells without a
    //   w_nb = # wells without b
    //   w_nc = # wells without c
    // The probabilities associated with each of these categories is easily expressible as a product of (1-f) values:
    //   p_o = (1-fa)*(1-fb)*(1-fc)*(1-fab)*(1-fac)*(1-fbc)
    //   p_nbc = (1-fb)*(1-fc)*(1-fab)*(1-fac)*(1-fbc)
    //   p_nac = (1-fa)*(1-fc)*(1-fab)*(1-fac)*(1-fbc)
    //   p_nab = (1-fa)*(1-fb)*(1-fab)*(1-fac)*(1-fbc)
    //   p_na = (1-fa)*(1-fab)*(1-fac)
    //   p_nb = (1-fb)*(1-fab)*(1-fbc)
    //   p_nc = (1-fc)*(1-fac)*(1-fbc)
    // If we take the logarithm of each equation, we can formulate these equations as a single matrix equation
    // in terms of the logarithms of the inverse cell frequencies as our unknowns. The coefficient matrix is:
    //  [1 1 1 1 1 1     [log(1-fa)        [log(p_o)
    //   0 1 1 1 1 1      log(1-fb)         log(p_nbc)
    //   1 0 1 1 1 1      log(1-fc)         log(p_nac)
    //   1 1 0 1 1 1  *   log(1-fab)   =    log(p_nab)
    //   1 0 0 1 1 0      log(1-fac)        log(p_na)
    //   0 1 0 1 0 1      log(1-fbc)]       log(p_nb)
    //   0 0 1 0 1 1]                       log(p_nc)]
    // Although this system of equations is overspecified, we can get the least squares estimate with the pseudoinverse
    // of the coefficient matrix, which yields
    //  [log(1-fa)             [ 5 -5  2  2 -2 -2 -2    [log(p_o)
    //   log(1-fb)               5  2 -5  2 -2 -2 -2     log(p_nbc)
    //   log(1-fc)         1     5  2  2 -5 -2 -2 -2     log(p_nac)
    //   log(1-fab)   =   --- * -3  3  3 -4  4  4 -3  *  log(p_nab)
    //   log(1-fac)        7    -3  3 -4  3  4 -3  4     log(p_na)
    //   log(1-fbc)]            -3 -4  3  3 -3  4  4]    log(p_nb)
    //                                                   log(p_nc)]
    // And then it is simple to determine the frequencies from these values.
    // Huge thanks to Sakul for a lot of this.

    int w_o = w_tot - w_a - w_b - w_c - w_ab - w_ac - w_bc - w_abc;

    double log_po, log_pnbc, log_pnac, log_pnab, log_pna, log_pnb, log_pnc;
    double log_na, log_nb, log_nc, log_nab, log_nac, log_nbc;
    
    log_po = log((double)w_o/w_tot);
    log_pnbc = log((w_o+w_a)/(double)w_tot);
    log_pnac = log((w_o+w_b)/(double)w_tot);
    log_pnab = log((w_o+w_c)/(double)w_tot);
    log_pna = log((w_o+w_b+w_c+w_bc)/(double)w_tot);
    log_pnb = log((w_o+w_a+w_c+w_ac)/(double)w_tot);
    log_pnc = log((w_o+w_a+w_b+w_ab)/(double)w_tot);

    log_na  = ( 5*log_po - 5*log_pnbc + 2*log_pnac + 2*log_pnab - 2*log_pna - 2*log_pnb - 2*log_pnc)/7;
    log_nb  = ( 5*log_po + 2*log_pnbc - 5*log_pnac + 2*log_pnab - 2*log_pna - 2*log_pnb - 2*log_pnc)/7;
    log_nc  = ( 5*log_po + 2*log_pnbc + 2*log_pnac - 5*log_pnab - 2*log_pna - 2*log_pnb - 2*log_pnc)/7;
    log_nab = (-3*log_po + 3*log_pnbc + 3*log_pnac - 4*log_pnab + 4*log_pna + 4*log_pnb - 3*log_pnc)/7;
    log_nac = (-3*log_po + 3*log_pnbc - 4*log_pnac + 3*log_pnab + 4*log_pna - 3*log_pnb + 4*log_pnc)/7;
    log_nbc = (-3*log_po - 4*log_pnbc + 3*log_pnac + 3*log_pnab - 3*log_pna + 4*log_pnb + 4*log_pnc)/7;

    f_a = max(0., 1 - exp(log_na));
    f_b = max(0., 1 - exp(log_nb));
    f_c = max(0., 1 - exp(log_nc));
    f_ab = max(0., 1 - exp(log_nab));
    f_ac = max(0., 1 - exp(log_nac));
    f_bc = max(0., 1 - exp(log_nbc));


//    if (w_a + w_o == 0)  f_a = 0;
//    else  f_a = (double)w_a / (w_a+w_o);
//    if (w_b + w_o == 0)  f_b = 0;
//    else  f_b = (double)w_b / (w_b+w_o);
//    if (w_c + w_o == 0)  f_c = 0;
//    else  f_c = (double)w_c / (w_c+w_o);
//
//    double temp1, temp2, temp3;
//    temp1 = (1 - (double)(w_a+w_ab+w_ac+w_abc)/w_tot)/(1 - f_a);
//    temp2 = (1 - (double)(w_b+w_ab+w_bc+w_abc)/w_tot)/(1 - f_b);
//    temp3 = (1 - (double)(w_c+w_ac+w_bc+w_abc)/w_tot)/(1 - f_c);
//
//    
//    f_ab = 1 - sqrt(temp1*temp2/temp3);
//    f_ac = 1 - sqrt(temp1*temp3/temp2);
//    f_bc = 1 - sqrt(temp2*temp3/temp1);
//    if (w_o + w_a + w_b + w_ab == 0)  f_ab = 0;
//    else  f_ab = max(0., ((double)w_ab / (w_o+w_a+w_b+w_ab) - f_a*f_b) / (1 - f_a*f_b));
//    if (w_o + w_a + w_c + w_ac == 0)  f_ac = 0;
//    else  f_ac = max(0., ((double)w_ac / (w_o+w_a+w_c+w_ac) - f_a*f_c) / (1 - f_a*f_c));
//    if (w_o + w_b + w_c + w_bc == 0)  f_bc = 0;
//    else  f_bc = max(0., ((double)w_bc / (w_o+w_b+w_c+w_bc) - f_b*f_c) / (1 - f_b*f_c));

    //if (w_ac==0 && w_bc==0 && w_b==0 && w_c==0)
//        cout << w_abc << " " << w_ab << " " << w_ac << " " << w_bc << " " << w_a << " " << w_b << " " << w_c << " " << f_a << " " << f_b << " " << f_c << " " << f_ab << " " << f_ac << " " << f_bc << endl;

}


// Match MLE estimator for f_ab,f_a,f_b

void match_frequency(int w_ab,int w_a,int w_b,int w_tot,double& f_ab,double& f_a,double& f_b)
{
    if (w_tot - w_ab - w_b == 0)  f_a = 0;
    else  f_a = (double)(w_a)/(double)(w_tot-w_ab-w_b);

    if (w_tot - w_ab - w_a == 0)  f_b = 0;
    else  f_b = (double)(w_b)/(double)(w_tot-w_ab-w_a);

    //f_ab = max(0.,1.-(1.-((double)(w_ab)/(double)(w_tot)))/(1.-f_a*f_b));
    f_ab = max(0.,(w_ab - f_a*f_b*w_tot) / ((1 - f_a*f_b)*w_tot));

    //cout << w_ab << " " << w_a << " " << w_b << " " << f_ab << " " << f_a << " " << f_b << endl;
}
void match_frequency_dual(int w_abc, int w_ab, int w_ac, int w_bc, int w_a, int w_b, int w_c, int w_tot, double &f_abc, double &f_ab, double &f_ac, double &f_bc, double &f_a, double &f_b, double &f_c) {
    int w_o = w_tot - w_a - w_b - w_c - w_ab - w_ac - w_bc - w_abc;

    if (w_a + w_o == 0)  f_a = 0;
    else  f_a = (double)w_a / (w_a+w_o);
    if (w_b + w_o == 0)  f_b = 0;
    else  f_b = (double)w_b / (w_b+w_o);
    if (w_c + w_o == 0)  f_c = 0;
    else  f_c = (double)w_c / (w_c+w_o);

    if (w_o + w_a + w_b + w_ab == 0)  f_ab = 0;
    else  f_ab = max(0., ((double)w_ab / (w_o+w_a+w_b+w_ab) - f_a*f_b) / (1 - f_a*f_b));
    //else  f_ab = ( ((double)w_ab / (w_o+w_a+w_b+w_ab) - f_a*f_b) / (1 - f_a*f_b));
    if (w_o + w_a + w_c + w_ac == 0)  f_ac = 0;
    else  f_ac = max(0., ((double)w_ac / (w_o+w_a+w_c+w_ac) - f_a*f_c) / (1 - f_a*f_c));
    //else  f_ac = ( ((double)w_ac / (w_o+w_a+w_c+w_ac) - f_a*f_c) / (1 - f_a*f_c));
    if (w_o + w_b + w_c + w_bc == 0)  f_bc = 0;
    else  f_bc = max(0., ((double)w_bc / (w_o+w_b+w_c+w_bc) - f_b*f_c) / (1 - f_b*f_c));
    //else  f_bc = ( ((double)w_bc / (w_o+w_b+w_c+w_bc) - f_b*f_c) / (1 - f_b*f_c));

    double f_false_abc = f_ab*(f_ac + (1-f_ac)*(f_bc + (1-f_bc)*f_c)) + (1-f_ab)*(f_ac*(f_bc + (1-f_bc)*f_b) + (1-f_ac)*(f_bc*f_a + (1-f_bc)*f_a*f_b*f_c));
    f_abc = max(0., ((double)w_abc / w_tot - f_false_abc) / (1 - f_false_abc));
    //f_abc = ( ((double)w_abc / w_tot - f_false_abc) / (1 - f_false_abc));
}

// Instantaneous probability for nonmatch instance

double nonmatch_instant_probability(int w_ab,int w_a,int w_b,int w_tot,double f_a,double f_b)
{
    double val =  multinomial_prob4(w_ab, w_a, w_b, w_tot-w_ab-w_a-w_b, f_a*f_b, f_a*(1-f_b), f_b*(1-f_a), (1-f_a)*(1-f_b));
    return val;
}


// Instantaneous probability for match instance

double match_instant_probability(int w_ab,int w_a,int w_b,int w_tot,double f_ab,double f_a,double f_b)
{
    double val =  multinomial_prob4(w_ab, w_a, w_b, w_tot-w_ab-w_a-w_b, f_a*f_b+f_ab-f_a*f_b*f_ab, f_a*(1-f_b)*(1-f_ab), f_b*(1-f_a)*(1-f_ab), (1-f_a)*(1-f_b)*(1-f_ab));
    return val;
}


// Non-match probability calculation

double nonmatch_probability(int w_ab,int w_a,int w_b,int w_tot)
{
    double f_a, f_b;
    nonmatch_frequency(w_ab,w_a,w_b,w_tot,f_a,f_b);
    double prob = nonmatch_instant_probability(w_ab,w_a,w_b,w_tot,f_a,f_b); 
    return prob;
}
double nonmatch_probability_dual(int w_abc, int w_ab, int w_ac, int w_bc, int w_a,int w_b, int w_c, int w_tot)
{
    double f_ab, f_ac, f_bc, f_a, f_b, f_c;
    nonmatch_frequency_dual(w_abc, w_ab, w_ac, w_bc, w_a, w_b, w_c, w_tot, f_ab, f_ac, f_bc, f_a, f_b, f_c);
    double n_ab, n_ac, n_bc, n_a, n_b, n_c;
    n_ab = 1-f_ab; n_ac = 1-f_ac; n_bc = 1-f_bc; n_a = 1-f_a; n_b = 1-f_b; n_c = 1-f_c;
    double p_abc, p_ab, p_ac, p_bc, p_a, p_b, p_c, p_o;
    p_abc = f_ab*(f_ac + n_ac*f_bc + n_ac*n_bc*f_c) + n_ab*(f_ac*(f_bc + n_bc*f_b) + n_ac*(f_bc*f_a + n_bc*f_a*f_b*f_c));
    p_ab = n_ac*n_bc*n_c*(f_ab + n_ab*f_a*f_b);
    p_ac = n_ab*n_bc*n_b*(f_ac + n_ac*f_a*f_c);
    p_bc = n_ab*n_ac*n_a*(f_bc + n_bc*f_b*f_c);
    p_a = n_ab*n_ac*n_bc*f_a*n_b*n_c;
    p_b = n_ab*n_ac*n_bc*n_a*f_b*n_c;
    p_c = n_ab*n_ac*n_bc*n_a*n_b*f_c;
    p_o = n_ab*n_ac*n_bc*n_a*n_b*n_c;

    //if (w_ac==0 && w_bc==0 && w_b==0 && w_c==0)
        //cout << w_abc << " " << w_ab << " " << w_ac << " " << w_bc << " " << w_a << " " << w_b << " " << w_c << " | " << f_a << " " << f_b << " " << f_c << " " << f_ab << " " << f_ac << " " << f_bc << " | " << p_abc << " " << p_ab << " " << p_ac << " " << p_bc << endl;
//        cout << w_abc << " " << w_ab << " " << w_ac << " " << w_bc << " " << w_a << " " << w_b << " " << w_c << " | " << p_abc*w_tot << " " << p_ab*w_tot  << " " << p_ac*w_tot  << " " << p_bc*w_tot << " " << p_a*w_tot << " " << p_b*w_tot << " " << p_c*w_tot  << " | " << f_a << " " << f_b << " " << f_c << " " << f_ab << " " << f_ac << " " << f_bc << endl;

    return multinomial_prob8(w_abc, w_ab, w_ac, w_bc, w_a, w_b, w_c, w_tot-w_abc-w_ab-w_bc-w_a-w_b-w_c, p_abc, p_ab, p_ac, p_bc, p_a, p_b, p_c, p_o);
}



// Match probability calculation

double match_probability(int w_ab,int w_a,int w_b,int w_tot)
{
    double f_a, f_b, f_ab;
    match_frequency(w_ab,w_a,w_b,w_tot,f_ab,f_a,f_b);
    if (f_ab==0)  return 0;
    double prob = match_instant_probability(w_ab,w_a,w_b,w_tot,f_ab,f_a,f_b); 
    return prob;   
}
double match_probability_dual(int w_abc,int w_ab, int w_ac, int w_bc, int w_a,int w_b, int w_c, int w_tot)
{ // TODO: check these probabilities
    double f_abc, f_ab, f_ac, f_bc, f_a, f_b, f_c;
    match_frequency_dual(w_abc, w_ab, w_ac, w_bc, w_a, w_b, w_c, w_tot, f_abc, f_ab, f_ac, f_bc, f_a, f_b, f_c);
    double n_abc, n_ab, n_ac, n_bc, n_a, n_b, n_c;
    n_abc = 1-f_abc; n_ab = 1-f_ab; n_ac = 1-f_ac; n_bc = 1-f_bc; n_a = 1-f_a; n_b = 1-f_b; n_c = 1-f_c;
    double p_abc, p_ab, p_ac, p_bc, p_a, p_b, p_c, p_o;
    p_abc = f_abc + n_abc*(f_ab*(f_ac + n_ac*f_bc + n_ac*n_bc*f_c) + n_ab*(f_ac*(f_bc + n_bc*f_b) + n_ac*(f_bc*f_a + n_bc*f_a*f_b*f_c)));
    p_ab = n_abc*n_ac*n_bc*n_c*(f_ab + f_a*f_b*n_ab);
    p_ac = n_abc*n_ab*n_bc*n_b*(f_ac + f_a*f_c*n_ac);
    p_bc = n_abc*n_ab*n_ac*n_a*(f_bc + f_b*f_c*n_bc);
    p_a = n_abc*n_ab*n_ac*n_bc*f_a*n_b*n_c;
    p_b = n_abc*n_ab*n_ac*n_bc*n_a*f_b*n_c;
    p_c = n_abc*n_ab*n_ac*n_bc*n_a*n_b*f_c;
    p_o = n_abc*n_ab*n_ac*n_bc*n_a*n_b*n_c;

    //cout << w_abc << " " << p_abc << " " << w_ab << " " << p_ab << " " << w_ac << " " << p_ac << " " << w_bc << " " << p_bc << " " << w_a << " " << p_a << " " << w_b << " " << p_b << " " << w_c << " " << p_c << endl;

    double prob = multinomial_prob8(w_abc, w_ab, w_ac, w_bc, w_a, w_b, w_c, w_tot-w_abc-w_ab-w_ac-w_bc-w_a-w_b-w_c, p_abc, p_ab, p_ac, p_bc, p_a, p_b, p_c, p_o);
    return prob;   
}






// Match score calculator

void match_score(int w_ab,int w_a,int w_b,int w_tot,double& score,double& freq)
{
    // If there are two or fewer matches, its unlikely to be real
    // JB - changed to 3 or fewer matches to be same as backup
    if ( w_ab <= 3 ){
        score = -numeric_limits<double>::infinity();
        freq = 0.f;
    }
    else {
        double f_a, f_b;
        double mp = match_probability(w_ab,w_a,w_b,w_tot);
        double nmp = nonmatch_probability(w_ab,w_a,w_b,w_tot);
        match_frequency(w_ab,w_a,w_b,w_tot,freq,f_a,f_b);
        score = log10(mp) - log10(nmp);
//        cout << w_ab << " " << w_a << " " << w_b << " " << mp << " " << nmp << " " << score << endl;
    }
}
void match_score_dual(int w_abc,int w_ab, int w_ac, int w_bc, int w_a,int w_b, int w_c, int w_tot,double& score,double& freq)
{
    // Calculates the match score for the chains a, b, and c
    // If there are three or fewer matches, its unlikely to be real
    if ( w_abc <= 3 ){
        score = 0.f;
        freq = 0.f;
    }
    else {
        double f_a, f_b, f_c;
        double f_ab, f_ac, f_bc;
        double mp = match_probability_dual(w_abc, w_ab, w_ac, w_bc,w_a,w_b, w_c,w_tot);
        double nmp = nonmatch_probability_dual(w_abc, w_ab, w_ac, w_bc,w_a,w_b, w_c,w_tot);
        match_frequency_dual(w_abc,w_ab,w_ac, w_bc,w_a,w_b,w_c,w_tot,freq,f_ab,f_ac,f_bc,f_a,f_b,f_c);
        score = log10(mp) - log10(nmp);
        //if (freq==0)  score = -numeric_limits<double>::infinity();
        //if (isinf(score))
        //  cout << w_abc << " " << w_ab << " " << w_ac << " " << w_bc << " " << w_a << " " << w_b << " " << w_c << " " << mp << " " << nmp << " " << score << endl;	
    }
}

// gross, there has got to be a better way :(
char* get_chain_data(chain_t &c, char **chain_data_a, char **chain_data_b) {
    if (c.is_alpha)
        return chain_data_a[c.id];
    else
        return chain_data_b[c.id];
}
int get_chain_count(chain_t &c, int *chain_count_a, int *chain_count_b) {
    if (c.is_alpha)
        return chain_count_a[c.id];
    else
        return chain_count_b[c.id];
}

int chainset_num_alphas(chainset_t &c) {
    return (c.alphas[0]!=-1) + (c.alphas[1]!=-1);
}
int chainset_num_betas(chainset_t &c) {
    return (c.betas[0]!=-1) + (c.betas[1]!=-1);
}
int chainset_size(chainset_t &c) {
    return (c.alphas[0]!=-1) + (c.alphas[1]!=-1) + (c.betas[0]!=-1) + (c.betas[1]!=-1);
}

chainset_t add_chain_to_chainset(chainset_t &cs, chain_t &c) {
    chainset_t new_cs = cs;
    if (c.is_alpha)
        new_cs.alphas[chainset_num_alphas(cs)] = c.id;
    else
        new_cs.betas[chainset_num_betas(cs)] = c.id;
    return new_cs;
}

string make_result_string(chainset_t cs, double freq, double score, vector<int> &uniques_a, vector<int> &uniques_b) {
    int a1, a2, b1, b2;
    a1 = cs.alphas[0]!=-1 ? uniques_a[cs.alphas[0]] : -1;
    a2 = cs.alphas[1]!=-1 ? uniques_a[cs.alphas[1]] : -1;
    b1 = cs.betas[0]!=-1 ? uniques_b[cs.betas[0]] : -1;
    b2 = cs.betas[1]!=-1 ? uniques_b[cs.betas[1]] : -1;
    stringstream ss;
    ss << score << "\t" << freq << "\t" << a1 << "\t" << a2 << "\t" << b1 << "\t" << b2;
    return ss.str(); 
}

// Attempt to add chains to an existing chainset
void try_chain_additions_nondual(chainset_t &chainset, vector<chain_t> &additions, float threshold, vector<int> &uniques_a, char **chain_data_a, int *chain_count_a, vector<int> &uniques_b, char **chain_data_b, int *chain_count_b, int w_tot, vector<string> &results) {
    // Calculate chainset well data
    char *chainset_welldata;
    int chainset_wellcount;
    if (chainset.alphas[0] != -1) {
        chainset_welldata = chain_data_a[chainset.alphas[0]];
        chainset_wellcount = chain_count_a[chainset.alphas[0]];
    } else {
        chainset_welldata = chain_data_b[chainset.betas[0]];
        chainset_wellcount = chain_count_b[chainset.betas[0]];
    }

    for (vector<chain_t>::iterator chain_it = additions.begin(); chain_it != additions.end(); chain_it++) {
        char *chain_welldata;
        int chain_wellcount;
        double score;
        double freq;

        // Retrieve data on chain being added
        chain_welldata = get_chain_data(*chain_it, chain_data_a, chain_data_b);
        chain_wellcount = get_chain_count(*chain_it, chain_count_a, chain_count_b);

        // Score probability of chainset+chain
        int w_ab = intersection_sum(chainset_welldata, chain_welldata, w_tot);
        int w_a = chainset_wellcount - w_ab;
        int w_b = chain_wellcount - w_ab;
        match_score(w_ab, w_a, w_b, w_tot, score, freq);
    
        // Store result if it's good enough
        //cout << "score: " << score << " | " << w_ab << " " << w_a << " " << w_b << endl;
        if (score > threshold) {
            chainset_t new_chainset = add_chain_to_chainset(chainset, *chain_it);
            results.push_back(make_result_string(new_chainset, freq, score, uniques_a, uniques_b)); 
        }
    }

}
void try_chain_additions_dual(chainset_t &chainset, vector<chain_t> &additions, float threshold, vector<int> &uniques_a, char **chain_data_a, int *chain_count_a, vector<int> &uniques_b, char **chain_data_b, int *chain_count_b, int w_tot, vector<string> &results) {
    // Calculate chainset well data
    char *chainset_a_welldata, *chainset_b_welldata, *chainset_ab_welldata;
    int chainset_a_wellcount, chainset_b_wellcount, chainset_ab_wellcount;

    //Assumes the chainset has exactly 1 alpha and 1 beta. We could generalize this without too much work.
    chainset_a_welldata = chain_data_a[chainset.alphas[0]];
    chainset_a_wellcount = chain_count_a[chainset.alphas[0]];
    chainset_b_welldata = chain_data_b[chainset.betas[0]];
    chainset_b_wellcount = chain_count_b[chainset.betas[0]];

    chainset_ab_welldata = new char[w_tot];
    intersection(chainset_a_welldata, chainset_b_welldata, chainset_ab_welldata, w_tot);
    chainset_ab_wellcount = sum(chainset_ab_welldata, w_tot);

    for (vector<chain_t>::iterator chain_it = additions.begin(); chain_it != additions.end(); chain_it++) {
        char *chain_welldata;
        int chain_wellcount;
        double score;
        double freq;

        // Retrieve data on chain being added
        chain_welldata = get_chain_data(*chain_it, chain_data_a, chain_data_b);
        chain_wellcount = get_chain_count(*chain_it, chain_count_a, chain_count_b);

        // Score probability of chainset+chain
        int w_abc, w_ab, w_ac, w_bc, w_a, w_b, w_c;
        w_abc = intersection_sum(chainset_ab_welldata, chain_welldata, w_tot);
        w_ab = chainset_ab_wellcount - w_abc;
        w_ac = intersection_sum(chainset_a_welldata, chain_welldata, w_tot) - w_abc;
        w_bc = intersection_sum(chainset_b_welldata, chain_welldata, w_tot) - w_abc;
        w_a = chainset_a_wellcount - w_ac - w_ab - w_abc;
        w_b = chainset_b_wellcount - w_bc - w_ab - w_abc;
        w_c = chain_wellcount - w_ac - w_bc - w_abc;
        match_score_dual(w_abc, w_ab, w_ac, w_bc, w_a, w_b, w_c, w_tot, score, freq);
    
        // Store result if it's good enough
        if (score > threshold) {
            chainset_t new_chainset = add_chain_to_chainset(chainset, *chain_it);
            results.push_back(make_result_string(new_chainset, freq, score, uniques_a, uniques_b)); 
        }
    }

    delete[] chainset_ab_welldata;

}
void try_chain_additions(chainset_t &chainset, vector<chain_t> &additions, float threshold, vector<int> &uniques_a, char **chain_data_a, int *chain_count_a, vector<int> &uniques_b, char **chain_data_b, int *chain_count_b, int w_tot, vector<string> &results) {
    switch (chainset_size(chainset)) {
        case 1:
            try_chain_additions_nondual(chainset, additions,threshold, uniques_a, chain_data_a, chain_count_a, uniques_b, chain_data_b, chain_count_b, w_tot, results);
            break;
        case 2:
            try_chain_additions_dual(chainset, additions,threshold, uniques_a, chain_data_a, chain_count_a, uniques_b, chain_data_b, chain_count_b, w_tot, results);
            break;
        case 3:
        break;
    }
}


int main(int argc, char *argv[])

{
    //cout << "starting... \n";
    //cout << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << endl;
    int w_tot = stoi(argv[1]);
    float threshold = stof(argv[2]);
    bool try_alphas = stoi(argv[3]);
    bool try_betas = stoi(argv[4]);
    int index = stoi(argv[5]);
    
    // Input parameters
    //cout << "Starting Process-" << index << " with parameters:\n";
    //cout << "  W_tot: " << w_tot << endl;
    //cout << "  Threshold: " << threshold << endl;
    //cout << "  try_alphas: " << try_alphas << endl;
    //cout << "  try_betas:  " << try_betas << endl;
    
    //const int w_tot = 96;
    //const float threshold = 4.0;
    
    // Initialize the storage variables
/*
    double* scores = new double[(w_tot+1)*(w_tot+1)*(w_tot+1)](); 
    double* freqs = new double[(w_tot+1)*(w_tot+1)*(w_tot+1)](); 
    
    // Fill storage variables
    for( int w_ab = 0; w_ab < w_tot+1; w_ab += 1 ){
        for( int w_a = 0; w_a < w_tot+1; w_a += 1 ){
            for( int w_b = 0; w_b < w_tot+1; w_b += 1 ){
                if ( w_ab+w_a+w_b <= w_tot ) {
                    double score;
                    double freq;  // Assign empty score/freq variables
                    match_score(w_ab,w_a,w_b,w_tot,score,freq);
                    scores[Index(w_ab,w_a,w_b,w_tot)] = score;
                    freqs[Index(w_ab,w_a,w_b,w_tot)] = freq;
                }
            }
        }
        cout << "Progress " << w_ab << "/" << w_tot << "    \r";
        cout.flush();
    }
    cout << endl;
 */   
    
    // Declare data variables
    string fname_init_chainsets = "./solver/initial_" + to_string(index) + ".txt";
    string fname_data_a = "./solver/chain_data_a.txt";
    string fname_data_b = "./solver/chain_data_b.txt";
    string fname_uniques_a = "./solver/uniques_a.txt";
    string fname_uniques_b = "./solver/uniques_b.txt";

    // Load unique chains for a/b
    //cout << "Loading unique chains..." << endl;
    vector<int> uniques_a = LoadUniques(fname_uniques_a);
    vector<int> uniques_b = LoadUniques(fname_uniques_b);
    //cout << uniques_a.size() << "/" << uniques_b.size() << " unique a/b chains loaded!" << endl;
    
    // Load well data for a/b 
    //cout << "Loading well data..." << endl;
    char** chain_data_a = LoadChainData(fname_data_a,uniques_a.size(),w_tot);
    //cout << "Finished loading chain data A!" << endl;
    char** chain_data_b = LoadChainData(fname_data_b,uniques_b.size(),w_tot);
    //cout << "Finished loading chain data B!" << endl;

    // Load well data for a/b 
    //cout << "Loading well data..." << endl;
    int *chain_count_a = LoadChainCount(fname_data_a,uniques_a.size());
    int *chain_count_b = LoadChainCount(fname_data_b,uniques_b.size());
    //cout << "Finished loading chain counts!" << endl;

    //cout << "Loading starting chainsets from " << fname_init_chainsets << ": ";
    vector<chainset_t> init_chainsets; // stores indices into uniques_a/uniques_b (or -1 if no chain)
    LoadInitChainsets(fname_init_chainsets, init_chainsets);
    //cout << init_chainsets.size() << " initial chain sets loaded!";

    vector<string> results;

    // Create vector of chain additions to be tested
    vector<chain_t> additions;
    if (try_alphas) {
        for (int i = 0; i < uniques_a.size(); i++) {
            additions.push_back((chain_t){.is_alpha=true, .id=i});
            //cout << i << " " << additions.back().is_alpha << " " << additions.back().id << endl;
        }
    }
    if (try_betas) {
        for (int i = 0; i < uniques_b.size(); i++) {
            additions.push_back((chain_t){.is_alpha=false, .id=i});
        }
    }
    //cout << "trying " << additions.size() << " additions";
    //cout << additions.back().is_alpha << " " << additions.back().id << endl;
    
    for (int i = 0; i < init_chainsets.size(); i++) {
        try_chain_additions(init_chainsets[i], additions, threshold, uniques_a, chain_data_a, chain_count_a, uniques_b, chain_data_b, chain_count_b, w_tot, results);
        //cout << "Finished " << i+1 << "/" << init_chainsets.size() << "      \r";
        //cout.flush();
    }
    //cout << endl;
    // Iterate through pairs 
    /*for (int i = 0; i < uniques_a.size(); i ++)
    {
        for (int j = 0; j < uniques_b.size(); j ++)
        {
            // Scores and freqs
            int w_ab = intersection_sum(chain_data_a[i],chain_data_b[j],w_tot);
            int w_a = chain_count_a[i] - w_ab; 
            int w_b = chain_count_b[j] - w_ab; 
            score = scores[Index(w_ab,w_a,w_b,w_tot)];
            freq = freqs[Index(w_ab,w_a,w_b,w_tot)];

            // Save results that have a score above threshold
            if (score > threshold) 
            {
                stringstream ss;
                ss << score << "\t" << freq << "\t" << uniques_a[i] << "\t" << uniques_b[j];
                results.push_back(ss.str()); 
            }
        } 
        cout << "Finished " << i+1 << "/" << uniques_a.size() << "      \r";
        //cout.flush();
    }
    cout << endl;
    */

    // Output results to txt file
    ofstream output_file("./solver/results_" + to_string(index) + ".txt");
    ostream_iterator<string> output_iterator(output_file, "\n");
    copy(results.begin(), results.end(), output_iterator); 

};




