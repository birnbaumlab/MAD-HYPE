
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
#include <algorithm>
#include <vector>
#include <array>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <unordered_set>

using namespace std;

///////////////////////
/// DATA PROCESSING ///
///      METHODS    ///
///////////////////////

// Print method for vectors
void Print(const vector<int>& v){
    for(unsigned i = 0; i< v.size(); ++i) {
        cout << v[i] << " ";
    }
    cout << endl;
}

int intersection(int *a, int *b, int size)
{
    int r = 0; 
    for (int i = 0; i < size; i++)
    {
        //cout << a[i] << "  " << b[i] << "  " << r << endl;
        r += a[i]*b[i];
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
vector<int> LoadChainCount(string fname, int unique_count)
{
    // Initialize values 
    vector<int> multilevel; 
    
    // Initialize data array
    vector<int> data(unique_count); 
    
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
int** LoadChainData(string fname,int unique_count,int w_tot) 
{
    // Initialize values 
    vector<int> multilevel; 
    
    // Initialize data array
    int** data = new int*[unique_count];
    for (int i = 0; i < unique_count; ++i)
    {
       data[i] = new int[w_tot];
       memset(data[i],0,sizeof(int)*w_tot);
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


// Non-match MLE estimator for f_ab,f_a,f_b

float nonmatch_frequency(int w_ab,int w_a,int w_b,int w_tot,float& f_a,float& f_b,float& f_ab)
{
    f_ab = (float)(w_ab)/(float)(w_tot);
    f_a = (float)(w_a + w_ab)/(float)(w_tot);
    f_b = (float)(w_b + w_ab)/(float)(w_tot);
}


// Match MLE estimator for f_ab,f_a,f_b

float match_frequency(int w_ab,int w_a,int w_b,int w_tot,float& f_a,float& f_b,float& f_ab)
{
    if ( w_tot-w_ab == 0 ){ 
        f_a = 0.f; 
        f_b = 0.f; 
    }
    else {
        f_a = (float)(w_a)/(float)(w_tot-w_ab);
        f_b = (float)(w_b)/(float)(w_tot-w_ab);
    }
    f_ab = max(0.f,1.f-(1.f-((float)(w_ab)/(float)(w_tot)))/(1.f-f_a*f_b));
}

// Instantaneous probability for nonmatch instance

float nonmatch_instant_probability(int w_ab,int w_a,int w_b,int w_tot,float f_ab,float f_a,float f_b)
{
    // We need to operate in logs unfortunately, since floating decimals seem to lose information
    float a = nCk(w_tot,w_a+w_ab)*pow(f_a,w_a+w_ab)*pow((1.-f_a),w_tot-(w_a+w_ab)); 
    float b = nCk(w_a+w_ab,w_ab)*pow(f_b,w_ab)*pow((1.-f_b),w_a); 
    float c = nCk(w_tot-(w_a+w_ab),w_b)*pow(f_b,w_b)*pow((1.-f_b),w_tot-(w_a+w_b+w_ab)); 

    //DELETE 
    if (a*b*c > 1)
    {
        cout << w_ab << "  "<< w_a << "  "<< w_b << "  " << w_tot << "  " << f_b << endl;
        cout << a << "  " << b << "  " << c << "  letters" << endl;
        cout << endl;
    }


    return a*b*c;
}


// Instantaneous probability for match instance

float match_instant_probability(int w_ab,int w_a,int w_b,int w_tot,float f_ab,float f_a,float f_b)
{
    float prob = 0.f;float int_prob;
    for( int w = 0.f; w < w_ab+1; w += 1 ){
        int_prob = nonmatch_instant_probability(w_ab-w,w_a,w_b,w_tot-w,f_ab,f_a,f_b);        
        prob += (nCk(w_tot,w)*pow(f_ab,w)*pow((1-f_ab),(w_tot-w)))*int_prob;
    }
    return prob;
}


// Non-match probability calculation

float nonmatch_probability(int w_ab,int w_a,int w_b,int w_tot)
{
    float f_a; float f_b; float f_ab;
    nonmatch_frequency(w_ab,w_a,w_b,w_tot,f_a,f_b,f_ab);
    float prob = nonmatch_instant_probability(w_ab,w_a,w_b,w_tot,f_ab,f_a,f_b); 
    return prob;
}


// Match probability calculation

float match_probability(int w_ab,int w_a,int w_b,int w_tot)
{
    float f_a; float f_b; float f_ab;
    match_frequency(w_ab,w_a,w_b,w_tot,f_a,f_b,f_ab);
    float prob = match_instant_probability(w_ab,w_a,w_b,w_tot,f_ab,f_a,f_b); 
    return prob;   
}





// Match score calculator

void match_score(int w_ab,int w_a,int w_b,int w_tot,float& score,float& freq)
{
    // If there are two or fewer matches, its unlikely to be real
    if ( w_ab <= 2 ){
        score = 0.f;
        freq = 0.f;
    }
    else {
        float f_a; float f_b;
        float mp = match_probability(w_ab,w_a,w_b,w_tot);
        float nmp = nonmatch_probability(w_ab,w_a,w_b,w_tot);
        match_frequency(w_ab,w_a,w_b,w_tot,f_a,f_b,freq);
        score = log10(mp) - log10(nmp);
    }
}

int main()

{
    
    const int w_tot = 96;
    const float threshold = 4.0;
    
    // Initialize the storage variables
    float scores[w_tot+1][w_tot+1][w_tot+1] = {}; 
    float freqs[w_tot+1][w_tot+1][w_tot+1] = {}; 

    // Fill storage variables
    for( int w_ab = 0; w_ab < w_tot+1; w_ab += 1 ){
        for( int w_a = 0; w_a < w_tot+1; w_a += 1 ){
            for( int w_b = 0; w_b < w_tot+1; w_b += 1 ){
                if ( w_ab+w_a+w_b <= w_tot ) {
                    float score, freq;  // Assign empty score/freq variables
                    match_score(w_ab,w_a,w_b,w_tot,score,freq);
                    scores[w_ab][w_a][w_b] = score;
                    freqs[w_ab][w_a][w_b] = freq;
                }
            }
        }
        cout << "Progress " << w_ab << "/" << w_tot << "    \r";
        cout.flush();
    }
    cout << endl;
    
    // Declare data variables
    string fname_data_a = "chain_data_a.txt";
    string fname_data_b = "chain_data_b.txt";
    string fname_uniques_a = "uniques_a.txt";
    string fname_uniques_b = "uniques_b.txt";

    // Load unique chains for a/b
    cout << "Loading unique chains..." << endl;
    vector<int> uniques_a = LoadUniques(fname_uniques_a);
    vector<int> uniques_b = LoadUniques(fname_uniques_b);
    cout << uniques_a.size() << "/" << uniques_b.size() << " unique a/b chains loaded!" << endl;
    
    // Load well data for a/b 
    cout << "Loading well data..." << endl;
    int** chain_data_a = LoadChainData(fname_data_a,uniques_a.size(),w_tot);
    cout << "Finished loading chain data A!" << endl;
    int** chain_data_b = LoadChainData(fname_data_b,uniques_b.size(),w_tot);
    cout << "Finished loading chain data B!" << endl;

    // Load well data for a/b 
    cout << "Loading well data..." << endl;
    vector<int> chain_count_a = LoadChainCount(fname_data_a,uniques_a.size());
    vector<int> chain_count_b = LoadChainCount(fname_data_b,uniques_b.size());
    cout << "Finished loading chain counts!" << endl;

    vector<string> results;

    // Iterate through pairs 
    for (int i = 0; i < uniques_a.size(); i ++)
    {
        for (int j = 0; j < uniques_b.size(); j ++)
        {
            
            // Scores and freqs
            int w_ab = intersection(chain_data_a[i],chain_data_b[j],w_tot);
            int w_a = w_ab - chain_count_a[i]; 
            int w_b = w_ab - chain_count_b[j]; 
            float score = scores[w_ab][w_a][w_b];
            float freq = freqs[w_ab][w_a][w_b];
            
            // Save results that have a score above threshold
            if (score > threshold) 
            {
                stringstream ss;
                ss << score << "\t" << freq << "\t" << uniques_a[i] << "\t" << uniques_b[j];
                results.push_back(ss.str()); 
            }

        } 
        cout << "Finished " << i+1 << "/" << uniques_a.size() << "      \r";
        cout.flush();
    }
    cout << endl;

    // Output results to txt file
    ofstream output_file("./results.txt");
    ostream_iterator<string> output_iterator(output_file, "\n");
    copy(results.begin(), results.end(), output_iterator); 
};




