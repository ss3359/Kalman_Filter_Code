//
//  main.cpp
//  NATIONAL_PARKS_XCODE
//
//  Created by Owner on 1/28/26.
//
#include<stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <algorithm>
#include "DENALI.hpp"

using namespace std;

void PrintMatrix(vector<vector<double>> A){
    int m=A.size(), n=A[0].size();
    
    for(int i =0; i<m; i++){
        for(int j=0; j<n; j++){
            cout<<A[i][j]<<"\t";
        }
        cout<<endl;
    }
}

int main(){
//    int r=3, c=3;
//    Matrix Mat1(r,c);
//    
//    Mat1.InputEntries();
//    
//    vector<vector<double>> A=Mat1.M, A_inv=Mat1.InverseMatrix(A);
//    
//    PrintMatrix(A_inv);
    
    Kalman KF;
    KF.KF_Algorithm();
    
    
    return 0;
}
