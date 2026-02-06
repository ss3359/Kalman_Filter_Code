//
//  DENALI.cpp
//  NATIONAL_PARKS_XCODE
//
//  Created by Owner on 1/28/26.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <algorithm>
#include "DENALI.hpp"

using namespace std;

void Matrix:: InputEntries(){
    cout<<"Enter The Values in the Matrix: "<<endl;
    for(int i=0; i<rows; i++){
        for(int j=0; j<columns; j++){
            cin>>M[i][j];
        }
        cout<<endl;
    }
}


vector<vector<double>> Matrix:: Transpose(const vector<vector<double>> &M){
    
    int m=M.size(),n=M[0].size();
    vector<vector<double>> M_T(n, vector<double>(m));
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            M_T[i][j]=M[j][i];
        }
    }
    return M_T;
}


vector<vector<double>> Matrix:: GetSubMatrix(const vector<vector<double>> &M,int rp, int cp){

    int n = M.size();
    vector<vector<double>> M_sub(n-1,vector<double>(n-1));
    
    int r=0;
    for(int i=0; i<n; i++){
        if(i==rp)
            continue;
        int c=0;
        for(int j=0; j<n; j++){
            if(j==cp)
                continue;
            M_sub[r][c]=M[i][j];
            c++;
        }
        r++;
    }
    return M_sub;
}

double Matrix::Determinant(const vector<vector<double>>&M_star){
    int n= M_star.size();
    
   
    if(n==1){
        return M_star[0][0];
    }
    if(n==2){
        return (M_star[0][0]*M_star[1][1])-(M_star[0][1]*M_star[1][0]);
    }
    double det =0.0;
    for(int j=0; j<n; j++){
        auto Sub_Matrix=GetSubMatrix(M_star,0,j);
      
        det+=pow(-1,j)*M_star[0][j]*Determinant(Sub_Matrix);
    }
    return det;
}


vector<vector<double>> Matrix:: Mij(const vector<vector<double>> M, int rp, int cp){
    int n= M.size();
    vector<vector<double>> Mij (n-1,vector<double>(n-1,0.0));
    
    int r=0;
    for(int i=0; i<n;i++){
        if(i==rp)
            continue;
        int c=0;
        for(int j=0; j<n; j++){
            if(j==cp)
                continue;
            Mij[r][c]=M[i][j];
            c++;
        }
        r++;
    }
    return Mij;
}

vector<vector<double>> Matrix::Adjugate(const vector<vector<double>> &M){
    int n=M.size();
    vector<vector<double>> Adj_M(n,vector<double>(n,0.0));

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            vector<vector<double>> mij=Mij(M,i,j);
            double Cji=pow(-1,i+j)*Determinant(mij);
            Adj_M[i][j]=Cji;
        }
    }
    return Transpose(Adj_M);
}

vector<vector<double>> Matrix::  InverseMatrix(const vector<vector<double>> &M){
    if (abs(Determinant(M))<1e-12){
        throw runtime_error("This is a Singular Matrix: det->0");
    }
    
    vector<vector<double>> Adj_M=Adjugate(M);
    
    return (1.0/Determinant(M))*Adj_M;
    
    
}


vector<vector<double>>Matrix::  InverseDiagonalMatrix(const vector<vector<double>> &M){
    
    int m= M.size();
    int n=M[0].size();
    
    vector<vector<double>> R(m, vector<double>(n,0.0));
    
    try{
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                if(i!=j && M[i][j]!=0.0){
                    throw runtime_error("The element is non-zero");
                }
                if(i==j&&abs(M[i][j])<1e-12){
                    throw runtime_error("Zero Diagonal Entry");
                }
                else{
                    R[i][i]= 1.0/M[i][i];
                }
            }
        }
    }catch(exception e){
        cout<<"Error Output"<<endl;
    }
    return R;
}



double Matrix::Trace(const vector<vector<double>> &M){
    
    double result=0.0;
    
    for(int i=0; i<M[0].size(); i++){
        result+=M[i][i];
    }
    return result;
}

//Matrix Operations

vector<vector<double>> operator*(double c, const vector<vector<double>> &M){
    
    int m=M.size();
    int n=M[0].size();
    
    vector<vector<double>> result(m,vector<double>(n,0.0));
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            result[i][j]=c*M[i][j];
        }
    }
    return result;
}


vector<vector<double>> operator*(const vector<vector<double>>& A, const vector<vector<double>> &B){
    int m=A.size(), pA=A[0].size();
    int pB= B.size(), n=B[0].size();

    vector<vector<double>> result(m,vector<double>(n,0.0));

    try{
        
        if(pA != pB){
            throw "Number of Columns and Number of Rows are not equal!";
        }
        
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                for(int k = 0; k<pA; k++){
                    result[i][j]+=A[i][k]*B[k][j];
                }
            }
        }
        
    }catch(runtime_error e){
        cout<<"Error! "<<e.what()<<endl;
    }
    return result;
}

vector<vector<double>> operator+(const vector<vector<double>> &A, const vector<vector<double>> &B){
    
    int m=A.size(),n=A[0].size();
    vector<vector<double>> result(m,vector<double>(n,0.0));
    
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

vector<vector<double>> operator-(const vector<vector<double>> &A, const vector<vector<double>> &B){
    
    int m=A.size(),n=A[0].size();
    vector<vector<double>> result(m,vector<double>(n,0.0));
    
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}


double norm(const vector<vector<double>> &V){
    double result=0.0;
    
    int m=V.size(),n=V[0].size();
    try{
        if(n>1){
            throw runtime_error("This is a Matrix, not a vector");
        }
        for(int i=0; i<m; i++){
            result+=V[i][0]*V[i][0];
        }
        
    }catch(exception e){
        cout<<"The vector has more than 1 columns"<<endl;
    }
    return sqrt(result);

}



//Statistics Operations:

double Variance(vector<double> data){
    int n= data.size();
    double v=0,sum=0,avg=0;
    
    for(int i=0; i<n; i++){
        sum+=(data[i]);
    }
    avg=(1/n)*sum;
    
    for(int i=0; i<n; i++){
        v+=pow((avg-data[i]),2);
    }
    return sqrt((1/n)*v);
}


vector<vector<double>> Covariance_Matrix(vector<vector<double>> M){
    
    int m=M.size(); // Number of Equations
    int n=M[0].size(); // Number of Unknowns (Variables);
    vector<double> Variances(n,0.0);
    vector<double> data_j(m,0.0);
    vector<vector<double>> COV(m,vector<double>(n,0.0));
    
    for(int j=0; j<n; j++){
        for(int i=0; i<m; i++){
            data_j[i]=M[i][j];
        }
        double sigma_j=Variance(data_j);
        Variances[j]= sigma_j;
    }
    
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            COV[i][j]=Variances[i]*Variances[j];
        }
    }
    return COV;
}


