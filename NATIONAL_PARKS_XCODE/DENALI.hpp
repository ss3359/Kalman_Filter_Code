//
//  DENALI.hpp
//  NATIONAL_PARKS_XCODE
//
//  Created by Owner on 1/28/26.
//

#ifndef DENALI_hpp
#define DENALI_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <algorithm>
#include<random>
#include<chrono>
using namespace std;

const vector<vector<double>> I ={{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},
    {0,0,0,1,0,0}, {0,0,0,0,1,0},{0,0,0,0,0,1}};


// Operations of Matrices

vector<vector<double>> operator*(double c, const vector<vector<double>> &M);
vector<vector<double>> operator*(const vector<vector<double>> &A, const vector<vector<double>> &B);
vector<vector<double>> operator+(const vector<vector<double>> &A, const vector<vector<double>> &B);
vector<vector<double>> operator-(const vector<vector<double>> &A, const vector<vector<double>> &B);

double norm(const vector<vector<double>> &V);

//Statistics Operations:

double Variance(vector<double> data);

vector<vector<double>> Covariance_Matrix(vector<vector<double>> M);



class Matrix{
    
private:
    int rows;
    int columns;
public:
    vector<vector<double>> M;
    Matrix(){}
    Matrix(int r, int c):rows(r),columns(c){
        M =vector<vector<double>> (rows,vector<double>(columns));
    }
    void InputEntries();
    vector<vector<double>> GetSubMatrix(const vector<vector<double>> &M,int rp, int cp);
    double Determinant(const vector<vector<double>>&M_star);
    
    vector<vector<double>> Mij(const vector<vector<double>> M, int rp, int cp);
    vector<vector<double>>  Adjugate(const vector<vector<double>> &M);
    vector<vector<double>>  InverseMatrix(const vector<vector<double>> &M);
    vector<vector<double>>  InverseDiagonalMatrix(const vector<vector<double>> &M);

    vector<vector<double>> Transpose(const vector<vector<double>> &M);
    double Trace(const vector<vector<double>> &M); 
};


class Kalman{
private:
    int m, n;
public:
    Matrix x;
    Matrix A;
    Matrix B;
    Matrix P;
    Matrix H;
    Matrix C;
    Matrix R;
    Matrix Q;
    
    double dt=0.1;
    
    Kalman():m(6),n(0),x(6,1),A(6,6),B(6,3),P(6,6),H(3,6),C(3,6),R(3,3),Q(6,6){}
    
    
    vector<vector<double>> GetX(double x0,double y0, double z0, double vx0=0, double vy0=0, double vz0=0){
        vector<vector<double>> x_vector=x.M;
        
        x_vector={{x0},{y0},{z0},{vx0},{vy0},{vz0}};
        
        cout<<"x.rows: "<<x_vector.size()<<" x_col:"<<x_vector[0].size()<<endl;
        return x_vector;
    }
    vector<vector<double>> GetU(double ax, double ay,double az){
        vector<vector<double>> u_vector={{ax},{ay},{az}};
        
        return u_vector;
    }
    vector<vector<double>> GetMatrixA(){
        vector<vector<double>> A_mat=A.M;
        A_mat={{1,0,0,dt,0,0},
               {0,1,0,0,dt,0},
               {0,0,1,0,0,dt},
               {0,0,0,1,0,0},
               {0,0,0,0,1,0},
               {0,0,0,0,0,1}};
        return A_mat;
    }
    
    vector<vector<double>> GetMatrixB(){
        vector<vector<double>> B_mat=B.M;

        B_mat={{0.5*dt*dt,0,0},
            {0,0.5*dt*dt,0},
            {0,0,0.5*dt*dt},
            {dt,0,0},
            {0,dt,0},
            {0,0,dt}};
        
        return B_mat;
    }
    
    
        vector<vector<double>> GetMatrixP(double x, double y, double z, double vx, double vy, double vz){
            vector<vector<double>> P_mat=P.M;
//            P_mat={{x*x,0,0,0,0,0},
//                    {0,y*y,0,0,0,0},
//                    {0,0,z*z,0,0,0},
//                    {0,0,0,vx*vx,0,0},
//                    {0,0,0,0,vy*vy,0},
//                    {0,0,0,0,0,vz*vz}};
            
            P_mat={{1e6,0,0,0,0,0},
                    {0,1e6,0,0,0,0},
                    {0,0,1e6,0,0,0},
                    {0,0,0,1e3,0,0},
                    {0,0,0,0,1e3,0},
                {0,0,0,0,0,1e3}};
            return P_mat;
        }

    vector<vector<double>> GetMatrixH(){
        vector<vector<double>> H_mat=H.M;
        
        H_mat={{1,0,0,0,0,0},
            {0,1,0,0,0,0},
            {0,0,1,0,0,0}};
        
        return H_mat;
    }
    vector<vector<double>> GetMatrixC(){
        vector<vector<double>> C_mat=C.M;
        
        C_mat={{1,0,0,0,0,0},
            {0,1,0,0,0,0},
            {0,0,1,0,0,0}};
        
        return C_mat;
    }
    
    vector<vector<double>> GetMatrixR(){
        vector<vector<double>> R_mat=R.M;
        
        R_mat={{100,0,0},{0,100,0},{0,0,100}};
        
        return R_mat;
    }
    
    vector<vector<double>> GetMatrixQ(){
        vector<vector<double>> Q_mat=Q.M;
        double sigma_a=0.01;
        
        Q_mat={{pow(dt, 4)/4,0,0,pow(dt, 3)/2,0,0},
               {0,pow(dt, 4)/4,0,0,pow(dt, 3)/2,0},
                {0,0,pow(dt, 4)/4,0,0,pow(dt, 3)/2},
                {pow(dt, 3)/2,0,0,pow(dt, 2),0,0},
                {0,pow(dt, 3)/2,0,0,pow(dt, 2),0},
                {0,0,pow(dt, 3)/2,0,0,pow(dt, 2)}};
        
        return pow(sigma_a,2)*Q_mat;
    }
    
    void KF_Algorithm(){
        Matrix MAT;
        vector<vector<double>> x_true= GetX(100,100,100,0.001,0.001,0.001);
        vector<vector<double>> x_est= x_true;

        vector<vector<double>> u=GetU(2, 2, 2);
        vector<vector<double>> A=GetMatrixA();
        vector<vector<double>> B=GetMatrixB();
        vector<vector<double>> P=GetMatrixP(100, 100, 100, 0.001, 0.001, 0.001);
        vector<vector<double>> H=GetMatrixH();
        vector<vector<double>> C=GetMatrixC();
        vector<vector<double>> R=GetMatrixR();
        vector<vector<double>> Q=GetMatrixQ();

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed);
//            // initialize a uniform distribution between 0 and 1
//            std::uniform_real_distribution<double> unif(0, 100);
        
        //Define the normal distribution
        double mean = 0.0;
        double std_dev = 10.0;
        std::normal_distribution<double> distribution(mean, std_dev);
        
        int n=40;
        cout<<"Initial vector: "<<endl;
        cout<<"<"<<x_true[0][0]<<","<<x_true[1][0]<<","<<x_true[2][0]<<">"<<endl;
        cout<<endl;
        for(int iter=0; iter<n; iter++){
            //Get The Next True Position (Including Noise)
            vector<vector<double>> w = {{(dt*dt/2)*distribution(generator)}, {(dt*dt/2)*distribution(generator)}, {(dt*dt/2)*distribution(generator)},{dt*distribution(generator)},{dt*distribution(generator)},{dt*distribution(generator)}};
            
            vector<vector<double>> v={{distribution(generator)}, {distribution(generator)},{distribution(generator)}}; //veclocity (noise v ~ N(0,R)
            
            //Get the True Position With Noise
            x_true=(A*x_true+B*u)+w;
            vector<vector<double>> z=H*x_true+v;
            
            
            //Get The Next Estimated Position (Without Noise)
            x_est=A*x_est+B*u;
            
            //Get Predicted State Covariance Matrix
            P=(A*P*MAT.Transpose(A))+Q;
            
            //Kalman Filter Update
            //  Update the innovation (suprise):
            vector<vector<double>> y=z-H*x_est;
            
            // Innovation Covariance Matrix:
            vector<vector<double>> S=H*P*MAT.Transpose(H)+R;
            
            // Determine The Kalman Gain:
            vector<vector<double>> K=P*MAT.Transpose(H)*MAT.InverseDiagonalMatrix(S);
            
            // Update the estimated vector
            x_est= x_est+(K*y);
            
            // Update the Covariance:
                P=(I-(K*H))*P*MAT.Transpose(I-(K*H))+K*R*MAT.Transpose(K);
            
           
//            //Print Out The Results:
            cout<<"New Vectors (iter "<<iter<<")"<<endl;
            cout<<"x_true = <"<<x_true[0][0]<<","<<x_true[1][0]<<","<<x_true[2][0]<<">"<<endl;
            cout<<"x_est = <"<<x_est[0][0]<<","<<x_est[1][0]<<","<<x_est[2][0]<<">"<<endl;
            
            cout<<"Error: "<<norm(x_true-x_est)<<endl;
            
            
            cout<<"Trace: "<<MAT.Trace(P)<<endl;
            cout<<endl;
        }
    }
};







#endif /* DENALI_hpp */
