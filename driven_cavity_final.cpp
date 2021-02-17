#include <iostream>
#include <eigen-3.3.7/Eigen/Sparse>
#include<algorithm> 
#include<fstream>

using namespace Eigen;

int n_x,n_y,i,j;
double dx,dy,Re;
double velocity=1.0;

double l2_norm_x=0.0;
double alpha_uv=1e-4;
double epsilon_uv=0.75;

double l2_norm_y=0.0;

double l2_norm_p=0.0;
double epsilon_p=1e-4;
double alpha_p=0.5;

void link_coefficients(MatrixXd &A_p,MatrixXd &A_e,MatrixXd &A_w,MatrixXd &A_n,MatrixXd &A_s,VectorXd &source_x,VectorXd &source_y,MatrixXd &u_face,MatrixXd &v_face,MatrixXd &p,MatrixXd &u_star,MatrixXd &v_star)

{
    int n;
    double D_e,D_w,D_n,D_s,F_e,F_n,F_s,F_w;
    
    D_e=dy/(dx*Re);
    D_w=dy/(dx*Re);
    D_n=dx/(dy*Re);
    D_s=dx/(dy*Re);
    
    //interior cells
    for(i=2;i<n_y;i++)
    {
        for (j=2;j<n_x;j++)
        {
            n=(i-1)*n_x + (j-1);

            F_e=dy*u_face(i,j);
            F_w=dy*u_face(i,j-1);
            F_n=dx*v_face(i-1,j);
            F_s=dx*v_face(i,j);

            A_e(i,j)=D_e + std::max(0.0,-F_e);
            A_w(i,j)=D_w + std::max(0.0,F_w);
            A_n(i,j)=D_n + std::max(0.0,-F_n);
            A_s(i,j)=D_s + std::max(0.0,F_s);
            A_p(i,j)=D_e + D_w + D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

            source_x[n]=0.5*alpha_uv*(p(i,j-1)-p(i,j+1))*dy + (1-alpha_uv)*A_p(i,j)*u_star(i,j);
            source_y[n]=0.5*alpha_uv*(p(i+1,j)-p(i-1,j))*dx + (1-alpha_uv)*A_p(i,j)*v_star(i,j);
        }
    }

    //left wall
    j=1;
    for (i=2;i<n_y;i++)
    { 
        n=(i-1)*n_x ;

        F_e=dy*u_face(i,j);
        F_w=dy*u_face(i,j - 1);  //left face velocity is initialised to zero
        F_n=dx*v_face(i - 1,j);
        F_s=dx*v_face(i,j);

        A_e(i,j)=D_e + std::max(0.0,-F_e);
        A_n(i,j)=D_n + std::max(0.0,-F_n);
        A_s(i,j)=D_s + std::max(0.0,F_s);
        A_p(i,j)=D_e + 2*D_w + D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

        source_x[n]=0.5*alpha_uv*(p(i,j) - p(i,j + 1))*dy + (1-alpha_uv)*A_p(i,j)*u_star(i,j);  // P_o - 0.5(P_o+P_e)
        source_y[n]=0.5*alpha_uv*(p(i + 1,j) - p(i - 1,j))*dx + (1-alpha_uv)*A_p(i,j)*v_star(i,j);
    }

    //bottom wall
    i=n_y;
    for(j=2;j<n_x;j++)
    {
        n=(n_y-1)*n_x + (j-1);

        F_e=dy*u_face(i,j);
        F_w=dy*u_face(i,j - 1);
        F_n=dx*v_face(i - 1,j);
        F_s=dx*v_face(i,j) ;        //bottom wall v-velocity is already initialised to zero

        A_e(i,j)=D_e + std::max(0.0,-F_e);
        A_w(i,j)=D_w + std::max(0.0,F_w);
        A_n(i,j)=D_n + std::max(0.0,-F_n);
        A_p(i,j)=D_e + D_w + D_n + 2*D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

        source_x[n]=0.5*alpha_uv*(p(i,j - 1) - p(i,j + 1))*dy + (1-alpha_uv)*A_p(i,j)*u_star(i,j);
        source_y[n]=0.5*alpha_uv*(p(i,j) - p(i - 1,j))*dx + (1-alpha_uv)*A_p(i,j)*v_star(i,j);    //P_o - 0.5(P_o+P_n)
    }
    //right wall
    j=n_x;
    for(i=2;i<n_y;i++)
    {
        n=i*n_x -1;    
     
        F_e=dy*u_face(i,j);
        F_w=dy*u_face(i,j - 1);   //right face velocity is initialised to zero
        F_n=dx*v_face(i - 1,j);
        F_s=dx*v_face(i,j);

        A_w(i,j)=D_w + std::max(0.0,F_w);
        A_n(i,j)=D_n + std::max(0.0,-F_n);
        A_s(i,j)=D_s + std::max(0.0,F_s);
        A_p(i,j)=2*D_e + D_w + D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

        source_x[n]=0.5*alpha_uv*(p(i,j-1) - p(i,j))*dy + (1-alpha_uv)*A_p(i,j)*u_star(i,j);   //0.5(P_w+P_o)-P_o
        source_y[n]=0.5*alpha_uv*(p(i + 1,j) - p(i - 1,j))*dx + (1-alpha_uv)*A_p(i,j)*v_star(i,j);
    }
    //top wall
    i=1;
    for(j=2;j<n_x;j++)
    {
        n=(j-1);

        F_e=dy*u_face(i,j);
        F_w=dy*u_face(i,j - 1);
        F_n=dx*v_face(i - 1,j);
        F_s=dx*v_face(i,j);

        A_e(i,j)=D_e + std::max(0.0,-F_e);
        A_w(i,j)=D_w + std::max(0.0,F_w);
        A_s(i,j)=D_s + std::max(0.0,F_s);
        A_p(i,j)=D_e + D_w + 2*D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n) ;

        source_x[n]=0.5*alpha_uv*(p(i,j - 1) - p(i,j + 1))*dy + (1-alpha_uv)*A_p(i,j)*u_star(i,j) + alpha_uv*velocity*(2*D_n + std::max(0.0,-F_n));
        source_y[n]=0.5*alpha_uv*(p(i + 1,j) - p(i,j))*dx  + (1-alpha_uv)*A_p(i,j)*v_star(i,j);   //0.5(P_s+P_o) - P_o
    }
    //top left corner
    i=1;
    j=1;
    n=0;

    F_e=dy*u_face(i,j);
    F_w=dy*u_face(i,j - 1);
    F_n=dx*v_face(i - 1,j);
    F_s=dx*v_face(i,j);

    A_e(i,j)=D_e + std::max(0.0,-F_e);
    A_s(i,j)=D_s + std::max(0.0,F_s);
    A_p(i,j)=D_e + 2*D_w + 2*D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n) ;

    source_x[n]=0.5*alpha_uv*(p(i,j) - p(i,j + 1))*dy + (1-alpha_uv)*A_p(i,j)*u_star(i,j) + alpha_uv*velocity*(2*D_n + std::max(0.0,-F_n));  //P_o - 0.5(P_o+P_e)
    source_y[n]=0.5*alpha_uv*(p(i + 1,j) - p(i,j))*dx + (1-alpha_uv)*A_p(i,j)*v_star(i,j);  //0.5(P_s+P_o) - P_o

    //top right corner
    i=1;
    j=n_x;
    n=n_x-1;

    F_e=dy*u_face(i,j);
    F_w=dy*u_face(i,j - 1);  //right face velocity is initialised to zero
    F_n=dx*v_face(i - 1,j);
    F_s=dx*v_face(i,j);

    A_w(i,j)=D_w + std::max(0.0,F_w);
    A_s(i,j)=D_s + std::max(0.0,F_s);
    A_p(i,j)=2*D_e + D_w + 2*D_n + D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

    source_x[n]=0.5*alpha_uv*(p(i,j - 1) - p(i,j))*dy + (1-alpha_uv)*A_p(i,j)*u_star(i,j) + alpha_uv*velocity*(2*D_n + std::max(0.0,-F_n));  //0.5(P_w+P_o)-P_o
    source_y[n]=0.5*alpha_uv*(p(i + 1,j) - p(i,j))*dx + (1-alpha_uv)*A_p(i,j)*v_star(i,j); //0.5(P_s+P_o) - P_o

    //bottom left corner
    i=n_y;
    j=1;
    n=(n_y-1)*n_x ;


    F_e=dy*u_face(i,j);
    F_w=dy*u_face(i,j - 1) ; //left face velocity is initialised to zero
    F_n=dx*v_face(i - 1,j);
    F_s=dx*v_face(i,j);

    A_e(i,j)=D_e + std::max(0.0,-F_e);
    A_n(i,j)=D_n + std::max(0.0,-F_n);
    A_p(i,j)=D_e + 2*D_w + D_n + 2*D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

    source_x[n]=0.5*alpha_uv*(p(i,j) - p(i,j + 1))*dy + (1-alpha_uv)*A_p(i,j)*u_star(i,j);  //P_o - 0.5(P_o+P_e)
    source_y[n]=0.5*alpha_uv*(p(i,j) - p(i - 1,j))*dx + (1-alpha_uv)*A_p(i,j)*v_star(i,j);  //P_o - 0.5(P_o+P_n)

    //bottom right corner
    i=n_y;
    j=n_x;
    n=n_x*n_y-1;

    F_e=dy*u_face(i,j);
    F_w=dy*u_face(i,j - 1) ; //right face velocity is initialised to zero
    F_n=dx*v_face(i - 1,j);
    F_s=dx*v_face(i,j);

    A_w(i,j)=D_w + std::max(0.0,F_w);
    A_n(i,j)=D_n + std::max(0.0,-F_n);
    A_p(i,j)=2*D_e + D_w + D_n + 2*D_s + std::max(0.0,F_e) + std::max(0.0,-F_w) + std::max(0.0,-F_s) + std::max(0.0,F_n);

    source_x[n]=0.5*alpha_uv*(p(i,j - 1)- p(i,j))*dy + (1-alpha_uv)*A_p(i,j)*u_star(i,j);  //0.5(P_w+P_o)-P_o
    source_y[n]=0.5*alpha_uv*(p(i,j) - p(i-1,j))*dx  + (1-alpha_uv)*A_p(i,j)*v_star(i,j); //P_o - 0.5(P_o+P_n)


    A_e=alpha_uv*A_e;
    A_w=alpha_uv*A_w;
    A_n=alpha_uv*A_n;
    A_s=alpha_uv*A_s;



}

void build_coefficient_matrix(SparseMatrix <double> &A,MatrixXd &ap,MatrixXd &ae,MatrixXd &aw,MatrixXd &an,MatrixXd &as){


   typedef Triplet<double> T;
   std::vector<T> tripletList;
   
   int n;

   //top wall
   int i=1;
   for (int j=2;j<n_x;j++)
   {
        n=(j-1);

        tripletList.push_back(T(n,n,ap(i,j)));
        tripletList.push_back(T(n,n-1,-aw(i,j)));
        tripletList.push_back(T(n,n+1,-ae(i,j)));
        tripletList.push_back(T(n,n+n_x,-as(i,j)));

   }
   
   //interior cells
   for(int i=2;i<n_y;i++)
   {
       for(int j=2;j<n_x;j++)
       {
           n=(i-1)*n_x + (j-1);

           tripletList.push_back(T(n,n,ap(i,j)));
           tripletList.push_back(T(n,n-1,-aw(i,j)));
           tripletList.push_back(T(n,n+1,-ae(i,j)));
           tripletList.push_back(T(n,n-n_x,-an(i,j)));
           tripletList.push_back(T(n,n+n_x,-as(i,j)));

           
        }
    }
    
    
    //left wall
    int j=1;
    for(int i=2;i<n_y;i++)
    {
           n=(i-1)*n_x ;

           tripletList.push_back(T(n,n,ap(i,j)));
           tripletList.push_back(T(n,n+1,-ae(i,j)));
           tripletList.push_back(T(n,n-n_x,-an(i,j)));
           tripletList.push_back(T(n,n+n_x,-as(i,j)));

    }
    
    
    
   //right  wall 
    j=n_x;
    for(int i=2;i<n_y;i++)
    {
        n=i*n_x -1;

        tripletList.push_back(T(n,n,ap(i,j)));
        tripletList.push_back(T(n,n-1,-aw(i,j)));
        tripletList.push_back(T(n,n-n_x,-an(i,j)));
        tripletList.push_back(T(n,n+n_x,-as(i,j)));   
    }
   
   
   //bottom wall
   i =n_y;
   for (int j=2;j<n_x;j++)
   {
       n=(n_y-1)*n_x + (j-1);

       tripletList.push_back(T(n,n,ap(i,j)));
       tripletList.push_back(T(n,n-1,-aw(i,j)));
       tripletList.push_back(T(n,n+1,-ae(i,j)));
       tripletList.push_back(T(n,n-n_x,-an(i,j)));

   }
    
   //top left corner
    n=0;
    i=1;
    j=1;
    tripletList.push_back(T(n,n,ap(i,j)));
    tripletList.push_back(T(n,n+1,-ae(i,j)));
    tripletList.push_back(T(n,n+n_x,-as(i,j)));
     
    
    //top right corner
    i=1;
    j=n_x;
    n=n_x-1;

    tripletList.push_back(T(n,n,ap(i,j)));
    tripletList.push_back(T(n,n-1,-aw(i,j)));
    tripletList.push_back(T(n,n+n_x,-as(i,j)));
    
    
    
    //bottom left corner
    n=(n_y-1)*n_x ;
    j=1;
    i=n_y;

    tripletList.push_back(T(n,n,ap(i,j)));
    tripletList.push_back(T(n,n+1,-ae(i,j)));
    tripletList.push_back(T(n,n-n_x,-an(i,j)));
    
    //bottom right corner
    n=(n_x*n_y)-1;
    i=n_y;
    j=n_x;

    tripletList.push_back(T(n,n,ap(i,j)));
    tripletList.push_back(T(n,n-1,-aw(i,j)));
    tripletList.push_back(T(n,n-n_x,-an(i,j)));
    
    
    A.setFromTriplets(tripletList.begin(),tripletList.end());

}

void solve(SparseMatrix <double> &A,VectorXd &b,MatrixXd &phi,double &l2_norm,double &epsilon)
{
    VectorXd x(n_x*n_y);

    int n=0;
    for(int i=1;i<n_y+1;i++)
    {
        for (int j=1;j<n_x+1;j++)
        {
         
         x[n]=phi(i,j);
         n+=1;

        }
    }
    l2_norm=(A*x-b).norm() ;

    BiCGSTAB<SparseMatrix<double> > solver;
    solver.compute(A);
    solver.setTolerance(epsilon);
    x = solver.solve(b);
    //std::cout << "L2 Norm : " <<  (A*x-b).norm()     << std::endl;
    //std::cout << "Number of iterations : " <<  solver.iterations()    << std::endl;

    n=0;
    for(int i=1;i<n_y+1;i++)
    {
        for (int j=1;j<n_x+1;j++)
        {
         
         phi(i,j)=x[n];
         n+=1;

        }
    }



}

void face_velocity(MatrixXd &u,MatrixXd &v,MatrixXd &u_face,MatrixXd &v_face,MatrixXd &p,MatrixXd &A_p)
{

    //u face velocity
    
    for(i=1;i<n_y+1;i++)
    {
      for (j=1;j<n_x;j++)
      {
        u_face(i,j)=0.5*(u(i,j) + u(i,j + 1)) + 0.25*alpha_uv*(p(i,j + 1) - p(i,j - 1))*dy/A_p(i,j) + 0.25*alpha_uv*(p(i,j + 2) - p(i,j))*dy/A_p(i,j + 1)- 0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(p(i,j + 1) - p(i,j))*dy;
          
      }  
    }



    //v face velocity
    for(i=2;i<n_y+1;i++)
    {
        for (j=1;j<n_x+1;j++)
        {

            v_face(i-1,j)=0.5*(v(i,j) + v(i - 1,j)) + 0.25*alpha_uv*(p(i - 1,j) - p(i + 1,j))*dy/A_p(i,j) + 0.25*alpha_uv*(p(i - 2,j) - p(i,j))*dy/A_p(i - 1,j)- 0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(p(i - 1,j) - p(i,j))*dy;
        }
    }

}

void pressure_correction_link_coefficients(MatrixXd &u,MatrixXd &u_face,MatrixXd &v_face,MatrixXd &Ap_p,MatrixXd &Ap_e,MatrixXd &Ap_w,MatrixXd &Ap_n,MatrixXd &Ap_s,VectorXd &source_p,MatrixXd &A_p)
{
    int n;
    //interior cells
    for(i=2;i<n_y;i++)
    {
        for(j=2;j<n_x;j++)
        {
            n=(i-1)*n_x + (j-1);

            Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
            Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
            Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
            Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
            Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

            source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;

        } 

    }
    
    //top
    i=1;
    for(j=2;j<n_x;j++)
    { 
        n=(j-1);
    
        Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
        Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
        Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
        Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

        source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    }
    
    //left
    j=1;
    for(i=2;i<n_y;i++)
    { 
        n=(i-1)*n_x ;

        Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
        Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
        Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
        Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

        source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    }

    //right
    j=n_x;
    for(i=2;i<n_y;i++)
    {
        n=i*n_x -1;    

        Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
        Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
        Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
        Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

        source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    }
    //bottom
    i=n_y;
    for(j=2;j<n_x;j++)
    { 
        n=(n_y-1)*n_x + (j-1);

        Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
        Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
        Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
        Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

        source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    }

    //top left corner
    i=1;
    j=1;
    n=0;

    Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
    Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
    Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

    source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    
    //top right corner
    i=1;
    j=n_x;
    n=n_x-1;

    Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
    Ap_s(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i + 1,j))*(dx*dx);
    Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

    source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    
    //bottom left corner
    i=n_y;
    j=1;
    n=(n_y-1)*n_x ;

    Ap_e(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j + 1))*(dy*dy);
    Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
    Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

    source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    
    //bottom right corner
    i=n_y;
    j=n_x;
    n=n_x*n_y-1;

    Ap_w(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i,j - 1))*(dy*dy);
    Ap_n(i,j)=0.5*alpha_uv*(1/A_p(i,j) + 1/A_p(i - 1,j))*(dx*dx);
    Ap_p(i,j)=Ap_e(i,j) + Ap_w(i,j) + Ap_n(i,j) + Ap_s(i,j);

    source_p[n]=-(u_face(i,j) - u_face(i,j - 1))*dy - (v_face(i - 1,j) - v_face(i,j))*dx;
    
}

void correct_pressure(MatrixXd &p_star,MatrixXd &p,MatrixXd &p_prime)
{
    
    //top wall
    p_prime.block(0,1,1,n_x)=p_prime.block(1,1,1,n_x);
    
    //left wall
    p_prime.block(1,0,n_y,1)=p_prime.block(1,1,n_y,1);
    
    //right wall
    p_prime.block(1,n_x+1,n_y,1)=p_prime.block(1,n_x,n_y,1);
    
    //bottom wall
    p_prime.block(n_y+1,1,1,n_x)=p_prime.block(n_y,1,1,n_x);
    
    //top left corner
    p_prime(0,0)=(p_prime(1,1)+p_prime(0,1)+p_prime(1,0))/3;

    //top right corner
    p_prime(0,n_x+1)=(p_prime(0,n_x)+p_prime(1,n_x)+p_prime(1,n_x+1))/3;

    //bottom left corner
    p_prime(n_y+1,0)=(p_prime(n_y,0)+p_prime(n_y,1)+p_prime(n_y+1,1))/3;

    //bottom right corner
    p_prime(n_y+1,n_x+1)=(p_prime(n_y,n_x+1)+p_prime(n_y+1,n_x)+p_prime(n_y,n_x))/3;

    MatrixXd p_ref=MatrixXd::Constant(n_y+2,n_x+2,p_prime(0,0));
    
    
    p_star=p+alpha_p*(p_prime);



}

void correct_cell_center_velocity(MatrixXd &u,MatrixXd &v,MatrixXd &u_star,MatrixXd &v_star,MatrixXd &p_prime,MatrixXd &A_p)
{

    //u velocity
    //interior cells
    for(i=1;i<n_y+1;i++) 
    {
        for(j=2;j<n_x;j++) 
        {
        
            u_star(i,j)= u(i,j) + 0.5*alpha_uv*(p_prime(i,j-1)-p_prime(i,j+1))*dy/A_p(i,j);
        }
    }


    //left
    j=1;
    for(i=1;i<n_y+1;i++)
    {
        u_star(i,j)=u(i,j) + 0.5*alpha_uv*(p_prime(i,j) - p_prime(i,j+1))*dy/A_p(i,j);
    }

    //right
    j=n_x;
    for(i=1;i<n_y+1;i++)
    {
        u_star(i,j)=u(i,j) + 0.5*alpha_uv*(p_prime(i,j-1) - p_prime(i,j))*dy/A_p(i,j);
    }

    //v velocity
    for(i=2;i<n_y;i++) 
    {
        for(j=1;j<n_x+1;j++) 
        {
            v_star(i,j)=v(i,j) + 0.5*alpha_uv*(p_prime(i+1,j)-p_prime(i-1,j))*dx/A_p(i,j);
        }
    }


    //top
    i=1;
    for(j=1;j<n_x+1;j++)
    { 
        v_star(i,j)=v(i,j) + 0.5*alpha_uv*(p_prime(i + 1,j) - p_prime(i,j))*dx/A_p(i,j);
    }

    //bottom
    i=n_y;
    for(j=1;j<n_x+1;j++)
    {
        v_star(i,j)=v(i,j) + 0.5*alpha_uv*(p_prime(i,j) - p_prime(i - 1,j))*dx/A_p(i,j);
    }
}

void correct_face_velocity(MatrixXd &u_face,MatrixXd &v_face,MatrixXd &p_prime,MatrixXd &A_p)
{
    for(i=1;i<n_y+1;i++)
    { 
        for(j=1;j<n_x;j++)
        {
            u_face(i,j)=u_face(i,j)+ 0.5*alpha_uv*(1/A_p(i,j)+1/A_p(i,j+1))*(p_prime(i,j)-p_prime(i,j+1))*dy;
        }
    }

    for(i=2;i<n_y+1;i++)
    {
        for(j=1;j<n_x+1;j++) 
        {
            v_face(i-1,j)=v_face(i-1,j) +  0.5*alpha_uv*(1/A_p(i,j)+1/A_p(i-1,j))*(p_prime(i,j)-p_prime(i-1,j))*dx;
        }
    }


}

void post_processing(MatrixXd u_star,MatrixXd v_star,MatrixXd p_star)
{
    VectorXd x(n_x+2),y(n_y+2);
    x << 0,VectorXd::LinSpaced(n_x,dx/2.0,1-dx/2.0),1;
    y << 0,VectorXd::LinSpaced(n_y,dy/2.0,1-dy/2.0),1;

    //Saving outputs for post-processing
     std::ofstream outFile;
     outFile.open("/Users/deep/Desktop/CFMHT final project/u.dat");
     outFile << u_star;
     outFile.close();

     outFile.open("/Users/deep/Desktop/CFMHT final project/v.dat");
     outFile << v_star;
     outFile.close();

     outFile.open("/Users/deep/Desktop/CFMHT final project/p.dat");
     outFile << p_star;
     outFile.close();

     outFile.open("/Users/deep/Desktop/CFMHT final project/x.dat");
     outFile << x;
     outFile.close();

     outFile.open("/Users/deep/Desktop/CFMHT final project/y.dat");
     outFile << y;
     outFile.close();




}

int main ()
{

    std::cout << "Enter number of divisions in x-direction:";
    std::cin >> n_x;
    
    std::cout << "Enter number of divisions in y-direction:";
    std::cin >> n_y;
    
    dx=1.0/n_x;
    dy=1.0/n_y;
    
    std::cout << "Enter Reynolds number:";
    std::cin >> Re;
    
    //Declaring primitive variables
    MatrixXd u(n_y+2,n_x+2),u_star(n_y+2,n_x+2);
    MatrixXd v(n_y+2,n_x+2),v_star(n_y+2,n_x+2);
    MatrixXd p(n_y+2,n_x+2),p_star(n_y+2,n_x+2),p_prime(n_y+2,n_x+2);
    
    //Declaring momentum-link coefficients
    MatrixXd A_p(n_y+2,n_x+2),A_e(n_y+2,n_x+2),A_w(n_y+2,n_x+2),A_n(n_y+2,n_x+2),A_s(n_y+2,n_x+2);
    
    //Declaring pressure-correction link coefficients
    MatrixXd Ap_p(n_y+2,n_x+2),Ap_e(n_y+2,n_x+2),Ap_w(n_y+2,n_x+2),Ap_n(n_y+2,n_x+2),Ap_s(n_y+2,n_x+2);
    
    //Declaring source-terms
    VectorXd source_x(n_y*n_x),source_y(n_y*n_x),source_p(n_y*n_x);
    
    //Declaring face velocities
    MatrixXd u_face(n_y+2,n_x+1),v_face(n_y+1,n_x+2);
    
    //Declaring coefficient matrix and vectors in Ax=b
    SparseMatrix<double> A_momentum(n_x*n_y,n_x*n_y),A_pressure_correction(n_x*n_y,n_x*n_y);
    
    //Boundary conditions    
    u.block(0,1,1,n_x)=MatrixXd::Constant(1,n_x,1.0);
    u_star.block(0,1,1,n_x)=MatrixXd::Constant(1,n_x,1.0);
    u_face.block(0,1,1,n_x-1)=MatrixXd::Constant(1,n_x-1,1.0);
    

    int max_outer_iterations=1000;
    
    for(int n=1;n<=max_outer_iterations;n++)
    {
        link_coefficients(A_p,A_e,A_w,A_n,A_s,source_x,source_y,u_face,v_face,p,u_star,v_star);
    
        build_coefficient_matrix(A_momentum,A_p,A_e,A_w,A_n,A_s);
      
        solve(A_momentum,source_x,u,l2_norm_x,epsilon_uv);
        
        solve(A_momentum,source_y,v,l2_norm_y,epsilon_uv);
        
        face_velocity(u,v,u_face,v_face,p,A_p);
        
        pressure_correction_link_coefficients(u,u_face,v_face,Ap_p,Ap_e,Ap_w,Ap_n,Ap_s,source_p,A_p);

        build_coefficient_matrix(A_pressure_correction,Ap_p,Ap_e,Ap_w,Ap_n,Ap_s);

        solve(A_pressure_correction,source_p,p_prime,l2_norm_p,epsilon_p);
        
        correct_pressure(p_star,p,p_prime);

        correct_cell_center_velocity(u,v,u_star,v_star,p_prime,A_p);

        correct_face_velocity(u_face,v_face,p_prime,A_p);
        
        p=p_star;
    
        std::cout << n << " " << l2_norm_x << " " << l2_norm_y << " " << l2_norm_p << "\n";
        
        if(l2_norm_x < 1e-5 & l2_norm_y < 1e-5 & l2_norm_p < 1e-6)
        { 
        std::cout << "Converged !";
        break;
        }

    }

    post_processing(u_star,v_star,p_star);

}