#include "feMetric.h"
#include "fullMatrix.h"
#include "feSimplex.h"

#include <assert.h>
#include <iostream>
#include <fstream>

void metricHechtKuate(int nbpoints, double* x, double* y, double &A, double &B, double &C, double epsilon, double* xNew, double* yNew);

feMetric::feMetric(feRecovery *recovery, feMetricOptions metricOptions) : _recovery(recovery), _options(metricOptions)
{	
  int degSol = recovery->getDegreeSolution();
  int dimRecovery = recovery->getDimRecovery();
  std::vector<int> &expX = recovery->getXExponentsRecovery(); 
  std::vector<int> &expY = recovery->getYExponentsRecovery();

  std::vector<int> &vertices = recovery->getVertices();

  std::map<int, std::vector<double>> &e = recovery->getErrorCoefficients();

  int nPhi = _options.nPhi;
  std::vector<double> phi(nPhi, 0.); // Error curve discretization
  for(int i = 0; i < nPhi; ++i){ phi[i] = i * 2.*M_PI/nPhi; }

  double* xC   = (double*) malloc(sizeof(double) * nPhi);
  double* yC   = (double*) malloc(sizeof(double) * nPhi);
  double* xNew = (double*) malloc(sizeof(double) * nPhi);
  double* yNew = (double*) malloc(sizeof(double) * nPhi);

  double err;
  int cnt;
  for(auto v : vertices){


    // printf("error at vertex %d (%3.4f,%3.4f,%3.4f) = %+-3.4f - %+-3.4f - %+-3.4f - %+-3.4f\n", v,
    //         recovery->_mesh->getVertex(v)->x(),
    //         recovery->_mesh->getVertex(v)->y(),
    //         recovery->_mesh->getVertex(v)->z(),
    //         e[v][0],
    //         e[v][1],
    //         e[v][2],
    //         e[v][3]);




    double A, B, C;
    // Construct the discretization of the error curve
    for(int iPhi = 0; iPhi < nPhi; ++iPhi){
      err = 0.0;
      cnt = 0;
      for(int n = 0; n < dimRecovery; ++n){
        if(expX[n] + expY[n] == degSol+1){ // Select only exponents of the homogeneous poynomial of order deg+1
          // std::cout<<"pos = "<<cnt<<" : "<<e[v][cnt]<<std::endl;
          err += e[v][cnt++] * pow(cos(phi[iPhi]), expX[n]) * pow(sin(phi[iPhi]), expY[n]);
        }
      }
      err = fmax(fabs(err), 1e-8);
      // std::cout<<"phi =  "<<phi[iPhi]<<std::endl;
      // std::cout<<"cos phi =  "<<cos(phi[iPhi])<<std::endl;
      // std::cout<<"sin phi =  "<<sin(phi[iPhi])<<std::endl;
      // std::cout<<"1 : "<<pow(1.0/fabs(err), 1.0/(degSol+1))* cos(phi[iPhi])<<std::endl;
      // std::cout<<"2 : "<<pow(1.0/fabs(err), 1.0/(degSol+1))* sin(phi[iPhi])<<std::endl;
      // std::cout<<"3 : "<<pow(1.0/fabs(err), 1.0/(degSol+1))<<std::endl;
      xC[iPhi] = pow(1.0/fabs(err), 1.0/(degSol+1)) * cos(phi[iPhi]);
      yC[iPhi] = pow(1.0/fabs(err), 1.0/(degSol+1)) * sin(phi[iPhi]);

    }

    // for(int i = 0; i < nPhi; ++i){
    //   std::cout<<"xC : "<<xC[i]<<std::endl;
    //   std::cout<<"yC : "<<yC[i]<<std::endl;
    // }
    // if(false){
      // Compute brute-force metric
      metricHechtKuate(nPhi, xC, yC, A, B, C, 1e-5, xNew, yNew);
      // std::cout<<"Computed metric : "<<A<<" - "<<B<<" - "<<C<<std::endl;
      SMetric3 M;
      M.set_m11(A);
      M.set_m21(C);
      M.set_m22(B);
       // M.print("Mavant");
      // TAKE ABS of eigenvalues
      fullMatrix<double> V(3, 3);
      fullVector<double> S(3);
      M.eig(V, S, false);
      // std::cout<<"Valeurs propres : "<<S(0)<<" - "<<S(1)<<" - "<<S(2)<<std::endl;
      SVector3 v0(V(0,0), V(0,1), V(0,2));  // Attention c'est peut-être la transposée
      SVector3 v1(V(1,0), V(1,1), V(1,2));
      SVector3 v2(V(2,0), V(2,1), V(2,2));
      // std::cout<<"Vecteurs propres : "<<std::endl;
      // v0.print();
      // v1.print();
      // v2.print();
      M = SMetric3(fabs(S(0)), fabs(S(1)), fabs(S(2)), v0, v1, v2);
      // M.print("M");
      // M = SMetric3();
      // M.set_m11(100);

      _metrics[v] = M;
    // } else{
      // metriqueSimplexe2D(nPhi, phi, e[v], A, B, C, 50);
      // std::cout<<"Computed metric simplexe : "<<A<<" - "<<B<<" - "<<C<<std::endl;
      // SMetric3 M;
      // M.set_m11(A);
      // M.set_m21(C);
      // M.set_m22(B);
      //  M.print("Mavant");
      // // TAKE ABS of eigenvalues
      // fullMatrix<double> V(3, 3);
      // fullVector<double> S(3);
      // M.eig(V, S, false);
      // // std::cout<<"Valeurs propres : "<<S(0)<<" - "<<S(1)<<" - "<<S(2)<<std::endl;
      // SVector3 v0(V(0,0), V(0,1), V(0,2));  // Attention c'est peut-être la transposée
      // SVector3 v1(V(1,0), V(1,1), V(1,2));
      // SVector3 v2(V(2,0), V(2,1), V(2,2));
      // // std::cout<<"Vecteurs propres : "<<std::endl;
      // // v0.print();
      // // v1.print();
      // // v2.print();
      // M = SMetric3(fabs(S(0)), fabs(S(1)), fabs(S(2)), v0, v1, v2);
      // M.print("M");
      // // M = SMetric3();
      // // M.set_m11(100);

      // _metrics[v] = M;
    // }
  }

  metricScaling();

  free(xNew);
  free(yNew);
  free(xC);
  free(yC);

}

void feMetric::metricScaling(){


  std::vector<double> &w = _recovery->_geoSpace->getQuadratureWeights();
  std::vector<double> &J = _recovery->_cnc->getJacobians();

  int nQuad = w.size();
  double N = (double) _options.nTargetVertices;
  double p = _options.LpNorm;
  double deg = (double) _recovery->getDegreeSolution();
  double dim = (double) _recovery->getDim();
  double exponent = p*(deg + 1.0) / (2.0 * (p*(deg + 1.0) + dim));

  // std::cout<<"exponent = "<<exponent<<std::endl;

  double I = 0.0;
  for(int iElm = 0; iElm < _recovery->_nElm; ++iElm){
    for(int k = 0; k < nQuad; ++k){
      // Interpolate det(Q) at quad nodes
      for(int iNode = 0; iNode < _recovery->_nNodePerElm; ++iNode){
        int v = _recovery->_cnc->getNodeConnectivity(iElm,iNode);
        I += J[nQuad*iElm+k] * w[k] * _recovery->_geoSpace->getFunctionAtQuadNode(iNode,k) * pow(_metrics[v].determinant(), exponent);   
      }
    }
  }

  // std::cout<<"Int det(Q) = "<<I<<std::endl;

  std::vector<int> &vertices = _recovery->getVertices();

  double hMin = _options.hMin, lMax = 1.0/(hMin*hMin);
  double hMax = _options.hMax, lMin = 1.0/(hMax*hMax);

  fullMatrix<double> V(3, 3);
  fullVector<double> S(3);
  SVector3 v0, v1, v2;
  for(auto v : vertices){

    SMetric3 M = _metrics[v];

    double factor = pow(N/I, 2.0/dim) * pow(M.determinant(), -1.0/(p*(deg + 1.0) + dim));

    M *= factor;

    M.eig(V, S, false);
    v0 = SVector3(V(0,0), V(0,1), V(0,2));
    v1 = SVector3(V(1,0), V(1,1), V(1,2));
    v2 = SVector3(V(2,0), V(2,1), V(2,2));
    M = SMetric3( fmin( lMax, fmax( lMin, S(0))),
                  fmin( lMax, fmax( lMin, S(1))),
                  fmin( lMax, fmax( lMin, S(2))), v0, v1, v2);

    _metrics[v] = M;
  }
}

void metricHechtKuate(int nbpoints, double* x, double* y, double &A, double &B, double &C, double epsilon, double* xNew, double* yNew)
{
  C = 0.0;

  int bool_assert = 1;
 
  double epsilon0 = 1e-5, precision = 1e-18, delta = 1e-10;
  double inf = 1e100;

  double Rmin = 1e20, Rmax = 0;

  int indiceX0 = 0;

  // double* xNew = new double[nbpoints];
  // double* yNew = new double[nbpoints];
  // double *PPoint = new double[nbpoints];

  for(int i=0; i<nbpoints; i++){
      Rmax = fmax(Rmax,sqrt(x[i]*x[i] + y[i]*y[i]));
    
      //---déplacement des points situées sur les axes--------------
      if(abs(x[i])<=precision){
          if(y[i]<0){
              x[i]=-delta;
              y[i]=-sqrt(pow(y[i],2)-pow(x[i],2));
          }
          if(y[i]>0){
              x[i]=delta;
              y[i]=sqrt(pow(y[i],2)-pow(x[i],2));
          }
      }
    
      if(abs(y[i])<=precision){
          if(x[i]<0){
              y[i]=-delta;
              x[i]=-sqrt(pow(x[i],2)-pow(y[i],2));
          }
          if(x[i]>0){
              y[i]=delta;
              x[i]=sqrt(pow(x[i],2)-pow(y[i],2));
          }
      }
      //-----------------------------------------------------------
      if(bool_assert){
        // std::cout<<"foo : "<<abs(x[i]*y[i])<<std::endl;
        // std::cout<<"foo : "<<pow(precision,2)<<std::endl;
          assert(abs(x[i]*y[i])>=pow(precision,2));
      }

      if(Rmin > sqrt(x[i]*x[i] + y[i]*y[i])){
          indiceX0 = i;
          Rmin = sqrt(x[i]*x[i] + y[i]*y[i]);
      }
  }

  //-------permutation des indices de la liste des points : 
  //ranger la liste en commençant par le point X0-------
  for(int k=0;k<nbpoints-indiceX0;k++){
      xNew[k] = x[k+indiceX0];
      yNew[k] = y[k+indiceX0];
  }
  for(int k=nbpoints-indiceX0;k<nbpoints;k++){
      xNew[k] = x[k-nbpoints+indiceX0];
      yNew[k] = y[k-nbpoints+indiceX0];
  }
  for(int i=0;i<nbpoints;i++){
      x[i] = xNew[i];
      y[i] = yNew[i];
  }

  //----------------------------------------------------------------

  int test = -1;

  double X0, Y0;
  double bmin= 0.,bmax=inf,b1,b2,aik=0.,bik=0.,cik=0.;
  double Xk=0.,Yk=0.,Ck=0.,Ak =0.,Bk=0.,Xi=0.,Yi=0.,ri,detXY=0.,Ri,R0,r0;


  X0 = x[0];
  Y0 = y[0];
  r0 = sqrt(x[0]*x[0] + y[0]*y[0]);
  if(bool_assert){
      // std::cout<<"Assert ON"<<std::endl;
      assert(r0 == Rmin);
  }

  //std::cout<<" Rmin = "<<Rmin<<" Rmax =  "<<Rmax<<std::endl;

  double EPS = 0.0;// pour recuperer la valeur de epsilon0 optimale
  double epsilonmax = r0*(1.-r0/Rmax)/20.;

  // double Tabepsilon[20];

  // int neps =4;
  // //--------- discertisation de epsilon0----------------------------------
  // if(epsilonmax>1e-2){
  //     neps=10;

  //     Tabepsilon[0]=1e-5;
  //     Tabepsilon[1]=1e-4;
  //     Tabepsilon[2]=1e-3;
  //     for(int i=3;i<neps;i++){
  //         Tabepsilon[i]=(i-3)*(epsilonmax-1e-2)/(neps-4.) +1e-2;
  //     }
  // }
  // else{
  //     Tabepsilon[0]=1e-5;
  //     Tabepsilon[1]=1e-4;
  //     Tabepsilon[2]=1e-3;
  //     Tabepsilon[3]=1e-2;
  // }
  //------------------------------------------------------------------------

  int condition = -1;

  if(r0 <= epsilon0)
      epsilon0 = r0*epsilon0;
  
  A = 1./((r0-epsilon0)*(r0-epsilon0));
  B = A;
  
  double epsilon0min = epsilon0;

  // if(abs(Rmin-Rmax) > 1e-5){
  if(abs(Rmin-Rmax) > precision){
      // for(int ee=0; ee<neps-1; ee++){ //boucle sur epsilon0---------------
          // epsilon0= Tabepsilon[ee];
          if(r0<=epsilon0) epsilon0 = r0*epsilon0;
          if(bool_assert){
              assert(r0>epsilon0);
          }
          R0=r0/(r0-epsilon0);
   
          for(int i=1;i<nbpoints;i++){     //boucle sur chaque noeud             
              Xi = x[i];
              Yi = y[i];
              ri = sqrt(x[i]*x[i] + y[i]*y[i]);

              if(ri <= epsilon)
                  epsilon = ri*epsilon;

              if(bool_assert){
                  assert(ri > epsilon);
              }

              Ri = ri/(ri-epsilon);

              detXY = Xi*Y0 - Yi*X0;



              //------deplacement des points alignés avec l'origine et X0-----------
              if(abs(detXY) <= precision){
                  printf("Point %d - \t x = %10.15f - \t y =%10.15f - \t x0 = %10.15f - \t y0 =%10.15f - \t detXY = %10.15f \n", i, x[i],y[i],X0,Y0,detXY);
                  Xi += delta;
           
                  if(Yi<0) Yi= -sqrt(pow(ri,2)- pow(Xi,2));
                  else Yi = sqrt(pow(ri,2)- pow(Xi,2));
                  x[i]=Xi;
                  y[i]=Yi;

                  std::cout<<ri<<std::endl;
                  ri=sqrt(x[i]*x[i] + y[i]*y[i]);

                  if(ri<=epsilon)
                      epsilon=ri*epsilon;

                  std::cout <<"xi ="<<x[i]<<" yi ="<<y[i]<<" ri = "<<ri<<" epsilon = "<<epsilon<<std::endl;
                  if(bool_assert){
                      assert(ri>epsilon);
                  }
                  Ri = ri/(ri-epsilon);
              }
           
              detXY=Xi*Y0-Yi*X0;
           
              if(bool_assert){
                  assert(abs(detXY)>=precision);   
              }
           
              //-----racines du polynome en b à minimiser----------------------------
              double bb1=(1./pow(detXY,2))*(pow(X0*Ri,2)+pow(Xi*R0,2)-2.*abs(Xi*X0)*sqrt(pow(R0*Ri,2)-pow(detXY/(Rmax*(r0-epsilon0)),2)));
              double bb2=(1./pow(detXY,2))*(pow(X0*Ri,2)+pow(Xi*R0,2)+2.*abs(Xi*X0)*sqrt(pow(R0*Ri,2)-pow(detXY/(Rmax*(r0-epsilon0)),2)));
              //--fin----racines du polynome en b à minimiser--------------------
           
              bmax = fmin(bb2, pow(Rmax/pow((r0),2),2) );
              bmin = fmax(1./(Rmax*Rmax),bb1);//minoration de b
              double Cte = fmax(1e-9,(bmax-bmin)*1e-9);
              bmin=bmin*(1.+Cte);
              bmax=bmax*(1.-Cte);

              //bornes de b-----------------------------------------------------------
           
              //cas:  majoration de c --------------------------------------------
              double Li=X0*Xi*(pow(Rmax/pow(r0-epsilon0min,2),2)-1./pow(Rmax,2))+(pow(Ri*X0,2)-pow(R0*Xi,2))/detXY;
              double LiXY=Xi*Y0+Yi*X0;
           
              if(abs(LiXY)>=precision){
                  condition=1;
           
                  if(Xi*X0>0){
                      if(LiXY>0) bmin=fmax(bmin,-Li/LiXY);
                      else bmax=fmin(bmax,-Li/LiXY);
                  }
                  else{
                      if(LiXY<0) bmin=fmax(bmin,-Li/LiXY);
                      else bmax=fmin(bmax,-Li/LiXY); 
                  }
              }
              else
              {
                  if(Li<0) condition = 0;
                  else condition =1;
              }

              //cas  minoration de c --------------------------------------------
              Li=X0*Xi*(-pow(Rmax/pow(r0-epsilon0min,2),2)+1./pow(Rmax,2))+(pow(Ri*X0,2)-pow(R0*Xi,2))/detXY;
              LiXY=Xi*Y0+Yi*X0;
           
              if(abs(LiXY)>=precision){
                  condition=1;
                  if(Xi*X0>0){
                      if(LiXY<0) bmin=fmax(bmin,-Li/LiXY);
                      else bmax=fmin(bmax,-Li/LiXY);
                  }
                  else{
                      if(LiXY>0) bmin=fmax(bmin,-Li/LiXY);
                      else bmax=fmin(bmax,-Li/LiXY); 
                  }
              }
              else{
                  if(Li>0) condition =0;
                  else condition =1;
              }
           
              if(condition){
           
                  //--cas : minoration de a-----------------------------------------------
               
                  double Gi=((Xi*Yi*R0*R0-X0*Y0*Ri*Ri)/detXY +Xi*X0/(Rmax*Rmax))/(Yi*Y0);
               
                  if(Xi*X0>0){
                      if(Yi*Y0>0) bmin=fmax(bmin,Gi);
                      else bmax=fmin(bmax,Gi);
                  }
                  else{
                      if(Yi*Y0<0) bmin=fmax(bmin,Gi);
                      else bmax=fmin(bmax,Gi);
                  }
               
                  //cas :majoration de a------------------------------------------------
                  double Hi=(Xi*X0*Rmax*Rmax/pow((r0-epsilon0min),4)+(Xi*Yi*R0*R0-X0*Y0*Ri*Ri)/detXY)/(Yi*Y0);
                  if(Xi*X0>0){
                      if(Yi*Y0>0) bmax=fmin(bmax,Hi);
                      else bmin=fmax(bmin,Hi);
                  }
                  else{
                      if(Yi*Y0<0) bmax=fmin(bmax,Hi);
                      else bmin=fmax(bmin,Hi); 
                  }
                  //------fin bornes de b------------------------------------------------
                  b2=bmax;
                  b1=bmin;
               
                  for(int k=1; k<nbpoints ;k++){   //on balaye les contraintes
                      Xk = x[k];
                      Yk = y[k];
                      Bk = (Yk*Yk*Xi*X0 +Xk*(Xk*Yi*Y0-Yk*(Yi*X0+Xi*Y0)))/(Xi*X0);
                      Ck = (X0*Xi*detXY-Xk*(Xi*R0*R0*(Yk*Xi-Yi*Xk) +X0*Ri*Ri*(-Yk*X0+Y0*Xk)))/(Xi*X0*detXY);
                   
                      if(bool_assert){
                          assert(abs(Xi*X0*Y0*Yi*Xk*Yk)>=pow(precision,5));
                      }
                      if(abs(Bk)>precision){  //non nul
                     
                          if(Bk<=0) bmax=fmin(bmax,Ck/Bk);
                          else  bmin=fmax(bmin,Ck/Bk);
                   
                          if((bmax<b1)||(bmin>b2)||(bmin>bmax)){
                              test=0;
                              break;  
                          }
                   
                          else
                              test=1;              
                      }
                      else{
                          if(Ck>precision){ 
                              test=0;
                              break;      
                          }
                          else //Ck<=0
                              test=-1;// 1 peut etre
                     }
                  }

                  if(test){ 
                      double a0= -pow((detXY/(Xi*X0)),2);
                      double a1= 2.*(pow(Ri/Xi,2)+pow(R0/X0,2));
                      if(((a0*bmax+a1)*bmax) < ((a0*bmin+a1)*bmin)) 
                          bik = bmax;
                      else 
                          bik = bmin;

                      aik=(Ri*Ri*Y0*X0 -R0*R0*Yi*Xi+bik*Yi*Y0*detXY)/(detXY*Xi*X0);
                      cik=( -Ri*Ri*X0*X0 + R0*R0*Xi*Xi-bik*(Yi*X0+Y0*Xi)*detXY)/(detXY*Xi*X0);
                  
                      if(bool_assert){
                          assert((4.*aik*bik-cik*cik)>=0.);// aire positive   
                          assert(abs((4.*aik*bik-cik*cik)-pow(2./(Rmax*(r0-epsilon0)),2))>0);// aire positive
                      }
                      if((4.*aik*bik-cik*cik) <= (4.*A*B-C*C)){
                          A=aik;
                          B=bik;
                          C=cik;
                          EPS=epsilon0;
                      }
                  } //if(test)
              } // if(condition)
          } // for(int i=1;i<nbpoints;i++)
      // } // for(int ee=0; ee<neps-1; ee++)
  } // if(abs(Rmin-Rmax)>1e-5)
  else{
      A = 1./(Rmin*Rmin);
      B = A;
      C = 0.;
  }
}

void feMetric::metriqueSimplexe2D(int nPhi, std::vector<double> phi, std::vector<double> erreur, double &A, double &B, double &C, int max_iter){
    
    double xi, yi, xj, yj, l1, l2, v11, v12, v21, v22, Q11, Q12, Q21, Q22;
    double err, normeXj;
    double L1, L2, L3, expL11, expL12, expL21, expL22;
    double Aprev, Bprev, Cprev;
    // Variables LAPACK
    char jobz='v', uplo='l';
    const int n=2, lwork=18;
    int info;
    double w[2], work[18];

    SVector3 v0, v1, v2;

    std::vector<int> &expX = _recovery->getXExponentsRecovery(); 
    std::vector<int> &expY = _recovery->getYExponentsRecovery();

    int nSimplex = 6, mSimplex = nPhi;

    std::vector<double> AConstraint(nPhi*nSimplex), bConstraint(nPhi), cObjective(nSimplex);
    cObjective[0] = -1.0;
    cObjective[1] =  1.0;
    cObjective[2] =  0.0;
    cObjective[3] =  0.0;
    cObjective[4] = -1.0;
    cObjective[5] =  1.0;

    A = 1.0;    // Candidat initial : Q = I
    C = 0.0;
    B = 1.0;

    std::pair<std::vector<double>, double> retSimplex;

    int i = 0;
    double residu = 1.0;
    while(i < max_iter && residu > 1e-8){
        Aprev = A;
        Bprev = B;
        Cprev = C;
        //==================================================================================
        // Calcul de Q^(-1/2)
        SMetric3 Q(1.0);
        Q.set_m11(A);
        Q.set_m21(C);
        Q.set_m22(B);
        fullMatrix<double> V(3, 3);
        fullVector<double> S(3);
        Q.eig(V, S, false);

        // S(0) = pow(S(0), -0.5);
        // S(1) = pow(S(1), -0.5);
        // S(2) = pow(S(2), -0.5);

        v0 = SVector3(V(0,0), V(0,1), V(0,2));  // Attention c'est peut-être la transposée
        v1 = SVector3(V(1,0), V(1,1), V(1,2));
        v2 = SVector3(V(2,0), V(2,1), V(2,2));

        Q = SMetric3(pow(S(0), -0.5), pow(S(1), -0.5), pow(S(2), -0.5), v0, v1, v2);
        
        // double aTest[4] = {A, C, C, B};
        // dsyev(&jobz, &uplo, &n, aTest, &n, w, work, &lwork, &info);

        // l1 = w[0];  v11 = aTest[0];  v21 = aTest[2];
        // l2 = w[1];  v12 = aTest[1];  v22 = aTest[3];

        // Q11 = v11*v11*pow(l1, -0.5) + v21*v21*pow(l2, -0.5);
        // Q12 = v12*v11*pow(l1, -0.5) + v21*v22*pow(l2, -0.5);
        // Q21 = v11*v12*pow(l1, -0.5) + v21*v22*pow(l2, -0.5);
        // Q22 = v12*v12*pow(l1, -0.5) + v22*v22*pow(l2, -0.5);
        //==================================================================================
        for(int iPhi = 0; iPhi < nPhi; ++iPhi){

            xi = cos(phi[iPhi]);
            yi = sin(phi[iPhi]);
            
            xj = Q(0,0)*xi + Q(0,1)*yi;
            yj = Q(1,0)*xi + Q(1,1)*yi;

            // err = erreur[0] * xj*xj + erreur[1] * xj*yj + erreur[2] * yj*yj;
            err = 0.0;
            int indice = 0;
            for(int iDeg = 0; iDeg < _recovery->getDimRecovery(); ++iDeg){
                if((expX[iDeg] + expY[iDeg]) == _recovery->getDegreeSolution()+1){
                    err += erreur[indice] * pow(xj, expX[iDeg]) * pow(yj, expY[iDeg]);
                    ++indice;
                }
            }

            xj = xi / fmax( pow( fabs(err) , 1.0/(_recovery->getDegreeSolution()+1)) ,1e-20);
            yj = yi / fmax( pow( fabs(err) , 1.0/(_recovery->getDegreeSolution()+1)) ,1e-20);

            // ==============================================================================
            // Contraintes lineaires
            AConstraint[iPhi*nSimplex  ] =        xj*xj;
            AConstraint[iPhi*nSimplex+1] =       -xj*xj;
            AConstraint[iPhi*nSimplex+2] =  2.0 * xj*yj;
            AConstraint[iPhi*nSimplex+3] = -2.0 * xj*yj;
            AConstraint[iPhi*nSimplex+4] =        yj*yj;
            AConstraint[iPhi*nSimplex+5] =       -yj*yj;

            normeXj = xj*xj + yj*yj;

            bConstraint[iPhi] = -(-normeXj * log(normeXj));
        }

        retSimplex = simplex(nSimplex,mSimplex,AConstraint,bConstraint,cObjective,0.0);
            
        if(isinf(retSimplex.second)){
            if (retSimplex.first[0] == -1) printf("Objective function unbounded!\n");
            else if (retSimplex.first[0] == -2) printf("Linear program infeasible!\n");
        }else{
            // printf("Solution: (");
            // for (int i=0;i<nSimplex+mSimplex;i++) printf("%lf%s", retSimplex.first[i], (i < nSimplex + mSimplex - 1) ? ", " : ")\n");
            // for (int i=0;i<nSimplex;i++) printf("%lf%s", retSimplex.first[i], (i < nSimplex - 1) ? ", " : ")\n");
            // printf("Optimal objective value: %lf\n", retSimplex.second);
            L1 = retSimplex.first[0]-retSimplex.first[1];
            L2 = retSimplex.first[2]-retSimplex.first[3];
            L3 = retSimplex.first[4]-retSimplex.first[5];
        }

        //==================================================================================
        // Calcul de Q^(1/2)
        // Q11 = v11*v11*pow(l1, 0.5) + v21*v21*pow(l2, 0.5);
        // Q12 = v12*v11*pow(l1, 0.5) + v21*v22*pow(l2, 0.5);
        // Q21 = v11*v12*pow(l1, 0.5) + v21*v22*pow(l2, 0.5);
        // Q22 = v12*v12*pow(l1, 0.5) + v22*v22*pow(l2, 0.5);
        Q = SMetric3(pow(S(0), 0.5), pow(S(1), 0.5), pow(S(2), 0.5), v0, v1, v2);
        Q11 = Q(0,0);
        Q12 = Q(0,1);
        Q21 = Q(1,0);
        Q22 = Q(1,1);
        //==================================================================================
        // Calcul de expm(L)

        SMetric3 L(1.0);
        L.set_m11(L1);
        L.set_m21(L2);
        L.set_m22(L3);
        // fullMatrix<double> V(3, 3);
        // fullVector<double> S(3);
        L.eig(V, S, false);

        // S(0) = pow(S(0), -0.5);
        // S(1) = pow(S(1), -0.5);
        // S(2) = pow(S(2), -0.5);

        v0 = SVector3(V(0,0), V(0,1), V(0,2));  // Attention c'est peut-être la transposée
        v1 = SVector3(V(1,0), V(1,1), V(1,2));
        v2 = SVector3(V(2,0), V(2,1), V(2,2));

        L = SMetric3(exp(S(0)), exp(S(1)), exp(S(2)), v0, v1, v2);

        expL11 = L(0,0);
        expL12 = L(0,1);
        expL21 = L(1,0);
        expL22 = L(1,1);

        // double aTest2[4] = {L1, L2, L2, L3};
        // dsyev(&jobz, &uplo, &n, aTest2, &n, w, work, &lwork, &info);

        // l1 = w[0];  v11 = aTest2[0];  v21 = aTest2[2];
        // l2 = w[1];  v12 = aTest2[1];  v22 = aTest2[3];

        // expL11 = v11*v11*exp(l1) + v21*v21*exp(l2);
        // expL12 = v12*v11*exp(l1) + v21*v22*exp(l2);
        // expL21 = v11*v12*exp(l1) + v21*v22*exp(l2);
        // expL22 = v12*v12*exp(l1) + v22*v22*exp(l2);
        //==================================================================================
        // Calcul de Q = Q^(1/2) * expm(L) * Q^(1/2)
        A = Q11*(Q11*expL11 + Q12*expL21) + Q21*(Q11*expL12 + Q12*expL22);
        C = Q12*(Q11*expL11 + Q12*expL21) + Q22*(Q11*expL12 + Q12*expL22);
        B = Q12*(Q21*expL11 + Q22*expL21) + Q22*(Q21*expL12 + Q22*expL22);

        residu = sqrt((A-Aprev)*(A-Aprev) + 2*(C-Cprev)*(C-Cprev) + (B-Bprev)*(B-Bprev));
        // cout<<"======== ITERATION "<<i+1<<" - RESIDU = "<<residu<<" ==========="<<endl;
        ++i;
    }
}

void feMetric::writeSizeFieldSol(std::string solFileName){
  // Write the size field to a .sol file
  int dim = _recovery->getDim();
  std::vector<int> &vertices = _recovery->getVertices();

  FILE* myfile = fopen(solFileName.c_str(),"w");
  fprintf(myfile, "MeshVersionFormatted 2\n\n");
  fprintf(myfile, "Dimension 3\n\n");
  fprintf(myfile, "SolAtVertices\n");
  fprintf(myfile, "%ld\n", _recovery->getVertices().size());
  fprintf(myfile, "1 3\n\n");

  for(auto v : vertices){
    SMetric3 M = _metrics[v];
    if(dim == 2) fprintf(myfile, "%+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f\n", 
      M(0,0), M(0,1), M(1,1), M(0,2), M(1,2), M(2,2));
    // if(dimelm == 2) fprintf(myfile, "%+-10.12f \t %+-10.12f \t %+-10.12f\n", fmax(A[iVar][i],B[iVar][i]), 0.0, fmax(A[iVar][i],B[iVar][i]));
    // if(dim == 3) fprintf(myfile, "%+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f\n", A[iVar][i], B[iVar][i], D[iVar][i], C[iVar][i], E[iVar][i], F[iVar][i]);
  }

  fprintf(myfile,"End");
  fclose(myfile);
}

void feMetric::writeSizeFieldGmsh(std::string meshName, std::string metricMeshName){
  // Write the size field as NodeData in a copy of the .msh file.
  // For now, the mesh file is simply copied and the NodeData is written at the end.
  feMesh2DP1 *mesh = dynamic_cast<feMesh2DP1*>(_recovery->_mesh);

  if(mesh->isMeshFileBinary()){
    printf("In readGmsh : Error - Only reading ASCII files.\n");
    return;
  }
  if(mesh->getGmshVersion() != 2.2){
    printf("In readGmsh : Error - Anisotropic size field can only be written to a msh 2.2 file for MMG compatibility.\n");
    return;
  }

  std::filebuf fbIn, fbOut;
  if(fbIn.open(meshName,std::ios::in)){
    if(fbOut.open(metricMeshName, std::ios::out)){
      std::istream input(&fbIn);
      std::ostream output(&fbOut);
      std::string buffer;
      // Copy .msh file except for the possible previous NodeData
      while(getline(input,buffer)){
        if(buffer == "$NodeData"){
          while(buffer != "$EndNodeData")
            getline(input,buffer);
          getline(input,buffer);
        }
        output << buffer << std::endl;
      }

      std::vector<int> &vertices = _recovery->getVertices();

      // Append the size field
      output << "$NodeData" << std::endl;
      output << 1 << std::endl; // number of string tag
      output << "\"sizeField:metric\"" << std::endl;
      output << 1 << std::endl; // number of real tag
      output << 0 << std::endl;
      output << 3 << std::endl; // number of integer tag
      output << 0 << std::endl; 
      output << 9 << std::endl; // metric type ( 1:scalar, 3:vector, 9:tensor)
      output << vertices.size() << std::endl;
      for(auto v : vertices){
        // std::cout<<"Writing metric at vertex "<<v;
        printf(" (%f, %f, %f)\n", _recovery->_mesh->getVertex(v)->x(),
            _recovery->_mesh->getVertex(v)->y(),
            _recovery->_mesh->getVertex(v)->z());
        SMetric3 M = _metrics[v];
        output << v+1 << " " << M(0,0) << " " << M(0,1) << " " << M(0,2) << " " << M(1,0) << " " << M(1,1) << " " << M(1,2) << " " << M(2,0)  << " " << M(2,1)  << " " << M(2,2) << std::endl;  
      }
      output << "$EndNodeData" << std::endl;

      fbOut.close();
    }
    fbIn.close();
  }
}