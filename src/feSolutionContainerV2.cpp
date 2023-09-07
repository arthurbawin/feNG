#include "feSolutionContainerV2.h"
#include "feLinearSystem.h"



void feSolutionContainerV2::rotate()
{
	for (int i=_nStorage-1; i>0; i--){
		_sol[i]       = _sol[i-1];
		_solDot[i]    = _solDot[i-1];
		_fResidual[i] = _fResidual[i-1];
	}
}


void feSolutionContainerV2::invRotate()
{
	for (int i=0; i>_nStorage-1; i++){
		_sol[i]       = _sol[i+1];
		_solDot[i]    = _solDot[i+1];
		_fResidual[i] = _fResidual[i+1];
	}
}




void feSolutionContainerV2::computeCoeffBDF(std::vector<double> &k)
{
	if (_orderBDF==0){
		_cBDF[0] = 0.0;
	}else if (_orderBDF==1){
		_cBDF[0] = 1.0/k[0];
		_cBDF[1] =-1.0/k[0];
	}else if (_orderBDF==2){
		_cBDF[0] = 1.0/k[0] + 1.0/(k[0]+k[1]);
    		_cBDF[1] =-1.0/k[0] - 1.0/(k[1]);
    		_cBDF[2] = k[0]/k[1]* 1.0/(k[0]+k[1]);
	}else {
		feErrorMsg(FE_STATUS_ERROR,"Unknown orderBDF");
	}

}






void feSolutionContainerV2::copyCurrentSolution(feSolution *currentSolution)
{
 _sol[0] = currentSolution->getSolutionCopy();
}

void feSolutionContainerV2::copyCurrentSolutionDot(feSolution *currentSolution)
{
 _solDot[0] = currentSolution->getSolutionDotCopy();
}



void feSolutionContainerV2::computeSolDot(feLinearSystem *linearSystem)
{
 std::vector<double> solDotBDF(_nbDDL,0.0);

 for (size_t i=0; i<_cBDF.size(); i++){
 	for (int j=0; j<_nbDDL; j++){
      solDotBDF[j] += _cBDF[i]*_sol[i][j];
    }
  }

 for (int j=0; j<_nbDDL; j++){
 	_solDot[0][j] = solDotBDF[j];
 }

 if (_correctionType=="SOLDOT"){
	for (int j=0; j<_nbDDL; j++){
		_solDot[0][j] += _d[j];
	}

 } else if (_correctionType=="RESIDUAL"){
  	// printf("coucou\n");
  	linearSystem->applyCorrectionToResidual(-1.0, _d);
 } 

}


void feSolutionContainerV2::printSolution(int i)
{
	for (int n=0; n<_nbDDL; n++){
		printf("%10.10f \n", _sol[i][n]);
	}
}


// ================================================
// Compute Coefficient
// ================================================
void feSolutionContainerV2DC2F::computeCorrectionCoeff(std::vector<double> &k)
{
	_cDC2[0] = 1.0/k[0];
	_cDC2[1] =-1.0/k[0];

	// feInfo("%f",_cDC2[0]);
	// feInfo("%f",_cDC2[1]);
}

void feSolutionContainerV2DC3F::computeCorrectionCoeff(std::vector<double> &k)
{
	_cDC3[0] = (2.*k[0] + k[1])/(k[0]*(k[0] + k[1]));
	_cDC3[1] = (-k[0] - k[1])/(k[0]*k[1]);
	_cDC3[2] = k[0]/(k[1]*(k[0] + k[1]));

	_bDC3[0] = (2./k[0]) * (1/(k[0] +k[1]));
	_bDC3[1] = -2./(k[0]*k[1]);
	_bDC3[2] = (2./k[1]) * (1/(k[0] +k[1]));


	// feInfo("%f",k[0]);
	// feInfo("%f",k[1]);
	// printf("\n");
	// feInfo("%f",_cDC3[0]);
	// feInfo("%f",_cDC3[1]);
	// feInfo("%f",_cDC3[2]);
	// printf("\n");
	// feInfo("%f",_bDC3[0]);
	// feInfo("%f",_bDC3[1]);
	// feInfo("%f",_bDC3[2]);
}

void feSolutionContainerV2DC4F::computeCorrectionCoeff(std::vector<double> &k)
{
  	_cDC4[0] = (k[0]*(k[0] + k[1]) + k[0]*(k[0] + k[1] + k[2]) + (k[0] + k[1])*(k[0] + k[1] + k[2]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]));
	_cDC4[1] = (-k[0]*(k[0] + k[1]) - k[0]*(k[1] + k[2]) - k[1]*(k[1] + k[2]))/(k[0]*k[1]*(k[1] + k[2]));
	_cDC4[2] = k[0]*(k[0] + k[1] + k[2])/(k[1]*k[2]*(k[0] + k[1]));
	_cDC4[3] = -k[0]*(k[0] + k[1])/(k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2]));

	_bDC4[0] = 2.*(3.*k[0]+2.*k[1]+k[2])/(k[0]*(k[0]+k[1])*(k[0]+k[1]+k[2]));
 	_bDC4[1] = -2.*(2.*k[0]+2*k[1]+k[2])/((k[0]*k[1])*(k[1]+k[2]));
	_bDC4[2] = 2.*(2.*k[0]+k[1]+k[2])/((k[1]*k[2])*(k[0]+k[1]));
	_bDC4[3] = -2.*(2.*k[0]+k[1])/(k[2]*(k[1]+k[2])*(k[0]+k[1]+k[2]));

	_aDC4[0] = 6./(k[0]*(k[0]+k[1])*(k[0]+k[1]+k[2]));
	_aDC4[1] = -6./(k[0]*k[1]*(k[1]+k[2]));
	_aDC4[2] = 6./(k[1]*k[2]*(k[0]+k[1]));
	_aDC4[3] = -6./(k[2]*(k[1]+k[2])*(k[0]+k[1]+k[2]));

 //  feInfo("%f",_aDC4[0]);
	// feInfo("%f",_aDC4[1]);
	// feInfo("%f",_aDC4[2]);
	// feInfo("%f",_aDC4[3]);
}


void feSolutionContainerV2DC5F::computeCorrectionCoeff(std::vector<double> &k)
{
	_cDC5[0] = (k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]) + k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2] + k[3]) + k[0]*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]) + (k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
  	_cDC5[1] = (-k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]) - k[0]*(k[0] + k[1])*(k[1] + k[2] + k[3]) - k[0]*(k[1] + k[2])*(k[1] + k[2] + k[3]) - k[1]*(k[1] + k[2])*(k[1] + k[2] + k[3]))/(k[0]*k[1]*(k[1] + k[2])*(k[1] + k[2] + k[3]));
	_cDC5[2] = k[0]*(k[2]*(k[2] + k[3]) + (k[0] + k[1])*(k[2] + k[3]) + (k[0] + k[1])*(k[0] + k[1] + k[2]))/(k[1]*k[2]*(k[0] + k[1])*(k[2] + k[3]));
	_cDC5[3] = k[0]*(-pow(k[0],2) - 2*k[0]*k[1] - k[0]*k[2] - k[0]*k[3] - pow(k[1],2) - k[1]*k[2] - k[1]*k[3])/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2*k[1]*k[2] + pow(k[2],2)));
	_cDC5[4] = k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));

	_bDC5[0] = 2.*(k[0]*(k[0] + k[1]) + k[0]*(k[0] + k[1] + k[2]) + (k[0] + k[1])*(k[0] + k[1] + k[2]) + (2*k[0] + k[1])*(k[0] + k[1] + k[2] + k[3]) + (k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
	_bDC5[1] = 2.*(-3.*pow(k[0],2) - 6.*k[0]*k[1] - 4.*k[0]*k[2] - 2.*k[0]*k[3] - 3.*pow(k[1],2) - 4.*k[1]*k[2] - 2.*k[1]*k[3] - pow(k[2],2) - k[2]*k[3])/(k[0]*k[1]*(pow(k[1],2) + 2.*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
	_bDC5[2] = 2.*(3.*pow(k[0],2) + 4.*k[0]*k[1] + 4.*k[0]*k[2] + 2.*k[0]*k[3] + pow(k[1],2) + 2.*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
	_bDC5[3] = 2.*(-3.*pow(k[0],2) - 4.*k[0]*k[1] - 2.*k[0]*k[2] - 2.*k[0]*k[3] - pow(k[1],2) - k[1]*k[2] - k[1]*k[3])/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2.*k[1]*k[2] + pow(k[2],2)));
	_bDC5[4] = 2.*(k[0]*(k[0] + k[1]) + k[0]*(k[0] + k[1] + k[2]) + (k[0] + k[1])*(k[0] + k[1] + k[2]))/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));

	_aDC5[0] = 6.*(4.*k[0] + 3.*k[1] + 2.*k[2] + k[3])/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
	_aDC5[1] = 6.*(-3.*k[0] - 3.*k[1] - 2.*k[2] - k[3])/(k[0]*k[1]*(pow(k[1],2) + 2.*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
	_aDC5[2] = 6.*(3.*k[0] + 2.*k[1] + 2.*k[2] + k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
	_aDC5[3] = 6.*(-3.*k[0] - 2.*k[1] - k[2] - k[3])/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2.*k[1]*k[2] + pow(k[2],2)));
	_aDC5[4] = 6.*(3.*k[0] + 2.*k[1] + k[2])/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));

  	_eDC5[0] = 24./(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
  	_eDC5[1] = -24./(k[0]*k[1]*(pow(k[1],2) + 2.*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
  	_eDC5[2] = 24./(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
  	_eDC5[3] = -24./(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2.*k[1]*k[2] + pow(k[2],2)));
  	_eDC5[4] = 24./(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));

 //  feInfo("%f",_eDC5[0]);
	// feInfo("%f",_eDC5[1]);
	// feInfo("%f",_eDC5[2]);
	// feInfo("%f",_eDC5[3]);
	// feInfo("%f",_eDC5[4]);
}


void feSolutionContainerV2DC3::computeCorrectionCoeff(std::vector<double> &k)
{
	_bDC3[0] = (2./k[0]) * (1/(k[0] +k[1]));
	_bDC3[1] = -2./(k[0]*k[1]);
	_bDC3[2] = (2./k[1]) * (1/(k[0] +k[1]));

	// feInfo("%f",_bDC3[0]);
	// feInfo("%f",_bDC3[1]);
	// feInfo("%f",_bDC3[2]);
}


void feSolutionContainerV2DC4::computeCorrectionCoeff(std::vector<double> &k)
{
	_bDC4[0] = 2.*(3.*k[0]+2.*k[1]+k[2])/(k[0]*(k[0]+k[1])*(k[0]+k[1]+k[2]));
  	_bDC4[1] = -2.*(2.*k[0]+2*k[1]+k[2])/((k[0]*k[1])*(k[1]+k[2]));
  	_bDC4[2] = 2.*(2.*k[0]+k[1]+k[2])/((k[1]*k[2])*(k[0]+k[1]));
  	_bDC4[3] = -2.*(2.*k[0]+k[1])/(k[2]*(k[1]+k[2])*(k[0]+k[1]+k[2]));

	_aDC4[0] = 6./(k[0]*(k[0]+k[1])*(k[0]+k[1]+k[2]));
 	 _aDC4[1] = -6./(k[0]*k[1]*(k[1]+k[2]));
 	 _aDC4[2] = 6./(k[1]*k[2]*(k[0]+k[1]));
 	 _aDC4[3] = -6./(k[2]*(k[1]+k[2])*(k[0]+k[1]+k[2]));

 //  feInfo("%f",_aDC4[0]);
	// feInfo("%f",_aDC4[1]);
	// feInfo("%f",_aDC4[2]);
	// feInfo("%f",_aDC4[3]);
}


void feSolutionContainerV2DC5::computeCorrectionCoeff(std::vector<double> &k)
{
	_bDC5[0] = 2.*(k[0]*(k[0] + k[1]) + k[0]*(k[0] + k[1] + k[2]) + (k[0] + k[1])*(k[0] + k[1] + k[2]) + (2*k[0] + k[1])*(k[0] + k[1] + k[2] + k[3]) + (k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
  	_bDC5[1] = 2.*(-3.*pow(k[0],2) - 6.*k[0]*k[1] - 4.*k[0]*k[2] - 2.*k[0]*k[3] - 3.*pow(k[1],2) - 4.*k[1]*k[2] - 2.*k[1]*k[3] - pow(k[2],2) - k[2]*k[3])/(k[0]*k[1]*(pow(k[1],2) + 2.*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
  	_bDC5[2] = 2.*(3.*pow(k[0],2) + 4.*k[0]*k[1] + 4.*k[0]*k[2] + 2.*k[0]*k[3] + pow(k[1],2) + 2.*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
  	_bDC5[3] = 2.*(-3.*pow(k[0],2) - 4.*k[0]*k[1] - 2.*k[0]*k[2] - 2.*k[0]*k[3] - pow(k[1],2) - k[1]*k[2] - k[1]*k[3])/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2.*k[1]*k[2] + pow(k[2],2)));
  	_bDC5[4] = 2.*(k[0]*(k[0] + k[1]) + k[0]*(k[0] + k[1] + k[2]) + (k[0] + k[1])*(k[0] + k[1] + k[2]))/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));

	_aDC5[0] = 6.*(4.*k[0] + 3.*k[1] + 2.*k[2] + k[3])/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
  	_aDC5[1] = 6.*(-3.*k[0] - 3.*k[1] - 2.*k[2] - k[3])/(k[0]*k[1]*(pow(k[1],2) + 2.*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
  	_aDC5[2] = 6.*(3.*k[0] + 2.*k[1] + 2.*k[2] + k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
  	_aDC5[3] = 6.*(-3.*k[0] - 2.*k[1] - k[2] - k[3])/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2.*k[1]*k[2] + pow(k[2],2)));
  	_aDC5[4] = 6.*(3.*k[0] + 2.*k[1] + k[2])/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));

 	_eDC5[0] = 24./(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
  	_eDC5[1] = -24./(k[0]*k[1]*(pow(k[1],2) + 2.*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
  	_eDC5[2] = 24./(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
  	_eDC5[3] = -24./(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2.*k[1]*k[2] + pow(k[2],2)));
  	_eDC5[4] = 24./(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));

 //  feInfo("%f",_eDC5[0]);
	// feInfo("%f",_eDC5[1]);
	// feInfo("%f",_eDC5[2]);
	// feInfo("%f",_eDC5[3]);
	// feInfo("%f",_eDC5[4]);
}



void feSolutionContainerV2DC6::computeCorrectionCoeff(std::vector<double> &k)
{
	_bDC6[0] = 2.*(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]) + (2.*k[0] + k[1])*(k[0] + k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3] + k[4]) + (k[0]*(k[0] + k[1]) + (2.*k[0] + k[1])*(k[0] + k[1] + k[2]))*(k[0] + k[1] + k[2] + k[3]) + (k[0]*(k[0] + k[1]) + (2.*k[0] + k[1])*(k[0] + k[1] + k[2]))*(k[0] + k[1] + k[2] + k[3] + k[4]) + (k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3] + k[4]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3] + k[4]));
	_bDC6[1] = 2.*(-k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]) - (2.*k[0] + k[1])*(k[1] + k[2] + k[3])*(k[1] + k[2] + k[3] + k[4]) - (k[1] + k[2])*(k[1] + k[2] + k[3])*(k[1] + k[2] + k[3] + k[4]) - (k[0]*(k[0] + k[1]) + (2.*k[0] + k[1])*(k[0] + k[1] + k[2]))*(k[0] + k[1] + k[2] + k[3]) - (k[0]*(k[0] + k[1]) + (2.*k[0] + k[1])*(k[0] + k[1] + k[2]))*(k[1] + k[2] + k[3] + k[4]))/(k[0]*k[1]*(k[1] + k[2])*(k[1] + k[2] + k[3])*(k[1] + k[2] + k[3] + k[4]));
	_bDC6[2] = 2.*(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]) + k[2]*(k[2] + k[3])*(k[2] + k[3] + k[4]) + (2.*k[0] + k[1])*(k[2] + k[3])*(k[2] + k[3] + k[4]) + (k[0]*(k[0] + k[1]) + (2.*k[0] + k[1])*(k[0] + k[1] + k[2]))*(k[2] + k[3] + k[4]) + (k[0]*(k[0] + k[1]) + (2.*k[0] + k[1])*(k[0] + k[1] + k[2]))*(k[0] + k[1] + k[2] + k[3]))/(k[1]*k[2]*(k[0] + k[1])*(k[2] + k[3])*(k[2] + k[3] + k[4])); 
	_bDC6[3] = 2.*(-4.*pow(k[0],3) - 9.*pow(k[0],2)*k[1] - 6.*pow(k[0],2)*k[2] - 6.*pow(k[0],2)*k[3] - 3.*pow(k[0],2)*k[4] - 6.*k[0]*pow(k[1],2) - 8.*k[0]*k[1]*k[2] - 8.*k[0]*k[1]*k[3] - 4.*k[0]*k[1]*k[4] - 2.*k[0]*pow(k[2],2) - 4.*k[0]*k[2]*k[3] - 2.*k[0]*k[2]*k[4] - 2.*k[0]*pow(k[3],2) - 2.*k[0]*k[3]*k[4] - pow(k[1],3) - 2.*pow(k[1],2)*k[2] - 2.*pow(k[1],2)*k[3] - pow(k[1],2)*k[4] - k[1]*pow(k[2],2) - 2.*k[1]*k[2]*k[3] - k[1]*k[2]*k[4] - k[1]*pow(k[3],2) - k[1]*k[3]*k[4])/(k[2]*k[3]*(k[0]*k[1]*k[3] + k[0]*k[1]*k[4] + k[0]*k[2]*k[3] + k[0]*k[2]*k[4] + pow(k[1],2)*k[3] + pow(k[1],2)*k[4] + 2.*k[1]*k[2]*k[3] + 2.*k[1]*k[2]*k[4] + pow(k[2],2)*k[3] + pow(k[2],2)*k[4]));
	_bDC6[4] = 2.*(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]) + k[4]*(k[0]*(k[0] + k[1]) + (2.*k[0] + k[1])*(k[0] + k[1] + k[2])) + (k[0]*(k[0] + k[1]) + (2.*k[0] + k[1])*(k[0] + k[1] + k[2]))*(k[0] + k[1] + k[2] + k[3]))/(k[3]*k[4]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
	_bDC6[5] = 2.*(-k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]) - (k[0]*(k[0] + k[1]) + (2.*k[0] + k[1])*(k[0] + k[1] + k[2]))*(k[0] + k[1] + k[2] + k[3]))/(k[4]*(k[3] + k[4])*(k[2] + k[3] + k[4])*(k[1] + k[2] + k[3] + k[4])*(k[0] + k[1] + k[2] + k[3] + k[4]));
	
	_aDC6[0] = 6.*(3.*k[0]*(2*k[0] + 3.*k[1] + 2.*k[2] + k[3]) + k[1]*(3.*k[1] + 4.*k[2] + 2.*k[3]) + k[2]*(k[2] + k[3]) + (3.*k[0] + 2.*k[1] + k[2])*(k[0] + k[1] + k[2] + k[3] + k[4]) + (k[0] + k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3] + k[4]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3] + k[4]));
  	_aDC6[1] = 6.*(-3.*k[0]*(2*k[0] + 3.*k[1] + 2.*k[2] + k[3]) - k[1]*(3.*k[1] + 4.*k[2] + 2.*k[3]) - k[2]*(k[2] + k[3]) - (3.*k[0] + 2.*k[1] + k[2])*(k[1] + k[2] + k[3] + k[4]) - (k[1] + k[2] + k[3])*(k[1] + k[2] + k[3] + k[4]))/(k[0]*k[1]*(k[1] + k[2])*(k[1] + k[2] + k[3])*(k[1] + k[2] + k[3] + k[4]));
  	_aDC6[2] = 6.*(3.*k[0]*(2*k[0] + 3.*k[1] + 2.*k[2] + k[3]) + k[1]*(3.*k[1] + 4.*k[2] + 2.*k[3]) + k[2]*(k[2] + k[3]) + (k[2] + k[3])*(k[2] + k[3] + k[4]) + (3.*k[0] + 2.*k[1] + k[2])*(k[2] + k[3] + k[4]))/(k[1]*k[2]*(k[0] + k[1])*(k[2] + k[3])*(k[2] + k[3] + k[4]));
  	_aDC6[3] = 6.*(-6.*pow(k[0],2) - 9.*k[0]*k[1] - 6.*k[0]*k[2] - 6.*k[0]*k[3] - 3.*k[0]*k[4] - 3.*pow(k[1],2) - 4.*k[1]*k[2] - 4.*k[1]*k[3] - 2.*k[1]*k[4] - pow(k[2],2) - 2.*k[2]*k[3] - k[2]*k[4] - pow(k[3],2) - k[3]*k[4])/(k[2]*k[3]*(k[0]*k[1]*k[3] + k[0]*k[1]*k[4] + k[0]*k[2]*k[3] + k[0]*k[2]*k[4] + pow(k[1],2)*k[3] + pow(k[1],2)*k[4] + 2.*k[1]*k[2]*k[3] + 2.*k[1]*k[2]*k[4] + pow(k[2],2)*k[3] + pow(k[2],2)*k[4]));
  	_aDC6[4] = 6.*(6.*pow(k[0],2) + 9.*k[0]*k[1] + 6.*k[0]*k[2] + 3.*k[0]*k[3] + 3.*k[0]*k[4] + 3.*pow(k[1],2) + 4.*k[1]*k[2] + 2.*k[1]*k[3] + 2.*k[1]*k[4] + pow(k[2],2) + k[2]*k[3] + k[2]*k[4])/(k[3]*k[4]*(k[0]*k[1]*k[2] + k[0]*k[1]*k[3] + k[0]*pow(k[2],2) + 2.*k[0]*k[2]*k[3] + k[0]*pow(k[3],2) + pow(k[1],2)*k[2] + pow(k[1],2)*k[3] + 2.*k[1]*pow(k[2],2) + 4.*k[1]*k[2]*k[3] + 2.*k[1]*pow(k[3],2) + pow(k[2],3) + 3.*pow(k[2],2)*k[3] + 3.*k[2]*pow(k[3],2) + pow(k[3],3)));
  	_aDC6[5] = 6.*(-3.*k[0]*(2*k[0] + 3.*k[1] + 2.*k[2] + k[3]) - k[1]*(3.*k[1] + 4.*k[2] + 2.*k[3]) - k[2]*(k[2] + k[3]))/(k[4]*(k[3] + k[4])*(k[2] + k[3] + k[4])*(k[1] + k[2] + k[3] + k[4])*(k[0] + k[1] + k[2] + k[3] + k[4]));
            
	_eDC6[0] = 24.*(5.*k[0] + 4.*k[1] + 3.*k[2] + 2.*k[3] + k[4])/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3] + k[4]));
  	_eDC6[1] = 24.*(-4.*k[0] - 4.*k[1] - 3.*k[2] - 2.*k[3] - k[4])/(k[0]*k[1]*(pow(k[1],3) + 3.*pow(k[1],2)*k[2] + 2.*pow(k[1],2)*k[3] + pow(k[1],2)*k[4] + 3*k[1]*pow(k[2],2) + 4.*k[1]*k[2]*k[3] + 2.*k[1]*k[2]*k[4] + k[1]*pow(k[3],2) + k[1]*k[3]*k[4] + pow(k[2],3) + 2.*pow(k[2],2)*k[3] + pow(k[2],2)*k[4] + k[2]*pow(k[3],2) + k[2]*k[3]*k[4]));
  	_eDC6[2] = 24.*(4.*k[0] + 3.*k[1] + 3.*k[2] + 2.*k[3] + k[4])/(k[1]*k[2]*(k[0]*pow(k[2],2) + 2.*k[0]*k[2]*k[3] + k[0]*k[2]*k[4] + k[0]*pow(k[3],2) + k[0]*k[3]*k[4] + k[1]*pow(k[2],2) + 2.*k[1]*k[2]*k[3] + k[1]*k[2]*k[4] + k[1]*pow(k[3],2) + k[1]*k[3]*k[4]));
  	_eDC6[3] = 24.*(-4.*k[0] - 3.*k[1] - 2.*k[2] - 2.*k[3] - k[4])/(k[2]*k[3]*(k[0]*k[1]*k[3] + k[0]*k[1]*k[4] + k[0]*k[2]*k[3] + k[0]*k[2]*k[4] + pow(k[1],2)*k[3] + pow(k[1],2)*k[4] + 2.*k[1]*k[2]*k[3] + 2.*k[1]*k[2]*k[4] + pow(k[2],2)*k[3] + pow(k[2],2)*k[4]));
  	_eDC6[4] = 24.*(4.*k[0] + 3.*k[1] + 2.*k[2] + k[3] + k[4])/(k[3]*k[4]*(k[0]*k[1]*k[2] + k[0]*k[1]*k[3] + k[0]*pow(k[2],2) + 2*k[0]*k[2]*k[3] + k[0]*pow(k[3],2) + pow(k[1],2)*k[2] + pow(k[1],2)*k[3] + 2.*k[1]*pow(k[2],2) + 4.*k[1]*k[2]*k[3] + 2.*k[1]*pow(k[3],2) + pow(k[2],3) + 3.*pow(k[2],2)*k[3] + 3.*k[2]*pow(k[3],2) + pow(k[3],3)));
  	_eDC6[5] = 24.*(-4.*k[0] - 3.*k[1] - 2.*k[2] - k[3])/(k[4]*(k[3] + k[4])*(k[2] + k[3] + k[4])*(k[1] + k[2] + k[3] + k[4])*(k[0] + k[1] + k[2] + k[3] + k[4]));

  	_fDC6[0] = 120./(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3] + k[4]));
  	_fDC6[1] = -120./(k[0]*k[1]*(k[1] + k[2])*(k[1] + k[2] + k[3])*(k[1] + k[2] + k[3] + k[4]));
  	_fDC6[2] = 120./(k[1]*k[2]*(k[0]*pow(k[2],2) + 2*k[0]*k[2]*k[3] + k[0]*k[2]*k[4] + k[0]*pow(k[3],2) + k[0]*k[3]*k[4] + k[1]*pow(k[2],2) + 2.*k[1]*k[2]*k[3] + k[1]*k[2]*k[4] + k[1]*pow(k[3],2) + k[1]*k[3]*k[4]));
  	_fDC6[3] = -120./(k[2]*k[3]*(k[0]*k[1]*k[3] + k[0]*k[1]*k[4] + k[0]*k[2]*k[3] + k[0]*k[2]*k[4] + pow(k[1],2)*k[3] + pow(k[1],2)*k[4] + 2.*k[1]*k[2]*k[3] + 2.*k[1]*k[2]*k[4] + pow(k[2],2)*k[3] + pow(k[2],2)*k[4]));
  	_fDC6[4] = 120./(k[3]*k[4]*(k[0]*k[1]*k[2] + k[0]*k[1]*k[3] + k[0]*pow(k[2],2) + 2.*k[0]*k[2]*k[3] + k[0]*pow(k[3],2) + pow(k[1],2)*k[2] + pow(k[1],2)*k[3] + 2.*k[1]*pow(k[2],2) + 4.*k[1]*k[2]*k[3] + 2.*k[1]*pow(k[3],2) + pow(k[2],3) + 3*pow(k[2],2)*k[3] + 3.*k[2]*pow(k[3],2) + pow(k[3],3)));
  	_fDC6[5] = -120./(k[4]*(k[3] + k[4])*(k[2] + k[3] + k[4])*(k[1] + k[2] + k[3] + k[4])*(k[0] + k[1] + k[2] + k[3] + k[4])); 
}


// =========================================================
// Start Correction Coefficient
// =========================================================

void feSolutionContainerV2DC3F::computeStartCorrectionCoeff(std::vector<double> &k)
{
	_cStartDC3_t1.resize(3);
    	_bStartDC3_t1.resize(3);

	_cStartDC3_t1[0] = k[1]/(k[0]*(k[0] + k[1]));
	_cStartDC3_t1[1] = (k[0] - k[1])/(k[0]*k[1]) ;
	_cStartDC3_t1[2] = -k[0]/(k[1]*(k[0] + k[1]));

	_bStartDC3_t1[0] = 2/(k[0]*(k[0] + k[1]));
	_bStartDC3_t1[1] = -2/(k[0]*k[1]);
	_bStartDC3_t1[2] = 2/(k[1]*(k[0] + k[1]));

}


void feSolutionContainerV2DC4F::computeStartCorrectionCoeff(std::vector<double> &k)
{
	_cStartDC4_t1.resize(4);
    	_bStartDC4_t1.resize(4);
    	_aStartDC4_t1.resize(4);

    	_cStartDC4_t2.resize(4);
    	_bStartDC4_t2.resize(4);
    	_aStartDC4_t2.resize(4);

    	_cStartDC4_t1[0] = -k[1]*k[2]/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]));
    	_cStartDC4_t1[1] = k[2]*(k[0] + k[1])/(k[0]*k[1]*(k[1] + k[2]));
    	_cStartDC4_t1[2] = (-k[1]*k[2] + k[1]*(k[0] + k[1]) - k[2]*(k[0] + k[1]))/(k[1]*k[2]*(k[0] + k[1]));
    	_cStartDC4_t1[3] = -k[1]*(k[0] + k[1])/(k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2]));

    	_cStartDC4_t2[0] = k[1]*(k[1] + k[2])/(k[0]*(pow(k[0],2) + 2*k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + k[1]*k[2]));                
     _cStartDC4_t2[1] = (k[0]*k[1] + k[0]*(k[1] + k[2]) - k[1]*(k[1] + k[2]))/(k[0]*k[1]*(k[1] + k[2]));
	_cStartDC4_t2[2] = -k[0]*(k[1] + k[2])/(k[1]*k[2]*(k[0] + k[1]));
	_cStartDC4_t2[3] = k[0]*k[1]/(k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
	  
	_bStartDC4_t1[0] = 2*(-k[1] + k[2])/(k[0]*(pow(k[0],2) + 2*k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + k[1]*k[2]));
	_bStartDC4_t1[1] = 2*(k[0] + k[1] - k[2])/(k[0]*k[1]*(k[1] + k[2]));
	_bStartDC4_t1[2] = 2*(-k[0] - 2*k[1] + k[2])/(k[1]*k[2]*(k[0] + k[1]));
	_bStartDC4_t1[3] = 2*(k[0] + 2*k[1])/(k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
	  
	_bStartDC4_t2[0] = 2*(2*k[1] + k[2])/(k[0]*(pow(k[0],2) + 2*k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + k[1]*k[2]));
	_bStartDC4_t2[1] = 2*(k[0] - 2*k[1] - k[2])/(k[0]*k[1]*(k[1] + k[2]));
	_bStartDC4_t2[2] = 2*(-k[0] + k[1] + k[2])/(k[1]*k[2]*(k[0] + k[1]));
	_bStartDC4_t2[3] = 2*(k[0] - k[1])/(k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
	  
	_aStartDC4_t1[0] = 6/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]));
	_aStartDC4_t1[1] = -6/(k[0]*k[1]*(k[1] + k[2]));
	_aStartDC4_t1[2] = 6/(k[1]*k[2]*(k[0] + k[1]));
	_aStartDC4_t1[3] = -6/(k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
	  
	_aStartDC4_t2[0] = 6/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2]));
	_aStartDC4_t2[1] = -6/(k[0]*k[1]*(k[1] + k[2]));
	_aStartDC4_t2[2] = 6/(k[1]*k[2]*(k[0] + k[1]));
	_aStartDC4_t2[3] = -6/(k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
}



void feSolutionContainerV2DC5F::computeStartCorrectionCoeff(std::vector<double> &k)
{
	_cStartDC5_t1.resize(5);
    	_bStartDC5_t1.resize(5);
    	_aStartDC5_t1.resize(5);
    	_eStartDC5_t1.resize(5);

    	_cStartDC5_t2.resize(5);
    	_bStartDC5_t2.resize(5);
    	_aStartDC5_t2.resize(5);
    	_eStartDC5_t2.resize(5);

    	_cStartDC5_t3.resize(5);
    	_bStartDC5_t3.resize(5);
    	_aStartDC5_t3.resize(5);
    	_eStartDC5_t3.resize(5);

    	_cStartDC5_t1[0] = (-k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2]) + k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
     _cStartDC5_t1[1] = -k[2]*k[3]*(k[0] + k[1] + k[2])/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3])) ;
     _cStartDC5_t1[2] = k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2*k[1]*k[2] + pow(k[2],2))/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _cStartDC5_t1[3] = (-k[2]*k[3]*(k[1] + k[2]) - k[2]*k[3]*(k[0] + k[1] + k[2]) + k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2]) - k[3]*(k[1] + k[2])*(k[0] + k[1] + k[2]))/(k[2]*k[3]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
     _cStartDC5_t1[4] = -k[2]*(k[1] + k[2])*(k[0] + k[1] + k[2])/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
     _cStartDC5_t2[0] = k[1]*(k[2]*(k[0] + k[1]) + (k[0] + k[1])*(k[0] + k[1] + k[2] + k[3]) - (k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));                                                          
     _cStartDC5_t2[1] = k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3])/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _cStartDC5_t2[2] = (2*k[0]*k[1]*k[2] + k[0]*k[1]*k[3] - k[0]*pow(k[2],2) - k[0]*k[2]*k[3] + 2*pow(k[1],2)*k[2] + pow(k[1],2)*k[3] - 2*k[1]*pow(k[2],2) - 2*k[1]*k[2]*k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _cStartDC5_t2[3] = -k[1]*(k[0] + k[1])*(k[2] + k[3])/(k[2]*k[3]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
     _cStartDC5_t2[4] = k[1]*k[2]*(k[0] + k[1])/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
     _cStartDC5_t3[0] = (-k[0]*k[1]*(k[1] + k[2]) - k[0]*k[1]*(k[0] + k[1] + k[2] + k[3]) - k[0]*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]) + (k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));                                                                                                               
     _cStartDC5_t3[1] = (3*k[0]*pow(k[1],2) + 4*k[0]*k[1]*k[2] + 2*k[0]*k[1]*k[3] + k[0]*pow(k[2],2) + k[0]*k[2]*k[3] - pow(k[1],3) - 2*pow(k[1],2)*k[2] - pow(k[1],2)*k[3] - k[1]*pow(k[2],2) - k[1]*k[2]*k[3])/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _cStartDC5_t3[2] = -k[0]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _cStartDC5_t3[3] = k[0]*k[1]*(k[1] + k[2] + k[3])/(k[2]*k[3]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
     _cStartDC5_t3[4] = -k[0]*k[1]*(k[1] + k[2])/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
            
     _bStartDC5_t1[0] = 2*(k[2]*(k[1] + k[2]) + (k[1] + 2*k[2])*(k[0] + k[1] + k[2]) - (k[1] + 2*k[2])*(k[0] + k[1] + k[2] + k[3]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
     _bStartDC5_t1[1] = 2*(-k[0]*k[2] + k[0]*k[3] - k[1]*k[2] + k[1]*k[3] - pow(k[2],2) + 2*k[2]*k[3])/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _bStartDC5_t1[2] = 2*(k[0]*k[1] + k[0]*k[2] - k[0]*k[3] + pow(k[1],2) + 2*k[1]*k[2] - 2*k[1]*k[3] + pow(k[2],2) - 2*k[2]*k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _bStartDC5_t1[3] = 2*(-k[0]*k[1] - 2*k[0]*k[2] + k[0]*k[3] - pow(k[1],2) - 4*k[1]*k[2] + 2*k[1]*k[3] - 3*pow(k[2],2) + 3*k[2]*k[3])/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2*k[1]*k[2] + pow(k[2],2)));
     _bStartDC5_t1[4] = 2*(k[2]*(k[1] + k[2]) + k[2]*(k[0] + k[1] + k[2]) + (k[1] + k[2])*(k[0] + k[1] + k[2]))/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
     _bStartDC5_t2[0] = 2*(-k[1]*k[2] - k[1]*(k[0] + k[1] + k[2] + k[3]) + (k[0] + k[1])*(k[1] - k[2]) - (k[0] + k[1])*(k[0] + k[1] + k[2] + k[3]) + (k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
     _bStartDC5_t2[1] = 2*(2*k[0]*k[2] + k[0]*k[3] + 2*k[1]*k[2] + k[1]*k[3] - pow(k[2],2) - k[2]*k[3])/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _bStartDC5_t2[2] = 2*(k[0]*k[1] - 2*k[0]*k[2] - k[0]*k[3] + pow(k[1],2) - 4*k[1]*k[2] - 2*k[1]*k[3] + pow(k[2],2) + k[2]*k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _bStartDC5_t2[3] = 2*(-k[1]*(k[0] + k[1]) + k[2]*(k[0] + 2*k[1]) + k[3]*(k[0] + 2*k[1]))/(k[2]*k[3]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
     _bStartDC5_t2[4] = 2*(-k[1]*k[2] + k[1]*(k[0] + k[1]) - k[2]*(k[0] + k[1]))/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
     _bStartDC5_t3[0] = 2*(-k[0]*(2*k[1] + k[2]) - k[0]*(k[0] + k[1] + k[2] + k[3]) + k[1]*(k[1] + k[2]) + k[1]*(k[0] + k[1] + k[2] + k[3]) + (k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]))/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
     _bStartDC5_t3[1] = 2*(3*k[0]*k[1] + 2*k[0]*k[2] + k[0]*k[3] - 3*pow(k[1],2) - 4*k[1]*k[2] - 2*k[1]*k[3] - pow(k[2],2) - k[2]*k[3])/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _bStartDC5_t3[2] = 2*(-2*k[0]*k[1] - 2*k[0]*k[2] - k[0]*k[3] + pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _bStartDC5_t3[3] = 2*(k[0]*k[1] + k[3]*(k[0] - k[1]) + (k[0] - k[1])*(k[1] + k[2]))/(k[2]*k[3]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
     _bStartDC5_t3[4] = 2*(-k[0]*k[1] - k[0]*(k[1] + k[2]) + k[1]*(k[1] + k[2]))/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
            
     _aStartDC5_t1[0] = 6*(-k[1] - 2*k[2] + k[3])/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
     _aStartDC5_t1[1] = 6*(k[0] + k[1] + 2*k[2] - k[3])/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _aStartDC5_t1[2] = 6*(-k[0] - 2*k[1] - 2*k[2] + k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _aStartDC5_t1[3] = 6*(k[0] + 2*k[1] + 3*k[2] - k[3])/(k[2]*k[3]*(k[1] + k[2])*(k[0] + k[1] + k[2]));
     _aStartDC5_t1[4] = -(6*k[0] + 12*k[1] + 18*k[2])/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
     _aStartDC5_t2[0] = 6*(-k[1] + 2*k[2] + k[3])/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
     _aStartDC5_t2[1] = 6*(k[0] + k[1] - 2*k[2] - k[3])/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _aStartDC5_t2[2] = 6*(-k[0] - 2*k[1] + 2*k[2] + k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _aStartDC5_t2[3] = 6*(k[0] + 2*k[1] - k[2] - k[3])/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2*k[1]*k[2] + pow(k[2],2)));
     _aStartDC5_t2[4] = 6*(-k[0] - 2*k[1] + k[2])/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
     _aStartDC5_t3[0] = 6*(3*k[1] + 2*k[2] + k[3])/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
     _aStartDC5_t3[1] = 6*(k[0] - 3*k[1] - 2*k[2] - k[3])/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _aStartDC5_t3[2] = 6*(-k[0] + 2*k[1] + 2*k[2] + k[3])/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _aStartDC5_t3[3] = 6*(k[0] - 2*k[1] - k[2] - k[3])/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2*k[1]*k[2] + pow(k[2],2)));
     _aStartDC5_t3[4] = 6*(-k[0] + 2*k[1] + k[2])/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
            
     _eStartDC5_t1[0] = 24/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
     _eStartDC5_t1[1] = -24/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _eStartDC5_t1[2] = 24/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _eStartDC5_t1[3] = -24/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2*k[1]*k[2] + pow(k[2],2)));
     _eStartDC5_t1[4] = 24/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
     _eStartDC5_t1[0] = 24/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
     _eStartDC5_t1[1] = -24/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _eStartDC5_t1[2] = 24/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _eStartDC5_t1[3] = -24/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2*k[1]*k[2] + pow(k[2],2)));
     _eStartDC5_t1[4] = 24/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
            
     _eStartDC5_t1[0] = 24/(k[0]*(k[0] + k[1])*(k[0] + k[1] + k[2])*(k[0] + k[1] + k[2] + k[3]));
     _eStartDC5_t1[1] = -24/(k[0]*k[1]*(pow(k[1],2) + 2*k[1]*k[2] + k[1]*k[3] + pow(k[2],2) + k[2]*k[3]));
     _eStartDC5_t1[2] = 24/(k[1]*k[2]*(k[0]*k[2] + k[0]*k[3] + k[1]*k[2] + k[1]*k[3]));
     _eStartDC5_t1[3] = -24/(k[2]*k[3]*(k[0]*k[1] + k[0]*k[2] + pow(k[1],2) + 2*k[1]*k[2] + pow(k[2],2)));
     _eStartDC5_t1[4] = 24/(k[3]*(k[2] + k[3])*(k[1] + k[2] + k[3])*(k[0] + k[1] + k[2] + k[3]));
}


// ================================================
// Compute Correction
// ================================================

void feSolutionContainerV2Stat::computeCorrection()
{
	for(int i=0; i<_nbInc; i++)
		_d[i]=0.0;
}

void feSolutionContainerV2BDF::computeCorrection()
{
	for(int i=0; i<_nbInc; i++)
		_d[i]=0.0;
}	

void feSolutionContainerV2DC2F::computeCorrection(feSolutionContainerV2BDF *solutionContainerBDF1, std::vector<double> &k) 
{
	this->computeCorrectionCoeff(k);

	int nRange = 2;
	std::vector<double> C(_nbDDL,0.0);

	if (_correctionType=="SOLDOT"){
		for (int i=0; i< nRange; i++){
			std::vector<double> &solDotBDF = solutionContainerBDF1->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				C[j] += _cDC2[i]*solDotBDF[j];
			}
		}
	
	} else if (_correctionType=="RESIDUAL"){
		for (int i=0; i< nRange; i++){
			std::vector<double> &solResidualBDF = solutionContainerBDF1->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				C[j] += _cDC2[i]*solResidualBDF[j];
			}
		}
	
	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}

	for (int j=0; j<_nbDDL; j++){
		_d[j] = C[j]*k[0]/2.0;
	}

}


void feSolutionContainerV2DC3F::computeCorrection(feSolutionContainerV2DC2F *solutionContainerDC2F, std::vector<double> &k)
{
	this->computeCorrectionCoeff(k);

	int nRange=3;
	std::vector<double> C(_nbDDL, 0.0);
	std::vector<double> B(_nbDDL, 0.0);

	if (_correctionType=="SOLDOT"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solDotDC2F = solutionContainerDC2F->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				C[j] += _cDC3[i]*solDotDC2F[j];
				B[j] -= _bDC3[i]*solDotDC2F[j];
			}
		}
	
	} else if (_correctionType=="RESIDUAL"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solResidualDC2F = solutionContainerDC2F->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				C[j] += _cDC3[i]*solResidualDC2F[j];
				B[j] -= _bDC3[i]*solResidualDC2F[j];
			}
		}

	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}

	for (int j=0; j<_nbDDL; j++){
		_d[j] = C[j]*k[0]/2.0 + B[j]*pow(k[0],2)/(2.*3.);
	}
	
}


void feSolutionContainerV2DC4F::computeCorrection(feSolutionContainerV2DC3F *solutionContainerDC3F, std::vector<double> &k)
{
	this->computeCorrectionCoeff(k);

	int nRange=4;
	std::vector<double> C(_nbDDL, 0.0);
	std::vector<double> B(_nbDDL, 0.0);
	std::vector<double> A(_nbDDL, 0.0);
	
	if (_correctionType=="SOLDOT"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solDotDC3F = solutionContainerDC3F->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				C[j] += _cDC4[i]*solDotDC3F[j];
				B[j] -= _bDC4[i]*solDotDC3F[j];
				A[j] += _aDC4[i]*solDotDC3F[j];
			}
		}
	
	} else if (_correctionType=="RESIDUAL"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &residualDC3F = solutionContainerDC3F->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				C[j] += _cDC4[i]*residualDC3F[j];
				B[j] -= _bDC4[i]*residualDC3F[j];
				A[j] += _aDC4[i]*residualDC3F[j];
			}
		}
	
	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}

	for (int j=0; j<_nbDDL; j++){
		_d[j] = C[j]*k[0]/2.0 + B[j]*pow(k[0],2)/(2.*3.) + A[j]*pow(k[0],3)/(2.*3.*4.);
	}
	
}


void feSolutionContainerV2DC5F::computeCorrection(feSolutionContainerV2DC4F *solutionContainerDC4F, std::vector<double> &k)
{
	this->computeCorrectionCoeff(k);

	int nRange=5;
	std::vector<double> C(_nbDDL, 0.0);
	std::vector<double> B(_nbDDL, 0.0);
	std::vector<double> A(_nbDDL, 0.0);
	std::vector<double> E(_nbDDL, 0.0);

	if (_correctionType=="SOLDOT"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solDotDC4F = solutionContainerDC4F->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				C[j] += _cDC5[i]*solDotDC4F[j];
				B[j] -= _bDC5[i]*solDotDC4F[j];
				A[j] += _aDC5[i]*solDotDC4F[j];
				E[j] -= _eDC5[i]*solDotDC4F[j];
			}
		}
	
	} else if (_correctionType=="RESIDUAL"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &residualDC4F = solutionContainerDC4F->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				C[j] += _cDC5[i]*residualDC4F[j];
				B[j] -= _bDC5[i]*residualDC4F[j];
				A[j] += _aDC5[i]*residualDC4F[j];
				E[j] -= _eDC5[i]*residualDC4F[j];
			}
		}

	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}

	for (int j=0; j<_nbDDL; j++){
		_d[j] = C[j]*k[0]/2.0 + B[j]*pow(k[0],2)/(2.*3.) + A[j]*pow(k[0],3)/(2.*3.*4.) + E[j]*pow(k[0],4)/(2.*3.*4.*5.);
	}
	
}


void feSolutionContainerV2DC3::computeCorrection(feSolutionContainerV2BDF *solutionContainerBDF2, std::vector<double> &k)
{
	int nRange=3;
	std::vector<double> B(_nbDDL, 0.0);
	if (_correctionType=="SOLDOT"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solDotBDF2 = solutionContainerBDF2->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				B[j] += _bDC3[i]*solDotBDF2[j];
			}
		}
	
	} else if (_correctionType=="RESIDUAL") {
		for (int i=0; i<nRange; i++){
			std::vector<double> &residualBDF2 = solutionContainerBDF2->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				B[j] += _bDC3[i]*residualBDF2[j];
			}
		}
	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}
	
	for (int j=0; j<_nbDDL; j++){
		_d[j] = B[j]*(_cBDF[1]*pow(k[0],3) + _cBDF[2]*pow(k[0]+k[1],3))/(2.*3.);
	}
}



void feSolutionContainerV2DC4::computeCorrection(feSolutionContainerV2DC3 *solutionContainerDC3, std::vector<double> &k)
{
	int nRange=4;
	std::vector<double> B(_nbDDL, 0.0);
	std::vector<double> A(_nbDDL, 0.0);

	if (_correctionType=="SOLDOT") {
		for (int i=0; i<nRange; i++){
			std::vector<double> &solDotDC3 = solutionContainerDC3->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				B[j] += _bDC4[i]*solDotDC3[j];
				A[j] -= _aDC4[i]*solDotDC3[j];
			}
		}

	} else if (_correctionType=="RESIDUAL") {
		for (int i=0; i<nRange; i++){
			std::vector<double> &residualDC3 = solutionContainerDC3->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				B[j] += _bDC4[i]*residualDC3[j];
				A[j] -= _aDC4[i]*residualDC3[j];
			}
		}

	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}
	
	for (int j=0; j<_nbDDL; j++){
		_d[j] = B[j]*(_cBDF[1]*pow(k[0],3) + _cBDF[2]*pow(k[0]+k[1],3))/(2.*3.) + A[j]*(_cBDF[1]*pow(k[0],4) + _cBDF[2]*pow(k[0]+k[1],4))/(2.*3.*4.);
	}
	
}

void feSolutionContainerV2DC5::computeCorrection(feSolutionContainerV2DC4 *solutionContainerDC4, std::vector<double> &k)
{
	int nRange=5;
	std::vector<double> B(_nbDDL, 0.0);
	std::vector<double> A(_nbDDL, 0.0);
	std::vector<double> E(_nbDDL, 0.0);

	if (_correctionType =="SOLDOT") {
		for (int i=0; i<nRange; i++){
			std::vector<double> &solDotDC4 = solutionContainerDC4->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				B[j] += _bDC5[i]*solDotDC4[j];
				A[j] -= _aDC5[i]*solDotDC4[j];
				E[j] += _eDC5[i]*solDotDC4[j];
			}
		}
	} else if (_correctionType =="RESIDUAL"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &residualDC4 = solutionContainerDC4->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				B[j] += _bDC5[i]*residualDC4[j];
				A[j] -= _aDC5[i]*residualDC4[j];
				E[j] += _eDC5[i]*residualDC4[j];
			}
		}

	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}
	
	for (int j=0; j<_nbDDL; j++){
		_d[j] = B[j]*(_cBDF[1]*pow(k[0],3) + _cBDF[2]*pow(k[0]+k[1],3))/(2.*3.) + A[j]*(_cBDF[1]*pow(k[0],4) + _cBDF[2]*pow(k[0]+k[1],4))/(2.*3.*4.) + E[j]*(_cBDF[1]*pow(k[0],5) + _cBDF[2]*pow(k[0]+k[1],5))/(2.*3.*4.*5.);
	}
	
}


void feSolutionContainerV2DC6::computeCorrection(feSolutionContainerV2DC5 *solutionContainerDC5, std::vector<double> &k)
{
	int nRange=6;
	std::vector<double> B(_nbDDL, 0.0);
	std::vector<double> A(_nbDDL, 0.0);
	std::vector<double> E(_nbDDL, 0.0);
	std::vector<double> F(_nbDDL, 0.0);

	if (_correctionType=="SOLDOT"){ 
		for (int i=0; i<nRange; i++){
			std::vector<double> &solDotDC5 = solutionContainerDC5->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				B[j] += _bDC6[i]*solDotDC5[j];
				A[j] -= _aDC6[i]*solDotDC5[j];
				E[j] += _eDC6[i]*solDotDC5[j];
				F[j] -= _fDC6[i]*solDotDC5[j];
			}
		}
	
	} else if (_correctionType=="RESIDUAL"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &residualDC5 = solutionContainerDC5->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				B[j] += _bDC6[i]*residualDC5[j];
				A[j] -= _aDC6[i]*residualDC5[j];
				E[j] += _eDC6[i]*residualDC5[j];
				F[j] -= _fDC6[i]*residualDC5[j];
			}
		}

	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}


	for (int j=0; j<_nbDDL; j++){
		_d[j] = B[j]*(_cBDF[1]*pow(k[0],3) + _cBDF[2]*pow(k[0]+k[1],3))/(2.*3.) + A[j]*(_cBDF[1]*pow(k[0],4) + _cBDF[2]*pow(k[0]+k[1],4))/(2.*3.*4.) + E[j]*(_cBDF[1]*pow(k[0],5) + _cBDF[2]*pow(k[0]+k[1],5))/(2.*3.*4.*5.) + F[j]*(_cBDF[1]*pow(k[0],6) + _cBDF[2]*pow(k[0]+k[1],6))/(2.*3.*4.*5.*6.);
	}
	
}


// =================================================
// Start Correction
// =================================================

void feSolutionContainerV2DC3F::computeStartCorrection(feSolutionContainerV2DC2F *solutionContainerDC2F, std::vector<double> &k)
{
	this->computeStartCorrectionCoeff(k);

	int nRange=3;
	std::vector<double> uDot2_t1(_nbDDL, 0.0);
	std::vector<double> uDot3_t1(_nbDDL, 0.0);

	if (_correctionType=="SOLDOT"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solDotDC2F = solutionContainerDC2F->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				uDot2_t1[j] += _cStartDC3_t1[i]*solDotDC2F[j];
				uDot3_t1[j] -= _bStartDC3_t1[i]*solDotDC2F[j];
			}
		}
	
	} else if (_correctionType=="RESIDUAL"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solResidualDC2F = solutionContainerDC2F->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				uDot2_t1[j] += _cStartDC3_t1[i]*solResidualDC2F[j];
				uDot3_t1[j] -= _bStartDC3_t1[i]*solResidualDC2F[j];
			}
		}

	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}

	for (int j=0; j<_nbDDL; j++){
		_d[j] = uDot2_t1[j]*k[0]/2.0 + uDot3_t1[j]*pow(k[0],2)/(2.*3.);
	}
	
}




void feSolutionContainerV2DC4F::computeStartCorrection(feSolutionContainerV2DC3F *solutionContainerDC3F, std::vector<double> &k)
{
	this->computeStartCorrectionCoeff(k);

	_dDC4_t1.resize(_nbDDL,0.);
	_dDC4_t2.resize(_nbDDL,0.);

	int nRange=4;
	std::vector<double> uDot2_t1(_nbDDL, 0.0);
	std::vector<double> uDot3_t1(_nbDDL, 0.0);
	std::vector<double> uDot4_t1(_nbDDL, 0.0);

	std::vector<double> uDot2_t2(_nbDDL, 0.0);
	std::vector<double> uDot3_t2(_nbDDL, 0.0);
	std::vector<double> uDot4_t2(_nbDDL, 0.0);

	if (_correctionType=="SOLDOT"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solDotDC3F = solutionContainerDC3F->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				uDot2_t1[j] += _cStartDC4_t1[i]*solDotDC3F[j];
				uDot3_t1[j] -= _bStartDC4_t1[i]*solDotDC3F[j];
				uDot4_t1[j] += _aStartDC4_t1[i]*solDotDC3F[j];

				uDot2_t2[j] += _cStartDC4_t2[i]*solDotDC3F[j];
				uDot3_t2[j] -= _bStartDC4_t2[i]*solDotDC3F[j];
				uDot4_t2[j] += _aStartDC4_t2[i]*solDotDC3F[j];
			}
		}
	
	} else if (_correctionType=="RESIDUAL"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solResidualDC3F = solutionContainerDC3F->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				uDot2_t1[j] += _cStartDC4_t1[i]*solResidualDC3F[j];
				uDot3_t1[j] -= _bStartDC4_t1[i]*solResidualDC3F[j];
				uDot4_t1[j] += _aStartDC4_t1[i]*solResidualDC3F[j];

				uDot2_t2[j] += _cStartDC4_t2[i]*solResidualDC3F[j];
				uDot3_t2[j] -= _bStartDC4_t2[i]*solResidualDC3F[j];
				uDot4_t2[j] += _aStartDC4_t2[i]*solResidualDC3F[j];
			}
		}

	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}

	for (int j=0; j<_nbDDL; j++){
		_dDC4_t1[j] = uDot2_t1[j]*k[0]/2.0 + uDot3_t1[j]*pow(k[0],2)/(2.*3.) + uDot4_t1[j]*pow(k[0],3)/(2.*3.*4.);
		_dDC4_t2[j] = uDot2_t2[j]*k[0]/2.0 + uDot3_t2[j]*pow(k[0],2)/(2.*3.) + uDot4_t2[j]*pow(k[0],3)/(2.*3.*4.);
	}
	
}



void feSolutionContainerV2DC5F::computeStartCorrection(feSolutionContainerV2DC4F *solutionContainerDC4F, std::vector<double> &k)
{
	this->computeStartCorrectionCoeff(k);

	_dDC5_t1.resize(_nbDDL,0.);
	_dDC5_t2.resize(_nbDDL,0.);
	_dDC5_t3.resize(_nbDDL,0.);

	int nRange=5;
	std::vector<double> uDot2_t1(_nbDDL, 0.0);
	std::vector<double> uDot3_t1(_nbDDL, 0.0);
	std::vector<double> uDot4_t1(_nbDDL, 0.0);
	std::vector<double> uDot5_t1(_nbDDL, 0.0);

	std::vector<double> uDot2_t2(_nbDDL, 0.0);
	std::vector<double> uDot3_t2(_nbDDL, 0.0);
	std::vector<double> uDot4_t2(_nbDDL, 0.0);
	std::vector<double> uDot5_t2(_nbDDL, 0.0);

	std::vector<double> uDot2_t3(_nbDDL, 0.0);
	std::vector<double> uDot3_t3(_nbDDL, 0.0);
	std::vector<double> uDot4_t3(_nbDDL, 0.0);
	std::vector<double> uDot5_t3(_nbDDL, 0.0);

	if (_correctionType=="SOLDOT"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solDotDC4F = solutionContainerDC4F->getSolDot(i);
			for (int j=0; j<_nbDDL; j++){
				uDot2_t1[j] += _cStartDC5_t1[i]*solDotDC4F[j];
				uDot3_t1[j] -= _bStartDC5_t1[i]*solDotDC4F[j];
				uDot4_t1[j] += _aStartDC5_t1[i]*solDotDC4F[j];
				uDot5_t1[j] -= _eStartDC5_t1[i]*solDotDC4F[j];

				uDot2_t2[j] += _cStartDC5_t2[i]*solDotDC4F[j];
				uDot3_t2[j] -= _bStartDC5_t2[i]*solDotDC4F[j];
				uDot4_t2[j] += _aStartDC5_t2[i]*solDotDC4F[j];
				uDot5_t2[j] -= _eStartDC5_t2[i]*solDotDC4F[j];

				uDot2_t3[j] += _cStartDC5_t3[i]*solDotDC4F[j];
				uDot3_t3[j] -= _bStartDC5_t3[i]*solDotDC4F[j];
				uDot4_t3[j] += _aStartDC5_t3[i]*solDotDC4F[j];
				uDot5_t3[j] -= _eStartDC5_t3[i]*solDotDC4F[j];
			}
		}
	
	} else if (_correctionType=="RESIDUAL"){
		for (int i=0; i<nRange; i++){
			std::vector<double> &solResidualDC4F = solutionContainerDC4F->getResidual(i);
			for (int j=0; j<_nbInc; j++){
				uDot2_t1[j] += _cStartDC5_t1[i]*solResidualDC4F[j];
				uDot3_t1[j] -= _bStartDC5_t1[i]*solResidualDC4F[j];
				uDot4_t1[j] += _aStartDC5_t1[i]*solResidualDC4F[j];
				uDot5_t1[j] -= _eStartDC5_t1[i]*solResidualDC4F[j];

				uDot2_t2[j] += _cStartDC5_t2[i]*solResidualDC4F[j];
				uDot3_t2[j] -= _bStartDC5_t2[i]*solResidualDC4F[j];
				uDot4_t2[j] += _aStartDC5_t2[i]*solResidualDC4F[j];
				uDot5_t2[j] -= _eStartDC5_t2[i]*solResidualDC4F[j];

				uDot2_t3[j] += _cStartDC5_t3[i]*solResidualDC4F[j];
				uDot3_t3[j] -= _bStartDC5_t3[i]*solResidualDC4F[j];
				uDot4_t3[j] += _aStartDC5_t3[i]*solResidualDC4F[j];
				uDot5_t3[j] -= _eStartDC5_t3[i]*solResidualDC4F[j];
			}
		}

	} else {
		feErrorMsg(FE_STATUS_ERROR,"Unknow solverID for computeCorrection");
	}

	for (int j=0; j<_nbDDL; j++){
		_dDC5_t1[j] = uDot2_t1[j]*k[0]/2.0 + uDot3_t1[j]*pow(k[0],2)/(2.*3.) + uDot4_t1[j]*pow(k[0],3)/(2.*3.*4.) + uDot5_t1[j]*pow(k[0],4)/(2.*3.*4.*5.);
		_dDC5_t2[j] = uDot2_t2[j]*k[0]/2.0 + uDot3_t2[j]*pow(k[0],2)/(2.*3.) + uDot4_t2[j]*pow(k[0],3)/(2.*3.*4.) + uDot5_t2[j]*pow(k[0],4)/(2.*3.*4.*5.);
		_dDC5_t3[j] = uDot2_t3[j]*k[0]/2.0 + uDot3_t3[j]*pow(k[0],2)/(2.*3.) + uDot4_t3[j]*pow(k[0],3)/(2.*3.*4.) + uDot5_t3[j]*pow(k[0],4)/(2.*3.*4.*5.);
	}
	
}






// ================================================
// initialize solution
// ================================================

void feSolutionContainerV2::initializeBDFSolution(feMesh *mesh, feMetaNumber *metaNumber, feSolution *currentSolution)
{
  for(int i = 0; i < _nbDDL; ++i) {
    currentSolution->setSolAtDOF(i, _sol[0][i]);
    // currentSolution->setSolDotAtDOF(i, _solDot[0][i]);
  }
  
  currentSolution->setC0(_cBDF[0]);
  currentSolution->initializeContainerEssentialBC(mesh, metaNumber, this);
}

void feSolutionContainerV2::initializeDCSolution(feMesh *mesh, feMetaNumber *metaNumber, feSolution *currentSolution, std::string codeIniDC)
{
	if (codeIniDC=="fromPreviousSolution"){
		feInfo("initializing DC with the previous solution");
	  	for(int i = 0; i < _nbDDL; ++i) {
	    		currentSolution->setSolAtDOF(i, _sol[0][i]);
	    		// currentSolution->setSolDotAtDOF(i, _solDot[0][i]);
	  	}

	} else if (codeIniDC=="fromDC-1"){
		feInfo("initializing DCn with the previous DC(n-1) solution");
		this->copyCurrentSolution(currentSolution);
	
	} else {
		feWarning("Unknow codeIniDC => by default 'fromPreviousSolution'");
		for(int i = 0; i < _nbDDL; ++i) {
	    		currentSolution->setSolAtDOF(i, _sol[0][i]);
	    		// currentSolution->setSolDotAtDOF(i, _solDot[0][i]);
	  	}
	}
 	currentSolution->initializeContainerEssentialBC(mesh, metaNumber, this);
}




// =========================================================
// Set Correction 
// =========================================================

void feSolutionContainerV2DC4F::setStartCorrection(int step)
{
	if (step == 0){
		_d = _dDC4_t1;
	}
	else if (step == 1){
		_d = _dDC4_t2;
	}
}


void feSolutionContainerV2DC5F::setStartCorrection(int step)
{
	if (step == 0){
		_d = _dDC5_t1;
	}
	else if (step == 1){
		_d = _dDC5_t2;
	}
	else if (step == 2){
		_d = _dDC5_t3;
	}
}






