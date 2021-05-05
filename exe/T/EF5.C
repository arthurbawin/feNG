#include "EF5.h"
#include "EF5_PRBSOL.h"
#include "EF5_PRBSTA.h"
#include "EF5_PRBBDF.h"
//#include "EF5_EF2TP.h"



int main( int argc, char** argv )
{
	int ret = 0;

	if(argc!=2) exit(-1);
 
// ===============================================================================
// 
// ===============================================================================
	EF5_PRBSOL*    PRBSOL     = new EF5_PRBSOL(argv[1]);
	EF5_PRBBDF*    PRBBDF     = new EF5_PRBBDF(PRBSOL);
	PRBBDF->RESOUDRE();
	delete PRBBDF;
	delete PRBSOL;

   	return ret;
}
