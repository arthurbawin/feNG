#include "EF5.h"
#include "EF5_PRBSOL.h"
#include "EF5_PRBSTA.h"
#include "EF5_PRBBDF.h"
#include "EF5_PRBVSU2VU.h"
#include "EF5_PRBVSU2TP.h"
#include "EF5_VTK.h"

typedef enum {VSUVTK, VSUTP, VSUVU} EF5_TYPEVSU;

EF5_TYPEVSU GETTYPEVSU(const char* Fichier)
{
	EF5_TYPEVSU typevsu=VSUTP;
	FILE* fp = fopen(Fichier,"r");
		char token[maxchar];

		fscanf(fp,"%s",token);
		const char* p = strstr(token,":");

		if(p!=NULL)
		{
			p++;

			if(strncmp(p,"EF2VTK",strlen("EF2VTK"))==0) typevsu=VSUVTK;
			if(strncmp(p,"EF2TP",strlen("EF2TP"))==0)   typevsu=VSUTP;
			if(strncmp(p,"EF2VU",strlen("EF2VU"))==0)   typevsu=VSUVU;
		}
		// printf("typevsu %ld\n",(EFint) typevsu);
	fclose(fp);
	return typevsu;
};

int main( int argc, char** argv )
{
	int ret = 0;
	assert (argc == 2);

	EF5_TYPEVSU typevsu = GETTYPEVSU(argv[1]);

	if(typevsu==VSUVTK) EF5_VTK VTK(argv[1]); 
	if(typevsu==VSUTP)  EF5_PRBVSU2TP TP(argv[1]);
	if(typevsu==VSUVU)  EF5_VTK VU(argv[1]);
	 
	return ret;
}
