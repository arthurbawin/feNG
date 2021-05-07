# feNG
TODO : Écrire un README

	feNG/src contient les fichiers sources (.cpp) et les headers (.h)
	feNG/exe contient les codes sources associés aux exécutables

Pour créer un makefile et compiler le code :

Dans feNG/, créer un répertoire build (qui est privé : il n'est pas partagé sur GitHub) :
 	mkdir build
 	cd build

 Lancer CMake en pointant vers le répertoire parent (..). CMake s'occupe de localiser MPI, mais
 il faut indiquer où se trouve petsc via les variables PETSC_DIR et PETSC_ARCH :
 	cmake -DPETSC_DIR=/chemin/vers/petsc -DPETSC_ARCH=architectureDetecteeParPETSc ..

 Par exemple :
 	cmake -DPETSC_DIR=/home/arthur/Code/petsc -DPETSC_ARCH=arch-linux-c-debug ..

 Pour la version graphique de CMake (pratique pour afficher les détails) :
 	ccmake ..
 (puis "c" pour configurer, "g" pour generate and exit)

 Compiler dans feNG/build :
 	make