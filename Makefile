CC = gcc -I/PH.D/water/flexible_water/INM_gromacs_MD/Utility -I. -L/PH.D/water/flexible_water/INM_gromacs_MD/Utility 
INM_RF.exe :  INM.o tred2.o tqli.o IOUtility.o pythag.o Potential.o
	$(CC) -o INM_RF.exe  INM.o tred2.o tqli.o IOUtility.o pythag.o Potential.o -lm 
INM.o : INM.c
	$(CC) -c INM.c
tred2.o : tred2.c
	$(CC) -c tred2.c
tqli.o : tqli.c
	$(CC) -c tqli.c
IOUtility.o : IOUtility.c
	$(CC) -c IOUtility.c
pythag.o : pythag.c
	$(CC) -c pythag.c
Potential.o : Potential.c
	$(CC) -c Potential.c
