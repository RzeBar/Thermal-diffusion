#include "lib/nrutil.c"
#include "lib/gaussj.c"

#include <iostream>
#include <cmath>
#include <fstream>

/////////////////////////////////////////////////////////////////
//                                                             //
//       Projekt MOFiT C18 - Bartłomiej Rzeszotarski FT3       //
//       g++ -Wall -std=c++0x main.C #lambda                   //
//                                                             //
/////////////////////////////////////////////////////////////////

using namespace std;

int main (){
	
//______ INITIAL PARAMETERS _______//
	const double L1 = 3.75 ;// Bok kwadratowego procesora [cm]
	const double L2 = 12.5 ; // Bok kwadratowego radiatora [cm]
	const double D = 0.6 ; // Wspolczynnik przewodnictwa cieplnego
	const double T = 30; // Czas symulacji
	const double Temp = 100; // Temperatura procesora [*C]
	const double TempAir = 18; // Temperatura otoczenia [*C]
	const double h = 0.1; // wspolczynnik konwekcyjnej transmisji ciepla

	const int Nt = 60 ; // Ilosc krokow czasowych
	const int Nx = 25 ; // Ilosc punktow siatki w 1D 
	const int Ny = 25 ; // Ilosc punktow siatki w 1D 
	// UWAGA:
	// przyjac Tp jako zrodlo S(T-Tp)= - alpha dT/dn

//______ VARIABLES _______//
	const double dt = T/Nt;
	const double dx = L2/(Nx-1);
	const double dy = L2/(Ny-1);

	const double iProcLeft(1+(L2-L1)/2./dx),iProcRight(1+(L2+L1)/2./dx);
	const double jProcUp(1+(L2-L1)/2./dy),jProcDown(1+(L2+L1)/2./dy);
//cout<< iProcLeft <<" " << iProcRight <<" " << jProcUp <<" " <<jProcDown<<endl;
	const double alpha = D*dt/2./dx/dx ;
	const double beta = D*dt/2./dy/dy ;
	const double gamma = 1. + 2.*alpha + 2.*beta;	

	double f00,fp0,f0p,fm0,f0m; // wartosci funkcji w chwili n o punktach w "fij" gdzie p=X+1, m=X-1, ponaddto X=i,j, a 0=i,j
	double value;
	double time(0);
	int row,col;
	//File variables
	ofstream file;
	char fname[50];

	float **A,**templA,**B, **tmpB;

	A=matrix(1,Nx*Ny,1,Nx*Ny); // Macierz glowna A - do obliczen
	templA=matrix(1,Nx*Ny,1,Nx*Ny); // Szablon macierzy A
	B=matrix(1,Nx*Ny,1,1); // Macierz rozwiazan B - do obliczen
	tmpB=matrix(1,Nx*Ny,1,1); // Macierz bufor B




	// MATRIX TEMPLATE PREPARATION
 	for(int k = 1 ; k <= Nx; k++)
	for(int l = 1 ; l <= Ny; l++)
	for(int i = 1 ; i <= Nx; i++){
	for(int j = 1 ; j <= Ny; j++){
		col=(i-1)*Ny+j;
		row=(k-1)*Ny+l;

		value = 0.;

		if( k==1 || l==1 || k==Nx || l==Ny ){// warunek brzegowy 
				// BOKI
				if(k==Nx && l!=Ny && l!=1){
					if( i==k && j==l){ value = h*dx/D+1; } 
					if( l==j && k-i==1){ value = -1; } 
				}
				if(l==Ny && k!=Nx && k!=1){
					if( i==k && j==l){ value = h*dy/D+1; } 
					if( k==i && l-j==1){ value = -1; } 
				}
				if(k==1 && l!=Ny && l!=1){
					if( i==k && j==l){ value = h*dx/D+1; } 
					if( l==j && k-i==-1){ value = -1; } 
				}
				if(l==1 && k!=Nx && k!=1){
					if( i==k && j==l){ value = h*dy/D+1; } 
					if( k==i && l-j==-1){ value = -1; }
				}
				// ROGI
				if(k==Nx && l==Ny && l!=1){
					if( i==k && j==l){ value = h*sqrt(dx*dx+dy*dy)/D/*+h*dx/D*TempAir+h*dy/D*TempAir+3*/; } 
					if( k-i==1 && l-j==1){ value = -1; } 
//					if( k-i==1 && l-j==2 || l-j==1 && k-i==2 ){ value = -1; } 
				}
				if(k==Nx && l!=Ny && l==1){
					if( i==k && j==l){ value = h*sqrt(dx*dx+dy*dy)/D/*+h*dx/D*TempAir+h*dy/D*TempAir+3*/; } 
					if( k-i==1 && l-j==-1){ value = -1; } 
//					if( k-i==1 && l-j==-2 || l-j==-1 && k-i==2 ){ value = -1; } 
				}
				if(k==1 && l!=Ny && l==1){
					if( i==k && j==l){ value = h*sqrt(dx*dx+dy*dy)/D/*+h*dx/D*TempAir+h*dy/D*TempAir+3*/; } 
					if( k-i==-1 && l-j==-1){ value = -1; } 
//					if( k-i==-1 && l-j==-2 || l-j==-1 && k-i==-2 ){ value = -1; } 
				}
				if(k==1 && l==Ny && l!=1){
					if( i==k && j==l){ value = h*sqrt(dx*dx+dy*dy)/D/*+h*dx/D*TempAir+h*dy/D*TempAir+3*/; } 
					if( k-i==-1 && l-j==1){ value = -1; }
//					if( k-i==-1 && l-j==2 || l-j==1 && k-i==-2 ){ value = -1; } 
				}
				

		}
		else{
			if( k>=iProcLeft && k<=iProcRight && l>=jProcUp && l<=jProcDown){// Styk
				if(i==k && j==l){value=1;}
				else{value=0.;}
			} 
			else{					// poza stykiem
				if( i==k && j==l ){ value = gamma; } // Diagonal
				if( k==i && abs(j-l)==1 ){ value = -beta; } // yDim
				if( l==j && abs(k-i)==1 ){ value = -alpha; } // xDim

			}
		}

		templA[row][col]=value;

	}	
	}
	// VECTOR B INITIATION
	for(int i = 1 ; i <= Nx*Ny; i++){ 
		B[i][1]=24; 
	} 

	//Print first
	if(Ny*Nx<=64)
	for(int i = 1 ; i <= Ny*Nx; i++){

		for(int j = 1 ; j <= Ny*Nx; j++){
			cout<<templA[i][j]<<" " ;

		}	
		//cout << B[i][1];
		cout << endl;
	}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^     MAIN LOOP  start  ^^^^
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	while(time<=T){
	sprintf(fname, "results/%f_out.dat",time);
	file.open(fname);	

	if(file.is_open()){
		// Zapis kroku
		for(int i = 1 ; i <= Nx; i++){
			for(int j = 1 ; j <= Ny; j++){
				row = (i-1)*Ny+j;
				file<<i*dx<<" "<<j*dy<<" "<<B[row][1]<<endl;
			}	
			file<<endl;
		}

		for(int i = 1 ; i <= Ny*Nx; i++){
			for(int j = 1 ; j <= Ny*Nx; j++){
				A[i][j]=templA[i][j];
			}	
		}
		// CORE
		gaussj(A,Nx*Ny,B,1);
		

		// 

		// Kopiowanie wektora odpowiedzi
		for(int i = 1 ; i <= Nx; i++){
			for(int j = 1 ; j <= Ny; j++){
				row = (i-1)*Ny+j;
				tmpB[row][1]=B[row][1];
			}	
		}
		// Przygotowanie wektora do dalszych obliczen
		for(int i = 1 ; i <= Nx; i++){
			for(int j = 1 ; j <= Ny; j++){
				row = (i-1)*Ny+j;


				if( i==1 || j==1 || i==Nx || j==Ny ){// warunek brzegowy
					value = sqrt(dx*dx+dy*dy)*h/D*TempAir/*+h*dx/D*TempAir+h*dy/D*TempAir*/; //rogi
					if( ( i==Nx || i==1 ) && j!=Ny && j!=1){  value=h*dx/D*TempAir;  } //brzegi
					if( ( j==Ny || j==1 ) && i!=Nx && i!=1){  value=h*dy/D*TempAir;  } //podstawy
				}
				else{
					if( i>=iProcLeft && i<=iProcRight && j>=jProcUp && j<=jProcDown){// Styk
						value = Temp;
					} 
					else{					// poza stykiem
						f00=tmpB[row][1];
						fm0=tmpB[(i-2)*Ny+j][1];
						fp0=tmpB[(i)*Ny+j][1];
						f0m=tmpB[(i-1)*Ny+j-1][1];
						f0p=tmpB[(i-1)*Ny+j+1][1];
						value = f00*(1.-2.*alpha-2.*beta) + alpha*(fp0+fm0) + beta*(f0p+f0m);
					}
				}




				B[row][1] = value;
			}	
		}

		file.close();
		time+=dt;

	}




	else{
		cout << "Error opening file" << endl;
		file.close();
		return 0;
	}
	}
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^     MAIN LOOP  end  ^^^^^^
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^









	free_matrix(A,1,Nx*Ny,1,Ny*Nx);
	free_matrix(templA,1,Nx*Ny,1,Ny*Nx);
	free_matrix(B,1,Nx*Ny,1,1);
	free_matrix(tmpB,1,Nx*Ny,1,1);

		cout<< "Symulacja zakończona powodzeniem" << endl; 

}
