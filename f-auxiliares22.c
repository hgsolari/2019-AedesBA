// Copyright by Victoria Romeo Aznar <vromeoaznar@gmail.com>
// Lucas Alonso <lucasalo28@gmail.com> and
// Hernán G Solari < hgsolari@gmail.com> under License
// GNU GENERAL PUBLIC LICENSE (see CPYING file)
//
/* save the data after a day */
/* SACAR EN savedata LA PARTE DE FW->SOBRA*/

void savedata(FILE *FilePobla[], FILE *FileEventos[], int Plist[], int Elist[],
	float time, float Temp, float lluvia, FoodAndWeight FW)	
{
int err=0,i,j;


  /*  for (i=0; i< LIMfilas ; i++)
       for (j=0; j< LIMcolumnas; j++){
           aux_l=0; aux_p=0; aux_ml=0;
           for(k=0; k<LarPobla; k++){
               aux_l+=pobla[i][j][k+3];
               aux_ml+=neventos[i][j][2*(k+3)];
           }
           mortalidades_L[i][j]=aux_ml;
           larvas[i][j]=aux_l;
           for(k=0; k<PupPobla; k++){
               aux_p+=pobla[i][j][LarPobla+3+k]; 
               aux_mp+=neventos[i][j][2*(k+LarPobla+3)];
            }
           mortalidades_P[i][j]=aux_mp;
           pupas[i][j]=aux_p;
       }*/
          
while ((Plist[err] >= 0) && (err < (POBLACIONES)) )
//  while ((Plist[err] >= 0) && (err < 8) )
	{ 
 	fprintf(FilePobla[err], "%f %f %f \n", time, Temp, lluvia);
	for (i=0; i< LIMfilas ; i++)
		{for (j=0; j< LIMcolumnas; j++)
                     //if((err<3) || (4 <err))
			fprintf(FilePobla[err], "%ld\t",pobla[i][j][err]);
                /*     else {
                       if(err==3)
                         fprintf(FilePobla[err], "%ld\t",larvas[i][j]);
                       if(err==4)
                         fprintf(FilePobla[err], "%ld\t",pupas[i][j]);
                     }*/

		fprintf(FilePobla[err], "\n"); 
		}
	err++ ; 
	}


 	fprintf(FilePobla[err], "%f %f %f \n", time, Temp, lluvia);
	for (i=0; i< LIMfilas ; i++)
		{for (j=0; j< LIMcolumnas; j++)
			fprintf(FilePobla[err], "%f\t",FW.sobra[i][j]);
		fprintf(FilePobla[err], "\n"); 
		}

 	fprintf(FilePobla[err+1], "%f %f %f \n", time, Temp, lluvia);
	for (i=0; i< LIMfilas ; i++)
		{for (j=0; j< LIMcolumnas; j++)
			fprintf(FilePobla[err+1], "%f\t",FW.produce[i][j]);
		fprintf(FilePobla[err+1], "\n"); 
		}

err=0;  
while ((Elist[err] >= 0) && (err < EVENTOS+7) )
//while ((Elist[err] >= 0) && (err < 24) )
	{
	fprintf(FileEventos[err], "%f %f %f \n", time, Temp, lluvia);
	for (i=0; i< LIMfilas ; i++)
		{for (j=0; j< LIMcolumnas; j++)
                     //if((err<6) || (err>9))
			fprintf(FileEventos[err], "%ld\t",neventos[i][j][err]);
                     /*else{
                       if(err==6)
                         fprintf(FileEventos[err], "%ld\t",mortalidades_L[i][j]);
                       if(err==7)
                         fprintf(FileEventos[err], "%ld\t",neventos[i][j][2*(LarPobla+2)+1]);
                       if(err==8)
                         fprintf(FileEventos[err], "%ld\t",mortalidades_P[i][j]);
                       if(err==9)
                         fprintf(FileEventos[err], "%ld\t",neventos[i][j][2*(LarPobla+PupPobla+2)+1]);
                     }*/
                          
		fprintf(FileEventos[err], "\n");
		}
	err++ ; 
	} 
}



int getdata(FILE *FilePobla[], FILE *FileEventos[], int Plist[], int Elist[],
	char *mastername, char *poblaname[], char *eventoname[],//
	int *Transitorio, int *Tsimulado, long *idum, float temperaturavec[], float lluviavec[], float hlsvec[])
{
int Npobla=0, Neventos=0;
int err=0,i;
char aux2[10]="NADA", aux1[50]="NADA";
FILE *inputfile;

inputfile=fopen(AEDESINPUT,"r");

if(inputfile == (FILE *) NULL)
	return -1; /* error of reading the file */

fscanf(inputfile,"%s",mastername);
//fprintf(stderr,"Lei %s",mastername);

fscanf(inputfile,"%d %d  %f", Transitorio, Tsimulado, &C.BS);

fscanf(inputfile,"%d %d",&Npobla,&Neventos);
fprintf(stderr,"Lei %s\n%d\t%d\t%f\n",mastername,Npobla,Neventos,C.BS);

/* Open the population files */
if ((Npobla < 0) || (Npobla > POBLACIONES) )
	return err=-2;
else
	{for (i=0; i< Npobla; i++)
		fscanf(inputfile,"%d", &Plist[i]);
	} 
	for (i=Npobla; i< POBLACIONES; i++)
		Plist[i]=-1;

for (i=0; i < Npobla ; i++) /* mientras este dentro de los leidos */
	{
	strcpy(aux1,mastername);/* copy the master name to the auxiliary */
	
	sprintf(aux2,"-P%d.dat",Plist[i]); /* I add the indicator */
	strcat(aux1,aux2); /* I add the indicator */
	strcpy(poblaname[i],aux1); /* save the name */
	FilePobla[i]=fopen(poblaname[i],"w"); /* pass the pointer to the filepointers vector */
	};
/* for the food */
	strcpy(aux1,mastername);
	strcat(aux1,"-Fsobra");
	strcpy(poblaname[Npobla],aux1);
	FilePobla[Npobla]=fopen(poblaname[Npobla],"w");
	strcpy(aux1,mastername);
	strcat(aux1,"-Fproduce");
	strcpy(poblaname[Npobla+1],aux1);
	FilePobla[Npobla+1]=fopen(poblaname[Npobla+1],"w");
/* End of food */


/* the event section */
if ((Neventos < 0) || (Neventos > EVENTOS+7) )
	return err=-3;
else
	{for (i=0; i< Neventos; i++)
		fscanf(inputfile,"%d", (Elist+i));
	}
	for (i=Neventos; i< EVENTOS+7; i++)
		Elist[i]=-1;

for (i=0; i < Neventos ; i++) /* mientras este dentro de los leidos */
	{
	strcpy(aux1,mastername);/* copy the master name to the auxiliary */
	
	sprintf(aux2,"-E%d.dat",Elist[i]);
	strcat(aux1, aux2); /* I add the indicator */
	strcpy(eventoname[i],aux1); /* save the name */
	FileEventos[i]=fopen(eventoname[i],"w"); /* pass the pointer to the filepointers vector*/
	};

fclose(inputfile);

if ((inputfile=fopen("semilla.inp","r") )!= NULL)
	fscanf(inputfile, "%ld", idum);
else
	{
	printf("preciso semilla ");
	scanf("%ld", idum);
	}


/* reading of Horas de luz solar */
int hls_i;
FILE *HorasLuzSolar;
HorasLuzSolar=fopen("HLS.txt","r");
//for(hls_i=0;hls_i<=364;hls_i++)fscanf(HorasLuzSolar,"%f",&hls[hls_i]);
for(hls_i=0;hls_i<=*Tsimulado + *Transitorio;hls_i++)fscanf(HorasLuzSolar,"%f",&hlsvec[hls_i]);
fclose(HorasLuzSolar);  


/* reading of rains and temperatures */
int filatemp;
char entrada[80];
FILE *tempylluv;  /*The file in two columns. The first column with the mean temperature and the other one with the daily rainfall.*/
tempylluv = fopen("TempyLL.dat", "r"); //TempyLL_TempBelg.dat","r");
if(tempylluv == (FILE *) NULL)
	return -4; /* error of reading the file */
for (filatemp=0; filatemp < *Tsimulado + *Transitorio; filatemp++) {
    if(NULL == fgets(entrada,80,tempylluv)){
		    fprintf(stderr,"NO MORE WEATHER DATA\n"); exit(-1);}
    if(2 != sscanf(entrada,"%f%f",&temperaturavec[filatemp], &lluviavec[filatemp]))
    {fprintf(stderr,"Error en Temperatura y lluvia linea %d\n",filatemp+1);
		    exit(-1);
    };
}

fclose(tempylluv); 
return err;
}

void save_transitorio(long p[][LIMcolumnas][POBLACIONES], FoodAndWeight FW) {
/* p is the pointer to de poblational vector */
   
   FILE *tedio;
   int i,j,k;

   if((tedio=fopen("CIniciales.dat","r")) != NULL) /* file exists */
	fprintf(stderr,"CIniciales existia\n");
	
   else { // start else
     fprintf(stderr, "Guardo transitorio  de %d días\n", TsaveCI);
     tedio=fopen("CIniciales.dat","w");
     for (i=0;i<LIMfilas; i++) { // start firs i
        for(j=0;j<LIMcolumnas;j++)
	        for(k=0; k< POBLACIONES; k++)
	     fprintf(tedio,"%ld ", p[i][j][k]);
        fprintf(tedio,"\n");
       }; // close first i
     for (i=0;i<LIMfilas; i++) { // start second i
        for(j=0;j<LIMcolumnas;j++)
          fprintf(tedio,"%f ", FW.produce[i][j]);
        fprintf(tedio,"\n");
     } // close second i
     for (i=0;i<LIMfilas; i++) { // start third i
        for(j=0;j<LIMcolumnas;j++)
          for(k=3; k<POBLACIONES-3; k++) /*I'm saving the adults weight too. Eggs not because they have always the same weight*/
            fprintf(tedio,"%f ", FW.peso[i][j][k]);
        fprintf(tedio,"\n");
     } // end of third i
     fclose(tedio); // I close the archive CIniciales.dat
   }; // end of else
} // end of the function save_transitorio



void nuevasemilla(long idum) {
   FILE *tedio;
   tedio=fopen("semilla.inp","w");
   fprintf(tedio,"%ld\n", idum);
   fclose(tedio);
}



void pobla_iniciales(int *Transitorio, int *tiempo_inicial, FoodAndWeight *FW) {
  FILE *tedio;
  int i,j,k;
  float aux2, aux3, alfa;  

  alfa=C.alfa; 

  if((tedio=fopen("CIniciales.dat","r")) == NULL) {
    fprintf(stderr,"Transitorio largo\n");
    *tiempo_inicial= TsaveCI;
    for (i=0;i<LIMfilas; i++)
      for(j=0;j < LIMcolumnas;j++) {
        for(k=1; k < POBLACIONES; k++)
          pobla[i][j][k]=0;
        pobla[i][j][0]= HInicialesLD; // huevos LD al inicio de la simulación
        pobla[i][j][95]= HInicialesSD;// Huevos SD al inicio de la simulación
        FW->produce[i][j]=1.0;
       
        for(k=3; k<POBLACIONES-3; k++)        
          FW->peso[i][j][k]=0.;
      }; /* initialize all in zero less the eggs */
    }
  /* reading the data of tedio */ 
  else {
    *tiempo_inicial=0.;
    *Transitorio= *Transitorio+TsaveCI;
    /* reading the initial population */
    for (i=0;i<LIMfilas; i++){
      for(j=0;j<LIMcolumnas;j++){
	    for(k=0; k< POBLACIONES; k++)
          fscanf(tedio,"%ld", pobla[i][j]+k);
      }
    };
    /* reading the initial food */  
    for(i=0; i<LIMfilas; i++) 
      for(j=0; j<LIMcolumnas; j++)
        fscanf(tedio,"%f",&(FW->produce[i][j]));
    /* reading the inicial pobla weight */
    for (i=0;i<LIMfilas; i++){
      for(j=0;j<LIMcolumnas;j++){
	    for(k=3; k< POBLACIONES-3; k++) 
          fscanf(tedio,"%f", &(FW->peso[i][j][k]));//podes poner &(FW->peso[i][j][k]) o FW->peso[i][j]+k, de otro modo me salta error, violacion de segmento
      }
    };

  }; // end of else
  
        
  if((tedio != NULL))
    fclose(tedio); /* closing tedio */

/* giving the eggs weight */
 /* for(i=0; i<LIMfilas; i++)
    for(j=0;j<LIMcolumnas; j++)
      for(k=0; k<3; k++)
        FW->peso[i][j][e]=C.peso_huevos;*/


/* Calculating the first condition for the term beta[i][j] */
  for(i=0; i<LIMfilas; i++)
     for(j=0; j<LIMcolumnas; j++) {
       for(k=0; k<3; k++)
         FW->peso[i][j][k]=C.peso_huevos; /* giving the eggs weight*/ /*ya está normalizado por Bmax en constantes*/
      // FW->peso[i][j][3]=C.peso_huevos;
 
       for(k=LarPobla+PupPobla+6; k<LarPobla+PupPobla+9; k++)
         FW->peso[i][j][k]=C.peso_huevos; /* giving the eggs with diapausa weight*/ /*ya está normalizado por Bmax en constantes*/

        aux2=FW->produce[i][j];
        FW->beta[i][j]=LaBeta(aux2);
        for(k=3; k<(LarPobla+3); k++){
            aux3=FW->peso[i][j][k];
            FW->peso_alfa[i][j][k]=aux3;
            FW->dBdt[i][j][k]= FW->beta[i][j]*log(log(FW->peso[i][j][k]/
	C.peso_huevos)/C.RateTop);
        }
        aux3=FW->peso[i][j][LarPobla+PupPobla+4];
        FW->Fecun[i][j]=fecundidad(aux3);
     }

}/* returning with the initial populations */   
