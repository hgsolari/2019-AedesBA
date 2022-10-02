// Copyright by Victoria Romeo Aznar <vromeoaznar@gmail.com>
// Lucas Alonso <lucasalo28@gmail.com> and
// Hernán G Solari < hgsolari@gmail.com> under License
// GNU GENERAL PUBLIC LICENSE (see CPYING file)
//
/* promedio weights*/

float Promedio(int poblacion_e, int ne_1, int ne_2, float peso_e, float peso_e_1, float peso_e_2){
  float peso_e_prom;
  if(poblacion_e>0)
     peso_e_prom=((poblacion_e-ne_1-ne_2)*peso_e + ne_1*peso_e_1 + ne_2*peso_e_2)/poblacion_e;
  else
     peso_e_prom=0.0;
  return peso_e_prom;
}

int fecundidad(float peso){

    int fecun;    
    if (peso<pm){
    fecun = 0;}
    else{
    fecun = (int) (((peso)-pm)*C.fecundidad)/(1-pm);}
    return fecun;
}



/* We make the allocation of population increases */

void reparto_y_promedio (int i, int j, double lam[],long  paux[][LIMcolumnas][POBLACIONES], int count[], long *idum, float PesoAux[][LIMcolumnas][POBLACIONES], int FecunAux[LIMfilas][LIMcolumnas], int hls_o)
{int im, jm, ip, jp, e, cnt[EVENTOS]={0}, N1=0, N2=0, n3, n2;
double prob;
float PesoSinProm[POBLACIONES];
int paso3=0, pasoV1, pasoV2, paso3d;
paso3d=0;
pasoV1=0;
pasoV2=0;
FecunAux[i][j]=0;

/* defino PesoSinProm */
for(e=0; e<POBLACIONES; e++)
   PesoSinProm[e]=PesoAux[i][j][e];

/*the allocation of the positive events */
			 /* immature eggs */
			if(lam[0] > 0) /* careful with the division by zero */
				{prob=Lbd[i][j][0]/lam[0];
				N1=binomial(prob, count[0], idum); //deaths
				n2=N2=count[0]-N1; //others
				prob=C.eclo;
				n3=binomial(prob, n2,idum);// a larva
				N2=n2-n3;
				paux[i][j][1] += N2 ; // adding to  the mature
				paux[i][j][3] += n3; // adding to  the larvae
				neventos[i][j][0]+= N1;
				neventos[i][j][1]+= N2;
				neventos[i][j][3] += n3;
                paso3 += n3; /* here I retain all the new eggs that pass to larvae */
/* Event 3 is the spontaneous eclosion */
/*Deaths in N1, transitions in N2 */
				};
               /*Mature eggs*/
			if(lam[1] > 0)
				{prob=Lbd[i][j][2]/lam[1];
				N1=binomial(prob, count[1], idum);
				N2=count[1]-N1;
				paux[i][j][3] += N2 ;
				neventos[i][j][2]+= N1;
				/*if(N2 != 0)
					{neventos[i][j][3]+= N2;
					fprintf(stderr,"esto debia ser cero N2=%d i=%d j=%d\n",N2,i,j);
					}*/
                paso3 += N2; /* here I retain all the new eggs that pass to larvae. At the moment N2 should be 0 */
				}
               /*rain eggs*/
			if(lam[2] > 0)
				{prob=Lbd[i][j][4]/lam[2];
				N1=binomial(prob, count[2], idum);
				N2=count[2]-N1;
				paux[i][j][3] += N2 ;
				neventos[i][j][4]+= N1;
				neventos[i][j][5]+= N2;
                paso3 += N2; /* here I retain all the new eggs that pass to larvae */
 //if(i==1)
                                    //if(j==1)
                                        //printf("PesoAux[3]=%f, pobla[1][1][3]=%i\n", PesoAux[1][1][3], paux[1][1][3]);
				}
        

			/* larvae and pupae*/
			
                        for(e=3; e<(LarPobla+PupPobla+3);e++) {
				if((PesoSinProm[e] < pm) &
						(e==LarPobla+PupPobla+2)) break;
// Si la pupa peso menos que el peso requerido para ser fertil muere al emerger.
                                if(lam[e] > 0)
				  {prob=Lbd[i][j][2*e]/lam[e];
				  N1=binomial(prob, count[e], idum); /* number of larvae o pupae death, the e-esima*/
				  N2=count[e]-N1; /*number of individuals from e to e+1*/
				  paux[i][j][e+1] += N2 ;/* increment de e+1 population*/
				  neventos[i][j][2*e]+= N1; /*increment of death event*/
                          	  neventos[i][j][2*e+1]+= N2; /* increment in the event of transfer from e to e+1*/
                                  /*if(e<(LarPobla+3))
                                      neventos[i][j][3]+=neventos[i][j][2*e];*/
                                  /* average weight for larvae, pupae and adults 1 */
                                  PesoAux[i][j][e+1]=Promedio(paux[i][j][e+1], N2, 0, PesoAux[i][j][e+1], PesoSinProm[e], 0.0);
                                  //if(i==1)
                                    //if(j==1)
                                      //  printf("PesoAux[%i]=%f, pobla[1][1][%i]=%i\n", e+1, PesoAux[1][1][e+1], e+1, paux[1][1][e+1]);
                                  /*if(paux[i][j][e+1]>0)
                                     PesoAux[i][j][e+1]=((paux[i][j][e+1]-N2)*PesoAux[i][j][e+1]+N2*PesoSinProm[e])/paux[i][j][e+1];
                                  else
                                     PesoAux[i][j][e+1]=0.0;*/
				}
                        }
			
			/* Adults 1 */
			if(lam[LarPobla+PupPobla+3] > 0)
				{prob=Lbd[i][j][2*(LarPobla+PupPobla+3)]/lam[LarPobla+PupPobla+3];
				N1=binomial(prob, count[LarPobla+PupPobla+3], idum);
				N2=count[LarPobla+PupPobla+3]-N1;
				paux[i][j][LarPobla+PupPobla+4] += N2 ;
				neventos[i][j][2*(LarPobla+PupPobla+3)]+= N1;
				neventos[i][j][2*(LarPobla+PupPobla+3)+1]+= N2;
                                pasoV1=N2; /* here I retain all the new adults that will fly*/
				}

			/* flying */

			count[LarPobla+PupPobla+9]=0; /*ex-count[8]*/
			if(lam[LarPobla+PupPobla+4] > 0)
				{prob=Lbd[i][j][2*(LarPobla+PupPobla+4)]/lam[LarPobla+PupPobla+4];
				N1=binomial(prob, count[LarPobla+PupPobla+4], idum); /* die */
				N2=count[LarPobla+PupPobla+4]-N1; /* fly or lay eggs */ 
				neventos[i][j][2*(LarPobla+PupPobla+4)]+= N1;
				if (N2 >= 0) /* some fly  or lay eggs */
					{prob=Lbd[i][j][2*(LarPobla+PupPobla+9)]/(Lbd[i][j][2*(LarPobla+PupPobla+4)+1]+Lbd[i][j][2*(LarPobla+PupPobla+9)]);
					count[LarPobla+PupPobla+9]=binomial(prob,N2,idum); /* fly */
					N2 -= count[LarPobla+PupPobla+9]; /* lay eggs */
                                        if(N2>0)
                                           FecunAux[i][j]=fecundidad(PesoAux[i][j][LarPobla+PupPobla+4]);
                                        else
                                           FecunAux[i][j]=0;
					}
				paux[i][j][LarPobla+PupPobla+5] += N2 ;
				neventos[i][j][2*(LarPobla+PupPobla+4)+1]+= N2; //laid eggs
	 			paux[i][j][0] += FecunAux[i][j]*N2*hls_o; /* huevos */ // FECUNDIDAD SE LO PASO COMO PARAMETRO
                paux[i][j][LarPobla+PupPobla+6]+=FecunAux[i][j]*N2*(1-hls_o); /* paux[i][j][0] += FecunAux[i][j]*N2; // huevos  FECUNDIDAD SE LO PASO COMO PARAMETRO
                        averaging Adult1's weight */
               PesoAux[i][j][LarPobla+PupPobla+5]=Promedio(paux[i][j][LarPobla+PupPobla+5], N2, 0,PesoAux[i][j][LarPobla+PupPobla+5], PesoSinProm[LarPobla+PupPobla+4], 0.0);
				}

			/* Adults 2 */
			if(lam[LarPobla+PupPobla+5] > 0)
				{prob=Lbd[i][j][2*(LarPobla+PupPobla+5)]/lam[LarPobla+PupPobla+5];
				N1=binomial(prob, count[LarPobla+PupPobla+5], idum);
				N2=count[LarPobla+PupPobla+5]-N1;
				paux[i][j][LarPobla+PupPobla+4] += N2 ; /* suma a voladoras */
				neventos[i][j][2*(LarPobla+PupPobla+5)]+= N1;
				neventos[i][j][2*(LarPobla+PupPobla+5)+1]+= N2;
                                /*here I retain all the new adults that will fly*/
                                pasoV2 =N2;
				}


/* immature eggs diapausa */
			if(lam[LarPobla+PupPobla+6] > 0) /* careful with the division by zero */
				{prob=Lbd[i][j][2*(LarPobla+PupPobla+6)]/lam[LarPobla+PupPobla+6];
				N1=binomial(prob, count[LarPobla+PupPobla+6], idum); //deaths
				n2=N2=count[LarPobla+PupPobla+6]-N1; //others
				prob=C.eclo;                                    //REVISAR
				n3=binomial(prob, n2,idum);// a larva
				N2=n2-n3;
				paux[i][j][LarPobla+PupPobla+7] += N2 ; // adding to  the mature
				paux[i][j][3] += n3; // adding to  the larvae
				neventos[i][j][2*(LarPobla+PupPobla+6)]+= N1;
				neventos[i][j][2*(LarPobla+PupPobla+6)+1]+= N2;
				neventos[i][j][2*(LarPobla+PupPobla+7)+1] += n3;
                paso3d += n3; /* here I retain all the new eggs that pass to larvae */
/* Event 3 is the spontaneous eclosion */
/*Deaths in N1, transitions in N2 */
				};
               /*Mature eggs diapausa*/
			if(lam[LarPobla+PupPobla+7] > 0)
				{prob=Lbd[i][j][2*(LarPobla+PupPobla+7)]/lam[LarPobla+PupPobla+7];
				N1=binomial(prob, count[LarPobla+PupPobla+7], idum);
				N2=count[LarPobla+PupPobla+7]-N1;
				paux[i][j][3] += N2 ;
				neventos[i][j][2*(LarPobla+PupPobla+7)]+= N1;
				/*if(N2 != 0)
					{neventos[i][j][3]+= N2;
					fprintf(stderr,"esto debia ser cero N2=%d i=%d j=%d\n",N2,i,j);
					}*/
                paso3d += N2; /* here I retain all the new eggs that pass to larvae. At the moment N2 should be 0 */
				}
               /*rain eggs diapausa*/
			if(lam[LarPobla+PupPobla+8] > 0)
				{prob=Lbd[i][j][2*(LarPobla+PupPobla+8)]/lam[LarPobla+PupPobla+8];
				N1=binomial(prob, count[LarPobla+PupPobla+8], idum);
				N2=count[LarPobla+PupPobla+8]-N1;
				paux[i][j][3] += N2 ;
				neventos[i][j][2*(LarPobla+PupPobla+8)]+= N1;
				neventos[i][j][2*(LarPobla+PupPobla+8)+1]+= N2;
                paso3d += N2; /* here I retain all the new eggs that pass to larvae */
                }

                                

                      /* average weight for larvae 1, and the flying *//* These populations have immigration from more than one subpopulation */
                      //PesoAux[i][j][3]= Promedio(paux[i][j][3], paso3, 0, PesoAux[i][j][3], PesoSinProm[2],0.0);
                      PesoAux[i][j][3]= Promedio(paux[i][j][3], paso3, paso3d, PesoAux[i][j][3], PesoSinProm[2],PesoSinProm[LarPobla+PupPobla+8]);
                      PesoAux[i][j][LarPobla+PupPobla+4]=Promedio(paux[i][j][LarPobla+PupPobla+4], pasoV1, pasoV2, PesoAux[i][j][LarPobla+PupPobla+4], PesoSinProm[LarPobla+PupPobla+3],PesoSinProm[LarPobla+PupPobla+5]);
                     //  if(i==1)
                       //  if(j==1)
                         //   printf("paux[1][1]= %i \t  PesoAux[1][1]=%f \n",paux[1][1][LarPobla+PupPobla+4], PesoAux[1][1][LarPobla+PupPobla+4]
			/* finally the flight */
			if(count[LarPobla+PupPobla+9] > 0)
				{multinomial(0,7,count[LarPobla+PupPobla+9], C.pvuel, cnt, idum); // VER QUE SON EL RESTO DE LOS PARAMETR, PARA VER SI HAY Q MODIFICAR
		                im=i-1; jm=j-1; ip=i+1;jp=j+1;
       			        if( (im <= 0) ) im=0;
                		if( (jm <= 0) ) jm=0;
                		if( (jp >= LIMcolumnas) ) jp=LIMcolumnas-1;
                		if( (ip >= LIMfilas) ) ip=LIMfilas-1;
				paux[im][jm][LarPobla+PupPobla+4]+= cnt[0];
				paux[im][j][LarPobla+PupPobla+4]+= cnt[1];
				paux[im][jp][LarPobla+PupPobla+4]+= cnt[2];
				paux[i][jp][LarPobla+PupPobla+4]+= cnt[3];
				paux[ip][jp][LarPobla+PupPobla+4]+= cnt[4];
				paux[ip][j][LarPobla+PupPobla+4]+= cnt[5];
				paux[ip][jm][LarPobla+PupPobla+4]+= cnt[6];
				paux[i][jm][LarPobla+PupPobla+4]+= cnt[7];                          
PesoAux[im][jm][LarPobla+PupPobla+4]=Promedio(paux[im][jm][LarPobla+PupPobla+4],cnt[0],0,PesoAux[im][jm][LarPobla+PupPobla+4],PesoSinProm[LarPobla+PupPobla+4],0.0);
PesoAux[im][j][LarPobla+PupPobla+4]=Promedio(paux[im][j][LarPobla+PupPobla+4],cnt[1],0,PesoAux[im][j][LarPobla+PupPobla+4],PesoSinProm[LarPobla+PupPobla+4],0.0);
PesoAux[im][jp][LarPobla+PupPobla+4]=Promedio(paux[im][jp][LarPobla+PupPobla+4],cnt[2],0,PesoAux[im][jp][LarPobla+PupPobla+4],PesoSinProm[LarPobla+PupPobla+4],0.0);
PesoAux[i][jp][LarPobla+PupPobla+4]=Promedio(paux[i][jp][LarPobla+PupPobla+4],cnt[3],0,PesoAux[i][jp][LarPobla+PupPobla+4],PesoSinProm[LarPobla+PupPobla+4],0.0);
PesoAux[ip][jp][LarPobla+PupPobla+4]=Promedio(paux[ip][jp][LarPobla+PupPobla+4],cnt[4],0,PesoAux[ip][jp][LarPobla+PupPobla+4],PesoSinProm[LarPobla+PupPobla+4],0.0);
PesoAux[ip][j][LarPobla+PupPobla+4]=Promedio(paux[ip][j][LarPobla+PupPobla+4],cnt[5],0,PesoAux[ip][j][LarPobla+PupPobla+4],PesoSinProm[LarPobla+PupPobla+4],0.0);
PesoAux[ip][jm][LarPobla+PupPobla+4]=Promedio(paux[ip][jm][LarPobla+PupPobla+4],cnt[6],0,PesoAux[ip][jm][LarPobla+PupPobla+4],PesoSinProm[LarPobla+PupPobla+4],0.0);
PesoAux[i][jm][LarPobla+PupPobla+4]=Promedio(paux[i][jm][LarPobla+PupPobla+4],cnt[7],0,PesoAux[i][jm][LarPobla+PupPobla+4],PesoSinProm[LarPobla+PupPobla+4],0.0);
		
                                for(e=0; e<8; e++) 
                                  neventos[i][j][2*(LarPobla+PupPobla+9)+e]+= cnt[e];
       // printf("celda %i  %i \t paux[1][1]=%i \t PesoAux[1][1]= %f \n", i, j, paux[1][1][LarPobla+PupPobla+4], PesoAux[1][1][LarPobla+PupPobla+4]);
       // printf("cnt: 0= %i  1= %i  2= %i  3= %i  4= %i  5= %i  6= %i  7= %i\n", cnt[0],cnt[1],cnt[2],cnt[3],cnt[4],cnt[5],cnt[6],cnt[7]);


				}; /* if(count) end */

 /* esto es porque si lam[k]<0 entonces no le aplico el promedio y por lo tanto si la poblacion se hizo cero no le hice cero su peso, en principio esto no afectaria pues el peso solo tiene relevancia cuando se lo multiplica por la poblacion y ahi no aportaria entonces*/
  for(e=3; e<POBLACIONES-3; e++)
      if (paux[i][j][e]<1)
        PesoAux[i][j][e]=0.0;
   //for(e=148; e<POBLACIONES; e++)
    // if(i==1)
        //  if(j==1)
          //  printf("pobla[1][1][%i]=%i  \t  PesoAux[1][1][%i]=%f\n", e, paux[1][1][e], e, PesoAux[1][1][e]);


} /* END OF REPARTO FUNCTION*/



/* Here, changing the rain */
void llovio(long *idum,FoodAndWeight FW) { 
  unsigned long n;
  double p;
  int i,j;
  float comida_norm, comida_ProbEcloMin;

/* Gillett effect with values that mirror the increment of mortality */
  for (i=0; i< LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++) {
        comida_norm=FW.produce[i][j];
        p=GILLETTM; /* default for c > 1 */
        if((comida_norm < 1.) && (comida_norm > Lstarving))
	        p=GILLETTM+(GILLETTm-p)*(1.-comida_norm)/(1.-Lstarving); /* Ramp of increasing eclosion */
        if (comida_norm < Lstarving)
            p=GILLETTm; /* No enough food for survival */
    
     /*if(i==1)
        if(j==1)
          printf("Antes Hm %ld, Hll %ld \n", pobla[i][j][1], pobla[i][j][2]); */
      n=binomial(p, pobla[i][j][1], idum);     
      pobla[i][j][1]-= n; 
      pobla[i][j][2]+= n;

      n=binomial(p, pobla[i][j][LarPobla+PupPobla+7], idum);     
      pobla[i][j][LarPobla+PupPobla+7]-= n; 
      pobla[i][j][LarPobla+PupPobla+8]+= n;
      //neventos[i][j][0] += n;

      /*if(i==1)
        if(j==1)
          printf("p=%f , Hm %ld, Hll %ld, n %lu \n", p, pobla[i][j][1], pobla[i][j][2], n); */
     
    }
}




void updatepobla (long *idum, FoodAndWeight *FW, int hls_o)
{ 
  int i,j,k, count[EVENTOS]={0};

  long Etotales=0;
  long resto[LIMfilas][LIMcolumnas][POBLACIONES], restototal=0, oldresto=-1;
  long paux[LIMfilas][LIMcolumnas][POBLACIONES];
  double Ltotal, MAXl=0.0, lam[EVENTOS];
  float PesoAux[LIMfilas][LIMcolumnas][POBLACIONES];
  int FecunAux[LIMfilas][LIMcolumnas];

// VER DONDE VOY A CALCULAR LA FECUNDIDAD, ME PARECE QUE EN comidas.c Y LUEGO LO PASO AQUI

/* hago una copia de trabajo de las poblaciones */
/* hago una copia de trabajo de los pesos */
  for (i=0; i < LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++)
      for(k=0; k<POBLACIONES ; k++){
        paux[i][j][k]=pobla[i][j][k];
        PesoAux[i][j][k]=FW->peso[i][j][k];
      }

  for (i=0; i< LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++) {
      Ltotal=0.0;
      for (k=0; k< POBLACIONES; k++) {
        lam[k] = Lbd[i][j][2*k]+Lbd[i][j][2*k+1];
		resto[i][j][k]=0;
		Ltotal += lam[k];
			/* lam[k] tiene la suma de los lambdas que restan a las poblacion k, tal como esta listado en el main */
		MAXl= (MAXl >= Ltotal? MAXl : Ltotal);
      }  
      
      lam[LarPobla+PupPobla+4] += Lbd[i][j][2*(LarPobla+PupPobla+9)]; /* vuelo */
      lam[LarPobla+PupPobla+9] = Lbd[i][j][2*(LarPobla+PupPobla+9)];
      Ltotal += Lbd[i][j][2*(LarPobla+PupPobla+9)];
      
      Etotales= poidev(Ltotal, idum); /*total de eventos en celda*/
      multinomial(0, POBLACIONES-1, Etotales, lam, count, idum); //OJO ACA... ESTO TRAE PROBLEMAS, le saque el p[30] pero igual jode

/* en count tengo los eventos negativos para cada poblacion */
      for (k=0; k< POBLACIONES; k++) {
		if (count[k] > pobla[i][j][k]) {/* exceso temporario de eventos */
	  		resto[i][j][k]= count[k]-pobla[i][j][k];
	  		count[k]= pobla[i][j][k];
	  		paux[i][j][k]=0; /* poblacion a cero */
	  		restototal += resto[i][j][k];
		}
		else /* reste no mas */ {
          paux[i][j][k] -= count[k];
		};
      };
      
      reparto_y_promedio(i, j, lam, paux, count, idum, PesoAux, FecunAux, hls_o);
   }; //fin for de j e i

/* fuera del loop en i y j actualizo las poblaciones */
  for (i=0; i< LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++)
      for (k=0; k< POBLACIONES; k++)
		pobla[i][j][k]=paux[i][j][k];

  while ((restototal != oldresto) && (restototal > 0)) {
    oldresto=restototal; restototal=0;
    for (i=0; i< LIMfilas; i++)
      for (j=0; j< LIMcolumnas; j++) {
        for(k=0; k < POBLACIONES ; k++) {
          lam[k] = Lbd[i][j][2*k]+Lbd[i][j][2*k+1];
          count[k]=resto[i][j][k];
        };
        lam[LarPobla+PupPobla+4] += Lbd[i][j][2*(LarPobla+PupPobla+9)];
        lam[LarPobla+PupPobla+9] = Lbd[i][j][2*(LarPobla+PupPobla+9)];	
        for(k=0; k < POBLACIONES ; k++) {
          paux[i][j][k]=pobla[i][j][k];
          if (count[k] > pobla[i][j][k]) {
            resto[i][j][k]= count[k]-pobla[i][j][k];
	    	count[k]= pobla[i][j][k];
	    	paux[i][j][k]=0;
	    	restototal += resto[i][j][k];
          }
          else {
            paux[i][j][k] -= count[k];
  		  };
        }; //fin del for k's 2

        reparto_y_promedio(i, j, lam, paux, count, idum, PesoAux, FecunAux, hls_o);

      };// fin for j

      for (i=0; i< LIMfilas; i++)
        for (j=0; j< LIMcolumnas; j++)
		  for (k=0; k< POBLACIONES; k++)
	    	pobla[i][j][k]=paux[i][j][k];
  }; /* while loop */

 for (i=0; i < LIMfilas; i++)
    for (j=0; j< LIMcolumnas; j++){
      for(k=3; k<POBLACIONES-3 ; k++)
          FW->peso[i][j][k]=PesoAux[i][j][k];
      FW->Fecun[i][j]=FecunAux[i][j];
      }

}//END UPDATE

