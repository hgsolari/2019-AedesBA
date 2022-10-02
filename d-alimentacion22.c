// Copyright by Victoria Romeo Aznar <vromeoaznar@gmail.com>
// Lucas Alonso <lucasalo28@gmail.com> and
// Hernán G Solari < hgsolari@gmail.com> under License
// GNU GENERAL PUBLIC LICENSE (see CPYING file)
//
/* Respecto de la version final en la tesis de Vico, aquí modificamos que: No
 * hay reciclado de cadáveres en comida. Los valores son inciertos y produce
 * inestabilidades a veces.
 * Estaba mal el manejo de la mortalidad por falta de comida.
 * Se incopora tolerancia en la busqueda del
 * intervalo de comida en algunos casos que faltaba (epsilon) MiroIntervalo no
 * retorna nunca 0, saco el switch que aislaba ese caso i=0 corresponde a no
 * crecimiento y mortalidad aumentada.
 * Se modifica la exprsión del peso según apunte "peso.lox/pdf" de Septiembre
 * 2022
 */

#define rateW(x) x/(x*0.47315+0.033377) // tasa característica de la comida

float LaBeta(float comida) {
  float betaNorm, aux,  ratio;
  static float betaHungry=0.;

  ratio=comida;
/* Calcula el factor de crecimiento corregido cuando la comida es óptima
 * referirse a Comida-UltimaVersion.lyx donde están las ecuaciones
 * y peso.lyx */


        if (betaHungry == 0.)
                betaHungry=rateW(comida);
// si venimos con valor de inicialización convertimos a tasa correspondiente


  if(ratio>1.0)  
    betaNorm=1.0;
  else {
       if(ratio >= Lhungry)
	       betaNorm=rateW(ratio);
       else
	       betaNorm=0.; // no engordan en el régimen de mortalidad aumentada
    }


  return betaNorm;
}


int MiroIntervalo(float comida) {
   int aux, j, i;
   aux=1;
   j=1;
 // los intervalos son abiertos a derecha
   if(comida>=C.a[7]){ // Comida mayor que 1
      i=7;
      aux=-1; }

   while(aux>0){
      if(comida>=C.a[j]) // empiezo la busqueda en intervalo 
// Como es mayor que el limite inferior que pertenece al intervalo
// anterior debe ser el siguiente
         j=j+1;
      else { // si no está a la derecha, entonces era este que se llama j-1
         i=j-1;
         aux=-1;}
   }
   return i;
}
    

float IntegroComida(float Comensales, float comida, float rate_food, float dt, int s, int l, float t){ 
   float m, n, am, aM, mnl;
   float aux1, aux2,  Dt1, raciones, C0, C1, epsilon=0.001;
   float auxaux1, auxaux2;
   int i,a, aux3;
  
//   raciones=VOL*C.BS; esto fue un error. Una hibridación de ideas.
   raciones=1;
   C0= comida;
   i=MiroIntervalo(C0);
   /* doy valores a los parámetros de beta_aprox: * m la pendiente, n la ordenada al origen * el intervalo de validez [am;aM) */
   m=C.m[i]; n=C.n[i]; am=C.a[i]; aM=C.a[i+1];
   
   aux3=1; // Indico sin progreso
   Dt1=0.;
      
   while(aux3>0) { // aux3 > 0 es sin terminar
   aux3=-1; // Soy optimista
      
	if(i <= 1)// es aproximadamente la region de mortalidad aumentada
	mnl=SCAV*C.mnl;

	else
		mnl=0.;
      aux1=DECAE+Comensales*(m/raciones+mnl/(Lhungry-Lstarving));
      aux2=(rate_food-Comensales*n+mnl/(1.-Lstarving/Lhungry))/raciones;
  
      C1=(C0-aux2/aux1)*exp(-(dt-Dt1)*aux1)+aux2/aux1;
         
      if((C1<am) || (C1>=aM)) {
         if(C1<am){ /* estoy bajando de intervalo*/
	aux3=2; // pase al intervalo a la izquierda
	if(i <= 1)// es aproximadamente la region de mortalidad aumentada
		mnl=SCAV*C.mnl;
	else
		mnl=0.;
      i=i-1; 
// indico que pase de intervalo la izquierda
      C1=am;
            }
         else {
            if(am==1) // indico terminado
               aux3=-1; // Estoy en el último intervalo. Termine.
            else {// estoy a la derecha
	aux3=3; // pase al intervalo a la derecha
	if(i <= 1)// es aproximadamente la region de mortalidad aumentada
		mnl=SCAV*C.mnl;
	else
		mnl=0.;
      i=i+1; 
               C1=aM;
            };
         };
	}; // End del if (C1< am) || (C1<aM))

	if(aux3 > 1) // Sali por borde de intervalo
	 { 
         if((C0-aux2/aux1)/(C1-aux2/aux1)<1.-epsilon)
            {printf("fallan los tiempos de integracion en la comida %f \n",(C0-aux2/aux1)/(C1-aux2/aux1));
 		exit(0);
	    } // Aborte el programa
         Dt1 += log((C0-aux2/aux1)/(C1-aux2/aux1))/aux1;
//prevent roudoff error
	if((dt-Dt1) < 0)
		if(Dt1-dt < epsilon) // fue por poco, tolerable
			aux3=-1; // Termine
		else
			fprintf(stderr,"Me pase %g %f %i\n", Dt1-dt, C0, i);

      auxaux1=DECAE+Comensales*(C.m[i]/raciones+mnl/(Lhungry-Lstarving));
      auxaux2=(rate_food-Comensales*C.n[i]+mnl/(1.-Lstarving/Lhungry))/raciones;

	auxaux2=auxaux2/auxaux1;
	switch (aux3){
		case	1: fprintf(stderr,"Unexpect case %d\n",aux3);// esto es un error
		case	2: {if((auxaux2 - C1) > 0.) aux3=-1; // sali por la izquierda y me plante en el borde
			break;};
		case	3: {if((auxaux2 - C1) < 0.) aux3=-1; // sali por la derecha y me plante en el borde
			break;};
		default: {aux3=1; // sigo integrando
			break;};
		}

// calculo los nuevos rates

        m=C.m[i];
	n=C.n[i];
	am=C.a[i];
	aM=C.a[i+1];
        C0=C1;
      }
       else
         aux3=-1; // termine
   } //while
   return C1;
}

void comensales_y_alimento(FoodAndWeight *FW, float dt, float t){
  float  Comensales, comida, rate_food,alfa;
  // C.a es la tasa de captura de peso
  int i,j,e,k;
  float aux1, aux3;
  rate_food=FW->r_food;
  alfa=C.alfa;
  int auxx;
/* update the leftover food */
for (i=0;i<LIMfilas;i++)
  for (j=0;j<LIMcolumnas;j++) {
       	comida=FW->produce[i][j];
        FW->beta[i][j]=LaBeta(comida);

	 Comensales=0.;
         for(e=3; e<(LarPobla+3); e++) { //solo importan las larvas, pues solo ellas comen
	   if(FW->peso[i][j][e] > 0){
//                FW->peso_alfa[i][j][e] = pow(FW->peso[i][j][e],alfa);
              	Comensales += (pobla[i][j][e]*FW->peso[i][j][e]); 
	   };
	}
	 Comensales=6.*Comensales;
// los comensales comen en proporción a la comida Con temperatura de 27ºC y
// lluvias diarias de 7.6mm El peso en el estado estacionario es bastante
// constante con el factor que estamos // buscando. Con un factor 6. hay
// aproximadamente 30 larvas que es lo que // debe haber con LOPT=6 y BS=5.
// CALIBRADO!

        FW->produce[i][j]= IntegroComida(Comensales, comida, rate_food, dt, i, j, t); 


/* salida con peso promediado de adultos
 */
            auxx=pobla[i][j][92]+pobla[i][j][93]+pobla[i][j][94];
           if(auxx ==0)
              FW->sobra[i][j]=0.0;
           else
              FW->sobra[i][j]=(pobla[i][j][92]*FW->peso[i][j][92]+pobla[i][j][93]*FW->peso[i][j][93]+pobla[i][j][94]*FW->peso[i][j][94])/(1.0*auxx);

/* Salida con peso de pupas
	   FW->sobra[i][j]=0; auxx=0;
	   for (k=LarPobla+3;k < LarPobla+PupPobla+3;k++){
              FW->sobra[i][j] += FW->peso[i][j][k]*pobla[i][j][k];
	      auxx += pobla[i][j][k];
	     }
	   if(auxx==0)
		   FW->sobra[i][j]=0;
	    else
		    FW->sobra[i][j]= (FW->sobra[i][j])/auxx ;
*/
  
/* acá evoluciono el peso *//* Si engordaron --> las engordo. Pero si adelgazan, n las hago adelgazar y las dejo con el mismo peso anterior*/   
        for(e=3; e<(LarPobla+3);e++) { //solo importan las larvas, pues solo a ellas les evoluciona el peso
          FW->dBdt[i][j][e]= FW->beta[i][j]*log(log(FW->peso[i][j][e]/
	C.peso_huevos)/C.RateTop);

	    aux1= FW->peso[i][j][e];
            FW->peso[i][j][e]=pow(aux1,exp(-FW->beta[i][j]*dt));
        } // end e
      } // end of j
} // end of function comensales_y_alimentos

void check_food(FoodAndWeight FW)
{int i,j;
for(i=0;i<LIMfilas;i++)
	for(j=0;j<LIMcolumnas;j++)
		if(FW.produce[i][j] <0)
			fprintf(stdout,"BAD check %f %d %d\n",FW.produce[i][j],i,j);

return;

}
