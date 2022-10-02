// Copyright by Victoria Romeo Aznar <vromeoaznar@gmail.com>
// Lucas Alonso <lucasalo28@gmail.com> and
// Hernán G Solari < hgsolari@gmail.com> under License
// GNU GENERAL PUBLIC LICENSE (see CPYING file)
//
/* Calculating the mean values ​​of the populations with border correction, include auxiliary functions. */

void	Gs(double Plus,double Minus,double *G0,double *G1,double *G2)
{/* calculo de las probabilidades de exeder por 0 1 y 2 un borde poblacional */
double z, y, yy;

if(Plus == 0){ *G0 = exp(-Minus);
	      *G1 = *G0*Minus;
	      *G2 = *G1*Minus/2.;
	      return;
	    }/* eventos que suman con tasa cero */

if(Minus == 0){*G0 = exp(-Plus);
		*G1=*G2=0.;
		return;
		} /* eventos que restan con tasa cero */

/* general case */
z=2.  * sqrt(Plus*Minus);
y= z/3.75;
yy= y*y;

if (y < 1) {
*G0=1.0+yy*(3.5156229+yy*(3.0899424+yy*(1.2067492
+yy*(0.2659732+yy*(0.360768e-1+yy*0.45813e-2)))));
*G0= *G0 * exp(-(Plus+Minus));
*G1= 2.*Minus*(0.5+yy*(0.87890594+yy*(0.51498869+yy*(0.15084934
+yy*(0.2658733e-1+yy*(0.301532e-2+yy*0.32411e-3))))));
*G1= *G1 * exp(-(Plus+Minus));
}
else {
y=1/y;
*G0=(exp(z-Plus-Minus)/sqrt(z))*(0.39894228+y*(0.1328592e-1
+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
+y*0.392377e-2))))))));
*G1=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
-y*0.420059e-2));
*G1=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
+y*(0.163801e-2+y*(-0.1031555e-1+y* *G1))));
*G1 =  sqrt(Minus/Plus)* *G1 * (exp(z-(Plus+Minus))/Plus);
}

*G2 = (- *G1 + Minus * *G0)/Plus;
return;
}

void fix(long n, double Plus, double Minus, double *Norm, double *correc,
double *G0, double *G1, double *G2)
/*calculo de la correccion por borde al valor medio devuelve Norma y correc*/
{
if((n < 0) ||  (n > 2)) {*Norm=1., *correc=0.; return;}
if ((Plus==0) && (Minus==0))
	{*correc=0.; *Norm=1; *G0=1.; *G1=*G2=0.; /* Evitar 0/0 */
	}
else
	{
Gs(Plus,Minus,G0,G1,G2);
switch (n){
     case 0: {*Norm = 1.-*G0-*G1-*G2;
	if(*Norm <= 0.)fprintf(stderr, "Case 0 Norma negativa  %g %g %g %g %g %g\n", *Norm,*G0, *G1,*G2, Plus, Minus);
	*correc = *G1+2* *G2;
        return; break;};
     case 1: {*Norm = 1- *G1- *G2;
	if(*Norm <= 0.)fprintf(stderr, "Case 1 Norma negativa %g\n", *Norm);
	*correc = *G2;
	return; break;};
     case 2: {*Norm = 1- *G2;
	if(*Norm <= 0.)fprintf(stderr, "Case 2 Norma negativa %g\n", *Norm);
	*correc = 0;
	   }
	  }
	}
}

/* This routine computes the mean values of the populations that are necessary for the rates derivate.*/

void VM(double vm[][LIMcolumnas][POBLACIONES], FoodAndWeight FW, int hls_o) //double Tnl[][LIMcolumnas][LarPobla+3], FoodAndWeight FW)
/* Lambdas, pobla y nevents son globales. Lbd los Lambdas es global */
/* Tnl es el termino no lineal  en larvas */
{
  int i,j,k, e,im,ip,jm,jp;
  long n;
  double LbdP, LbdM, G0, G1, G2, Norm, correc, correc2, hold;
  double aux, peso_alfa;

  for (i=0 ; i< LIMfilas; i++)
    for (j=0; j < LIMcolumnas; j++) {

/* fix retorna la correccion por borde de poblacion */

/* huevos no maduros*/
      n=pobla[i][j][0];
      LbdP=Lbd[i][j][2*(LarPobla+PupPobla+4)+1]*hls_o; //oviposicion del volador
      LbdM=Lbd[i][j][0]+Lbd[i][j][1];
/* El caso de los huevos el fix es diferente y lo hacemos directamente */
      if(n >= 3)
        vm[i][j][0]= n+FW.Fecun[i][j] * LbdP-LbdM;
	
      else {
        switch (n){
	  case 0: {G0=exp(-(LbdP+LbdM)); G1=G0*LbdM; G2=G1*LbdM/2.;
		Norm=1.-G0-G1-G2;
		correc= G1+2*G2;
		break;}
	  case 1: {G0=exp(-(LbdP+LbdM))*LbdM; G1=G0*LbdM/2.; G2=0.;
		Norm=1.-G1-G2;
		correc= G2;
		break;}
	  case 2: {G0=exp(-(LbdP+LbdM))*LbdM*LbdM/2.; G1=0.; G2=0.;
		Norm=1.-G2;
		correc= 0;;
		break;}
	  default: {fprintf(stderr,"error, poblacion negativa  i=%d j=%d poblacion %d\n", i,j,0); break;}
	};
	if((LbdP==0.) && (LbdM==0.)) Norm=1;
	vm[i][j][0]= n+(FW.Fecun[i][j] * LbdP-LbdM + correc)/Norm;
      };

/* Huevos Maduros */
      n=pobla[i][j][1];
      LbdP=Lbd[i][j][1]*(1.-C.eclo); /* descuento los que van a larva */
      LbdM=Lbd[i][j][2]+Lbd[i][j][3];
      fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
      vm[i][j][1]= n+(LbdP-LbdM + correc)/Norm;


/* Huevos lluvia */
      n=pobla[i][j][2];
      LbdP=0.0;
      LbdM=Lbd[i][j][4]+Lbd[i][j][5];
      fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
      vm[i][j][2]= n+(LbdP-LbdM + correc)/Norm;

/* algo que necesito para los terminos no lineales */
/*aux=Lbd[i][j][3]+Lbd[i][j][1]*C.eclo;
for(e=3;e<(LarPobla+3);e++) {
	LbdP=Lbd[i][j][2*e-1];
	LbdM=Lbd[i][j][2*e]+Lbd[i][j][2*e+1];
	aux += FW.peso_alfa[i][j][e]*(LbdP-LbdM)/C.Bmax_alfa;
}*/

/* larvae, first class of laevae */
      n=pobla[i][j][3];
      //LbdP=Lbd[i][j][3]+Lbd[i][j][5]+Lbd[i][j][1]*C.eclo;
	LbdP=Lbd[i][j][3]+Lbd[i][j][5]+Lbd[i][j][1]*C.eclo+Lbd[i][j][2*(LarPobla+PupPobla+7)+1]+Lbd[i][j][2*(LarPobla+PupPobla+8)+1]+Lbd[i][j][2*(LarPobla+PupPobla+6)+1]*C.eclo;;
/* Sumo la eclosion espontanea. Lbd[i][j][3] deberia ser 0 */
      LbdM=Lbd[i][j][6]+Lbd[i][j][7];
      fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
      hold = LbdP-LbdM;
      vm[i][j][3]= n+(LbdP-LbdM + correc)/Norm;
	
	/*switch (n){
		case 0: {correc2= -(2.* G1+6.* G2);break;}
		case 1: {correc2= -(2.*G2); break;}
		default: {correc2=0; break;}
		  }
        Tnl[i][j][3] = (1.*n +hold)*aux+ (FW.peso_alfa[i][j][3]*(LbdP + LbdM + Lbd[i][j][1]*C.eclo*(C.eclo-1))-Lbd[i][j][7]*FW.peso_alfa[i][j][4])/C.Bmax_alfa; */

/* Larvae up to pupae */ //lo pongo separado de las pupas por si llegan a tener algun termino no lineal
      if(LarPobla>1)
        for(e=4; e<(LarPobla+3);e++){
          n=pobla[i][j][e];
		  LbdP=Lbd[i][j][2*(e-1)+1];
	  	  LbdM=Lbd[i][j][2*e]+Lbd[i][j][2*e+1];
	  	  fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
          hold = LbdP-LbdM;
	  	  vm[i][j][e]= n+(LbdP-LbdM + correc)/Norm;
           /* switch (n){
		case 0: {correc2= -(2.* G1+6.* G2);break;}
		case 1: {correc2= -(2.*G2); break;}
		default: {correc2=0; break;}
		  }
            Tnl[i][j][e]=(1.*n+hold)*aux+(FW.peso_alfa[i][j][e]*(LbdP+LbdM)-FW.peso_alfa[i][j][e-1]*Lbd[i][j][2*e-1]-FW.peso_alfa[i][j][e+1]*Lbd[i][j][2*e+1])/C.Bmax_alfa;
            if((e==LarPobla+2))
               Tnl[i][j][e] += FW.peso_alfa[i][j][e+1]*Lbd[i][j][2*e+1]/C.Bmax_alfa;*/ 
        }
        
/*Pupae */
      for(e=(LarPobla+3); e<(LarPobla+PupPobla+3); e++) {
	n=pobla[i][j][e];
	LbdP=Lbd[i][j][2*(e-1)+1];
	LbdM=Lbd[i][j][2*e]+Lbd[i][j][2*e+1];
	fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	vm[i][j][e]= n+(LbdP-LbdM + correc)/Norm;
      }

/* Adultos 1 */
      n=pobla[i][j][LarPobla+PupPobla+3];
      LbdP=Lbd[i][j][2*(LarPobla+PupPobla+2)+1];
      LbdM=Lbd[i][j][2*(LarPobla+PupPobla+3)]+Lbd[i][j][2*(LarPobla+PupPobla+3)+1];
      fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
      vm[i][j][LarPobla+PupPobla+3]= n+(LbdP-LbdM + correc)/Norm;

/* Voladoras */
/* arreglo de indices para los vuelos */
      im=i-1; jm=j-1; ip=i+1;jp=j+1;
      if( (im < 0) ) im=0;
      if( (jm < 0) ) jm=0;
      if( (jp >= LIMcolumnas) ) jp=LIMcolumnas-1;
      if( (ip >= LIMfilas) ) ip=LIMfilas-1;
	n=pobla[i][j][LarPobla+PupPobla+4];
	LbdP=Lbd[i][j][2*(LarPobla+PupPobla+3)+1]+Lbd[i][j][2*(LarPobla+PupPobla+5)+1]+ /* mas los que llegan */
	(Lbd[ip][j][2*(LarPobla+PupPobla+6)] + Lbd[i][jm][2*(LarPobla+PupPobla+6)]+ Lbd[im][j][2*(LarPobla+PupPobla+6)] + Lbd[i][jp][2*(LarPobla+PupPobla+6)]) *C.pvuel[1]+
	(Lbd[ip][jp][2*(LarPobla+PupPobla+6)]+Lbd[im][jm][2*(LarPobla+PupPobla+6)]+ Lbd[ip][jm][2*(LarPobla+PupPobla+6)]+ Lbd[im][jp][2*(LarPobla+PupPobla+6)])*C.pvuel[0];
	LbdM=Lbd[i][j][2*(LarPobla+PupPobla+4)]+Lbd[i][j][2*(LarPobla+PupPobla+4)+1]+ Lbd[i][j][2*(LarPobla+PupPobla+6)];
		/* y los que se van estan solo en 12 */
	fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	vm[i][j][LarPobla+PupPobla+4]= n+(LbdP-LbdM + correc)/Norm;

/* Adultos 2 */
	n=pobla[i][j][LarPobla+PupPobla+5];
	LbdP=Lbd[i][j][2*(LarPobla+PupPobla+4)+1];
	LbdM=Lbd[i][j][2*(LarPobla+PupPobla+5)]+Lbd[i][j][2*(LarPobla+PupPobla+5)+1];
	fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
	vm[i][j][LarPobla+PupPobla+5]= n+(LbdP-LbdM + correc)/Norm;

/* DIAPAUSA */
/* huevos no maduros*/
      n=pobla[i][j][LarPobla+PupPobla+6];
      LbdP=Lbd[i][j][2*(LarPobla+PupPobla+4)+1]*(1-hls_o); //oviposicion del volador
      LbdM=Lbd[i][j][2*(LarPobla+PupPobla+6)]+Lbd[i][j][2*(LarPobla+PupPobla+6)+1];
/* El caso de los huevos el fix es diferente y lo hacemos directamente */
      if(n >= 3)
        vm[i][j][(LarPobla+PupPobla+6)]= n+FW.Fecun[i][j] * LbdP-LbdM;
	
      else {
        switch (n){
	  case 0: {G0=exp(-(LbdP+LbdM)); G1=G0*LbdM; G2=G1*LbdM/2.;
		Norm=1.-G0-G1-G2;
		correc= G1+2*G2;
		break;}
	  case 1: {G0=exp(-(LbdP+LbdM))*LbdM; G1=G0*LbdM/2.; G2=0.;
		Norm=1.-G1-G2;
		correc= G2;
		break;}
	  case 2: {G0=exp(-(LbdP+LbdM))*LbdM*LbdM/2.; G1=0.; G2=0.;
		Norm=1.-G2;
		correc= 0;;
		break;}
	  default: {fprintf(stderr,"error, poblacion negativa  i=%d j=%d poblacion %d\n", i,j,0); break;}
	};
	if((LbdP==0.) && (LbdM==0.)) Norm=1;
	vm[i][j][(LarPobla+PupPobla+6)]= n+(FW.Fecun[i][j] * LbdP-LbdM + correc)/Norm;
      };

/* Huevos Maduros */
      n=pobla[i][j][(LarPobla+PupPobla+7)];
      LbdP=Lbd[i][j][2*(LarPobla+PupPobla+6)+1]*(1.-C.eclo); /* descuento los que van a larva */
      LbdM=Lbd[i][j][2*(LarPobla+PupPobla+7)]+Lbd[i][j][2*(LarPobla+PupPobla+7)+1];
      fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
      vm[i][j][LarPobla+PupPobla+7]= n+(LbdP-LbdM + correc)/Norm;

/* Huevos lluvia */
      n=pobla[i][j][LarPobla+PupPobla+8];
      LbdP=0.0;
      LbdM=Lbd[i][j][2*(LarPobla+PupPobla+8)]+Lbd[i][j][2*(LarPobla+PupPobla+8)+1];
      fix(n, LbdP, LbdM, &Norm, &correc, &G0, &G1, &G2);
      vm[i][j][LarPobla+PupPobla+8]= n+(LbdP-LbdM + correc)/Norm;

	for (k=0; k < POBLACIONES; k++)
		if(vm[i][j][k] < 0.){
		fprintf(stderr, "Valor medio negativo (%d,%d,%d) %e\n",i,j,k, vm[i][j][k]);
		vm[i][j][k]=0. ;
		}

    };//fin del j
}//fin de la funcion
