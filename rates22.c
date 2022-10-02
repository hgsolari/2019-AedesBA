// Copyright by Victoria Romeo Aznar <vromeoaznar@gmail.com>
// Lucas Alonso <lucasalo28@gmail.com> and
// Hernán G Solari < hgsolari@gmail.com> under License
// GNU GENERAL PUBLIC LICENSE (see CPYING file)
//
/* Developmental cycle. The values for 25ºC have been regionalized using
data for Cordoba strain published by Grech et al.,  
J Vector Ecol, 2010, 35, 277--285
Except for the development of larvae where we use our own data,
V Romeo Aznar et al.  Journal of Theoretical Biology 365 (2015) 311­324
*/

float thh(float grados) /* Coeficiente de huevo-huevo */
{float load;
load= (-0.0167105+0.03866*exp(((grados)+7.26)/16.47906));  //vale 0.2 a 21 grados
return load;
  /*without inhibition*/
};

float lp(float grados) /* Coeficiente de larva-pupa */
	/* ajuste segun DeMajo 2016, Romeo-Aznar 2018 y Bar-Zeev 1958
	 * con 89 pasos
	 */
{float load;
if(grados <= 10) load = 0.00314159712157291;
if (grados > 10) load= -0.0503046 + 0.00534462 * grados;
if (grados > 21) load=0.00884845*(grados-21)+0.0619324080398323;

load = load*(LarPobla*PupPobla); /* steps of development */
return load;
};

float mp(float grados) /* coeficiente mortalidad larva o pupa*/
{float load;
	/* Ajuste sobre datos de DeMajo-2016
	 * De Majo da la fracción de supervivientes, k, esta es el resultado
	 * de N=(89+7) pasos de desarrollo. k = (lp/(lp+mp))^N 
	 * se resuelve en mp=(k^(-1/N)-1)*lp ajustamos
	 * linealmente con (-0.00454284 +grados * 0.000430817) de lo cual
	 * resulta la función programada.
	 *
	 * Notamos que el rate aumenta con la temperatura aunque la mortalidad
	 * decrece ya que la mortalidad es el resultado de la competencia
	 * entre desarrollo y muerte */
if(grados > 10.)
	load=  (-0.00415 + grados* 0.0003935);
else
	load= 0.0000001; /* hay una ligera incosistencia, más que esperable,
			    entre los dos fits. Hay un límite cero sobre cero.
			    La arreglamos dentro del error de fiteo cambiando
			    la ordenda en el origen*/
return load;
};


/*float lp(float grados)  Coeficiente de larva-pupa */
/*{float load;
 load = (((0.005502-0.37266)/(1+exp(((grados)-25.189)/4.6456)))+0.37266);
return load;
}; */


/*float pa(float grados)  Coeficiente de pupa-adulto */
/*{float load
load= (0.06202*exp(((grados)-4.7655)/11.1062));
return load;
}; */

/*float mp(float grados) // coeficiente mortalidad larva o pupa
{float load;
 load= (0.06+0.97248*exp(-((grados)-5)/2.70346));  R.Aznar 
 load= (0.01+0.97248*exp(-((grados)-5)/2.70346));  Focks 
return load; 
}; */

float ovi1(float grados) /* first gonotrophic-cycle coefficient Focks */
{float load;
load=(0.03154*exp(((grados)-4.7511)/10.580));
return load;
};

float ovi2(float grados) /* gonotrophic-cycle coefficient Focks */
{float load;
load=((0.05427931*exp(((grados)-4.7511)/10.580)));
return load;
};

float rate_food(float grados) /*production rate of food for larvaes*/
/* APPLIED AND ENVIRONMENTAL MICROBIOLOGY, Apr. 1991, p. 1094-1101
modelo Ratkowsky-2  */
{float load, dif1, dif2,dif3;
 double LoptTot, DEGREE=1.;// tol=0.00000001;
 static float norma=0.,tmax=27.0; 
 dif1=TOPT-TMIN;
 dif2=grados-TOPT;
 dif3=grados-TMIN;
 if((grados > TMIN) & (grados < TOPT))
	{LoptTot=LOPT*C.BS;
/* codigo para hallar el máximo, normalizo ahora a la produccion de 27C
	if((tmax <= TMIN)|| (tmax > TOPT))
		tmax=0.75*dif1+TMIN;
	 while( abs(memo1-tmax) > tol)
	 	{memo1=tmax;
		tmax= TOPT - log((tmax-TMIN)*DEGREE+1) /DEGREE;
		}
*/
	if (norma == 0) /* remember the norma value */
		norma=(tmax-TMIN)*(1-exp(DEGREE*(tmax-TOPT)));
	load=dif3*(1-exp(DEGREE*dif2))/norma;
	load= LoptTot*load*load; 
	}
 else
    load = 0.0;
 return load;
};

/* sinusoidal variation of the temperature */

float varterm(int tp) {
  return TPROMEDIO+AMPLITUDT*cos(2*M_PI*tp/365.25+9.1305);   
}

void coeficientes(float Temp, FoodAndWeight  *FW) /* Produce the coefficients that varying per day */
{
	C.mpu=mp(Temp);
    C.mla=C.mpu;
 	C.transferhh=C.hh*thh(Temp);
 	C.transferhhd=C.hhd*thh(Temp);
    C.transferhml=C.eclo; // algun día rate de riego C.eclo*thh(Temp);
    C.transferhlll=C.eclo_ll*thh(Temp);
//    C.transferlp=lp(Temp);
    C.transferlp_no_comida=lp(Temp); /* tda la nmormalización en lp() */
//    C.transferlp_no_comida=lp(Temp)*81.27;
/* este número 139, ver el archivo
 *     Actualizacion_aedes06_10-4-2014.dat de Vico. */
	C.transferpa=C.transferlp_no_comida;
	C.ov1=ovi1(Temp);
	C.ov2=ovi2(Temp);
    FW->r_food=rate_food(Temp);
    FW->b=82.32; /* en el paper de JTB 2014, está en unidades de 1/hora, aqui en el programa las tasas van en 1/día. Lo pongo aquí porque en algun momento dependerá de la temperatura. Es un factor de adelgazamiento, dBdt= beta*(B/B_max)^alfa - b*B */
/*   C.RateTop=FW->b*6.629; 6.629 es el peso maximo del mosquito en mg, tomado como el peso en condiciones óptimas segun el paper (JTB 2014) en el grafico de peso vrs tiempo, con los datos de christophers*/

    if(FW->r_food < 0.) fprintf(stdout,"COEF %f\n",FW->r_food);
}

/*Calculating the "total" larvae mortality */
/* Si no adelgazan no aumenta la mortalidad. Pero cuando adelgazan, la mortalidad crece a los niveles de Barrera, es decir a la maxima mortalidad que teniamos en aedesBA06 */
float m_lar(float comida)
{
   float mlaTotal;
   mlaTotal=C.mla;

   if (comida < Lhungry)
	if (comida > Lstarving)
	mlaTotal+= C.mnl*(Lhungry-comida)/(Lhungry-Lstarving);
         else
           mlaTotal+= C.mnl;
	 
//Esta es el tipo de mortalidad 5   
   
   return mlaTotal;
}

float trfLP_comida(float comida){ // 29-12 Ojo que ahora cambio la comida ... fijarse si esto sigue asi. Donde dice comida voy a tener que usar comida/(V/Vf)
   float tlpTotal, aux, c_corte;
   c_corte=0.446; //¿de donde salió? Lhungry?
   if(comida/c_corte>=1.0)
      aux=1.0;
   else
      aux= comida/c_corte;
    
   tlpTotal=C.transferlp_no_comida*aux; /*fijarse que tlpTotal<=C.transferlp_no_comida pues beta2<=1*/  
    
   return tlpTotal;
}   



