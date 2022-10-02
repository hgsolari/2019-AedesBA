// Copyright by Victoria Romeo Aznar <vromeoaznar@gmail.com>
// Lucas Alonso <lucasalo28@gmail.com> and
// Hernán G Solari < hgsolari@gmail.com> under License
// GNU GENERAL PUBLIC LICENSE (see CPYING file)
//
/* Codigo Aedes para el mosquito Aedes aegypti
Esta version incorpora la comida como variable dinámica y los rates de
desarrollo y mortalidad denso-dependientes medidos por Macia 2009.
(ISSN 0373-5680 Rev. Soc. Entomol. Argent. 68 (1-2): 107-114, 2009)
El rate de produccion  de comida es regulado por la ley de Ratkowsky-2
(JOURNAL OF BACTERiOLOGY, Jan. 1982, p. 1-5 Vol. 149, No. 1)
La eclosion de los huevos puede deberse a tres factores:
Lluvia, espontanea (inmediata a maduracion)
*/ 

/* POPULATIONS REFERENCES
   0 ---------------------------------------------- huevos no maduros 
   1 ---------------------------------------------- huevos maduros 
   2 ---------------------------------------------- huevos lluvia 
   3 a LarPobla+2 --------------------------------- larvas
   LarPobla+3 a LarPobla+PupPobla+2 --------------- pupas (from 3 to N_LP_pobl+3)
   LarPobla+PupPobla+3 ---------------------------- Adultos 1 154
   LarPobla+PupPobla+4 ---------------------------- Voladoras
   LarPobla+PupPobla+5 ---------------------------- Adultos 2 

   LarPobla+PupPobla+6----------------------------- huevos no maduros con diapausa
   LarPobla+PupPobla+7----------------------------- huevos maduros con diapausa
   LarPobla+PupPobla+8----------------------------- huevos lluvia con diapausa  */

/* EVENTS REFERENCES
   0 ---------------------------------------------- mortalidad de huevos no maduros
   1 ---------------------------------------------- pasaje de huevo no maduros a huevos maduros
   2 ---------------------------------------------- mortalidad de huevos maduros
   3 ---------------------------------------------- pasaje de huevos inmaduros a larva espontanea.
   4 ---------------------------------------------- mortalidad de huevos lluvia
   5 ---------------------------------------------- pasaje de huevos lluvia a larva
   2*i con 3<i<LarPobla+2 ------------------------- mortalidad de larvas --> pobla[i]
   2*i+1 con 3<i<LarPobla+2 ----------------------- pasaje de larva a larva ---> de pobla[i] a pobla[i+1]
   2*(LarPobla+2)+1-------------------------------- pasaje de larva a pupa
   2*i con (LarPobla+3)<i<(LarPobla+PupPobla+2) --- mortalidad de pupas --> pobla[i]
   2*i+1 con (LarPobla+2)<i<(LarPobla+PupPobla+2) - pasaje de pupa a pupa ---> de pobla[i] a pobla[i+1]
   2*(LarPobla+PupPobla+2)+1----------------------- pasaje de pupa a adultas tipo 1
   2*(LarPobla+PupPobla+3) ------------------------ mortalidad de adultos tipo 1
   2*(LarPobla+PupPobla+3)+1----------------------- pasaje de adultos tipo 1 a volador
   2*(LarPobla+PupPobla+4)------------------------- mortalidad de volador
   2*(LarPobla+PupPobla+4)+1----------------------- oviposicion por volador y pasaje a adultos 2
   2*(LarPobla+PupPobla+5)------------------------- mortalidad de adulto 2
   2*(LarPobla+PupPobla+5)+1----------------------- pasaje de adulto 2 a volador

   2*(LarPobla+PupPobla+6)------------------------- mortalidad de huevos no maduros con diapausa
   2*(LarPobla+PupPobla+6)+1----------------------- pasaje de huevo no maduros con diapausa a huevos maduros con diapausa 
   2*(LarPobla+PupPobla+7) ------------------------ mortalidad de huevos maduros con diapausa 
   2*(LarPobla+PupPobla+7)+1----------------------- pasaje de huevos INmaduros con diapausa a larva espontanea.
   2*(LarPobla+PupPobla+8)------------------------- mortalidad de huevos lluvia con diapausa
   2*(LarPobla+PupPobla+8)+1----------------------- pasaje de huevos lluvia con diapausa a larva

   2*(LarPobla+PupPobla+9)+j------------------------ vuelo con 0<= j <= 7

	Distribuido con C.pvuel[]
		0 vuelo a -1 -1
		1 vuelo a -1 0
		2 vuelo a -1 1
		3 vuelo a  0 1
		4 vuelo a  1 1
		5 vuelo a  1 0
		6 vuelo a  1 -1
		7 vuelo a  0 -1
*/

#include "aedes22.h"
#include "algoritmos.c" /* generacion de numeros aleatorios */

/* defino las variables globales */
common C;
ciudades Cd;

long pobla[LIMfilas][LIMcolumnas][POBLACIONES];
long poblaA[LIMfilas][LIMcolumnas][POBLACIONES];
long neventos[LIMfilas][LIMcolumnas][EVENTOS+7];
/* El vuelo se distribuye en 8 eventos */
double Lbd[LIMfilas][LIMcolumnas][EVENTOS];

#include "f-auxiliares22.c" /*calculo de coeficientes y rutinas de entrada */
#include "rates22.c"//"rates06_ciudades.c", se utiliza para utilizar parametros de rates segun las ciudades
#include "d-alimentacion22.c" /*calculo de la produccion de alimento y lo que comen las larvas*/
#include "vm22.c" /* rutina de calculo de valores medios de poblaciones */
#include "rk2.c"
#include "deriv22.c"
#include "reparto22_b.c"


/* Cuerpo principal del programa */

int main()
{ 

   int DP=1;  /*1 para que pongan huevos con Diapausa, 0 para que no*/

   long idum= 3480;  /*Semilla para la generacion de numeros pseudoaleatorios*/

   double t=0.0, *pointerLbd;

   FoodAndWeight FW;



   char *poblaname[POBLACIONES+2], mastername[40], *eventoname[EVENTOS+7];
   char poblaN[POBLACIONES+2][48], eventoN[EVENTOS+7][48];
   FILE *FilePobla[POBLACIONES+2], *FileEventos[EVENTOS+7];
    /*lectura de datos e impresion de resultados*/

   int i,l,o,Toffset;
   long *pointereventos;
    /* tiempo del transitorio y de las simulacion (en días)  a ser leidos.
    fecha de inicio queda fija en el 1 de julio */
   int Transitorio=0, Tsimulado=0, TransitorioLeido=-1, repite;
   float Dt=1./PASO; /* paso de tiempo en fraccion dias */

   int Plist[POBLACIONES], Elist[EVENTOS+7];
    /* listado de eventos y poblaciones para imprimir -Print- */

   float temperaturavec[15000],lluviavec[15000],hlsvec[15000];                              
    /*listado de temperaturas y lluvias*/

   float Temp=0.0;   /*Temperatura del día*/                                                  
   float lluvia=0.0; /*Lluvia del día*/      
   float hls=0.0; /*horas d eluz solar del día*/ 
   int hls_o; /*función escalón para poner huevos originales o con diapausa*/ 
                                             
    /* coeficiente ambiental */

   int tiempo_inicial=0; /* tiempo en dias al comenzar el cálculo */

   int k,j,e;

   pointerLbd=Lbd[0][0]; /* pointer para recorrer los lambdas */
   pointereventos=neventos[0][0]; /* pointer para recorrer eventos */

    /* asigno los pointers de los nombres (POCO ELEGANTE) */
   for (i=0;i<POBLACIONES+2; i++)
      poblaname[i]=poblaN[i];
   for (i=0;i<EVENTOS+7; i++)
      eventoname[i]=eventoN[i];

    /*lee los datos de entrada y abre los files, todos los argumentos son salida */
   if((i=getdata(FilePobla, FileEventos, Plist, Elist, mastername, poblaname, eventoname, &Transitorio, &Tsimulado, &idum, temperaturavec, lluviavec, hlsvec))<0){
      fprintf(stderr,"Error al leer %s es %d\n", AEDESINPUT,i);
   };

   #include "CONSTANTES22"
   fprintf(stderr,"Transitorio %d  Simulado %d semilla %ld mnl %f\n", Transitorio, Tsimulado, idum, C.mnl);

   TransitorioLeido=Transitorio;

   /* inicializo las constantes de evolucion que no dependen de la Temperatura */
 
/* REPETICION para la estadistica */

for (repite=1; repite <= REPITE; repite++) {

    fprintf(stderr, "Repetición %d\n",repite); 

    pobla_iniciales(&Transitorio, &tiempo_inicial, &FW);     /* inicializo las poblaciones */ 


    for(i=0; i< (EVENTOS+7) * LIMfilas * LIMcolumnas; i++)
	    *(pointereventos+i)=0; /* pongo en cero los eventos del dia */

    /* inicio el transitorio */
    for (l=tiempo_inicial; l< Transitorio+tiempo_inicial ; l++) {  
	    Temp=temperaturavec[l-TsaveCI];  /*Marce se fija la temperatura en el vector*/
        lluvia=lluviavec[l-TsaveCI];     /*Marce se fija la lluvia en el vector*/ 
        hls=hlsvec[l-TsaveCI];

        if(hls > UmbralDiaPausa)
            hls_o=1;
        else
            hls_o=1-DP;
 
   	coeficientes(Temp,&FW);

        if((l == 0) && (tiempo_inicial != 0)) save_transitorio(pobla, FW);

         /* corre un día */
        /* Hernan pone a cero los eventos */
        for(i=0; i< (EVENTOS+7) * LIMfilas * LIMcolumnas; i++)
            *(pointereventos+i)=0; /* pongo en cero los eventos del dia */
        for(o=0; o< PASO; o++){ 
            t=l+(o+1.)*Dt; 
	        for(i=0; i< TOTALEVENTOS; i++)
	            *(pointerLbd+i)=0; /* pongo a cero el lambda */         

	        rk2(deri, Lbd[0][0], TOTALEVENTOS, t, Dt, FW, hls_o); 
            comensales_y_alimento(&FW,Dt,t);

            updatepobla(&idum,&FW,hls_o); /* aqui tambien promedio los pesos */
         
            if((lluvia > UmbralLluvia) && (PASO/2 == o)){
                /* AQUI camnio la lluvia */
                llovio(&idum,FW);
	        };     
       
	    }; //fin for (o) 
    }; //fin for (l)

    /* aseguro que el valor de la Temp y Lluvia impresos sean correctos */
    if(Transitorio+tiempo_inicial-TsaveCI > 0) {
        Temp=temperaturavec[Transitorio+tiempo_inicial-TsaveCI-1];  /*miro la temperatura en el vector al final del transitorio*/
        lluvia=lluviavec[Transitorio+tiempo_inicial-TsaveCI-1];     /*miro la lluvia en el vector al final del transitorio*/    
        hls=hlsvec[Transitorio+tiempo_inicial-TsaveCI-1]; 
    }
    else {
        Temp=temperaturavec[0];
        lluvia=lluviavec[0];
        hls=hlsvec[0];
    }
    /* calculo la temperatura para los datos iniciales */

    for (i=0; (i<POBLACIONES) && (Plist[i] >=0 ); i++)
         fprintf(FilePobla[i],"%d\n", repite); 
    for (i=0; (i<EVENTOS) && (Elist[i] >=0 ); i++)
  	     fprintf(FileEventos[i],"%d\n",repite);
    /* Pongo un cartel en los files */


  /* SIMULACION */

    if((l == 0) && (tiempo_inicial != 0) ) save_transitorio(pobla,FW);  /*guardo las cond iniciales*/ /* FALTA EDITAR */

    /* Hernan corrige la salida de Y */
    Toffset= l - (TransitorioLeido+TsaveCI);  /* para poner el día de simulación en files */
    savedata(FilePobla,FileEventos,Plist,Elist, t=(float)Toffset, Temp, lluvia,FW); /* guardo la salida */ /* FALTA EDITAR */

    /* loop de tiempo despues de guardar la CI */
    for (l=tiempo_inicial+Transitorio; l< Tsimulado+Transitorio+tiempo_inicial-1 ; l++) {
        Temp=temperaturavec[l-TsaveCI];
        lluvia=lluviavec[l-TsaveCI]; /* lluvia del dia*/
        hls=hlsvec[l-TsaveCI];       

        if(hls > UmbralDiaPausa)
            hls_o=1;
        else
            hls_o=1-DP;

	coeficientes(Temp,&FW);

        for(i=0; i< (EVENTOS+7) * LIMfilas * LIMcolumnas; i++)
	        *(pointereventos+i)=0; /* pongo en cero los eventos del dia */
         
	    for(o=0; o< PASO; o++) { 
 	        t=l+o*Dt;
	        for(i=0; i< TOTALEVENTOS; i++)
	            *(pointerLbd+i)=0; /* pongo a cero el lambda */
            
	        rk2(deri, Lbd[0][0], TOTALEVENTOS, t, Dt, FW, hls_o);
 	        comensales_y_alimento(&FW,Dt,t);

	        updatepobla(&idum,&FW,hls_o);  
           
            /* el tiempo va en dias de la parte visible de la corrida. Controla en cada paso, de lo contrario no hay acumulacion de comida */

            if((lluvia > UmbralLluvia) && (PASO/2 == o)){
                llovio(&idum,FW);
	        };
            
        }; /*fin del for del indice o*/
      
        lluvia=lluviavec[l-TsaveCI];     /*mira la lluvia en el vector*/  
        Toffset= l - (TransitorioLeido+TsaveCI) + 1;

        savedata(FilePobla,FileEventos,Plist,Elist, t=(float)Toffset, Temp, lluvia,FW); /*guardo las pob, eventos, temp y lluvia*/ /* FALTA EDITAR */
            
    };  /*fin del loop del tiempo, indice l*/      

    Transitorio= TransitorioLeido; /* Transitorio va a ser decrementado en TsaveCI al leer el archivo de CI generado */

} /* FIN DEL LOOP DE REPETICION */

nuevasemilla(-idum); /* guarde la semilla el menos para ran1 de nrc con + para ran33 */

exit(0);
} /* CIERRO MAIN */
