#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define AEDESINPUT "aedes19.inp" /* input file for the code */
#define TPROMEDIO 18.0
#define AMPLITUDT 6.7389
#define LIMfilas 10 /*total number of rows considered */ //5x5 
#define LIMcolumnas 6  /*total number of columns considered */
#define PASO 12 //
#define REPITE 20 /* total repetitions  */ //20
#define UmbralDiaPausa 11 // umbral para la diapausa
#define UmbralLluvia 7.5 /* rainfall threshold for the eggs eclosion*/
#define TsaveCI -800 /* -670 Transition time (time before the simulation). -2863 PARA MDQ, AZUL, TAND */
/* larvae's parameters  alimentation */
#define TOPT 40.0  /* Optimum temperature for the food production */ 
#define TMIN 11.0  /* Minimum temperature for the food production */ 
#define DECAE 0.3  /* food decay rate */
#define LOPT  6  /* Optimum number of larvae per breeding site */
#define VOL   6 /* Volumen of breading site computed as the number of
max-size larvae needed to filter it in one day */
#define GILLETTM 0.8 /* Max probabilidad of eclosion in the presence of food */
#define GILLETTm 0.05 /* Probability of eclosion in the abscence of food */
#define Lhungry 1./70.921986 /* Food level to start to increase mortality */
#define Lstarving 0.0 /* Food level to complete the increase in mortality */
#define SCAV 0.01 /* Dead larvae return to food with this coefficient */

#define pm 0.4 /* peso mínimo del mosquito fértil */


typedef struct Common
	{float mh,mhd,mla,mpu,madul,eclo,eclo_ll,hh,hhd;
	 float transferhh, transferhml,transferhlll, transferhhd, transferlp, transferpa, ov1, ov2, transferlp_no_comida;
	 int fecundidad;
	 float postura;
	 float migr;
	 float mnl;
	 double pvuel[8];
	 float BS;
         /*  son constantes de food and weight */
         float alfa, RateTop, peso_huevos, Bmax, Vf, V; /*Bmax^alfa, rateTop: rate de comida a temperatura optima*/
         float a[8];
         float m[8];
         float n[8];
	} common;  /* define the common type for the Common structure that have the pointers to variable that are used in internal subroutines like the calculation of derivated */ /* Some of these values are in CONSTANTES */

/*#define JULIO 184*/          /* off set del 1 de julio a 1 de enero */
#define PupPobla  22      /* total number of pupae populations. This is for the case femmale and male together */ //61
#define LarPobla 67      //90 /* total number of larvae populations, there are 7 that depends on food and the rest that don't. This is for the case femmale and male together */
#define EVENTOS 204        /*total number of events per cell. Two for each population and one for the flight */ //315
#define POBLACIONES 98   /* total number of populations per cell. There are,three populations of eggs, three for adults and the rest are larva-pupa ones (90 for larvae and 65 for pupae) */ //157
#define TOTALPOBLA LIMcolumnas * LIMfilas * POBLACIONES /* total populations */
#define TOTALEVENTOS   LIMcolumnas * LIMfilas * EVENTOS /* total events */

typedef struct foodandweight /*variables que no sean constates en el tiempo*/
	{float produce[LIMfilas][LIMcolumnas]; /* disponible= sobrante + producido. Takes values in producimos function at comidas.c file */
         float sobra[LIMfilas][LIMcolumnas]; /* sobrante por celda. Its place it is equal that produce */
         float r_food; /* production rate, depends on the temperature. You can found it in rates.c */
	 /*float decae;  decay factor = exp(-DECAE/PASOS) */
         /*float Bmax_alfa  It is defined in CONSTANTES*/
         /* Now for the weight 
            First the parameters of the weight equation; dB/dt = a*B^alpha - b*B */
         /*float alpha;  It is defined in CONSTANTES*/
         float b; /* It is defined in rates.c because in principle depends on the temperature*/
         float beta[LIMfilas][LIMcolumnas]; /* This one is more complex, depend on the food and larvae poblation, it is definded in comidas.c in LaBeta function */
         /* and here the wheight and its derivated */
         float peso[LIMfilas][LIMcolumnas][POBLACIONES]; // fijarse si esta bien la dimension de esto porque tambien tengo que tener en cuenta los huevos
         float peso_alfa[LIMfilas][LIMcolumnas][LarPobla+3];// SOLO LARVAS
         float dBdt[LIMfilas][LIMcolumnas][LarPobla+3]; // SOLO LARVAS
         int Fecun[LIMfilas][LIMcolumnas];
}FoodAndWeight; /* define the FoodAndWeight type for the foodandweight structure, then I point them with a pointer and take them to deriv. In this way, I avoid to declare more globals variables and change the structure of rk2 too*/


/* prototypes functions */
float ran1(long *);         /* generates a random number with uniform distribution */
float gammln(float );       /* Ln gamma distribution, it is necessary for Poisson distribution*/
int poidev(float , long *); /* Poisson distribution*/
float thh(float );          /* Eclosion coefficient */
float lp(float );           /* Pupation coefficient */
float pa(float );           /* Coefficient of emerging adult  */
float ovi1(float );         /* Coefficient of gonotrophic cycle, type 1 */
float ovi2(float );         /* Coefficient of gonotrophic cycle, type 2 */
float mp(float );           /* Mortality coefficient for larvae and pupae */
float gasdev(long *);       /* Generation of numbers with normal N(0,1) distribution */
float varterm(int);         /* Function to define the temperature*/
float LaBeta(float);
int fecundidad(float);

/* simple multinomial */
unsigned long BINV(double, long ,long *);
unsigned long BTRD(double ,long ,long *);
unsigned long binomial(double ,long ,  long *);
void multinomialS(int , int , int ,double [],int [],long *);
void multinomialR(int ,int ,int ,double [],int [],long *);


typedef struct Ciudades
 { float a_lp;
   float a_pa;
   float a_mp;
   float transferlp_26;
   float a_ml;
   float a_ovi1;
   float a_ovi2;
 } ciudades; /* These values are in CONSTANTES */


