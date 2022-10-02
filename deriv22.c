// Copyright by Victoria Romeo Aznar <vromeoaznar@gmail.com>
// Lucas Alonso <lucasalo28@gmail.com> and
// Hernán G Solari < hgsolari@gmail.com> under License
// GNU GENERAL PUBLIC LICENSE (see CPYING file)
//
/* Calculating the rates derivates */

void deri(int dimension, double l[], double dl[], double t, FoodAndWeight FW, int hls_o)
{   /* l is the lambda in a vector way, dl is the derivate in a vector way and t the time that here we don't used. */

double vm[LIMfilas][LIMcolumnas][POBLACIONES]; //Tnl[LIMfilas][LIMcolumnas][LarPobla+3], vm[LIMfilas][LIMcolumnas][POBLACIONES];
double hold1,hold2,hold3=1.;
int i,j,e,m;
float transfLP_comida, aux; //transLP_NL;
//long LarTot;

VM(vm, FW, hls_o);//VM(vm, Tnl, FW);
//trasfPP=trfLP(0,1)*PupPobla/LarPobla;
hold1= (0.585*C.transferpa + C.mpu);
hold2= (0.415 * C.transferpa);

for (i=0; i < LIMfilas; i++)
	for (j=0; j < LIMcolumnas; j++){        
        aux=1.;
        transfLP_comida = trfLP_comida(FW.produce[i][j]);                        
        m=(LIMcolumnas *i+ j)*EVENTOS;

        /* eggs */
        for(e=0;e<3;e++)
		    (dl+m)[2*e] = C.mh * vm[i][j][e];              /* egg mortality */
		    
		(dl+m)[1] = C.transferhh * vm[i][j][0]; 	       /* Transition from immature to mature egg */
		(dl+m)[3] = C.transferhml * vm[i][j][1]; 	       /* Transition of mature egg to larva */
		(dl+m)[5] = hold3 * C.transferhlll * vm[i][j][2];  /* Transition of rain egg to larva */
                
        /* Larvae */                  
        for(e=3;e<(LarPobla+3); e++) {
            (dl+m)[2*e] = m_lar(FW.produce[i][j])*vm[i][j][e]; /* mortalidad larvas*/
            if(e>(LarPobla-5)) /*larvae stage that depend on food */
                (dl+m)[2*e+1] = transfLP_comida*vm[i][j][e] ;//+ transLP_NL*Tnl[i][j][e];
            else
                (dl+m)[2*e+1] = C.transferlp_no_comida*vm[i][j][e]; /* larvae stage doesn't depend on food*/

        }
        /* Pupae */
        if(PupPobla>1)
            for(e = (LarPobla+3);e<(LarPobla+PupPobla+2); e++) {
                (dl+m)[2*e] = C.mpu * vm[i][j][e]; /*pupa mortality*/              
                (dl+m)[2*e+1] = C.transferpa * vm[i][j][e]; //c.transferpa * vm[i][j][e]/PupPobla;/* transition inter pupa*/
            }

		(dl+m)[2*(LarPobla+PupPobla+2)]= hold1 * vm[i][j][LarPobla+PupPobla+2];        /* pupae mortality,  discount for males and emergence mortality. */
		(dl+m)[2*(LarPobla+PupPobla+2)+1]= hold2 * vm[i][j][LarPobla+PupPobla+2];      /* Transition of  pupa to adult-type 1 */

        /* Adultos */ 
		(dl+m)[2*(LarPobla+PupPobla+3)]= C.madul * vm[i][j][LarPobla+PupPobla+3];                /* Adult mortality 1 */
		(dl+m)[2*(LarPobla+PupPobla+3)+1]= C.ov1 * vm[i][j][LarPobla+PupPobla+3];                /* Transition of  adult-type 1 to flying */
		(dl+m)[2*(LarPobla+PupPobla+4)]= C.madul * vm[i][j][LarPobla+PupPobla+4];                /* Flying mortality */
		(dl+m)[2*(LarPobla+PupPobla+4)+1]= C.postura * vm[i][j][LarPobla+PupPobla+4];            /* oviposition of  flying and transition to adult-2 */
		(dl+m)[2*(LarPobla+PupPobla+5)]= C.madul * vm[i][j][LarPobla+PupPobla+5];                /* Adult mortality 2 */
		(dl+m)[2*(LarPobla+PupPobla+5)+1]= C.ov2 * vm[i][j][LarPobla+PupPobla+5];                /* transition of adult-2 to flying */
		(dl+m)[2*(LarPobla+PupPobla+9)]= C.migr * vm[i][j][LarPobla+PupPobla+4];                 /* flying */

        /* eggs diapausa */
        for(e=(LarPobla+PupPobla+6);e<(LarPobla+PupPobla+9);e++)
		    (dl+m)[2*e] = C.mhd * vm[i][j][e];              /* egg mortality */
		(dl+m)[2*(LarPobla+PupPobla+6)+1] = C.transferhhd * vm[i][j][LarPobla+PupPobla+6];     /* Transition from immature to mature egg */
		(dl+m)[2*(LarPobla+PupPobla+7)+1] = C.transferhml * vm[i][j][LarPobla+PupPobla+7];    /* Transition of mature egg to larva */
		(dl+m)[2*(LarPobla+PupPobla+8)+1] = hold3 * C.transferhlll * vm[i][j][LarPobla+PupPobla+8]; /* Transition of rain egg to larva */

    };
}

