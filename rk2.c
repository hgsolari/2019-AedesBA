// Copyright by Victoria Romeo Aznar <vromeoaznar@gmail.com>
// Lucas Alonso <lucasalo28@gmail.com> and
// Hernán G Solari < hgsolari@gmail.com> under License
// GNU GENERAL PUBLIC LICENSE (see CPYING file)
//
/* Runge Kutta integrator of order 2, error in dx^4 */

void rk2(void deri(int , double [], double [], double, FoodAndWeight, int), \
double h[], int n, double t, double dt, FoodAndWeight FW, int hls_o)
{

int i;
double k1[TOTALEVENTOS],k2[TOTALEVENTOS],h0[TOTALEVENTOS];
double dt2;

dt2=2.*dt/3.;

for (i = 0 ; i<n; i++)
	h0[i] = h[i];

deri(n,h0,k1,t,FW, hls_o);
for (i = 0 ; i<n; i++)
	h0[i] = h[i]+dt2*k1[i];
deri(n,h0,k2,t+dt2,FW, hls_o);

for (i =0; i<n ; i++)
	h[i]=h[i]+dt*((k1[i]+3*k2[i])/4.);

return;
}
