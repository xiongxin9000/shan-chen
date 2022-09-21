#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath> // for fabs()
int const static n=9,mx=200,my=200; //number of latttice nodes
int c=1;//different cases
int const nfluids=1;//number of component
double w[n],ux[mx][my],uy[mx][my],x[mx],y[my];
double f[nfluids][n][mx][my],f2[nfluids][n][mx][my],feq[nfluids][n],rho[nfluids][mx][my];
const int cx[n]={0, 1, 0, -1, 0, 1, -1, -1, 1};
const int cy[n]={0, 0, 1, 0, -1, 1, 1, -1, -1};
double Fx[nfluids][mx][my]; // x-component of Shan-Chen force
double Fy[nfluids][mx][my]; // y-component of Shan-Chen force
double forcing[nfluids][n];
double press[nfluids][mx][my];
double rho_l=1.95,rho_g=0.15;
double radius=25;
double a=mx/10,b=my/5;
int i,j;
int dx=1,dy=1; //space and time step
// double const alpha=0.4;
// double omega=1.0/(3.*alpha+0.5);
double tauA=1.0;
double tauB=1.0;
int mstep=10; // The total number of time steps
int freq=1;
double rho0=1.0;//reference density 
const double gA = -4.7;//interaction strength of fluid A.
const double gB = 0;//interaction strength of fluid B.
const double gAB = 6;//interaction strength between fluid A and fluid B.

double g[nfluids][nfluids]; // interaction strengths
double tau[nfluids]; // relaxation times

void result(std::string filename,int time)
{

std::ofstream out(filename);
  out << "TITLE=\"LBM output\"" << std::endl;
  out << "VARIABLES = \"X\", \"Y\", \"Ux\",\"Uy\",\"density\",\"density1\",\"f0\",\"f1\",\"f2\",\"f3\",\"f4\",\"f5\",\"f6\",\"f7\",\"f8\"" << std::endl;
  out << "ZONE T = \"fluid\", I=" << mx << ", J=" << my << ", F=POINT" << std::endl;
  //out << "SOLUTIONTIME="<< time << std::endl;
  for (unsigned j = 0; j < my; ++j)
    for (unsigned i = 0; i < mx; ++i)
    {
        out << std::scientific << std::setprecision(5) << std::setw(15) << x[i];
        out << std::scientific << std::setprecision(5) << std::setw(15) << y[j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << ux[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << uy[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << rho[0][i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << rho[1][i][j];
      for(int k=0;k<9;k++)
      {
          out << std::scientific << std::setprecision(5) << std::setw(15) << f[0][k][i][j];
      }
      out << std::endl;
    }
  out.close();

}

void ComputeDesnity(double ftem[nfluids][n][mx][my])
{
    for(i=0;i<mx;i++)
    {
        for(j=0;j<my;j++)
        {
           double ssum[nfluids]={0};
            for (int s = 0; s < nfluids; s++)
            {
                for(int k=0;k<9;k++)
                {
                    ssum[s]=ssum[s]+ftem[s][k][i][j];
                }
           }
           for (int s = 0; s < nfluids; s++)
           {
                rho[s][i][j]=ssum[s];
           }
        }
    }
}
double ComputeTotalDensity(int i,int j)
{
    double dens;
  if (nfluids == 1) {
    dens = rho[0][i][j];
  }
  else if (nfluids == 2) {
    dens = rho[0][i][j] + rho[1][i][j];
  }
  return dens;
}
double psi(double dense)
{
    return rho0 * (1. - exp(-dense / rho0));
}
void ComputeSCForce() {
for (int i = 0; i < mx; i++) {
for (int j = 0; j < my; j++) {
      for (int s = 0; s < nfluids; s++) {
        Fx[s][i][j] = 0.;
        Fy[s][i][j] = 0.;
        for (int ss = 0; ss < nfluids; ss++) {
          double fxtemp = 0.;
          double fytemp = 0.;

          for (int k = 1; k < 9; k++) {
            const int i2 = (i + cx[k] + mx) % mx;
            const int j2 = (j + cy[k] + my) % my;
            double dense=rho[ss][i2][j2];
            double psinb = psi(dense);
            fxtemp += w[k] * cx[k] * psinb;
            fytemp += w[k] * cy[k] * psinb;
          }

          double psiloc = psi(rho[s][i][j]);
          fxtemp *= (-g[s][ss] * psiloc);
          fytemp *= (-g[s][ss] * psiloc);
          Fx[s][i][j] += fxtemp;
          Fy[s][i][j] += fytemp;
        }
      }
    }
  }

  return;
}
void ComputeVelocity(double ftem[nfluids][n][mx][my])//issue is here
{   
    for(int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
            ux[i][j] = 0.;
            uy[i][j] = 0.;
            for (int s = 0; s < nfluids; s++)
            {
            for(int k=0;k<9;k++)
            {
                    ux[i][j] += ftem[s][k][i][j]*cx[k];
                    uy[i][j] += ftem[s][k][i][j]*cy[k];
            }
            ux[i][j] += (0.5*Fx[s][i][j]);
            uy[i][j] += (0.5*Fy[s][i][j]);
            }
            double dens = ComputeTotalDensity(i,j);
            ux[i][j] /= dens;
            uy[i][j] /= dens;
        }
    }
}
void ComputeEquilibrium(int i,int j,int s)
{
    double dens = rho[s][i][j];
    double vx = ux[i][j];
    double vy = uy[i][j];
    double usq = vx * vx + vy * vy;
    // std::cout<<"dens= "<<dens<<std::endl;
    // std::cout<<"vx= "<<vx<<std::endl;
    // std::cout<<"vy= "<<vy<<std::endl;
    // std::cout<<"usq= "<<usq<<std::endl;
    feq[s][0] = w[0] * dens * (1. - 1.5 * usq);
    feq[s][1] = w[1] * dens * (1. + 3. * vx + 4.5 * vx * vx - 1.5 * usq);
    feq[s][2] = w[2] * dens * (1. + 3. * vy + 4.5 * vy * vy - 1.5 * usq);
    feq[s][3] = w[3] * dens * (1. - 3. * vx + 4.5 * vx * vx - 1.5 * usq);
    feq[s][4] = w[4] * dens * (1. - 3. * vy + 4.5 * vy * vy - 1.5 * usq);
    feq[s][5] = w[5] * dens * (1. + 3. * ( vx + vy) + 4.5 * ( vx + vy) * ( vx + vy) - 1.5 * usq);
    feq[s][6] = w[6] * dens * (1. + 3. * (-vx + vy) + 4.5 * (-vx + vy) * (-vx + vy) - 1.5 * usq);
    feq[s][7] = w[7] * dens * (1. + 3. * (-vx - vy) + 4.5 * ( vx + vy) * ( vx + vy) - 1.5 * usq);
    feq[s][8] = w[8] * dens * (1. + 3. * ( vx - vy) + 4.5 * ( vx - vy) * ( vx - vy) - 1.5 * usq);
}
void ComputeGuoForce(int i,int j,int k)
{
    for (int s = 0; s < nfluids; s++)
    {
        double temp=cx[k] * ux[i][j] + cy[k] * uy[i][j];
        const double omega = 1. / tau[s];
        forcing[s][k] = w[k] * (1. - 0.5 * omega) * 
        ((3. * (cx[k] - ux[i][j]) + 9. * cx[k] * temp) * Fx[s][i][j] + 
        (3. * (cy[k] - uy[i][j]) + 9. * cy[k] * temp) * Fy[s][i][j]);
    }
}
void Collision(double ftem[nfluids][n][mx][my],double ftem2[nfluids][n][mx][my])
{
    for (int s = 0; s < nfluids; s++)
    {
    const double omega = 1. / tau[s];
        for(int i=0;i<mx;i++)
        {
            for(int j=0;j<my;j++)
            {
                for (int k=0;k<9;k++)
                {
                    ComputeEquilibrium(i,j,s);
                    ComputeGuoForce(i,j,k);
                    int i2 = (i + cx[k] + mx) % mx;
                    int j2 = (j + cy[k] + my) % my;
                    ftem2[s][k][i2][j2] = ftem[s][k][i][j] * (1. - omega) + feq[s][k] * omega + forcing[s][k];
                }
            }
        }
    }
}
void Streaming(double ftem[nfluids][n][mx][my])
{
    double f_hlp[nfluids][9][mx][my];
    for (int s = 0; s < nfluids; s++)
    {
        for(int j=0;j<my;j++)
        {
            for (int i=0;i<mx;i++)
            {
                int y_n = (1+j)%my;
                int x_e = (1+i)%mx;
                int y_s = my - 1 - (my- j)%my;
                int x_w = mx - 1 - (mx- i)%mx;
                f_hlp[s][1][x_e][j] = ftem[s][1][i][j];
                f_hlp[s][2][i][y_n] = ftem[s][2][i][j];
                f_hlp[s][3][x_w][j] = ftem[s][3][i][j];
                f_hlp[s][4][i][y_s] = ftem[s][4][i][j];
                f_hlp[s][5][x_e][y_n] = ftem[s][5][i][j];
                f_hlp[s][6][x_w][y_n] = ftem[s][6][i][j];
                f_hlp[s][7][x_w][y_s] = ftem[s][7][i][j];
                f_hlp[s][8][x_e][y_s] = ftem[s][8][i][j];
            }
        }
    }
    for (int s = 0; s < nfluids; s++)
    {
        for(int j=0;j<my;j++)
        {
            for (int i=0;i<mx;i++)
            {
                for(int k=1;k<9;k++)
                {
                    ftem[s][k][i][j]=f_hlp[s][k][i][j];
                }
            }
        }
    }
}
void BoundaryCondition(double ftem[nfluids][n][mx][my])
{
    for (int s = 0; s < nfluids; s++)
    {
        for(int j=0;j<my;j++)
        {
            ftem[s][1][0][j]=ftem[s][3][0][j];//left
            ftem[s][5][0][j]=ftem[s][7][0][j];
            ftem[s][8][0][j]=ftem[s][6][0][j];

            ftem[s][3][mx-1][j]=ftem[s][1][mx-1][j];//right
            ftem[s][7][mx-1][j]=ftem[s][5][mx-1][j];
            ftem[s][6][mx-1][j]=ftem[s][8][mx-1][j];
        }

        for (int i=0;i<mx;i++)
        {
            ftem[s][2][i][0]=ftem[s][4][i][0];//bottom
            ftem[s][5][i][0]=ftem[s][7][i][0];
            ftem[s][6][i][0]=ftem[s][8][i][0];

            ftem[s][4][i][my-1]=ftem[s][2][i][my-1];//top
            ftem[s][7][i][my-1]=ftem[s][5][i][my-1];
            ftem[s][8][i][my-1]=ftem[s][6][i][my-1];
        }
    }
}
void initialize()
{
x[0] =0.0; 
y[0] =0.0;
//coordinate
for (i=1;i<mx;i++)
{
    x[i]=x[i-1]+dx;
}
for (j=1;j<my;j++)
{
    y[j]=y[j-1]+dy;
}
// Set viscosity
tau[0] = tauA;

if (nfluids == 2) {
tau[1] = tauB;
}

// Set interaction strength
g[0][0] = gA;

if (nfluids == 2) {
g[1][1] = gB;
g[0][1] = g[1][0] = gAB;
}    

std::cout<<"omegaA= "<<1/tauA<<std::endl;
std::cout<<"omegaB= "<<1/tauB<<std::endl;
/*-------weight factor ------*/
w[0]=4./9;
    for(i=1;i<5;i++)
    {
        w[i]=1./9;
    }
    for(i=5;i<9;i++)
    {
        w[i]=1./36;
    }
/*---------weight factor ----------*/
/*initial condition--------------*/
    for(int i=0;i<mx;i++)
    {
        for (int j=0;j<my;j++)
        {
            switch (c)
            {
            case 1:
            {
                if((x[i]-mx/2)*(x[i]-my/2)+(y[j]-mx/2)*(y[j]-my/2)<radius*radius)
            {
                rho[0][i][j]=rho_l;
            }
            else
            rho[0][i][j]=rho_g;
                break;
            }
            case 2:
            {
                if((x[i]-mx/2+radius-1)*(x[i]-mx/2+radius-1)+(y[j]-my/2)*(y[j]-my/2)<radius*radius
            ||(x[i]-mx/2-radius+1)*(x[i]-mx/2-radius+1)+(y[j]-my/2)*(y[j]-my/2)<radius*radius)
            {
                rho[0][i][j]=rho_l;
            }
            else
            rho[0][i][j]=rho_g;
                break;
            }
            case 3:
            {
                if((x[i]-mx/2)*(x[i]-my/2)/(a*a)+(y[j]-mx/2)*(y[j]-my/2)/(b*b)<1)
            {
                rho[0][i][j]=rho_l;
            }
            else
            rho[0][i][j]=rho_g;
                break;
            }
            }
            //rho[1][i][j]=0;
        }
    }
//initial condition of distribution function 
    for (int s = 0; s < nfluids; s++)
    {
        for(int i=0;i<mx;i++)
        {
            for(int j=0;j<my;j++)
            {
                for(int k=0;k<9;k++)
                {
                    ux[i][j]=0;
                    uy[i][j]=0;
                    ComputeEquilibrium(i,j,s);
                    f[s][k][i][j]=feq[s][k];
                    f2[s][k][i][j]=feq[s][k];
                }
            }
        }
    }
/*initial condition--------------*/
}
void VerifyLapalaceLaw()
{
//calculate pressure
for (int j = 0; j < my; j++) {
    for (int i = 0; i < mx; i++) {
        press[0][i][j] = 0.;
        press[0][i][j] += rho[0][i][j] / 3.;
        press[0][i][j] += (g[0][0] * psi(rho[0][i][j]) * psi(rho[0][i][j]) / 6.);
    }
}
//calculate deltapressure
double rho_gas = 0.;
double rho_liq = 0.;
double press_gas = 0.;
double press_liq = 0.;

// Average gas pressure in a square box with 4 lattice nodes
for (int i = 0; i < 4; i++) 
{
    for (int j = 0; j < 4; j++) {
        rho_gas += rho[0][i][j];
        press_gas += press[0][i][j];
    }
}

  // Average liquid pressure in a square box with 4 lattice nodes
  for (int i = mx / 2 - 2; i < mx / 2 + 2; i++) {
    for (int j = my / 2 - 2; j < my / 2 + 2; j++) {
      rho_liq += rho[0][i][j];
      press_liq += press[0][i][j];
    }
  }

  rho_gas /= 16.;
  rho_liq /= 16.;
  press_gas /= 16.;
  press_liq /= 16.;
  
  double delta_press = press_liq - press_gas;
  double rho_av = 0.5 * (rho_gas + rho_liq);
  const int y = my / 2;
//calculate radius
    double rad;
    for (int i = 0; i < mx / 2; i++) {
        if (rho[0][i][y] > rho_av) {
            const double drho = rho[i][y] - rho[i-1][y];
            const double dx = (rho_av - rho[0][i-1][y]) / drho;
            rad = (mx / 2. - i) + (1. - dx);
            break;
        }
    }
    //write data
    std::ofstream out;
    out.open("pressure-rad.dat",std::ios::app);
    
        out << std::scientific << std::setprecision(5) << std::setw(15) << delta_press;
        out << std::scientific << std::setprecision(5) << std::setw(15) << 1./rad;
        out << std::scientific << std::setprecision(5) << std::setw(15) << delta_press*rad;
        out << std::endl;
    out.close();
}
int main()
{
//main loop
initialize();
int time=0;
std::string filename = std::string("animation") +std::to_string(time)+std::string(".dat");
result(filename,time);
for (time=1;time<mstep+1;++time)
{
ComputeDesnity(f);
ComputeSCForce();//only novelty which contribute to equilibrium velocity and macroscopic velocity 
ComputeVelocity(f);
Collision(f,f2);//plus GuoForce 
std::swap(f,f2);
//Streaming(f);
BoundaryCondition(f);
std::string filename = std::string("animation") +std::to_string(time)+std::string(".dat");
if(time%freq==0) result(filename,time);
if(time%freq==0) {
    std::cout << "Iteration: " << time << " / " << mstep << " (" << 100.0*double(time)/double(mstep) << " %)" <<  std::endl;
}
VerifyLapalaceLaw();
}//main loop end
return 0;
}