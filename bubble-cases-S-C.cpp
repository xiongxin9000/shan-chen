#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath> // for fabs()
int const static n=9,mx=200,my=200; //number of latttice nodes
int c=2;//different cases
double f[n][mx][my],f2[n][mx][my],feq[n],rho[mx][my],w[n],ux[mx][my],uy[mx][my],x[mx],y[my];
const int cx[n]={0, 1, 0, -1, 0, 1, -1, -1, 1};
const int cy[n]={0, 0, 1, 0, -1, 1, 1, -1, -1};
double Fx[mx][my]; // x-component of Shan-Chen force
double Fy[mx][my]; // y-component of Shan-Chen force
double forcing[n];
double rho_l=0.15,rho_g=1.95;
double radius=mx/8;
double a=mx/10,b=my/5;
int i,j;
int dx=1,dy=1; //space and time step
// double const alpha=0.4;
// double omega=1.0/(3.*alpha+0.5);
double omega=1.0;
int mstep=50000; // The total number of time steps
int freq=1000;
double rho0=1.0;//reference density 
const double g = -4.7;//interaction strength between particles.
void result(std::string filename,int time)
{

std::ofstream out(filename);
  out << "TITLE=\"LBM output\"" << std::endl;
  out << "VARIABLES = \"X\", \"Y\", \"Ux\",\"Uy\",\"density\",\"f0\",\"f1\",\"f2\",\"f3\",\"f4\",\"f5\",\"f6\",\"f7\",\"f8\",\"feq0\",\"feq1\",\"feq2\",\"feq3\",\"feq4\",\"feq5\",\"feq6\",\"feq7\",\"feq8\"," << std::endl;
  out << "ZONE T = \"fluid\", I=" << mx << ", J=" << my << ", F=POINT" << std::endl;
  //out << "SOLUTIONTIME="<< time << std::endl;
  for (unsigned j = 0; j < my; ++j)
    for (unsigned i = 0; i < mx; ++i)
    {
        out << std::scientific << std::setprecision(5) << std::setw(15) << x[i];
        out << std::scientific << std::setprecision(5) << std::setw(15) << y[j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << ux[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << uy[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << rho[i][j];
      for(int k=0;k<9;k++)
      {
          out << std::scientific << std::setprecision(5) << std::setw(15) << f[k][i][j];
      }
      for(int k=0;k<9;k++)
      {
          out << std::scientific << std::setprecision(5) << std::setw(15) << feq[k];
      }
      out << std::endl;
    }
  out.close();

}

void ComputeDesnity(double ftem[n][mx][my])
{
    for(i=0;i<mx;i++)
    {
        for(j=0;j<my;j++)
        {
           double ssum=0;
           for(int k=0;k<9;k++)
           {
               ssum=ssum+ftem[k][i][j];
           }
           rho[i][j]=ssum;
        }
    }
}
double ComputeTotalDensity(int i,int j)
{
    double dense = rho[i][j];
    return dense;
}
double psi(double dense)
{
    return rho0 * (1. - exp(-dense / rho0));
}
void ComputeSCForce()
{
    for(int i=0;i<mx;i++)
    {
        for (int j=0;j<my;j++)
        {
            double fxtem=0;
            double fytem=0;
            for(int k=1;k<9;k++)
            {
                int i1=(i+cx[k]+mx)%mx;
                int j1=(j+cy[k]+my)%my;
                // std::cout<<"i1="<<i1<<std::endl;
                // std::cout<<"j1="<<i1<<std::endl;
                double dense= rho[i1][j1];
                // std::cout<<"dense="<<dense<<std::endl;
                double psinb=psi(dense);//psi(x+c*deltat)
                // std::cout<<"psinb="<<psinb<<std::endl;
                fxtem += w[k] * cx[k] * psinb;//psi(x+c*deltat)*w*c
                fytem += w[k] * cy[k] * psinb;
                // std::cout<<"fxtem="<<fxtem<<std::endl;
                // std::cout<<"fytem="<<fytem<<std::endl;
            }
            double psiloc = psi(rho[i][j]);//psi(x)
            fxtem *= (-g * psiloc);//g*psi(x)*psi(x+c*deltat)*w*c
            fytem *= (-g * psiloc);
            Fx[i][j] = fxtem;
            Fy[i][j] = fytem;
        }
    }
}
void ComputeVelocity(double ftem[n][mx][my])
{   
    for(int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
            ux[i][j] = 0.;
            uy[i][j] = 0.;
            for(int k=0;k<9;k++)
            {
                ux[i][j] += ftem[k][i][j]*cx[k];
                uy[i][j] += ftem[k][i][j]*cy[k];
            }
                ux[i][j] += (0.5 * Fx[i][j]);
                uy[i][j] += (0.5 * Fy[i][j]);
            double dens = ComputeTotalDensity(i,j);
            ux[i][j] /= dens;
            uy[i][j] /= dens;
        }
    }
}
void ComputeEquilibrium(int i,int j)
{
    double dens = rho[i][j];
    double vx = ux[i][j];
    double vy = uy[i][j];
    double usq = vx * vx + vy * vy;
    // std::cout<<"dens= "<<dens<<std::endl;
    // std::cout<<"vx= "<<vx<<std::endl;
    // std::cout<<"vy= "<<vy<<std::endl;
    // std::cout<<"usq= "<<usq<<std::endl;
    feq[0] = w[0] * dens * (1. - 1.5 * usq);
    feq[1] = w[1] * dens * (1. + 3. * vx + 4.5 * vx * vx - 1.5 * usq);
    feq[2] = w[2] * dens * (1. + 3. * vy + 4.5 * vy * vy - 1.5 * usq);
    feq[3] = w[3] * dens * (1. - 3. * vx + 4.5 * vx * vx - 1.5 * usq);
    feq[4] = w[4] * dens * (1. - 3. * vy + 4.5 * vy * vy - 1.5 * usq);
    feq[5] = w[5] * dens * (1. + 3. * ( vx + vy) + 4.5 * ( vx + vy) * ( vx + vy) - 1.5 * usq);
    feq[6] = w[6] * dens * (1. + 3. * (-vx + vy) + 4.5 * (-vx + vy) * (-vx + vy) - 1.5 * usq);
    feq[7] = w[7] * dens * (1. + 3. * (-vx - vy) + 4.5 * ( vx + vy) * ( vx + vy) - 1.5 * usq);
    feq[8] = w[8] * dens * (1. + 3. * ( vx - vy) + 4.5 * ( vx - vy) * ( vx - vy) - 1.5 * usq);
}
void ComputeGuoForce(int i,int j,int k)
{
        double temp=cx[k] * ux[i][j] + cy[k] * uy[i][j];
        forcing[k] = w[k] * (1. - 0.5 * omega) * 
        ((3. * (cx[k] - ux[i][j]) + 9. * cx[k] * temp) * Fx[i][j] + 
        (3. * (cy[k] - uy[i][j]) + 9. * cy[k] * temp) * Fy[i][j]);
}
void Collision(double ftem[n][mx][my],double ftem2[n][mx][my])
{
    for(int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
            for (int k=0;k<9;k++)
            {
                ComputeEquilibrium(i,j);
                ComputeGuoForce(i,j,k);
                int i2 = (i + cx[k] + mx) % mx;
                int j2 = (j + cy[k] + my) % my;
                ftem2[k][i2][j2] = ftem[k][i][j] * (1. - omega) + feq[k] * omega + forcing[k] ;
            }
        }
    }
}
void Streaming(double ftem[n][mx][my])
{
    double f_hlp[9][mx][my];
    for(int j=0;j<my;j++)
    {
        for (int i=0;i<mx;i++)
        {
            int y_n = (1+j)%my;
            int x_e = (1+i)%mx;
            int y_s = my - 1 - (my- j)%my;
            int x_w = mx - 1 - (mx- i)%mx;
            f_hlp[1][x_e][j] = ftem[1][i][j];
            f_hlp[2][i][y_n] = ftem[2][i][j];
            f_hlp[3][x_w][j] = ftem[3][i][j];
            f_hlp[4][i][y_s] = ftem[4][i][j];
            f_hlp[5][x_e][y_n] = ftem[5][i][j];
            f_hlp[6][x_w][y_n] = ftem[6][i][j];
            f_hlp[7][x_w][y_s] = ftem[7][i][j];
            f_hlp[8][x_e][y_s] = ftem[8][i][j];
        }
    }
    for(int j=0;j<my;j++)
    {
        for (int i=0;i<mx;i++)
        {
            for(int k=1;k<9;k++)
            {
                ftem[k][i][j]=f_hlp[k][i][j];
            }
        }
    }
}
void BoundaryCondition(double ftem[n][mx][my])
{
    for(int j=0;j<my;j++)
    {
        ftem[1][0][j]=ftem[3][0][j];//left
        ftem[5][0][j]=ftem[7][0][j];
        ftem[8][0][j]=ftem[6][0][j];

        ftem[3][mx-1][j]=ftem[1][mx-1][j];//right
        ftem[7][mx-1][j]=ftem[5][mx-1][j];
        ftem[6][mx-1][j]=ftem[8][mx-1][j];
    }

    for (int i=0;i<mx;i++)
    {
        ftem[2][i][0]=ftem[4][i][0];//bottom
        ftem[5][i][0]=ftem[7][i][0];
        ftem[6][i][0]=ftem[8][i][0];

        ftem[4][i][my-1]=ftem[2][i][my-1];//top
        ftem[7][i][my-1]=ftem[5][i][my-1];
        ftem[8][i][my-1]=ftem[6][i][my-1];
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

std::cout<<"omega= "<<omega<<std::endl;
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
                rho[i][j]=rho_l;
            }
            else
            rho[i][j]=rho_g;
                break;
            }
            case 2:
            {
                if((x[i]-mx/2+radius-1)*(x[i]-mx/2+radius-1)+(y[j]-my/2)*(y[j]-my/2)<radius*radius
            ||(x[i]-mx/2-radius+1)*(x[i]-mx/2-radius+1)+(y[j]-my/2)*(y[j]-my/2)<radius*radius)
            {
                rho[i][j]=rho_g;
            }
            else
            rho[i][j]=rho_l;
                break;
            }
            case 3:
            {
                if((x[i]-mx/2)*(x[i]-my/2)/(a*a)+(y[j]-mx/2)*(y[j]-my/2)/(b*b)<1)
            {
                rho[i][j]=rho_l;
            }
            else
            rho[i][j]=rho_g;
                break;
            }
            }
        }
    }
//initial condition of distribution function 
    for(int j=0;j<my;j++)
    {
       for(int i=0;i<mx;i++)
       {
           for(int k=0;k<9;k++)
           {
               ux[i][j]=0;
               uy[i][j]=0;
               ComputeEquilibrium(i,j);
               f[k][i][j]=feq[k];
               f2[k][i][j]=feq[k];
           }
       }
    }
/*initial condition--------------*/
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
//BoundaryCondition(f);
std::string filename = std::string("animation") +std::to_string(time)+std::string(".dat");
if(time%freq==0) result(filename,time);
if(time%freq==0) {
    std::cout << "Iteration: " << time << " / " << mstep << " (" << 100.0*double(time)/double(mstep) << " %)" <<  std::endl;
}
}//main loop end
return 0;
}