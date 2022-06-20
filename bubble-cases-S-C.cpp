#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath> // for fabs()
int const static n=9,mx=200,my=200; //number of latttice nodes
int c=2;//different cases
double f[n][mx][my],feq[n],rho[mx][my],cx[n],cy[n],w[n],ux[mx][my],uy[mx][my],x[mx],y[my];
double Fx[mx][my]; // x-component of Shan-Chen force
double Fy[mx][my]; // y-component of Shan-Chen force
double forcing[mx][my];
double rho_l=1.95,rho_g=0.15;
double radius=mx/8;
double a=mx/10,b=my/5;
int i,j;
int dx=1,dy=1; //space and time step
double const alpha=0.02;
double omega=1.0/(3.*alpha+0.5);
int mstep=10; // The total number of time steps
double rho0=1.0;//reference density 
const double g = -4.7;//interaction strength between particles.
void result(double ux[mx][my],double uy[mx][my],double rho[mx][my],double x[mx], double y[my],std::string filename,double f[9][mx][my],double feq[9],int time)
{

std::ofstream out(filename);
  out << "TITLE=\"LBM output\"" << std::endl;
  out << "VARIABLES = \"X\", \"Y\", \"Ux\",\"Uy\",\"density\",\"f0\",\"f1\",\"f2\",\"f3\",\"f4\",\"f5\",\"f6\",\"f7\",\"f8\",\"feq0\",\"feq1\",\"feq2\",\"feq3\",\"feq4\",\"feq5\",\"feq6\",\"feq7\",\"feq8\"," << std::endl;
  out << "ZONE T = \"fluid\", I=" << mx << ", J=" << my << ", F=POINT" << std::endl;
  out << "SOLUTIONTIME="<< time << std::endl;
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
    double fxtem=0;
    double fytem=0;
    for(int i=1;i<mx-1;i++)
    {
        for (int j=1;j<my-1;j++)
        {
            for(int k=0;k<9;k++)
            {

                double dense= rho[i][j];
                double psinb=psi(dense);//psi(x+c*deltat)
                fxtem += w[k] * cx[k] * psinb;//psi(x+c*deltat)*w*c
                fytem += w[k] * cy[k] * psinb;
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
void ComputeGuoForce(int i,int j)
{
    for (i = 0; i < mx; i++) 
    {
        for(j=0;j<my;j++)
        {
            forcing[i][j] = w[i] * (1. - 0.5 * omega) * ((3. * (cx[i] - ux[i][j]) + 9. * cx[i] * (cx[i] * ux[i][j] + cy[i] * uy[i][j])) * Fx[i][j] + (3. * (cy[i] - uy[i][j]) + 9. * cy[i] * (cx[i] * ux[i][j] + cy[i] * uy[i][j])) * Fy[i][j]);
        }
    }
}
void Collision(double ftem[n][mx][my])
{
    for(int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
            for (int k=0;k<9;k++)
            {
                ftem[k][i][j] = ftem[k][i][j] * (1. - omega) + feq[k] * omega + forcing[i][j] ;
            }
        }
    }
}
void Streaming(double ftem[n][mx][my])
{
    for (int j=0;j<my;j++)
    {
        for(int i=mx-1;i>0;i--)
        {
            f[1][i][j]=f[1][i-1][j];
        }
        for(int i=0;i<mx-1;i++)
        {
            f[3][i][j]=f[3][i+1][j];
        }
    }

    for(int j=my-1;j>0;j--)
    {
        for(int i=0;i<mx;i++)
        {
            f[2][i][j]=f[2][i][j-1];
        }
        for(int i=mx-1;i>0;i--)
        {
            f[5][i][j]=f[5][i-1][j-1];
        }
        for(int i=0;i<mx-1;i++)
        {
            f[6][i][j]=f[6][i+1][j-1];
        }
    }

    for(int j=0;j<my;j++)
    {
        for(int i=0;i<mx;i++)
        {
            f[4][i][j]=f[4][i][j+1];
        }
        for(int i=0;i<mx-1;i++)
        {
            f[7][i][j]=f[7][i+1][j+1];
        }
        for(int i=mx-1;i>0;i--)
        {
            f[8][i][j]=f[8][i-1][j+1];
        }
    }
}
void BoundaryCondition(double ftem[n][mx][my])
{
    for(int j=0;j<my;j++)
    {
        f[1][0][j]=f[3][0][j];//left
        f[5][0][j]=f[7][0][j];
        f[8][0][j]=f[6][0][j];

        f[3][mx-1][j]=f[1][mx-1][j];//right
        f[7][mx-1][j]=f[5][mx-1][j];
        f[6][mx-1][j]=f[8][mx-1][j];
    }

    for (int i=0;i<mx;i++)
    {
        f[2][i][0]=f[4][i][0];//bottom
        f[5][i][0]=f[7][i][0];
        f[6][i][0]=f[8][i][0];

        f[4][i][my-1]=f[2][i][my-1];//top
        f[7][i][my-1]=f[5][i][my-1];
        f[8][i][my-1]=f[6][i][my-1];
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

/*---------streaming vector--------*/
cx[0]=0,cx[1]=1,cx[2]=0,cx[3]=-1,cx[4]=0,cx[5]=1,cx[6]=-1,cx[7]=-1,cx[8]=1;
cy[0]=0,cy[1]=0,cy[2]=1,cy[3]=0,cy[4]=-1,cy[5]=1,cy[6]=1,cy[7]=-1,cy[8]=-1;
/*---------streaming vector--------*/

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
               ComputeEquilibrium(i,j);
               f[k][i][j]=feq[k];
               ux[i][j]=0;
               uy[i][j]=0;
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
result(ux,uy,rho,x, y,filename,f,feq,time);
for (time=1;time<mstep+1;++time)
{
ComputeSCForce();//only novelty which contribute to equilibrium velocity and macroscopic velocity 
ComputeVelocity(f);
ComputeEquilibrium(i,j);
ComputeGuoForce(i,j);
Collision(f);//plus GuoForce 
Streaming(f);
BoundaryCondition(f);
ComputeDesnity(f);
ComputeVelocity(f);
std::string filename = std::string("animation") +std::to_string(time)+std::string(".dat");
result(ux,uy,rho,x, y,filename,f,feq,time);
}//main loop end
return 0;
}




