#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath> // for fabs()
int const static n=9,mx=200,my=200; //number of latttice nodes
int c=1;//different cases
double f[n][mx][my],f2[n][mx][my],feq[n],rho[mx][my],w[n],ux[mx][my],uy[mx][my],x[mx],y[my];
const int cx[n]={0, 1, 0, -1, 0, 1, -1, -1, 1};
const int cy[n]={0, 0, 1, 0, -1, 1, 1, -1, -1};
double press[mx][my];
double Fx[mx][my]; // x-component of Shan-Chen force
double Fy[mx][my]; // y-component of Shan-Chen force
double forcing[n];
double rho_l= 1.5,rho_g=0.25;
double radius=50;
double a=mx/10,b=my/5;
int i,j;
int dx=1,dy=1; //space and time step
// double const alpha=0.4;
// double omega=1.0/(3.*alpha+0.5);
double omega=1.0;
int mstep=6000; // The total number of time steps
int freq=100;
double rho0=1.0;//reference density 
const double g = -4.5;//interaction strength between particles.
void result(std::string filename,int time)
{

std::ofstream out(filename);
  out << "TITLE=\"LBM output\"" << std::endl;
  out << "VARIABLES = \"X\", \"Y\", \"Ux\",\"Uy\",\"density\",\"pressure\",\"f0\",\"f1\",\"f2\",\"f3\",\"f4\",\"f5\",\"f6\",\"f7\",\"f8\",\"feq0\",\"feq1\",\"feq2\",\"feq3\",\"feq4\",\"feq5\",\"feq6\",\"feq7\",\"feq8\"," << std::endl;
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
        out << std::scientific << std::setprecision(5) << std::setw(15) << press[i][j];
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
void Non_Equilibrium_extrapolation(double ftemp[n][mx][my])
{
    
    //west
    for(int j=0;j<my;j++)
    {
        ux[0][j]=ux[1][j];
        uy[0][j]=ux[1][j];
        rho[0][j]=1.6;
        for (int k = 0; k < 9; k++)
        {
            ComputeEquilibrium(0,j);
            ftemp[k][0][j]=feq[k]+ftemp[k][1][j];
            ComputeEquilibrium(1,j);
            ftemp[k][0][j]-=feq[k];
        }
    //east
        ux[mx-1][j]=ux[mx-2][j];
        uy[mx-1][j]=ux[mx-2][j];
        rho[mx-1][j]=1.6;
        for (int k = 0; k < 9; k++)
        {
            ComputeEquilibrium(mx-1,j);
            ftemp[k][mx-1][j]=feq[k]+ftemp[k][mx-2][j];
            ComputeEquilibrium(mx-2,j);
            ftemp[k][mx-1][j]-=feq[k];
        }
    }
    //north
    // for (int i=1;i<mx-1;i++)
    // {
    //     // ux[i][my-1]=0.0;
    //     // uy[i][my-1]=0.0;
    //     rho[i][my-1]=1.6;
    //     for (int k = 0; k < 9; k++)
    //     {
    //         ComputeEquilibrium(i,my-1);
    //         ftemp[k][i][my-1]=feq[k]+ftemp[k][i][my-2];
    //         ComputeEquilibrium(i,my-2);
    //         ftemp[k][i][my-1]-=feq[k];
    //     }
    //     //south
    //     // ux[i][0]=0.0;
    //     // uy[i][0]=0.0;
    //     rho[i][0]=1.6;
    //     for (int k = 0; k < 9; k++)
    //     {
    //         ComputeEquilibrium(i,0);
    //         ftemp[k][i][0]=feq[k]+ftemp[k][i][1];
    //         ComputeEquilibrium(i,1);
    //         ftemp[k][i][0]-=feq[k];
    //     }
    // }
    // corner
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
            double w=10;
            double factor1 = 0.5 * (rho_l-rho_g);
            double temp=0.5 * (rho_l+rho_g);
            double factor2 = sqrt(pow((i-0.5*mx), 2) + pow((j-0.5*my), 2));
            factor2 = 2.0 * (factor2 - radius) / w;
            switch (c)
            {
            case 1:
            {
            //smooth the interface 
            //droplet
            //rho[i][j] = temp-factor1 * tanh(factor2);
            //bubble 
            rho[i][j] = temp+factor1 * tanh(factor2);
            break;
            //  if((x[i]-mx/2)*(x[i]-my/2)+(y[j]-mx/2)*(y[j]-my/2)<radius*radius)
            // {
            //     rho[i][j]=rho_l;
            // }
            // else
            // rho[i][j]=rho_g;
            //     break;
            }
            case 2:
            {
                //1.Wetting boundary conditions for multicomponent pseudopotential lattice Boltzmann
                // double factor3=pow((i-0.5*mx-radius), 2) + pow((j-0.5*my), 2)-radius*radius;
                // double factor4=pow((i-0.5*mx+radius), 2) + pow((j-0.5*my), 2)-radius*radius;
                // rho[i][j]=rho_g+factor1*(1-tanh(factor3/(w*w)))+factor1*(1-tanh(factor4/(w*w)));
                // break;
                //2.initialized in two axisymmetric parts based on case 1. 
                double factor3 = sqrt(pow((i-0.5*mx-radius), 2) + pow((j-0.5*my), 2));
                factor3 = 2.0 * (factor3 - radius) / w;
                double factor4 = sqrt(pow((i-0.5*mx+radius), 2) + pow((j-0.5*my), 2));
                factor4 = 2.0 * (factor4 - radius) / w;
                if (x[i]>mx/2) rho[i][j] = temp+factor1 * tanh(factor3);
                else rho[i][j] = temp+factor1 * tanh(factor4);
            break;
            }
            case 3:
            {
                //Force approach for the pseudopotential lattice Boltzmann method
                double e=sqrt(1-(a/b)*(a/b));
                double theta;
                double R0;
                if (i!=0.5*mx&&j!=0.5*my)
                {
                    theta=atan((j-0.5*my)/(i-0.5*mx));
                    R0=a/sqrt(1-pow(e*cos(theta),2));
                }
                else if(j==0.5*my&&i<0.5*mx) R0=a/sqrt(1-pow(-e,2));
                else if(j==0.5*my&&i>0.5*mx) R0=a/sqrt(1-pow(e,2));
                else if(i==0.5*mx) R0=a/sqrt(1-pow(0,2));
                factor2=sqrt(pow((i-0.5*mx), 2) + pow((j-0.5*my), 2));
                factor2 = 2*(factor2 - R0) / w;
                rho[i][j] = temp-factor1 * tanh(factor2);
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
void Non_Equilibrium_Bounce_Back(double ftemp[n][mx][my])
{
//west
    for (int j=0;j<my;j++)
    {
        uy[0][j]=0.0;
        rho[0][j]=1.4;
        ux[0][j]=1.0-(ftemp[0][0][j]+ftemp[2][0][j]+ftemp[4][0][j]+2.0*(ftemp[3][0][j]+ftemp[6][0][j]+ftemp[7][0][j]))/rho[0][j];
        ftemp[1][0][j]=ftemp[3][0][j]+2.0/3.0*rho[0][j]*ux[0][j];
        ftemp[5][0][j]=ftemp[7][0][j]-1.0/2.0*(ftemp[2][0][j]-ftemp[4][0][j])+1.0/6.0*rho[0][j]*ux[0][j];
        ftemp[8][0][j]=ftemp[6][0][j]-1.0/2.0*(ftemp[2][0][j]-ftemp[4][0][j])+1.0/6.0*rho[0][j]*ux[0][j];
    }
//east
for (int j=0;j<my;j++)
    {
        uy[mx-1][j]=0.0;
        rho[mx-1][j]=1.4;
        ux[mx-1][j]=(ftemp[0][mx-1][j]+ftemp[2][mx-1][j]+ftemp[4][mx-1][j]+2.0*(ftemp[1][mx-1][j]+ftemp[5][mx-1][j]+ftemp[8][mx-1][j]))/rho[mx-1][j]-1;
        ftemp[3][mx-1][j]=ftemp[1][mx-1][j]-2.0/3.0*rho[mx-1][j]*ux[mx-1][j];
        ftemp[6][mx-1][j]=ftemp[8][mx-1][j]-1.0/2.0*(ftemp[2][mx-1][j]-ftemp[4][mx-1][j])-1.0/6.0*rho[mx-1][j]*ux[mx-1][j];
        ftemp[7][mx-1][j]=ftemp[5][mx-1][j]+1.0/2.0*(ftemp[2][mx-1][j]-ftemp[4][mx-1][j])-1.0/6.0*rho[mx-1][j]*ux[mx-1][j];
    }
//north
for (int i=0;i<mx;i++)
    {
        ux[i][my-1]=0.0;
        rho[i][my-1]=1.4;
        uy[i][my-1]=(ftemp[0][i][my-1]+ftemp[1][i][my-1]+ftemp[3][i][my-1]+2.0*(ftemp[2][i][my-1]+ftemp[5][i][my-1]+ftemp[6][i][my-1]))/rho[i][my-1]-1;
        ftemp[4][i][my-1]=ftemp[2][i][my-1]-2.0/3.0*rho[i][my-1]*uy[i][my-1];
        ftemp[7][i][my-1]=ftemp[5][i][my-1]+1.0/2.0*(ftemp[1][i][my-1]-ftemp[3][i][my-1])-1.0/6.0*rho[i][my-1]*uy[i][my-1];
        ftemp[8][i][my-1]=ftemp[6][i][my-1]+1.0/2.0*(ftemp[3][i][my-1]-ftemp[1][i][my-1])-1.0/6.0*rho[i][my-1]*uy[i][my-1];
    }
//south
for (int i=0;i<mx;i++)
    {
        ux[i][0]=0.0;
        rho[i][0]=1.4;
        uy[i][0]=1-(ftemp[0][i][0]+ftemp[1][i][0]+ftemp[3][i][0]+2.0*(ftemp[4][i][0]+ftemp[7][i][0]+ftemp[8][i][0]))/rho[i][0];
        ftemp[2][i][0]=ftemp[4][i][0]+2.0/3.0*rho[i][0]*uy[i][0];
        ftemp[5][i][0]=ftemp[7][i][0]-1.0/2.0*(ftemp[1][i][0]-ftemp[3][i][0])+1.0/6.0*rho[i][0]*uy[i][0];
        ftemp[6][i][0]=ftemp[8][i][0]+1.0/2.0*(ftemp[1][i][0]-ftemp[3][i][0])+1.0/6.0*rho[i][0]*uy[i][0];
    }
}

void VerifyLapalaceLaw()
{
//calculate pressure
for (int j = 0; j < my; j++) {
    for (int i = 0; i < mx; i++) {
        press[i][j] = 0.;
        press[i][j] += rho[i][j] / 3.;
        press[i][j] += (g * psi(rho[i][j]) * psi(rho[i][j]) / 6.);
    }
}
//calculate deltapressure
double rho_gas = 0.;
double rho_liq = 0.;
double press_gas = 0.;
double press_liq = 0.;

// Average liquid pressure in a square box with 4 lattice nodes
for (int i = 0; i < 4; i++) 
{
    for (int j = 0; j < 4; j++) {
        rho_liq += rho[i][j];
        press_liq += press[i][j];
    }
}

  // Average gas pressure in a square box with 4 lattice nodes
  for (int i = mx / 2 - 2; i < mx / 2 + 2; i++) {
    for (int j = my / 2 - 2; j < my / 2 + 2; j++) {
      rho_gas += rho[i][j];
      press_gas += press[i][j];
    }
  }

  rho_gas /= 16.;
  rho_liq /= 16.;
  press_gas /= 16.;
  press_liq /= 16.;
  
  double delta_press = press_gas - press_liq;
  double rho_av = 0.5 * (rho_gas + rho_liq);
  const int y = my / 2;
//calculate radius
    double rad;
    for (int i = 0; i < mx / 2; i++) {
        if (rho[i][y] < rho_av) {
            const double drho = rho[i][y] - rho[i-1][y];
            const double dx = (rho_av - rho[i-1][y]) / drho;
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
    //BoundaryCondition(f);
    VerifyLapalaceLaw();
    //Non_Equilibrium_extrapolation(f);
    Non_Equilibrium_Bounce_Back(f);
    std::string filename = std::string("animation") +std::to_string(time)+std::string(".dat");
    if(time%freq==0) result(filename,time);
    if(time%freq==0)
    {
        std::cout << "Iteration: " << time << " / " << mstep << " (" << 100.0*double(time)/double(mstep) << " %)" <<  std::endl;
    }
}//main loop end
return 0;
}