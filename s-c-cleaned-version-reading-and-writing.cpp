#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath> // for fabs()
char setupfile[]="setupfile.txt";
int n,mx,my; //number of latttice nodes
int cas;//different cases
int bc;//1:periodic boundary condition 2:zou-he boundary condition
int eos;//different equation of state :: eos=1 and default standard eos. eos=2 C-S equation of state. 
double R;//constant in the c-s eos
double A;//constant in the c-s eos
double B;//constant in the c-s eos
double T;//temperature in c-s eos 
double lx,ly;//lattice length 
double dx,dy; //space
double dt;//time step
int i11=0;
int j11=0;
int c;//lattice speed
double rho_l,rho_g;//liquid and gas density where liquid density changes around equilibrium density to trigger growth and collpase
//equilibrium density rho_l=0.33 and rho_g=0.011 
double radius;
double a=mx/10,b=my/5;//pararmeter for elliptical case
// double const alpha=0.4;
// double omega=1.0/(3.*alpha+0.5);
double omega;
int mstep;// The total number of time steps
int freq;//frequency for output
double rho0;//reference density 
double g;//interaction strength between particles.
double if_th;//interface thickness
std::string to_string_with_precision(const double a_value, const int n = 2)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}
void Readinitialfile(char* filename,int *n,int *mx,int *my,double *omega,int *mstep
,int* cas,int *bc,int* eos,double* R,double* A,double* B,double* T,double *lx,double* ly,int*c,
double* rho_l,double* rho_g,double* radius,int* freq,double* rho0,double* g,double* if_th)
{ 
    FILE *file=NULL;
    file=fopen(filename,"r");
    if(!file)
    {
        std::cout<<"cannot find file"<<std::endl;
    }
    fscanf(file,"%d",n);
    fscanf(file,"%d",mx);
    fscanf(file,"%d",my);
    fscanf(file,"%lf",omega);
    fscanf(file,"%d",mstep);
    fscanf(file,"%d",cas);
    fscanf(file,"%d",bc);
    fscanf(file,"%d",eos);
    fscanf(file,"%lf",R);
    fscanf(file,"%lf",A);
    fscanf(file,"%lf",B);
    fscanf(file,"%lf",T);
    fscanf(file,"%lf",lx);
    fscanf(file,"%lf",ly);
    fscanf(file,"%d",c);
    fscanf(file,"%lf",rho_l);
    fscanf(file,"%lf",rho_g);
    fscanf(file,"%lf",radius);
    fscanf(file,"%d",freq);
    fscanf(file,"%lf",rho0);
    fscanf(file,"%lf",g);
    fscanf(file,"%lf",if_th);
    fclose(file);
}
void result(std::string filename,int time,double *x,double *y,double** ux,double **uy,double** rho,double **press
,double*** f, double* feq)
{

std::ofstream out(filename);
  out << "TITLE=\"LBM output\"" << std::endl;
  out << "VARIABLES = \"X\", \"Y\", \"Ux\",\"Uy\",\"density\",\"pressure\",\"f0\",\"f1\",\"f2\",\"f3\",\"f4\",\"f5\",\"f6\",\"f7\",\"f8\",\"feq0\",\"feq1\",\"feq2\",\"feq3\",\"feq4\",\"feq5\",\"feq6\",\"feq7\",\"feq8\"," << std::endl;
  out << "ZONE T = \"fluid\", I=" << mx << ", J=" << my << ", F=POINT" << std::endl;
  //out << "SOLUTIONTIME="<< time << std::endl;
  for (int j = 0; j < my; ++j)
  {
    for (int i = 0; i < mx; ++i)
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
  }
  out.close();
}

void ComputeDesnity(double ***f,double **rho)
{
    for(int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
           double ssum=0;
           for(int k=0;k<9;k++)
           {
               ssum=ssum+f[k][i][j];
           }
           rho[i][j]=ssum;
        }
    }
}
double ComputeTotalDensity(int i,int j,double **rho)
{
    double dense = rho[i][j];
    return dense;
}
double psi(double dense)
{
    switch (eos)
    {
        case 2://C-S EOS
            {
                double p=0;
                p=dense*R*T*(1+B*dense/4+(B*dense/4)*(B*dense/4)-pow(B*dense/4,3.0))/pow((1-B*dense/4),3.0)-A*dense*dense;
                return sqrt(2/(g*1/3)*(p-dense*1/3));
                break;
            }
        default://original shan-chen psi 
        {
            return rho0 * (1. - exp(-dense / rho0));
            break;
        }
    }
}
void ComputeSCForce(int *cx,int *cy,double **rho,double *w,double **Fx,double **Fy)
{
  for(int i=0;i<mx;i++)
    {
        for (int j=0;j<my;j++)
        {
            double fxtem=0;
            double fytem=0;
            // int i11=0;
            // int j11=0;
            for(int k=0;k<9;k++)
            {
                if(bc==1)
                {
                //lattice speed cx,cy=1 and i11,j11 is just an index so it is no need to use double 
                i11 = (i + cx[k] + mx) % (mx);
                j11 = (j + cy[k] + my) % (my);
                }
                else
                {
                    i11 = i + cx[k];
                    if (i11 < 0) {
                        i11 = 0;
                    } else if (i11 >= mx) {
                        i11 = mx - 1;
                    }
                    j11 = j + cy[k];
                    if (j11 < 0) {
                        j11 = 0;
                    } else if (j11 >= my) {
                        j11 = my - 1;
                    }
                    // std::cout<<"i11 before"<<i11<<std::endl;
                }
                // std::cout<<"i11 after"<<i11<<std::endl;
                double dense= rho[i11][j11];
                //neighbor psi 
                double psinb=psi(dense);//psi(x+c*deltat)
                fxtem += w[k] * cx[k] * psinb;//psi(x+c*deltat)*w*c
                fytem += w[k] * cy[k] * psinb;
            }
            //local psi
            double psiloc = psi(rho[i][j]);//psi(x)
            fxtem *= (-g * psiloc);//g*psi(x)*psi(x+c*deltat)*w*c
            fytem *= (-g * psiloc);
            Fx[i][j] = fxtem;
            Fy[i][j] = fytem;
        }
    }
}
void ComputeVelocity(double ***f,double** ux,double **uy,int* cx,int* cy,double **rho)
{   
    for(int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
            ux[i][j] = 0.;
            uy[i][j] = 0.;
            for(int k=0;k<9;k++)
            {
                ux[i][j] += f[k][i][j]*cx[k];
                uy[i][j] += f[k][i][j]*cy[k];
            }
            double dens = ComputeTotalDensity(i,j,rho);
            ux[i][j] /= dens;
            uy[i][j] /= dens;
        }
    }
}
void ComputeEquilibrium(int i,int j,double **ux,double** uy,double ** rho,double **Fx,double **Fy,int* cx
,int *cy, double *w,double *feq)
{
    double dens = rho[i][j];
    //implement the vsm force scheme here 
    double vx = ux[i][j]+Fx[i][j]/(omega*rho[i][j]);
    double vy = uy[i][j]+Fy[i][j]/(omega*rho[i][j]);
    double usq = vx * vx + vy * vy;
    for(int k=0;k<9;k++)
    {
        feq[k]=w[k]*dens*(1.0+3.0*(cx[k]*vx+cy[k]*vy)/(c*c)+4.5*(cx[k]*vx+cy[k]*vy)*(cx[k]*vx+cy[k]*vy)/(c*c*c*c)-
        1.5*usq/(c*c));
    }
}
void Collision(double ***f,double **ux,double** uy,double ** rho,double **Fx,double **Fy,int* cx
,int *cy, double *w,double *feq)
{
    for(int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
            ComputeEquilibrium(i,j,ux,uy,rho,Fx,Fy,cx,cy,w,feq);
            for (int k=0;k<9;k++)
            {
                f[k][i][j] = f[k][i][j] * (1. - omega) + feq[k] * omega;
            }
        }
    }
}
void Streaming(double ***f)
{
    for (int j=0;j<my;j++)
    {
        for(int i=mx-1;i>0;i--)
        {
            f[1][i][j]=f[1][i-1][j];
        }
        for(int i=0;i<mx-1;i++)
        {
            f[3][i][j]=f[3][i+1][j];//?
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

    for(int j=0;j<my-1;j++)
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
void initialize(int *cx,int *cy,double *w,double *x,double *y,double **rho,double **ux,double **uy,
double ***f,double *feq,double** Fx,double** Fy,double if_th)
{
dx=lx/(mx-1);
dy=ly/(my-1); //space
dt=dx;//time step
mstep=mstep/dt;
std::cout<<"dx "<<dx<<std::endl;
std::cout<<"dy "<<dy<<std::endl;
std::cout<<"dt "<<dt<<std::endl;
std::cout<<"mstep "<<mstep<<std::endl;
std::string filename = std::string("animation0") +std::string(".dat");
/*---------streaming vector--------*/
cx[0]=0,cx[1]=1,cx[2]=0,cx[3]=-1,cx[4]=0,cx[5]=1,cx[6]=-1,cx[7]=-1,cx[8]=1;
cy[0]=0,cy[1]=0,cy[2]=1,cy[3]=0,cy[4]=-1,cy[5]=1,cy[6]=1,cy[7]=-1,cy[8]=-1;
/*---------streaming vector--------*/
/*-------weight factor ------*/
w[0]=4./9;
    for(int i=1;i<5;i++)
    {
        w[i]=1./9;
    }
    for(int i=5;i<9;i++)
    {
        w[i]=1./36;
    }
/*---------weight factor ----------*/
//coordinate
x[0] =0.0; 
y[0] =0.0;
for (int i=1;i<mx;i++)
{
    x[i]=dx*i;
}
for (int j=1;j<my;j++)
{
    y[j]=dy*j;
}
/*initial condition--------------*/
// std::cout<<"n "<<n<<std::endl;
// std::cout<<"mx "<<mx<<std::endl;
// std::cout<<"my "<<my<<std::endl;
// std::cout<<"omega= "<<omega<<std::endl;
// std::cout<<"mstep= "<<mstep<<std::endl;
// std::cout<<"cas= "<<cas<<std::endl;
// std::cout<<"bc= "<<bc<<std::endl;
// std::cout<<"eos= "<<eos<<std::endl;
// std::cout<<"R= "<<R<<std::endl;
// std::cout<<"A= "<<A<<std::endl;
// std::cout<<"B= "<<B<<std::endl;
// std::cout<<"T= "<<T<<std::endl;
// std::cout<<"lx= "<<lx<<std::endl;
// std::cout<<"ly= "<<ly<<std::endl;
// std::cout<<"c= "<<c<<std::endl;
// std::cout<<"rho_l= "<<rho_l<<std::endl;
// std::cout<<"rho_g= "<<rho_g<<std::endl;
// std::cout<<"radius "<<radius<<std::endl;
// std::cout<<"freq "<<freq<<std::endl;
// std::cout<<"rho0 "<<rho0<<std::endl;
// std::cout<<"g "<<g<<std::endl;
// std::cout<<"if_th "<<if_th<<std::endl;
/*initial condition--------------*/
    for(int i=0;i<mx;i++)
    {
        for (int j=0;j<my;j++)
        {
            double factor1 = 0.5 * (rho_l-rho_g);
            double temp=0.5 * (rho_l+rho_g);
            double factor2 = sqrt(pow((x[i]-0.5*lx), 2) + pow((y[j]-0.5*ly), 2));
            factor2 = 2.0 * (factor2 - radius) / if_th;
            // std::cout<<"factor1= "<<rho_g<<std::endl;
            // std::cout<<"temp "<<radius<<std::endl;
            // std::cout<<"factor2 "<<freq<<std::endl;
            switch (cas)
            {
            case 1://single bubble or droplet with smooth interface
            {
                //smooth the interface 
                //droplet
                //rho[i][j] = temp-factor1 * tanh(factor2);
                //bubble 
                rho[i][j] = temp+factor1 * tanh(factor2);
                break;
            }
            case 2: //two droplets or bubbles merge with smmoth interface 
            {
                //2.initialized in two axisymmetric parts based on case 1. 
                double factor3 = sqrt(pow((x[i]-0.5*lx-radius), 2) + pow((y[j]-0.5*ly), 2));
                factor3 = 2.0 * (factor3 - radius) / if_th;
                double factor4 = sqrt(pow((x[i]-0.5*lx+radius), 2) + pow((y[j]-0.5*ly), 2));
                factor4 = 2.0 * (factor4 - radius) / if_th;
                //bubble
                rho[i][j] = temp+factor1 * tanh(factor3)*tanh(factor4);
                //droplet
                //rho[i][j] = temp-factor1 * tanh(factor3)*tanh(factor4);
            break;
            }
            case 3: //elliptic droplet or bubble with smooth interface
            {
                //Force approach for the pseudopotential lattice Boltzmann method
                double e=sqrt(1-(a/b)*(a/b));
                double theta;
                double R0;
                if (x[i]!=0.5*lx&&y[j]!=0.5*ly)
                {
                    theta=atan((y[j]-0.5*ly)/(x[i]-0.5*lx));
                    R0=a/sqrt(1-pow(e*cos(theta),2));
                }
                else if(y[j]==0.5*ly&&x[i]<0.5*lx) R0=a/sqrt(1-pow(-e,2));
                else if(y[j]==0.5*ly&&x[i]>0.5*lx) R0=a/sqrt(1-pow(e,2));
                else if(x[i]==0.5*lx) R0=a/sqrt(1-pow(0,2));
                factor2=sqrt(pow((x[i]-0.5*lx), 2) + pow((y[j]-0.5*ly), 2));
                factor2 = 2*(factor2 - R0) / if_th;
                //droplet
                // rho[i][j] = temp-factor1 * tanh(factor2);
                //bubble
                rho[i][j] = temp+factor1 * tanh(factor2);
                break;
            }
            case 4://single bubble or droplet without smooth interface
            {
                 if((x[i]-lx/2)*(x[i]-ly/2)+(y[j]-lx/2)*(y[j]-ly/2)<radius*radius)
                {
                    rho[i][j]=rho_l;
                }
                else
                rho[i][j]=rho_g;
                    break;
            }
            case 5: //two bubbles or droplets merge without smooth interface
            {
                if((x[i]-lx/2+radius-1)*(x[i]-lx/2+radius-1)+(y[j]-ly/2)*(y[j]-ly/2)<radius*radius
            ||(x[i]-lx/2-radius+1)*(x[i]-lx/2-radius+1)+(y[j]-ly/2)*(y[j]-ly/2)<radius*radius)
            {
                rho[i][j]=rho_g;
            }
            else
            rho[i][j]=rho_l;
                break;
            }
            case 6: //elliptical bubble or droplet without smooth interface
            {
                if((x[i]-lx/2)*(x[i]-ly/2)/(a*a)+(y[j]-lx/2)*(y[j]-ly/2)/(b*b)<1)
            {
                rho[i][j]=rho_g;
            }
            else
            rho[i][j]=rho_l;
                break;
            }
            case 7: //bubble cluster mx=101
            {
                if(//first part
                (x[i]-30)*(x[i]-30)+(y[j]-93)*(y[j]-93)<radius*radius||(x[i]-43)*(x[i]-43)+(y[j]-93)*(y[j]-93)<radius*radius||(x[i]-56)*(x[i]-56)+(y[j]-93)*(y[j]-93)<radius*radius||(x[i]-69)*(x[i]-69)+(y[j]-93)*(y[j]-93)<radius*radius
                ||(x[i]-21)*(x[i]-21)+(y[j]-78)*(y[j]-78)<radius*radius||(x[i]-35.5)*(x[i]-35.5)+(y[j]-78)*(y[j]-78)<radius*radius||(x[i]-50)*(x[i]-50)+(y[j]-78)*(y[j]-78)<radius*radius||(x[i]-64.5)*(x[i]-64.5)+(y[j]-78)*(y[j]-78)<radius*radius||(x[i]-79)*(x[i]-79)+(y[j]-78)*(y[j]-78)<radius*radius
                ||(x[i]-13)*(x[i]-13)+(y[j]-64)*(y[j]-64)<radius*radius||(x[i]-27.8)*(x[i]-27.8)+(y[j]-64)*(y[j]-64)<radius*radius||(x[i]-42.6)*(x[i]-42.6)+(y[j]-64)*(y[j]-64)<radius*radius||(x[i]-57.4)*(x[i]-57.4)+(y[j]-64)*(y[j]-64)<radius*radius||(x[i]-72.2)*(x[i]-72.2)+(y[j]-64)*(y[j]-64)<radius*radius||(x[i]-87)*(x[i]-87)+(y[j]-64)*(y[j]-64)<radius*radius
                ||(x[i]-5)*(x[i]-5)+(y[j]-50)*(y[j]-50)<radius*radius||(x[i]-20)*(x[i]-20)+(y[j]-50)*(y[j]-50)<radius*radius||(x[i]-35)*(x[i]-35)+(y[j]-50)*(y[j]-50)<radius*radius||(x[i]-50)*(x[i]-50)+(y[j]-50)*(y[j]-50)<radius*radius||(x[i]-65)*(x[i]-65)+(y[j]-50)*(y[j]-50)<radius*radius||(x[i]-80)*(x[i]-80)+(y[j]-50)*(y[j]-50)<radius*radius||(x[i]-94)*(x[i]-94)+(y[j]-50)*(y[j]-50)<radius*radius
                //second part
                ||(x[i]-13)*(x[i]-13)+(y[j]-35)*(y[j]-35)<radius*radius||(x[i]-27.8)*(x[i]-27.8)+(y[j]-35)*(y[j]-35)<radius*radius||(x[i]-42.6)*(x[i]-42.6)+(y[j]-35)*(y[j]-35)<radius*radius||(x[i]-57.4)*(x[i]-57.4)+(y[j]-35)*(y[j]-35)<radius*radius||(x[i]-72.2)*(x[i]-72.2)+(y[j]-35)*(y[j]-35)<radius*radius||(x[i]-87)*(x[i]-87)+(y[j]-35)*(y[j]-35)<radius*radius
                ||(x[i]-21)*(x[i]-21)+(y[j]-21)*(y[j]-21)<radius*radius||(x[i]-35.5)*(x[i]-35.5)+(y[j]-21)*(y[j]-21)<radius*radius||(x[i]-50)*(x[i]-50)+(y[j]-21)*(y[j]-21)<radius*radius||(x[i]-64.5)*(x[i]-64.5)+(y[j]-21)*(y[j]-21)<radius*radius||(x[i]-79)*(x[i]-79)+(y[j]-21)*(y[j]-21)<radius*radius
                ||(x[i]-30)*(x[i]-30)+(y[j]-7)*(y[j]-7)<radius*radius||(x[i]-43)*(x[i]-43)+(y[j]-7)*(y[j]-7)<radius*radius||(x[i]-56)*(x[i]-56)+(y[j]-7)*(y[j]-7)<radius*radius||(x[i]-69)*(x[i]-69)+(y[j]-7)*(y[j]-7)<radius*radius)
            {
                rho[i][j]=rho_g;
            }
            else
            rho[i][j]=rho_l;
                break;
            }
            case 8: //bubble cluster with smooth interface
            {   
                //first row
                double factor5 = sqrt(pow((i-30), 2) + pow((j-93), 2));
                factor5 = 2.0 * (factor5 - radius) / if_th;
                double factor6 = sqrt(pow((i-43), 2) + pow((j-93), 2));
                factor6 = 2.0 * (factor6 - radius) / if_th;
                double factor7 = sqrt(pow((i-56), 2) + pow((j-93), 2));
                factor7 = 2.0 * (factor7 - radius) / if_th;
                double factor8 = sqrt(pow((i-69), 2) + pow((j-93), 2));
                factor8 = 2.0 * (factor8 - radius) / if_th;
                //second row
                double factor9 = sqrt(pow((i-21), 2) + pow((j-78), 2));
                factor9 = 2.0 * (factor9 - radius) / if_th;
                double factor10 = sqrt(pow((i-35.5), 2) + pow((j-78), 2));
                factor10 = 2.0 * (factor10 - radius) / if_th;
                double factor11 = sqrt(pow((i-50), 2) + pow((j-78), 2));
                factor11 = 2.0 * (factor11 - radius) / if_th;
                double factor12 = sqrt(pow((i-64.5), 2) + pow((j-78), 2));
                factor12 = 2.0 * (factor12 - radius) / if_th;
                double factor13 = sqrt(pow((i-79), 2) + pow((j-78), 2));
                factor13 = 2.0 * (factor13 - radius) / if_th;
                //third row
                double factor14 = sqrt(pow((i-13), 2) + pow((j-64), 2));
                factor14 = 2.0 * (factor14 - radius) / if_th;
                double factor15 = sqrt(pow((i-27.8), 2) + pow((j-64), 2));
                factor15 = 2.0 * (factor15 - radius) / if_th;
                double factor16 = sqrt(pow((i-42.6), 2) + pow((j-64), 2));
                factor16 = 2.0 * (factor16 - radius) / if_th;
                double factor17 = sqrt(pow((i-57.4), 2) + pow((j-64), 2));
                factor17 = 2.0 * (factor17 - radius) / if_th;
                double factor18 = sqrt(pow((i-72.2), 2) + pow((j-64), 2));
                factor18 = 2.0 * (factor18 - radius) / if_th;
                double factor19 = sqrt(pow((i-87), 2) + pow((j-64), 2));
                factor19 = 2.0 * (factor19 - radius) / if_th;
                //fourth row
                double factor20 = sqrt(pow((i-5), 2) + pow((j-50), 2));
                factor20 = 2.0 * (factor20 - radius) / if_th;
                double factor21 = sqrt(pow((i-20), 2) + pow((j-50), 2));
                factor21 = 2.0 * (factor21 - radius) / if_th;
                double factor22 = sqrt(pow((i-35), 2) + pow((j-50), 2));
                factor22 = 2.0 * (factor22 - radius) / if_th;
                double factor23 = sqrt(pow((i-50), 2) + pow((j-50), 2));
                factor23 = 2.0 * (factor23 - radius) / if_th;
                double factor24 = sqrt(pow((i-65), 2) + pow((j-50), 2));
                factor24 = 2.0 * (factor24 - radius) / if_th;
                double factor25 = sqrt(pow((i-80), 2) + pow((j-50), 2));
                factor25 = 2.0 * (factor25 - radius) / if_th;
                double factor26 = sqrt(pow((i-94), 2) + pow((j-50), 2));
                factor26 = 2.0 * (factor26 - radius) / if_th;
                //second part
                //fifth row
                double factor27 = sqrt(pow((i-13), 2) + pow((j-35), 2));
                factor27 = 2.0 * (factor27 - radius) / if_th;
                double factor28 = sqrt(pow((i-27.8), 2) + pow((j-35), 2));
                factor28 = 2.0 * (factor28 - radius) / if_th;
                double factor29 = sqrt(pow((i-42.6), 2) + pow((j-35), 2));
                factor29 = 2.0 * (factor29 - radius) / if_th;
                double factor30 = sqrt(pow((i-57.4), 2) + pow((j-35), 2));
                factor30 = 2.0 * (factor30 - radius) / if_th;
                double factor31 = sqrt(pow((i-72.2), 2) + pow((j-35), 2));
                factor31 = 2.0 * (factor31 - radius) / if_th;
                double factor32 = sqrt(pow((i-87), 2) + pow((j-35), 2));
                factor32 = 2.0 * (factor32 - radius) / if_th;
                //sixth row
                double factor33 = sqrt(pow((i-21), 2) + pow((j-21), 2));
                factor33 = 2.0 * (factor33 - radius) / if_th;
                double factor34 = sqrt(pow((i-35.5), 2) + pow((j-21), 2));
                factor34 = 2.0 * (factor34 - radius) / if_th;
                double factor35 = sqrt(pow((i-50), 2) + pow((j-21), 2));
                factor35 = 2.0 * (factor35 - radius) / if_th;
                double factor36 = sqrt(pow((i-64.5), 2) + pow((j-21), 2));
                factor36 = 2.0 * (factor36 - radius) / if_th;
                double factor37 = sqrt(pow((i-79), 2) + pow((j-21), 2));
                factor37 = 2.0 * (factor37 - radius) / if_th;
                //seventh row
                double factor38 = sqrt(pow((i-30), 2) + pow((j-7), 2));
                factor38 = 2.0 * (factor38 - radius) / if_th;
                double factor39 = sqrt(pow((i-43), 2) + pow((j-7), 2));
                factor39 = 2.0 * (factor39 - radius) / if_th;
                double factor40 = sqrt(pow((i-56), 2) + pow((j-7), 2));
                factor40 = 2.0 * (factor40 - radius) / if_th;
                double factor41 = sqrt(pow((i-69), 2) + pow((j-7), 2));
                factor41 = 2.0 * (factor41 - radius) / if_th;
                rho[i][j] = temp+factor1 * tanh(factor5)*tanh(factor6)*tanh(factor7)*tanh(factor8)*
                tanh(factor9)*tanh(factor10)*tanh(factor11)*tanh(factor12)*tanh(factor13)
                *tanh(factor14)*tanh(factor15)*tanh(factor16)*tanh(factor17)*tanh(factor18)*tanh(factor19)
                *tanh(factor20)*tanh(factor21)*tanh(factor22)*tanh(factor23)*tanh(factor24)*tanh(factor25)*tanh(factor26)
                *tanh(factor27)*tanh(factor28)*tanh(factor29)*tanh(factor30)*tanh(factor31)*tanh(factor32)
                *tanh(factor33)*tanh(factor34)*tanh(factor35)*tanh(factor36)*tanh(factor37)
                *tanh(factor38)*tanh(factor39)*tanh(factor40)*tanh(factor41);
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
            ComputeEquilibrium(i,j,ux,uy,rho,Fx,Fy,cx,cy,w,feq);
            f[k][i][j]=feq[k];
        }
    }
}
/*initial condition--------------*/
}
void Non_Equilibrium_Bounce_Back(double ***f,double **ux,double **uy,double **Fx,double **Fy,double **rho)
{
//west
    for (int j=1;j<my-1;j++)
    {
        //set boundary density the same as liquid density
        rho[0][j]=rho_l;
        uy[0][j]=-Fy[0][j]*dt/(2*rho[0][j]);
        ux[0][j]=1.0-(f[0][0][j]+f[2][0][j]+f[4][0][j]+2.0*(f[3][0][j]+f[6][0][j]+f[7][0][j]))/rho[0][j];
        f[1][0][j]=f[3][0][j]+2.0/3.0*rho[0][j]*(ux[0][j]+Fx[0][j]/(omega*rho[0][j]));
        f[5][0][j]=f[7][0][j]-1.0/2.0*(f[2][0][j]-f[4][0][j])+1.0/6.0*rho[0][j]*ux[0][j]-Fy[0][j]*dt/4-1/3*Fx[0][j]/(omega*rho[0][j]);
        f[8][0][j]=f[6][0][j]-1.0/2.0*(f[2][0][j]-f[4][0][j])+1.0/6.0*rho[0][j]*ux[0][j]+Fy[0][j]*dt/4-1/3*Fx[0][j]/(omega*rho[0][j]);
    }
//east
for (int j=1;j<my-1;j++)
    {
        rho[mx-1][j]=rho_l;
        uy[mx-1][j]=-Fy[mx-1][j]*dt/(2*rho[mx-1][j]);
        ux[mx-1][j]=(f[0][mx-1][j]+f[2][mx-1][j]+f[4][mx-1][j]+2.0*(f[1][mx-1][j]+f[5][mx-1][j]+f[8][mx-1][j]))/rho[mx-1][j]-1;
        f[3][mx-1][j]=f[1][mx-1][j]-2.0/3.0*rho[mx-1][j]*(ux[mx-1][j]+Fx[mx-1][j]/(omega*rho[mx-1][j]));
        f[6][mx-1][j]=f[8][mx-1][j]-1.0/2.0*(f[2][mx-1][j]-f[4][mx-1][j])-1.0/6.0*rho[mx-1][j]*ux[mx-1][j]-1/4*Fy[mx-1][j]*dt+1/3*Fx[mx-1][j]/omega;
        f[7][mx-1][j]=f[5][mx-1][j]+1.0/2.0*(f[2][mx-1][j]-f[4][mx-1][j])-1.0/6.0*rho[mx-1][j]*ux[mx-1][j]+1/4*Fy[mx-1][j]*dt+1/3*Fx[mx-1][j]/omega;
    }
//north
for (int i=1;i<mx-1;i++)
    {
        rho[i][my-1]=rho_l;
        ux[i][my-1]=-Fx[i][my-1]*dt/(2*rho[i][my-1]);
        uy[i][my-1]=(f[0][i][my-1]+f[1][i][my-1]+f[3][i][my-1]+2.0*(f[2][i][my-1]+f[5][i][my-1]+f[6][i][my-1]))/rho[i][my-1]-1;
        f[4][i][my-1]=f[2][i][my-1]-2.0/3.0*rho[i][my-1]*(uy[i][my-1]+Fy[i][my-1]/(omega*rho[i][my-1]));
        f[7][i][my-1]=f[5][i][my-1]+1.0/2.0*(f[1][i][my-1]-f[3][i][my-1])-1.0/6.0*rho[i][my-1]*uy[i][my-1]+1/3*Fy[i][my-1]/omega+1/4*Fx[i][my-1]*dt;
        f[8][i][my-1]=f[6][i][my-1]+1.0/2.0*(f[3][i][my-1]-f[1][i][my-1])-1.0/6.0*rho[i][my-1]*uy[i][my-1]+1/3*Fy[i][my-1]/omega-1/4*Fx[i][my-1]*dt;
    }
//south
for (int i=1;i<mx-1;i++)
    {
        rho[i][0]=rho_l;
        ux[i][0]=-Fx[i][0]*dt/(2*rho[i][0]);
        uy[i][0]=1-(f[0][i][0]+f[1][i][0]+f[3][i][0]+2.0*(f[4][i][0]+f[7][i][0]+f[8][i][0]))/rho[i][0];
        f[2][i][0]=f[4][i][0]+2.0/3.0*rho[i][0]*(uy[i][0]+Fy[i][0]/(omega*rho[i][0]));
        f[5][i][0]=f[7][i][0]-1.0/2.0*(f[1][i][0]-f[3][i][0])+1.0/6.0*rho[i][0]*uy[i][0]-1/4*Fx[i][0]*dt-1/3*Fy[i][0]/omega;
        f[6][i][0]=f[8][i][0]+1.0/2.0*(f[1][i][0]-f[3][i][0])+1.0/6.0*rho[i][0]*uy[i][0]+1/4*Fx[i][0]*dt-1/3*Fy[i][0]/omega;
    }
//corner
// Bottom left (inlet)
rho[0][0]=rho_l;
// rho[0][1]=rho_l;
// rho[1][1]=rho_l;
// rho[2][1]=rho_l;
// rho[3][1]=rho_l;
ux[0][0]=0.0;
uy[0][0]=0.0;
// std::cout<<rho[0][0]<<std::endl;
// std::cout<<rho[0][1]<<std::endl;
f[1][0][0] =f[3][0][0];
f[2][0][0] =f[4][0][0];
f[5][0][0] =f[7][0][0];
f[6][0][0]=1/2*(rho[0][0]-(f[0][0][0]+f[1][0][0]+f[2][0][0]+f[3][0][0]+f[4][0][0]+f[5][0][0]+f[7][0][0]));
f[8][0][0]=f[6][0][0];

//Bottom right
int i=mx-1; int j=0;
rho[i][j]=rho_l;
ux[i][j]=0.0;
uy[i][j]=0.0;
f[2][i][j]=f[4][i][j];
f[3][i][j] =f[1][i][j];
f[6][i][j]=f[8][i][j];
f[7][i][j] =1/2*(rho[i][j]-(f[0][i][j]+f[1][i][j]+f[2][i][j]+f[3][i][j]+f[4][i][j]+f[6][i][j]+f[8][i][j]));
f[5][i][j] =f[7][i][j];
//Top left
i=0; j=my-1;
rho[i][j]=rho_l;
ux[i][j]=0.0;
uy[i][j]=0.0;
f[1][i][j] =f[3][i][j];
f[4][i][j] =f[2][i][j];
f[8][i][j] =f[6][i][j];
f[7][i][j] =1/2*(rho[i][j]-(f[0][i][j]+f[1][i][j]+f[2][i][j]+f[3][i][j]+f[4][i][j]+f[6][i][j]+f[8][i][j]));
f[5][i][j] =f[7][i][j];
//Top right
i=mx-1; j=my-1;
rho[i][j]=rho_l;
ux[i][j]=0.0;
uy[i][j]=0.0;
f[3][i][j] =f[1][i][j];
f[4][i][j] =f[2][i][j];
f[7][i][j] =f[5][i][j];
f[6][i][j]=1/2*(rho[i][j]-(f[0][i][j]+f[1][i][j]+f[2][i][j]+f[3][i][j]+f[4][i][j]+f[5][i][j]+f[7][i][j]));
f[8][i][j]=f[6][i][j];
}
void VerifyLapalaceLaw(double **press,double **rho,double *x)
{
//calculate pressure according to different eos
for (int j = 0; j < my; j++) {
    for (int i = 0; i < mx; i++) {
        switch (eos)
        {
        case 2:
        {
            press[i][j]=rho[i][j]*R*T*(1+B*rho[i][j]/4+(B*rho[i][j]/4)*(B*rho[i][j]/4)-pow(B*rho[i][j]/4,3.0))/pow((1-B*rho[i][j]/4),3.0)-A*rho[i][j]*rho[i][j];
            break;
        }
        default:
        {
            press[i][j] = 0.;
            press[i][j] += rho[i][j] / 3.;
            press[i][j] += (g * psi(rho[i][j]) * psi(rho[i][j]) / 6.);
            break;
        }
        }
    }
}
//calculate deltapressure
double rho_gas = 0.;
double rho_liq = 0.;
double press_gas = 0.;
double press_liq = 0.;

// static int iteration = 0;
// iteration++;
// Average liquid pressure in a square box with 4 lattice nodes
for (int i = 0; i < 4; i++) 
{
    for (int j = 0; j < 4; j++) {
        // if(iteration==2) std::cout << "   " << i << ", " << j << ", " << rho[i][j] << std::endl;
        rho_liq += rho[i][j];
        press_liq += press[i][j];
    }
}

    double rad;
  // Average gas pressure in a square box with 4 lattice nodes
  //Originally it is i=mx/2-2,i<mx/2+2;i++ j = my / 2 - 2; j < my / 2 + 2; j++
  //but there will be high density in the center with high mesh resolution so I change it to what it looks like now.
  for (int i = mx / 2 - 10; i < mx / 2 - 8; i++) {
    for (int j = my / 2 - 10; j < my / 2 - 8; j++) {
      rho_gas += rho[i][j];
      press_gas += press[i][j];
    }
  }

rho_gas /= 16.;
rho_liq /= 16.;
press_gas /= 16.;
press_liq /= 16.;

double delta_press = press_gas - press_liq;
double rho_av = 0.5 * (rho_gas + rho_liq);//average density of gas and liquid 
const int y = my / 2;
// std::cout << rho_liq << ", " << rho_gas << " -> ";
//calculate radius
    for (int i = 1; i < mx / 2; i++) {
        if (rho[i][y] < rho_av) {
            const double drho = rho[i-1][y] - rho[i][y];
            const double dx = (rho_av - rho[i][y]) / drho;
            rad = (lx / 2. - x[i]) + dx;
            break;
        }
    }

    //std::cout << iteration << ": " << rad << ", " << delta_press << ", " << rho_av << std::endl; 
    //iteration += 1;

    //write data for radius,pressure and surface tension 
    std::ofstream out;
    std::string outputfile="pressure_"+to_string_with_precision(lx)+"_"+std::to_string(mx)+"_"+to_string_with_precision(rho_l)+".dat";
    out.open(outputfile,std::ios::app);
        out << std::scientific << std::setprecision(5) << std::setw(15) << delta_press;
        out << std::scientific << std::setprecision(5) << std::setw(15) << 1./rad;
        out << std::scientific << std::setprecision(5) << std::setw(15) << delta_press*rad;
        out << std::endl;
    out.close();
}
void Periodic_Boundary_Condition(double ***f)
{
//west
for (int j = 0; j < my; j++)
{
    f[1][0][j]=f[1][mx-1][j];
    f[5][0][j]=f[5][mx-1][j];
    f[8][0][j]=f[8][mx-1][j];
}
//east
for (int j = 0; j < my; j++)
{
    f[3][mx-1][j]=f[3][0][j];
    f[6][mx-1][j]=f[6][0][j];
    f[7][mx-1][j]=f[7][0][j];
}
//north
for (int i = 0; i < mx; i++)
{
    f[4][i][my-1]=f[4][i][0];
    f[7][i][my-1]=f[7][i][0];
    f[8][i][my-1]=f[8][i][0];
}
//south
for (int i = 0; i < mx; i++)
{
    f[2][i][0]=f[2][i][my-1];
    f[5][i][0]=f[5][i][my-1];
    f[6][i][0]=f[6][i][my-1];
}
}
int main()
{
//main loop
Readinitialfile(setupfile,&n,&mx,&my,&omega,&mstep,
&cas,&bc,&eos,&R,&A,&B,&T,&lx,&ly,&c,
&rho_l,&rho_g,&radius,&freq,&rho0,&g,&if_th);
int *cx = new int[n]; 
int *cy = new int[n];
double *w = new double[n]; 
double *x = new double[mx]; 
double *y = new double[my];
double *feq = new double[n];

double **rho = new double *[mx];
double **Fx = new double *[mx];
double **Fy = new double *[mx];
double **press = new double *[mx];
double **ux = new double *[mx];
double **uy = new double *[mx];
for(int s=0;s<mx;s++)
{
    rho[s] = new double [my];
    Fx[s] = new double [my];
    Fy[s] = new double [my];
    press[s] = new double [my];
    ux[s] = new double [my];
    uy[s] = new double [my];
}
double ***f = new double**[n];
for(int s=0;s<n;s++)
{
    f[s]= new double *[mx];
    for(int ss=0;ss<mx;ss++)
    {
        f[s][ss]=new double [my];
    }
}
initialize(cx,cy,w,x,y,rho,ux,uy,f,feq,Fx,Fy,if_th);
int time=0;
std::string filename = std::string("animation") +std::to_string(time)+std::string(".dat");
result(filename,time,x,y,ux,uy,rho,press,f,feq);
for (time=0;time<mstep+1;time=time+1)
{
    ComputeDesnity(f,rho);
    ComputeSCForce(cx,cy,rho,w,Fx,Fy);//only novelty which contribute to equilibrium velocity and macroscopic velocity 
    ComputeVelocity(f,ux,uy,cx,cy,rho);
    Collision(f,ux,uy,rho,Fx,Fy,cx,cy,w,feq);
    Streaming(f);
    switch (bc)
    {
    case 1:
        Periodic_Boundary_Condition(f);
        break;
    case 2: 
        Non_Equilibrium_Bounce_Back(f,ux,uy,Fx,Fy,rho);
        break;
    }
    VerifyLapalaceLaw(press,rho,x);//output radius,pressure and surface tension
    std::string filename = std::string("animation") +std::to_string(time)+std::string(".dat");
    if(time%freq==0) result(filename,time,x,y,ux,uy,rho,press,f,feq);
    if(time%freq==0)
    {
        std::cout << "Iteration: " << time << " / " << mstep << " (" << 100.0*double(time)/double(mstep) << " %)" <<  std::endl;
    }
}//main loop end
    //delete memory
    for(int s=0;s<mx;s++)
    {
        delete [] rho[s];
        delete [] Fx[s];
        delete [] Fy[s];
        delete [] press[s];
        delete [] ux[s];
        delete [] uy[s];
    }
    delete [] rho;
    delete [] Fx;
    delete [] Fy;
    delete [] press;
    delete [] ux;
    delete [] uy;

    for(int s=0; s<n; s++)
    {
        for(int ss=0; ss<mx; ss++)
        {
            delete [] f[s][ss];
        }
        delete [] f[s];
    }
    delete [] f;

    delete [] cx;
    delete [] cy;
    delete [] w;
    delete [] x;
    delete [] y;
    delete [] feq;
return 0;
}