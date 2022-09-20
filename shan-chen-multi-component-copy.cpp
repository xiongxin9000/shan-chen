#include <iostream>
#include <cmath>
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
using namespace std;

// Lattice parameters
const int nx = 200; // number of nodes along x-axis
const int ny = 200; // number of nodes along y-axis
const int nsteps = 10; // number of time steps
const int noutput = 1; // data output interval (data written to "/data")
const int nfluids = 1; // number of fluid components (choose 1 or 2)
const double tauA = 1; // relaxation time for fluid A
const double tauB = 1; // relaxation time for fluid B (will only be used if nfluids = 2)

// Shan-Chen parameters
const double gA = -4.7; // self-interaction strength of fluid A
const double gB = 0.; // self-interaction strength of fluid B (will only be used if nfluids = 2)
const double gAB = 6.; // interaction strength of fluids A and B (will only be used if nfluids = 2)
const double rho0 = 1.; // reference density for pseudopotential (usually set to 1)
const double rho_l = 1.95; // initial liquid density
const double rho_g = 0.15; // initial gas density
const double radius = 25.; // initial droplet radius

// Fixed parameters; DO NOT CHANGE
const int npop = 9; // number of populations
const int cx[] = {0, 1, 0, -1, 0, 1, -1, -1, 1}; // x-components of lattice vectors
const int cy[] = {0, 0, 1,  0,-1, 1,  1, -1,-1}; // y-components of lattice vectors
const double w[]={4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};

// Arrays
double g[nfluids][nfluids]; // interaction strengths
double rho[nfluids][nx][ny]; // density
double ux[nx][ny]; // x-component of fluid velocity
double uy[nx][ny]; // y-component of fluid velocity
double Fx[nfluids][nx][ny]; // x-component of Shan-Chen force
double Fy[nfluids][nx][ny]; // y-component of Shan-Chen force
double feq[nfluids][npop]; // equilibrium populations
double forcing[nfluids][npop]; // force populations
double tau[nfluids]; // relaxation times
double f1[nfluids][npop][nx][ny]; // populations (old)
double f2[nfluids][npop][nx][ny]; // populations (new)
double press[nx][ny];
// Compute density everywhere.
void computeDensity(double pop[nfluids][npop][nx][ny]) {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j <ny; j++)
    {
        for (int s = 0; s < nfluids; s++) {
        rho[s][i][j] = pop[s][0][i][j] + pop[s][1][i][j] + pop[s][2][i][j] + pop[s][3][i][j] + pop[s][4][i][j] + pop[s][5][i][j] + pop[s][6][i][j] + pop[s][7][i][j] + pop[s][8][i][j];
        }
    }
  }
  return;
}

// Compute pseudopotential for given density value.
double psi(double dens) {
  return rho0 * (1. - exp(-dens / rho0));
}

// Compute Shan-Chen forces.
void ComputeSCForce() {
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      for (int s = 0; s < nfluids; s++) {
        Fx[s][i][j] = 0.;
        Fy[s][i][j] = 0.;
        for (int ss = 0; ss < nfluids; ss++) {
          double fxtemp = 0.;
          double fytemp = 0.;

          for (int k = 1; k < 9; k++) {
            const int i2 = (i + cx[k] + nx) % nx;
            const int j2 = (j + cy[k] + ny) % ny;
            const double psinb = psi(rho[ss][i2][j2]);
            fxtemp += w[k] * cx[k] * psinb;
            fytemp += w[k] * cy[k] * psinb;
          }

          const double psiloc = psi(rho[s][i][j]);
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

// Calculate total density at given point.
double computeTotalDensity(int i,int j) {
  double dens;
  
  if (nfluids == 1) {
    dens = rho[0][i][j];
  }
  else if (nfluids == 2) {
    dens = rho[0][i][j] + rho[1][i][j];
  }

  return dens;
}
// Calculate barycentric velocity everywhere.
void computeVelocity(double pop[nfluids][npop][nx][ny]) {
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      ux[i][j] = 0.;
      uy[i][j] = 0.;
      for (int s = 0; s < nfluids; s++) {
        ux[i][j] += (pop[s][1][i][j] - pop[s][3][i][j] + pop[s][5][i][j] - pop[s][6][i][j] - pop[s][7][i][j] + pop[s][8][i][j]);
        uy[i][j] += (pop[s][2][i][j] - pop[s][4][i][j] + pop[s][5][i][j] + pop[s][6][i][j] - pop[s][7][i][j] - pop[s][8][i][j]);
        ux[i][j] += (0.5 * Fx[s][i][j]);
        uy[i][j] += (0.5 * Fy[s][i][j]);
      }
      const double dens = computeTotalDensity(i,j);
      ux[i][j] /= dens;
      uy[i][j] /= dens;
    }
  }

  return;
}

// Compute equilibrium distributions for fluid s at point k.
void equilibrium(int s, int i,int j) {
  const double dens = rho[s][i][j];
  const double vx = ux[i][j];
  const double vy = uy[i][j];
  const double usq = vx * vx + vy * vy;
  feq[s][0] = w[0] * dens * (1. - 1.5 * usq);
  feq[s][1] = w[1] * dens * (1. + 3. * vx + 4.5 * vx * vx - 1.5 * usq);
  feq[s][2] = w[2] * dens * (1. + 3. * vy + 4.5 * vy * vy - 1.5 * usq);
  feq[s][3] = w[3] * dens * (1. - 3. * vx + 4.5 * vx * vx - 1.5 * usq);
  feq[s][4] = w[4] * dens * (1. - 3. * vy + 4.5 * vy * vy - 1.5 * usq);
  feq[s][5] = w[5] * dens * (1. + 3. * ( vx + vy) + 4.5 * ( vx + vy) * ( vx + vy) - 1.5 * usq);
  feq[s][6] = w[6] * dens * (1. + 3. * (-vx + vy) + 4.5 * (-vx + vy) * (-vx + vy) - 1.5 * usq);
  feq[s][7] = w[7] * dens * (1. + 3. * (-vx - vy) + 4.5 * ( vx + vy) * ( vx + vy) - 1.5 * usq);
  feq[s][8] = w[8] * dens * (1. + 3. * ( vx - vy) + 4.5 * ( vx - vy) * ( vx - vy) - 1.5 * usq);

  return;
}

// Initialise simulation
void initialisation(double f1[nfluids][npop][nx][ny], double f2[nfluids][npop][nx][ny]) {
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

  // Initialise 1-component system.
  // In this case a liquid droplet is created in a gas.
  if (nfluids == 1) {
    for(int i=0;i<nx;i++)
    {
        for (int j=0;j<ny;j++)
        {
              const int k = j * nx + i;
              if ((k / nx - ny / 2.) * (k / nx - ny / 2.) + (k % nx - nx / 2.) * (k % nx - nx / 2.) <= radius * radius)
              //if((i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)<radius*radius)
              {
                rho[0][i][j]=rho_l;
              }
            else
            {
            rho[0][i][j]=rho_g;
            }
        }
    }
  }
  // Initialise populations.
  for (int s = 0; s < nfluids; s++) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j <ny; j++)
        {
            // Set initial velocity field (zero velocity everywhere).
            ux[i][j] = uy[i][j] = 0.;

            // Compute equilibrium distributions
            equilibrium(s,i,j);

            // Set populations to equilibrium
            for (int k = 0; k < npop; k++) {
                f1[s][k][i][j] = feq[s][k];
                f2[s][k][i][j] = feq[s][k];
            }
        }
    }
  }
  
  return;
}

// LBM: collide and propagate (push).
void push(double f1[nfluids][npop][nx][ny], double f2[nfluids][npop][nx][ny]) {
  for (int s = 0; s < nfluids; s++) {
    const double omega = 1. / tau[s];

    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        // Compute Guo's forcing terms
        for (int k = 0; k < npop; k++) {
          forcing[s][k] = w[k] * (1. - 0.5 * omega) * ((3. * (cx[k] - ux[i][j]) + 9. * cx[k] * (cx[k] * ux[i][j] + cy[k] * uy[i][j])) * Fx[s][i][j] + (3. * (cy[k] - uy[i][j]) + 9. * cy[k] * (cx[k] * ux[i][j] + cy[k] * uy[i][j])) * Fy[s][i][j]);
        }

        // Compute equilibrium distributions
        equilibrium(s, i,j);

        // Collide and propagate
        for (int k = 0; k < npop; k++)
        {
          int i2 = (i + cx[k] + nx) % nx;
          int j2 = (j + cy[k] + ny) % ny;
          f2[s][k][i2][j2] = f1[s][k][i][j] * (1. - omega) + feq[s][k] * omega + forcing[s][k];
        }
      }
    }
  }
  
  return;
}
void result(std::string filename,int time)
{

std::ofstream out(filename);
  out << "TITLE=\"LBM output\"" << std::endl;
  out << "VARIABLES = \"X\", \"Y\", \"Ux\",\"Uy\",\"density\",\"f0\",\"f1\",\"f2\",\"f3\",\"f4\",\"f5\",\"f6\",\"f7\",\"f8\"" << std::endl;
  out << "ZONE T = \"fluid\", I=" << nx << ", J=" << ny << ", F=POINT" << std::endl;
  //out << "SOLUTIONTIME="<< time << std::endl;
  for (unsigned j = 0; j < ny; ++j)
    for (unsigned i = 0; i < nx; ++i)
    {
        out << std::scientific << std::setprecision(5) << std::setw(15) << i;
        out << std::scientific << std::setprecision(5) << std::setw(15) << j;
        out << std::scientific << std::setprecision(5) << std::setw(15) << ux[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << uy[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << rho[0][i][j]+rho[1][i][j];
      for(int k=0;k<9;k++)
      {
          out << std::scientific << std::setprecision(5) << std::setw(15) << f1[0][k][i][j]+f1[1][k][i][j];
      }
      out << std::endl;
    }
  out.close();

}
void VerifyLapalaceLaw()
{
//calculate pressure
for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
        press[i][j] = 0.;
        press[i][j] += rho[0][i][j] / 3.;
        press[i][j] += (g[0][0] * psi(rho[0][i][j]) * psi(rho[0][i][j]) / 6.);
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
        press_gas += press[i][j];
    }
}

  // Average liquid pressure in a square box with 4 lattice nodes
  for (int i = nx / 2 - 2; i < nx / 2 + 2; i++) {
    for (int j = ny / 2 - 2; j < ny / 2 + 2; j++) {
      rho_liq += rho[0][i][j];
      press_liq += press[i][j];
    }
  }

  rho_gas /= 16.;
  rho_liq /= 16.;
  press_gas /= 16.;
  press_liq /= 16.;
  
  double delta_press = press_liq - press_gas;
  double rho_av = 0.5 * (rho_gas + rho_liq);
  const int y = ny / 2;
//calculate radius
    double rad;
    for (int i = 0; i < nx / 2; i++) {
        if (rho[0][i][y] > rho_av) {
            const double drho = rho[0][i][y] - rho[0][i-1][y];
            const double dx = (rho_av - rho[0][i-1][y]) / drho;
            rad = (nx / 2. - i) + (1. - dx);
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
// Main function.
int main(int argc, char** argv) {

// Initialise domain
initialisation(f1, f2);
int step=0;
std::string filename = std::string("animation") +std::to_string(step)+std::string(".dat");
result(filename,step);

  // Main loop
  for (int step = 1; step <= nsteps; ++step)
   {
    // Calculate density field.
    computeDensity(f1);
    
    // Calculate SC forces.
    ComputeSCForce();//????
    
    // Calculate barycentric velocity.
    computeVelocity(f1);
    
    // Collide and propagate.
    push(f1, f2);

    // Swap old and new populations.
    std::swap(f1, f2);
    
    // Write data to disk/console
    std::string filename = std::string("animation") +std::to_string(step)+std::string(".dat");
    if (step % noutput == 0) 
    {
      cout << "Running time step " << step << endl;    
      result(filename,step);
    }
  VerifyLapalaceLaw();
  }
  return 0;
}
