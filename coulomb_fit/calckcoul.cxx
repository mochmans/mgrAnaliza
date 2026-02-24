#include <iostream>
#include <math.h>
#include <TGraph.h>
#include <TMath.h>
#include <TFile.h>

using namespace std;

#define NSTEPS 200

long double ac;

struct _dcomplex 
{
  long double re;
  long double im;
};

typedef struct _dcomplex dcomplex;

dcomplex mult(dcomplex arga, dcomplex argb)
{
  dcomplex res;
  
  res.re = arga.re * argb.re - arga.im * argb.im;
  res.im = arga.re * argb.im + argb.re * arga.im;

  return res;
}

dcomplex mult(dcomplex arga, long double argb)
{
  dcomplex res;
  
  res.re = arga.re * argb;
  res.im = arga.im * argb;

  return res;
}

void disp(dcomplex arg)
{
  cout << arg.re;
  if (arg.im > 0) cout << "+";
  cout << arg.im << "i";
}

long double Gamow(long double arg)
{
  return (2.0*TMath::Pi()/(ac * arg)) * 1.0/(exp(2.0*TMath::Pi()/(ac * arg)) - 1.0);
}

TGraph *getgr(int nsteps, double *xy, long double *yk, const char *name)
{
  double *yki = (double *) malloc(sizeof(double *) * nsteps);
  
  for (int iter=0; iter<nsteps; iter++) {
    yki[iter] = yk[iter]; 
  }
  
  TGraph *grbuf = new TGraph(nsteps, xy, yki);
  grbuf->SetName(name);

  free(yk);

  return grbuf;
}

int main(int argc, char **argv) 
{
  int ptype = 0;
  float radinp = 5.0;

  if (argc < 2) {
    cout << "Program calculating the squared Coulomb wave-function" << endl
	 << "integrated over a gaussian source with equal sizes in all dimensions" << endl
	 << "for pion, kaon or antiproton-xi system" << endl
	 <<"Usage" << endl << "calckcoul <particle type> <radius>" << endl
	 << "where particle type is 0 for pions 1 for kaons and 2 for antiproton-xi" << endl
	 << "and radius is single-particle Rinv size" << endl;
    ptype = 0;
    radinp = 5.0;
  }
  else {
    ptype = atoi(argv[1]);
    radinp = atof(argv[2]);
  }
  
  if (ptype==1) {
    ac = 109.55 / 0.197327;
    cout << "Calculation for kaons" << endl;
  }
  if (ptype==2) {
    // calucate bohr radius ac:
    // ac = (h*c)/(mu*c*alpha)
    // mu - reduced mass
    // 1 GeV^(-1) = 0.197 fm
    
    ac = 49.28 / 0.197327;
    cout << "Calculation for antiproton-xi" << endl;
  }
  if(ptype == 3){
    ac = 44.9 / 0.197327;
    cout << "Calculation for antiproton-Omega" << endl;
  }
  else {
    ac = 387.5 / 0.197327;
    cout << "Calculation for pions" << endl;
  }
  cout << "Single particle radius is " << radinp << endl;

  long double radius = TMath::Sqrt(2.0) * radinp / 0.197327; 
  long double eta, ksi, rad, radx, rady, radz, coskr;

  int nstep = 150;

  long double **compre;
  long double **compim;

  long double **compaccrep;
  long double **compaccimp;
  long double **compaccrem;
  long double **compaccimm;

  TGraph **grcompre;
  TGraph **grcompim;
  TGraph **grcompaccre;
  TGraph **grcompaccim;
  TGraph **grall;
  TGraph **gralla;

  long double **funallt;
  long double **funalla;

  double xy[NSTEPS];

  grcompre = (TGraph **) malloc(sizeof(TGraph *) * nstep);
  grcompim = (TGraph **) malloc(sizeof(TGraph *) * nstep);
  grcompaccre = (TGraph **) malloc(sizeof(TGraph *) * nstep);
  grcompaccim = (TGraph **) malloc(sizeof(TGraph *) * nstep);
  grall  = (TGraph **) malloc(sizeof(TGraph *) * nstep);
  gralla = (TGraph **) malloc(sizeof(TGraph *) * nstep);
  
  compre = (long double **) malloc(sizeof(long double *) * nstep);
  compim = (long double **)malloc(sizeof(long double *) * nstep);
  compaccrep = (long double **)malloc(sizeof(long double *) * nstep);
  compaccimp = (long double **)malloc(sizeof(long double *) * nstep);
  compaccrem = (long double **)malloc(sizeof(long double *) * nstep);
  compaccimm = (long double **)malloc(sizeof(long double *) * nstep);

  funallt = (long double **) malloc(sizeof(long double *) * nstep);
  funalla = (long double **) malloc(sizeof(long double *) * nstep);

  for (int iter=0; iter<nstep; iter++) {
    compre[iter] =    (long double *) malloc(sizeof(long double) * NSTEPS);
    compim[iter] =    (long double *) malloc(sizeof(long double) * NSTEPS);
    compaccrep[iter] = (long double *) malloc(sizeof(long double) * NSTEPS);
    compaccimp[iter] = (long double *) malloc(sizeof(long double) * NSTEPS);
    compaccrem[iter] = (long double *) malloc(sizeof(long double) * NSTEPS);
    compaccimm[iter] = (long double *) malloc(sizeof(long double) * NSTEPS);
    funallt[iter] =   (long double *) malloc(sizeof(long double) * NSTEPS);
    funalla[iter] =   (long double *) malloc(sizeof(long double) * NSTEPS);
  }

  dcomplex alfa, zetp, zetm;
  long double stepwidth = 0.4/NSTEPS;

  for (int kiter = 0; kiter<NSTEPS; kiter++) {
    long double kstar = kiter*stepwidth+stepwidth/2.0;
    xy[kiter] = kstar;

    if (!((kiter+1) % 10))
      cout << "Iter: " << kiter << endl;

    eta = 1.0/(kstar * ac);
    alfa.re = 0.0;
    alfa.im = -eta;

    dcomplex fcomp, scompp, scompm;
    long double tcomp;
    dcomplex summ, sump;
    dcomplex fcmult;
    long double wgt;
    long double norm = 0.0;

    long double pcp, mcm, pcm;

    for (int istep=0; istep<nstep; istep++) {
      funallt[istep][kiter] = 0.0;
      funalla[istep][kiter] = 0.0;
    }
    
    for (int rxiter=0; rxiter<50; rxiter++) {
      for (int ryiter=0; ryiter<50; ryiter++) {
	for (int rziter=0; rziter<50; rziter++) {
	  radx = (rxiter - 25 + 0.5) * radius * 3.0 / 25.0;
	  rady = (ryiter - 25 + 0.5) * radius * 3.0 / 25.0;
	  radz = (rziter - 25 + 0.5) * radius * 3.0 / 25.0;
	  rad  = sqrt(radx*radx + rady*rady + radz*radz);
	  coskr = radx/rad;
	  wgt = exp(-(radx*radx + rady*rady + radz*radz)/(2.0*radius*radius));

	  ksi = kstar*rad*(1+coskr);
    
	  zetp.re = 0.0;
	  zetp.im = ksi;
      
	  ksi = kstar*rad*(1-coskr);
	  
	  if (ksi<15.0) {
	    zetm.re = 0.0;
	    zetm.im = ksi;
	    
	    fcomp.re = 1.0;
	    fcomp.im = 0.0;
	    scompp.re = 1.0; 
	    scompp.im = 0.0;
	    scompm.re = 1.0; 
	    scompm.im = 0.0;
	    tcomp = 1.0;
	    
	    for (int istep=0; istep<nstep; istep++) {
	      sump = mult(fcomp, scompp);
	      sump = mult(sump, 1.0/(tcomp*tcomp));
	      
	      summ = mult(fcomp, scompm);
	      summ = mult(summ, 1.0/(tcomp*tcomp));
	      
	      compre[istep][kiter] = sump.re;
	      compim[istep][kiter] = sump.im;
	      
	      if (istep == 0) {
		compaccrep[istep][kiter] = sump.re;
		compaccimp[istep][kiter] = sump.im;
		compaccrem[istep][kiter] = summ.re;
		compaccimm[istep][kiter] = summ.im;
	      }
	      else {
		compaccrep[istep][kiter] = compaccrep[istep-1][kiter] + sump.re;
		compaccimp[istep][kiter] = compaccimp[istep-1][kiter] + sump.im;
		compaccrem[istep][kiter] = compaccrem[istep-1][kiter] + summ.re;
		compaccimm[istep][kiter] = compaccimm[istep-1][kiter] + summ.im;
	      }
	      
	      fcmult.re = alfa.re + istep;
	      fcmult.im = alfa.im;
	      
	      fcomp = mult(fcomp, fcmult);
	      scompp = mult(scompp, zetp);
	      scompm = mult(scompm, zetm);
	      tcomp *= (istep+1);
	      funallt[istep][kiter] += 
		(compaccrep[istep][kiter]*compaccrep[istep][kiter] + 
		 compaccimp[istep][kiter]*compaccimp[istep][kiter]) * 
		wgt;
	      
	      pcp = (compaccrep[istep][kiter]*compaccrep[istep][kiter] + 
		     compaccimp[istep][kiter]*compaccimp[istep][kiter]);
	      mcm = (compaccrem[istep][kiter]*compaccrem[istep][kiter] + 
		     compaccimm[istep][kiter]*compaccimm[istep][kiter]);
	      pcm = (compaccrep[istep][kiter]*compaccrem[istep][kiter] + 
		     compaccimp[istep][kiter]*compaccimm[istep][kiter]);
	      
	      funalla[istep][kiter] += 
		((pcp + mcm - 2.0*pcm)/(1+cos(2*kstar*radx)) + 2.0*pcm) * wgt * 0.5;  
	    }
	    norm += wgt;
	  }
	  else {
	    funalla[nstep-1][kiter] += (1.0 - 1.0/(rad*(1+coskr)*ac*kstar*kstar))*wgt/Gamow(kstar);
	    funallt[nstep-1][kiter] += (1.0 - 1.0/(rad*(1+coskr)*ac*kstar*kstar))*wgt/Gamow(kstar);
	    norm += wgt;
	  }
	}
      }
    }
    for (int istep=0; istep<nstep; istep++) {
      funallt[istep][kiter] /= norm;
      funalla[istep][kiter] /= norm;
    }
  }

  char gname[100];


  TFile *ofile = new TFile(Form("ffcomp_%f.root", radinp),"RECREATE");
  ofile->cd();

  for (int iter=0; iter<nstep; iter++) {
    sprintf(gname, "compre%i", iter);
    grcompre[iter] = getgr(NSTEPS, xy, compre[iter], gname);
    sprintf(gname, "compim%i", iter);
    grcompim[iter] = getgr(NSTEPS, xy, compim[iter], gname);
    sprintf(gname, "compaccre%i", iter);
    grcompaccre[iter] = getgr(NSTEPS, xy, compaccrep[iter], gname);
    sprintf(gname, "compaccim%i", iter);
    grcompaccim[iter] = getgr(NSTEPS, xy, compaccimp[iter], gname);
    sprintf(gname, "Full%i", iter);
    grall[iter] = getgr(NSTEPS, xy, funallt[iter], gname);
    sprintf(gname, "Full%ia", iter);
    gralla[iter] = getgr(NSTEPS, xy, funalla[iter], gname);

    grcompre[iter]->Write();
    grcompim[iter]->Write();
    grcompaccre[iter]->Write();
    grcompaccim[iter]->Write();
    grall[iter]->Write();
  }
  
  TGraph *grkcoul = new TGraph(*(grall[nstep-1]));
  grkcoul->SetName("KCoulomb");
  grkcoul->SetTitle("Coulomb wave-function squared over the source");
  for (int ibin=0; ibin<grkcoul->GetN(); ibin++) {
    grkcoul->GetY()[ibin] = grall[nstep-1]->GetY()[ibin]*Gamow(grkcoul->GetX()[ibin]);
  }
  grkcoul->Write();

  TGraph *grkcoula = new TGraph(*(grall[nstep-1]));
  grkcoula->SetName("KCoulombAdvanced");
  grkcoula->SetTitle("Pair wave function divided by Quantum Statistics");
  for (int ibin=0; ibin<grkcoul->GetN(); ibin++) {
    grkcoula->GetY()[ibin] = gralla[nstep-1]->GetY()[ibin]*Gamow(grkcoul->GetX()[ibin]);
  }
  grkcoula->Write();

  delete ofile;
}
