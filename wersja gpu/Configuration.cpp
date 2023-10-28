#pragma once
#include "Configuration.h"

//---------------------------------------------------------------------------------------------------------------------------
void geoinp(int nnmax, int nwmax, int ntmax, int &nnod, float *xnod, float *ynod, float *znod, \
	        int &nw, int *iwn, float *wrad, int *nsw, int &nt, int *iwt)
{
	fstream geofile;
	string fname;
	string linia;
		
	while (1)
	{
		fname = "dipol.geo";
		cout << " Nazwa pliku geo: " << fname << endl;
		//cin >> fname;
		geofile.open(fname, ios::in);		
		if (geofile.good() == false) cout << " Nie udalo sie otworzyc pliku!" << endl;
		else break;
	}
	
	// wczytanie i wyswietlenie komentarzy w pliku geo
	getline(geofile, linia);	
	getline(geofile, linia);
	
	//sekcja wezlow
	geofile >> nnod;
	if (nnod > nnmax)
	{
		cout << " Przekroczona liczba wezlow w pliku geo" << endl;
		cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
		getchar();getchar();
		exit(0);
	}
	for (int i = 1; i <= nnod; i++)
	{
		geofile >> xnod[i]; 
		geofile >> ynod[i]; 
		geofile >> znod[i]; 
	}
	// sekcja drucikow
	geofile >> nw;
	if (nw > nwmax)
	{
		cout << " Przekroczona liczba drucikow w pliku geo (nw > nwmax)" << endl;
		cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
		getchar(); getchar();
		exit(0);
	}
	for (int i = 1; i <= nw; i++)
	{
		geofile >> iwn[i];
		geofile >> iwn[i + nw];
		geofile >> wrad[i];
		geofile >> nsw[i];
	}
	// sekcja platkow
	geofile >> nt;
	if (nt > ntmax)
	{
		cout << " Przekroczona liczba platkow w pliku geo (nt > ntmax)" << endl;
		cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
		getchar(); getchar();
		exit(0);
	}
	for (int i = 0; i < nt; i++)
	{
		geofile >> iwt[i];
		geofile >> iwt[i + nt];
		geofile >> iwt[i + 2 * nt];
	}


	geofile.close();
}

//---------------------------------------------------------------------------------------------------------------------------
void config(int &npnt, int nnod, float *xnod, float *ynod, float *znod, float *x, float *y, float *z, int nw, \
	int *iwn, float *wrad, int *nsw, int *ise, int *isp, int *ips, int *ipc, int *iwe, float *rad, \
	int &npls, int &nplsw, int &nseg, int nsmax, int nwmax, int ienv)

	
{
	configwire(npnt, nnod, xnod, ynod, znod, x, y, z, nw, iwn, wrad, nsw, ise, isp, ips, ipc, iwe, rad, nplsw, nseg, nsmax, nwmax, ienv);

	//int npntw = 0;
	//int nplsbp=0;
	//configbody(npntw, nnod, xnod, ynod, znod, x, y, z, nt, iwt, ipt, itp, iti, nplsbp, npbmax);
	//nplsb = nplsbp;
	npls = nplsw; //nplsw + nplsb + nplsj;
	cout << endl;
	cout << " ##################################################################\n";
	cout << " +-----------------------------------------------------------------+\n";
	cout << " | Inspect file configwire.txt for a list of current pulses        |\n";
	cout << " | assciated with wires, and coordinates of their nodes            |\n";
	cout << " +-----------------------------------------------------------------+\n";
	cout << " ##################################################################\n";
    //cout << " +-----------------------------------------------------------------+\n";
	//cout << " | Inspect file configbody.txt for a list of current pulses        |\n";
	//cout << " | asociated with patches, and coordinates of their nodes and edges|\n";
    //cout << " +-----------------------------------------------------------------+\n";
	//cout << " ##################################################################\n";
	//cout << " +-----------------------------------------------------------------+\n";
	//cout << " | Inspect file configjun.txt for a list of current pulses         |\n";
	//cout << " | associated with junctions, and coordinates of their nodes       |\n";
	//cout << " +-----------------------------------------------------------------+\n";
	cout << endl;
}


//---------------------------------------------------------------------------------------------------------------------------
void configwire(int &npnt, int nnod, float *xnod, float *ynod, float *znod, float *x, float *y, float *z, int nw, \
	int *iwn, float *wrad, int *nsw, int *ise, int *isp, int *ips, int *ipc, int *iwe, float *rad, \
	int &nplsw, int &nseg, int nsmax, int nwmax, int ienv)
	
{
	fstream wirefile;
	wirefile.open("configwire.txt", ios::out);
	if (nw == 0)
	{
		wirefile << endl << " -== NO WIRES IN GEOMETRY ==- " << endl;
		wirefile.close();
		return;
	}
	// brak sprawdzenia czy nw > nwmax -> to ju¿ jest w geoinp
	nseg = 0;
	nplsw = 0;
	int npntw = 0;
	
	//sprawdzenie, czy przy opcji PEC ground nie na ¿adnego wêz³a o wspó³rzêdnej z < 0
	if (ienv == 2)
	{
		for (int i = 1; i <= nnod; i++)
		{
			if (znod[i] < 0.)
			{
				cout << endl << "ERROR: z-coordinate of any node can not be negative" << endl;
				cout << "         when environment is perfectly conducting ground!" << endl;
				cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
				getchar(); getchar();
				exit(0);
			}
		}
	}

	int nvnod = 0; // liczba virtualnych wezlow
	int	nvw = 0; // liczba wirtualnych drucikow
	if (ienv == 2)
	{
		for (int i = 1; i <= nnod; i++)
		{
			if (abs(znod[i]) < 1.e-6)	// oryginalnie tu by³o znod[i]=0. tutaj zmienione bo nie wiem jak z porównywaniem liczb float do 0!
			{
				nvnod = nvnod + 1;
				nvw = nvw + 1;
			}
		}
	}
	int nadnod = nvnod;
	nvnod = nvnod + nnod;
	nvw = nvw + nw;

	float *xnodv, *ynodv, *znodv, *wradv;
	int *iwnv, *nswv;
	xnodv = new float[nvnod+1]();
	ynodv = new float[nvnod+1]();
	znodv = new float[nvnod+1]();
	wradv = new float[nvw+1]();
	iwnv = new int[2 * nvw+1]();
	nswv = new int[nvw+1]();

	int ivn = 0;
	int ivw = 0;
	float dx, dy, dz, dl;
	if (ienv == 2)
	{
		for (int i = 1; i <= nnod; i++)
		{
			if (abs(znod[i]) < 1.e-6)	// oryginalnie tu by³o znod[i]=0. tutaj zmienione bo nie wiem jak z porównywaniem liczb float do 0!
			{
				for (int j = 1; j <= nw; j++)
				{
					if ((iwn[j] == i) || (iwn[j + nw] == i))
					{
						ivn = ivn + 1;
						ivw = ivw + 1;
						dx = xnod[iwn[j + nw]] - xnod[iwn[j]];
						dy = ynod[iwn[j + nw]] - ynod[iwn[j]];
						dz = znod[iwn[j + nw]] - znod[iwn[j]];
						dl = sqrt(dx*dx + dy * dy + dz * dz) / nsw[j];						
						xnodv[ivn] = xnod[i] + dl;
						ynodv[ivn] = ynod[i];
						znodv[ivn] = znod[i];
						iwnv[ivw] = i + nadnod;
						iwnv[ivw + nvw] = ivn;
						wradv[ivw] = wrad[j];
						nswv[ivw] = 1;												
					}
				}
			}
		}
	}
	

	// przepisanie (dopelnienie) tablic wirtualnych
	int ivnpi;
	for (int i = 1; i <= nnod; i++)
	{
		ivnpi = ivn + i;
		xnodv[ivnpi] = xnod[i];
		ynodv[ivnpi] = ynod[i];
		znodv[ivnpi] = znod[i];		
	}

	int ivwpi;
	for (int i = 1; i <= nw; i++)
	{
		ivwpi = ivw + i;
		iwnv[ivwpi] = iwn[i] + nadnod;
		iwnv[ivwpi + nvw] = iwn[i + nw] + nadnod;
		wradv[ivwpi] = wrad[i];
		nswv[ivwpi] = nsw[i];		
	}		
	
	int i = 0; 
	int np = 0;  
	int ns;
	bool conect;	

	for (int m = 1; m <= nvw; m++)
	{
		conect = false;
		ns = nswv[m];
		dx = xnodv[iwnv[m + nvw]] - xnodv[iwnv[m]];
		dy = ynodv[iwnv[m + nvw]] - ynodv[iwnv[m]];
		dz = znodv[iwnv[m + nvw]] - znodv[iwnv[m]];		
		
		if (ns > 1)
		{
			dx = dx / ns;
			dy = dy / ns;
			dz = dz / ns;
		}

		if (m != 1)
		{
			for (int n = 1; n <= m - 1; n++)
			{
				if (iwnv[m] == iwnv[n]) // connection 1
				{
					conect = true;
					nseg = nseg + 1;
					if (nseg > nsmax)
					{
						cout << endl << " Number of segments exceeds nsmax" << endl;
						cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
						getchar(); getchar();
						exit(0);
					}
					ise[nseg] = ise[iwe[n] + nsmax];
					ise[nseg + nsmax] = ise[iwe[n]];
					rad[nseg] = wradv[n];
					np = np + 1;
					if (np > (nsmax - nwmax))   // npwmax=nsmax-nwmax
					{
						cout << endl << " Number of basis functions on wires exceeds nwmax" << endl;
						cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
						getchar(); getchar();
						exit(0);
					}
					ipc[np] = -ise[nseg + nsmax];
					isp[nseg] = 0;
					isp[nseg + nsmax] = np;
					ips[np] = nseg;					
 					break;
				}
				if (iwnv[m] == iwnv[n+nvw]) // connection 2
				{
					conect = true;
					nseg = nseg + 1;
					if (nseg > nsmax)
					{
						cout << endl << " Number of segments exceeds nsmax" << endl;
						cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
						getchar(); getchar();
						exit(0);
					}
					ise[nseg] = ise[iwe[n+nsmax]];
					ise[nseg + nsmax] = ise[iwe[n+nsmax]+nsmax];
					rad[nseg] = wradv[n];
					np = np + 1;
					if (np > (nsmax - nwmax))   // npwmax=nsmax-nwmax
					{
						cout << endl << " Number of basis functions on wires exceeds nwmax" << endl;
						cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
						getchar(); getchar();
						exit(0);
					}
					ipc[np] = -ise[nseg + nsmax];
					isp[nseg] = 0;
					isp[nseg + nsmax] = np;
					ips[np] = nseg;					
					break;
				}				
			}			
		}
		
		if (!conect)
		{
			i = i + 1;
			x[i] = xnodv[iwnv[m]];
			y[i] = ynodv[iwnv[m]];
			z[i] = znodv[iwnv[m]];						
		}

		for (int j = 1; j <= ns; j++)
		{
			i = i + 1;
			if (j == ns)
			{				
				x[i] = xnodv[iwnv[m + nvw]];
				y[i] = ynodv[iwnv[m + nvw]];
				z[i] = znodv[iwnv[m + nvw]];				
			}
			else
			{
				x[i] = xnodv[iwnv[m]] + j * dx;
				y[i] = ynodv[iwnv[m]] + j * dy;
				z[i] = znodv[iwnv[m]] + j * dz;											
			}
			
			nseg = nseg + 1;
			if (nseg > nsmax)
			{
				cout << endl << " Number of segments exceeds nsmax" << endl;
				cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
				getchar(); getchar();
				exit(0);
			}
			rad[nseg] = wradv[m];			
			if (j == 1)
			{
				iwe[m] = nseg;
				if (conect) iwe[m] = nseg - 1;				
			}
			ise[nseg] = i - 1;  
			if (conect) ise[nseg] = ise[nseg - 1 + nsmax]; 
			ise[nseg + nsmax] = i;   
			isp[nseg] = 0;
			isp[nseg + nsmax] = 0;			
		}

		int mfrst = iwe[m];

		int mlast = mfrst + ns - 2;
		if (conect)
		{
			mlast = mlast + 1;
			np = np - 1;
		}
				
		for (int n = mfrst; n <= mlast; n++)
		{
			np = np + 1;
			if (np > (nsmax - nwmax))   // npwmax=nsmax-nwmax
			{
				cout << endl << " Number of basis functions on wires exceeds nwmax" << endl;
				cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
				getchar(); getchar();
				exit(0);
			}
			isp[n + nsmax] = np;
			isp[n + 1] = np;
			ips[np] = n;
			ips[np + nsmax] = n + 1;
			if (!(conect  && (n==mfrst)))
			{
				ipc[np] = ise[n + nsmax];
			}			
		}

		iwe[m + nsmax] = nseg;  //Wiersz 478

		if (m != 1)
		{
			for (int n = 1; n <= m-1; n++)
			{
				if (iwnv[m+nvw] == iwnv[n]) // connection 3
				{
					nseg = nseg + 1;
					if (nseg > nsmax)
					{
						cout << endl << " Number of segments exceeds nsmax" << endl;
						cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
						getchar(); getchar();
						exit(0);
					}
					i = i - 1;
					ise[nseg] = ise[iwe[n]];
					ise[nseg + nsmax] = ise[iwe[n] + nsmax];
					rad[nseg] = wradv[n];
					np = np + 1;
					if (np > (nsmax - nwmax))   // npwmax=nsmax-nwmax
					{
						cout << endl << " Number of basis functions on wires exceeds nwmax" << endl;
						cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
						getchar(); getchar();
						exit(0);
					}
					ipc[np] = -ise[nseg];
					isp[iwe[m + nsmax] + nsmax] = np;
					isp[nseg] = np;
					isp[nseg + nsmax] = 0;
					ips[np] = iwe[m + nsmax];
					ips[np + nsmax] = nseg;
					ise[iwe[m + nsmax] + nsmax] = abs(ipc[np]);
					iwe[m + nsmax] = nseg;					
					break;
				}
				if (iwnv[m + nvw] == iwnv[n + nvw]) // connection 4
				{
					nseg = nseg + 1;
					if (nseg > nsmax)
					{
						cout << endl << " Number of segments exceeds nsmax" << endl;
						cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
						getchar(); getchar();
						exit(0);
					}
					i = i - 1;
					ise[nseg] = ise[iwe[n + nsmax] + nsmax];
					ise[nseg + nsmax] = ise[iwe[n + nsmax]];
					rad[nseg] = wradv[n];
					np = np + 1;
					if (np > (nsmax - nwmax))   // npwmax=nsmax-nwmax
					{
						cout << endl << " Number of basis functions on wires exceeds nwmax" << endl;
						cout << " Nacisnij dowolny klawisz zeby zakonczyc program";
						getchar(); getchar();
						exit(0);
					}
					ipc[np] = -ise[nseg];
					isp[iwe[m + nsmax] + nsmax] = np;
					isp[nseg] = np;
					isp[nseg + nsmax] = 0;
					ips[np] = iwe[m + nsmax];
					ips[np + nsmax] = nseg;
					ise[iwe[m + nsmax] + nsmax] = abs(ipc[np]);
					iwe[m + nsmax] = nseg;					
					break;
				}				
			}
		}			
	}

	npntw = i;
	npnt = npntw;
	nplsw = np;	

	// zapisywanie wyników do pliku
	int iwrk, nb, ne;
	wirefile << endl << " Number of wires and basis functions: " << nw << '\t' << nplsw << endl;
	for (int ii = 1; ii <= nw; ii++)
	{
		wirefile << " Drucik nr " << ii << endl;
		iwrk = ii + ivw;
		nb = isp[iwe[iwrk] + nsmax];
		ne = isp[iwe[iwrk + nsmax]];
		for (int j = nb; j <= ne; j++)
		{
			if (j == 0)
			{
				wirefile << " no basis function(s)" << endl;
			}
			else
			{
				np = abs(ipc[j]);
				wirefile << j << '\t' << x[np] << '\t' << y[np] << '\t' << z[np] << endl;
			}
		}
	}
    
	//przeskalowanie tablic ¿eby index zaczyna³ siê od 0
	
	for (int i = 1; i <= nnod; i++)
	{
		xnod[i-1] = xnod[i];
		ynod[i-1] = ynod[i];
		znod[i-1] = znod[i];		
	}
	for (int i = 1; i <= nw; i++)
	{
		iwn[i-1]=iwn[i];
		iwn[i-1 + nw]=iwn[i+nw];
		wrad[i-1]=wrad[i];
		nsw[i-1]=nsw[i];
	}

	for (int i = 1; i <= nsmax; i++)
	{
		x[i-1] = x[i];
		y[i-1] = y[i];
		z[i-1] = z[i];
		ise[i-1] = ise[i];
		ise[i - 1 + nsmax] = ise[i + nsmax];
		isp[i-1] = isp[i];
		isp[i - 1 + nsmax] = isp[i + nsmax];
		iwe[i - 1] = iwe[i];
		iwe[i - 1 + nsmax] = iwe[i + nsmax];
		rad[i - 1] = rad[i];
		ips[i - 1] = ips[i];
		ips[i - 1 + nsmax] = ips[i + nsmax];
		ipc[i - 1] = ipc[i];
	}
	
	delete[] xnodv, ynodv, znodv, wradv;
	delete[] iwnv, nswv;

	wirefile.close();
}


//---------------------------------------------------------------------------------------------------------------------------
void tablice(int npnt, int nseg, int nsmax, float *x, float *y, float *z, int nplsw, int *ise, int *isp, int *ips, int *ipc)
{
	fstream tablicefile;
	tablicefile.open("tablice.txt", ios::out);

	//--------------------------------wyniki konfiguracji modelu
	tablicefile << " ------------------------------------------------" << endl;;
	tablicefile << " Coordinates of points: " << endl;
	tablicefile << " Number of nodes: " << npnt << endl;
	for (int i = 0; i < npnt; i++)
	{
		tablicefile << " Point: " << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << endl;
	}
	tablicefile << " ------------------------------------------------" << endl;
	tablicefile << " ---== WIRES ==---" << endl;
	tablicefile << " ipc: number of node for i-th basis function for wires" << endl;
	for (int i = 0; i < nplsw; i++)
	{
		tablicefile << " Basis function: " << i << "\t" << "Node: " << ipc[i] << endl;
	}
	tablicefile << " Number of 1 and 2 nodes of i-th segment" << endl;
	for (int i = 0; i < nseg; i++)
	{
		tablicefile << " segment: " << i << "\t" << "nodes: " << ise[i] << "\t" << ise[i + nsmax] << endl;
	}
	tablicefile << " Number of 1 and 2 segment of i-th basis function: " << endl;
	for (int i = 0; i < nplsw; i++)
	{
		tablicefile << " Basis function: " << i << "\t" << "segments: " << ips[i] << "\t" << ips[i + nsmax] << endl;
	}

	tablicefile.close();
}

/*void tablice(int nnod, int nw, int nt, int npls, int nplsb, int npbmax, float *x, float *y, float *z, \
	int *ipt, int *itp, int *iti)
{
	fstream tablicefile;
	tablicefile.open("tablice.txt", ios::out);

	tablicefile << " Number of nodes: " << nnod << endl;
	tablicefile << " Number of wires: " << nw << endl;
	tablicefile << " Number of pathes: " << nt << endl;
	tablicefile << " Number of basis functions associated with patches: " << nplsb << endl;
	tablicefile << " Total number of basis functions (unknowns): " << npls << endl;
	tablicefile << " ------------------------------------------------" << endl;
	tablicefile << " Coordinates of points : " << endl;
	tablicefile << " Number of points: \t" << nnod << endl;
	tablicefile.precision(8);
	tablicefile << scientific;
	for (int i = 0; i < nnod; i++)
	{
		tablicefile << " Point: \t" << i + 1 << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << endl;
	}
	tablicefile << " ------------------------------------------------" << endl;
	tablicefile << " ------------------------------------------------" << endl;
	tablicefile << " ---== PATCHES ==---" << endl;
	tablicefile << " ipt - Number of 1 and 2 patch of i-th basis function:" << endl;
	for (int i = 0; i < nplsb; i++)
	{
		tablicefile << " Basis function: " << i + 1 << "\t Patches: \t" << ipt[i] << "\t" << ipt[i + npbmax] << endl;
	}
	tablicefile << " itp - Number of basis function associated with 1, 2, 3 patch edge:" << endl;
	for (int i = 0; i < nt; i++)
	{
		tablicefile << "\t" << i + 1 << "\t" << itp[i] << "\t" << itp[i + nt] << "\t" << itp[i + 2 * nt] << endl;
	}
	tablicefile << " iti - Number of 1, 2, 3 node of i-th patch:" << endl;
	for (int i = 0; i < nt; i++)
	{
		tablicefile << "\t" << i + 1 << "\t" << iti[i] << "\t" << iti[i + nt] << "\t" << iti[i + 2 * nt] << endl;
	}

	tablicefile.close();
}*/