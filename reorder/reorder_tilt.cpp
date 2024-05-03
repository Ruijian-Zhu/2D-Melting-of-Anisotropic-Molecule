// C++ code for arranging COM file calculated by LAMMPS into a LAMMPS style trajectory file
// Modified for triclinic box
// Written by Ruijian Zhu (ITP-CAS) in Nov 2022

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

const int frame = 5;
const int molecule = 2496;      // Number of COMs
const int sample_gap = 2000;    // Sample interval, ensuring the timestep in new traj in agreement with the original one
std::string fileHead ("com_1.2");      // file name of COM of each snapshot calculated by LAMMPS (fileHead.step.dat)
std::string dumpFile ("dump_npt_h.1.2.lammpstrj");
std::string writeFile ("dump_1.2_reorder.lammpstrj");
const int dumpHeadLine = 5;
const int comHeadLine = 4;
const int atom = molecule * 6;  // 6 is the case for hexagon, can be modified for other cases

using namespace std;
ifstream infile;   
ofstream write;

struct coordinate
{
    double x;
    double y;
};

string filename(int a)
{
    stringstream ss;
    ss << fileHead << "." << a << ".dat";
    return ss.str();
}

// Note that the COM calculated by LAMMPS does not consider PBC (if the original traj file has ix iy iz)
coordinate pbc(double xlo, double xhi, double ylo, double yhi, double xy, coordinate atom)
{
    coordinate basis1;
    basis1.x = xhi - xlo;
    basis1.y = 0;
    coordinate basis2;
    basis2.x = xy;
    basis2.y = yhi - ylo;
    while(atom.y > yhi)
    {
        atom.x -= xy;
        atom.y -= basis2.y;
    }
    while(atom.y < ylo)
    {
        atom.x += xy;
        atom.y += basis2.y;
    }
    double xleft, xright;
    double dx;
    dx = (atom.y-ylo)*(xy)/(yhi-ylo);
    xleft = xlo + dx;
    xright = xhi + dx;
    while(atom.x > xright)
    {
        atom.x -= basis1.x;
    }
    while(atom.x < xleft)
    {
        atom.x += basis1.x;
    }
    return atom;
}

int main()
{
    double box_size[(frame+1)*9]={};
    double xlo, xhi, xy, ylo, yhi;
    infile.open(dumpFile,ios::in);
	while (infile.bad())
	{
		cerr << "cannot open the file" << endl;
		return 0;
	}
    int datalen = 0;
    int j = 0;

    // Read the box information from original traj file
    // Note that, for example, run 40000 and sample per 2000, then you get 21 snapshot in the original output traj file
    // So we need to give up the first snapshot (timestep = 0)
    for(j = 0; j < frame+1; j++)
    {
        for(int i = 0; i < dumpHeadLine; )
        {
            if(infile.get()=='\n') i++; 
        }
        for(int i = 0; i < 9; i++)
        {
            infile >> box_size[datalen++];
        }
        for(int i = 0; i < 2+atom; )
        {
            if(infile.get()=='\n') i++; 
        }
    }
    infile.close();

    // Read in coordinate of COMS and write to the new traj file spontaneously
    write.open(writeFile, ios::app); 
    double rdinate[4*molecule];
    coordinate particle[molecule];
    int molIndex;
    for(j=0;j<frame;j++)
    {
        datalen = 0;
        int index = (j+1) * sample_gap;
        string comFileName = filename(index);
        infile.open(comFileName,ios::in);
        while (infile.bad())
	    {
		    cerr << "cannot open the file" << endl;
		    return 0;
	    }
        for(int i = 0; i < comHeadLine; )
        {
            if(infile.get()=='\n') i++; 
        }
        for(int i = 0; i < molecule*4; i++)
        {
            infile >> rdinate[datalen++];
        }
        // in LAMMPS, xlo = x_left-bound (at y = 0) + xy if xy < 0, xhi = x_right-bound (at y = 0) +xy if xy > 0
        // Here we set xlo = x_lef-bound (at y = 0) and xhi = x_right-bound (at y = 0)
        // So we have to modify this before PBC calculation
        if(box_size[j*9+11] > 0)
        {
            xlo = box_size[j*9+9];
            xhi = box_size[j*9+10] - box_size[j*9+11];
            ylo = box_size[j*9+12];
            yhi = box_size[j*9+13];
            xy = box_size[j*9+11];
        }
        else
        {
            xlo = box_size[j*9+9] - box_size[j*9+11];
            xhi = box_size[j*9+10];
            ylo = box_size[j*9+12];
            yhi = box_size[j*9+13];
            xy = box_size[j*9+11];
        }
        for(int i = 0; i < molecule; i++)
        {
            particle[i].x = rdinate[4*i+1];
            particle[i].y = rdinate[4*i+2];
            particle[i] = pbc(xlo,xhi,ylo,yhi,xy,particle[i]);
        }

        // Write to the new traj file in LAMMPS style
        write<<"ITEM: TIMESTEP"<<endl;
	    write << index <<endl;
	    write<<"ITEM: NUMBER OF ATOMS"<<endl;
	    write<< molecule <<endl;
	    write<<"ITEM: BOX BOUNDS xy xz yz pp pp pp"<<endl;
	    write<< box_size[j*9+9] << ' ' << box_size[j*9+10] << ' ' << box_size[j*9+11] <<endl;
	    write<< box_size[j*9+12] << ' ' << box_size[j*9+13] << ' ' << box_size[j*9+14] <<endl;
	    write<< box_size[j*9+15] << ' ' << box_size[j*9+16] << ' ' << box_size[j*9+17] <<endl;
	    write<<"ITEM: ATOMS id mol type x y z"<<endl;
        for(molIndex=0;molIndex<molecule;molIndex++)
        {
            write << rdinate[molIndex*4] << ' ' << rdinate[molIndex*4] << ' ' << "1" <<\
            ' ' << particle[molIndex].x << ' ' << particle[molIndex].y << ' ' << "0" << endl;
        }
        infile.close();
    }
    return 0;
}
