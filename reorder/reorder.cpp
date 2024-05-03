// C++ code for arranging COM file calculated by LAMMPS into a LAMMPS style trajectory file
// Written by Ruijian Zhu (ITP-CAS) in Feb 2023
// 2496 molecules, 2000 snapshots, 1 CPU, requires about 2 min

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

const int frame = 5;
const int molecule = 2496;          // Number of COMs
const int sample_gap = 2000;        // Sample interval, ensuring the timestep in new traj in agreement with the original one
std::string fileHead ("com_1");   // file name of COM of each snapshot calculated by LAMMPS (fileHead.step.dat)
std::string dumpFile ("dump_npt_o.1.lammpstrj");
std::string writeFile ("dump_1_reorder.lammpstrj");
const int dumpHeadLine = 5;
const int comHeadLine = 4;
const int atom = molecule * 5;      // 5 is the case for pentagon, can be modified for other cases

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
// Periodic boundary condition for rectangular box
coordinate pbc(double xlo, double xhi, double ylo, double yhi, coordinate atom)
{
    coordinate basis1;
    basis1.x = xhi - xlo;
    basis1.y = 0;
    coordinate basis2;
    basis2.x = 0;
    basis2.y = yhi - ylo;
    while(atom.y > yhi)
    {
        atom.x -= 0;
        atom.y -= basis2.y;
    }
    while(atom.y < ylo)
    {
        atom.x += 0;
        atom.y += basis2.y;
    }
    while(atom.x > xhi)
    {
        atom.x -= basis1.x;
    }
    while(atom.x < xlo)
    {
        atom.x += basis1.x;
    }
    return atom;
}

int main()
{
    double box_size[(frame+1)*6]={};
    double xlo, xhi, ylo, yhi;
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
        for(int i = 0; i < 6; i++)
        {
            infile >> box_size[datalen++];
        }
        for(int i = 0; i < 2+atom; )
        {
            if(infile.get()=='\n') i++; 
        }
    }
    infile.close();
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
        xlo = box_size[j*6+6];
        xhi = box_size[j*6+7];
        ylo = box_size[j*6+8];
        yhi = box_size[j*6+9];
        for(int i = 0; i < molecule; i++)
        {
            particle[i].x = rdinate[4*i+1];
            particle[i].y = rdinate[4*i+2];
            particle[i] = pbc(xlo,xhi,ylo,yhi,particle[i]);
        }
        write<<"ITEM: TIMESTEP"<<endl;
	    write << index <<endl;
	    write<<"ITEM: NUMBER OF ATOMS"<<endl;
	    write<< molecule <<endl;
	    write<<"ITEM: BOX BOUNDS pp pp pp"<<endl;
	    write<< box_size[j*6+6] << ' ' << box_size[j*6+7] <<endl;
	    write<< box_size[j*6+8] << ' ' << box_size[j*6+9] <<endl;
	    write<< box_size[j*6+10] << ' ' << box_size[j*6+11] <<endl;
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
