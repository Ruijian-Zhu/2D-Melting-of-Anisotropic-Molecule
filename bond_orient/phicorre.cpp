// C++ code for calculate bond-orientational correlation function along r direction
// Assuming typical lammpstrj output of BOOP as input (see example)
// Written by Ruijian Zhu (ITP-CAS) in April 2023

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include <string>
#include<sstream>

using namespace std;
ifstream infile1;  
ifstream infile2; 
ofstream write;

const int headLine = 4;
const int unitCell = 2496;          // number of unit-cell representative points
const int frame = 5;
const int sample_gap = 2000;        // used to calculate step in the file name below
const int dumpHeadLine = 9;
const int hexHeadLine1 = 5;
const int hexHeadLine2 = 2;
const double scale = 60;           // max r
const int interval = 400;           // number of bins
std::string fileHead ("hex_1");  // BOOP of each monomer each snapshot is named as filehead_step.lammpstrj
std::string dumpFile ("dump_1_reorder.lammpstrj"); // COM traj file name

struct coordinate
{
    double x;
    double y;
};

struct hexorder
{
    int id;
    coordinate complex;
};

string hexname(int a)
{
    stringstream ss;
    ss << fileHead << "_" << (a+1)*sample_gap << ".lammpstrj";
    return ss.str();
}

coordinate movevector(coordinate pR, coordinate m, bool judge)
{
    if(judge == true)
    {
        pR.x += m.x;
        pR.y += m.y;
    }
    else
    {
        pR.x -= m.x;
        pR.y -= m.y;
    }
    return pR;
}

double distance(coordinate p1, coordinate p2)
{
    coordinate delta;
    delta.x = p1.x - p2.x;
    delta.y = p1.y - p2.y;
    double d;
    d = (delta.x)*(delta.x) + (delta.y)*(delta.y);
    return sqrt(d);
}

// This function is now used for rectangular box but can be modified to triclinic box easily by adopting a non-zero basis2.x
// For large tilting box, it is difficult to judge the nearest image, so we compare all the 9 neighbors
// This is slightly time-consuming for rectangular box or slightly tilting box
double triclinicPBC(double yhi, double ylo, double xhi, double xlo, coordinate pR, coordinate pC)
{
    coordinate basis1;
    basis1.x = xhi - xlo;
    basis1.y = 0;
    coordinate basis2;
    basis2.x = 0;
    basis2.y = yhi - ylo;
    coordinate PCR[9];  // 9 neighboring images
    PCR[0] = pC;
    PCR[1] = movevector(pC,basis1,true);
    PCR[2] = movevector(pC,basis2,true);
    PCR[3] = movevector(pC,basis1,false);
    PCR[4] = movevector(pC,basis2,false);
    PCR[5] = movevector(PCR[1],basis2,true);
    PCR[6] = movevector(PCR[1],basis2,false);
    PCR[7] = movevector(PCR[3],basis2,true);
    PCR[8] = movevector(PCR[3],basis2,false);
    double dist = 100000;
    double test;
    for(int i=0; i < 9; i++)
    {
        test = distance(pR,PCR[i]);
        if(test<dist)
        {
            dist = test;
        }
    }
    return dist;
}

// assign to each bin
int assign(double dist)
{
    double a;
    a = dist * interval / scale;
    int ass = int(a);
    return ass;
}

// decide the middle value of each bin
double assignInverse(int i)
{
    double d;
    d = i * (scale/interval) + 0.5 * (scale/interval);
    return d;
}

int main()
{
    double box[6];
    coordinate atom[unitCell];  // To store the coordinate of each representative com
    hexorder hex[unitCell];     // To store the hexagonal order parameter of each com with its id together
    double d[unitCell*6];       // To store the data read in from dump file (id, type, x, y, z)
    coordinate intermediateCoordinate1;
    coordinate intermediateCoordinate2;
    coordinate intermediateHex1;
    coordinate intermediateHex2;
    int count[interval]={};        //Count how many numbers of pairs drop into each interval
    coordinate statistic[interval]; // the value on each interval
    for(int j=0; j<interval;j++)
    {
        statistic[j].x = 0;
        statistic[j].y = 0;
    }
    double twoDist; // just to calculate two particle distance
    infile1.open(dumpFile,ios::in);
    int i = 0;
    for(i = 0; i < frame; i++)
    {
        // We first read in the dump file, this file would not be closed until finish
        // skip all the head lines, directly go to coordinates
        for(int j = 0; j < dumpHeadLine; )
        {
            if(infile1.get()=='\n') j++; 
        }
        // store the data in d temprorally
        for(int j = 0; j < unitCell*6; j++)
        {
            infile1 >> d[j];
        }
        // gap the last '\n'
        for(int j = 0; j < 1; )
        {
            if(infile1.get()=='\n') j++; 
        }
        // assign the coordinates to atom array
        for(int j = 0; j < unitCell; j++)
        {
            atom[j].x = d[6*j+3];
            atom[j].y = d[6*j+4];
        }
        //
        // open the hexagonal order parameter file
        string hexFile = hexname(i);
        infile2.open(hexFile,ios::in);
        // skip the head lines to the box size
        for(int j = 0; j < hexHeadLine1; )
        {
            if(infile2.get()=='\n') j++; 
        }
        // read in the box size
        for(int j = 0; j < 6; j++)
        {
            infile2 >> box[j];
        }
        for(int j = 0; j < hexHeadLine2; )
        {
            if(infile2.get()=='\n') j++;
        }
        // read in hexorder of each atoms
        for(int j = 0; j < unitCell; j++)
        {
            //infile2 >> d[j];
            infile2 >> hex[j].id >> hex[j].complex.x >> hex[j].complex.y;
        }
        //Now we have to use a double-loop and calculate their correlation
        for(int j = 0; j < unitCell-1; j++)
        {
            intermediateCoordinate1 = atom[hex[j].id];
            intermediateHex1 = hex[j].complex;
            for(int k = j + 1; k < unitCell; k++)
            {
                intermediateCoordinate2 = atom[hex[k].id];
                intermediateHex2 = hex[k].complex;
                twoDist = triclinicPBC(box[3],box[2],box[1],box[0],intermediateCoordinate1, intermediateCoordinate2);
                if(twoDist < scale)
                {
                    statistic[assign(twoDist)].x += intermediateHex1.x * intermediateHex2.x + intermediateHex1.y\
                    * intermediateHex2.y;
                    statistic[assign(twoDist)].y += intermediateHex1.x * intermediateHex2.y - intermediateHex2.x\
                    * intermediateHex1.y;
                    count[assign(twoDist)]++;
                }
            }
        }
        infile2.close(); 
    }
    write.open("correlation.txt",ios::app);      // Output file name may be modified
    for(i = 0; i < interval; i++)
    {
        if(count[i] > 0)
        {
            statistic[i].x /= count[i];
            statistic[i].y /= count[i];
            write << assignInverse(i) << '\t' << sqrt((statistic[i].x * statistic[i].x) + (statistic[i].y\
            * statistic[i].y)) << endl;
        }
    }
    return 0;
}
