// C++ code for determine the crystal axis (scanning at fixed r)
// Note that tan\theta can be calculated directly (y/x) and what we need is also tan\theta, so we never convert it to theta
// For triclinic box
// Written by Ruijian Zhu (ITP-CAS) in Oct 2022
// 2496 molecules, 2000 snapshots, 1 CPU, requires about 30 min

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include <string>
#include<sstream>

const int unitCell = 2496;
const int frame = 5;
const int dumpHeadLine1 = 5;
const int dumpHeadLine2 = 2;
const double scale = 0.3;
const int interval = 300;
const double r0 = 72.375;
const double error = 4e-1;
std::string fileHead ("theta_1.2_all.txt"); 
std::string dumpFile ("dump_1.2_reorder.lammpstrj");

using namespace std;
ifstream infile;   
ofstream write;

struct coordinate
{
    double x;
    double y;
};

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

coordinate triclincPBC(double xy, double yhi, double ylo, double xhi, double xlo, coordinate pR, coordinate pC)
{
    coordinate basis1;
    basis1.x = xhi - xlo;
    basis1.y = 0;
    coordinate basis2;
    basis2.x = xy;
    basis2.y = yhi - ylo;
    coordinate PCR[9];
    PCR[0] = pC;
    PCR[1] = movevector(pC,basis1,true);
    PCR[2] = movevector(pC,basis2,true);
    PCR[3] = movevector(pC,basis1,false);
    PCR[4] = movevector(pC,basis2,false);
    PCR[5] = movevector(PCR[2],basis1,true);
    PCR[6] = movevector(PCR[2],basis1,false);
    PCR[7] = movevector(PCR[4],basis1,true);
    PCR[8] = movevector(PCR[4],basis1,false);
    double dist = distance(pR,PCR[0]);
    double test;
    int index = 0;
    for(int i=1; i < 9; i++)
    {
        test = distance(pR,PCR[i]);
        if(test<dist)
        {
            dist = test;
            index = i;
        }
    }
    coordinate delta;
    delta.x = PCR[index].x - pR.x;
    delta.y = PCR[index].y - pR.y;
    return delta;
}

bool judge(double xy, double yhi, double ylo, double xhi, double xlo, coordinate atom1, coordinate atom2)
{
    coordinate deltar;
    deltar = triclincPBC(xy,yhi,ylo,xhi,xlo,atom1,atom2);
    double r = deltar.x;
    if(fabs(r-r0) < error)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int assign(double dist)
{
    double a;
    a = (dist + scale) * interval / (2 * scale);
    int ass = int(a);
    return ass;
}

double assignInverse(int i)
{
    double d;
    d = i * (2 * scale/interval) + 0.5 * (2 * scale/interval) - scale;
    return d;
}

int main()
{
    double box[9];
    coordinate atom[unitCell];
    double id;
    double type;
    double mol;
    double z;
    double count[interval] = {};
    double xlo, xhi, xy, ylo, yhi;
    double tanTheta;
    int in;
    infile.open(dumpFile,ios::in);
    for(int k = 0; k < frame; k++)
    {
        for(int i = 0; i < dumpHeadLine1; )
        {
            if(infile.get()=='\n') i++; 
        }
        for(int i = 0; i < 9; i++)
        {
            infile >> box[i];
        }
        if(box[2] > 0)
        {
            xlo = box[0];
            xhi = box[1] - box[2];
            ylo = box[3];
            yhi = box[4];
            xy = box[2];
        }
        else
        {
            xlo = box[0] - box[2];
            xhi = box[1];
            ylo = box[3];
            yhi = box[4];
            xy = box[2];
        }
        for(int i = 0; i < dumpHeadLine2; )
        {
            if(infile.get()=='\n') i++; 
        }
        for(int i = 0; i < unitCell; i++)
        {
            infile >> id >> z >> type >> atom[i].x >> atom[i].y >> z;
        }
        for(int i = 0; i < 1; )
        {
            if(infile.get()=='\n') i++; 
        }
        for(int i = 0; i < unitCell; i++)
        {
            for(int j = i + 1; j < unitCell; j++)
            {
                bool countin = judge(xy,yhi,ylo,xhi,xlo,atom[i],atom[j]);
                coordinate dr = triclincPBC(xy,yhi,ylo,xhi,xlo,atom[i],atom[j]);
                tanTheta = dr.y / dr.x;
                if(countin && -scale <= tanTheta && tanTheta <= scale)
                {
                    //cout << "here" << endl;
                    in = assign(tanTheta);
                    count[in] ++;
                }
            }
        }
    }
    infile.close();
    write.open(fileHead,ios::app);
    for(int i = 0; i < interval; i++)
    {
        write << assignInverse(i) << '\t' << count[i]/frame << endl;
    }
    write.close();
    return 0;
}
