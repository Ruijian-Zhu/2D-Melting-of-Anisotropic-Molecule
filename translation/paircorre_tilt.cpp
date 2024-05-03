// C++ code for the calculation of translational function
// Note that the crystal axis orientation has to be determined at first
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
const double scale = 100;
const int interval = 400;
const double theta = 0.063;
const double error = 0.2;
std::string fileHead ("py_1.2_all.txt"); 
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
    return d;
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
    double deltay = r * (100*theta);
    double real = deltar.y*100;
    if(fabs(real-deltay) < error*100)
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
    a = dist * interval / scale;
    int ass = int(a);
    return ass;
}

double assignInverse(int i)
{
    double d;
    d = i * (scale/interval) + 0.5 * (scale/interval);
    return d;
}

int main()
{
    double box[9];
    coordinate atom[unitCell];
    double id;
    double mol;
    double type;
    double z;
    double count[interval] = {};
    double xlo, xhi, xy, ylo, yhi;
    double r;
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
        for(int i = 0; i < dumpHeadLine2; )
        {
            if(infile.get()=='\n') i++; 
        }
        for(int i = 0; i < unitCell; i++)
        {
            infile >> id >> mol >> type >> atom[i].x >> atom[i].y >> z;
        }
        for(int i = 0; i < 1; )
        {
            if(infile.get()=='\n') i++; 
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
        for(int i = 0; i < unitCell ; i++)
        {
            for(int j = i; j < unitCell; j++)
            {
                bool countin = judge(xy,yhi,ylo,xhi,xlo,atom[i],atom[j]);
                if(countin)
                {
                    coordinate dr = triclincPBC(xy,yhi,ylo,xhi,xlo,atom[i],atom[j]);
                    r = dr.x;
                    if (r < 0)
                    {
                        r *= -1;
                    }
                    if(r <= scale)
                    {
                        in = assign(r);
                        count[in]++;
                    }
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
