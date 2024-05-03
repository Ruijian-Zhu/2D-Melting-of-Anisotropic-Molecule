// C++ code for determine the crystal axis (scanning at fixed r)
// Note that tan\theta can be calculated directly (y/x) and what we need is also tan\theta, so we never convert it to theta
// Written by Ruijian Zhu (ITP-CAS) in Jan 2023
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
const double scale = 0.3;       // Range of tan\theta
const int interval = 300;
const double r0 = 64.125;       // The fixed r
const double error = 4e-1;      // Error allowed in r 
std::string fileHead ("theta_1_all.txt"); // Output file
std::string dumpFile ("dump_1_reorder.lammpstrj");

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
    d = sqrt(d);
    return d;
}

// Modified from triclinic box, so still decide the nearest image by comparing distance directly

coordinate triclincPBC(double yhi, double ylo, double xhi, double xlo, coordinate pR, coordinate pC)
{
    coordinate basis1;
    basis1.x = xhi - xlo;
    basis1.y = 0;
    coordinate basis2;
    basis2.x = 0;
    basis2.y = yhi - ylo;
    coordinate PCR[5];
    PCR[0] = pC;
    PCR[1] = movevector(pC,basis1,true);
    PCR[2] = movevector(pC,basis2,true);
    PCR[3] = movevector(pC,basis1,false);
    PCR[4] = movevector(pC,basis2,false);
    double dist = distance(pR,PCR[0]);
    double test;
    int index = 0;
    for(int i=1; i < 5; i++)
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

// Judge whether the distance between i and j is the fixed r0 (with a tolerance error)

bool judge(double yhi, double ylo, double xhi, double xlo, coordinate atom1, coordinate atom2)
{
    coordinate deltar;
    deltar = triclincPBC(yhi,ylo,xhi,xlo,atom1,atom2);
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
    a = (dist + scale) * interval / (2 * scale); // 2 comes from tha fact that tan\theta is in the range \pm scale
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
    double box[6];
    coordinate atom[unitCell];
    double id;    // since reorder traj file is ordered as atom index, so this is useless
    double type;
    double z;
    double mol;
    double count[interval] = {};
    double tanTheta;
    int in;
    infile.open(dumpFile,ios::in);
    for(int k = 0; k < frame; k++)
    {
        for(int i = 0; i < dumpHeadLine1; )
        {
            if(infile.get()=='\n') i++; 
        }
        for(int i = 0; i < 6; i++)
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
        for(int i = 0; i < unitCell; i++)
        {
            for(int j = i + 1; j < unitCell; j++)
            {
                bool countin = judge(box[3],box[2],box[1],box[0],atom[i],atom[j]);  // judge whether they are separated at r0
                coordinate dr = triclincPBC(box[3],box[2],box[1],box[0],atom[i],atom[j]);
                tanTheta = dr.y / dr.x;         // Calculate the tan\theta
                if(countin && -scale <= tanTheta && tanTheta <= scale)  // if at r0 && tan\theta in the range, assign it into appropriate bin
                {
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
