// C++ code for the calculation of translational function
// Note that the crystal axis orientation has to be determined at first
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
const double scale = 80;    // The correlation function would be calculated from 0 to 'scale' 
const int interval = 320;   // Number of bins
const double theta = 0.003; // Crystal axis
const double error = 0.2;   // Tolerance on r, note that this should be independent from r
std::string fileHead ("trans_1_all.txt"); 
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
    return sqrt(d);
}

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

bool judge(double yhi, double ylo, double xhi, double xlo, coordinate atom1, coordinate atom2)
{
    coordinate deltar;
    deltar = triclincPBC(yhi,ylo,xhi,xlo,atom1,atom2);
    double r = deltar.x;    // 'Automatically' $\delta(|x_{i}-x_{j}|-r)$
    double deltay = r * (100*theta);    // 100 is used to amplify the error, avoiding truncated error in computer
    double real = deltar.y*100;         // deltay: expected, real: real value
    if(fabs(real-deltay) < error*100)   // $\delta(|y_{i}-y_{j}|-x tan\theta)$
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
    double box[6];
    coordinate atom[unitCell];
    double id;
    double type;
    double mol;
    double z;
    double count[interval] = {};
    double r;
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
        for(int i = 0; i < unitCell ; i++)
        {
            for(int j = i; j < unitCell; j++)   // Here we count in i itself, but it would better be deleted when plotting
            {
                bool countin = judge(box[3],box[2],box[1],box[0],atom[i],atom[j]);
                if(countin)
                {
                    coordinate dr = triclincPBC(box[3],box[2],box[1],box[0],atom[i],atom[j]);
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
