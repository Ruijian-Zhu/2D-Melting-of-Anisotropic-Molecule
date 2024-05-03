// C++ code for calculate body-orientational (theta) of each monomer
// Written by Ruijian Zhu (ITP-CAS) in March 2023

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

const int molecule = 2496;          // number of molecules
const int frame = 5;             // number of snapshots
const int dumpHeadLine = 5;         // skip the first several lines in the original traj until box information
const int comHeadLine = 9;          // skip the first several lines in the COM traj until coordinate information
double dc = 1.5*1.5;
const int n = 5;                // polygon edge number
std::string comFile ("dump_1_reorder.lammpstrj");    // original traj file
std::string dumpFile ("dump_npt_o.1.lammpstrj");   // COM traj file

using namespace std;
ifstream infile1;
ifstream infile2;
ofstream write;

string writefilename(int j)
{
    stringstream ss;
    ss << "theta_" << j << ".dat";
    return ss.str();
}

struct CartesianComplex
{
    double costheta;
    double sintheta;
};

struct coordinate
{
    double x;
    double y;
};

double modLength(coordinate a)
{
    double dist = a.x * a.x + a.y * a.y;
    return sqrtf(dist);
}

double distance(coordinate a1, coordinate a2)
{
    coordinate delta;
    delta.x = a1.x - a2.x;
    delta.y = a1.y - a2.y;
    return modLength(delta);
}

double bodyOrientation(coordinate com, coordinate a1, coordinate a2, coordinate a3, coordinate a4, coordinate a5)
{
    double theta;
    coordinate deltaVector;
    bool alright = 1;
    double cost , sint;
    if(distance(com,a1)<dc)
    {
        deltaVector.x = a1.x - com.x;
        deltaVector.y = a1.y - com.y;
    }
    else if(distance(com,a2)<dc)
    {
        deltaVector.x = a2.x - com.x;
        deltaVector.y = a2.y - com.y;
    }
    else if(distance(com,a3)<dc)
    {
        deltaVector.x = a3.x - com.x;
        deltaVector.y = a3.y - com.y;
    }
    else if(distance(com,a4)<dc)
    {
        deltaVector.x = a4.x - com.x;
        deltaVector.y = a4.y - com.y;
    }
    else if(distance(com,a5)<dc)
    {
        deltaVector.x = a5.x - com.x;
        deltaVector.y = a5.y - com.y;
    }
    else
    {
        theta = 1000;
        alright = 0;
    }
    if(alright)
    {
        theta = atan2(deltaVector.y,deltaVector.x);
        while(theta>2*M_PI/n)       // Taking n-fold symmetry into account
        {
            theta -= 2*M_PI/n;
        }
        while(theta<0)
        {
            theta += 2*M_PI/n;
        }
    }
    return theta;
}

int main()
{
    double box[6];
    coordinate atom[molecule*5];
    coordinate com[molecule];
    double orient[molecule];
    int id[molecule*5];
    int mol[molecule*5];
    int type[molecule*5];
    double z;
    double ix;
    double iy;
    int i;
    double d[10*molecule];
    infile1.open(dumpFile,ios::in);
    while (infile1.bad())
	{
		cout << "cannot open the file" << endl;
		return 0;
	}
    infile2.open(comFile,ios::in);
    while (infile2.bad())
	{
		cout << "cannot open the file" << endl;
		return 0;
	}
    // Open two files successfully, below start to read
    for(int j = 0; j < dumpHeadLine; )
    {
        if(infile1.get()=='\n') j++; 
    }
    for(int j = 0; j < 6; j++)
    {
        infile1 >> box[j];
    }
    for(int j = 0; j < 2; )
    {
        if(infile1.get()=='\n') j++; 
    }
    for(int j = 0; j < 5*molecule; j++)
    {
        infile1 >> id[j] >> mol[j] >> type[j] >> d[2*j] >> d[2*j+1] >> z >> ix >> iy;  
    }
    for(int j = 0; j < 1; )
    {
        if(infile1.get()=='\n') j++; 
    }
    //Start to read useful START START START START START
    for(i = 0; i < frame; i++)
    {     
        for(int j = 0; j < dumpHeadLine; )
        {
            if(infile1.get()=='\n') j++; 
        }
        for(int j = 0; j < 6; j++)
        {
            infile1 >> box[j];
        }
        for(int j = 0; j < 2; )
        {
            if(infile1.get()=='\n') j++; 
        }
        for(int j = 0; j < 5*molecule; j++)
        {
            infile1 >> id[j] >> mol[j] >> type[j] >> d[2*j] >> d[2*j+1] >> z >> ix >> iy; 
        }
        for(int j = 0; j < 5*molecule; j++)
        {
            int index = id[j];
            atom[index-1].x = d[2*j];
            atom[index-1].y = d[2*j+1];
        }
        for(int j = 0; j < 1; )
        {
            if(infile1.get()=='\n') j++; 
        }
        //read the second file
        for(int j = 0; j < 9; )
        {
            if(infile2.get()=='\n') j++;
        }
        for(int j = 0; j < molecule; j++)
        {
            infile2 >> id[j] >> type[j] >> type[j] >> com[j].x >> com[j].y >> z; 
        }
        for(int j = 0; j < 1; )
        {
            if(infile2.get()=='\n') j++; 
        }
        //Calculate and write
        string openWrite = writefilename(i+1);
        write.open(openWrite,ios::app);
        for(int j = 0; j < molecule; j++)
        {
            orient[j] = bodyOrientation(com[j],atom[5*j],atom[5*j+1],atom[5*j+2],atom[5*j+3],atom[5*j+4]);
            write << j+1 << ' ' << orient[j] << endl;
        }
        write.close();
    }
    infile1.close();
    infile2.close();
    return 0;
}
