// C++ code for calculate body-orientational order parameter
// Written by Ruijian Zhu (ITP-CAS) in March 2023

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include<string>
#include<sstream>

const int molecule = 2496;      // number of molecules
const int frame = 5;         // number of snapshots
const int dumpHeadLine = 5;     // skip the first several lines in the original traj until box information
const int comHeadLine = 9;      // skip the first several lines in the COM traj until coordinate information
double dc = 1.5*1.5;            // cutoff distance for atom-COM distance to judge whether they are in the same image
                                // This have to be changed for polygons with many edges and it also depends on bond length 
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
    ss << "bo_" << j << ".dat";
    return ss.str();
}

struct CartesianComplex     // A struct used to save complex number
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

// This program is designed for pentagon -> cos(10*theta), sin(10*theta)
double cosTheta(double cost)    // This function and next function have to be replaced if used for other polygons
{
    double cos2 = cost * cost;
    double cos4 = cos2 * cos2;
    return ((-1+2*cos2)*(256*cos4*cos4-512*cos4*cos2+304*cos4-48*cos2+1));
}

double sinTheta(double cost,double sint)
{
    double a1 = cost * sint * (4 * sint * sint + 2 * sint - 1);
    double a2 = 4 * sint * sint - 2 * sint - 1;
    double a3 = -20 * sint * sint + 5 + 16 * sint * sint * sint * sint;
    return 2 * a1 * a2 * a3;
}

CartesianComplex bodyOrientation(coordinate com, coordinate a1, coordinate a2, coordinate a3, coordinate a4, coordinate a5)
{
    // There is no guarantee that an arbitrary atom and the COM locate in the same image
    // But at least one atom has to locate in the same image with the COM
    // This is designed for pentagon, 5 atoms, which should be modified for other cases
    CartesianComplex bodyOrient;
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
    // All the atoms have a distance larger than dc with the COM
    else    // Two possible reasons: (1) dc is too small (2) something wrong in the trajectory file
    {
        bodyOrient.costheta = 1000;
        bodyOrient.sintheta = 1000;
        alright = 0;
    }
    if(alright)
    {
        cost = deltaVector.x / modLength(deltaVector);
        sint = deltaVector.y / modLength(deltaVector);
        bodyOrient.costheta = cosTheta(cost);
        bodyOrient.sintheta = sinTheta(cost,sint);
    }
    return bodyOrient;
}

int main()
{
    double box[6];      // For triclinic box, there should be 9 elements
    coordinate atom[molecule*n];
    coordinate com[molecule];
    CartesianComplex orient[molecule];
    int id[molecule*n];
    int mol[molecule*n];
    int type[molecule*n];
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
    // Now we could read in information from traj
    // But jump the first snapshot in original traj file 
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
    // Typical dump file: id mol type x y z ix iy
    for(int j = 0; j < n*molecule; j++)
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
        for(int j = 0; j < n*molecule; j++)
        {
            infile1 >> id[j] >> mol[j] >> type[j] >> d[2*j] >> d[2*j+1] >> z >> ix >> iy; 
        }
        // Note that the line index is not the atom index, so we first read in then assign
        for(int j = 0; j < n*molecule; j++)
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
        for(int j = 0; j < comHeadLine; )
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
        string openWrite = writefilename(i);
        write.open(openWrite,ios::app);
        for(int j = 0; j < molecule; j++)
        {
            orient[j] = bodyOrientation(com[j],atom[n*j],atom[n*j+1],atom[n*j+2],atom[n*j+3],atom[n*j+4]);
            write << j+1 << ' ' << orient[j].costheta << ' ' << orient[j].sintheta << endl;
        }
        write.close();
    }
    infile1.close();
    infile2.close();
    return 0;
}
