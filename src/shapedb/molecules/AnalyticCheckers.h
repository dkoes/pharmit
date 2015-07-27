/*
 * AnalyticCheckers.h
 *
 *  Created on: Jun 2, 2014
 *      Author: dkoes
 *
 *  Various classes for supporting analytic calcultation of molecular surfaces.
 */


#ifndef ANALYTICCHECKERS_H_
#define ANALYTICCHECKERS_H_


#include <iostream>
#include <string>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/atom.h>
#include <cmath>
#include "MolSphere.h"
#include "Cube.h"
#include "MGrid.h"
#include "PMol.h"


/* given three close spheres and probe radius, does the necessary
 * calculations to determine if a point is inaccessible by the probe
 */
class TriChecker
{
	double probeSq; //probe radius squared
	double h_ijk; //height of probe
	vector3 pc1, pc2; //probe centers

	vector3 normals[6]; //plane normals
	vector3 points[6]; //points on the plane - the centers of the atoms
	static float distSq(float x1, float y1, float z1, float x2, float y2,
			float z2)
	{
		float X = (x2 - x1);
		X *= X;
		float Y = (y2 - y1);
		Y *= Y;
		float Z = (z2 - z1);
		Z *= Z;

		return X + Y + Z;
	}

public:
	double height() const
	{
		return h_ijk;
	}
	TriChecker(const MolSphere& a, const MolSphere& b, const MolSphere& c,
			float p)
	{
		probeSq = p * p;
		//calculate position of two spheres that touch a, b and c
		//take these calculations from connely
		vector3 a_i(a.x, a.y, a.z);
		vector3 a_j(b.x, b.y, b.z);
		vector3 a_k(c.x, c.y, c.z);

		//distances
		double d_ij = (a_j - a_i).length();
		double d_ik = (a_k - a_i).length();

		//torus axis
		vector3 u_ij = (a_j - a_i) / d_ij;
		vector3 u_ik = (a_k - a_i) / d_ik;

		//torus center
		vector3 t_ij = .5 * (a_i + a_j)
				+ .5 * (a_j - a_i)
						* ((a.r + p) * (a.r + p) - (b.r + p) * (b.r + p))
						/ (d_ij * d_ij);
		vector3 t_ik = .5 * (a_i + a_k)
				+ .5 * (a_k - a_i)
						* ((a.r + p) * (a.r + p) - (c.r + p) * (c.r + p))
						/ (d_ik * d_ik);

		//base triangle angle
		double w_ijk = acos(dot(u_ij, u_ik));

		//base plane normal vector
		vector3 u_ijk = cross(u_ij, u_ik) / sin(w_ijk);

		//torus base point unit vector
		vector3 u_tb = cross(u_ijk, u_ij);

		//base point
		vector3 b_ijk = t_ij + u_tb * (dot(u_ik, t_ik - t_ij)) / sin(w_ijk);

		//probe height
		h_ijk = sqrt((a.r + p) * (a.r + p) - (b_ijk - a_i).length_2());

		pc1 = b_ijk + h_ijk * u_ijk;
		pc2 = b_ijk - h_ijk * u_ijk;

		//calculate six planes that are defined by the centers of a,b,c and where
		//they intersect the two probe spheres.
		//normals all point in
		points[0] = a_i;
		points[1] = a_i;
		points[2] = a_j;
		points[3] = a_j;
		points[4] = a_k;
		points[5] = a_k;

		normals[0] = cross(a_j - a_i, pc1 - a_i);
		normals[1] = cross(a_j - a_i, pc2 - a_i);
		normals[2] = cross(a_k - a_j, pc1 - a_j);
		normals[3] = cross(a_k - a_j, pc2 - a_j);
		normals[4] = cross(a_i - a_k, pc1 - a_k);
		normals[5] = cross(a_i - a_k, pc2 - a_k);

		//make absolutely sure signs are correct
		for (unsigned i = 0; i < 6; i++)
		{
			unsigned index = (i + 4) % 6; //opposite
			vector3 opppt = points[index];
			double val = dot(opppt - points[i], normals[i]);
			if (val < 0)
			{
				normals[i] = -normals[i];
			}
		}
	}

	bool isValid() const
	{
		if (::isfinite(h_ijk))
		{
			//ignore cusps where the sasa isn't connected
			if (height() * height() < probeSq)
				return false;
			return true;
		}
		return false;
	}

	//return true if checking is necessary
	static bool isNonTrivial(const MolSphere& a, const MolSphere& b,
			const MolSphere& c, float probe)
	{
		float d = distSq(a.x, a.y, a.z, b.x, b.y, b.z);
		float r = a.r + b.r + probe * 2;
		r *= r;

		if (d >= r)
			return false;
		d = distSq(a.x, a.y, a.z, c.x, c.y, c.z);
		r = a.r + c.r + probe * 2;
		r *= r;

		if (d >= r)
			return false;

		d = distSq(b.x, b.y, b.z, c.x, c.y, c.z);
		r = b.r + c.r + probe * 2;
		r *= r;

		if (d >= r)
			return false;

		return true;
	}

	bool containsPoint(float x, float y, float z) const
	{
		//can't be within the probe
		float d1 = distSq(x, y, z, pc1.x(), pc1.y(), pc1.z());
		if (d1 < probeSq)
			return false;
		float d2 = distSq(x, y, z, pc2.x(), pc2.y(), pc2.z());
		if (d2 < probeSq)
			return false;

		//if it is within the double pyramid, then contained
		vector3 pt(x, y, z);
		for (unsigned i = 0; i < 6; i++)
		{
			double val = dot(pt - points[i], normals[i]);
			if (val < 0)
				return false;
		}

		return true;
	}

};

/* given two close spheres and a probe radius, this class does the necessary
 * calculations to determine if a point is accessible by the probe
 */
class ToroidChecker
{
	MolSphere a;
	MolSphere b;
	float probe;
	float smallSq; //shortest distance to probe
	float largeSq; //longest, from line between a and b

	float AminxSq; //shortest distance from A allowed
	float AmaxxSq; //longest from A that is acceptable
	float BminxSq; //shortest distance from B allowed
	float BmaxxSq; //longest from B that is acceptable
	float probecY; //distance to probe center from line
	float probecX; //x-component of probe center from a

	float asq;
	float bsq;

	float AcheckR;

	vector<TriChecker> tricheckers;

	static float distSq(float x1, float y1, float z1, float x2, float y2,
			float z2)
	{
		float X = (x2 - x1);
		X *= X;
		float Y = (y2 - y1);
		Y *= Y;
		float Z = (z2 - z1);
		Z *= Z;

		return X + Y + Z;
	}

	//return true of x,y,z is within the gap between a and b
	bool toroidContainsPoint(float x, float y, float z) const
	{
		//find the closest point on the line between a and b
		float Ax = a.x, Ay = a.y, Az = a.z;
		float Bx = b.x, By = b.y, Bz = b.z;
		float Cx = x, Cy = y, Cz = z;

		//must be within SA radius
		float Bsq = distSq(Bx, By, Bz, Cx, Cy, Cz);
		if (Bsq > bsq)
			return false;
		float Asq = distSq(Ax, Ay, Az, Cx, Cy, Cz);
		if (Asq > asq)
			return false;

		float t1 = Ax * Ax;
		float t3 = Bx * Ax;
		float t5 = Ay * Ay;
		float t7 = By * Ay;
		float t9 = Az * Az;
		float t11 = Bz * Az;
		float t13 = t1 - Ax * Cx - t3 + Bx * Cx + t5 - Ay * Cy - t7 + By * Cy
				+ t9 - Az * Cz - t11 + Bz * Cz;
		float t15 = Bx * Bx;
		float t17 = By * By;
		float t19 = Bz * Bz;
		float u = t13
				/ (t1 - 2.0 * t3 + t15 + t5 - 0.2e1 * t7 + t17 + t9 - 2.0 * t11
						+ t19);

		//calculate actual point on line
		float Px = Ax + u * (Bx - Ax);
		float Py = Ay + u * (By - Ay);
		float Pz = Az + u * (Bz - Az);

		//distance with A
		float Ad = distSq(Ax, Ay, Az, Px, Py, Pz);
		if (Ad < AminxSq || Ad > AmaxxSq)
			return false;
		float Bd = distSq(Bx, By, Bz, Px, Py, Pz);
		if (Bd < BminxSq || Bd > BmaxxSq)
			return false;

		//distance between C and P
		float d = distSq(Px, Py, Pz, Cx, Cy, Cz);
		if (d < smallSq) //definitely fits
		{
			return true;
		}
		if (d > largeSq)
			return false;

		//falls somewhere in the curvature region
		//compute distance from probe
		float xd = probecX - Ad;
		xd *= xd;
		float yd = probecY - d;
		yd *= yd;
		if (xd + yd < probe * probe)
			return false;
		return true;
	}
public:
	//return true if x and y are close enough to need a checker
	static bool isNonTrivial(const MolSphere& a, const MolSphere& b,
			float probe)
	{
		float d = distSq(a.x, a.y, a.z, b.x, b.y, b.z);
		float r = a.r + b.r + probe * 2;
		r *= r;

		return d < r;
	}

	friend void swap(ToroidChecker& a, ToroidChecker& b)
	{
		using std::swap;
		// bring in swap for built-in types
		swap(a.a, b.a);
		swap(a.b, b.b);
		swap(a.probe, b.probe);
		swap(a.smallSq, b.smallSq);
		swap(a.largeSq, b.largeSq);

		swap(a.AminxSq, b.AminxSq);
		swap(a.AmaxxSq, b.AmaxxSq);
		swap(a.BminxSq, b.BminxSq);
		swap(a.BmaxxSq, b.BmaxxSq);
		swap(a.probecY, b.probecY);
		swap(a.probecX, b.probecX);

		swap(a.asq, b.asq);
		swap(a.bsq, b.bsq);

		swap(a.AcheckR, b.AcheckR);

		swap(a.tricheckers, b.tricheckers);
	}

	ToroidChecker(const MolSphere& x, const MolSphere& y, float p) :
			a(x), b(y), probe(p)
	{
		//precalculate distances from the line between a and b that define the curved region
		//distances
		double ap = a.r + probe;
		double bp = b.r + probe;
		double ab = sqrt(distSq(a.x, a.y, a.z, b.x, b.y, b.z));
		//angles
		double A = acos((ap * ap + ab * ab - bp * bp) / (2 * ap * ab));
		double B = acos((bp * bp + ab * ab - ap * ap) / (2 * bp * ab));

		//bound both A and B distances since distance is unsigned
		double mina = a.r * cos(A);
		double minb = b.r * cos(B);

		AminxSq = mina * mina;
		AmaxxSq = ab - minb;
		AmaxxSq *= AmaxxSq;

		BminxSq = minb * minb;
		BmaxxSq = ab - mina;
		BmaxxSq *= BmaxxSq;

		double maxa = a.r * sin(A);
		double maxb = b.r * sin(B);
		largeSq = max(maxa, maxb);
		largeSq *= largeSq;
		//also the projection of the probe center
		probecY = (a.r + probe) * sin(A);
		probecX = (a.r + probe) * cos(A);
		smallSq = max(probecY - probe, 0.0f);
		smallSq *= smallSq;

		//calculate radius from A that includes the full toroid
		float axr = ab - minb;
		float ayr = maxb;
		asq = axr * axr + ayr * ayr;

		float bxr = ab - mina;
		float byr = maxa;
		bsq = bxr * bxr + byr * byr;
	}

	float getACheckRadius() const
	{
		return sqrt(asq);
	}

	bool containsPoint(float x, float y, float z) const
	{
		if (toroidContainsPoint(x, y, z))
			return true;

		//check tri-points
		for (unsigned i = 0, n = tricheckers.size(); i < n; i++)
		{
			if (tricheckers[i].containsPoint(x, y, z))
			{
				return true;
			}
		}
		return false;
	}

	//create tri checkers
	void addNeighbors(const vector<MolSphere>& spheres, unsigned start,
			float probe)
	{
		for (unsigned i = start, n = spheres.size(); i < n; i++)
		{
			if (TriChecker::isNonTrivial(a, b, spheres[i], probe))
			{
				TriChecker check(a, b, spheres[i], probe);
				if (check.isValid())
				{
					tricheckers.push_back(check);
				}
			}
		}

	}
};

/* This class performs intersection testing of cubes and points within the
 * context of a set of connected spheres.  It will check the sphere radius and
 * any connected toroids and triploids.
 */
class SphereChecker
{
	MolSphere atom;
	MolSphere sphere; //sphere with adjusted check radius
	vector<ToroidChecker> toroids;
public:
	SphereChecker(const MolSphere& a) :
			atom(a), sphere(a)
	{
	}

	friend void swap(SphereChecker& a, SphereChecker& b)
	{
		using std::swap;
		// bring in swap for built-in types

		swap(a.atom, b.atom);
		swap(a.sphere, b.sphere);
		swap(a.toroids, b.toroids);
	}

	//calculate toroids and triploid if necessary
	//update radius of sphere to max check distance
	void addNeighbors(const vector<MolSphere>& spheres, unsigned start,
			float probe)
	{
		for (unsigned i = start, n = spheres.size(); i < n; i++)
		{
			MolSphere b = spheres[i];
			if (probe > 0 && ToroidChecker::isNonTrivial(atom, b, probe))
			{
				toroids.push_back(ToroidChecker(atom, b, probe));
				float checkR = toroids.back().getACheckRadius();
				if (checkR > sphere.r)
					sphere.r = checkR;

				toroids.back().addNeighbors(spheres, i + 1, probe);
			}
		}
	}

	//cube possibly intersects the represented volume
	bool intersectsCube(const Cube& cube) const
	{
		return sphere.intersectsCube(cube);
	}

	//point is definitely contained in represented volume
	bool containsPoint(float x, float y, float z) const
	{
		if (atom.containsPoint(x, y, z))
			return true;

		if (sphere.containsPoint(x, y, z))
		{
			for (unsigned i = 0, n = toroids.size(); i < n; i++)
			{
				if (toroids[i].containsPoint(x, y, z))
					return true;
			}
		}
		return false;
	}

};
#endif
