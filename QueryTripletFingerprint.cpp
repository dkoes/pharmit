/*
Pharmer: Efficient and Exact 3D Pharmacophore Search
Copyright (C) 2011  David Ryan Koes and the University of Pittsburgh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/*
 * QueryTripletFingerprint.cpp
 *
 *  Created on: Nov 2, 2010
 *      Author: dkoes
 */

#include "QueryTripletFingerprint.h"
#include "CommandLine2/CommandLine.h"
#include "BoundingBox.h"
#include "Triplet.h"
#include "PharmerQuery.h"
#include "BitSetTree.h"
#include "basis.h"

extern cl::opt<bool> Quiet;
cl::opt<bool> SkipFingers("skip-fingers", cl::desc("don't check triplet fingerprints"), cl::Hidden);
cl::opt<bool> SkipRangedFingers("skip-ranged-fingers", cl::desc("don't check triplet fingerprint ranges"), cl::Hidden);

struct Plane
{
	vector3 normal;
	vector3 point;

	Plane(const vector3& n, const vector3& p): normal(n), point(p)
	{
		assert(isfinite(n.length_2()));
		assert(isfinite(p.length_2()));
	}

	double distanceToPoint(const vector3& p)
	{
		return dot(normal,p-point);
	}
};

//compute the two similarity points between a and b (inside then outside)
//return false if a and b are at the same place
//assuming points in the z-plane
//if inner point doesn't exist, set to -inf
static bool computeSimilarityPoints(vector3 a, vector3 b, double ar, double br, vector<vector3>& simAB)
{
	simAB.clear();
	simAB.resize(2);
	if(a == b)
		return false;

	//inside point - circles must not overlap
	vector3 d = a-b;
	if(d.length() < ar+br)
	{
		simAB[0] = vector3(HUGE_VAL,HUGE_VAL,HUGE_VAL);
	}
	else
	{
		//top of a with bottom of b; intersect line segment with z plane
		vector3 a2 = a; a2.z()+=ar;
		vector3 b2 = b; b2.z()-=br;
		vector3 l = a2-b2;
		l.normalize();
		vector3 p = (-a2.z()/l.z())*l+a2;
		simAB[0] = p;
	}

	//outside
	if(ar == br) //parallel to z
	{
		simAB[1] = vector3(HUGE_VAL,HUGE_VAL,HUGE_VAL);
	}
	else
	{
		vector3 a2 = a; a2.z()+=ar;
		vector3 b2 = b; b2.z()+=br;
		vector3 l = a2-b2;
		l.normalize();
		vector3 p = (-a2.z()/l.z())*l+a2;
		simAB[1] = p;
	}
	return true;
}

//compute inserection of A with sphere S and plane N (which has A as a point)
static vector3 computeSSPIntersection(vector3 A, double ar, vector3 S, double sr, vector3 N, bool aup)
{
	/* the following is derived from the maple:
with(codegen,optimize,makeproc,cost);
C1 := (x-x1)^2 + (y-y1)^2 + (z)^2 = r1^2;
C2 := (x-x2)^2 + (y-y2)^2 + (z)^2 = r2^2;
L :=  (x-x1)*nx+(y-y1)*ny = 0;
SOL := [op([allvalues(solve({C1,C2,L},{x,y,z}))][1])];

CodeGeneration[C]([optimize(SOL,tryhard)],coercetypes=false);

	 *
	 */
	assert(fabs(A.z()) < THRESHOLD  && fabs(S.z()) < THRESHOLD && fabs(N.z()) < THRESHOLD);
	double x1 = A.x();
	double y1 = A.y();
	double x2 = S.x();
	double y2 = S.y();
	double r1 = ar;
	double r2 = sr;
	double nx = N.x();
	double ny = N.y();

	double t4 = y2 * y2;
	double t47 = y1 * y2;
	double t7 = y1 * y1;
	double t45 = -2 * t47 + t7;
	double t34 = t4 + t45;
	double t40 = x2 * x1;
	double t48 = 4 * t40;
	double t50 = t48 - 2 * t34;
	double t39 = -x2 + x1;
	double t35 = (y1 - y2) * nx;
	double t46 = t39 * ny;
	double t10 = x2 * x2;
	double t17 = r1 * r1;
	double t13 = x1 * x1;
	double t33 = -t13 + t34;
	double t44 = 2 * (2 * t40 - t10 + t33) * t17;
	double t43 = t10 + t17;
	double t42 = 4 * t47;
	double t1 = 0.1e1 / (-t46 + t35);
	double t41 = t1 / 2;
	double t32 = t13 + t4 - 2 * t40 + t43;
	double t15 = r2 * r2;
	double t30 = t34 * t48 + (-t13 + t50) * t13 - t15 * t15 - t17 * t17 + (-6
			* t7 + t42 - t4) * t4 + 2 * (t32 + t45) * t15 + (t42 - t7) * t7
			+ (-6 * t13 - t10 + t50) * t10;
	double x = (2 * x1 * t35 + (-t15 + t33 + t43) * ny) * t41;
	double y = -(2 * y1 * t46 + (-t7 - t15 + t32) * nx) * t1 / 2;
	double z = sqrt((t30 + t44) * nx * nx + (-8 * t39 * t17 * t35 + (t30 - t44) * ny)
			* ny) * t41;

	if((z>= 0) != aup)
		z = -z;

	return vector3(x,y,z);
}

//compute intersection of 3 spheres, choose which point based on aup
static vector3 computeSSSIntersection(vector3 A, double r1, vector3 SAB, double r2, vector3 SCA, double r3, bool aup)
{
/* this comes from the maple code:
with(codegen);
C1 := (x-x1)^2 + (y-y1)^2 + (z)^2 = r1^2;
C2 := (x-x2)^2 + (y-y2)^2 + (z)^2 = r2^2;
C3 := (x-x3)^2 + (y-y3)^2 + (z)^2 = r3^2;

SOL := [op([allvalues(solve({C1,C2,C3},{x,y,z}))][1])]:

CodeGeneration[C]([optimize(SOL,tryhard)],coercetypes=false);

 *
 */
	assert(fabs(A.z()) < THRESHOLD && fabs(SAB.z()) < THRESHOLD && fabs(SCA.z()) < THRESHOLD);
	double x1 = A.x();
	double y1 = A.y();
	double x2 = SAB.x();
	double y2 = SAB.y();
	double x3 = SCA.x();
	double y3 = SCA.y();

	double t149 = y2 * y2;
	double t165 = r2 * r2;
	double t246 = t149 - t165;
	double t146 = y3 * y3;
	double t163 = r3 * r3;
	double t281 = t163 - t146;
	double t145 = y3 * t146;
	double t268 = 2 * t145;
	double t233 = 2 * t146;
	double t148 = y2 * t149;
	double t267 = 2 * t148;
	double t280 = 2 * t149;
	double t152 = y1 * y1;
	double t151 = y1 * t152;
	double t279 = 2 * t151;
	double t278 = 2 * t152;
	double t155 = x3 * x3;
	double t154 = x3 * t155;
	double t264 = 2 * t154;
	double t158 = x2 * x2;
	double t157 = x2 * t158;
	double t263 = 2 * t157;
	double t161 = x1 * x1;
	double t160 = x1 * t161;
	double t262 = 2 * t160;
	double t228 = 2 * t161;
	double t277 = 2 * t163;
	double t167 = r1 * r1;
	double t243 = t163 - t167;
	double t203 = 2 * t243;
	double t259 = -2 * t165;
	double t275 = t163 + t259;
	double t274 = t146 * t146 + t163 * t163 + t155 * t155;
	double t273 = t149 * t149 + t165 * t165 + t158 * t158;
	double t272 = t161 * t161 + t167 * t167 + t152 * t152;
	double t258 = -2 * t167;
	double t211 = t278 + t258 + 4 * y2 * y3;
	double t261 = -2 * t163;
	double t271 = -(8 * y3 - 4 * y2) * y1 + t261 + t233 - 4 * t149 + t211;
	double t270 = -(8 * y2 - 4 * y3) * y1 + t259 + t280 - 4 * t146 + t211;
	double t269 = -2 * y3;
	double t257 = x2 * x3;
	double t256 = y1 * y2;
	double t255 = y1 * y3;
	double t143 = 1 / ((-y2 + y1) * x3 + (-y1 + y3) * x2 + (y2 - y3) * x1);
	double t252 = t143 / 2;
	double t251 = x3 * t157;
	double t250 = t167 * t165;
	double t249 = -2 * t257;
	double t248 = 2 * t257;
	double t247 = -2 * t255;
	double t245 = t163 - t155;
	double t244 = t163 + t165;
	double t242 = -t165 + t163;
	double t241 = t165 - t158;
	double t240 = -t167 + t165;
	double t239 = t167 - t152;
	double t238 = t167 - t158;
	double t237 = 2 * x2;
	double t236 = 2 * x3;
	double t235 = 4 * y1;
	double t234 = 2 * y3;
	double t230 = 2 * t155;
	double t229 = 2 * t158;
	double t226 = 2 * t167;
	double t225 = t158 + t249;
	double t223 = 2 * t274;
	double t221 = 2 * t273;
	double t219 = t155 - t243;
	double t218 = 2 * t272;
	double t217 = t258 + t244;
	double t216 = -t161 + t241;
	double t215 = -t163 + t241;
	double t214 = t163 + t241;
	double t213 = -t146 + t239;
	double t212 = -t155 + t238;
	double t209 = t155 - t214;
	double t208 = -t152 + t219;
	double t207 = -t146 - t215;
	double t206 = t239 - t241;
	double t205 = -t167 - t216;
	double t202 = 2 * t240;
	double t200 = 4 * t239;
	double t199 = y1 * t267 - 2 * t160 * x2;
	double t198 = t148 * t269 + 2 * t251;
	double t197 = -t225 + t245;
	double t196 = t255 + t213;
	double t195 = -2 * t158 - t167 + t248 - t245;
	double t194 = -4 * t257 + t209;
	double t193 = t257 - t155 + t216;
	double t192 = -t161 + t146 + t208;
	double t191 = t161 - t149 - t206;
	double t190 = t149 - t155 + t207;
	double t189 = t248 - 2 * t155 + t205;
	double t186 = t262 + t263 - t240 * t237;
	double t185 = t165 * t277 - t273 - t274;
	double t184 = t167 * t277 - t272 - t274;
	double t183 = 2 * t250 - t272 - t273;
	double x = -(t191 * y3 + t192 * y2 + t190 * y1) * t143 / 2;
	double y = (t191 * x3 + t192 * x2 + t190 * x1) * t252;
	double z = sqrt((t263 + t242 * t237) * t160 + (t279 + y1 * t203) * t148 + t186 * t154 + (t267 + t279 + y1 * t202) * t145 + (y1 * t268 + t184) * t158 + (-t157 * t203 + t215 * t262 + (-2 * t250 + t163 * t202 + t218) * x2) * x3 + (t217 * t229 + t183 + t199) * t155 + ((-t155 + t244) * t229 + t186 * x3 + t183 - t199) * t146 + ((t165 - t197) * t279 + t219 * t267 + (-4 * t251 - t245 * t229 + t197 * t226 + t221 + 2 * t195 * t165) * y1) * y3 + ((t165 + t167) * t229 + (t163 + t238) * t230 + (t264 + (t258 - t244) * t236) * x2 + (t257 + t212 + t275) * t233 + t185 + t198) * t152 + ((t163 + t239) * t233 + (t267 + t268 - t242 * t234) * y1 + (t261 + t165 + t196) * t229 + (-2 * t154 + (-t146 + t278 + t247 + t217) * t236) * x2 + (3 * t158 + t196 + t275) * t230 + t185 - t198) * t161 + (x3 * t262 + t151 * t269 + t214 * t230 + (t264 + t262 - x3 * t203) * x2 + (-t257 + t165 + t212) * t228 + (t261 + t167 + t193) * t278 + (-2 * t145 + (-t161 + t259 - t195) * t234) * y1 + (3 * t152 + t258 + t163 + t193) * t233 + t184) * t149 + ((t161 + t206) * t268 + (t155 + t249 + t207) * t279 + ((t258 + 4 * t257 + t209) * t228 + 4 * (-x2 - x3) * t160 + t218 + t194 * t226 + (t228 + t258 - t194) * t278 + (t277 - t230) * t241) * y3 + (-4 * t154 * x2 + (t261 - t189) * t233 + t223 + t205 * t230 + t189 * t277 + (t226 - t228) * (t165 - t225)) * y1) * y2 + (((4 * t256 - t239 + t246) * t233 + t270 * t163 + 4 * (-y1 - y2) * t145 + (y2 * t200 - t246 * t235) * y3 + (-4 * t163 - t270) * t155 + t223) * x2 + (-2 * t256 + t149 - t206) * t264 + (t247 + t146 - t208) * t263 + (t271 * t165 + 4 * (-y3 - y1) * t148 + t221 + (y3 * t200 + t281 * t235) * y2 + (-4 * t165 - t271) * t158 + (4 * t255 - t163 - t213) * t280) * x3 + 2 * (-t246 * x2 + t281 * x3) * t239) * x1) * t252;

	if((z>= 0) != aup)
		z = -z;

	return vector3(x,y,z);

}


//compute the tangent point by finding the interection of three spheres centered
//at A, SAB, and SCA; A has radius ar, SAB and SCA have a radius equal to the distance to A
//choose the point above or below the z-axis using aup
//special case: an S point at infinity - match to plane
static vector3 computeTangentPoint(vector3 A, double ar, vector3 B, vector3 C, vector3 SAB, vector3 SCA, bool aup)
{
	if(isfinite(SAB.x()) && isfinite(SCA.x()))
	{
		double ab = (SAB-A).length();
		double ac = (SCA-A).length();
		return computeSSSIntersection(A, ar, SAB, ab, SCA, ac, aup);
	}
	else if(isfinite(SAB.x())) //sphere plane
	{
		double ab = (SAB-A).length();
		return computeSSPIntersection(A, ar, SAB, ab, (C-A).normalize(), aup);
	}
	else if(isfinite(SCA.x())) //sphere plane
	{
		double ac = (SCA-A).length();
		return computeSSPIntersection(A, ar, SCA, ac, (B-A).normalize(), aup);
	}
	else //plane-plane - just top of a
	{
		if(aup)
			return vector3(A.x(), A.y(), A.z()+ar);
		else
			return vector3(A.x(), A.y(), A.z()-ar);
	}
}

//compute a tangent plane for A,B,C spheres of radius ar,br,cr
//the location of the plane relative to the z-plane at each sphere
//is given by aup,bup,cup
//simularity points (inside-outside) are given in sim*
static void computeTangentPlane(vector3 A, vector3 B, vector3 C, double ar, double br,
		double cr, bool aup, bool bup, bool cup, const vector<vector3>& simAB,
		const vector<vector3>& simBC, const vector<vector3>& simCA,
		vector<Plane>& planes)
{
	vector3 SAB, SBC,SCA;
	//compute point tangent to A
	if(aup == bup) //simab outside
		SAB = simAB[1];
	else
		SAB = simAB[0];

	if(aup == cup)
		SCA = simCA[1];
	else
		SCA = simCA[0];

	if(bup == cup)
		SBC = simBC[1];
	else
		SBC = simBC[0];

	//check for non-existent similarity points
	if(!isfinite(SCA.x()) && SCA.x() < 0)
		return;
	if(!isfinite(SAB.x()) && SAB.x() < 0)
		return;
	if(!isfinite(SBC.x()) && SBC.x() < 0)
		return;

	vector3 tanA = computeTangentPoint(A, ar, B, C, SAB, SCA, aup);
	if(!isfinite(tanA.length_2()))
		return; //not all up/down combinations may be possible
	vector3 tanB = computeTangentPoint(B, br, A, C, SAB, SBC, bup);
	if(!isfinite(tanB.length_2()))
		return;
	vector3 tanC = computeTangentPoint(C, cr, B, A, SBC, SCA, cup);
	if(!isfinite(tanC.length_2()))
		return;

	planes.push_back(Plane(cross(tanB-tanA,tanC-tanA).normalize(), tanA));
}

//compute tangent planes to three spheres
//planes are in the passed coordinate basis which is assumed to have a at the orgin and b on x-axis
//if no plane exists, the normal is returned as zero and false is returned
static void computeTangentPlanes(const CoordinateBasis& basis, const PharmaPoint& a, const PharmaPoint& b, const PharmaPoint& c,
		vector< Plane >& planes)
{
	planes.clear();

	vector3 A, B, C;
	basis.replot(a.x,a.y,a.z, A.x(),A.y(),A.z());
	basis.replot(b.x,b.y,b.z, B.x(),B.y(),B.z());
	basis.replot(c.x,c.y,c.z, C.x(),C.y(),C.z());


	//calculate the similarity points
	vector<vector3> simAB;
	vector<vector3> simBC;
	vector<vector3> simCA;

	if(!computeSimilarityPoints(A, B, a.radius, b.radius, simAB))
		return;
	if(!computeSimilarityPoints(B, C, b.radius, c.radius, simBC))
		return;
	if(!computeSimilarityPoints(C, A, c.radius, a.radius, simCA))
		return;

	//this is a little redundant.. really just need to compute 4 planes than flip them over the z-axis
	for(unsigned i = 0; i < 2; i++)
	{
		for(unsigned j = 0; j < 2; j++)
		{
			for(unsigned k = 0; k < 2; k++)
			{
				computeTangentPlane(A,B,C,a.radius, b.radius, c.radius,i,j,k,simAB,simBC,simCA, planes);
			}
		}
	}
}

void QueryTripletFingerprint::set(unsigned i, unsigned j, unsigned k,
		const PharmerQuery& query)
{
	const vector<PharmaPoint>& points = query.getPoints();
	//generate ALL possible points
	bigqfingers.clear();
	smallqfingers.clear();
	bigqfingers.resize(points.size());
	smallqfingers.resize(points.size());

	//figure out tangent planes as bounds on the chirality
	CoordinateBasis basis(points[i],points[j],points[k],true);
	vector<Plane> tangentPlanes;
	computeTangentPlanes(basis, points[i],points[j],points[k], tangentPlanes);

	unsigned cnt = 0;
	for (unsigned p = 0, n = points.size(); p < n; p++)
	{
		if (p == i)
			continue;
		if (p == j)
			continue;
		if (p == k)
			continue;
		if(points[p].requirements != PharmaPoint::Required)
			continue;

		//precompute possibiliy of multiple chirality
		const int UP = 1;
		const int DOWN = 2;
		unsigned chiral = 0;

		if(tangentPlanes.size() == 0) //co-linear
		{
			chiral = UP|DOWN;
		}
		else
		{
			for (unsigned t = 0, nt = tangentPlanes.size(); t < nt; t++)
			{
				vector3 P;
				basis.replot(points[p].x, points[p].y, points[p].z, P.x(),
						P.y(), P.z());
				double d = tangentPlanes[t].distanceToPoint(P);
				if (fabs(d) <= points[p].radius)
					chiral |= (DOWN | UP);
				else if (d < 0)
					chiral |= DOWN;
				else
					chiral |= UP;
			}
		}

		vector<TripletFingerprint> smallvals[2]; smallvals[0].reserve(4*4*4); smallvals[1].reserve(4*4*4);
		//compute a val that is a combination of distances, type
		unsigned l1min = trunc((PharmaPoint::pharmaDist(points[i], points[p])
				- points[i].radius - points[p].radius) / TripletFingerprint::DISTSPACES[0]);
		unsigned l1max = trunc((PharmaPoint::pharmaDist(points[i], points[p])
				+ points[i].radius + points[p].radius) / TripletFingerprint::DISTSPACES[0]);

		for (unsigned l1 = l1min; l1 <= l1max; l1++)
		{
			unsigned l2min = trunc((PharmaPoint::pharmaDist(points[j],
					points[p]) - points[j].radius - points[p].radius)
					/ TripletFingerprint::DISTSPACES[0]);
			unsigned l2max = trunc((PharmaPoint::pharmaDist(points[j],
					points[p]) + points[j].radius + points[p].radius)
					/ TripletFingerprint::DISTSPACES[0]);
			double minip = l1*TripletFingerprint::DISTSPACES[0];
			double maxip = (l1+1)*TripletFingerprint::DISTSPACES[0];

			for (unsigned l2 = l2min; l2 <= l2max; l2++)
			{
				double minjp = l2*TripletFingerprint::DISTSPACES[0];
				double maxjp = (l2+1)*TripletFingerprint::DISTSPACES[0];

				if(!SkipRangedFingers && !query.inRange(i, j, p, minip, maxip, minjp, maxjp))
					continue;

				unsigned l3min = trunc((PharmaPoint::pharmaDist(points[k],
						points[p]) - points[k].radius - points[p].radius)
						/ TripletFingerprint::DISTSPACES[0]);
				unsigned l3max = trunc((PharmaPoint::pharmaDist(points[k],
						points[p]) + points[k].radius + points[p].radius)
						/ TripletFingerprint::DISTSPACES[0]);

				for (unsigned l3 = l3min; l3 <= l3max; l3++)
				{
					double minkp = l3*TripletFingerprint::DISTSPACES[0];
					double maxkp = (l3+1)*TripletFingerprint::DISTSPACES[0];

					if(!SkipRangedFingers && !query.inRange(j, k, p, minjp, maxjp, minkp, maxkp))
						continue;
					if(!SkipRangedFingers && !query.inRange(k, i, p, minkp, maxkp, minip, maxip))
						continue;

					unsigned long val = 0;
					val += l1;
					val *= (TripletFingerprint::MAXDISTS[0] + 1);

					val += l2;
					val *= (TripletFingerprint::MAXDISTS[0] + 1);

					val += l3;
					val *= TripletFingerprint::MAXPHARMA;
					val += points[p].pharma->index;

					val *= 2;

					smallvals[0].clear(); smallvals[1].clear();
					//now all the same but with higher resolution! this should really really be refactored..
					unsigned sml1min = minip/ TripletFingerprint::DISTSPACES[1];
					unsigned sml1max = maxip / TripletFingerprint::DISTSPACES[1];

					for (unsigned sml1 = sml1min; sml1 < sml1max; sml1++)
					{
						unsigned sml2min = minjp/ TripletFingerprint::DISTSPACES[1];
						unsigned sml2max = maxjp / TripletFingerprint::DISTSPACES[1];

						double sminip = sml1*TripletFingerprint::DISTSPACES[1];
						double smaxip = (sml1+1)*TripletFingerprint::DISTSPACES[1];

						for (unsigned sml2 = sml2min; sml2 < sml2max; sml2++)
						{
							double sminjp = sml2*TripletFingerprint::DISTSPACES[1];
							double smaxjp = (sml2+1)*TripletFingerprint::DISTSPACES[1];

							if(!SkipRangedFingers && !query.inRange(i, j, p, sminip, smaxip, sminjp, smaxjp))
								continue;

							unsigned sml3min = minkp/ TripletFingerprint::DISTSPACES[1];
							unsigned sml3max = maxkp / TripletFingerprint::DISTSPACES[1];


							for (unsigned sml3 = sml3min; sml3 < sml3max; sml3++)
							{
								double sminkp = sml3*TripletFingerprint::DISTSPACES[1];
								double smaxkp = (sml3+1)*TripletFingerprint::DISTSPACES[1];

								if(!SkipRangedFingers && !query.inRange(j, k, p, sminjp, smaxjp, sminkp, smaxkp))
									continue;
								if(!SkipRangedFingers && !query.inRange(k, i, p, sminkp, smaxkp, sminip, smaxip))
									continue;

								unsigned long val = 0;
								val += sml1;
								val *= (TripletFingerprint::MAXDISTS[1] + 1);

								val += sml2;
								val *= (TripletFingerprint::MAXDISTS[1] + 1);

								val += sml3;
								val *= TripletFingerprint::MAXPHARMA;
								val += points[p].pharma->index;

								val *= 2;

								//seperate smallvals for each chirality
								if (chiral & UP)
								{
									TripletFingerprint finger;
									TripletFingerprint::bloom(val, finger.f1, finger.f2, BloomBitsSmall);
									smallvals[0].push_back(finger);
								}
								if (chiral & DOWN)
								{
									val++;
									TripletFingerprint finger;
									TripletFingerprint::bloom(val, finger.f1, finger.f2, BloomBitsSmall);
									smallvals[1].push_back(finger);
								}
							}
						}
					}

					if (chiral & UP)
					{
						TripletFingerprint finger;
						TripletFingerprint::bloom(val, finger.f1, finger.f2, BloomBitsLarge);
						bigqfingers[p].push_back(finger);
						smallqfingers[p].push_back(smallvals[0]);
					}
					if (chiral & DOWN)
					{
						val++;
						TripletFingerprint finger;
						TripletFingerprint::bloom(val, finger.f1, finger.f2,BloomBitsLarge);
						bigqfingers[p].push_back(finger);

						smallqfingers[p].push_back(smallvals[1]);
					}
				}
			}
		}

		cnt += bigqfingers[p].size();
	}

	if(!Quiet && false)
	{
		cout << "QFinger " << cnt;
		for (unsigned p = 0, n = bigqfingers.size(); p < n; p++)
		{
			cout << " " << bigqfingers[p].size() << "(";
			unsigned smcnt = 0;
			BOOST_FOREACH(vector<TripletFingerprint>& v, smallqfingers[p])
			{
				smcnt += v.size();
			}
			cout << smcnt << ") ";
		}
		cout << "\n";
	}
}

bool QueryTripletFingerprint::isValid(const TripletFingerprint& f) const
{
	if(SkipFingers)
		return true;

	//there must be a hit in every one of the points
	//unfortunately, more clever methods (bitsettree) aren't actually faster
	for(unsigned p = 0, np = bigqfingers.size(); p < np; p++)
	{
		if(bigqfingers[p].size() > 0)
		{
			unsigned i;
			unsigned n = bigqfingers[p].size();
			for(i = 0; i < n; i++)
			{
				if(f.contains(bigqfingers[p][i]))
				{
					//check small
					unsigned j;
					unsigned m = smallqfingers[p][i].size();
					for(j = 0; j < m; j++)
					{
						if(f.contains(smallqfingers[p][i][j]))
							break;
					}
					if(j != m)
						break; //small fingers were also a match, otherwise keep going
				}
			}
			if(i == n)
			{
				return false;
			}
		}
	}

	return true;
}

