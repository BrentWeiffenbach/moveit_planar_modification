/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2013, Ioan A. Sucan
 *  Copyright (c) 2008-2013, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Willow Garage nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

/* Author: Ioan Sucan */

#include <moveit/robot_model/planar_joint_model.h>
#include <boost/math/constants/constants.hpp>
#include <limits>
#include <cmath>

namespace
 {
     // The comments, variable names, etc. use the nomenclature from the Reeds & Shepp paper.
     const double pi = boost::math::constants::pi<double>();
     const double twopi = 2. * pi;
 #ifndef NDEBUG
     const double RS_EPS = 1e-6;
 #endif
     const double ZERO = 10 * std::numeric_limits<double>::epsilon();

     inline double mod2pi(double x)
     {
         double v = fmod(x, twopi);
         if (v < -pi)
             v += twopi;
         else if (v > pi)
             v -= twopi;
         return v;
     }
     inline void polar(double x, double y, double &r, double &theta)
     {
         r = sqrt(x * x + y * y);
         theta = atan2(y, x);
     }
     inline void tauOmega(double u, double v, double xi, double eta, double phi, double &tau, double &omega)
     {
         double delta = mod2pi(u - v), A = sin(u) - sin(delta), B = cos(u) - cos(delta) - 1.;
         double t1 = atan2(eta * A - xi * B, xi * A + eta * B), t2 = 2. * (cos(delta) - cos(v) - cos(u)) + 3;
         tau = (t2 < 0) ? mod2pi(t1 + pi) : mod2pi(t1);
         omega = mod2pi(tau - u + v - phi);
     }

     // formula 8.1 in Reeds-Shepp paper
     inline bool LpSpLp(double x, double y, double phi, double &t, double &u, double &v)
     {
         polar(x - sin(phi), y - 1. + cos(phi), u, t);
         if (t >= -ZERO)
         {
             v = mod2pi(phi - t);
             if (v >= -ZERO)
             {
                 assert(fabs(u * cos(t) + sin(phi) - x) < RS_EPS);
                 assert(fabs(u * sin(t) - cos(phi) + 1 - y) < RS_EPS);
                 assert(fabs(mod2pi(t + v - phi)) < RS_EPS);
                 return true;
             }
         }
         return false;
     }
     // formula 8.2
     inline bool LpSpRp(double x, double y, double phi, double &t, double &u, double &v)
     {
         double t1, u1;
         polar(x + sin(phi), y - 1. - cos(phi), u1, t1);
         u1 = u1 * u1;
         if (u1 >= 4.)
         {
             double theta;
             u = sqrt(u1 - 4.);
             theta = atan2(2., u);
             t = mod2pi(t1 + theta);
             v = mod2pi(t - phi);
             assert(fabs(2 * sin(t) + u * cos(t) - sin(phi) - x) < RS_EPS);
             assert(fabs(-2 * cos(t) + u * sin(t) + cos(phi) + 1 - y) < RS_EPS);
             assert(fabs(mod2pi(t - v - phi)) < RS_EPS);
             return t >= -ZERO && v >= -ZERO;
         }
         return false;
     }
     void CSC(double x, double y, double phi, moveit::core::PlanarJointModel::ReedsSheppPath &path)
     {
         double t, u, v, Lmin = path.length(), L;
         if (LpSpLp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[14], t, u, v);
             Lmin = L;
         }
         if (LpSpLp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[14], -t, -u, -v);
             Lmin = L;
         }
         if (LpSpLp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[15], t, u, v);
             Lmin = L;
         }
         if (LpSpLp(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[15], -t, -u, -v);
             Lmin = L;
         }
         if (LpSpRp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[12], t, u, v);
             Lmin = L;
         }
         if (LpSpRp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[12], -t, -u, -v);
             Lmin = L;
         }
         if (LpSpRp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[13], t, u, v);
             Lmin = L;
         }
         if (LpSpRp(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[13], -t, -u, -v);
     }
     // formula 8.3 / 8.4  *** TYPO IN PAPER ***
     inline bool LpRmL(double x, double y, double phi, double &t, double &u, double &v)
     {
         double xi = x - sin(phi), eta = y - 1. + cos(phi), u1, theta;
         polar(xi, eta, u1, theta);
         if (u1 <= 4.)
         {
             u = -2. * asin(.25 * u1);
             t = mod2pi(theta + .5 * u + pi);
             v = mod2pi(phi - t + u);
             assert(fabs(2 * (sin(t) - sin(t - u)) + sin(phi) - x) < RS_EPS);
             assert(fabs(2 * (-cos(t) + cos(t - u)) - cos(phi) + 1 - y) < RS_EPS);
             assert(fabs(mod2pi(t - u + v - phi)) < RS_EPS);
             return t >= -ZERO && u <= ZERO;
         }
         return false;
     }
     void CCC(double x, double y, double phi, moveit::core::PlanarJointModel::ReedsSheppPath &path)
     {
         double t, u, v, Lmin = path.length(), L;
         if (LpRmL(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[0], t, u, v);
             Lmin = L;
         }
         if (LpRmL(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[0], -t, -u, -v);
             Lmin = L;
         }
         if (LpRmL(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[1], t, u, v);
             Lmin = L;
         }
         if (LpRmL(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[1], -t, -u, -v);
             Lmin = L;
         }

         // backwards
         double xb = x * cos(phi) + y * sin(phi), yb = x * sin(phi) - y * cos(phi);
         if (LpRmL(xb, yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[0], v, u, t);
             Lmin = L;
         }
         if (LpRmL(-xb, yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[0], -v, -u, -t);
             Lmin = L;
         }
         if (LpRmL(xb, -yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[1], v, u, t);
             Lmin = L;
         }
         if (LpRmL(-xb, -yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[1], -v, -u, -t);
     }
     // formula 8.7
     inline bool LpRupLumRm(double x, double y, double phi, double &t, double &u, double &v)
     {
         double xi = x + sin(phi), eta = y - 1. - cos(phi), rho = .25 * (2. + sqrt(xi * xi + eta * eta));
         if (rho <= 1.)
         {
             u = acos(rho);
             tauOmega(u, -u, xi, eta, phi, t, v);
             assert(fabs(2 * (sin(t) - sin(t - u) + sin(t - 2 * u)) - sin(phi) - x) < RS_EPS);
             assert(fabs(2 * (-cos(t) + cos(t - u) - cos(t - 2 * u)) + cos(phi) + 1 - y) < RS_EPS);
             assert(fabs(mod2pi(t - 2 * u - v - phi)) < RS_EPS);
             return t >= -ZERO && v <= ZERO;
         }
         return false;
     }
     // formula 8.8
     inline bool LpRumLumRp(double x, double y, double phi, double &t, double &u, double &v)
     {
         double xi = x + sin(phi), eta = y - 1. - cos(phi), rho = (20. - xi * xi - eta * eta) / 16.;
         if (rho >= 0 && rho <= 1)
         {
             u = -acos(rho);
             if (u >= -.5 * pi)
             {
                 tauOmega(u, u, xi, eta, phi, t, v);
                 assert(fabs(4 * sin(t) - 2 * sin(t - u) - sin(phi) - x) < RS_EPS);
                 assert(fabs(-4 * cos(t) + 2 * cos(t - u) + cos(phi) + 1 - y) < RS_EPS);
                 assert(fabs(mod2pi(t - v - phi)) < RS_EPS);
                 return t >= -ZERO && v >= -ZERO;
             }
         }
         return false;
     }
     void CCCC(double x, double y, double phi, moveit::core::PlanarJointModel::ReedsSheppPath &path)
     {
         double t, u, v, Lmin = path.length(), L;
         if (LpRupLumRm(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[2], t, u, -u, v);
             Lmin = L;
         }
         if (LpRupLumRm(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // timeflip
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[2], -t, -u, u, -v);
             Lmin = L;
         }
         if (LpRupLumRm(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[3], t, u, -u, v);
             Lmin = L;
         }
         if (LpRupLumRm(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // timeflip + reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[3], -t, -u, u, -v);
             Lmin = L;
         }

         if (LpRumLumRp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[2], t, u, u, v);
             Lmin = L;
         }
         if (LpRumLumRp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // timeflip
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[2], -t, -u, -u, -v);
             Lmin = L;
         }
         if (LpRumLumRp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[3], t, u, u, v);
             Lmin = L;
         }
         if (LpRumLumRp(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + 2. * fabs(u) + fabs(v)))  // timeflip + reflect
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[3], -t, -u, -u, -v);
     }
     // formula 8.9
     inline bool LpRmSmLm(double x, double y, double phi, double &t, double &u, double &v)
     {
         double xi = x - sin(phi), eta = y - 1. + cos(phi), rho, theta;
         polar(xi, eta, rho, theta);
         if (rho >= 2.)
         {
             double r = sqrt(rho * rho - 4.);
             u = 2. - r;
             t = mod2pi(theta + atan2(r, -2.));
             v = mod2pi(phi - .5 * pi - t);
             assert(fabs(2 * (sin(t) - cos(t)) - u * sin(t) + sin(phi) - x) < RS_EPS);
             assert(fabs(-2 * (sin(t) + cos(t)) + u * cos(t) - cos(phi) + 1 - y) < RS_EPS);
             assert(fabs(mod2pi(t + pi / 2 + v - phi)) < RS_EPS);
             return t >= -ZERO && u <= ZERO && v <= ZERO;
         }
         return false;
     }
     // formula 8.10
     inline bool LpRmSmRm(double x, double y, double phi, double &t, double &u, double &v)
     {
         double xi = x + sin(phi), eta = y - 1. - cos(phi), rho, theta;
         polar(-eta, xi, rho, theta);
         if (rho >= 2.)
         {
             t = theta;
             u = 2. - rho;
             v = mod2pi(t + .5 * pi - phi);
             assert(fabs(2 * sin(t) - cos(t - v) - u * sin(t) - x) < RS_EPS);
             assert(fabs(-2 * cos(t) - sin(t - v) + u * cos(t) + 1 - y) < RS_EPS);
             assert(fabs(mod2pi(t + pi / 2 - v - phi)) < RS_EPS);
             return t >= -ZERO && u <= ZERO && v <= ZERO;
         }
         return false;
     }
     void CCSC(double x, double y, double phi, moveit::core::PlanarJointModel::ReedsSheppPath &path)
     {
         double t, u, v, Lmin = path.length() - .5 * pi, L;
         if (LpRmSmLm(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[4], t, -.5 * pi, u, v);
             Lmin = L;
         }
         if (LpRmSmLm(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
         {
             path =
                 moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[4], -t, .5 * pi, -u, -v);
             Lmin = L;
         }
         if (LpRmSmLm(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[5], t, -.5 * pi, u, v);
             Lmin = L;
         }
         if (LpRmSmLm(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
         {
             path =
                 moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[5], -t, .5 * pi, -u, -v);
             Lmin = L;
         }

         if (LpRmSmRm(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[8], t, -.5 * pi, u, v);
             Lmin = L;
         }
         if (LpRmSmRm(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
         {
             path =
                 moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[8], -t, .5 * pi, -u, -v);
             Lmin = L;
         }
         if (LpRmSmRm(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[9], t, -.5 * pi, u, v);
             Lmin = L;
         }
         if (LpRmSmRm(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
         {
             path =
                 moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[9], -t, .5 * pi, -u, -v);
             Lmin = L;
         }

         // backwards
         double xb = x * cos(phi) + y * sin(phi), yb = x * sin(phi) - y * cos(phi);
         if (LpRmSmLm(xb, yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[6], v, u, -.5 * pi, t);
             Lmin = L;
         }
         if (LpRmSmLm(-xb, yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
         {
             path =
                 moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[6], -v, -u, .5 * pi, -t);
             Lmin = L;
         }
         if (LpRmSmLm(xb, -yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[7], v, u, -.5 * pi, t);
             Lmin = L;
         }
         if (LpRmSmLm(-xb, -yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
         {
             path =
                 moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[7], -v, -u, .5 * pi, -t);
             Lmin = L;
         }

         if (LpRmSmRm(xb, yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
         {
             path =
                 moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[10], v, u, -.5 * pi, t);
             Lmin = L;
         }
         if (LpRmSmRm(-xb, yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
         {
             path =
                 moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[10], -v, -u, .5 * pi, -t);
             Lmin = L;
         }
         if (LpRmSmRm(xb, -yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
         {
             path =
                 moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[11], v, u, -.5 * pi, t);
             Lmin = L;
         }
         if (LpRmSmRm(-xb, -yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
             path =
                 moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[11], -v, -u, .5 * pi, -t);
     }
     // formula 8.11 *** TYPO IN PAPER ***
     inline bool LpRmSLmRp(double x, double y, double phi, double &t, double &u, double &v)
     {
         double xi = x + sin(phi), eta = y - 1. - cos(phi), rho, theta;
         polar(xi, eta, rho, theta);
         if (rho >= 2.)
         {
             u = 4. - sqrt(rho * rho - 4.);
             if (u <= ZERO)
             {
                 t = mod2pi(atan2((4 - u) * xi - 2 * eta, -2 * xi + (u - 4) * eta));
                 v = mod2pi(t - phi);
                 assert(fabs(4 * sin(t) - 2 * cos(t) - u * sin(t) - sin(phi) - x) < RS_EPS);
                 assert(fabs(-4 * cos(t) - 2 * sin(t) + u * cos(t) + cos(phi) + 1 - y) < RS_EPS);
                 assert(fabs(mod2pi(t - v - phi)) < RS_EPS);
                 return t >= -ZERO && v >= -ZERO;
             }
         }
         return false;
     }
     void CCSCC(double x, double y, double phi, moveit::core::PlanarJointModel::ReedsSheppPath &path)
     {
         double t, u, v, Lmin = path.length() - pi, L;
         if (LpRmSLmRp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[16], t, -.5 * pi, u,
                                                         -.5 * pi, v);
             Lmin = L;
         }
         if (LpRmSLmRp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[16], -t, .5 * pi, -u,
                                                         .5 * pi, -v);
             Lmin = L;
         }
         if (LpRmSLmRp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // reflect
         {
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[17], t, -.5 * pi, u,
                                                         -.5 * pi, v);
             Lmin = L;
         }
         if (LpRmSLmRp(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))  // timeflip + reflect
             path = moveit::core::PlanarJointModel::ReedsSheppPath(moveit::core::PlanarJointModel::reedsSheppPathType[17], -t, .5 * pi, -u,
                                                         .5 * pi, -v);
     }
     moveit::core::PlanarJointModel::ReedsSheppPath reedsShepp(double x, double y, double phi)
     {
         moveit::core::PlanarJointModel::ReedsSheppPath path;
         CSC(x, y, phi, path);
         CCC(x, y, phi, path);
         CCCC(x, y, phi, path);
         CCSC(x, y, phi, path);
         CCSCC(x, y, phi, path);
         return path;
     }
 }
namespace moveit
{
namespace core
{
PlanarJointModel::PlanarJointModel(const std::string& name) : JointModel(name), angular_distance_weight_(1.0), motion_model_(HOLONOMIC), min_translational_distance_(0.01)
{
  type_ = PLANAR;

  local_variable_names_.push_back("x");
  local_variable_names_.push_back("y");
  local_variable_names_.push_back("theta");
  for (int i = 0; i < 3; ++i)
  {
    variable_names_.push_back(name_ + "/" + local_variable_names_[i]);
    variable_index_map_[variable_names_.back()] = i;
  }

  variable_bounds_.resize(3);
  variable_bounds_[0].position_bounded_ = true;
  variable_bounds_[1].position_bounded_ = true;
  variable_bounds_[2].position_bounded_ = false;

  variable_bounds_[0].min_position_ = -std::numeric_limits<double>::infinity();
  variable_bounds_[0].max_position_ = std::numeric_limits<double>::infinity();
  variable_bounds_[1].min_position_ = -std::numeric_limits<double>::infinity();
  variable_bounds_[1].max_position_ = std::numeric_limits<double>::infinity();
  variable_bounds_[2].min_position_ = -boost::math::constants::pi<double>();
  variable_bounds_[2].max_position_ = boost::math::constants::pi<double>();

  computeVariableBoundsMsg();
}

unsigned int PlanarJointModel::getStateSpaceDimension() const
{
  return 3;
}

double PlanarJointModel::getMaximumExtent(const Bounds& other_bounds) const
{
  double dx = other_bounds[0].max_position_ - other_bounds[0].min_position_;
  double dy = other_bounds[1].max_position_ - other_bounds[1].min_position_;
  return sqrt(dx * dx + dy * dy) + boost::math::constants::pi<double>() * angular_distance_weight_;
}

void PlanarJointModel::getVariableDefaultPositions(double* values, const Bounds& bounds) const
{
  for (unsigned int i = 0; i < 2; ++i)
  {
    // if zero is a valid value
    if (bounds[i].min_position_ <= 0.0 && bounds[i].max_position_ >= 0.0)
      values[i] = 0.0;
    else
      values[i] = (bounds[i].min_position_ + bounds[i].max_position_) / 2.0;
  }
  values[2] = 0.0;
}

void PlanarJointModel::getVariableRandomPositions(random_numbers::RandomNumberGenerator& rng, double* values,
                                                  const Bounds& bounds) const
{
  if (bounds[0].max_position_ >= std::numeric_limits<double>::infinity() ||
      bounds[0].min_position_ <= -std::numeric_limits<double>::infinity())
    values[0] = 0.0;
  else
    values[0] = rng.uniformReal(bounds[0].min_position_, bounds[0].max_position_);
  if (bounds[1].max_position_ >= std::numeric_limits<double>::infinity() ||
      bounds[1].min_position_ <= -std::numeric_limits<double>::infinity())
    values[1] = 0.0;
  else
    values[1] = rng.uniformReal(bounds[1].min_position_, bounds[1].max_position_);
  values[2] = rng.uniformReal(bounds[2].min_position_, bounds[2].max_position_);
}

void PlanarJointModel::getVariableRandomPositionsNearBy(random_numbers::RandomNumberGenerator& rng, double* values,
                                                        const Bounds& bounds, const double* near,
                                                        const double distance) const
{
  if (bounds[0].max_position_ >= std::numeric_limits<double>::infinity() ||
      bounds[0].min_position_ <= -std::numeric_limits<double>::infinity())
    values[0] = 0.0;
  else
    values[0] = rng.uniformReal(std::max(bounds[0].min_position_, near[0] - distance),
                                std::min(bounds[0].max_position_, near[0] + distance));
  if (bounds[1].max_position_ >= std::numeric_limits<double>::infinity() ||
      bounds[1].min_position_ <= -std::numeric_limits<double>::infinity())
    values[1] = 0.0;
  else
    values[1] = rng.uniformReal(std::max(bounds[1].min_position_, near[1] - distance),
                                std::min(bounds[1].max_position_, near[1] + distance));

  double da = angular_distance_weight_ * distance;
  // limit the sampling range to 2pi to work correctly even if the distance is very large
  if (da > boost::math::constants::pi<double>())
    da = boost::math::constants::pi<double>();
  values[2] = rng.uniformReal(near[2] - da, near[2] + da);
  normalizeRotation(values);
}

void PlanarJointModel::interpolate(const double* from, const double* to, const double t, double* state) const
{
    double dx = to[0] - from[0];
    double dy = to[1] - from[1];
    double c = cos(from[2]);
    double s = sin(from[2]);
    double x = c * dx + s * dy, y = -s * dx + c * dy, phi = to[2] - from[2];
    ReedsSheppPath path = reedsShepp(x / turning_radius_, y / turning_radius_, phi);

    double st[3];
    double seg = t * path.length(), phinew, v;
    st[0] = 0;
    st[1] = 0;
    st[2] = from[2];

    for (unsigned int i = 0; i < 5 && seg > 0; ++i)
    {
        if (path.length_[i] < 0)
        {
            v = std::max(-seg, path.length_[i]);
            seg += v;
        }
        else
        {
            v = std::min(seg, path.length_[i]);
            seg -= v;
        }
        phinew = st[2];
        switch (path.type_[i])
        {
            case RS_LEFT:
                st[0] = st[0] + sin(phinew + v) - sin(phinew);
                st[1] = st[1] -cos(phinew + v) + cos(phinew);
                st[2] = phinew + v;
                break;
            case RS_RIGHT:
                st[0] = st[0] - sin(phinew - v) + sin(phinew);
                st[1] = st[1] + cos(phinew - v) - cos(phinew);
                st[2] = phinew - v;
                break;
            case RS_STRAIGHT:
                st[0] = st[0] + v * cos(phinew);
                st[1] = st[1] + v * sin(phinew);
                break;
            case RS_NOP:
                break;
        }
    }
    state[0] = st[0] * turning_radius_ + from[0];
    state[1] = st[1] * turning_radius_ + from[1];
    state[2] = st[2];
    // Need to enforce bounds?
}

double PlanarJointModel::distance(const double* values1, const double* values2) const
{
  double dx = values2[0] - values1[0];
  double dy = values2[1] - values1[1];
  double c = cos(values1[2]);
  double s = sin(values1[2]);
  double x = c * dx + s * dy, y = -s * dx + c * dy, phi = values2[2] - values1[2];

  return turning_radius_ * reedsShepp(x / turning_radius_, y / turning_radius_, phi).length();
}

bool PlanarJointModel::satisfiesPositionBounds(const double* values, const Bounds& bounds, double margin) const
{
  for (unsigned int i = 0; i < 3; ++i)
    if (values[0] < bounds[0].min_position_ - margin || values[0] > bounds[0].max_position_ + margin)
      return false;
  return true;
}

bool PlanarJointModel::normalizeRotation(double* values) const
{
  double& v = values[2];
  if (v >= -boost::math::constants::pi<double>() && v <= boost::math::constants::pi<double>())
    return false;
  v = fmod(v, 2.0 * boost::math::constants::pi<double>());
  if (v < -boost::math::constants::pi<double>())
    v += 2.0 * boost::math::constants::pi<double>();
  else if (v > boost::math::constants::pi<double>())
    v -= 2.0 * boost::math::constants::pi<double>();
  return true;
}

bool PlanarJointModel::enforcePositionBounds(double* values, const Bounds& bounds) const
{
  bool result = normalizeRotation(values);
  for (unsigned int i = 0; i < 2; ++i)
  {
    if (values[i] < bounds[i].min_position_)
    {
      values[i] = bounds[i].min_position_;
      result = true;
    }
    else if (values[i] > bounds[i].max_position_)
    {
      values[i] = bounds[i].max_position_;
      result = true;
    }
  }
  return result;
}

void PlanarJointModel::computeTransform(const double* joint_values, Eigen::Isometry3d& transf) const
{
  transf = Eigen::Isometry3d(Eigen::Translation3d(joint_values[0], joint_values[1], 0.0) *
                             Eigen::AngleAxisd(joint_values[2], Eigen::Vector3d::UnitZ()));
}

void PlanarJointModel::computeVariablePositions(const Eigen::Isometry3d& transf, double* joint_values) const
{
  joint_values[0] = transf.translation().x();
  joint_values[1] = transf.translation().y();

  Eigen::Quaterniond q(transf.rotation());
  // taken from Bullet
  double s_squared = 1.0 - (q.w() * q.w());
  if (s_squared < 10.0 * std::numeric_limits<double>::epsilon())
    joint_values[2] = 0.0;
  else
  {
    double s = 1.0 / sqrt(s_squared);
    joint_values[2] = (acos(q.w()) * 2.0f) * (q.z() * s);
  }
}

const PlanarJointModel::ReedsSheppPathSegmentType PlanarJointModel::reedsSheppPathType[18][5] = {
    {RS_LEFT, RS_RIGHT, RS_LEFT, RS_NOP, RS_NOP},         // 0
    {RS_RIGHT, RS_LEFT, RS_RIGHT, RS_NOP, RS_NOP},        // 1
    {RS_LEFT, RS_RIGHT, RS_LEFT, RS_RIGHT, RS_NOP},       // 2
    {RS_RIGHT, RS_LEFT, RS_RIGHT, RS_LEFT, RS_NOP},       // 3
    {RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_NOP},    // 4
    {RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_NOP},   // 5
    {RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_LEFT, RS_NOP},    // 6
    {RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_RIGHT, RS_NOP},   // 7
    {RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP},   // 8
    {RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_NOP},    // 9
    {RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_LEFT, RS_NOP},   // 10
    {RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_RIGHT, RS_NOP},    // 11
    {RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_NOP, RS_NOP},     // 12
    {RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_NOP, RS_NOP},     // 13
    {RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_NOP, RS_NOP},      // 14
    {RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP, RS_NOP},    // 15
    {RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_RIGHT},  // 16
    {RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_LEFT}   // 17
};

PlanarJointModel::ReedsSheppPath::ReedsSheppPath(const ReedsSheppPathSegmentType *type, double t,
                                                                  double u, double v, double w, double x)
   : type_(type)
 {
     length_[0] = t;
     length_[1] = u;
     length_[2] = v;
     length_[3] = w;
     length_[4] = x;
     totalLength_ = fabs(t) + fabs(u) + fabs(v) + fabs(w) + fabs(x);
 }

}  // end of namespace core
}  // end of namespace moveit