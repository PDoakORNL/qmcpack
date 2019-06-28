//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Eric Neuscamman, eneuscamman@berkeley.edu, University of California, Berkeley
//
// File created by: Eric Neuscamman, eneuscamman@berkeley.edu, University of California, Berkeley
//////////////////////////////////////////////////////////////////////////////////////


/** @file ShortRangeCuspFunctor.h
 * @brief Functor designed to encode short-ranged part of nuclear cusp
 */
#ifndef QMCPLUSPLUS_SHORTRANGECUSPFUNCTOR_H
#define QMCPLUSPLUS_SHORTRANGECUSPFUNCTOR_H
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>
#include <stdexcept>
// #include <vector>
#include "OhmmsPETE/TinyVector.h"


namespace qmcplusplus
{

/*********************************************************************************************************
*
*  A functor intended for helping with the electron-nuclear (en) cusp at short range.
*
*  In QMCPACK, an en radial Jastrow is parameterized as exp(-U(r)), in which U(r) is the functor.
*  This class implements the functor
*
*    U(r) = -1.0 * exp(-r/R0) * ( A * R0 + sum_{k=0}^{N-1} B_k * (r/R0)^{k+2} / ( 1.0 + (r/R0)^{k+2} ) )
*
*  in which A, R0, and B_k are the variational parameters, although A is likely to be held fixed
*  as it enforces the en cusp condition ( dU/dr -> A as r -> 0 ).
*
*  The exp(-r/R0) keeps the functor very short range.  For example, initial testing on the LiH molecule
*  suggests that in that case, an R0 value of about 0.03 Bohr is desirable.  With this R0 value, r values
*  beyond just 0.5 Bohr already incur an exponential factor smaller than 10^{-7} and so the functor is
*  quite short ranged.
*
*  The B_k variables add structure to how the functor decays to zero with increasing r by expanding
*  in a set of sigmoidal functions.  The sigmoidal form (as opposed to a bare polynomial expansion)
*
*                          (r/R0)^{k+2} / ( 1.0 + (r/R0)^{k+2} )
*
*  was chosen so that they do not compete with the exponential decay and so help keep things short ranged.
*  Note that the lowest power of (r/R0) in this expansion is 2, so they do not affect the cusp condition.
*
**********************************************************************************************************/
template<class T>
struct ShortRangeCuspFunctor : public OptimizableFunctorBase
{

  //*****************************************************************************//
  //*******************************  Member Data ********************************//
  //*****************************************************************************//

  ///true, if A is optimizable
  bool Opt_A;
  ///true, if R0 is optimizable
  bool Opt_R0;
  ///true, if B variables are optimizable
  bool Opt_B;
  ///variable that controls the cusp
  real_type A;
  ///the soft cutoff distance that controls how short ranged the functor is
  real_type R0;
  ///variables that add detail through an expansion in sigmoidal functions
  std::vector<real_type> B;
  ///id of A
  std::string ID_A;
  ///id of R0
  std::string ID_R0;
  ///id of B
  std::string ID_B;

  //*****************************************************************************//
  //*****************************  Member Functions *****************************//
  //*****************************************************************************//

  ///default constructor
  ShortRangeCuspFunctor() : Opt_A(true), Opt_R0(true), Opt_B(true), A(1.0), R0(0.1), ID_A("A"), ID_R0("R0"), ID_B("B") { reset(); }

  ///default constructor
  void setCusp(real_type cusp)
  {
    //throw std::runtime_error("ShortRangeCuspFunctor::setCusp was called");
    A     = cusp;
    Opt_A = false;
    reset();
  }

  ///clone the functor
  OptimizableFunctorBase* makeClone() const { return new ShortRangeCuspFunctor(*this); }

  ///Implement the reset function, which was pure virtual in OptimizableFunctorBase, even though we don't need it
  void reset()
  {
    //cutoff_radius = 1.0e4; //some big range
  }

  ///compute U(r) at a particular value of r
  inline real_type evaluate(real_type r) const
  {

    // get the ratio of the distance and the soft cutoff distance
    const real_type s = r / R0;

    // sum up the sigmoidal function expansion
    real_type sig_sum = 0.0;
    real_type n = 2.0;
    real_type sn = s * s; // s^n
    for (int i = 0; i < B.size(); i++) {
      sig_sum += B[i] * sn / ( 1.0 + sn );
      sn *= s;  // update s^n
      n += 1.0; // update n
    }

    // return the functor's value
    return -1.0 * std::exp(-s) * ( A * R0 + sig_sum );

  }

  ///compute U(r), dU/dr, and d^2U/dr^2 at a particular value of r
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
  {

    // get the ratio of the distance and the soft cutoff distance
    const real_type s = r / R0;

    // get the exponential factor
    const real_type ex = std::exp(-s);

    // sum the terms related to the sigmoidal functions that are needed for U, dU, and d^2U
    real_type sig_sum_0 = 0.0; // contributes to U(r)
    real_type sig_sum_1 = 0.0; // contributes to dU/dr
    real_type sig_sum_2 = 0.0; // contributes to d2U/dr2
    real_type n = 2.0;
    real_type sn = s * s; // s^{n  }
    real_type snm1 = s;   // s^{n-1}
    real_type snm2 = 1.0; // s^{n-2}
    for (int i = 0; i < B.size(); i++) {
      const real_type d = 1.0 + sn;
      const real_type d2 = d * d;
      const real_type soverd = s / d;
      const real_type noverd2 = n / d2;
      sig_sum_0 += B[i] * sn / d;
      sig_sum_1 += B[i] * ( sn * sn + sn - n * snm1 ) / d2;
      sig_sum_2 += B[i] * snm2 * ( d * noverd2 * noverd2 * ( sn - 1.0 ) + noverd2 * ( 1.0 + 2.0 * s ) - soverd * s );
      snm2 = snm1; // update s^{n-2}
      snm1 = sn;   // update s^{n-1}
      sn *= s;     // update s^{n  }
      n += 1.0;    // update n
    }

    // get some useful ratios ( ex / R0 and ex / R0^2 )
    const real_type exoverR0  = ex / R0;
    const real_type exoverR02 = exoverR0 / R0;

    // evaluate dU/dr
    dudr = A * ex + exoverR0 * sig_sum_1;

    // evaluate d^2U/dr^2
    d2udr2 = -A * exoverR0 + exoverR02 * sig_sum_2;

    // return U(r)
    return -1.0 * ex * ( A * R0 + sig_sum_0 );

  }

  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
  {
    throw std::runtime_error("evaluate d3udr3 not implemented for ShortRangeCuspFunctor");
    return 0.0;
  }

  inline real_type evaluateV(const int iat,
                             const int iStart,
                             const int iEnd,
                             const T* restrict _distArray,
                             T* restrict distArrayCompressed) const
  {
    real_type sum(0);
    for (int idx = iStart; idx < iEnd; idx++)
      if (idx != iat)
        sum += evaluate(_distArray[idx]);
    return sum;
  }

  inline void evaluateVGL(const int iat,
                          const int iStart,
                          const int iEnd,
                          const T* distArray,
                          T* restrict valArray,
                          T* restrict gradArray,
                          T* restrict laplArray,
                          T* restrict distArrayCompressed,
                          int* restrict distIndices) const
  {
    for (int idx = iStart; idx < iEnd; idx++)
    {
      valArray[idx] = evaluate(distArray[idx], gradArray[idx], laplArray[idx]);
      gradArray[idx] /= distArray[idx];
    }
    if (iat >= iStart && iat < iEnd)
      valArray[iat] = gradArray[iat] = laplArray[iat] = T(0);
  }

  inline real_type f(real_type r)
  {
    if (r >= cutoff_radius)
      return 0.0;
    return evaluate(r);
  }

  inline real_type df(real_type r)
  {
    if (r >= cutoff_radius)
      return 0.0;
    real_type du, d2u;
    evaluate(r, du, d2u);
    return du;
  }

  /// compute derivatives with respect to the variational parameters
  inline bool evaluateDerivatives(real_type r, std::vector<TinyVector<real_type, 3>>& derivs)
  {

    const real_type s = r / R0;
    const real_type ex = std::exp(-s);

    int i = 0;

    if (Opt_A)
    {
      derivs[i][0] = -ex * R0;   //dU/da
      derivs[i][1] = ex;         //d(dU/da)/dr
      derivs[i][2] = -ex / R0;   //d^2 (dU/da)/dr^2
      ++i;
    }

    if (Opt_R0)
    {
      real_type n = 2.0;
      real_type sn = s * s;
      real_type snm1 = s;
      real_type snm2 = 1.0;
      real_type bpart_0 = 0.0; // contributes to dUdr0
      real_type bpart_1 = 0.0; // contributes to d(dU/dr0)/dr
      real_type bpart_2 = 0.0; // contributes to d^2(dU/dr0)/dr^2
      for (int j = 0; j < B.size(); j++) {
        const real_type d = 1.0 + sn;
        const real_type d2 = d * d;
        const real_type soverd = s / d;
        const real_type noverd2 = n / d2;
        bpart_0 += B[j] * sn * ( noverd2 - soverd );
        bpart_1 += B[j] * snm1 * ( d * noverd2 * noverd2 * ( 1.0 - sn ) - 2.0 * noverd2 * s + soverd * ( s - 1.0 ) );
        bpart_2 += B[j] * snm2 * (   s * ( 3.0 * s - 1.0 ) * noverd2 - s * ( s - 2.0 ) * soverd
                                   + ( n + n * sn * ( sn - 4.0 ) + ( 3.0 * s + 1.0 ) * ( sn * sn - 1.0 ) ) * noverd2 * noverd2 );
        snm2 = snm1;
        snm1 = sn;
        sn *= s;
        n += 1.0;
      }
      const real_type exoverR0 = ex / R0;
      const real_type exoverR02 = exoverR0 / R0;
      const real_type exoverR03 = exoverR02 / R0;
      const real_type apart_0 = -A * ex        * ( s + 1.0 );
      const real_type apart_1 =  A * exoverR0  * s;
      const real_type apart_2 = -A * exoverR02 * ( s - 1.0 );
      derivs[i][0] = apart_0 + exoverR0  * bpart_0;
      derivs[i][1] = apart_1 + exoverR02 * bpart_1;
      derivs[i][2] = apart_2 + exoverR03 * bpart_2;
      ++i;
    }

    if (Opt_B)
    {
      const real_type exoverR0 = ex / R0;
      const real_type exoverR02 = exoverR0 / R0;
      real_type n = 2.0;
      real_type sn = s * s;
      real_type snm1 = s;
      real_type snm2 = 1.0;
      for (int j = 0; j < B.size(); j++) {
        const real_type d = 1.0 + sn;
        const real_type d2 = d * d;
        const real_type soverd = s / d;
        const real_type noverd2 = n / d2;
        derivs[i][0] = -ex * sn / d;
        derivs[i][1] = exoverR0 * ( sn * sn + sn - n * snm1 ) / d2;
        derivs[i][2] = exoverR02 * snm2 * ( d * noverd2 * noverd2 * ( sn - 1.0 ) + noverd2 * ( 1.0 + 2.0 * s ) - soverd * s );
        snm2 = snm1;
        snm1 = sn;
        sn *= s;
        n += 1.0;
        ++i;
      }
    }

    return true;
  }

  /// compute derivatives with respect to the variational parameters
  inline bool evaluateDerivatives(real_type r, std::vector<real_type>& derivs)
  {

    const real_type s = r / R0;
    const real_type ex = std::exp(-s);

    int i = 0;

    if (Opt_A)
    {
      derivs[i] = -ex * R0;   //dU/da
      ++i;
    }

    if (Opt_R0)
    {
      real_type n = 2.0;
      real_type sn = s * s;
      //real_type snm1 = s;
      //real_type snm2 = 1.0;
      real_type bpart_0 = 0.0; // contributes to dUdr0
      //real_type bpart_1 = 0.0; // contributes to d(dU/dr0)/dr
      //real_type bpart_2 = 0.0; // contributes to d^2(dU/dr0)/dr^2
      for (int j = 0; j < B.size(); j++) {
        const real_type d = 1.0 + sn;
        const real_type d2 = d * d;
        const real_type soverd = s / d;
        const real_type noverd2 = n / d2;
        bpart_0 += B[j] * sn * ( noverd2 - soverd );
        //bpart_1 += B[j] * snm1 * ( d * noverd2 * noverd2 * ( 1.0 - sn ) - 2.0 * noverd2 * s + soverd * ( s - 1.0 ) );
        //bpart_2 += B[j] * snm2 * (   s * ( 3.0 * s - 1.0 ) * noverd2 - s * ( s - 2.0 ) * soverd
        //                           + ( n + n * sn * ( sn - 4.0 ) + ( 3.0 * s + 1.0 ) * ( sn * sn - 1.0 ) ) * noverd2 * noverd2 );
        //snm2 = snm1;
        //snm1 = sn;
        sn *= s;
        n += 1.0;
      }
      const real_type exoverR0 = ex / R0;
      //const real_type exoverR02 = exoverR0 / R0;
      //const real_type exoverR03 = exoverR02 / R0;
      const real_type apart_0 = -A * ex        * ( s + 1.0 );
      //const real_type apart_1 =  A * exoverR0  * s;
      //const real_type apart_2 = -A * exoverR02 * ( s - 1.0 );
      derivs[i] = apart_0 + exoverR0  * bpart_0;
      //derivs[i][1] = apart_1 + exoverR02 * bpart_1;
      //derivs[i][2] = apart_2 + exoverR03 * bpart_2;
      ++i;
    }

    if (Opt_B)
    {
      //const real_type exoverR0 = ex / R0;
      //const real_type exoverR02 = exoverR0 / R0;
      real_type n = 2.0;
      real_type sn = s * s;
      //real_type snm1 = s;
      //real_type snm2 = 1.0;
      for (int j = 0; j < B.size(); j++) {
        const real_type d = 1.0 + sn;
        //const real_type d2 = d * d;
        //const real_type soverd = s / d;
        //const real_type noverd2 = n / d2;
        derivs[i] = -ex * sn / d;
        //derivs[i][1] = exoverR0 * ( sn * sn + sn - n * snm1 ) / d2;
        //derivs[i][2] = exoverR02 * snm2 * ( d * noverd2 * noverd2 * ( sn - 1.0 ) + noverd2 * ( 1.0 + 2.0 * s ) - soverd * s );
        //snm2 = snm1;
        //snm1 = sn;
        sn *= s;
        n += 1.0;
        ++i;
      }
    }

//    int i       = 0;
//    real_type u = 1.0 / (1.0 + B * r);
//    if (Opt_A)
//    {
//      derivs[i] = r * u - 1 / B; //du/da
//      ++i;
//    }
//    if (Opt_B)
//    {
//      derivs[i] = -A * r * r * u * u + AoverB / B; //du/db
//    }

    return true;
  }

  bool put(xmlNodePtr cur)
  {
    ReportEngine PRE("ShortRangeCuspFunctor", "put(xmlNodePtr)");
    int nB = -1;
    real_type radius = -1.0;
    std::string name("0");
    std::string optA("yes");
    std::string optR0("yes");
    OhmmsAttributeSet rAttrib;
    rAttrib.add(nB, "size");
    rAttrib.add(radius, "rcut");
    rAttrib.add(radius, "cutoff");
    rAttrib.add(A, "cusp");
    rAttrib.add(R0, "r0");
    rAttrib.add(optA, "optA");
    rAttrib.add(optR0, "optR0");
    rAttrib.add(name, "elementType");
    rAttrib.put(cur);
    if (radius > 0.0)
      cutoff_radius = radius;
    if (nB < 0)
      PRE.error("Number of B paramters must be non-negative.  Did you forget to specify \"size\"?", true);
    if (name == "0")
      PRE.error("ShortRangeCuspFunctor did not find an acceptable name.", true);
    ID_A  = name + "_SRC_A";
    ID_B  = name + "_SRC_B";
    ID_R0 = name + "_SRC_R0";
    tolower(optA);
    tolower(optR0);
    Opt_A  = ( optA  == "yes" );
    Opt_R0 = ( optR0 == "yes" );
    if ( Opt_A )
      myVars.insert(ID_A, A, Opt_A, optimize::OTHER_P);
    if ( Opt_R0 )
      myVars.insert(ID_R0, R0, Opt_R0, optimize::OTHER_P);
    // Now read B coefficents
    xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
    while (xmlCoefs != NULL)
    {
      std::string cname((const char*)xmlCoefs->name);
      if (cname == "coefficients")
      {
        std::string type("0");
        //std::string id("0");
        std::string optB("yes");
        OhmmsAttributeSet cAttrib;
        //cAttrib.add(id, "id");
        cAttrib.add(type, "type");
        cAttrib.add(optB, "optimize");
        cAttrib.put(xmlCoefs);
        if (type != "Array")
          PRE.error("ShortRangeCuspFunctor expected B paramter array type to be \"Array\"", true);
        std::vector<real_type> params;
        putContent(params, xmlCoefs);
        if (params.size() != nB)
          PRE.error("ShortRangeCuspFunctor encountered a B parameter array that is the wrong length.", true);
        B = params;
        tolower(optB);
        Opt_B = ( optB == "yes" );
        if ( Opt_B ) {
          for (int i = 0; i < B.size(); i++)
          {
            std::stringstream sstr;
            //sstr << id << "_" << i;
            sstr << ID_B << "_" << i;
            myVars.insert(sstr.str(), B.at(i), Opt_B, optimize::OTHER_P);
          }
        }
        //int left_pad_space = 5;
        //app_log() << std::endl;
        //myVars.print(app_log(), left_pad_space, true);
        break;
      }
      xmlCoefs = xmlCoefs->next;
    }
    //app_summary() << "     Number of B parameters: " << nB << std::endl;
    app_summary() << "                   Cusp (A): " << A << std::endl;
    app_summary() << "                         R0: " << R0 << std::endl;
    app_summary() << "             B coefficients:";
    for (int i = 0; i < B.size(); i++)
      app_summary() << " " << B.at(i);
    app_summary() << std::endl;
    app_summary() << "              Cutoff radius: " << cutoff_radius << std::endl;
    app_summary() << "                      Opt_A: " << Opt_A << std::endl;
    app_summary() << "                     Opt_R0: " << Opt_R0 << std::endl;
    app_summary() << "                      Opt_B: " << Opt_B << std::endl;
    reset();
    return true;
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
    //myVars.print(std::cout);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    //myVars.print(std::cout);
  }

  void resetParameters(const opt_variables_type& active)
  {
    if (myVars.size())
    {
      int ia = myVars.where(0);
      if (ia > -1)
      {
        int i = 0;
        if (Opt_A)
          A = std::real( myVars[i++] = active[ia++] );
        if (Opt_R0)
          R0 = std::real( myVars[i++] = active[ia++] );
        for (int j = 0; Opt_B && j < B.size(); j++)
          B.at(j) = std::real( myVars[i++] = active[ia++] );
      }
      reset();
    }
  }

};

} // namespace qmcplusplus
#endif
