//
// This source code resides at www.jaeckel.org/LetsBeRational.7z .
//
=========================================================================================
 Copyright © 2013-2024 Peter Jäckel.
 
 Permission to use, copy, modify, and distribute this software is freely granted,
 provided that this notice is preserved.

 WARRANTY DISCLAIMER
 The Software is provided "as is" without warranty of any kind, either express or implied,
 including without limitation any implied warranties of condition, uninterrupted use,
 merchantability, fitness for a particular purpose, or non-infringement.
==========================================================================================

 /////////////////////////////////////////////////////////////////////////////////////////
 //
 // This is the reference implementation of the implied Black volatility computation algorithm
 // published in "Let's Be Rational" by Peter Jaeckel, © 2013-2024.
 //
 // See www.jaeckel.org/LetsBeRational.pdf for details of the mathematics
 // published in November 2013; Wilmott, pp 40-53, January 2015.
 //
 /////////////////////////////////////////////////////////////////////////////////////////

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

 Minor Update (2024-02-15)
 =========================

 + Implementation of erfinv(𝑒) using the internal branches of Φ⁻¹(·) avoiding catastrophic subtractive cancellation for small arguments 𝑒.
 + Consolidation of the exactly-at-the-money implied volatility calculation into √8·erfinv(𝛽) where 𝛽 is the normalised Black price.

 No changes in speed or accuracy (within the respective uncertainties, e.g., ±0.5nsec).

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

 NEW VERSION (2024-02-09)
 ========================

 Approximately 30% faster:
 ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
 My Windows test run
	letsberational_timing.exe       1048576 0 1e-08 0.0001 0.001 0.01 0.1 0.25 0.5 1 2 4 8 16 32 64 128 256 512
 (compiled with MSVC 1938) previously recorded ~326 nanoseconds, now shows ~219 nanoseconds on average (note that
 these numbers vary and we should round them to the leading 2 digits). The Linux (WSL) test on my machine used to
 show ~266 nanoseconds, now shows ~182 nanoseconds. On the new test invocation command line
    letsberational_timing	2097152	0 1e-14 1e-08 0.0001 0.001 0.01 0.05 0.1 0.25 0.5 1 2 4 8 16 32 64 128 256 512
 it even goes down to an average of 219 / 180 nsec for x64 / Linux-WSL.
 
 Improved accuracy. For some parameter combinations (|x/s|≈10, s/2>≈0.21), the accuracy of the Black function would
 perhaps have been only 1E-13 (in relative terms) whereas we now get (within the expected noise range) the attainable
 accuracy (around 1E-14 in this example). This was the actual intention of this version, the above mentioned speedup
 was more or less a side effect.
 
 The Black function evaluation partitioning has changed.
 We now have, with h:=x/s and t:=s/2 (s=σ·√T), for non-positive x values (i.e., all h<0) :-
 
  η := -13; τ = 2·DBL_EPSILON^(1/16) ≈ 0.21
  
  Region I. [evaluated 
    h < η && t < (τ+½)+η-h
	
  Region II.
    η ≤ h ≤ 0 && t < τ + (½/|η|)·|h|

  Region III. otherwise.
  
 See Figure-6-rev-1517.pdf for reference.
  
 Region I. is evaluated with and asymptotic expansion of the 'scaled normalised Black' function defined in analogy to
 the 'scaled' complementary error function as
   bx  :=  b / ( ∂(b(x,s)/∂s )
 as before, only that now the expansion order is adjust according to the distance from its upper (anti-diagonally aligned)
 border. This was the result of extensive tests as to what expansion order gives the target accuracy of DBL_EPSILON/2≈1.1E-16.
 See the logic in function asymptotic_expansion_of_scaled_normalised_black() for details - a lookup in a static threshold table
 determines the number of expansion terms used. This provided significant speedup in this region since the band where a higher
 order is needed is in fact quite slim.
 
 Region II. is as before evaluated via a small-t expansion to sixth/seventh order in t². Only odd order terms appear in the
 expansion and a total of seven terms are used. The value for Y'(h) = 1+h·Y(h) that was the dominant limit for accuracy is
 now evaluated via Nonlinear-Remez optimized minimax rational function approximations in two branches over h which guarantees
 accuracy better than DBL_EPSILON throughout - no more subtractive cancellation for h->∞.
 
 In Region III. (which is the union of previous regions III. and IV.), we use Cody's erfc() and erfcx() functions with
 judicious choice in each term to minimise the number of exponential function invocations, see function
 normalised_black_with_optimal_use_of_codys_functions() for details.
 
 The corner area where all of these regions meet (previously the three-region meeting point between regions I., II., and IV.)
 was my concern for accuracy of the internal normalised_black() function. This new version has no issues.
 Please note that all accuracy comparisons (with multi-precision packages such as mpfr or arbitrary precision software such
 as Maxima or Matlab) should be done under consideration of what is attainable in the evaluation of b(x,s). Even if we ignore
 its dependency on x, and just consider it a function of s, a first order error propagation analysis still suggests that we
 cannot expect relative accuracy much beyond
   ( 1+|(s·∂b(x,s)/∂s)/b(x,s)| ) · ε
 where ε is the machine's floating point precision.

 Another speedup was the analytical contraction of b_c defined in the main document LetsBeRational.pdf in equation (4.2) to
    b_c  =  exp(x/2) / 2 · [ 1 - erfcx(√|x|) ]
 where we remind the reader that we always have x≤0. To avoid loss of precision due to subtractive cancellation, the term
   1 - erfcx(√|x|)
 is also evaluated via a Nonlinear-Remez optimized minimax rational function approximation where appropriate. Here, again,
 the intention was preservation of accuracy and the speed was the side effect.

 Further, a small number of redundant intermediate function call levels have been removed, plus the removal of the enforcement
 of bracket preservation in the main two-iteration loop in lets_be_rational(). These changes gave very small speed improvements.
 
 Additional speedups were achieved by plain analytical simplification of some expressions, and by specialisation of the univariate
 functions b_l(x) and b_u(x) defined in equations (4.7) and (4.8), adding Remez-optimized minimax rational function approximations
 in branch sections.
 
 Also running: specialisation of the exact case F≡K, i.e., x=ln(F/K)≡0 since this is in fact a not-so-uncommon real world use case.
 The analytical solution here is σ=-2·Φ⁻¹(½-𝛽/2) where 𝛽 is the normalised Black price. The crux of this solution is that for small
 Black prices 𝛽, we evaluate Φ⁻¹(y) for y≈½ which is precisely where the inverse cumulative normal function suffers catastrophic
 loss of relative accuracy. Obviously, one can bypass this issue by the aid of expansions, and for this purpose, I implemented a
 Remez-oprimised rational function approximation of order (6,6) in 𝛽² that is accurate to 3.5E-17 on σ ∈ [0,2] (β ∈ [0,0.6826894921370859]).
 The implied volatility exactly ATM (x=0) is then evaluated on average in 13 nanoseconds (in Linux-WSL on my development laptop).
 
 LAST BUT NOT LEAST: I designed a faster but equally 'perfectly' accurate inverse cumulative normal function algorithm that I
 call "PJ-2024-Inverse-Normal". It is about 20%-30% faster in my tests than AS241 (my previous favourite) - feel free to report back
 to me what you see. I say 'perfectly' in the sense of 'within what is attainable'. For example, for x≈0, we have Φ(x)≈½. As a result,
 the relative accuracy of any round-loop calculation computed as Φ⁻¹(Φ(x))/x-1 for some small (but non-zero) value x is limited by the
 finite precision of p=Φ(x). As x → 0, p → ½ and the closest p can be to ½, but still be distinct from ½, is p = ½·(1±ε) where ε is the
 floating point granularity (DBL_EPSILON). Near x≈0, we have Φ(x) = ½ + x/√(2π) + O(x³). Hence, for any nonzero |x| < √(π/2)·DBL_EPSILON,
 [√(π/2)=1.2533..] we obtain Φ(x)=0.5, and thus for any such x we have Φ⁻¹(Φ(x))/x-1 = Φ⁻¹(0.5)/x-1 = 0/x-1 = -1 : 100% accuracy loss.
 The attainable accuracy can be computed by error propagation analysis. To first order, if x has relative accuracy ε, the attainable
 relative accuracy of f(x) should be (1+|x·f'(x)/f(x)|)·ε. As for the inverse function of f, for any x-y value pair where y=f(x),
 given y with relative accuracy ε, to first order, the inverse function f⁻¹(y) may attain the relative accuracy (1+|f(x)/(x·f'(x))|)·ε.
 
 I compared the accuracy of "PJ-2024-Inverse-Normal" against Algorithm AS241, 'The Percentage Points of the Normal Distribution', by M.J. Wichura
 see http://lib.stat.cmu.edu/apstat/241 or https://csg.sph.umich.edu/abecasis/gas_power_calculator/algorithm-as-241-the-percentage-points-of-the-normal-distribution.pdf
 and found both to be almost identical in their accuracy - perfect within attainable limits - the sole difference being that
 PJ-2024-Inverse-Normal is about 20%-30% faster on my hardware.
 
 For the source code of "PJ-2024-Inverse-Normal", please see normaldistribution.cpp. Feel free to use in other applications but please
 retain / copy to your target usage the copyright and warranty disclaimer just above the functions inverse_norm_cdf_for_low_probabilities()
 and inverse_norm_cdf().
 
 I tested the new rational function approximations for b_l(x)/b_max(x) and b_u(x)/b_max(x) and Φ⁻¹(p) in separate standalone tests based
 on high-precision floating point numerics (using mpreal/mpfr/mpir).
 
──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

 NEW VERSION (2023-12-18)
 ========================

 For a quick usage overview and accuracy illustration, please see spreadsheet LetsBeRational.xlsx.

 Changes:
 
   1) More careful handling of underflow and overflow, which enabled the removal of any bracketing checks (speed improvement).
   
   2) Better handling of small and very small prices (in the lowest branch where the objective function is ~ 1/ln(b)), by the aid of
      monolithic computation of b(s)/b'(s) [instead of both parts individually]. This improves speed and accuracy.
   
   3) Better handling of the limit σ → ∞ [b → b_max=exp(θ·x/2)] in the highest branch where the objective function is ~ ln(b_max-b)),
      by the aid of a direct evaluation of the 'complementary normalised Black function' b̄(x,s) :=  bₘₐₓ - b(x,s) which reduces to
                b̄(x,s)       =      eˣ𝄍²·Φ(-x/s-s/2)  +  e⁻ˣ𝄍²·Φ(x/s-s/2)
      and involves no subtractive cancellation. This improves speed and accuracy.
   
   4) The above also enabled the use of -Ofast with g++ 11.4.0.
      Under Linux [launched on Windows 11 via Microsoft's "Windows Subsystem for Linux" (WSL) layer], on a
	  12th Gen Intel(R) Core(TM) i5-12500H CPU (Asus Zenbook 14 - a plain laptop),
	  
	  *the calculation of a single implied volatility is now down to just under 270 nanoseconds.* [on average]
	  
	  That's more than 3 million implied volatilities *per second*.
	  =============================================================

      MSVC (Visual studio 2022, MSVC version 1936) speed is similar, though g++ (on Linux) outperforms MSVC.

   5) Prebuilt binaries for Win32 anc x64 (Excel usage, but also as a DLL in other applications - see lets_be_rational.h etc.),
      and Linux [built with Ubuntu 22.04.3 under Windows 11 via Microsoft's "Windows Subsystem for Linux" (WSL) layer].

   6) Source code and prebuilt packages (provided "as is" without warranty of any kind, etc.), i.e., add-ins (for x64 and Linux) for:-

      + GNU Octave (8.3.0)  [open source alternative to Matlab]: see files x64/letsberational.oct and Linux/letsberational.oct.
        Note that the *.oct files each depend on LetsBeRational.xll (or LetsBeRational.so) to be present in the same directory.


      + Python 3            both for direct usage directly in x64/ via 
                         letsberational.py → _letsberational.pyd → LetsBeRational.xll
        or in Linux via
                         letsberational.py → _letsberational.so → LetsBeRational.so
        or via the pre-built binary package files ('wheel's)
             x64/letsberational-1.0.1515-cp311-cp311-win_amd64.whl
        or
             Linux/letsberational-1.0.1515-cp310-cp310-linux_x86_64.whl

      + Gnuplot      supports a plugin API as of version 5.0. We give access via the source file
          lets_be_rational_in_gnuplot.cpp  (gnuplot API shim layer)
        and the front end load script
	      letsberational.gnuplot.
        There is also a demo script letsberational_demo.gnuplot which generates an example plot.
        See file LetsBeRational.xlsx for its graphical output.

   7) Better handling of extreme numbers for |x| in the lowest branch (0<b<b_l): when |x| > 190, we now use the Householder(4) method which
      is of 5th order (instead of 'only' Householder(3) which is of 4th order). The computational overhead is small (in relative terms).
      Likewise for |x| > 580 in the highest branch. This vanquishes the small residual inaccuracy I observed for such extreme values.
      Note that both of these changes have no impact on commercial usage since moneyness ratios of F/K < exp(-180) = 6.7E-79 are very
      rare, even in expansion (or otherwise analytically-induced) computations.

   NOTE: DBL_MAX is about 1.7977E+308, DBL_MIN is 2.2251E-308, and therefore |ln(DBL_MAX)|≈710 and |ln(DBL_MIN)|≈708.
   Hence, a natural limit for x = ln(F/K) is about 700 - above that, numbers become nigh impossible to be represented in
   standard IEEE 754 64 bit floating point numbers (53 bit mantissa).

==========================================================================================
 WARRANTY DISCLAIMER
 The Software is provided "as is" without warranty of any kind, either express or implied,
 including without limitation any implied warranties of condition, uninterrupted use,
 merchantability, fitness for a particular purpose, or non-infringement.
==========================================================================================

 The main function is

     extern "C" double ImpliedBlackVolatility(double price, double forward, double strike, double T, double q /* q=±1 for calls/puts */)

 in file lets_be_rational.cpp. Since the Black price function incurs massive subtractive cancellation errors when evaluated naively, we
 also provide a precise implementation for it. Its signature is

     extern "C" double Black(double forward, double strike, double sigma, double T, double q /* q=±1 for calls/puts */)

 and it is also implemented in file lets_be_rational.cpp.
 
 These functions can be directly called from any C/C++ code by the aid of the usual shared library inclusion/linking/loading techniques.

 In addition, in line with the text in "Let's Be Rational", we also provide

    // As defined in equation (2.3) in "Let's Be Rational".
    extern "C" double NormalisedBlack(double x, double s, double q /* q=±1 */);

 and

   // The input 'beta' is taken as a "normalised Black" price as defined in equation (2.3) in "Let's Be Rational".
   extern "C" double NormalisedImpliedBlackVolatility(double beta, double x, double q /* q=±1 */);

 This archive comes with pre-built binaries in the subdirectories Win32/ (32 bit Windows), x64/ (64 bit Windows), and Linux (64 bit).
 The Windows binaries were produced with Visual Studio 2022. The Linux binaries were generated with g++ 11.4.0 under Ubuntu 22.04.3
 launched on Windows 11 via Microsoft's "Windows Subsystem for Linux" (WSL) layer. I found WSL under Windows 11 to be surprisingly
 convenient since you can run all of genuine Linux, 64 bit Windows, and 32 bit Windows binaries directly from the command line within
 the WSL - no wine layer or similar required, and all Windows binaries run as native Windows (even though they were launched within the WSL).

 ------------------------------------------------------------------------------------------
  WARRANTY DISCLAIMER
  The Software is provided "as is" without warranty of any kind, either express or implied,
  including without limitation any implied warranties of condition, uninterrupted use,
  merchantability, fitness for a particular purpose, or non-infringement.
 ------------------------------------------------------------------------------------------


 To build with Microsoft Visual Studio 2022:
 ===========================================
  Load "LetsBeRational (Visual Studio 2022).sln"
 
 To build with GNU make/g++/etc:
 ===============================
 Use the provided Makefile, i.e.,
 
  make

 The Makefile is designed to build, by default, the underlying shared library LetsBeRational.xll/so, and,
 depending on what is found on your system, also:-
 
  -  a GNU Octave addin

  -  a Python 'wheel'
  
  in the subdirectory Win32/x64/Linux (whatever your platform).
  
 Cross-compilation on Linux for Win32/x64 is also somewhat supported.

 Kindly peruse the Makefile for other options and further details, if you would, please.


 EXCEL SUPPORT
 =============
 
 The library, as it is, compiles directly with a standard Excel 'XLL' interface. The above functions appear under the following synonyms in Excel:-

    Black   ---   Black
            ---   LetsBeRational.Black

    NormalisedBlack   ---   NormalisedBlack
                      ---   LetsBeRational.NormalisedBlack

    ImpliedBlackVolatility   ---   ImpliedBlackVolatility
                             ---   ImpliedVolatility
                             ---   LetsBeRational.ImpliedBlackVolatility
                             ---   LetsBeRational.ImpliedVolatility

    NormalisedImpliedBlackVolatility   ---   NormalisedImpliedBlackVolatility
                                       ---   NormalisedImpliedVolatility
                                       ---   LetsBeRational.NormalisedImpliedBlackVolatility

 There are some more functions, see LetsBeRational.xlsx, for example, or search the source code for all occurences of 'DECLARE_XL_FUNCTION'.
 All of the Excel worksheet functions will appear in the category "LetsBeRational".

 The XLL wrapper code automatically detects if the respectively used version of Excel provides the Excel-12 API or only the earlier Excel-4 API, and uses the higher available version.

 The XLL wrapper code (comprised by {main,XLFunctions,XLOper}.{h,cpp}) is fairly generic and can be used to build other XLLs, 32 bit or 64 bit. By default, without the definition of
 the preprocessor macro XL_CATEGORY_NAME, the code detects the name of the XLL at run time, and registers all provided Excel worksheet functions under that base name as the Excel
 "category" name. You can override this by giving a C/C++ preprocessor macro definition such as -DXL_CATEGORY_NAME=LetsBeRational or similar. Caution: if you *do not* specify
 -DXL_CATEGORY_NAME=<any_name_you_choose> at compile time, and then rename the XLL file, the Excel category of the provided functions changes in Excel, and any functions defined with
 the category as a prefix [e.g., LetsBeRational.Black(forward,strike,sigma,T,q)], then appear with a different effective name since the prefix changes, too! This, however, may be
 desirable since you may want to prevent anyone renaming the XLL anyway.


 MATLAB SUPPORT
 ===============
 According to https://de.mathworks.com/help/matlab/ref/loadlibrary.html, Matlab can load 'extern "C"' functions directly, e.g., by the instruction

   loadlibrary("LetsBeRational.xll","lets_be_rational.h")

 though I have not tested this due to not having access to Matlab.
 
 CAUTION: Matlab has its own implementation of "Let's Be Rational" according to https://uk.mathworks.com/help/finance/blsimpv.html:
	Volatility = blsimpv(Spot,Strike,Rate,Time,Value)
 The function blsimpv() offers two different calculation algorithms. Its documentation says that is uses "Let's Be Rational" by default.
 It also assumes that the option is a call option by default, unless specified otherwise via the addition of the optional arguments 'Class',{'put'}, e.g.,
    blsimpv(Spot,Strike,Rate,Time,Value,'Class',{'put'})
 To translate this into the use of ImpliedBlackVolatility(double price, double forward, double strike, double T, double q /* q=±1 for calls/puts */),
 we should have the equivalence
   blsimpv(spot,strike,r,T,value)  =  ImpliedBlackVolatility(value·exp(r·T),spot·exp(r·T),strike,T,1)
 and
   blsimpv(spot,strike,r,T,value,'Class',{'put'})  =  ImpliedBlackVolatility(value·exp(r·T),spot·exp(r·T),strike,T,-1)
 HOWEVER, there has recently been a report submitted at
   https://www.researchgate.net/publication/373160753_Tighter_bounds_for_implied_volatility_based_on_the_Dirac_delta_family_method
 in which the authors conclude that "Let's Be Rational" returns very large errors for certain (perfectly sensible) data sets via the
 Matlab built-in function blsimpv(). The alleged failures of "Let's Be Rational" are most easily identified in their figures 9 and 10
 on page 26. The authors also announced their paper on https://www.linkedin.com/feed/update/urn:li:activity:7097647500831903744 and,
 there, mention that their source code is available at https://github.com/cuizhyu/DiracDeltaIV. From figure 9 and their Matlab source
 code, we glean that the furthest out-of-the-money test case is for a call option with S0=100,r=0.03,sigma=0.2,T=0.1,K=180.
 This gives us the forward F=S0·exp(r·T) = 100.3004505, Black(F,K,sigma,T,+1) = 1.04392E-20 and the call option value is
	Black(F,K,sigma,T,+1)·exp(-r·T) = 1.04079E-20.
 "Let's Be Rational"'s ImpliedBlackVolatility(value·exp(r·T),S0·exp(r·T),K,T,1) returns *exactly* 0.2 for these data.
 Not surpringly, it also produces exact results for the various other alleged failures.
 
 Not having access to Matlab, and Octave (the GNU equivalent system) not having a version of blsimpv() that goes back to "Let's Be Rational",
 I have been unable to reproduce the miscalculations reported in the mentioned article.
 
 The only explanation is that there is a mishap in the Matlab invocation of "Let's Be Rational" inside their function blsimpv().

 Hence, anyone wishing to use "Let's Be Rational" from within Matlab may be well advised to check if blsimpv() is reliable,
 else load ImpliedBlackVolatility() via

   loadlibrary("LetsBeRational.xll","lets_be_rational.h")

 as mentioned above.

 GNU OCTAVE SUPPORT
 ==================
 GNU Octave can load DLL-provided functions from so called '*.oct' files if they are in the running Octave session's load path,
 including all their dependencies. Octave comes with a compilation utility called 'mkoctfile' to generate the respective oct-file.
 
 We provide two ways to have access to the "Let's Be Rational" algorithm in GNU Octave:
 
   1) A handwritten c++ file lets_be_rational_in_octave.cc. See its source code for details as to its compilation with 'mkoctfile',
      or use the precompiled Windows binary Release/x64/lets_be_rational_in_octave.oct. Note that you need to copy that file and 
	  its dependency LetsBeRational.xll into your Octave load path. Then, load by the simple Octave command line instruction
	     lets_be_rational_in_octave
	  This will load lets_be_rational_in_octave.oct which in turn will register Black(...) and ImpliedBlackVolatility(...) into your running Octave session.
	  
   2) A SWIG (see https://www.swig.org/) interface generation file lets_be_rational_in_octave_via_swig.i.
      See the file lets_be_rational_in_octave_via_swig.i how to generate lets_be_rational_in_octave_via_swig.cc from it,
	  and then how to compile lets_be_rational_in_octave_via_swig.oct from that using 'mkoctfile', or use the precompiled
	  Windows binary Release/x64/lets_be_rational_in_octave_via_swig.oct.  Note that you need to copy that file and 
	  its dependency LetsBeRational.xll into your Octave load path. Then, load by the simple Octave command line instruction
	     lets_be_rational_in_octave_via_swig
	  This will load lets_be_rational_in_octave.oct which in turn will register Black(...) and ImpliedBlackVolatility(...) into
	  your running Octave session.
	  
	  Unfortunately, due to Octave continuously evolving in its functionality and its API, recent versions of Octave require SWIG 4.1.1 or higher,
	  and even then the two SWIG octave-specific library files files octruntime.swg and octrun.swg (usually residing in Lib\octave\ from the swig
	  installation directory) may have to be manually updated according to https://github.com/swig/swig/pull/2512/files.
      You can find accordingly patched versions of octruntime.swg and octrun.swg in the directory "patched SWIG files".
	  These extra hurdles swayed me not to include the SWIG generated Octave add-in lets_be_rational_in_octave_via_swig.oct in the Makefile's
	  Default target.

	I have been unable to find any performance difference between the two approaches. In general, I like using SWIG since it enables
	me to update API layers without the need for any manual (API, aka 'shim') code work, but in this simple case here it was educational
	to find out how to write Octave API code.

 All provided octave addin files were tested with GNU Octave 8.3.0 under 64 bit Windows, using the binaries in x64/, and with GNU Octave 6.4.0
 under Ubuntu (Windows 11 WSL as mentioned above). 


 PYTHON SUPPORT
 ==============
 
 Prebuilt binary wheels for 'letsberational' (as a wrapper around LetsBeRational.xll/so) are provided in the Win32/, x64/, and Linux/ subdirectories.
 
 You should also be able to test (or even use) the python 'letsberational' module without wheel installation via, e.g.,

   cd Linux; python -c 'import letsberational; print(dir(letsberational))'
 
 This should show you something like
 
    ['Black', 'CPUName', 'ImpliedBlackVolatility', 'NormalisedBlack', 'NormalisedImpliedBlackVolatility' ...]


==========================================================================================
 WARRANTY DISCLAIMER
 The Software is provided "as is" without warranty of any kind, either express or implied,
 including without limitation any implied warranties of condition, uninterrupted use,
 merchantability, fitness for a particular purpose, or non-infringement.
==========================================================================================
