//  boost auto_timers_construction.cpp  ------------------------------------------------//

//  Copyright Beman Dawes 2007, 2011

//  Distributed under the Boost Software License, Version 1.0.
//  See http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/timer for documentation.

//--------------------------------------------------------------------------------------//

//  These constructors are in a separate file so that this translation unit will
//  not be linked in except when one of the constructors is actually used. This
//  is important since header <iostream> is required, and it incurs the cost of
//  the standard stream objects even if they are not used.

//--------------------------------------------------------------------------------------//

// define BOOST_TIMER_SOURCE so that <boost/timer/config.hpp> knows
// the library is being built (possibly exporting rather than importing code)
#ifndef BOOST_TIMER_SOURCE
# define BOOST_TIMER_SOURCE
#endif

#include <boost/timer/timer.hpp>
#include <iostream>

namespace
{
  // CAUTION: must be identical to same constant in cpu_timer.cpp
  const std::string default_fmt(" %ws wall, %us user + %ss system = %ts CPU (%p%)\n");
}

namespace boost
{
  namespace timer
  {

    // GEH
    // modify to remove warning in R package
    // since using std::cout is not allowed
    std::ostream out( NULL );

    auto_cpu_timer::auto_cpu_timer(short places)                                  // #1
      : m_places(places), m_os(&out), m_format(default_fmt) { start(); }

    auto_cpu_timer::auto_cpu_timer(short places, const std::string& format)       // #2
      : m_places(places), m_os(&out), m_format(format) { start(); }

    auto_cpu_timer::auto_cpu_timer(const std::string& format)                     // #3
      : m_places(default_places), m_os(&out), m_format(format) { start(); }

  } // namespace timer
} // namespace boost
