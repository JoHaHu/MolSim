// Copyright (c) 2005-2023 Code Synthesis.
//
// This program was generated by CodeSynthesis XSD, an XML Schema to
// C++ data binding compiler.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 2 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
//
// In addition, as a special exception, Code Synthesis gives permission
// to link this program with the Xerces-C++ library (or with modified
// versions of Xerces-C++ that use the same license as Xerces-C++), and
// distribute linked combinations including the two. You must obey the GNU
// General Public License version 2 in all respects for all of the code
// used other than Xerces-C++. If you modify this copy of the program, you
// may extend this exception to your version of the program, but you are
// not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// Furthermore, Code Synthesis makes a special exception for the Free/Libre
// and Open Source Software (FLOSS) which is described in the accompanying
// FLOSSE file.
//

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/pre.hxx>

#include "SimulationInputSchema.hxx"

// cuboidType
//

const cuboidType::lower_left_front_coordinate_type& cuboidType::
lower_left_front_coordinate () const
{
  return this->lower_left_front_coordinate_.get ();
}

cuboidType::lower_left_front_coordinate_type& cuboidType::
lower_left_front_coordinate ()
{
  return this->lower_left_front_coordinate_.get ();
}

void cuboidType::
lower_left_front_coordinate (const lower_left_front_coordinate_type& x)
{
  this->lower_left_front_coordinate_.set (x);
}

void cuboidType::
lower_left_front_coordinate (::std::unique_ptr< lower_left_front_coordinate_type > x)
{
  this->lower_left_front_coordinate_.set (std::move (x));
}

const cuboidType::dimensional_particle_numbers_type& cuboidType::
dimensional_particle_numbers () const
{
  return this->dimensional_particle_numbers_.get ();
}

cuboidType::dimensional_particle_numbers_type& cuboidType::
dimensional_particle_numbers ()
{
  return this->dimensional_particle_numbers_.get ();
}

void cuboidType::
dimensional_particle_numbers (const dimensional_particle_numbers_type& x)
{
  this->dimensional_particle_numbers_.set (x);
}

void cuboidType::
dimensional_particle_numbers (::std::unique_ptr< dimensional_particle_numbers_type > x)
{
  this->dimensional_particle_numbers_.set (std::move (x));
}

const cuboidType::distance_h_type& cuboidType::
distance_h () const
{
  return this->distance_h_.get ();
}

cuboidType::distance_h_type& cuboidType::
distance_h ()
{
  return this->distance_h_.get ();
}

void cuboidType::
distance_h (const distance_h_type& x)
{
  this->distance_h_.set (x);
}

const cuboidType::mass_m_type& cuboidType::
mass_m () const
{
  return this->mass_m_.get ();
}

cuboidType::mass_m_type& cuboidType::
mass_m ()
{
  return this->mass_m_.get ();
}

void cuboidType::
mass_m (const mass_m_type& x)
{
  this->mass_m_.set (x);
}

const cuboidType::initial_velocity_type& cuboidType::
initial_velocity () const
{
  return this->initial_velocity_.get ();
}

cuboidType::initial_velocity_type& cuboidType::
initial_velocity ()
{
  return this->initial_velocity_.get ();
}

void cuboidType::
initial_velocity (const initial_velocity_type& x)
{
  this->initial_velocity_.set (x);
}

void cuboidType::
initial_velocity (::std::unique_ptr< initial_velocity_type > x)
{
  this->initial_velocity_.set (std::move (x));
}


// intArray
//

const intArray::value_sequence& intArray::
value () const
{
  return this->value_;
}

intArray::value_sequence& intArray::
value ()
{
  return this->value_;
}

void intArray::
value (const value_sequence& s)
{
  this->value_ = s;
}


// Data
//

const Data::base_name_type& Data::
base_name () const
{
  return this->base_name_.get ();
}

Data::base_name_type& Data::
base_name ()
{
  return this->base_name_.get ();
}

void Data::
base_name (const base_name_type& x)
{
  this->base_name_.set (x);
}

void Data::
base_name (::std::unique_ptr< base_name_type > x)
{
  this->base_name_.set (std::move (x));
}

const Data::output_write_frequency_type& Data::
output_write_frequency () const
{
  return this->output_write_frequency_.get ();
}

Data::output_write_frequency_type& Data::
output_write_frequency ()
{
  return this->output_write_frequency_.get ();
}

void Data::
output_write_frequency (const output_write_frequency_type& x)
{
  this->output_write_frequency_.set (x);
}

void Data::
output_write_frequency (::std::unique_ptr< output_write_frequency_type > x)
{
  this->output_write_frequency_.set (std::move (x));
}

const Data::t_end_type& Data::
t_end () const
{
  return this->t_end_.get ();
}

Data::t_end_type& Data::
t_end ()
{
  return this->t_end_.get ();
}

void Data::
t_end (const t_end_type& x)
{
  this->t_end_.set (x);
}

const Data::delta_t_type& Data::
delta_t () const
{
  return this->delta_t_.get ();
}

Data::delta_t_type& Data::
delta_t ()
{
  return this->delta_t_.get ();
}

void Data::
delta_t (const delta_t_type& x)
{
  this->delta_t_.set (x);
}

const Data::epsilon_type& Data::
epsilon () const
{
  return this->epsilon_.get ();
}

Data::epsilon_type& Data::
epsilon ()
{
  return this->epsilon_.get ();
}

void Data::
epsilon (const epsilon_type& x)
{
  this->epsilon_.set (x);
}

const Data::sigma_type& Data::
sigma () const
{
  return this->sigma_.get ();
}

Data::sigma_type& Data::
sigma ()
{
  return this->sigma_.get ();
}

void Data::
sigma (const sigma_type& x)
{
  this->sigma_.set (x);
}

const Data::average_brownian_motion_type& Data::
average_brownian_motion () const
{
  return this->average_brownian_motion_.get ();
}

Data::average_brownian_motion_type& Data::
average_brownian_motion ()
{
  return this->average_brownian_motion_.get ();
}

void Data::
average_brownian_motion (const average_brownian_motion_type& x)
{
  this->average_brownian_motion_.set (x);
}

const Data::Cuboid_sequence& Data::
Cuboid () const
{
  return this->Cuboid_;
}

Data::Cuboid_sequence& Data::
Cuboid ()
{
  return this->Cuboid_;
}

void Data::
Cuboid (const Cuboid_sequence& s)
{
  this->Cuboid_ = s;
}


// output_write_frequency
//


#include <xsd/cxx/xml/dom/parsing-source.hxx>

// cuboidType
//

cuboidType::
cuboidType (const lower_left_front_coordinate_type& lower_left_front_coordinate,
            const dimensional_particle_numbers_type& dimensional_particle_numbers,
            const distance_h_type& distance_h,
            const mass_m_type& mass_m,
            const initial_velocity_type& initial_velocity)
: ::xml_schema::type (),
  lower_left_front_coordinate_ (lower_left_front_coordinate, this),
  dimensional_particle_numbers_ (dimensional_particle_numbers, this),
  distance_h_ (distance_h, this),
  mass_m_ (mass_m, this),
  initial_velocity_ (initial_velocity, this)
{
}

cuboidType::
cuboidType (::std::unique_ptr< lower_left_front_coordinate_type > lower_left_front_coordinate,
            ::std::unique_ptr< dimensional_particle_numbers_type > dimensional_particle_numbers,
            const distance_h_type& distance_h,
            const mass_m_type& mass_m,
            ::std::unique_ptr< initial_velocity_type > initial_velocity)
: ::xml_schema::type (),
  lower_left_front_coordinate_ (std::move (lower_left_front_coordinate), this),
  dimensional_particle_numbers_ (std::move (dimensional_particle_numbers), this),
  distance_h_ (distance_h, this),
  mass_m_ (mass_m, this),
  initial_velocity_ (std::move (initial_velocity), this)
{
}

cuboidType::
cuboidType (const cuboidType& x,
            ::xml_schema::flags f,
            ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  lower_left_front_coordinate_ (x.lower_left_front_coordinate_, f, this),
  dimensional_particle_numbers_ (x.dimensional_particle_numbers_, f, this),
  distance_h_ (x.distance_h_, f, this),
  mass_m_ (x.mass_m_, f, this),
  initial_velocity_ (x.initial_velocity_, f, this)
{
}

cuboidType::
cuboidType (const ::xercesc::DOMElement& e,
            ::xml_schema::flags f,
            ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  lower_left_front_coordinate_ (this),
  dimensional_particle_numbers_ (this),
  distance_h_ (this),
  mass_m_ (this),
  initial_velocity_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void cuboidType::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // lower_left_front_coordinate
    //
    if (n.name () == "lower_left_front_coordinate" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< lower_left_front_coordinate_type > r (
        lower_left_front_coordinate_traits::create (i, f, this));

      if (!lower_left_front_coordinate_.present ())
      {
        this->lower_left_front_coordinate_.set (::std::move (r));
        continue;
      }
    }

    // dimensional_particle_numbers
    //
    if (n.name () == "dimensional_particle_numbers" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< dimensional_particle_numbers_type > r (
        dimensional_particle_numbers_traits::create (i, f, this));

      if (!dimensional_particle_numbers_.present ())
      {
        this->dimensional_particle_numbers_.set (::std::move (r));
        continue;
      }
    }

    // distance_h
    //
    if (n.name () == "distance_h" && n.namespace_ ().empty ())
    {
      if (!distance_h_.present ())
      {
        this->distance_h_.set (distance_h_traits::create (i, f, this));
        continue;
      }
    }

    // mass_m
    //
    if (n.name () == "mass_m" && n.namespace_ ().empty ())
    {
      if (!mass_m_.present ())
      {
        this->mass_m_.set (mass_m_traits::create (i, f, this));
        continue;
      }
    }

    // initial_velocity
    //
    if (n.name () == "initial_velocity" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< initial_velocity_type > r (
        initial_velocity_traits::create (i, f, this));

      if (!initial_velocity_.present ())
      {
        this->initial_velocity_.set (::std::move (r));
        continue;
      }
    }

    break;
  }

  if (!lower_left_front_coordinate_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "lower_left_front_coordinate",
      "");
  }

  if (!dimensional_particle_numbers_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "dimensional_particle_numbers",
      "");
  }

  if (!distance_h_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "distance_h",
      "");
  }

  if (!mass_m_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "mass_m",
      "");
  }

  if (!initial_velocity_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "initial_velocity",
      "");
  }
}

cuboidType* cuboidType::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class cuboidType (*this, f, c);
}

cuboidType& cuboidType::
operator= (const cuboidType& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->lower_left_front_coordinate_ = x.lower_left_front_coordinate_;
    this->dimensional_particle_numbers_ = x.dimensional_particle_numbers_;
    this->distance_h_ = x.distance_h_;
    this->mass_m_ = x.mass_m_;
    this->initial_velocity_ = x.initial_velocity_;
  }

  return *this;
}

cuboidType::
~cuboidType ()
{
}

// intArray
//

intArray::
intArray ()
: ::xml_schema::type (),
  value_ (this)
{
}

intArray::
intArray (const intArray& x,
          ::xml_schema::flags f,
          ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  value_ (x.value_, f, this)
{
}

intArray::
intArray (const ::xercesc::DOMElement& e,
          ::xml_schema::flags f,
          ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  value_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void intArray::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // value
    //
    if (n.name () == "value" && n.namespace_ ().empty ())
    {
      this->value_.push_back (value_traits::create (i, f, this));
      continue;
    }

    break;
  }
}

intArray* intArray::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class intArray (*this, f, c);
}

intArray& intArray::
operator= (const intArray& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->value_ = x.value_;
  }

  return *this;
}

intArray::
~intArray ()
{
}

// Data
//

Data::
Data (const base_name_type& base_name,
      const output_write_frequency_type& output_write_frequency,
      const t_end_type& t_end,
      const delta_t_type& delta_t,
      const epsilon_type& epsilon,
      const sigma_type& sigma,
      const average_brownian_motion_type& average_brownian_motion)
: ::xml_schema::type (),
  base_name_ (base_name, this),
  output_write_frequency_ (output_write_frequency, this),
  t_end_ (t_end, this),
  delta_t_ (delta_t, this),
  epsilon_ (epsilon, this),
  sigma_ (sigma, this),
  average_brownian_motion_ (average_brownian_motion, this),
  Cuboid_ (this)
{
}

Data::
Data (const Data& x,
      ::xml_schema::flags f,
      ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  base_name_ (x.base_name_, f, this),
  output_write_frequency_ (x.output_write_frequency_, f, this),
  t_end_ (x.t_end_, f, this),
  delta_t_ (x.delta_t_, f, this),
  epsilon_ (x.epsilon_, f, this),
  sigma_ (x.sigma_, f, this),
  average_brownian_motion_ (x.average_brownian_motion_, f, this),
  Cuboid_ (x.Cuboid_, f, this)
{
}

Data::
Data (const ::xercesc::DOMElement& e,
      ::xml_schema::flags f,
      ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  base_name_ (this),
  output_write_frequency_ (this),
  t_end_ (this),
  delta_t_ (this),
  epsilon_ (this),
  sigma_ (this),
  average_brownian_motion_ (this),
  Cuboid_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void Data::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // base_name
    //
    if (n.name () == "base_name" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< base_name_type > r (
        base_name_traits::create (i, f, this));

      if (!base_name_.present ())
      {
        this->base_name_.set (::std::move (r));
        continue;
      }
    }

    // output_write_frequency
    //
    if (n.name () == "output_write_frequency" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< output_write_frequency_type > r (
        output_write_frequency_traits::create (i, f, this));

      if (!output_write_frequency_.present ())
      {
        this->output_write_frequency_.set (::std::move (r));
        continue;
      }
    }

    // t_end
    //
    if (n.name () == "t_end" && n.namespace_ ().empty ())
    {
      if (!t_end_.present ())
      {
        this->t_end_.set (t_end_traits::create (i, f, this));
        continue;
      }
    }

    // delta_t
    //
    if (n.name () == "delta_t" && n.namespace_ ().empty ())
    {
      if (!delta_t_.present ())
      {
        this->delta_t_.set (delta_t_traits::create (i, f, this));
        continue;
      }
    }

    // epsilon
    //
    if (n.name () == "epsilon" && n.namespace_ ().empty ())
    {
      if (!epsilon_.present ())
      {
        this->epsilon_.set (epsilon_traits::create (i, f, this));
        continue;
      }
    }

    // sigma
    //
    if (n.name () == "sigma" && n.namespace_ ().empty ())
    {
      if (!sigma_.present ())
      {
        this->sigma_.set (sigma_traits::create (i, f, this));
        continue;
      }
    }

    // average_brownian_motion
    //
    if (n.name () == "average_brownian_motion" && n.namespace_ ().empty ())
    {
      if (!average_brownian_motion_.present ())
      {
        this->average_brownian_motion_.set (average_brownian_motion_traits::create (i, f, this));
        continue;
      }
    }

    // Cuboid
    //
    if (n.name () == "Cuboid" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< Cuboid_type > r (
        Cuboid_traits::create (i, f, this));

      this->Cuboid_.push_back (::std::move (r));
      continue;
    }

    break;
  }

  if (!base_name_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "base_name",
      "");
  }

  if (!output_write_frequency_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "output_write_frequency",
      "");
  }

  if (!t_end_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "t_end",
      "");
  }

  if (!delta_t_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "delta_t",
      "");
  }

  if (!epsilon_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "epsilon",
      "");
  }

  if (!sigma_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "sigma",
      "");
  }

  if (!average_brownian_motion_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "average_brownian_motion",
      "");
  }
}

Data* Data::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class Data (*this, f, c);
}

Data& Data::
operator= (const Data& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->base_name_ = x.base_name_;
    this->output_write_frequency_ = x.output_write_frequency_;
    this->t_end_ = x.t_end_;
    this->delta_t_ = x.delta_t_;
    this->epsilon_ = x.epsilon_;
    this->sigma_ = x.sigma_;
    this->average_brownian_motion_ = x.average_brownian_motion_;
    this->Cuboid_ = x.Cuboid_;
  }

  return *this;
}

Data::
~Data ()
{
}

// output_write_frequency
//

output_write_frequency::
output_write_frequency (const ::xml_schema::positive_integer& _xsd_positive_integer_base)
: ::xsd::cxx::tree::fundamental_base< ::xml_schema::positive_integer, char, ::xml_schema::simple_type > (_xsd_positive_integer_base)
{
}

output_write_frequency::
output_write_frequency (const output_write_frequency& x,
                        ::xml_schema::flags f,
                        ::xml_schema::container* c)
: ::xsd::cxx::tree::fundamental_base< ::xml_schema::positive_integer, char, ::xml_schema::simple_type > (x, f, c)
{
}

output_write_frequency::
output_write_frequency (const ::xercesc::DOMElement& e,
                        ::xml_schema::flags f,
                        ::xml_schema::container* c)
: ::xsd::cxx::tree::fundamental_base< ::xml_schema::positive_integer, char, ::xml_schema::simple_type > (e, f, c)
{
}

output_write_frequency::
output_write_frequency (const ::xercesc::DOMAttr& a,
                        ::xml_schema::flags f,
                        ::xml_schema::container* c)
: ::xsd::cxx::tree::fundamental_base< ::xml_schema::positive_integer, char, ::xml_schema::simple_type > (a, f, c)
{
}

output_write_frequency::
output_write_frequency (const ::std::string& s,
                        const ::xercesc::DOMElement* e,
                        ::xml_schema::flags f,
                        ::xml_schema::container* c)
: ::xsd::cxx::tree::fundamental_base< ::xml_schema::positive_integer, char, ::xml_schema::simple_type > (s, e, f, c)
{
}

output_write_frequency* output_write_frequency::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class output_write_frequency (*this, f, c);
}

output_write_frequency::
~output_write_frequency ()
{
}

#include <istream>
#include <xsd/cxx/xml/sax/std-input-source.hxx>
#include <xsd/cxx/tree/error-handler.hxx>

::std::unique_ptr< ::Data >
Data_ (const ::std::string& u,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::tree::error_handler< char > h;

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

  return ::std::unique_ptr< ::Data > (
    ::Data_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::Data >
Data_ (const ::std::string& u,
       ::xml_schema::error_handler& h,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::Data > (
    ::Data_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::Data >
Data_ (const ::std::string& u,
       ::xercesc::DOMErrorHandler& h,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::Data > (
    ::Data_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::Data_ (isrc, f, p);
}

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       ::xml_schema::error_handler& h,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::Data_ (isrc, h, f, p);
}

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       ::xercesc::DOMErrorHandler& h,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::Data_ (isrc, h, f, p);
}

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       const ::std::string& sid,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::Data_ (isrc, f, p);
}

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       const ::std::string& sid,
       ::xml_schema::error_handler& h,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::Data_ (isrc, h, f, p);
}

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       const ::std::string& sid,
       ::xercesc::DOMErrorHandler& h,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::Data_ (isrc, h, f, p);
}

::std::unique_ptr< ::Data >
Data_ (::xercesc::InputSource& i,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xsd::cxx::tree::error_handler< char > h;

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

  return ::std::unique_ptr< ::Data > (
    ::Data_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::Data >
Data_ (::xercesc::InputSource& i,
       ::xml_schema::error_handler& h,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::Data > (
    ::Data_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::Data >
Data_ (::xercesc::InputSource& i,
       ::xercesc::DOMErrorHandler& h,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::Data > (
    ::Data_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::Data >
Data_ (const ::xercesc::DOMDocument& doc,
       ::xml_schema::flags f,
       const ::xml_schema::properties& p)
{
  if (f & ::xml_schema::flags::keep_dom)
  {
    ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
      static_cast< ::xercesc::DOMDocument* > (doc.cloneNode (true)));

    return ::std::unique_ptr< ::Data > (
      ::Data_ (
        std::move (d), f | ::xml_schema::flags::own_dom, p));
  }

  const ::xercesc::DOMElement& e (*doc.getDocumentElement ());
  const ::xsd::cxx::xml::qualified_name< char > n (
    ::xsd::cxx::xml::dom::name< char > (e));

  if (n.name () == "Data" &&
      n.namespace_ () == "")
  {
    ::std::unique_ptr< ::Data > r (
      ::xsd::cxx::tree::traits< ::Data, char >::create (
        e, f, 0));
    return r;
  }

  throw ::xsd::cxx::tree::unexpected_element < char > (
    n.name (),
    n.namespace_ (),
    "Data",
    "");
}

::std::unique_ptr< ::Data >
Data_ (::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d,
       ::xml_schema::flags f,
       const ::xml_schema::properties&)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > c (
    ((f & ::xml_schema::flags::keep_dom) &&
     !(f & ::xml_schema::flags::own_dom))
    ? static_cast< ::xercesc::DOMDocument* > (d->cloneNode (true))
    : 0);

  ::xercesc::DOMDocument& doc (c.get () ? *c : *d);
  const ::xercesc::DOMElement& e (*doc.getDocumentElement ());

  const ::xsd::cxx::xml::qualified_name< char > n (
    ::xsd::cxx::xml::dom::name< char > (e));

  if (f & ::xml_schema::flags::keep_dom)
    doc.setUserData (::xml_schema::dom::tree_node_key,
                     (c.get () ? &c : &d),
                     0);

  if (n.name () == "Data" &&
      n.namespace_ () == "")
  {
    ::std::unique_ptr< ::Data > r (
      ::xsd::cxx::tree::traits< ::Data, char >::create (
        e, f, 0));
    return r;
  }

  throw ::xsd::cxx::tree::unexpected_element < char > (
    n.name (),
    n.namespace_ (),
    "Data",
    "");
}

#include <ostream>
#include <xsd/cxx/tree/error-handler.hxx>
#include <xsd/cxx/xml/dom/serialization-source.hxx>

void
Data_ (::std::ostream& o,
       const ::Data& s,
       const ::xml_schema::namespace_infomap& m,
       const ::std::string& e,
       ::xml_schema::flags f)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0);

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::Data_ (s, m, f));

  ::xsd::cxx::tree::error_handler< char > h;

  ::xsd::cxx::xml::dom::ostream_format_target t (o);
  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    h.throw_if_failed< ::xsd::cxx::tree::serialization< char > > ();
  }
}

void
Data_ (::std::ostream& o,
       const ::Data& s,
       ::xml_schema::error_handler& h,
       const ::xml_schema::namespace_infomap& m,
       const ::std::string& e,
       ::xml_schema::flags f)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0);

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::Data_ (s, m, f));
  ::xsd::cxx::xml::dom::ostream_format_target t (o);
  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    throw ::xsd::cxx::tree::serialization< char > ();
  }
}

void
Data_ (::std::ostream& o,
       const ::Data& s,
       ::xercesc::DOMErrorHandler& h,
       const ::xml_schema::namespace_infomap& m,
       const ::std::string& e,
       ::xml_schema::flags f)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::Data_ (s, m, f));
  ::xsd::cxx::xml::dom::ostream_format_target t (o);
  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    throw ::xsd::cxx::tree::serialization< char > ();
  }
}

void
Data_ (::xercesc::XMLFormatTarget& t,
       const ::Data& s,
       const ::xml_schema::namespace_infomap& m,
       const ::std::string& e,
       ::xml_schema::flags f)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::Data_ (s, m, f));

  ::xsd::cxx::tree::error_handler< char > h;

  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    h.throw_if_failed< ::xsd::cxx::tree::serialization< char > > ();
  }
}

void
Data_ (::xercesc::XMLFormatTarget& t,
       const ::Data& s,
       ::xml_schema::error_handler& h,
       const ::xml_schema::namespace_infomap& m,
       const ::std::string& e,
       ::xml_schema::flags f)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::Data_ (s, m, f));
  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    throw ::xsd::cxx::tree::serialization< char > ();
  }
}

void
Data_ (::xercesc::XMLFormatTarget& t,
       const ::Data& s,
       ::xercesc::DOMErrorHandler& h,
       const ::xml_schema::namespace_infomap& m,
       const ::std::string& e,
       ::xml_schema::flags f)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::Data_ (s, m, f));
  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    throw ::xsd::cxx::tree::serialization< char > ();
  }
}

void
Data_ (::xercesc::DOMDocument& d,
       const ::Data& s,
       ::xml_schema::flags)
{
  ::xercesc::DOMElement& e (*d.getDocumentElement ());
  const ::xsd::cxx::xml::qualified_name< char > n (
    ::xsd::cxx::xml::dom::name< char > (e));

  if (n.name () == "Data" &&
      n.namespace_ () == "")
  {
    e << s;
  }
  else
  {
    throw ::xsd::cxx::tree::unexpected_element < char > (
      n.name (),
      n.namespace_ (),
      "Data",
      "");
  }
}

::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument >
Data_ (const ::Data& s,
       const ::xml_schema::namespace_infomap& m,
       ::xml_schema::flags f)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::serialize< char > (
      "Data",
      "",
      m, f));

  ::Data_ (*d, s, f);
  return d;
}

void
operator<< (::xercesc::DOMElement& e, const cuboidType& i)
{
  e << static_cast< const ::xml_schema::type& > (i);

  // lower_left_front_coordinate
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "lower_left_front_coordinate",
        e));

    s << i.lower_left_front_coordinate ();
  }

  // dimensional_particle_numbers
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "dimensional_particle_numbers",
        e));

    s << i.dimensional_particle_numbers ();
  }

  // distance_h
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "distance_h",
        e));

    s << ::xml_schema::as_double(i.distance_h ());
  }

  // mass_m
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "mass_m",
        e));

    s << ::xml_schema::as_double(i.mass_m ());
  }

  // initial_velocity
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "initial_velocity",
        e));

    s << i.initial_velocity ();
  }
}

void
operator<< (::xercesc::DOMElement& e, const intArray& i)
{
  e << static_cast< const ::xml_schema::type& > (i);

  // value
  //
  for (intArray::value_const_iterator
       b (i.value ().begin ()), n (i.value ().end ());
       b != n; ++b)
  {
    const intArray::value_type& x (*b);

    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "value",
        e));

    s << x;
  }
}

void
operator<< (::xercesc::DOMElement& e, const Data& i)
{
  e << static_cast< const ::xml_schema::type& > (i);

  // base_name
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "base_name",
        e));

    s << i.base_name ();
  }

  // output_write_frequency
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "output_write_frequency",
        e));

    s << i.output_write_frequency ();
  }

  // t_end
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "t_end",
        e));

    s << ::xml_schema::as_double(i.t_end ());
  }

  // delta_t
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "delta_t",
        e));

    s << ::xml_schema::as_double(i.delta_t ());
  }

  // epsilon
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "epsilon",
        e));

    s << ::xml_schema::as_double(i.epsilon ());
  }

  // sigma
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "sigma",
        e));

    s << ::xml_schema::as_double(i.sigma ());
  }

  // average_brownian_motion
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "average_brownian_motion",
        e));

    s << ::xml_schema::as_double(i.average_brownian_motion ());
  }

  // Cuboid
  //
  for (Data::Cuboid_const_iterator
       b (i.Cuboid ().begin ()), n (i.Cuboid ().end ());
       b != n; ++b)
  {
    const Data::Cuboid_type& x (*b);

    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "Cuboid",
        e));

    s << x;
  }
}

void
operator<< (::xercesc::DOMElement& e, const output_write_frequency& i)
{
  e << static_cast< const ::xsd::cxx::tree::fundamental_base< ::xml_schema::positive_integer, char, ::xml_schema::simple_type >& > (i);
}

void
operator<< (::xercesc::DOMAttr& a, const output_write_frequency& i)
{
  a << static_cast< const ::xsd::cxx::tree::fundamental_base< ::xml_schema::positive_integer, char, ::xml_schema::simple_type >& > (i);
}

void
operator<< (::xml_schema::list_stream& l,
            const output_write_frequency& i)
{
  l << static_cast< const ::xsd::cxx::tree::fundamental_base< ::xml_schema::positive_integer, char, ::xml_schema::simple_type >& > (i);
}

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

