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

#ifndef SIMULATION_INPUT_SCHEMA_HXX
#define SIMULATION_INPUT_SCHEMA_HXX

#ifndef XSD_CXX11
#define XSD_CXX11
#endif

#ifndef XSD_USE_CHAR
#define XSD_USE_CHAR
#endif

#ifndef XSD_CXX_TREE_USE_CHAR
#define XSD_CXX_TREE_USE_CHAR
#endif

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/config.hxx>

/**
#if (LIBXSD_VERSION != 400002000000000L)
#error XSD runtime version mismatch
#endif
**/

#include <xsd/cxx/pre.hxx>

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/types.hxx>

#include <xsd/cxx/xml/error-handler.hxx>

#include <xsd/cxx/xml/dom/auto-ptr.hxx>

#include <xsd/cxx/tree/parsing.hxx>
#include <xsd/cxx/tree/parsing/byte.hxx>
#include <xsd/cxx/tree/parsing/unsigned-byte.hxx>
#include <xsd/cxx/tree/parsing/short.hxx>
#include <xsd/cxx/tree/parsing/unsigned-short.hxx>
#include <xsd/cxx/tree/parsing/int.hxx>
#include <xsd/cxx/tree/parsing/unsigned-int.hxx>
#include <xsd/cxx/tree/parsing/long.hxx>
#include <xsd/cxx/tree/parsing/unsigned-long.hxx>
#include <xsd/cxx/tree/parsing/boolean.hxx>
#include <xsd/cxx/tree/parsing/float.hxx>
#include <xsd/cxx/tree/parsing/double.hxx>
#include <xsd/cxx/tree/parsing/decimal.hxx>

#include <xsd/cxx/xml/dom/serialization-header.hxx>
#include <xsd/cxx/tree/serialization.hxx>
#include <xsd/cxx/tree/serialization/byte.hxx>
#include <xsd/cxx/tree/serialization/unsigned-byte.hxx>
#include <xsd/cxx/tree/serialization/short.hxx>
#include <xsd/cxx/tree/serialization/unsigned-short.hxx>
#include <xsd/cxx/tree/serialization/int.hxx>
#include <xsd/cxx/tree/serialization/unsigned-int.hxx>
#include <xsd/cxx/tree/serialization/long.hxx>
#include <xsd/cxx/tree/serialization/unsigned-long.hxx>
#include <xsd/cxx/tree/serialization/boolean.hxx>
#include <xsd/cxx/tree/serialization/float.hxx>
#include <xsd/cxx/tree/serialization/double.hxx>
#include <xsd/cxx/tree/serialization/decimal.hxx>

namespace xml_schema
{
  // anyType and anySimpleType.
  //
  typedef ::xsd::cxx::tree::type type;
  typedef ::xsd::cxx::tree::simple_type< char, type > simple_type;
  typedef ::xsd::cxx::tree::type container;

  // 8-bit
  //
  typedef signed char byte;
  typedef unsigned char unsigned_byte;

  // 16-bit
  //
  typedef short short_;
  typedef unsigned short unsigned_short;

  // 32-bit
  //
  typedef int int_;
  typedef unsigned int unsigned_int;

  // 64-bit
  //
  typedef long long long_;
  typedef unsigned long long unsigned_long;

  // Supposed to be arbitrary-length integral types.
  //
  typedef long long integer;
  typedef long long non_positive_integer;
  typedef unsigned long long non_negative_integer;
  typedef unsigned long long positive_integer;
  typedef long long negative_integer;

  // Boolean.
  //
  typedef bool boolean;

  // Floating-point types.
  //
  typedef float float_;
  typedef double double_;
  typedef double decimal;

  // String types.
  //
  typedef ::xsd::cxx::tree::string< char, simple_type > string;
  typedef ::xsd::cxx::tree::normalized_string< char, string > normalized_string;
  typedef ::xsd::cxx::tree::token< char, normalized_string > token;
  typedef ::xsd::cxx::tree::name< char, token > name;
  typedef ::xsd::cxx::tree::nmtoken< char, token > nmtoken;
  typedef ::xsd::cxx::tree::nmtokens< char, simple_type, nmtoken > nmtokens;
  typedef ::xsd::cxx::tree::ncname< char, name > ncname;
  typedef ::xsd::cxx::tree::language< char, token > language;

  // ID/IDREF.
  //
  typedef ::xsd::cxx::tree::id< char, ncname > id;
  typedef ::xsd::cxx::tree::idref< char, ncname, type > idref;
  typedef ::xsd::cxx::tree::idrefs< char, simple_type, idref > idrefs;

  // URI.
  //
  typedef ::xsd::cxx::tree::uri< char, simple_type > uri;

  // Qualified name.
  //
  typedef ::xsd::cxx::tree::qname< char, simple_type, uri, ncname > qname;

  // Binary.
  //
  typedef ::xsd::cxx::tree::buffer< char > buffer;
  typedef ::xsd::cxx::tree::base64_binary< char, simple_type > base64_binary;
  typedef ::xsd::cxx::tree::hex_binary< char, simple_type > hex_binary;

  // Date/time.
  //
  typedef ::xsd::cxx::tree::time_zone time_zone;
  typedef ::xsd::cxx::tree::date< char, simple_type > date;
  typedef ::xsd::cxx::tree::date_time< char, simple_type > date_time;
  typedef ::xsd::cxx::tree::duration< char, simple_type > duration;
  typedef ::xsd::cxx::tree::gday< char, simple_type > gday;
  typedef ::xsd::cxx::tree::gmonth< char, simple_type > gmonth;
  typedef ::xsd::cxx::tree::gmonth_day< char, simple_type > gmonth_day;
  typedef ::xsd::cxx::tree::gyear< char, simple_type > gyear;
  typedef ::xsd::cxx::tree::gyear_month< char, simple_type > gyear_month;
  typedef ::xsd::cxx::tree::time< char, simple_type > time;

  // Entity.
  //
  typedef ::xsd::cxx::tree::entity< char, ncname > entity;
  typedef ::xsd::cxx::tree::entities< char, simple_type, entity > entities;

  typedef ::xsd::cxx::tree::content_order content_order;
  // Namespace information and list stream. Used in
  // serialization functions.
  //
  typedef ::xsd::cxx::xml::dom::namespace_info< char > namespace_info;
  typedef ::xsd::cxx::xml::dom::namespace_infomap< char > namespace_infomap;
  typedef ::xsd::cxx::tree::list_stream< char > list_stream;
  typedef ::xsd::cxx::tree::as_double< double_ > as_double;
  typedef ::xsd::cxx::tree::as_decimal< decimal > as_decimal;
  typedef ::xsd::cxx::tree::facet facet;

  // Flags and properties.
  //
  typedef ::xsd::cxx::tree::flags flags;
  typedef ::xsd::cxx::tree::properties< char > properties;

  // Parsing/serialization diagnostics.
  //
  typedef ::xsd::cxx::tree::severity severity;
  typedef ::xsd::cxx::tree::error< char > error;
  typedef ::xsd::cxx::tree::diagnostics< char > diagnostics;

  // Exceptions.
  //
  typedef ::xsd::cxx::tree::exception< char > exception;
  typedef ::xsd::cxx::tree::bounds< char > bounds;
  typedef ::xsd::cxx::tree::duplicate_id< char > duplicate_id;
  typedef ::xsd::cxx::tree::parsing< char > parsing;
  typedef ::xsd::cxx::tree::expected_element< char > expected_element;
  typedef ::xsd::cxx::tree::unexpected_element< char > unexpected_element;
  typedef ::xsd::cxx::tree::expected_attribute< char > expected_attribute;
  typedef ::xsd::cxx::tree::unexpected_enumerator< char > unexpected_enumerator;
  typedef ::xsd::cxx::tree::expected_text_content< char > expected_text_content;
  typedef ::xsd::cxx::tree::no_prefix_mapping< char > no_prefix_mapping;
  typedef ::xsd::cxx::tree::serialization< char > serialization;

  // Error handler callback interface.
  //
  typedef ::xsd::cxx::xml::error_handler< char > error_handler;

  // DOM interaction.
  //
  namespace dom
  {
    // Automatic pointer for DOMDocument.
    //
    using ::xsd::cxx::xml::dom::unique_ptr;

#ifndef XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
#define XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
    // DOM user data key for back pointers to tree nodes.
    //
    const XMLCh* const tree_node_key = ::xsd::cxx::tree::user_data_keys::node;
#endif
  }
}

// Forward declarations.
//
class double_array;
class celestial_body;
class cuboid;
class disc;
class Data;
class header;
class gravity;
class lennard_jones;
class settings;
class cuboids;
class discs;

#include <memory>    // ::std::unique_ptr
#include <limits>    // std::numeric_limits
#include <algorithm> // std::binary_search
#include <utility>   // std::move

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/containers.hxx>
#include <xsd/cxx/tree/list.hxx>

#include <xsd/cxx/xml/dom/parsing-header.hxx>

class double_array: public ::xml_schema::type
{
  public:
  // value
  //
  typedef ::xml_schema::double_ value_type;
  typedef ::xsd::cxx::tree::sequence< value_type > value_sequence;
  typedef value_sequence::iterator value_iterator;
  typedef value_sequence::const_iterator value_const_iterator;
  typedef ::xsd::cxx::tree::traits< value_type, char, ::xsd::cxx::tree::schema_type::double_ > value_traits;

  const value_sequence&
  value () const;

  value_sequence&
  value ();

  void
  value (const value_sequence& s);

  // Constructors.
  //
  double_array ();

  double_array (const ::xercesc::DOMElement& e,
                ::xml_schema::flags f = 0,
                ::xml_schema::container* c = 0);

  double_array (const double_array& x,
                ::xml_schema::flags f = 0,
                ::xml_schema::container* c = 0);

  virtual double_array*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  double_array&
  operator= (const double_array& x);

  virtual 
  ~double_array ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  value_sequence value_;
};

class celestial_body: public ::xml_schema::type
{
  public:
  // coordinate
  //
  typedef ::double_array coordinate_type;
  typedef ::xsd::cxx::tree::traits< coordinate_type, char > coordinate_traits;

  const coordinate_type&
  coordinate () const;

  coordinate_type&
  coordinate ();

  void
  coordinate (const coordinate_type& x);

  void
  coordinate (::std::unique_ptr< coordinate_type > p);

  // velocity
  //
  typedef ::double_array velocity_type;
  typedef ::xsd::cxx::tree::traits< velocity_type, char > velocity_traits;

  const velocity_type&
  velocity () const;

  velocity_type&
  velocity ();

  void
  velocity (const velocity_type& x);

  void
  velocity (::std::unique_ptr< velocity_type > p);

  // mass
  //
  typedef ::xml_schema::double_ mass_type;
  typedef ::xsd::cxx::tree::traits< mass_type, char, ::xsd::cxx::tree::schema_type::double_ > mass_traits;

  const mass_type&
  mass () const;

  mass_type&
  mass ();

  void
  mass (const mass_type& x);

  // Constructors.
  //
  celestial_body (const coordinate_type&,
                  const velocity_type&,
                  const mass_type&);

  celestial_body (::std::unique_ptr< coordinate_type >,
                  ::std::unique_ptr< velocity_type >,
                  const mass_type&);

  celestial_body (const ::xercesc::DOMElement& e,
                  ::xml_schema::flags f = 0,
                  ::xml_schema::container* c = 0);

  celestial_body (const celestial_body& x,
                  ::xml_schema::flags f = 0,
                  ::xml_schema::container* c = 0);

  virtual celestial_body*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  celestial_body&
  operator= (const celestial_body& x);

  virtual 
  ~celestial_body ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< coordinate_type > coordinate_;
  ::xsd::cxx::tree::one< velocity_type > velocity_;
  ::xsd::cxx::tree::one< mass_type > mass_;
};

class cuboid: public ::xml_schema::type
{
  public:
  // coordinate
  //
  typedef ::double_array coordinate_type;
  typedef ::xsd::cxx::tree::traits< coordinate_type, char > coordinate_traits;

  const coordinate_type&
  coordinate () const;

  coordinate_type&
  coordinate ();

  void
  coordinate (const coordinate_type& x);

  void
  coordinate (::std::unique_ptr< coordinate_type > p);

  // particle_counts
  //
  typedef ::double_array particle_counts_type;
  typedef ::xsd::cxx::tree::traits< particle_counts_type, char > particle_counts_traits;

  const particle_counts_type&
  particle_counts () const;

  particle_counts_type&
  particle_counts ();

  void
  particle_counts (const particle_counts_type& x);

  void
  particle_counts (::std::unique_ptr< particle_counts_type > p);

  // velocity
  //
  typedef ::double_array velocity_type;
  typedef ::xsd::cxx::tree::traits< velocity_type, char > velocity_traits;

  const velocity_type&
  velocity () const;

  velocity_type&
  velocity ();

  void
  velocity (const velocity_type& x);

  void
  velocity (::std::unique_ptr< velocity_type > p);

  // Constructors.
  //
  cuboid (const coordinate_type&,
          const particle_counts_type&,
          const velocity_type&);

  cuboid (::std::unique_ptr< coordinate_type >,
          ::std::unique_ptr< particle_counts_type >,
          ::std::unique_ptr< velocity_type >);

  cuboid (const ::xercesc::DOMElement& e,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  cuboid (const cuboid& x,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  virtual cuboid*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  cuboid&
  operator= (const cuboid& x);

  virtual 
  ~cuboid ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< coordinate_type > coordinate_;
  ::xsd::cxx::tree::one< particle_counts_type > particle_counts_;
  ::xsd::cxx::tree::one< velocity_type > velocity_;
};

class disc: public ::xml_schema::type
{
  public:
  // coordinate
  //
  typedef ::double_array coordinate_type;
  typedef ::xsd::cxx::tree::traits< coordinate_type, char > coordinate_traits;

  const coordinate_type&
  coordinate () const;

  coordinate_type&
  coordinate ();

  void
  coordinate (const coordinate_type& x);

  void
  coordinate (::std::unique_ptr< coordinate_type > p);

  // velocity
  //
  typedef ::double_array velocity_type;
  typedef ::xsd::cxx::tree::traits< velocity_type, char > velocity_traits;

  const velocity_type&
  velocity () const;

  velocity_type&
  velocity ();

  void
  velocity (const velocity_type& x);

  void
  velocity (::std::unique_ptr< velocity_type > p);

  // radius
  //
  typedef ::xml_schema::int_ radius_type;
  typedef ::xsd::cxx::tree::traits< radius_type, char > radius_traits;

  const radius_type&
  radius () const;

  radius_type&
  radius ();

  void
  radius (const radius_type& x);

  // Constructors.
  //
  disc (const coordinate_type&,
        const velocity_type&,
        const radius_type&);

  disc (::std::unique_ptr< coordinate_type >,
        ::std::unique_ptr< velocity_type >,
        const radius_type&);

  disc (const ::xercesc::DOMElement& e,
        ::xml_schema::flags f = 0,
        ::xml_schema::container* c = 0);

  disc (const disc& x,
        ::xml_schema::flags f = 0,
        ::xml_schema::container* c = 0);

  virtual disc*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  disc&
  operator= (const disc& x);

  virtual 
  ~disc ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< coordinate_type > coordinate_;
  ::xsd::cxx::tree::one< velocity_type > velocity_;
  ::xsd::cxx::tree::one< radius_type > radius_;
};

class Data: public ::xml_schema::type
{
  public:
  // header
  //
  typedef ::header header_type;
  typedef ::xsd::cxx::tree::traits< header_type, char > header_traits;

  const header_type&
  header () const;

  header_type&
  header ();

  void
  header (const header_type& x);

  void
  header (::std::unique_ptr< header_type > p);

  // gravity
  //
  typedef ::gravity gravity_type;
  typedef ::xsd::cxx::tree::optional< gravity_type > gravity_optional;
  typedef ::xsd::cxx::tree::traits< gravity_type, char > gravity_traits;

  const gravity_optional&
  gravity () const;

  gravity_optional&
  gravity ();

  void
  gravity (const gravity_type& x);

  void
  gravity (const gravity_optional& x);

  void
  gravity (::std::unique_ptr< gravity_type > p);

  // lennard_jones
  //
  typedef ::lennard_jones lennard_jones_type;
  typedef ::xsd::cxx::tree::optional< lennard_jones_type > lennard_jones_optional;
  typedef ::xsd::cxx::tree::traits< lennard_jones_type, char > lennard_jones_traits;

  const lennard_jones_optional&
  lennard_jones () const;

  lennard_jones_optional&
  lennard_jones ();

  void
  lennard_jones (const lennard_jones_type& x);

  void
  lennard_jones (const lennard_jones_optional& x);

  void
  lennard_jones (::std::unique_ptr< lennard_jones_type > p);

  // Constructors.
  //
  Data (const header_type&);

  Data (::std::unique_ptr< header_type >);

  Data (const ::xercesc::DOMElement& e,
        ::xml_schema::flags f = 0,
        ::xml_schema::container* c = 0);

  Data (const Data& x,
        ::xml_schema::flags f = 0,
        ::xml_schema::container* c = 0);

  virtual Data*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  Data&
  operator= (const Data& x);

  virtual 
  ~Data ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< header_type > header_;
  gravity_optional gravity_;
  lennard_jones_optional lennard_jones_;
};

class header: public ::xml_schema::type
{
  public:
  // base_name
  //
  typedef ::xml_schema::string base_name_type;
  typedef ::xsd::cxx::tree::traits< base_name_type, char > base_name_traits;

  const base_name_type&
  base_name () const;

  base_name_type&
  base_name ();

  void
  base_name (const base_name_type& x);

  void
  base_name (::std::unique_ptr< base_name_type > p);

  // t_end
  //
  typedef ::xml_schema::double_ t_end_type;
  typedef ::xsd::cxx::tree::traits< t_end_type, char, ::xsd::cxx::tree::schema_type::double_ > t_end_traits;

  const t_end_type&
  t_end () const;

  t_end_type&
  t_end ();

  void
  t_end (const t_end_type& x);

  // output_frequency
  //
  typedef ::xml_schema::double_ output_frequency_type;
  typedef ::xsd::cxx::tree::traits< output_frequency_type, char, ::xsd::cxx::tree::schema_type::double_ > output_frequency_traits;

  const output_frequency_type&
  output_frequency () const;

  output_frequency_type&
  output_frequency ();

  void
  output_frequency (const output_frequency_type& x);

  // output_file_name
  //
  typedef ::xml_schema::string output_file_name_type;
  typedef ::xsd::cxx::tree::traits< output_file_name_type, char > output_file_name_traits;

  const output_file_name_type&
  output_file_name () const;

  output_file_name_type&
  output_file_name ();

  void
  output_file_name (const output_file_name_type& x);

  void
  output_file_name (::std::unique_ptr< output_file_name_type > p);

  // seed
  //
  typedef ::xml_schema::int_ seed_type;
  typedef ::xsd::cxx::tree::traits< seed_type, char > seed_traits;

  const seed_type&
  seed () const;

  seed_type&
  seed ();

  void
  seed (const seed_type& x);

  // Constructors.
  //
  header (const base_name_type&,
          const t_end_type&,
          const output_frequency_type&,
          const output_file_name_type&,
          const seed_type&);

  header (const ::xercesc::DOMElement& e,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  header (const header& x,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  virtual header*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  header&
  operator= (const header& x);

  virtual 
  ~header ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< base_name_type > base_name_;
  ::xsd::cxx::tree::one< t_end_type > t_end_;
  ::xsd::cxx::tree::one< output_frequency_type > output_frequency_;
  ::xsd::cxx::tree::one< output_file_name_type > output_file_name_;
  ::xsd::cxx::tree::one< seed_type > seed_;
};

class gravity: public ::xml_schema::type
{
  public:
  // total_bodies
  //
  typedef ::xml_schema::double_ total_bodies_type;
  typedef ::xsd::cxx::tree::traits< total_bodies_type, char, ::xsd::cxx::tree::schema_type::double_ > total_bodies_traits;

  const total_bodies_type&
  total_bodies () const;

  total_bodies_type&
  total_bodies ();

  void
  total_bodies (const total_bodies_type& x);

  // celestial_body
  //
  typedef ::celestial_body celestial_body_type;
  typedef ::xsd::cxx::tree::sequence< celestial_body_type > celestial_body_sequence;
  typedef celestial_body_sequence::iterator celestial_body_iterator;
  typedef celestial_body_sequence::const_iterator celestial_body_const_iterator;
  typedef ::xsd::cxx::tree::traits< celestial_body_type, char > celestial_body_traits;

  const celestial_body_sequence&
  celestial_body () const;

  celestial_body_sequence&
  celestial_body ();

  void
  celestial_body (const celestial_body_sequence& s);

  // Constructors.
  //
  gravity (const total_bodies_type&);

  gravity (const ::xercesc::DOMElement& e,
           ::xml_schema::flags f = 0,
           ::xml_schema::container* c = 0);

  gravity (const gravity& x,
           ::xml_schema::flags f = 0,
           ::xml_schema::container* c = 0);

  virtual gravity*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  gravity&
  operator= (const gravity& x);

  virtual 
  ~gravity ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< total_bodies_type > total_bodies_;
  celestial_body_sequence celestial_body_;
};

class lennard_jones: public ::xml_schema::type
{
  public:
  // settings
  //
  typedef ::settings settings_type;
  typedef ::xsd::cxx::tree::traits< settings_type, char > settings_traits;

  const settings_type&
  settings () const;

  settings_type&
  settings ();

  void
  settings (const settings_type& x);

  void
  settings (::std::unique_ptr< settings_type > p);

  // cuboids
  //
  typedef ::cuboids cuboids_type;
  typedef ::xsd::cxx::tree::optional< cuboids_type > cuboids_optional;
  typedef ::xsd::cxx::tree::traits< cuboids_type, char > cuboids_traits;

  const cuboids_optional&
  cuboids () const;

  cuboids_optional&
  cuboids ();

  void
  cuboids (const cuboids_type& x);

  void
  cuboids (const cuboids_optional& x);

  void
  cuboids (::std::unique_ptr< cuboids_type > p);

  // discs
  //
  typedef ::discs discs_type;
  typedef ::xsd::cxx::tree::optional< discs_type > discs_optional;
  typedef ::xsd::cxx::tree::traits< discs_type, char > discs_traits;

  const discs_optional&
  discs () const;

  discs_optional&
  discs ();

  void
  discs (const discs_type& x);

  void
  discs (const discs_optional& x);

  void
  discs (::std::unique_ptr< discs_type > p);

  // Constructors.
  //
  lennard_jones (const settings_type&);

  lennard_jones (::std::unique_ptr< settings_type >);

  lennard_jones (const ::xercesc::DOMElement& e,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

  lennard_jones (const lennard_jones& x,
                 ::xml_schema::flags f = 0,
                 ::xml_schema::container* c = 0);

  virtual lennard_jones*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  lennard_jones&
  operator= (const lennard_jones& x);

  virtual 
  ~lennard_jones ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< settings_type > settings_;
  cuboids_optional cuboids_;
  discs_optional discs_;
};

class settings: public ::xml_schema::type
{
  public:
  // delta_t
  //
  typedef ::xml_schema::double_ delta_t_type;
  typedef ::xsd::cxx::tree::traits< delta_t_type, char, ::xsd::cxx::tree::schema_type::double_ > delta_t_traits;

  const delta_t_type&
  delta_t () const;

  delta_t_type&
  delta_t ();

  void
  delta_t (const delta_t_type& x);

  // sigma
  //
  typedef ::xml_schema::double_ sigma_type;
  typedef ::xsd::cxx::tree::traits< sigma_type, char, ::xsd::cxx::tree::schema_type::double_ > sigma_traits;

  const sigma_type&
  sigma () const;

  sigma_type&
  sigma ();

  void
  sigma (const sigma_type& x);

  // epsilon
  //
  typedef ::xml_schema::double_ epsilon_type;
  typedef ::xsd::cxx::tree::traits< epsilon_type, char, ::xsd::cxx::tree::schema_type::double_ > epsilon_traits;

  const epsilon_type&
  epsilon () const;

  epsilon_type&
  epsilon ();

  void
  epsilon (const epsilon_type& x);

  // mass_m
  //
  typedef ::xml_schema::double_ mass_m_type;
  typedef ::xsd::cxx::tree::traits< mass_m_type, char, ::xsd::cxx::tree::schema_type::double_ > mass_m_traits;

  const mass_m_type&
  mass_m () const;

  mass_m_type&
  mass_m ();

  void
  mass_m (const mass_m_type& x);

  // distance_h
  //
  typedef ::xml_schema::double_ distance_h_type;
  typedef ::xsd::cxx::tree::traits< distance_h_type, char, ::xsd::cxx::tree::schema_type::double_ > distance_h_traits;

  const distance_h_type&
  distance_h () const;

  distance_h_type&
  distance_h ();

  void
  distance_h (const distance_h_type& x);

  // brown_motion
  //
  typedef ::xml_schema::double_ brown_motion_type;
  typedef ::xsd::cxx::tree::traits< brown_motion_type, char, ::xsd::cxx::tree::schema_type::double_ > brown_motion_traits;

  const brown_motion_type&
  brown_motion () const;

  brown_motion_type&
  brown_motion ();

  void
  brown_motion (const brown_motion_type& x);

  // Constructors.
  //
  settings (const delta_t_type&,
            const sigma_type&,
            const epsilon_type&,
            const mass_m_type&,
            const distance_h_type&,
            const brown_motion_type&);

  settings (const ::xercesc::DOMElement& e,
            ::xml_schema::flags f = 0,
            ::xml_schema::container* c = 0);

  settings (const settings& x,
            ::xml_schema::flags f = 0,
            ::xml_schema::container* c = 0);

  virtual settings*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  settings&
  operator= (const settings& x);

  virtual 
  ~settings ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< delta_t_type > delta_t_;
  ::xsd::cxx::tree::one< sigma_type > sigma_;
  ::xsd::cxx::tree::one< epsilon_type > epsilon_;
  ::xsd::cxx::tree::one< mass_m_type > mass_m_;
  ::xsd::cxx::tree::one< distance_h_type > distance_h_;
  ::xsd::cxx::tree::one< brown_motion_type > brown_motion_;
};

class cuboids: public ::xml_schema::type
{
  public:
  // Cuboid
  //
  typedef ::cuboid Cuboid_type;
  typedef ::xsd::cxx::tree::sequence< Cuboid_type > Cuboid_sequence;
  typedef Cuboid_sequence::iterator Cuboid_iterator;
  typedef Cuboid_sequence::const_iterator Cuboid_const_iterator;
  typedef ::xsd::cxx::tree::traits< Cuboid_type, char > Cuboid_traits;

  const Cuboid_sequence&
  Cuboid () const;

  Cuboid_sequence&
  Cuboid ();

  void
  Cuboid (const Cuboid_sequence& s);

  // Constructors.
  //
  cuboids ();

  cuboids (const ::xercesc::DOMElement& e,
           ::xml_schema::flags f = 0,
           ::xml_schema::container* c = 0);

  cuboids (const cuboids& x,
           ::xml_schema::flags f = 0,
           ::xml_schema::container* c = 0);

  virtual cuboids*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  cuboids&
  operator= (const cuboids& x);

  virtual 
  ~cuboids ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  Cuboid_sequence Cuboid_;
};

class discs: public ::xml_schema::type
{
  public:
  // Disc
  //
  typedef ::disc Disc_type;
  typedef ::xsd::cxx::tree::sequence< Disc_type > Disc_sequence;
  typedef Disc_sequence::iterator Disc_iterator;
  typedef Disc_sequence::const_iterator Disc_const_iterator;
  typedef ::xsd::cxx::tree::traits< Disc_type, char > Disc_traits;

  const Disc_sequence&
  Disc () const;

  Disc_sequence&
  Disc ();

  void
  Disc (const Disc_sequence& s);

  // Constructors.
  //
  discs ();

  discs (const ::xercesc::DOMElement& e,
         ::xml_schema::flags f = 0,
         ::xml_schema::container* c = 0);

  discs (const discs& x,
         ::xml_schema::flags f = 0,
         ::xml_schema::container* c = 0);

  virtual discs*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  discs&
  operator= (const discs& x);

  virtual 
  ~discs ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  Disc_sequence Disc_;
};

#include <iosfwd>

#include <xercesc/sax/InputSource.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>

// Parse a URI or a local file.
//

::std::unique_ptr< ::Data >
Data_ (const ::std::string& uri,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::Data >
Data_ (const ::std::string& uri,
       ::xml_schema::error_handler& eh,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::Data >
Data_ (const ::std::string& uri,
       ::xercesc::DOMErrorHandler& eh,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse std::istream.
//

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       ::xml_schema::error_handler& eh,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       ::xercesc::DOMErrorHandler& eh,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       const ::std::string& id,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       const ::std::string& id,
       ::xml_schema::error_handler& eh,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::Data >
Data_ (::std::istream& is,
       const ::std::string& id,
       ::xercesc::DOMErrorHandler& eh,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse xercesc::InputSource.
//

::std::unique_ptr< ::Data >
Data_ (::xercesc::InputSource& is,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::Data >
Data_ (::xercesc::InputSource& is,
       ::xml_schema::error_handler& eh,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::Data >
Data_ (::xercesc::InputSource& is,
       ::xercesc::DOMErrorHandler& eh,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse xercesc::DOMDocument.
//

::std::unique_ptr< ::Data >
Data_ (const ::xercesc::DOMDocument& d,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::unique_ptr< ::Data >
Data_ (::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d,
       ::xml_schema::flags f = 0,
       const ::xml_schema::properties& p = ::xml_schema::properties ());

#include <iosfwd>

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>
#include <xercesc/framework/XMLFormatter.hpp>

#include <xsd/cxx/xml/dom/auto-ptr.hxx>

// Serialize to std::ostream.
//

void
Data_ (::std::ostream& os,
       const ::Data& x, 
       const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
       const ::std::string& e = "UTF-8",
       ::xml_schema::flags f = 0);

void
Data_ (::std::ostream& os,
       const ::Data& x, 
       ::xml_schema::error_handler& eh,
       const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
       const ::std::string& e = "UTF-8",
       ::xml_schema::flags f = 0);

void
Data_ (::std::ostream& os,
       const ::Data& x, 
       ::xercesc::DOMErrorHandler& eh,
       const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
       const ::std::string& e = "UTF-8",
       ::xml_schema::flags f = 0);

// Serialize to xercesc::XMLFormatTarget.
//

void
Data_ (::xercesc::XMLFormatTarget& ft,
       const ::Data& x, 
       const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
       const ::std::string& e = "UTF-8",
       ::xml_schema::flags f = 0);

void
Data_ (::xercesc::XMLFormatTarget& ft,
       const ::Data& x, 
       ::xml_schema::error_handler& eh,
       const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
       const ::std::string& e = "UTF-8",
       ::xml_schema::flags f = 0);

void
Data_ (::xercesc::XMLFormatTarget& ft,
       const ::Data& x, 
       ::xercesc::DOMErrorHandler& eh,
       const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
       const ::std::string& e = "UTF-8",
       ::xml_schema::flags f = 0);

// Serialize to an existing xercesc::DOMDocument.
//

void
Data_ (::xercesc::DOMDocument& d,
       const ::Data& x,
       ::xml_schema::flags f = 0);

// Serialize to a new xercesc::DOMDocument.
//

::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument >
Data_ (const ::Data& x, 
       const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
       ::xml_schema::flags f = 0);

void
operator<< (::xercesc::DOMElement&, const double_array&);

void
operator<< (::xercesc::DOMElement&, const celestial_body&);

void
operator<< (::xercesc::DOMElement&, const cuboid&);

void
operator<< (::xercesc::DOMElement&, const disc&);

void
operator<< (::xercesc::DOMElement&, const Data&);

void
operator<< (::xercesc::DOMElement&, const header&);

void
operator<< (::xercesc::DOMElement&, const gravity&);

void
operator<< (::xercesc::DOMElement&, const lennard_jones&);

void
operator<< (::xercesc::DOMElement&, const settings&);

void
operator<< (::xercesc::DOMElement&, const cuboids&);

void
operator<< (::xercesc::DOMElement&, const discs&);

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

#endif // SIMULATION_INPUT_SCHEMA_HXX
