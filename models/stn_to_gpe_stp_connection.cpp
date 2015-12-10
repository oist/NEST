/*
 *  stn_to_gpe_stp_connection.cpp
 *
 *	Created by Osamu Shouno	2013/09/25
 *	Copyright 2011 Honda Research Institute Japan Co., Ltd. All rights reserved.
 *
 *
 *	Permission is granted to compile and modify
 *	this file for non-commercial use.
 *	See the file LICENSE for details.
 *
 */

#include "stn_to_gpe_stp_connection.h"
#include "network.h"
#include "connector_model.h"

namespace nest
{

  STNtoGPeSTPConnection::STNtoGPeSTPConnection() :
	ConnectionHetWD(),
	tau_fac(148.0),
	inc_fac(1.64),
	f0_(1.0),
    f_(1.0),
	f_bound(5.0),
	tau_dep(764.0),
    inc_dep(0.55),
    d0_(1.0),
    d_(1.0)
  { }

  void STNtoGPeSTPConnection::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);

	def<double_t>(d, "tau_fac", tau_fac);
	def<double_t>(d, "inc_fac", inc_fac);
	def<double_t>(d, "f0", f0_);
	def<double_t>(d, "f_bound", f_bound);
	def<double_t>(d, "f", f_);
	def<double_t>(d, "tau_dep", tau_dep);
	def<double_t>(d, "inc_dep", d_);
	def<double_t>(d, "d0", d0_);
	def<double_t>(d, "d", d_);	  
  }
  
  void STNtoGPeSTPConnection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);

	updateValue<double_t>(d, "tau_fac", tau_fac);
	updateValue<double_t>(d, "inc_fac", inc_fac);
	updateValue<double_t>(d, "f0", f0_);
	updateValue<double_t>(d, "f_bound", f_bound);
	updateValue<double_t>(d, "f", f_);
	updateValue<double_t>(d, "tau_dep", tau_dep);
	updateValue<double_t>(d, "inc_dep", inc_dep);
	updateValue<double_t>(d, "d0", d0_);
	updateValue<double_t>(d, "d", d_);
  }

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void STNtoGPeSTPConnection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);

	set_property<double_t>(d, "tau_facs", p, tau_fac);
    set_property<double_t>(d, "inc_facs", p, inc_fac);
    set_property<double_t>(d, "f0s", p, f0_);
	set_property<double_t>(d, "f_bounds", p, f_bound);
    set_property<double_t>(d, "fs", p, f_);
	  
    set_property<double_t>(d, "tau_deps", p, tau_dep);
    set_property<double_t>(d, "inc_deps", p, inc_dep);
    set_property<double_t>(d, "d0s", p, d0_);
    set_property<double_t>(d, "ds", p, d_);	  
  }

  void STNtoGPeSTPConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, "tau_facs");
    initialize_property_array(d, "inc_facs");  
    initialize_property_array(d, "f0s");
	initialize_property_array(d, "f_bounds");  
    initialize_property_array(d, "fs");
	initialize_property_array(d, "tau_deps");
	initialize_property_array(d, "inc_deps");  
	initialize_property_array(d, "d0s"); 
	initialize_property_array(d, "ds");  
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void STNtoGPeSTPConnection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<double_t>(d, "tau_facs", tau_fac);
    append_property<double_t>(d, "inc_facs", inc_fac);  
    append_property<double_t>(d, "f0s", f0_); 
	append_property<double_t>(d, "f_bounds", f_bound);  
    append_property<double_t>(d, "fs", f_);
	append_property<double_t>(d, "tau_deps", tau_dep);
	append_property<double_t>(d, "inc_deps", inc_dep);  
	append_property<double_t>(d, "d0s", d0_); 
	append_property<double_t>(d, "ds", d_);
  }

} 
