/*
 *  gpe_stn_gabaa_connection3.cpp
 *
 *	Created by Osamu Shouno	2013/07/01
 *	Copyright 2013 Honda Research Institute Japan Co., Ltd. All rights reserved.
 *
 *
 *	Permission is granted to compile and modify
 *	this file for non-commercial use.
 *	See the file LICENSE for details.
 *
 */

#include "gpe_stn_gabaa_connection3.h"
#include "network.h"
#include "connector_model.h"

namespace nest
{

  GPeSTN_GABAA_Connection3::GPeSTN_GABAA_Connection3() :
    ConnectionHetWD(),
	n_contacts_(1.),
	n_contacts_mean(15.),
    d_(1.0),
	d0_(1.0),
	d_bound(0.0),
	tau_dep(20000.0),
    fd_(0.97),
	n_bound(1.0),
	myu_(1.0)
  { 
	
  }

  void GPeSTN_GABAA_Connection3::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);

	def<double_t>(d, "n_contacts", n_contacts_);
	def<double_t>(d, "n_contacts_mean", n_contacts_mean);  
	def<double_t>(d, "d", d_);
	def<double_t>(d, "d0", d0_);
	def<double_t>(d, "d_bound", d_bound); 
	def<double_t>(d, "n_bound", n_bound);
	def<double_t>(d, "tau_dep", tau_dep);
    def<double_t>(d, "fd", fd_);
	  
	def<double_t>(d, "myu", myu_);  
	
  }
  
  void GPeSTN_GABAA_Connection3::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
	updateValue<double_t>(d, "n_contacts", n_contacts_);
	updateValue<double_t>(d, "n_contacts_mean", n_contacts_mean);  
	updateValue<double_t>(d, "d", d_);
	updateValue<double_t>(d, "d0", d0_);
	updateValue<double_t>(d, "d_bound", d_bound);
	updateValue<double_t>(d, "n_bound", n_bound);
	updateValue<double_t>(d, "tau_dep", tau_dep);
    updateValue<double_t>(d, "fd", fd_);

	updateValue<double_t>(d, "myu", myu_);
  }

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void GPeSTN_GABAA_Connection3::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);

	set_property<double_t>(d, "n_contactss", p, n_contacts_);
	set_property<double_t>(d, "n_contacts_means", p, n_contacts_mean);  
	set_property<double_t>(d, "ds", p, d_);
	set_property<double_t>(d, "d0s", p, d0_);  
	set_property<double_t>(d, "d_bounds", p, d_bound);
	set_property<double_t>(d, "n_bounds", p, n_bound);  
	set_property<double_t>(d, "tau_deps", p, tau_dep);
    set_property<double_t>(d, "fds", p, fd_);
	set_property<double_t>(d, "myus", p, myu_);
  }

  void GPeSTN_GABAA_Connection3::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

	initialize_property_array(d, "n_contactss");
	initialize_property_array(d, "n_contacts_means");
	initialize_property_array(d, "ds"); 
	initialize_property_array(d, "d0s");
	initialize_property_array(d, "d_bounds");
	initialize_property_array(d, "n_bounds");  
	initialize_property_array(d, "tau_deps");
    initialize_property_array(d, "fds");  
	initialize_property_array(d, "myus");
	
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void GPeSTN_GABAA_Connection3::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

	append_property<double_t>(d, "n_contactss", n_contacts_);
	append_property<double_t>(d, "n_contacts_means", n_contacts_mean);
	  
	append_property<double_t>(d, "ds", d_);  
	append_property<double_t>(d, "d0s", d0_); 
	append_property<double_t>(d, "d_bounds", d_bound); 
	append_property<double_t>(d, "n_bounds", n_bound);  
    append_property<double_t>(d, "tau_deps", tau_dep);
	append_property<double_t>(d, "fds", fd_);  

	append_property<double_t>(d, "myus", myu_);
  }

} // of namespace bgnest
