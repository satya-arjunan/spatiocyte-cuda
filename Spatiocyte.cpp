//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of the Spatiocyte package
//
//        Copyright (C) 2006-2009 Keio University
//        Copyright (C) 2010-2014 RIKEN
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// Spatiocyte is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// Spatiocyte is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with Spatiocyte -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
//

#include <iostream> 
#include <boost/date_time/posix_time/posix_time.hpp>
#include <Spatiocyte.hpp>
#include <Model.hpp>
#include <VisualLogger.hpp>

int main() {
  Model model(10000);
  Species A("A", 16, 1e-12, model, model.get_compartment(),
            model.get_compartment().get_volume_species());
  model.initialize();
  A.populate();
  VisualLogger visual_logger(model);
  model.get_stepper().set_diffuser(A.get_diffuser());
  model.get_stepper().set_visual_logger(visual_logger);
  visual_logger.push_species(A);
  //visual_logger.push_species(model.get_compartment().get_surface_species());
  //visual_logger.push_species(model.get_compartment().get_volume_species());
  visual_logger.initialize();

  model.run(0.0001);
  boost::posix_time::ptime start(
      boost::posix_time::microsec_clock::universal_time()); 
  //model.run(0.1);
  model.run(0.01);
  boost::posix_time::ptime end(
      boost::posix_time::microsec_clock::universal_time());
  std::cout << "duration:" << (end-start)/1.0 << std::endl;
}
