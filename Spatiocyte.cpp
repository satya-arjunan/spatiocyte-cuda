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

#define NUM_BOX 1
#define BOX_MOL 16000000
int main() {
  std::cout << "start" << std::endl;
  Model model(NUM_BOX);
  std::cout << "1" << std::endl;
  Species A("A", BOX_MOL, 1e-12, model, model.get_compartment(),
            model.get_compartment().get_volume_species());
  std::cout << "2" << std::endl;
  model.initialize();
  std::cout << "3" << std::endl;
  A.populate();
  std::cout << "4" << std::endl;
  VisualLogger visual_logger(model);
  std::cout << "5" << std::endl;
  model.get_stepper().set_diffuser(A.get_diffuser());
  std::cout << "6" << std::endl;
  model.get_stepper().set_visual_logger(visual_logger);
  std::cout << "7" << std::endl;
  visual_logger.push_species(A);
  std::cout << "8" << std::endl;
  //visual_logger.push_species(model.get_compartment().get_surface_species());
  //visual_logger.push_species(model.get_compartment().get_volume_species());
  visual_logger.initialize();
  std::cout << "9" << std::endl;

  model.run(0.0001);
  std::cout << "10" << std::endl;
  boost::posix_time::ptime start(
      boost::posix_time::microsec_clock::universal_time()); 
  std::cout << "11" << std::endl;
  //model.run(0.1);
  model.run(0.005);
  std::cout << "12" << std::endl;
  boost::posix_time::ptime end(
      boost::posix_time::microsec_clock::universal_time());
  std::cout << "13" << std::endl;
  std::cout << "duration:" << (end-start)/1.0 << std::endl;
}
