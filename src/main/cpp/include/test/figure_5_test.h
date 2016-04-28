#ifndef TEST_FIGURE_5_TEST_H_
#define TEST_FIGURE_5_TEST_H_

#include <array>
#include <memory>

#include "param.h"
#include "java_util.h"

#include "cells/cell.h"
#include "cells/cell_factory.h"
#include "local_biology/cell_element.h"
#include "local_biology/neurite_element.h"
#include "local_biology/local_biology_module.h"
#include "simulation/ecm.h"
#include "simulation/scheduler.h"
#include "synapse/physical_bouton.h"
#include "synapse/biological_bouton.h"
#include "synapse/physical_spine.h"
#include "synapse/biological_spine.h"

namespace cx3d {

using cells::Cell;
using cells::CellFactory;
using local_biology::CellElement;
using local_biology::NeuriteElement;
using local_biology::LocalBiologyModule;
using simulation::ECM;
using simulation::Scheduler;
using synapse::PhysicalBouton;
using synapse::BiologicalBouton;
using synapse::PhysicalSpine;
using synapse::BiologicalSpine;

/**
 *  This class was used to produce Figure 5 of the paper
 * "A framework for modeling the growth and development of neurons and networks", Zubler & Douglas 2009.
 *
 * This class demonstrates the mechanical interactions between a chain of
 * PhysicalCylinders and a three PhysicalSpheres. The simulation starts by growing an axon,
 * then adds three overlapping spheres and pause. If un-paused, the mechanical interactions
 * push the objects apart.
 */
class Figure5Test {
 public:
  Figure5Test() {
  }

  void simulate(const std::shared_ptr<ECM>& ecm, const std::shared_ptr<JavaUtil2>& java) {
    // 1) Prepare the environment :
    //    eight extra PhysicalNodes :
    auto scheduler = Scheduler::getInstance(ecm);

    bool init_physical_node_movement_listener = true;  // hack to emulate original initialization order in Java (important for random number sequence)
    for (int i = 0; i < 18; i++) {
      double angle = 2 * Param::kPi * java->getRandomDouble1();
      std::array<double, 3> loc { 200 * java->sin(angle), 200 * java->cos(angle), -20 + 300 * java->getRandomDouble1() };
      if (init_physical_node_movement_listener) {
        init_physical_node_movement_listener = false;
        java->initPhysicalNodeMovementListener();
      }
      ecm->getPhysicalNodeInstance(loc);
    }

    // 2) creating a first cell, with a neurite going straight up.
    //    creating a 4-uple Cell-SomaElement-PhysicalSphere-SpatialOrganizerNode
    auto cell_a = CellFactory::getCellInstance( { 0, 0, -100 }, ecm);
    cell_a->setColorForAllPhysicalObjects(Param::kRedSolid);
    //    creating a single neurite
    auto ne = cell_a->getSomaElement()->extendNewNeurite(2.0, 0, 0);
    //    ne.setHasCytoskeletonMotor(false);
    //    elongating the neurite :
    std::array<double, 3> direction_up = { 0, 0, 1 };
    auto pc = ne->getPhysicalCylinder();
    for (int i = 0; i < 103; i++) {
      pc->movePointMass(300, direction_up);
      scheduler->simulateOneStep();
    }
    //    ecm.pause(3000);

    // 3) creating three additional spheres:
    auto cell_b = CellFactory::getCellInstance( { 10, 0, 0 }, ecm);
    auto ps_b = cell_b->getSomaElement()->getPhysicalSphere();

    ps_b->setMass(3);
    ps_b->setColor(Param::kYellowSolid);
    ps_b->setColor(Param::kYellow);
    auto cellC = CellFactory::getCellInstance( { -10, 0, 100 }, ecm);
    auto ps_c = cellC->getSomaElement()->getPhysicalSphere();

    ps_c->setMass(3);
    ps_c->setColor(Param::kYellowSolid);
    ps_c->setColor(Param::kYellow);
    auto cell_d = CellFactory::getCellInstance( { 10, 0, 160 }, ecm);
    auto ps_d = cell_d->getSomaElement()->getPhysicalSphere();

    ps_d->setColor(Param::kYellowSolid);
    ps_d->setColor(Param::kYellow);
    ps_d->setMass(2);

    // 4) setting a large diameter OR letting them grow
    bool growing = true;
    if (growing) {
      for (int i = 0; i < 30; i++) {
        ps_b->changeDiameter(400);
        ps_c->changeDiameter(300);
        ps_d->changeDiameter(200);
        scheduler->simulateOneStep();
      }
    } else {
      ps_b->setDiameter(140);
      ps_c->setDiameter(100);
      ps_d->setDiameter(50);
    }

    // 5) running the simulation slowly
    for (int i = 0; i < 1000; i++) {
      scheduler->simulateOneStep();
    }
  }
};

}  // namespace cx3d

#endif  // TEST_FIGURE_5_TEST_H_
