/**
 *  \file SequenceElement.cpp  \brief A special type of IMP::atom::Fragment
 *  particle that contains decorators for assigning sequence transformations
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/StructureElementMolecule.h>
#include <IMP/statistics/internal/random_generator.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Selection.h>
#include <IMP/core/XYZ.h>
#include <IMP/algebra/vector_generators.h>

IMPTHREADING_BEGIN_NAMESPACE

FloatKey StructureElementMolecule::get_molname_index_key() {
  static FloatKey k("molname_index");
  return k;
}

FloatKey StructureElementMolecule::get_copy_key() {
  static FloatKey k("copy");
  return k;
}


void StructureElementMolecule::show(std::ostream &out) const {
  out << get_molname_index_key();
}

IMPTHREADING_END_NAMESPACE
