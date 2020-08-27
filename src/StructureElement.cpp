/**
 *  \file SequenceElement.cpp  \brief A special type of IMP::atom::Fragment
 *  particle that contains decorators for assigning sequence transformations
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/threading/StructureElement.h>
#include <IMP/statistics/internal/random_generator.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Selection.h>
#include <IMP/core/XYZ.h>
#include <IMP/algebra/vector_generators.h>

IMPTHREADING_BEGIN_NAMESPACE

FloatKey StructureElement::get_start_res_key() {
  static FloatKey k("start_res");
  return k;
}
FloatKey StructureElement::get_polarity_key() {
  static FloatKey k("polarity");
  return k;
}
FloatKey StructureElement::get_length_key() {
  static FloatKey k("length");
  return k;
}
FloatKey StructureElement::get_offset_key() {
  static FloatKey k("offset");
  return k;
}
StringKey StructureElement::get_chain_key() {
  StringKey k("chain_id");
  return k;
}

void StructureElement::show(std::ostream &out) const {
  out << get_offset_key();
}

IMPTHREADING_END_NAMESPACE
