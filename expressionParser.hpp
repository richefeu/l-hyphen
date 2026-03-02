#pragma once

#include "toofus-gate/exprtk/exprtk.hpp"

#include <iostream>

class ExpressionParser {
public:
  exprtk::symbol_table<double> symbol_table;
  exprtk::expression<double> expression;
  exprtk::parser<double> parser;

  ExpressionParser();
  void addConstant(std::string &name, double &value);
  void getValue(std::istream &is, double &value);
  void getValue(std::istream &is, int &value);
  void getValue(std::istream &is, size_t &value);
};