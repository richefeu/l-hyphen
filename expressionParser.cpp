#include "expressionParser.hpp"

ExpressionParser::ExpressionParser() {
  symbol_table.add_constants();
  expression.register_symbol_table(symbol_table);
}

void ExpressionParser::addConstant(std::string &name, double &value) {
  symbol_table.add_constant(name, value);
  expression.register_symbol_table(symbol_table);
}

void ExpressionParser::getValue(std::istream &is, double &value) {
  char C;
  is >> C; // Read the first character

  if (C == '$') {
    std::string expr_string;
    // Read until the closing '$'
    while (is >> std::noskipws >> C && C != '$') {
      expr_string += C;
    }

    // Now, expr_string contains the expression (without the closing '$').
    // We reset the skip-word-space behavior
    is >> std::skipws;

    // Parse expr_string and set 'value'
    parser.compile(expr_string, expression);
    value = static_cast<double>(expression.value());

    // Vérification si value est NaN
    if (std::isnan(value)) {
      std::cout << "Warning: Expression $" << expr_string << "$ resulted in NaN." << std::endl;
    }

  } else {
    // Put the character back and read coming token as a double
    is.putback(C);
    is >> value;
  }
}

void ExpressionParser::getValue(std::istream &is, int &value) {
  double dvalue;
  getValue(is, dvalue);
  value = static_cast<int>(std::round(dvalue));
}

void ExpressionParser::getValue(std::istream &is, size_t &value) {
  double dvalue;
  getValue(is, dvalue);
  value = static_cast<size_t>(std::round(dvalue));
}