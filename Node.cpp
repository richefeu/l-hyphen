#include "Node.hpp"
#include "null_size_t.hpp"

Node::Node(double x, double y)
    : pos(x, y), vel(), acc(), force(), mass(1.0), ictrl(null_size_t), prevNode(null_size_t), nextNode(null_size_t),
      kr(0.0), mz(0.0), mz_max(0.0) {}

void Node::init(double kr_, double mz_max_, size_t prev, size_t next) {
  prevNode = prev;
  nextNode = next;
  kr = kr_;
  mz_max = mz_max_;
}