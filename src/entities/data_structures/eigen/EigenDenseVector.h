#ifndef SPINNER_EIGENDENSEVECTOR_H
#define SPINNER_EIGENDENSEVECTOR_H

#include <Eigen/Dense>

#include "src/entities/data_structures/AbstractDenseVector.h"

namespace quantum::linear_algebra {
template <typename T>
class EigenDenseVector: public AbstractDenseVector {
  public:
    void resize(uint32_t new_size);
    void concatenate_with(const std::unique_ptr<AbstractDenseVector>& rhs) override;
    void add_identical_values(size_t number, double value) override;
    void subtract_minimum() override;
    void wise_exp() override;
    std::unique_ptr<AbstractDenseVector> multiply_by(double multiplier) const override;
    double dot(const std::unique_ptr<AbstractDenseVector>& rhs) const override;
    std::unique_ptr<AbstractDenseVector>
    element_wise_multiplication(const std::unique_ptr<AbstractDenseVector>& rhs) const override;
    uint32_t size() const override;
    double at(uint32_t i) const override;
    void print(std::ostream& os) const override;

    Eigen::Vector<T, -1>& modifyDenseVector();
    const Eigen::Vector<T, -1>& getDenseVector();

  private:
    Eigen::Vector<T, -1> vector_;

    // c-like pointers are necessary to avoid double-free error
    static const EigenDenseVector* downcast_ptr(const std::unique_ptr<AbstractDenseVector>& ptr);
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENDENSEVECTOR_H
