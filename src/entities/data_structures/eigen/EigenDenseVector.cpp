#include "EigenDenseVector.h"

namespace quantum::linear_algebra {

template <typename T>
void EigenDenseVector<T>::resize(uint32_t new_size) {
    vector_.resize(new_size);
}

template <typename T>
void EigenDenseVector<T>::concatenate_with(const std::unique_ptr<AbstractDenseVector>& rhs) {
    auto rhs_ = downcast_ptr(rhs);

    Eigen::Vector<T, -1> tmp(vector_.size() + rhs_->vector_.size());
    tmp << vector_, rhs_->vector_;

    vector_ = tmp;
}

template <typename T>
void EigenDenseVector<T>::add_identical_values(size_t number, double value) {
    Eigen::Vector<T, -1> tmp(number);
    tmp.fill(value);

    Eigen::Vector<T, -1> tmptmp(vector_.size() + number);
    tmptmp << vector_, tmp;
    vector_ = tmptmp;
}

template <typename T>
void EigenDenseVector<T>::subtract_minimum() {
    double minimum = vector_.minCoeff();
    vector_.array() -= minimum;
}

template <typename T>
std::unique_ptr<AbstractDenseVector> EigenDenseVector<T>::multiply_by(double multiplier) const {
    auto answer = std::make_unique<EigenDenseVector>();
    answer->resize(vector_.size());

    answer->vector_.array() = vector_.array() * multiplier;

    return answer;
}

template <typename T>
void EigenDenseVector<T>::wise_exp() {
    vector_.array() = exp(vector_.array());
}

template <typename T>
double EigenDenseVector<T>::dot(const std::unique_ptr<AbstractDenseVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    return vector_.dot(rhs_->vector_);
}

template <typename T>
std::unique_ptr<AbstractDenseVector> EigenDenseVector<T>::element_wise_multiplication(
    const std::unique_ptr<AbstractDenseVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    auto answer = std::make_unique<EigenDenseVector>();
    answer->resize(vector_.size());

    answer->vector_.array() = vector_.array() * rhs_->vector_.array();

    return answer;
}

template <typename T>
uint32_t EigenDenseVector<T>::size() const {
    return vector_.size();
}

template <typename T>
double EigenDenseVector<T>::at(uint32_t i) const {
    return vector_(i);
}

template <typename T>
void EigenDenseVector<T>::print(std::ostream& os) const {
    os << vector_ << std::endl;
}

template <typename T>
Eigen::Vector<T, -1>& EigenDenseVector<T>::modifyDenseVector() {
    return vector_;
}

template <typename T>
const Eigen::Vector<T, -1>& EigenDenseVector<T>::getDenseVector() {
    return vector_;
}

template <typename T>
const EigenDenseVector<T>*
EigenDenseVector<T>::downcast_ptr(const std::unique_ptr<AbstractDenseVector>& ptr) {
    auto answer = dynamic_cast<const EigenDenseVector<T>*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

template class EigenDenseVector<double>;
template class EigenDenseVector<float>;
}  // namespace quantum::linear_algebra