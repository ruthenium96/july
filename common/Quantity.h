#ifndef JULY_QUANTITYENUM_H
#define JULY_QUANTITYENUM_H

#include <string>
#include <ostream>

namespace common {
enum QuantityEnum {
    Energy,
    S_total_squared,
    };

[[nodiscard]] std::string get_quantity_name(common::QuantityEnum);
}

std::ostream &operator<<(std::ostream &os, const common::QuantityEnum &quantity);

#endif //JULY_QUANTITYENUM_H
