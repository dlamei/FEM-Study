#include <iostream>
#include "utils.h"

void print(int i) {
    std::cout << i << std::endl;
}

int main() {


    defer(print(0));
    defer(print(1));
    defer(print(2));
    defer(print(3));

    assert(1 == 2, "this is false");

    return 0;
}
