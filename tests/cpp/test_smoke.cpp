// Build-plumbing check: proves the GoogleTest wiring (fetch, link, discovery) works.
#include <gtest/gtest.h>

TEST(Smoke, PointerIs64Bit) {
    EXPECT_EQ(sizeof(void*) * 8, 64u);
}
