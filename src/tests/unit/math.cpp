#include <time.h>

#include <algorithm>
#include "libgcm/Math.hpp"
#include "libgcm/node/CalcNode.hpp"

#include <gtest/gtest.h>

using namespace gcm;

#define ITERATIONS 1000

#define MIN_ROOT 0.1
#define MAX_ROOT 1000

bool solveThirdOrderTest(gcm::real a1, gcm::real a2, gcm::real a3)
{
    gcm::real e[3];
    e[0] = a1;
    e[1] = a2;
    e[2] = a3;
    
    gcm::real a = - (e[0] + e[1] + e[2]);
    gcm::real b = e[0]*e[1] + e[0]*e[2] + e[1]*e[2];
    gcm::real c = - e[0] * e[1] * e[2];
    
    gcm::real r[3];
    solvePolynomialThirdOrder(a, b, c, r[0], r[1], r[2]);
    
    std::sort(std::begin(e), std::end(e));
    std::sort(std::begin(r), std::end(r));
    bool match = true;
    for( int i = 0; i < 3; i++ )
        if( fabs(e[i] - r[i]) > 0.025 * fabs(e[i] + r[i]) || isnan(r[i]) )
            match = false;
    return match;
};

TEST(Math, solvePolynomialThirdOrderThreeRoots)
{
    srand(0);

    gcm::real e[3];
    
    for (int k = 0; k < ITERATIONS; k++)
    {
        e[0] = MIN_ROOT + (MAX_ROOT-MIN_ROOT) * (double) rand() / RAND_MAX;
        e[1] = MIN_ROOT + (MAX_ROOT-MIN_ROOT) * (double) rand() / RAND_MAX;
        e[2] = MIN_ROOT + (MAX_ROOT-MIN_ROOT) * (double) rand() / RAND_MAX;
        
        ASSERT_TRUE( solveThirdOrderTest(e[0], e[1], e[2]) );
    }
};

TEST(Math, solvePolynomialThirdOrderTwoRoots)
{
    srand(0);

    gcm::real e[3];
    
    for (int k = 0; k < ITERATIONS; k++)
    {
        e[0] = MIN_ROOT + (MAX_ROOT-MIN_ROOT) * (double) rand() / RAND_MAX;
        e[1] = MIN_ROOT + (MAX_ROOT-MIN_ROOT) * (double) rand() / RAND_MAX;
        e[2] = e[1];
        
        ASSERT_TRUE( solveThirdOrderTest(e[0], e[1], e[2]) );
    }
}

TEST(Math, solvePolynomialThirdOrderOneRoot)
{
    srand(0);

    gcm::real e[3];
    
    for (int k = 0; k < ITERATIONS; k++)
    {
        e[0] = MIN_ROOT + (MAX_ROOT-MIN_ROOT) * (double) rand() / RAND_MAX;
        e[1] = e[0];
        e[2] = e[0];
        
        ASSERT_TRUE( solveThirdOrderTest(e[0], e[1], e[2]) );
    }
};

TEST(Math, solvePolynomialThirdOrderHydrostaticStress)
{
    srand(0);

    CalcNode node;
    
    for (int k = 0; k < ITERATIONS; k++)
    {
        gcm::real p = MIN_ROOT + (MAX_ROOT-MIN_ROOT) * (double) rand() / RAND_MAX;
        node.sxx = node.syy = node.szz = p;
        node.sxy = node.sxz = node.syz = 0;
        
        gcm::real a = - node.getJ1();
        gcm::real b = node.getJ2();
        gcm::real c = - node.getJ3();
        
        gcm::real r[3];
        solvePolynomialThirdOrder(a, b, c, r[0], r[1], r[2]);

        bool match = true;
        for( int i = 0; i < 3; i++ )
            if( fabs(p - r[i]) > 0.025 * fabs(p + r[i]) || isnan(r[i]) )
                match = false;
        
        ASSERT_TRUE( match );
    }
};

TEST(Math, vectorAdd)
{
    vector3r v1 = { 1, 2, 3 };
    vector3r v2 = { 3, 2, 1 };

    vector3r v = v1 + v2;

    ASSERT_FLOAT_EQ(v.x, 4);
    ASSERT_FLOAT_EQ(v.y, 4);
    ASSERT_FLOAT_EQ(v.z, 4);
}

TEST(Math, vectorSubtract)
{
    vector3r v1 = { 1, 2, 3 };
    vector3r v2 = { 3, 2, 1 };

    vector3r v = v1 - v2;

    ASSERT_FLOAT_EQ(v.x, -2);
    ASSERT_FLOAT_EQ(v.y, 0);
    ASSERT_FLOAT_EQ(v.z, 2);
}

TEST(Math, vectorScalarProduct)
{
    vector3r v1 = {1, 2, 3};
    vector3r v2 = {3, 2, 1};

    vector3r v = v1 - v2;

    ASSERT_FLOAT_EQ(v.x, -2);
    ASSERT_FLOAT_EQ(v.y, 0);
    ASSERT_FLOAT_EQ(v.z, 2);
}

TEST(Math, isPointInNormalDirection)
{
    ASSERT_TRUE(isPointInNormalDirection({0, 0, 0}, {0, 0, 1}, {0, 0, 0.5}));
    ASSERT_FALSE(isPointInNormalDirection({0, 0, 0}, {0, 0, 1}, {0, 0, -0.5}));
    
    ASSERT_TRUE(isPointInNormalDirection({1, 0, 0}, {1, 1, 1}, {3, 2, 1}));
    ASSERT_FALSE(isPointInNormalDirection({1, 0, 0}, {1, 1, 1}, {-3, -2, -1}));
    
    ASSERT_TRUE(isPointInNormalDirection({1, 1, 1}, {0, 0, 1}, {1, 1, 1}));
}