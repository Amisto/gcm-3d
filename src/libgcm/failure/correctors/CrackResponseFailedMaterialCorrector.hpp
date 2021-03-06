#ifndef GCM_CrackResponseFailedMaterialCorrector_H
#define GCM_CrackResponseFailedMaterialCorrector_H

#include "libgcm/failure/correctors/FailedMaterialCorrector.hpp"


namespace gcm {

    class CrackResponseFailedMaterialCorrector : public FailedMaterialCorrector {
    public:
		CrackResponseFailedMaterialCorrector();
		void applyCorrection(ICalcNode& node, const float tau);
    private:
        USE_LOGGER;
    };
}

#endif

