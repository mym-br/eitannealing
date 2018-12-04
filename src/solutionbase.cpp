#include "solutionbase.h"

void shuffler::addShufflerFeedback(const shuffleData &data, bool pos)
{
	if (pos) { // positive feedback
		if (data.swap)
			this->swapshuffleconsts[data.ncoef] /= 2;
		else
			this->shuffleConsts[data.ncoef] /= 2;
	}
	else {
		if (data.swap)
			this->swapshuffleconsts[data.ncoef]++;
		else
			this->shuffleConsts[data.ncoef]++;
	}
}